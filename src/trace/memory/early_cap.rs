//! Streaming top-k caps for occurrence and per-contig work.
//!
//! These helpers retain at most the configured number of records while an
//! occurrence stream is decoded.  They never materialize a complete posting
//! and never sort a complete posting before applying its cap.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap};

/// A deterministic ranking where larger scores win and smaller ties win.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct HeapRank {
    score: i64,
    tie_break: u64,
}

impl HeapRank {
    #[must_use]
    pub const fn new(score: i64, tie_break: u64) -> Self {
        Self { score, tie_break }
    }

    #[must_use]
    pub const fn score(self) -> i64 {
        self.score
    }

    #[must_use]
    pub const fn tie_break(self) -> u64 {
        self.tie_break
    }
}

#[derive(Debug)]
struct Ranked<T> {
    rank: HeapRank,
    item: T,
}

impl<T> PartialEq for Ranked<T> {
    fn eq(&self, other: &Self) -> bool {
        self.rank == other.rank
    }
}

impl<T> Eq for Ranked<T> {}

impl<T> Ord for Ranked<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        // BinaryHeap's root is the "worst" retained record: lower score,
        // then larger tie-break.  This lets a new record replace the root in
        // O(log k) while keeping the retained set bounded.
        other
            .rank
            .score
            .cmp(&self.rank.score)
            .then_with(|| self.rank.tie_break.cmp(&other.rank.tie_break))
    }
}

impl<T> PartialOrd for Ranked<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// A fixed-size max-quality set backed by a binary heap.
#[derive(Debug)]
pub struct FixedTopK<T> {
    capacity: usize,
    heap: BinaryHeap<Ranked<T>>,
}

impl<T> FixedTopK<T> {
    /// Creates a top-k set.  A zero capacity is valid and retains nothing.
    #[must_use]
    pub fn new(capacity: usize) -> Self {
        Self {
            capacity,
            heap: BinaryHeap::with_capacity(capacity),
        }
    }

    #[must_use]
    pub const fn capacity(&self) -> usize {
        self.capacity
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.heap.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.heap.is_empty()
    }

    #[must_use]
    pub fn is_full(&self) -> bool {
        self.len() == self.capacity
    }

    /// Offers one item and returns whether it is retained after the offer.
    pub fn push(&mut self, item: T, rank: HeapRank) -> bool {
        if self.capacity == 0 {
            return false;
        }
        if self.heap.len() < self.capacity {
            self.heap.push(Ranked { rank, item });
            return true;
        }

        let should_replace = self
            .heap
            .peek()
            .is_some_and(|worst| better(rank, worst.rank));
        if !should_replace {
            return false;
        }
        let _ = self.heap.pop();
        self.heap.push(Ranked { rank, item });
        true
    }

    /// Removes all retained records in best-first deterministic order.
    pub fn into_sorted(self) -> Vec<T> {
        let mut records: Vec<_> = self.heap.into_vec();
        records.sort_unstable_by(|left, right| compare_best_first(left.rank, right.rank));
        records.into_iter().map(|record| record.item).collect()
    }

    /// Removes all retained records with their ranks in best-first order.
    pub fn into_ranked_sorted(self) -> Vec<(HeapRank, T)> {
        let mut records: Vec<_> = self.heap.into_vec();
        records.sort_unstable_by(|left, right| compare_best_first(left.rank, right.rank));
        records
            .into_iter()
            .map(|record| (record.rank, record.item))
            .collect()
    }
}

/// A bounded collection of fixed-size heaps keyed by contig ID.
///
/// `max_contigs` bounds the number of maps allocated.  New contigs after that
/// bound are ignored deterministically, while existing contigs continue to
/// retain their best records.  This is useful when the caller already has a
/// bounded chain/contig admission policy.
#[derive(Debug)]
pub struct PerContigTopK<T> {
    per_contig_capacity: usize,
    max_contigs: usize,
    heaps: BTreeMap<u32, FixedTopK<T>>,
}

impl<T> PerContigTopK<T> {
    #[must_use]
    pub fn new(per_contig_capacity: usize, max_contigs: usize) -> Self {
        Self {
            per_contig_capacity,
            max_contigs,
            heaps: BTreeMap::new(),
        }
    }

    #[must_use]
    pub const fn per_contig_capacity(&self) -> usize {
        self.per_contig_capacity
    }

    #[must_use]
    pub const fn max_contigs(&self) -> usize {
        self.max_contigs
    }

    #[must_use]
    pub fn contig_count(&self) -> usize {
        self.heaps.len()
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.heaps.values().map(FixedTopK::len).sum()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.heaps.values().all(FixedTopK::is_empty)
    }

    /// Offers one occurrence to its contig heap.
    pub fn push(&mut self, contig_id: u32, item: T, rank: HeapRank) -> bool {
        if let Some(heap) = self.heaps.get_mut(&contig_id) {
            return heap.push(item, rank);
        }
        if self.heaps.len() >= self.max_contigs || self.per_contig_capacity == 0 {
            return false;
        }
        let mut heap = FixedTopK::new(self.per_contig_capacity);
        let retained = heap.push(item, rank);
        self.heaps.insert(contig_id, heap);
        retained
    }

    /// Removes all records in ascending contig-ID and best-first rank order.
    pub fn into_sorted(self) -> BTreeMap<u32, Vec<T>> {
        self.heaps
            .into_iter()
            .map(|(contig, heap)| (contig, heap.into_sorted()))
            .collect()
    }

    /// Removes all ranked records in ascending contig-ID and best-first order.
    pub fn into_ranked_sorted(self) -> BTreeMap<u32, Vec<(HeapRank, T)>> {
        self.heaps
            .into_iter()
            .map(|(contig, heap)| (contig, heap.into_ranked_sorted()))
            .collect()
    }
}

/// Retains a bounded top-k set from an occurrence iterator.
pub fn retain_top_k<I, T, F>(occurrences: I, capacity: usize, mut rank: F) -> Vec<T>
where
    I: IntoIterator<Item = T>,
    F: FnMut(&T) -> HeapRank,
{
    let mut retained = FixedTopK::new(capacity);
    for occurrence in occurrences {
        let occurrence_rank = rank(&occurrence);
        let _ = retained.push(occurrence, occurrence_rank);
    }
    retained.into_sorted()
}

/// Retains bounded top-k records independently for each contig.
pub fn retain_per_contig_top_k<I, T, C, F>(
    occurrences: I,
    per_contig_capacity: usize,
    max_contigs: usize,
    mut contig: C,
    mut rank: F,
) -> BTreeMap<u32, Vec<T>>
where
    I: IntoIterator<Item = T>,
    C: FnMut(&T) -> u32,
    F: FnMut(&T) -> HeapRank,
{
    let mut retained = PerContigTopK::new(per_contig_capacity, max_contigs);
    for occurrence in occurrences {
        let contig_id = contig(&occurrence);
        let occurrence_rank = rank(&occurrence);
        let _ = retained.push(contig_id, occurrence, occurrence_rank);
    }
    retained.into_sorted()
}

fn better(left: HeapRank, right: HeapRank) -> bool {
    left.score > right.score || (left.score == right.score && left.tie_break < right.tie_break)
}

fn compare_best_first(left: HeapRank, right: HeapRank) -> Ordering {
    right
        .score
        .cmp(&left.score)
        .then_with(|| left.tie_break.cmp(&right.tie_break))
}
