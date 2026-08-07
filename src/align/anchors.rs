//! Exact anchors between the query and an alignment window, found through the query table.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

// Every `stride`-th window position is looked up, so an identical stretch of a + stride - 1
// bases always gives an anchor. A circular query also holds the k-mers across its origin.

use super::chains::ChainConfig;
use super::{AlignmentError, Strand};
use crate::seeds::{base_code, complement, kmers};
use jamhash::jamhash_u64;

/// How anchors are found.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnchorConfig {
    /// The anchor length a.
    pub a: u8,
    /// Window positions looked up: every `stride`-th, a power of two.
    pub stride: u8,
    /// Query k-mers whose base entropy (bits) is lower give no anchors; dense walks cover them.
    pub min_entropy: f64,
}

impl Default for AnchorConfig {
    fn default() -> Self {
        Self {
            a: 15,
            stride: 4,
            min_entropy: 1.0,
        }
    }
}

impl AnchorConfig {
    /// Checks a (8 to 31), a power-of-two stride and a finite entropy limit.
    pub fn validate(self) -> Result<(), AlignmentError> {
        if !(8..=31).contains(&self.a)
            || self.stride == 0
            || !self.stride.is_power_of_two()
            || !self.min_entropy.is_finite()
        {
            return Err(AlignmentError::InvalidAnchorConfig);
        }
        Ok(())
    }
}

/// A shared k-mer. `contig` is the oriented window position: the k-mer start in the window for
/// the forward strand, and in the reverse complement of the window for the reverse strand.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct Anchor {
    /// The contig strand of the match.
    pub strand: Strand,
    /// The oriented window position.
    pub contig: u32,
    /// The query position.
    pub query: u32,
}

/// Every k-mer of the query on both strands, keyed by its 2-bit code. The forward occurrences of
/// a k-mer come first, ordered by the code of the query base before them (4 for none or another
/// byte), then its reverse occurrences.
pub struct QueryTable {
    anchors: AnchorConfig,
    chains: ChainConfig,
    query_length: u32,
    circular: bool,
    period: Option<usize>,
    slots: Vec<Slot>,
    /// One bit for each key, so that most window k-mers that the query lacks are rejected without
    /// a probe into `slots`, which is many times larger. It never rejects a key of the table.
    filter: Vec<u64>,
    /// The most occurrences of any one k-mer.
    max_occurrences: usize,
    occurrences: Vec<u32>,
    /// For each occurrence: the code of the base before a forward occurrence, 5 for a reverse one.
    left: Vec<u8>,
    query: Vec<u8>,
}

impl QueryTable {
    /// The table of the query k-mers that pass the entropy limit.
    pub fn new(
        query: &[u8],
        circular: bool,
        anchors: AnchorConfig,
        chains: ChainConfig,
    ) -> Result<Self, AlignmentError> {
        Self::with_kmers(query, circular, anchors, chains, |forward| {
            base_entropy(forward, anchors.a) >= anchors.min_entropy
        })
    }

    /// The table of the k-mers that `new` leaves out for their low entropy, all left to dense
    /// walks. Short conserved blocks can rebuild a longer copy across a clustered difference.
    pub fn low_complexity(
        query: &[u8],
        circular: bool,
        anchors: AnchorConfig,
        chains: ChainConfig,
    ) -> Result<Self, AlignmentError> {
        Self::with_kmers(query, circular, anchors, chains, |forward| {
            base_entropy(forward, anchors.a) < anchors.min_entropy
        })
    }

    fn with_kmers(
        query: &[u8],
        circular: bool,
        anchors: AnchorConfig,
        chains: ChainConfig,
        keep: impl Fn(u64) -> bool,
    ) -> Result<Self, AlignmentError> {
        anchors.validate()?;
        chains.validate()?;
        let query_length = u32::try_from(query.len())
            .ok()
            .filter(|&length| length < REVERSE)
            .ok_or(AlignmentError::QueryTooLong)?;
        let a = usize::from(anchors.a);
        let period = period(query, circular, a);
        let mut scanned = query.to_vec();
        if circular && query.len() >= a {
            scanned.extend_from_slice(&query[..a - 1]);
        }
        let before = |start: usize| match start.checked_sub(1) {
            Some(at) => query[at],
            None if circular => query[query.len() - 1],
            None => b'N',
        };
        // (k-mer, left code, occurrence): sorting puts the forward occurrences of a k-mer first,
        // by left code, and its reverse ones (code 5) last.
        let mut pairs: Vec<(u64, u8, u32)> = Vec::new();
        kmers(&scanned, anchors.a, |kmer| {
            let (start, forward, reverse) = (kmer.position, kmer.forward, kmer.reverse);
            if start < query.len() && keep(forward) {
                let left = base_code(before(start)).unwrap_or(4);
                pairs.push((forward, left, start as u32));
                pairs.push((reverse, 5, start as u32 | REVERSE));
            }
        });
        pairs.sort_unstable();
        let distinct = pairs.chunk_by(|a, b| a.0 == b.0).count();
        let max_occurrences = pairs
            .chunk_by(|a, b| a.0 == b.0)
            .map(<[_]>::len)
            .max()
            .unwrap_or(0);
        let slots = (distinct * 2).next_power_of_two().max(16);
        let empty = Slot {
            key: EMPTY,
            start: 0,
            count: 0,
        };
        let mut table = vec![empty; slots];
        let mut filter = vec![0u64; ((distinct * 8).next_power_of_two().max(1024)) / 64];
        let mut occurrences = Vec::with_capacity(pairs.len());
        let mut left = Vec::with_capacity(pairs.len());
        for group in pairs.chunk_by(|a, b| a.0 == b.0) {
            let hash = jamhash_u64(group[0].0);
            let mut slot = (hash as usize) & (slots - 1);
            while table[slot].key != EMPTY {
                slot = (slot + 1) & (slots - 1);
            }
            table[slot] = Slot {
                key: group[0].0,
                start: occurrences.len() as u32,
                count: group.len() as u32,
            };
            let bit = filter_bit(hash, filter.len());
            filter[bit / 64] |= 1 << (bit % 64);
            occurrences.extend(group.iter().map(|pair| pair.2));
            left.extend(group.iter().map(|pair| pair.1));
        }
        Ok(Self {
            anchors,
            chains,
            query_length,
            circular,
            period,
            slots: table,
            filter,
            max_occurrences,
            occurrences,
            left,
            query: scanned,
        })
    }

    /// The occurrences of `key`: query positions, with `REVERSE` set for the reverse strand.
    fn get(&self, key: u64) -> &[u32] {
        &self.occurrences[self.occurrence_range(key)]
    }

    /// The range of the occurrences of `key`, empty when the query lacks it.
    fn occurrence_range(&self, key: u64) -> std::ops::Range<usize> {
        let hash = jamhash_u64(key);
        let bit = filter_bit(hash, self.filter.len());
        if self.filter[bit / 64] & (1 << (bit % 64)) == 0 {
            return 0..0;
        }
        let mask = self.slots.len() - 1;
        let mut slot = (hash as usize) & mask;
        loop {
            let found = self.slots[slot];
            if found.key == key {
                let start = found.start as usize;
                return start..start + found.count as usize;
            }
            if found.key == EMPTY {
                return 0..0;
            }
            slot = (slot + 1) & mask;
        }
    }

    /// Occurrences of `key` on both strands.
    pub(super) fn multiplicity(&self, key: u64) -> usize {
        self.occurrence_range(key).len()
    }

    /// Forward query positions of `key` that do not follow a base of code `left`, so a match
    /// with a window position whose preceding base has that code cannot extend to the left.
    /// Code 4 (none or another byte) excludes nothing.
    pub(super) fn match_starts(&self, key: u64, left: u8) -> impl Iterator<Item = u32> + '_ {
        let range = self.occurrence_range(key);
        let codes = &self.left[range.clone()];
        let forward = codes.partition_point(|code| *code < 5);
        let (low, high) = if left < 4 {
            (
                codes.partition_point(|code| *code < left),
                codes.partition_point(|code| *code <= left),
            )
        } else {
            (forward, forward)
        };
        let occurrences = &self.occurrences[range];
        occurrences[..low]
            .iter()
            .chain(&occurrences[high..forward])
            .copied()
    }

    /// The query as scanned: for a circular query, followed by its first a - 1 bases.
    pub(super) fn scanned(&self) -> &[u8] {
        &self.query
    }

    /// Appends `span` query bases from `start` on to `bases`; on a circular query they run on
    /// over the origin.
    pub(super) fn copy_bases(
        &self,
        start: u64,
        span: u64,
        bases: &mut Vec<u8>,
    ) -> Result<(), AlignmentError> {
        let query = self.query();
        let (start, span) = (start as usize, span as usize);
        if !self.circular {
            let piece = query.get(start..start + span);
            bases.extend_from_slice(piece.ok_or(AlignmentError::EndpointOutsideWindow)?);
            return Ok(());
        }
        if query.is_empty() {
            return Err(AlignmentError::EmptyQuery);
        }
        let (mut position, mut left) = (start % query.len(), span);
        while left > 0 {
            let count = left.min(query.len() - position);
            bases.extend_from_slice(&query[position..position + count]);
            (position, left) = (0, left - count);
        }
        Ok(())
    }

    /// The query bases, without the bases a circular query repeats.
    pub fn query(&self) -> &[u8] {
        &self.query[..self.query_length as usize]
    }

    /// The number of query bases.
    pub fn query_length(&self) -> u32 {
        self.query_length
    }

    /// Whether the query is circular.
    pub fn circular(&self) -> bool {
        self.circular
    }

    /// The anchor settings the table was built with.
    pub fn anchors(&self) -> AnchorConfig {
        self.anchors
    }

    /// The chaining limits the table was built with.
    pub fn chains(&self) -> ChainConfig {
        self.chains
    }

    /// The most occurrences of any one k-mer.
    pub fn max_occurrences(&self) -> usize {
        self.max_occurrences
    }

    /// When the query is a repeat of a period of at most a bases and `window` is a copy of the
    /// same repeat: the period and, per strand, the diagonal of one exact copy.
    pub(super) fn periodic_pair(&self, window: &[u8]) -> Option<Periodic> {
        if self.occurrences.is_empty() || self.query.is_empty() || window.is_empty() {
            return None;
        }
        let period = self.period?;
        let query = self.query();
        if window.iter().any(|base| base_code(*base).is_none()) {
            return None;
        }
        let mut phases = Vec::new();
        for strand in [Strand::Forward, Strand::Reverse] {
            let base = |at: usize| {
                if strand.is_reverse() {
                    complement(window[window.len() - 1 - at]).to_ascii_uppercase()
                } else {
                    window[at].to_ascii_uppercase()
                }
            };
            for shift in 0..period {
                if (0..window.len())
                    .all(|at| base(at) == query[(at + shift) % period].to_ascii_uppercase())
                {
                    let diagonal = -(shift as i64);
                    phases.push(Phase { strand, diagonal });
                    break;
                }
            }
        }
        (!phases.is_empty()).then_some(Periodic { period, phases })
    }
}

/// Appends the anchors of `window`, sorted by strand, contig and query position. K-mers with
/// more than `above` query occurrences give none, when it is set; `expected` anchors are
/// reserved first.
pub(super) fn find_anchors(
    table: &QueryTable,
    window: &[u8],
    anchors: &mut Vec<Anchor>,
    above: Option<usize>,
    expected: usize,
) {
    let a = u32::from(table.anchors.a);
    let stride = u32::from(table.anchors.stride);
    let length = window.len() as u32;
    let first = anchors.len();
    if anchors.capacity() < first + expected {
        anchors.reserve_exact(expected);
    }
    kmers(window, table.anchors.a, |kmer| {
        let (start, forward) = (kmer.position, kmer.forward);
        let start = start as u32;
        if !start.is_multiple_of(stride) {
            return;
        }
        let occurrences = table.get(forward);
        if above.is_some_and(|threshold| occurrences.len() > threshold) {
            return;
        }
        for &occurrence in occurrences {
            let (strand, contig) = match occurrence & REVERSE {
                0 => (Strand::Forward, start),
                _ => (Strand::Reverse, length - a - start),
            };
            anchors.push(Anchor {
                strand,
                contig,
                query: occurrence & !REVERSE,
            });
        }
    });
    anchors[first..].sort_unstable();
}

/// The anchors that the looked-up positions of `window` give.
pub(super) fn count_anchors(table: &QueryTable, window: &[u8]) -> usize {
    let stride = usize::from(table.anchors.stride);
    let mut count = 0usize;
    kmers(window, table.anchors.a, |kmer| {
        let (start, forward) = (kmer.position, kmer.forward);
        if start % stride == 0 {
            count = count.saturating_add(table.multiplicity(forward));
        }
    });
    count
}

/// How many looked-up positions of `window` have k-mers of each query occurrence count.
pub(super) fn occurrence_counts(
    table: &QueryTable,
    window: &[u8],
) -> std::collections::BTreeMap<usize, usize> {
    let stride = usize::from(table.anchors.stride);
    let mut positions = std::collections::BTreeMap::<usize, usize>::new();
    kmers(window, table.anchors.a, |kmer| {
        let (start, forward) = (kmer.position, kmer.forward);
        if start % stride == 0 {
            let occurrences = table.multiplicity(forward);
            if occurrences > 0 {
                *positions.entry(occurrences).or_default() += 1;
            }
        }
    });
    positions
}

/// A query that repeats a short period, and an alignment window that copies the repeat.
pub(super) struct Periodic {
    /// The period in bases.
    pub(super) period: usize,
    /// The diagonal of one exact copy per strand that has one.
    pub(super) phases: Vec<Phase>,
}

/// The diagonal of one exact copy of a periodic query on one strand.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct Phase {
    /// The contig strand of the copy.
    pub(super) strand: Strand,
    /// Oriented window position minus query position, at most zero.
    pub(super) diagonal: i64,
}

/// One slot of the open-addressing table: a k-mer and where its occurrences lie.
#[derive(Clone, Copy)]
struct Slot {
    key: u64,
    start: u32,
    count: u32,
}

/// The key of an empty slot; no k-mer of at most 31 bases has it.
const EMPTY: u64 = u64::MAX;
/// The bit of an occurrence that marks the reverse strand; a query stays shorter than it.
const REVERSE: u32 = 1 << 31;

/// The bit of a key with hash `hash` in a filter of `words` 64-bit words (a power of two).
#[inline(always)]
fn filter_bit(hash: u64, words: usize) -> usize {
    (hash >> (64 - (words * 64).trailing_zeros())) as usize
}

/// The Shannon entropy in bits of the bases of a 2-bit k-mer; k-mers below the entropy limit
/// are left to dense walks.
fn base_entropy(kmer: u64, k: u8) -> f64 {
    let mut counts = [0u8; 4];
    for at in 0..k {
        counts[((kmer >> (2 * at)) & 3) as usize] += 1;
    }
    counts
        .iter()
        .filter(|&&count| count > 0)
        .map(|&count| {
            let share = f64::from(count) / f64::from(k);
            -share * share.log2()
        })
        .sum()
}

/// The shortest period of at most a bases of a query of plain bases; a circular query must hold
/// a whole number of periods.
fn period(query: &[u8], circular: bool, a: usize) -> Option<usize> {
    if query.iter().any(|base| base_code(*base).is_none()) {
        return None;
    }
    (1..=a.min(query.len())).find(|&candidate| {
        (!circular || query.len().is_multiple_of(candidate))
            && query
                .iter()
                .enumerate()
                .all(|(at, base)| base.eq_ignore_ascii_case(&query[at % candidate]))
    })
}
