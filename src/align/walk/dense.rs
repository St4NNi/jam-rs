//! Finds walks: ungapped copies of k-mers too dense or too simple for anchors.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

// Only matches whose previous bases differ start a walk, grouped by strand and diagonal and
// extended with X-drop, so every copy of a repeat becomes one walk on its own diagonal.

use crate::align::BandConfig;
use crate::align::anchors::QueryTable;
use crate::align::{AlignmentError, Budgeted, Strand};
use crate::seeds::{base_code, kmers, reverse_complement};
use std::mem::size_of;

/// The best-scoring ungapped stretch around one or more match starts of a diagonal. `contig` is
/// the start in the oriented window (the reverse complement for the reverse strand) and `query`
/// the start in the query, below its length.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Walk {
    /// The contig strand of the walk.
    pub strand: Strand,
    /// The start in the oriented window.
    pub contig: u32,
    /// The start in the query.
    pub query: u32,
    /// The bases of the walk.
    pub length: u32,
}

/// The walks of an alignment window, and the reverse complement of the window that reverse
/// walks are read from (empty when no walk needs it).
#[derive(Debug, Default)]
pub struct DenseWalks {
    /// The walks of both strands.
    pub walks: Vec<Walk>,
    /// The reverse complement of the window.
    pub reverse: Vec<u8>,
}

/// A match start: a window k-mer that equals a query k-mer whose previous bases differ.
#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
struct Start {
    /// Oriented window position minus query position; modulo the query length when circular.
    diagonal: i64,
    contig: u32,
    query: u32,
}

/// The walks of every match start of a k-mer with more than `above` query occurrences, on both
/// strands, and the reverse complement of the window that reverse walks are read from. Starts,
/// that copy and the walks are priced against `budget` bytes before they are allocated.
pub fn dense_walks(
    table: &QueryTable,
    window: &[u8],
    above: usize,
    config: BandConfig,
    budget: usize,
) -> Result<DenseWalks, AlignmentError> {
    let a = usize::from(table.anchors().a);
    if window.len() < a || window.len() >= 1 << 32 {
        return Ok(DenseWalks::default());
    }
    if table.circular()
        && let Some(walks) = periodic_walks(table, window, budget)?
    {
        return Ok(walks);
    }
    // Count, price and only then hold the reverse complement and the starts.
    let copies = if table.circular() {
        window.len().div_ceil(table.query_length() as usize)
    } else {
        1
    };
    let price = |count: usize| {
        count
            .saturating_mul(size_of::<Start>() + copies.saturating_mul(size_of::<Walk>()))
            .saturating_add(window.len())
    };
    let mut count = 0usize;
    match_starts(table, window, above, &mut |_, _| count += 1);
    if price(count) > budget {
        return Err(AlignmentError::ChainingBudget(Budgeted::MatchStarts));
    }
    let reverse = reverse_complement(window);
    match_starts(table, &reverse, above, &mut |_, _| count += 1);
    if price(count) > budget {
        return Err(AlignmentError::ChainingBudget(Budgeted::MatchStarts));
    }
    let length = i64::from(table.query_length());
    let diagonal = |contig: u32, query: u32| {
        let diagonal = i64::from(contig) - i64::from(query);
        if table.circular() {
            diagonal.rem_euclid(length)
        } else {
            diagonal
        }
    };
    let mut walks = Vec::with_capacity(count.saturating_mul(copies));
    let mut starts = Vec::with_capacity(count);
    for (strand, sequence) in [(Strand::Forward, window), (Strand::Reverse, &reverse[..])] {
        starts.clear();
        match_starts(table, sequence, above, &mut |contig, query| {
            starts.push(Start {
                diagonal: diagonal(contig, query),
                contig,
                query,
            })
        });
        starts.sort_unstable();
        for group in starts.chunk_by(|a, b| a.diagonal == b.diagonal) {
            walk_diagonal(table, sequence, strand, group, config, &mut walks);
        }
    }
    Ok(DenseWalks { walks, reverse })
}

/// Whether a walk is long enough to count, as a chain must reach the minimum chain length.
pub fn long_enough(walk: &Walk, table: &QueryTable) -> bool {
    u64::from(walk.length) >= minimum_length(table)
}

/// The walks of a circular periodic query in a window of its repeat: one equivalent phase per
/// period, from every window position.
fn periodic_walks(
    table: &QueryTable,
    window: &[u8],
    budget: usize,
) -> Result<Option<DenseWalks>, AlignmentError> {
    let Some(periodic) = table.periodic_pair(window) else {
        return Ok(None);
    };
    let (period, phases) = (periodic.period, periodic.phases);
    let minimum = minimum_length(table) as usize;
    let positions = window.len().saturating_sub(minimum - 1);
    let count = positions.saturating_mul(phases.len());
    if count
        .saturating_mul(size_of::<Walk>())
        .saturating_add(window.len())
        > budget
    {
        return Err(AlignmentError::ChainingBudget(Budgeted::PeriodicWalks));
    }
    let mut walks = Vec::with_capacity(count);
    for phase in phases {
        for contig in 0..positions {
            walks.push(Walk {
                strand: phase.strand,
                contig: contig as u32,
                query: (contig as i64 - phase.diagonal).rem_euclid(period as i64) as u32,
                length: (window.len() - contig).min(table.query_length() as usize) as u32,
            });
        }
    }
    let reverse = reverse_complement(window);
    Ok(Some(DenseWalks { walks, reverse }))
}

/// Calls `visit(contig, query)` for every match start of an oriented window: a k-mer with more
/// than `above` query occurrences, at a query position whose preceding base differs.
fn match_starts(
    table: &QueryTable,
    oriented: &[u8],
    above: usize,
    visit: &mut dyn FnMut(u32, u32),
) {
    kmers(oriented, table.anchors().a, |kmer| {
        let (start, code) = (kmer.position, kmer.forward);
        if table.multiplicity(code) <= above {
            return;
        }
        let left = start
            .checked_sub(1)
            .and_then(|at| base_code(oriented[at]))
            .unwrap_or(4);
        for query in table.match_starts(code, left) {
            visit(start as u32, query);
        }
    });
}

/// Walks the starts of one diagonal; a start inside the previous walk gives none. A circular
/// walk over the whole query continues with the next copy when it is exact there.
fn walk_diagonal(
    table: &QueryTable,
    sequence: &[u8],
    strand: Strand,
    group: &[Start],
    config: BandConfig,
    walks: &mut Vec<Walk>,
) {
    let a = usize::from(table.anchors().a);
    let mut end = 0u32;
    for start in group {
        if start.contig < end {
            continue;
        }
        let mut next = *start;
        loop {
            let found = walk(table, sequence, strand, next, end, config);
            end = found.contig + found.length;
            // Short walks are kept too: joins link them across divergent pieces.
            walks.push(found);
            if !table.circular() || found.length != table.query_length() {
                break;
            }
            let (contig, query) = (end as usize, found.query as usize);
            if sequence.get(contig..contig + a) != table.scanned().get(query..query + a) {
                break;
            }
            next = Start {
                diagonal: start.diagonal,
                contig: end,
                query: found.query,
            };
        }
    }
}

/// Running state of one ungapped extension.
struct Extension {
    /// The scores of an equal and of a different base.
    matched: i64,
    mismatched: i64,
    /// The drop below the best score that ends a walk: a chaining gap of matches, so a walk
    /// crosses what a chain may bridge.
    drop: i64,
    /// The score after `length` bases, and the best score with the length that reached it.
    score: i64,
    best: i64,
    best_length: usize,
    length: usize,
}

impl Extension {
    fn new(table: &QueryTable, config: BandConfig) -> Self {
        let matched = i64::from(config.alignment.match_score);
        Self {
            matched,
            mismatched: i64::from(config.alignment.mismatch_score),
            drop: i64::from(table.chains().max_gap) * matched,
            score: 0,
            best: 0,
            best_length: 0,
            length: 0,
        }
    }

    /// Adds one base; false once the score dropped too far below the best.
    fn step(&mut self, same: bool) -> bool {
        self.score += if same { self.matched } else { self.mismatched };
        self.length += 1;
        if self.score > self.best {
            (self.best, self.best_length) = (self.score, self.length);
        }
        self.best - self.score <= self.drop
    }

    /// Adds a block of equal bases. Callers only do so at the best score, so the block end
    /// becomes the new best.
    fn equal_block(&mut self, length: usize) {
        self.score += self.matched * length as i64;
        self.length += length;
        (self.best, self.best_length) = (self.score, self.length);
    }

    /// Extends over equally long `contig` and `query` from their starts; false once the walk ends.
    fn forward(&mut self, contig: &[u8], query: &[u8]) -> bool {
        let mut at = 0;
        while at < contig.len() {
            if at + 8 <= contig.len()
                && self.score >= self.best
                && contig[at..at + 8] == query[at..at + 8]
            {
                self.equal_block(8);
                at += 8;
            } else if self.step(contig[at] == query[at]) {
                at += 1;
            } else {
                return false;
            }
        }
        true
    }

    /// Extends over equally long `contig` and `query` from their ends; false once the walk ends.
    fn backward(&mut self, contig: &[u8], query: &[u8]) -> bool {
        let mut at = contig.len();
        while at > 0 {
            if at >= 8 && self.score >= self.best && contig[at - 8..at] == query[at - 8..at] {
                self.equal_block(8);
                at -= 8;
            } else if self.step(contig[at - 1] == query[at - 1]) {
                at -= 1;
            } else {
                return false;
            }
        }
        true
    }
}

/// The walk around one match start, not reaching back before `previous_end`.
fn walk(
    table: &QueryTable,
    oriented: &[u8],
    strand: Strand,
    start: Start,
    previous_end: u32,
    config: BandConfig,
) -> Walk {
    let query = table.query();
    let length = query.len();
    let circular = table.circular();
    let (contig, from) = (start.contig as usize, start.query as usize);
    // A circular walk covers the query at most once; a linear one stays inside it. The query
    // is compared in its contiguous pieces: up to its end, then from its start when circular.
    let limit = (oriented.len() - contig).min(if circular { length } else { length - from });
    let mut right = Extension::new(table, config);
    let first = limit.min(length - from);
    if right.forward(
        &oriented[contig..contig + first],
        &query[from..from + first],
    ) {
        let rest = limit - first;
        right.forward(&oriented[contig + first..contig + limit], &query[..rest]);
    }
    let right = right.best_length;
    let room = contig - (previous_end as usize).min(contig);
    let limit = room.min(if circular { length - right } else { from });
    let mut left = Extension::new(table, config);
    let first = limit.min(from);
    if left.backward(
        &oriented[contig - first..contig],
        &query[from - first..from],
    ) {
        let rest = limit - first;
        left.backward(
            &oriented[contig - limit..contig - first],
            &query[length - rest..],
        );
    }
    let left = left.best_length;
    Walk {
        strand,
        contig: (contig - left) as u32,
        query: ((from + length - left) % length) as u32,
        length: (left + right) as u32,
    }
}

fn minimum_length(table: &QueryTable) -> u64 {
    table.chains().min_length(table.anchors().a)
}
