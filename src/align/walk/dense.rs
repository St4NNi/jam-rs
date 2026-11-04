//! Finds walks: ungapped copies of k-mers too dense or too simple for anchors.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

// Only matches whose previous bases differ start a walk, grouped by strand and diagonal and
// extended with X-drop, so every copy of a repeat becomes one walk on its own diagonal.

use crate::align::BandConfig;
use crate::align::anchors::QueryTable;
use crate::align::{AlignmentError, Budgeted, Strand, bases_match, known_bases};
use crate::seeds::{base_code, kmers, reverse_complement};
use rayon::prelude::*;
use std::mem::size_of;
use std::ops::Range;

/// Window bases from which the starts and walks of a window are found on all threads.
#[cfg(not(test))]
const PARALLEL_FROM: usize = 1 << 20;
/// As in builds, but tests may lower it.
#[cfg(test)]
pub(crate) static PARALLEL_FROM: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(1 << 20);

/// Window positions scanned for match starts by one task.
const SCAN_POSITIONS: usize = 1 << 16;

/// Match starts walked by one task, extended to whole diagonals.
const WALK_STARTS: usize = 1 << 14;

/// Whether a window of `bases` is scanned and walked on the threads of the current pool.
fn parallel(bases: usize) -> bool {
    #[cfg(test)]
    let from = PARALLEL_FROM.load(std::sync::atomic::Ordering::Relaxed);
    #[cfg(not(test))]
    let from = PARALLEL_FROM;
    bases >= from && rayon::current_num_threads() > 1
}

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
/// that copy and the walks are priced against `budget` bytes as the starts are found.
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
    let length = i64::from(table.query_length());
    let diagonal = |contig: u32, query: u32| {
        let diagonal = i64::from(contig) - i64::from(query);
        if table.circular() {
            diagonal.rem_euclid(length)
        } else {
            diagonal
        }
    };
    let over = |count: usize| {
        let needed = Some(price(count) as u64);
        AlignmentError::ChainingBudget(Budgeted::MatchStarts, needed)
    };
    let parallel = parallel(window.len());
    let (count, mut forward, mut backward) = if parallel {
        // Counted on all threads first, so the starts are held only once they fit.
        let positions = window.len() + 1 - a;
        let parts: Vec<Range<usize>> = (0..positions)
            .step_by(SCAN_POSITIONS)
            .map(|start| start..(start + SCAN_POSITIONS).min(positions))
            .collect();
        let counts = count_starts(table, window, above, &parts);
        let count = counts
            .iter()
            .map(|(forward, reverse)| forward + reverse)
            .sum();
        if price(count) > budget {
            return Err(over(count));
        }
        let (forward, backward) = fill_starts(table, window, above, &parts, &counts, &diagonal);
        (count, forward, backward)
    } else {
        // One scan finds the starts of both strands; the starts are priced as they come, and
        // none are kept once they pass the budget.
        let (mut count, mut forward, mut backward) = (0usize, Vec::new(), Vec::new());
        let positions = 0..window.len() + 1 - a;
        match_starts(
            table,
            window,
            positions,
            above,
            &mut |strand, contig, query| {
                count += 1;
                if price(count) <= budget {
                    let starts = if strand.is_reverse() {
                        &mut backward
                    } else {
                        &mut forward
                    };
                    starts.push(Start {
                        diagonal: diagonal(contig, query),
                        contig,
                        query,
                    });
                }
            },
        );
        (count, forward, backward)
    };
    if price(count) > budget {
        return Err(over(count));
    }
    let reverse = if backward.is_empty() {
        Vec::new()
    } else {
        reverse_complement(window)
    };
    let mut walks = Vec::with_capacity(count.saturating_mul(copies));
    for (strand, sequence, starts) in [
        (Strand::Forward, window, &mut forward),
        (Strand::Reverse, &reverse[..], &mut backward),
    ] {
        if !parallel {
            starts.sort_unstable();
            for group in starts.chunk_by(|a, b| a.diagonal == b.diagonal) {
                walk_diagonal(table, sequence, strand, group, config, &mut walks);
            }
            continue;
        }
        // Equal starts are equal in every field, so the unstable sort gives one order.
        starts.par_sort_unstable();
        let walked = WalkedStarts {
            table,
            sequence,
            strand,
            config,
            copies,
        };
        walked.walk_all(starts, &mut walks);
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
    let bytes = count
        .saturating_mul(size_of::<Walk>())
        .saturating_add(window.len());
    if bytes > budget {
        let needed = Some(bytes as u64);
        return Err(AlignmentError::ChainingBudget(
            Budgeted::PeriodicWalks,
            needed,
        ));
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

/// Calls `visit(strand, contig, query)` for every match start of the window on both strands: a
/// k-mer with more than `above` query occurrences, at a query position whose preceding base
/// differs. Reverse starts count positions in the reverse complement of the window.
fn match_starts(
    table: &QueryTable,
    window: &[u8],
    positions: Range<usize>,
    above: usize,
    visit: &mut dyn FnMut(Strand, u32, u32),
) {
    let a = usize::from(table.anchors().a);
    let end = (positions.end + a - 1).min(window.len());
    kmers(&window[positions.start..end], table.anchors().a, |kmer| {
        let start = positions.start + kmer.position;
        if !table.may_hold(kmer.forward) {
            return;
        }
        if table.multiplicity(kmer.forward) > above {
            let left = start
                .checked_sub(1)
                .and_then(|at| base_code(window[at]))
                .unwrap_or(4);
            for query in table.match_starts(kmer.forward, left) {
                visit(Strand::Forward, start as u32, query);
            }
        }
        if table.multiplicity(kmer.reverse) > above {
            // In the reverse complement this k-mer starts at `len - start - a`, after the
            // complement of the base that follows it here.
            let left = window
                .get(start + a)
                .and_then(|&base| base_code(base))
                .map_or(4, |code| 3 - code);
            let contig = (window.len() - start - a) as u32;
            for query in table.match_starts(kmer.reverse, left) {
                visit(Strand::Reverse, contig, query);
            }
        }
    });
}

/// The forward and reverse match starts of each part of the window positions, counted on all
/// threads.
fn count_starts(
    table: &QueryTable,
    window: &[u8],
    above: usize,
    parts: &[Range<usize>],
) -> Vec<(usize, usize)> {
    parts
        .par_iter()
        .map(|part| {
            let (mut forward, mut reverse) = (0, 0);
            match_starts(
                table,
                window,
                part.clone(),
                above,
                &mut |strand, _, _| match strand.is_reverse() {
                    true => reverse += 1,
                    false => forward += 1,
                },
            );
            (forward, reverse)
        })
        .collect()
}

/// The forward and reverse match starts of the window, found on all threads in the parts that
/// `counts` counted, each written to its own place, in the order one scan finds them.
fn fill_starts(
    table: &QueryTable,
    window: &[u8],
    above: usize,
    parts: &[Range<usize>],
    counts: &[(usize, usize)],
    diagonal: &(dyn Fn(u32, u32) -> i64 + Sync),
) -> (Vec<Start>, Vec<Start>) {
    let empty = Start {
        diagonal: 0,
        contig: 0,
        query: 0,
    };
    let (forwards, reverses): (usize, usize) =
        counts.iter().fold((0, 0), |(f, r), &(a, b)| (f + a, r + b));
    let (mut forward, mut backward) = (vec![empty; forwards], vec![empty; reverses]);
    let (mut forward_rest, mut backward_rest) = (&mut forward[..], &mut backward[..]);
    let mut places = Vec::with_capacity(parts.len());
    for (part, &(forwards, reverses)) in parts.iter().zip(counts) {
        let (own_forward, rest) = std::mem::take(&mut forward_rest).split_at_mut(forwards);
        forward_rest = rest;
        let (own_backward, rest) = std::mem::take(&mut backward_rest).split_at_mut(reverses);
        backward_rest = rest;
        places.push((part.clone(), own_forward, own_backward));
    }
    places
        .into_par_iter()
        .for_each(|(part, forward, backward)| {
            let (mut forwards, mut reverses) = (0, 0);
            match_starts(table, window, part, above, &mut |strand, contig, query| {
                let start = Start {
                    diagonal: diagonal(contig, query),
                    contig,
                    query,
                };
                if strand.is_reverse() {
                    backward[reverses] = start;
                    reverses += 1;
                } else {
                    forward[forwards] = start;
                    forwards += 1;
                }
            });
        });
    (forward, backward)
}

/// What the walks of the sorted starts of one strand read, for walking them on all threads.
struct WalkedStarts<'a> {
    table: &'a QueryTable,
    sequence: &'a [u8],
    strand: Strand,
    config: BandConfig,
    /// The most walks one start gives.
    copies: usize,
}

impl WalkedStarts<'_> {
    /// Walks the sorted `starts` in tasks of whole diagonals and adds the walks to `walks` in the
    /// order one loop adds them. Each task writes into the places of its starts; a task with more
    /// walks keeps the rest apart.
    fn walk_all(&self, starts: &[Start], walks: &mut Vec<Walk>) {
        let mut tasks = Vec::new();
        let mut rest = starts;
        while !rest.is_empty() {
            let mut end = WALK_STARTS.min(rest.len());
            while end < rest.len() && rest[end - 1].diagonal == rest[end].diagonal {
                end += 1;
            }
            let (task, after) = rest.split_at(end);
            tasks.push(task);
            rest = after;
        }
        let first = walks.len();
        let places = starts.len().saturating_mul(self.copies);
        let empty = Walk {
            strand: self.strand,
            contig: 0,
            query: 0,
            length: 0,
        };
        walks.resize(first + places, empty);
        let mut free = &mut walks[first..];
        let mut jobs = Vec::with_capacity(tasks.len());
        for task in tasks {
            let (own, after) = std::mem::take(&mut free).split_at_mut(task.len() * self.copies);
            free = after;
            jobs.push((task, own));
        }
        let ends: Vec<(usize, usize, Vec<Walk>)> = jobs
            .into_par_iter()
            .map(|(task, own)| {
                let (mut used, mut more) = (0, Vec::new());
                let (table, sequence, strand, config) =
                    (self.table, self.sequence, self.strand, self.config);
                let mut found = Vec::new();
                for group in task.chunk_by(|a, b| a.diagonal == b.diagonal) {
                    walk_diagonal(table, sequence, strand, group, config, &mut found);
                    for walk in found.drain(..) {
                        match own.get_mut(used) {
                            Some(place) => {
                                *place = walk;
                                used += 1;
                            }
                            None => more.push(walk),
                        }
                    }
                }
                (own.len(), used, more)
            })
            .collect();
        // Moves the walks of each task after those of the tasks before it.
        let (mut read, mut write) = (first, first);
        let spilled = ends.iter().any(|(_, _, more)| !more.is_empty());
        let mut kept = Vec::new();
        for (places, used, more) in ends {
            match spilled {
                true => {
                    kept.extend_from_slice(&walks[read..read + used]);
                    kept.extend(more);
                }
                false => walks.copy_within(read..read + used, write),
            }
            (read, write) = (read + places, write + used);
        }
        walks.truncate(match spilled {
            true => first,
            false => write,
        });
        walks.extend(kept);
    }
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

/// The bases an extension compares at once.
const BLOCK: usize = 8;

/// The number of equal bytes in two blocks of `BLOCK` bytes.
fn equal_bytes(contig: &[u8], query: &[u8]) -> u32 {
    let word = |bytes: &[u8]| u64::from_le_bytes(bytes.try_into().unwrap_or([0; BLOCK]));
    let differ = word(contig) ^ word(query);
    // The high bit of each byte is set when the byte is not zero.
    const LOW: u64 = 0x7f7f_7f7f_7f7f_7f7f;
    let nonzero = ((differ & LOW) + LOW) | differ;
    BLOCK as u32 - (nonzero & !LOW).count_ones()
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

    /// Adds a block of `BLOCK` bases when no base of it can reach a new best score or drop too
    /// far, as then only the number of equal bases counts. Returns whether it did.
    fn far_block(&mut self, contig: &[u8], query: &[u8]) -> bool {
        let bases = BLOCK as i64;
        if self.score + bases * self.matched > self.best
            || self.best - self.score - bases * self.mismatched > self.drop
            || !known_bases(query)
        {
            return false;
        }
        let equal = i64::from(equal_bytes(contig, query));
        self.score += equal * self.matched + (bases - equal) * self.mismatched;
        self.length += BLOCK;
        true
    }

    /// Extends over equally long `contig` and `query` from their starts; false once the walk ends.
    fn forward(&mut self, contig: &[u8], query: &[u8]) -> bool {
        let mut at = 0;
        while at < contig.len() {
            if at + BLOCK <= contig.len()
                && self.score >= self.best
                && contig[at..at + BLOCK] == query[at..at + BLOCK]
                && known_bases(&query[at..at + BLOCK])
            {
                self.equal_block(BLOCK);
                at += BLOCK;
            } else if at + BLOCK <= contig.len()
                && self.far_block(&contig[at..at + BLOCK], &query[at..at + BLOCK])
            {
                at += BLOCK;
            } else if self.step(bases_match(contig[at], query[at])) {
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
            let block = at.saturating_sub(BLOCK)..at;
            if at >= BLOCK
                && self.score >= self.best
                && contig[block.clone()] == query[block.clone()]
                && known_bases(&query[block.clone()])
            {
                self.equal_block(BLOCK);
                at -= BLOCK;
            } else if at >= BLOCK && self.far_block(&contig[block.clone()], &query[block]) {
                at -= BLOCK;
            } else if self.step(bases_match(contig[at - 1], query[at - 1])) {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{random_bases, rng};
    use rand::RngExt;

    /// The extension of one base at a time, without blocks.
    fn stepped(extension: &mut Extension, pairs: impl Iterator<Item = (u8, u8)>) -> bool {
        pairs
            .into_iter()
            .all(|(contig, query)| extension.step(contig == query))
    }

    #[test]
    fn blocks_extend_stepwise() {
        let mut random = rng(11);
        let table = QueryTable::new(
            &random_bases(&mut random, 100),
            false,
            crate::align::AnchorConfig::default(),
            crate::align::ChainConfig {
                max_gap: 200,
                ..crate::align::ChainConfig::default()
            },
        )
        .unwrap();
        let config = BandConfig::default();
        for round in 0..400 {
            let length = random.random_range(0..1500);
            let contig = random_bases(&mut random, length);
            // From a copy to unrelated bases, so walks end anywhere.
            let changes = [0.0, 0.05, 0.3, 0.75][round % 4];
            let query: Vec<u8> = contig
                .iter()
                .map(|&base| {
                    if random.random_bool(changes) {
                        b"ACGT"[random.random_range(0..4)]
                    } else {
                        base
                    }
                })
                .collect();
            let state =
                |extension: &Extension| (extension.score, extension.best, extension.best_length);
            let (mut fast, mut slow) = (
                Extension::new(&table, config),
                Extension::new(&table, config),
            );
            let pairs = contig.iter().copied().zip(query.iter().copied());
            let went = fast.forward(&contig, &query);
            assert_eq!(went, stepped(&mut slow, pairs.clone()));
            assert_eq!(fast.best_length, slow.best_length);
            if went {
                assert_eq!(state(&fast), state(&slow));
            }
            let (mut fast, mut slow) = (
                Extension::new(&table, config),
                Extension::new(&table, config),
            );
            let went = fast.backward(&contig, &query);
            assert_eq!(went, stepped(&mut slow, pairs.rev()));
            assert_eq!(
                (went, fast.best, fast.best_length),
                (went, slow.best, slow.best_length)
            );
        }
    }

    #[test]
    fn equal_bytes_counted() {
        assert_eq!(equal_bytes(b"ACGTACGT", b"ACGTACGT"), 8);
        assert_eq!(equal_bytes(b"ACGTACGT", b"TGCATGCA"), 0);
        assert_eq!(equal_bytes(b"AAAAAAAA", b"AxAAAyAA"), 6);
        assert_eq!(
            equal_bytes(&[0, 255, 128, 1, 0, 0, 0, 0], &[0, 255, 0, 1, 0, 0, 0, 1]),
            6
        );
    }
}
