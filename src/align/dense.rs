//! Copies in dense repeats, found without expanding all their anchors. A window k-mer with more
//! query occurrences than the chaining budget allows gives no anchor per occurrence. Instead only
//! maximal exact matches are started: at every oriented window position, occurrences whose
//! preceding query base equals the preceding window base continue a match that started earlier
//! and are skipped. The starts are grouped by strand and diagonal and walked with ungapped X-drop
//! extension, so every copy of a repeat becomes one ungapped alignment on its own diagonal.

use super::AlignmentError;
use super::anchors::{QueryTable, base_code, kmers};
use super::chain_alignment::ChainConfig;
use super::complement;
use std::mem::size_of;

/// The best-scoring ungapped stretch around one or more match starts of a diagonal. `contig` is
/// the start in the oriented window (the reverse complement for the reverse strand) and `query`
/// the start in the query, below its length.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Walk {
    pub reverse: bool,
    pub contig: u32,
    pub query: u32,
    pub length: u32,
}

#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
struct Start {
    /// Oriented window position minus query position; modulo the query length when circular.
    line: i64,
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
    config: ChainConfig,
    budget: usize,
) -> Result<(Vec<Walk>, Vec<u8>), AlignmentError> {
    let k = usize::from(table.config().k);
    if window.len() < k || window.len() >= 1 << 32 {
        return Ok((Vec::new(), Vec::new()));
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
        return Err(AlignmentError::ChainingBudget("dense match starts"));
    }
    let reverse: Vec<u8> = window.iter().rev().map(|base| complement(*base)).collect();
    match_starts(table, &reverse, above, &mut |_, _| count += 1);
    if price(count) > budget {
        return Err(AlignmentError::ChainingBudget("dense match starts"));
    }
    let length = i64::from(table.query_length());
    let line = |contig: u32, query: u32| {
        let line = i64::from(contig) - i64::from(query);
        if table.circular() {
            line.rem_euclid(length)
        } else {
            line
        }
    };
    let mut walks = Vec::with_capacity(count.saturating_mul(copies));
    let mut starts = Vec::with_capacity(count);
    for (reverse_strand, sequence) in [(false, window), (true, &reverse[..])] {
        starts.clear();
        match_starts(table, sequence, above, &mut |contig, query| {
            starts.push(Start {
                line: line(contig, query),
                contig,
                query,
            })
        });
        starts.sort_unstable();
        for group in starts.chunk_by(|a, b| a.line == b.line) {
            walk_line(table, sequence, reverse_strand, group, config, &mut walks);
        }
    }
    Ok((walks, reverse))
}

/// The walks of a circular periodic query in a window of its repeat: one equivalent phase per
/// period, from every window position.
#[allow(clippy::type_complexity)]
fn periodic_walks(
    table: &QueryTable,
    window: &[u8],
    budget: usize,
) -> Result<Option<(Vec<Walk>, Vec<u8>)>, AlignmentError> {
    let Some((period, phases)) = table.periodic_pair(window) else {
        return Ok(None);
    };
    let minimum = minimum_length(table) as usize;
    let positions = window.len().saturating_sub(minimum - 1);
    let count = positions.saturating_mul(phases.len());
    if count
        .saturating_mul(size_of::<Walk>())
        .saturating_add(window.len())
        > budget
    {
        return Err(AlignmentError::ChainingBudget("circular periodic walks"));
    }
    let mut walks = Vec::with_capacity(count);
    for (reverse, phase) in phases {
        for contig in 0..positions {
            walks.push(Walk {
                reverse,
                contig: contig as u32,
                query: (contig as i64 - phase).rem_euclid(period as i64) as u32,
                length: (window.len() - contig).min(table.query_length() as usize) as u32,
            });
        }
    }
    let reverse = window.iter().rev().map(|base| complement(*base)).collect();
    Ok(Some((walks, reverse)))
}

/// Calls `visit(contig, query)` for every match start of an oriented window: a k-mer with more
/// than `above` query occurrences, at a query position whose preceding base differs.
fn match_starts(
    table: &QueryTable,
    oriented: &[u8],
    above: usize,
    visit: &mut dyn FnMut(u32, u32),
) {
    kmers(oriented, table.config().k, |start, code, _| {
        if table.multiplicity(code) <= above {
            return;
        }
        let left = start
            .checked_sub(1)
            .and_then(|at| base_code(oriented[at]))
            .map_or(4, |code| code as u8);
        for query in table.match_starts(code, left) {
            visit(start as u32, query);
        }
    });
}

/// Walks the starts of one diagonal; a start inside the previous walk gives none. A circular
/// walk over the whole query continues with the next copy when it is exact there.
fn walk_line(
    table: &QueryTable,
    sequence: &[u8],
    reverse: bool,
    group: &[Start],
    config: ChainConfig,
    walks: &mut Vec<Walk>,
) {
    let k = usize::from(table.config().k);
    let mut end = 0u32;
    for start in group {
        if start.contig < end {
            continue;
        }
        let mut next = *start;
        loop {
            let found = walk(table, sequence, reverse, next, end, config);
            end = found.contig + found.length;
            // Shorter walks are kept as seeds of gapped joins.
            walks.push(found);
            if !table.circular() || found.length != table.query_length() {
                break;
            }
            let (contig, query) = (end as usize, found.query as usize);
            if sequence.get(contig..contig + k) != table.scanned().get(query..query + k) {
                break;
            }
            next = Start {
                line: start.line,
                contig: end,
                query: found.query,
            };
        }
    }
}

/// Running state of one ungapped extension.
struct Extension {
    matched: i64,
    mismatched: i64,
    /// The drop below the best score that ends a walk: a chaining gap of matches, so a walk
    /// crosses what a chain may bridge.
    drop: i64,
    score: i64,
    best: i64,
    best_length: usize,
    length: usize,
}

impl Extension {
    fn new(table: &QueryTable, config: ChainConfig) -> Self {
        let matched = i64::from(config.scores.match_score);
        Self {
            matched,
            mismatched: i64::from(config.scores.mismatch_score),
            drop: i64::from(table.config().max_gap) * matched,
            score: 0,
            best: 0,
            best_length: 0,
            length: 0,
        }
    }

    fn step(&mut self, same: bool) -> bool {
        self.score += if same { self.matched } else { self.mismatched };
        self.length += 1;
        if self.score > self.best {
            (self.best, self.best_length) = (self.score, self.length);
        }
        self.best - self.score <= self.drop
    }

    fn equal_block(&mut self, length: usize) {
        self.score += self.matched * length as i64;
        self.length += length;
        (self.best, self.best_length) = (self.score, self.length);
    }

    /// Extends over equally long `a` and `b` from their starts; false once the walk ends.
    fn forward(&mut self, a: &[u8], b: &[u8]) -> bool {
        let mut at = 0;
        while at < a.len() {
            if at + 8 <= a.len() && self.score >= self.best && a[at..at + 8] == b[at..at + 8] {
                self.equal_block(8);
                at += 8;
            } else if self.step(a[at] == b[at]) {
                at += 1;
            } else {
                return false;
            }
        }
        true
    }

    /// Extends over equally long `a` and `b` from their ends; false once the walk ends.
    fn backward(&mut self, a: &[u8], b: &[u8]) -> bool {
        let mut at = a.len();
        while at > 0 {
            if at >= 8 && self.score >= self.best && a[at - 8..at] == b[at - 8..at] {
                self.equal_block(8);
                at -= 8;
            } else if self.step(a[at - 1] == b[at - 1]) {
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
    reverse: bool,
    start: Start,
    previous_end: u32,
    config: ChainConfig,
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
        reverse,
        contig: (contig - left) as u32,
        query: ((from + length - left) % length) as u32,
        length: (left + right) as u32,
    }
}

/// Whether a walk is long enough to count, as a chain must reach the minimum chain length.
pub fn long_enough(walk: &Walk, table: &QueryTable) -> bool {
    u64::from(walk.length) >= minimum_length(table)
}

fn minimum_length(table: &QueryTable) -> u64 {
    let anchors = table.config();
    u64::from(anchors.k)
        .max(anchors.min_score.max(1) as u64)
        .max(u64::from(anchors.k) + u64::from(anchors.min_anchors) - 1)
}
