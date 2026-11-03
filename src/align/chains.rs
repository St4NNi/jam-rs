//! Chains: anchors of one strand in order along nearly the same diagonal.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::anchors::{Anchor, QueryTable, count_anchors, find_anchors, occurrence_counts};
use super::{AlignmentError, Budgeted, Strand};

/// Colinear anchors of one strand, in increasing oriented contig position.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Chain {
    /// The contig strand of all anchors.
    pub strand: Strand,
    /// The anchors in order.
    pub anchors: Vec<Anchor>,
    /// Query offset of every anchor from the first, counted across the origin of a circular query.
    pub offsets: Vec<u32>,
    /// The chain score: matched anchor bases less the cost of drifting links.
    pub score: i64,
}

/// How anchors are chained.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ChainConfig {
    /// Largest contig and query distance between consecutive anchors of a chain.
    pub max_gap: u32,
    /// The fewest anchors of a chain.
    pub min_anchors: u32,
    /// The shortest exact stretch a chain must reach. Its score may fall 2 (stride - 1) below,
    /// as the anchors of such a stretch can miss bases at both ends.
    pub min_score: i64,
}

impl Default for ChainConfig {
    fn default() -> Self {
        Self {
            max_gap: 5_000,
            min_anchors: 3,
            min_score: 40,
        }
    }
}

impl ChainConfig {
    /// Checks that the gap and the anchor count are positive.
    pub fn validate(self) -> Result<(), AlignmentError> {
        if self.max_gap == 0 || self.min_anchors == 0 {
            return Err(AlignmentError::InvalidChainConfig);
        }
        Ok(())
    }

    /// The shortest stretch a chain of anchors of length `a` can span; walks must reach it too.
    pub fn min_length(self, a: u8) -> u64 {
        let a = u64::from(a);
        a.max(self.min_score.max(1) as u64)
            .max(a + u64::from(self.min_anchors) - 1)
    }
}

/// Bytes of anchors and chaining arrays of one window.
pub const CHAINING_BYTES: usize = 64 << 20;

/// The chains of an alignment window, and the k-mers it leaves to dense walks.
pub struct WindowChains {
    /// The chains, best score first per strand.
    pub chains: Vec<Chain>,
    /// K-mers with more query occurrences than this gave no anchors, so dense walks recover
    /// their copies (0: every k-mer). `None` when every k-mer gave anchors.
    pub dense_above: Option<usize>,
}

/// The chains of a window within the chaining budget. The chains of a periodic window may hold
/// at most `periodic_bytes`, and chaining may examine at most `max_examined` predecessors.
/// `stopped` is asked now and then while chaining; when it says yes, chaining ends with
/// `AlignmentError::Stopped`.
pub fn window_chains(
    table: &QueryTable,
    window: &[u8],
    scratch: &mut Vec<Anchor>,
    periodic_bytes: usize,
    max_examined: u64,
    stopped: &dyn Fn() -> bool,
) -> Result<WindowChains, AlignmentError> {
    scratch.clear();
    // Every copy of a circular repeat is on many diagonals at once; dense walks take them all.
    if table.circular() && table.periodic_pair(window).is_some() {
        let chains = Vec::new();
        return Ok(WindowChains {
            chains,
            dense_above: Some(0),
        });
    }
    if let Some(chains) = periodic_chains(table, window, periodic_bytes)? {
        return Ok(WindowChains {
            chains,
            dense_above: None,
        });
    }
    // The cut depends on this window alone; a larger buffer kept from an earlier window is
    // shrunk, so it neither changes the cut nor holds memory beyond the budget.
    let budget = anchor_budget(table, window);
    scratch.shrink_to(budget.expected);
    find_anchors(table, window, scratch, budget.above, budget.expected);
    let mut chains = Vec::new();
    let split = scratch.partition_point(|anchor| anchor.strand == Strand::Forward);
    let mut examined = Examined::new(max_examined, stopped);
    examined.window_anchors = scratch.len();
    for strand in [&scratch[..split], &scratch[split..]] {
        chain_anchors(table, strand, &mut chains, &mut examined)?;
    }
    Ok(WindowChains {
        chains,
        dense_above: budget.above,
    })
}

/// Forward window positions `start..end` covered by the anchors of a chain.
pub fn forward_extent(chain: &Chain, window: u64, a: u64) -> (u64, u64) {
    let first = u64::from(chain.anchors[0].contig);
    let last = u64::from(chain.anchors[chain.anchors.len() - 1].contig) + a;
    if chain.strand.is_reverse() {
        (window - last, window - first)
    } else {
        (first, last)
    }
}

/// Splits a chain where its anchors drift more than `max_drift` off the first diagonal; parts of
/// fewer than `min_anchors` anchors are dropped.
pub fn split_by_drift(chain: &Chain, max_drift: u32, min_anchors: u32) -> Vec<Chain> {
    let mut parts: Vec<Chain> = Vec::new();
    let mut start = 0;
    for index in 0..=chain.anchors.len() {
        let drift = |index: usize| {
            let contig =
                i64::from(chain.anchors[index].contig) - i64::from(chain.anchors[start].contig);
            let query = i64::from(chain.offsets[index]) - i64::from(chain.offsets[start]);
            (contig - query).unsigned_abs()
        };
        if index < chain.anchors.len() && drift(index) <= u64::from(max_drift) {
            continue;
        }
        if index - start >= min_anchors as usize {
            let base = chain.offsets[start];
            parts.push(Chain {
                strand: chain.strand,
                anchors: chain.anchors[start..index].to_vec(),
                offsets: chain.offsets[start..index]
                    .iter()
                    .map(|offset| offset - base)
                    .collect(),
                score: chain.score,
            });
        }
        start = index;
    }
    parts
}

/// The anchors a window may give within `CHAINING_BYTES` and the query occurrence count above
/// which k-mers give none (`None`: all fit). When all anchors would not fit, only k-mers with the
/// fewest occurrences give anchors, those that fit half of the budget.
fn anchor_budget(table: &QueryTable, window: &[u8]) -> AnchorBudget {
    let price = |count: usize| anchor_bytes(count, table);
    // The most anchors a window can give, before its k-mers are looked up.
    let stride = usize::from(table.anchors().stride);
    let bound = (window.len() / stride + 1).saturating_mul(table.max_occurrences());
    if price(bound) <= CHAINING_BYTES {
        // Even a window of the most frequent k-mers fits, so the lookups need not be counted.
        return AnchorBudget {
            expected: bound,
            above: None,
        };
    }
    let count = count_anchors(table, window);
    if price(count) <= CHAINING_BYTES {
        return AnchorBudget {
            expected: count,
            above: None,
        };
    }
    let (mut sparse, mut threshold) = (0usize, 0usize);
    for (occurrences, count) in occurrence_counts(table, window) {
        let next = sparse.saturating_add(occurrences.saturating_mul(count));
        if price(next) > CHAINING_BYTES / 2 {
            break;
        }
        (sparse, threshold) = (next, occurrences);
    }
    AnchorBudget {
        expected: sparse,
        above: Some(threshold),
    }
}

/// The anchors a window gives at most, and the query occurrence count above which its k-mers
/// give none (`None`: all fit).
struct AnchorBudget {
    expected: usize,
    above: Option<usize>,
}

/// Bytes of `count` anchors and of the chaining arrays, traceback and chain headers that can
/// coexist with them.
fn anchor_bytes(count: usize, table: &QueryTable) -> usize {
    use std::mem::size_of;
    count
        .saturating_mul(
            2 * size_of::<Anchor>()
                + 2 * size_of::<i64>()
                + 4 * size_of::<usize>()
                + size_of::<bool>()
                + size_of::<u32>(),
        )
        .saturating_add(
            (2 * count / table.chains().min_anchors as usize + 4)
                .saturating_mul(size_of::<Chain>()),
        )
}

/// The chains of a linear periodic query in a window that is a copy of its repeat: one chain on
/// every diagonal of an exact copy, without storing the dense anchors, within `budget` bytes.
fn periodic_chains(
    table: &QueryTable,
    window: &[u8],
    budget: usize,
) -> Result<Option<Vec<Chain>>, AlignmentError> {
    if table.circular() {
        return Ok(None);
    }
    let Some(periodic) = table.periodic_pair(window) else {
        return Ok(None);
    };
    let (period, phases) = (periodic.period, periodic.phases);
    let a = u64::from(table.anchors().a);
    let anchors = table.chains().min_anchors as usize;
    let minimum = table.chains().min_length(table.anchors().a);
    let low = minimum as i64 - i64::from(table.query_length());
    let high = window.len() as i64 - minimum as i64;
    let mut total = 0usize;
    for phase in &phases {
        let first = low + (phase.diagonal - low).rem_euclid(period as i64);
        if first <= high {
            total = total.saturating_add(((high - first) / period as i64 + 1) as usize);
        }
    }
    let bytes = total.saturating_mul(
        2 * std::mem::size_of::<Chain>()
            + anchors.saturating_mul(std::mem::size_of::<Anchor>() + std::mem::size_of::<u32>()),
    );
    if bytes > budget {
        let needed = Some(bytes as u64);
        return Err(AlignmentError::ChainingBudget(
            Budgeted::PeriodicChains,
            needed,
        ));
    }
    let mut chains = Vec::with_capacity(total);
    for phase in phases {
        let mut diagonal = low + (phase.diagonal - low).rem_euclid(period as i64);
        while diagonal <= high {
            let (query_start, contig_start) = ((-diagonal).max(0) as u64, diagonal.max(0) as u64);
            let length = (u64::from(table.query_length()) - query_start)
                .min(window.len() as u64 - contig_start);
            let mut positions = Vec::with_capacity(anchors);
            let mut offsets = Vec::with_capacity(anchors);
            for at in 0..anchors {
                let step = (length - a) * at as u64 / anchors.saturating_sub(1).max(1) as u64;
                positions.push(Anchor {
                    strand: phase.strand,
                    contig: (contig_start + step) as u32,
                    query: (query_start + step) as u32,
                });
                offsets.push(step as u32);
            }
            chains.push(Chain {
                strand: phase.strand,
                anchors: positions,
                offsets,
                score: length as i64,
            });
            diagonal += period as i64;
        }
    }
    chains.sort_unstable_by_key(|chain| std::cmp::Reverse(chain.score));
    Ok(Some(chains))
}

/// Chains the anchors of one strand (sorted by contig, then query position), best score first.
/// Every anchor joins at most one chain. A chain of a circular query may cross the origin but
/// never covers more than the query length.
pub(super) fn chain_anchors(
    table: &QueryTable,
    anchors: &[Anchor],
    out: &mut Vec<Chain>,
    examined: &mut Examined<'_>,
) -> Result<(), AlignmentError> {
    let config = table.chains();
    // Lookups every `stride` positions can miss up to stride - 1 bases at each end of an exact
    // stretch, so a stretch of `min_score` bases may score 2 (stride - 1) less.
    let mut lowest = config.min_score - 2 * (i64::from(table.anchors().stride) - 1);
    // A chain of a short query may be one anchor, which scores a.
    let fewest = table.fewest_anchors();
    if fewest < config.min_anchors {
        lowest = lowest.min(i64::from(table.anchors().a));
    }
    let links = chain_links(table, anchors, examined)?;
    let mut order: Vec<usize> = (0..anchors.len()).collect();
    order.sort_unstable_by_key(|&index| (std::cmp::Reverse(links[index].score), index));
    let mut used = vec![false; anchors.len()];
    for end in order {
        if used[end] || links[end].score < lowest {
            continue;
        }
        let (mut path, mut base, mut at) = (Vec::new(), 0, end);
        loop {
            if used[at] {
                base = links[at].score;
                break;
            }
            used[at] = true;
            path.push(at);
            if links[at].previous == usize::MAX {
                break;
            }
            at = links[at].previous;
        }
        if path.len() < fewest as usize || links[end].score - base < lowest {
            continue;
        }
        path.reverse();
        let origin = links[path[0]].span;
        out.push(Chain {
            strand: anchors[end].strand,
            anchors: path.iter().map(|&index| anchors[index]).collect(),
            offsets: path
                .iter()
                .map(|&index| (links[index].span - origin) as u32)
                .collect(),
            score: links[end].score - base,
        });
    }
    Ok(())
}

/// Predecessors with a valid link that one anchor considers at most.
const LOOKBACK: usize = 64;

/// Diagonals of one group of the diagonal index.
const DIAGONALS: i64 = 1 << 10;

/// The anchors of one strand grouped by diagonal (contig less query position), each group in
/// anchor order: 4 bytes per anchor and per group.
struct DiagonalIndex {
    lowest: i64,
    starts: Vec<u32>,
    order: Vec<u32>,
}

impl DiagonalIndex {
    /// The index of `anchors`, whose window gives `window_anchors` anchors, for links that drift
    /// at most `drift`. `None` when one lookup could need more groups than `Near` holds, or when
    /// the index and what `anchor_bytes` prices for the window pass `CHAINING_BYTES`.
    fn new(
        table: &QueryTable,
        anchors: &[Anchor],
        window_anchors: usize,
        drift: i64,
    ) -> Option<Self> {
        let centres = if table.circular() { 2 } else { 1 };
        if centres * ((2 * drift) / DIAGONALS + 2) > NEAR_GROUPS as i64 {
            return None;
        }
        let diagonal = |anchor: &Anchor| i64::from(anchor.contig) - i64::from(anchor.query);
        let lowest = anchors.iter().map(diagonal).min().unwrap_or(0);
        let highest = anchors.iter().map(diagonal).max().unwrap_or(0);
        let groups = ((highest - lowest) / DIAGONALS) as usize + 1;
        let bytes = (anchors.len() + groups + 2) * std::mem::size_of::<u32>();
        let priced = anchor_bytes(window_anchors.max(anchors.len()), table);
        if priced.saturating_add(bytes) > CHAINING_BYTES {
            return None;
        }
        let group = |anchor: &Anchor| ((diagonal(anchor) - lowest) / DIAGONALS) as usize;
        // Counts at `group + 2`, summed: `starts[group + 1]` is where the group starts, and it
        // moves to the group's end as its anchors are placed.
        let mut starts = vec![0u32; groups + 2];
        for anchor in anchors {
            starts[group(anchor) + 2] += 1;
        }
        for at in 1..starts.len() {
            starts[at] += starts[at - 1];
        }
        let mut order = vec![0u32; anchors.len()];
        for (index, anchor) in anchors.iter().enumerate() {
            let at = &mut starts[group(anchor) + 1];
            order[*at as usize] = index as u32;
            *at += 1;
        }
        starts.pop();
        Some(Self {
            lowest,
            starts,
            order,
        })
    }

    /// The anchors before index `end` whose diagonal is within `drift` of the diagonal of
    /// `anchor` or, on a circular query, of that diagonal less the query length, from the last
    /// to the first. Their groups may hold other anchors too.
    fn near(&self, table: &QueryTable, anchor: &Anchor, end: usize, drift: i64) -> Near<'_> {
        let diagonal = i64::from(anchor.contig) - i64::from(anchor.query);
        let length = i64::from(table.query_length());
        let last = self.starts.len() as i64 - 2;
        let mut near = Near {
            heads: [&[]; NEAR_GROUPS],
            count: 0,
        };
        let mut previous = -1;
        let centres = [diagonal - length, diagonal];
        for centre in &centres[usize::from(!table.circular())..] {
            let first = (centre - drift - self.lowest)
                .div_euclid(DIAGONALS)
                .max(previous + 1);
            let high = (centre + drift - self.lowest)
                .div_euclid(DIAGONALS)
                .min(last);
            for group in first.max(0)..=high {
                let group = group as usize;
                let members =
                    &self.order[self.starts[group] as usize..self.starts[group + 1] as usize];
                let cut = members.partition_point(|&index| (index as usize) < end);
                near.heads[near.count] = &members[..cut];
                near.count += 1;
                previous = group as i64;
            }
        }
        near
    }
}

/// The predecessors of an anchor to examine, from the last to the first.
enum Predecessors<'a> {
    All(std::iter::Rev<std::ops::Range<usize>>),
    Near(Near<'a>),
}

impl Iterator for Predecessors<'_> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        match self {
            Self::All(all) => all.next(),
            Self::Near(near) => near.next(),
        }
    }
}

/// Whether a strand of `anchors` anchors that examined `total` predecessors so far looks up
/// the next ones by diagonal: once it passed the lookback of every anchor.
fn index_from(total: u64, anchors: usize) -> bool {
    total > (LOOKBACK * anchors) as u64
}

/// Groups of the diagonal index that one lookup reads at most.
const NEAR_GROUPS: usize = 8;

/// Anchor indices from up to `NEAR_GROUPS` groups of the diagonal index, highest first.
struct Near<'a> {
    heads: [&'a [u32]; NEAR_GROUPS],
    count: usize,
}

impl Iterator for Near<'_> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        let heads = &mut self.heads[..self.count];
        let (best, &index) = heads
            .iter()
            .enumerate()
            .filter_map(|(at, head)| head.last().map(|index| (at, index)))
            .max_by_key(|&(_, index)| index)?;
        heads[best] = &heads[best][..heads[best].len() - 1];
        Some(index as usize)
    }
}

/// Predecessors that the chaining of one window may examine by default, valid or not: seconds
/// of work. A window past it fails with `ChainingBudget` rather than chaining fewer anchors.
pub const MAX_EXAMINED: u64 = 1 << 32;

/// Predecessors examined between two questions whether to stop.
const STOP_CHECK: u64 = 1 << 20;

/// The predecessors examined so far in one window, their limit and the stop question.
pub(super) struct Examined<'a> {
    count: u64,
    limit: u64,
    stopped: &'a dyn Fn() -> bool,
    /// The anchors of the window, priced by `anchor_bytes` (0: those of one strand).
    window_anchors: usize,
}

impl<'a> Examined<'a> {
    pub(super) fn new(limit: u64, stopped: &'a dyn Fn() -> bool) -> Self {
        Self {
            count: 0,
            limit,
            stopped,
            window_anchors: 0,
        }
    }

    /// Adds `more` examined predecessors; fails past the limit or when asked to stop.
    fn add(&mut self, more: u64) -> Result<(), AlignmentError> {
        let before = self.count;
        self.count = self.count.saturating_add(more);
        if self.count > self.limit {
            return Err(AlignmentError::ChainingBudget(Budgeted::ChainingWork, None));
        }
        if before / STOP_CHECK != self.count / STOP_CHECK && (self.stopped)() {
            return Err(AlignmentError::Stopped);
        }
        Ok(())
    }
}

/// A link may leave the diagonal by this many bases plus one part in `DRIFT_SHARE` of its length.
const DRIFT_BASES: i64 = 50;
const DRIFT_SHARE: i64 = 5;

/// A link that leaves the diagonal costs one per this many drifted bases, plus one, plus half
/// the bits of the drift.
const DRIFT_COST: i64 = 10;

/// The best chain that ends in an anchor: its score, the anchor before (`usize::MAX` for none)
/// and the query span from the chain start.
#[derive(Clone, Copy)]
struct Link {
    score: i64,
    previous: usize,
    span: i64,
}

/// The best link of every anchor.
fn chain_links(
    table: &QueryTable,
    anchors: &[Anchor],
    examined: &mut Examined<'_>,
) -> Result<Vec<Link>, AlignmentError> {
    let a = i64::from(table.anchors().a);
    let max_gap = i64::from(table.chains().max_gap);
    let length = i64::from(table.query_length());
    // A link drifts at most this far from the diagonal, as its query step is at most the gap.
    let band = DRIFT_BASES + max_gap / DRIFT_SHARE;
    let (mut index, mut tried, mut total) = (None, false, 0u64);
    let mut links: Vec<Link> = Vec::with_capacity(anchors.len());
    for (current, anchor) in anchors.iter().enumerate() {
        let (mut best, mut from, mut best_span) = (a, usize::MAX, 0);
        let end = anchors[..current].partition_point(|earlier| earlier.contig < anchor.contig);
        let (mut eligible, mut scanned) = (0, 0u64);
        // Once the predecessors examined pass the lookback of every anchor, as most are far
        // off the diagonal, only those near it are examined, in the same order.
        if !tried && index_from(total, anchors.len()) {
            tried = true;
            index = DiagonalIndex::new(table, anchors, examined.window_anchors, band);
        }
        let predecessors = match &index {
            Some(index) => Predecessors::Near(index.near(table, anchor, end, band)),
            None => Predecessors::All((0..end).rev()),
        };
        for earlier in predecessors {
            scanned += 1;
            let contig_step = i64::from(anchor.contig) - i64::from(anchors[earlier].contig);
            if contig_step > max_gap {
                break;
            }
            let mut query_step = i64::from(anchor.query) - i64::from(anchors[earlier].query);
            if table.circular() && query_step <= 0 {
                query_step += length;
            }
            if contig_step <= 0
                || query_step <= 0
                || query_step > max_gap
                || links[earlier].span + query_step + a > length
            {
                continue;
            }
            let drift = (contig_step - query_step).abs();
            if drift > DRIFT_BASES + contig_step.min(query_step) / DRIFT_SHARE {
                continue;
            }
            eligible += 1;
            let gap_cost = if drift == 0 {
                0
            } else {
                drift / DRIFT_COST + 1 + i64::from(64 - (drift as u64).leading_zeros()) / 2
            };
            let gained = links[earlier].score + contig_step.min(query_step).min(a) - gap_cost;
            if gained > best {
                (best, from, best_span) = (gained, earlier, links[earlier].span + query_step);
            }
            if eligible == LOOKBACK {
                break;
            }
        }
        total += scanned;
        examined.add(scanned)?;
        links.push(Link {
            score: best,
            previous: from,
            span: best_span,
        });
    }
    Ok(links)
}

#[cfg(test)]
#[path = "chains_tests.rs"]
mod tests;
