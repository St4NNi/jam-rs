//! Chains: anchors of one strand in order along nearly the same diagonal.

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
    /// The lowest chain score.
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

/// Predecessors with a valid link that one anchor considers at most.
const LOOKBACK: usize = 64;

/// A link may leave the diagonal by this many bases plus one part in `DRIFT_SHARE` of its length.
const DRIFT_BASES: i64 = 50;
const DRIFT_SHARE: i64 = 5;

/// A link that leaves the diagonal costs one per this many drifted bases, plus one, plus half
/// the bits of the drift.
const DRIFT_COST: i64 = 10;

/// The chains of an alignment window, and the k-mers it leaves to dense walks.
pub struct WindowChains {
    /// The chains, best score first per strand.
    pub chains: Vec<Chain>,
    /// K-mers with more query occurrences than this gave no anchors, so dense walks recover
    /// their copies (0: every k-mer). `None` when every k-mer gave anchors.
    pub dense_above: Option<usize>,
}

/// The chains of a window within the chaining budget.
pub fn window_chains(
    table: &QueryTable,
    window: &[u8],
    scratch: &mut Vec<Anchor>,
) -> Result<WindowChains, AlignmentError> {
    scratch.clear();
    if table.circular() && table.periodic_pair(window).is_some() {
        let chains = Vec::new();
        return Ok(WindowChains {
            chains,
            dense_above: Some(0),
        });
    }
    if let Some(chains) = periodic_chains(table, window)? {
        return Ok(WindowChains {
            chains,
            dense_above: None,
        });
    }
    let budget = anchor_budget(table, window, scratch.capacity());
    find_anchors(table, window, scratch, budget.above, budget.expected);
    let mut chains = Vec::new();
    let split = scratch.partition_point(|anchor| anchor.strand == Strand::Forward);
    for strand in [&scratch[..split], &scratch[split..]] {
        chain_anchors(table, strand, &mut chains);
    }
    Ok(WindowChains {
        chains,
        dense_above: budget.above,
    })
}

/// The anchors a window may give within `CHAINING_BYTES` and the query occurrence count above
/// which k-mers give none (`None`: all fit). When all anchors would not fit, only k-mers with the
/// fewest occurrences give anchors, those that fit half of the budget.
fn anchor_budget(table: &QueryTable, window: &[u8], capacity: usize) -> AnchorBudget {
    let price = |count: usize| anchor_bytes(capacity, count, table);
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
fn anchor_bytes(capacity: usize, count: usize, table: &QueryTable) -> usize {
    use std::mem::size_of;
    capacity
        .max(count)
        .saturating_mul(size_of::<Anchor>())
        .saturating_add(count.saturating_mul(
            2 * size_of::<i64>()
                + 4 * size_of::<usize>()
                + size_of::<bool>()
                + size_of::<Anchor>()
                + size_of::<u32>(),
        ))
        .saturating_add(
            (2 * count / table.chains().min_anchors as usize + 4)
                .saturating_mul(size_of::<Chain>()),
        )
}

/// The chains of a linear periodic query in a window that is a copy of its repeat: one chain on
/// every diagonal of an exact copy, without storing the dense anchors.
fn periodic_chains(
    table: &QueryTable,
    window: &[u8],
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
    if bytes > CHAINING_BYTES {
        return Err(AlignmentError::ChainingBudget(Budgeted::PeriodicChains));
    }
    let mut chains = Vec::with_capacity(total);
    for phase in phases {
        let mut diagonal = low + (phase.diagonal - low).rem_euclid(period as i64);
        while diagonal <= high {
            let (q, t) = ((-diagonal).max(0) as u64, diagonal.max(0) as u64);
            let length = (u64::from(table.query_length()) - q).min(window.len() as u64 - t);
            let mut positions = Vec::with_capacity(anchors);
            let mut offsets = Vec::with_capacity(anchors);
            for at in 0..anchors {
                let step = (length - a) * at as u64 / anchors.saturating_sub(1).max(1) as u64;
                positions.push(Anchor {
                    strand: phase.strand,
                    contig: (t + step) as u32,
                    query: (q + step) as u32,
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
pub(super) fn chain_anchors(table: &QueryTable, anchors: &[Anchor], out: &mut Vec<Chain>) {
    let config = table.chains();
    let links = chain_links(table, anchors);
    let n = anchors.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by_key(|&index| (std::cmp::Reverse(links[index].score), index));
    let mut used = vec![false; n];
    for end in order {
        if used[end] || links[end].score < config.min_score {
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
        if path.len() < config.min_anchors as usize || links[end].score - base < config.min_score {
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
}

/// The best chain that ends in an anchor: its score, the anchor before (`usize::MAX` for none)
/// and the query span from the chain start.
#[derive(Clone, Copy)]
struct Link {
    score: i64,
    previous: usize,
    span: i64,
}

/// The best link of every anchor.
fn chain_links(table: &QueryTable, anchors: &[Anchor]) -> Vec<Link> {
    let a = i64::from(table.anchors().a);
    let max_gap = i64::from(table.chains().max_gap);
    let length = i64::from(table.query_length());
    let n = anchors.len();
    let mut links: Vec<Link> = Vec::with_capacity(n);
    for j in 0..n {
        let (mut best, mut from, mut best_span) = (a, usize::MAX, 0);
        let end = anchors[..j].partition_point(|anchor| anchor.contig < anchors[j].contig);
        let mut eligible = 0;
        for i in (0..end).rev() {
            let dt = i64::from(anchors[j].contig) - i64::from(anchors[i].contig);
            if dt > max_gap {
                break;
            }
            let mut dq = i64::from(anchors[j].query) - i64::from(anchors[i].query);
            if table.circular() && dq <= 0 {
                dq += length;
            }
            if dt <= 0 || dq <= 0 || dq > max_gap || links[i].span + dq + a > length {
                continue;
            }
            let drift = (dt - dq).abs();
            if drift > DRIFT_BASES + dt.min(dq) / DRIFT_SHARE {
                continue;
            }
            eligible += 1;
            let gap_cost = if drift == 0 {
                0
            } else {
                drift / DRIFT_COST + 1 + i64::from(64 - (drift as u64).leading_zeros()) / 2
            };
            let gained = links[i].score + dt.min(dq).min(a) - gap_cost;
            if gained > best {
                (best, from, best_span) = (gained, i, links[i].span + dq);
            }
            if eligible == LOOKBACK {
                break;
            }
        }
        links.push(Link {
            score: best,
            previous: from,
            span: best_span,
        });
    }
    links
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

#[cfg(test)]
#[path = "chains_tests.rs"]
pub(super) mod tests;
