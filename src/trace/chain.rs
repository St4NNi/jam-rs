//! Bounded, deterministic chaining of positional trace anchors.
//!
//! Chaining is evidence selection only.  It groups anchors by contig, strand,
//! and k-mer level, scores monotone paths with a bounded predecessor window,
//! and exposes circular query segments for a path that crosses the plasmid
//! origin.  It does not perform alignment or make a presence call.

use crate::trace::anchors::Anchor;
use crate::trace::config::SensitivityConfig;
use crate::trace::model::{BaseInterval, CoordinateError, Strand};
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Limits and scoring constants for one bounded chain search.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct ChainConfig {
    pub max_chains: u32,
    pub min_anchors: u32,
    pub max_predecessors: u32,
    pub max_query_gap: u64,
    pub max_target_gap: u64,
    pub gap_penalty: i64,
}

impl ChainConfig {
    #[must_use]
    pub fn from_sensitivity(config: &SensitivityConfig) -> Self {
        Self {
            max_chains: config.max_chains_per_candidate,
            min_anchors: config.min_chain_anchors,
            max_predecessors: config.max_anchors_per_candidate.min(256),
            max_query_gap: config.max_alignment_window_bases,
            max_target_gap: config.max_alignment_window_bases,
            gap_penalty: 1,
        }
    }
}

impl Default for ChainConfig {
    fn default() -> Self {
        Self {
            max_chains: 8,
            min_anchors: 3,
            max_predecessors: 256,
            max_query_gap: 1_000_000,
            max_target_gap: 1_000_000,
            gap_penalty: 1,
        }
    }
}

/// One retained monotone anchor path.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct AnchorChain {
    pub contig_id: u32,
    pub strand: Strand,
    pub k: u8,
    pub anchors: Vec<Anchor>,
    pub score: i64,
    /// Normalized first query segment.  For origin-crossing chains,
    /// `query_segments` is authoritative and contains the second segment too.
    pub query_interval: BaseInterval,
    pub query_segments: Vec<BaseInterval>,
    pub target_interval: BaseInterval,
    pub origin_crossing: bool,
    /// Unwrapped coordinates used by the chain DP; useful for diagnostics.
    pub linear_query_start: u64,
    pub linear_query_end: u64,
}

impl AnchorChain {
    #[must_use]
    pub fn anchor_count(&self) -> usize {
        self.anchors.len()
    }
}

#[derive(Debug, Error)]
pub enum ChainError {
    #[error("chain query length must be greater than zero")]
    ZeroQueryLength,
    #[error("chain limit {field} must be greater than zero")]
    ZeroLimit { field: &'static str },
    #[error("chain gap penalty must not be negative")]
    NegativeGapPenalty,
    #[error("anchor k must be greater than zero")]
    ZeroK,
    #[error(
        "anchor query interval is outside the circular query: position={position}, k={k}, length={length}"
    )]
    AnchorOutsideQuery { position: u64, k: u8, length: u64 },
    #[error("chain coordinate overflow")]
    CoordinateOverflow,
    #[error("invalid chain interval: {0}")]
    InvalidInterval(#[from] CoordinateError),
}

#[derive(Clone, Copy)]
struct ExpandedAnchor {
    anchor: Anchor,
    query_position: u64,
    source_index: usize,
}

/// Build up to `config.max_chains` chains per contig/strand/k group.
///
/// Query positions are duplicated one plasmid length apart, allowing one
/// origin crossing while rejecting paths that use the same source anchor
/// twice.  The predecessor scan is capped by `max_predecessors`, so runtime
/// and intermediate state remain bounded by the caller's anchor limit.
pub fn chain_anchors(
    anchors: &[Anchor],
    query_length: u64,
    config: ChainConfig,
) -> Result<Vec<AnchorChain>, ChainError> {
    validate_config(config, query_length)?;

    let mut validated = Vec::with_capacity(anchors.len());
    for (source_index, &anchor) in anchors.iter().enumerate() {
        if anchor.k == 0 {
            return Err(ChainError::ZeroK);
        }
        let end = anchor
            .query_position
            .checked_add(u64::from(anchor.k))
            .ok_or(ChainError::CoordinateOverflow)?;
        if end > query_length {
            return Err(ChainError::AnchorOutsideQuery {
                position: anchor.query_position,
                k: anchor.k,
                length: query_length,
            });
        }
        anchor
            .target_position
            .checked_add(u64::from(anchor.k))
            .ok_or(ChainError::CoordinateOverflow)?;
        validated.push((source_index, anchor));
    }

    // The grouping key is explicit, keeping primary and rescue coordinates
    // from being mixed when both seed levels are present in one candidate.
    validated.sort_by_key(|(_, anchor)| anchor_sort_key(anchor));
    let max_chains = usize::try_from(config.max_chains).unwrap_or(usize::MAX);
    // Keep the cross-group result bounded.  Retaining only the best global
    // `max_chains` candidates after each insertion is sufficient because the
    // final ordering is total and deterministic.
    let mut all_chains = Vec::with_capacity(validated.len().min(max_chains));
    let mut group_start = 0usize;
    while group_start < validated.len() {
        let current_group_key = group_key(&validated[group_start].1);
        let mut group_end = group_start + 1;
        while group_end < validated.len() && group_key(&validated[group_end].1) == current_group_key
        {
            group_end += 1;
        }
        let group = &validated[group_start..group_end];
        let mut used = vec![false; group.len()];
        for _ in 0..config.max_chains {
            let Some(chain) = best_chain(group, &used, query_length, config)? else {
                break;
            };
            for source_index in &chain.source_indices {
                used[*source_index] = true;
            }
            all_chains.push(chain.chain);
            all_chains.sort_by(chain_sort_key);
            all_chains.truncate(max_chains);
        }
        group_start = group_end;
    }

    all_chains.sort_by(chain_sort_key);
    all_chains.truncate(max_chains);
    Ok(all_chains)
}

/// Alias for callers that use the shorter operation name.
pub fn build_chains(
    anchors: &[Anchor],
    query_length: u64,
    config: ChainConfig,
) -> Result<Vec<AnchorChain>, ChainError> {
    chain_anchors(anchors, query_length, config)
}

struct BestChain {
    chain: AnchorChain,
    source_indices: Vec<usize>,
}

fn best_chain(
    group: &[(usize, Anchor)],
    used: &[bool],
    query_length: u64,
    config: ChainConfig,
) -> Result<Option<BestChain>, ChainError> {
    if group.is_empty() {
        return Ok(None);
    }

    let expanded_capacity = group
        .len()
        .checked_mul(2)
        .ok_or(ChainError::CoordinateOverflow)?;
    let mut expanded = Vec::with_capacity(expanded_capacity);
    for (group_index, &(source_index, anchor)) in group.iter().enumerate() {
        if used[group_index] {
            continue;
        }
        expanded.push(ExpandedAnchor {
            anchor,
            query_position: anchor.query_position,
            source_index: group_index,
        });
        let wrapped = anchor
            .query_position
            .checked_add(query_length)
            .ok_or(ChainError::CoordinateOverflow)?;
        // Keep `source_index` tied to the group position, not the caller's
        // original slice index, so duplicate circular copies can be blocked.
        let _ = source_index;
        expanded.push(ExpandedAnchor {
            anchor,
            query_position: wrapped,
            source_index: group_index,
        });
    }
    expanded.sort_by(|left, right| expanded_cmp(*left, *right));

    let length = expanded.len();
    let mut scores = vec![i64::MIN; length];
    let mut counts = vec![0u32; length];
    let mut predecessors: Vec<Option<usize>> = vec![None; length];

    for i in 0..length {
        let current = expanded[i];
        scores[i] = anchor_base_score(current.anchor);
        counts[i] = 1;
        let mut examined = 0u32;
        let mut j = i;
        while j > 0 && examined < config.max_predecessors {
            j -= 1;
            // Count every candidate examined, including same-source and
            // non-monotone candidates.  Otherwise a repeat-heavy group could
            // scan the whole preceding vector despite the predecessor cap.
            examined += 1;
            let previous = expanded[j];
            if current.query_position <= previous.query_position {
                continue;
            }
            let query_gap = current.query_position - previous.query_position;
            if query_gap > config.max_query_gap {
                break;
            }
            if current.source_index == previous.source_index {
                continue;
            }
            let target_gap = match current.anchor.strand {
                Strand::Forward => {
                    if current.anchor.target_position <= previous.anchor.target_position {
                        continue;
                    }
                    current.anchor.target_position - previous.anchor.target_position
                }
                Strand::Reverse => {
                    if current.anchor.target_position >= previous.anchor.target_position {
                        continue;
                    }
                    previous.anchor.target_position - current.anchor.target_position
                }
            };
            if target_gap == 0 || target_gap > config.max_target_gap {
                continue;
            }
            if scores[j] == i64::MIN {
                continue;
            }
            let gap = query_gap.abs_diff(target_gap);
            let extension = gap_extension(gap, config.gap_penalty);
            let candidate_score = scores[j].saturating_add(extension);
            let candidate_count = counts[j].saturating_add(1);
            if better_predecessor(
                candidate_score,
                candidate_count,
                j,
                scores[i],
                counts[i],
                predecessors[i],
                &expanded,
            ) {
                scores[i] = candidate_score;
                counts[i] = candidate_count;
                predecessors[i] = Some(j);
            }
        }
    }

    let mut endpoint: Option<usize> = None;
    for index in 0..length {
        if counts[index] < config.min_anchors {
            continue;
        }
        let replace = endpoint
            .is_none_or(|current| better_endpoint(index, current, &scores, &counts, &expanded));
        if replace {
            endpoint = Some(index);
        }
    }
    let Some(endpoint) = endpoint else {
        return Ok(None);
    };

    let mut path = Vec::with_capacity(counts[endpoint] as usize);
    let mut cursor = Some(endpoint);
    while let Some(index) = cursor {
        path.push(index);
        cursor = predecessors[index];
    }
    path.reverse();

    let first = expanded[path[0]];
    let last = expanded[*path.last().expect("non-empty chain")];
    let linear_start = first.query_position;
    let linear_end = last
        .query_position
        .checked_add(u64::from(last.anchor.k))
        .ok_or(ChainError::CoordinateOverflow)?;
    let span = linear_end
        .checked_sub(linear_start)
        .ok_or(ChainError::CoordinateOverflow)?;
    if span > query_length {
        return Ok(None);
    }

    let query_segments = circular_segments(linear_start, span, query_length)?;
    let query_interval = *query_segments
        .first()
        .expect("non-empty chain has a query segment");
    let origin_crossing = query_segments.len() == 2;

    let (target_start, target_end) =
        path.iter()
            .try_fold((u64::MAX, 0u64), |(start, end), &index| {
                let anchor = expanded[index].anchor;
                let anchor_end = anchor
                    .target_position
                    .checked_add(u64::from(anchor.k))
                    .ok_or(ChainError::CoordinateOverflow)?;
                Ok::<_, ChainError>((start.min(anchor.target_position), end.max(anchor_end)))
            })?;
    let target_interval = BaseInterval::new(target_start, target_end)?;
    let anchors = path.iter().map(|&index| expanded[index].anchor).collect();
    let source_indices = path
        .iter()
        .map(|&index| expanded[index].source_index)
        .collect();
    Ok(Some(BestChain {
        chain: AnchorChain {
            contig_id: first.anchor.contig_id,
            strand: first.anchor.strand,
            k: first.anchor.k,
            anchors,
            score: scores[endpoint],
            query_interval,
            query_segments,
            target_interval,
            origin_crossing,
            linear_query_start: linear_start,
            linear_query_end: linear_end,
        },
        source_indices,
    }))
}

fn validate_config(config: ChainConfig, query_length: u64) -> Result<(), ChainError> {
    if query_length == 0 {
        return Err(ChainError::ZeroQueryLength);
    }
    if config.max_chains == 0 {
        return Err(ChainError::ZeroLimit {
            field: "max_chains",
        });
    }
    if config.min_anchors == 0 {
        return Err(ChainError::ZeroLimit {
            field: "min_anchors",
        });
    }
    if config.max_predecessors == 0 {
        return Err(ChainError::ZeroLimit {
            field: "max_predecessors",
        });
    }
    if config.max_query_gap == 0 {
        return Err(ChainError::ZeroLimit {
            field: "max_query_gap",
        });
    }
    if config.max_target_gap == 0 {
        return Err(ChainError::ZeroLimit {
            field: "max_target_gap",
        });
    }
    if config.gap_penalty < 0 {
        return Err(ChainError::NegativeGapPenalty);
    }
    Ok(())
}

fn anchor_base_score(_anchor: Anchor) -> i64 {
    100
}

fn gap_extension(gap: u64, gap_penalty: i64) -> i64 {
    let gap = i64::try_from(gap).unwrap_or(i64::MAX);
    100i64.saturating_sub(gap.saturating_mul(gap_penalty))
}

fn better_predecessor(
    candidate_score: i64,
    candidate_count: u32,
    candidate_index: usize,
    current_score: i64,
    current_count: u32,
    current_predecessor: Option<usize>,
    expanded: &[ExpandedAnchor],
) -> bool {
    if candidate_score != current_score {
        return candidate_score > current_score;
    }
    if candidate_count != current_count {
        return candidate_count > current_count;
    }
    match current_predecessor {
        None => true,
        Some(current_index) if candidate_index != current_index => {
            expanded[candidate_index].query_position < expanded[current_index].query_position
        }
        Some(_) => false,
    }
}

fn better_endpoint(
    candidate: usize,
    current: usize,
    scores: &[i64],
    counts: &[u32],
    expanded: &[ExpandedAnchor],
) -> bool {
    if scores[candidate] != scores[current] {
        return scores[candidate] > scores[current];
    }
    if counts[candidate] != counts[current] {
        return counts[candidate] > counts[current];
    }
    let candidate_anchor = expanded[candidate].anchor;
    let current_anchor = expanded[current].anchor;
    (
        expanded[candidate].query_position,
        candidate_anchor.target_position,
        candidate_anchor.hash,
    ) < (
        expanded[current].query_position,
        current_anchor.target_position,
        current_anchor.hash,
    )
}

fn circular_segments(
    linear_start: u64,
    span: u64,
    query_length: u64,
) -> Result<Vec<BaseInterval>, ChainError> {
    if query_length == 0 || span == 0 || span > query_length {
        return Err(ChainError::CoordinateOverflow);
    }
    let start = linear_start % query_length;
    let end = start
        .checked_add(span)
        .ok_or(ChainError::CoordinateOverflow)?;
    if end <= query_length {
        Ok(vec![BaseInterval::new(start, end)?])
    } else {
        Ok(vec![
            BaseInterval::new(start, query_length)?,
            BaseInterval::new(0, end - query_length)?,
        ])
    }
}

fn group_key(anchor: &Anchor) -> (u32, u8, u8) {
    (anchor.contig_id, strand_rank(anchor.strand), anchor.k)
}

fn anchor_sort_key(anchor: &Anchor) -> (u32, u8, u8, u64, u64, u64, u64) {
    (
        anchor.contig_id,
        strand_rank(anchor.strand),
        anchor.k,
        anchor.query_position,
        anchor.target_position,
        anchor.hash,
        anchor.canonical_kmer,
    )
}

fn expanded_cmp(left: ExpandedAnchor, right: ExpandedAnchor) -> std::cmp::Ordering {
    left.query_position
        .cmp(&right.query_position)
        .then_with(|| match left.anchor.strand {
            Strand::Forward => left
                .anchor
                .target_position
                .cmp(&right.anchor.target_position),
            Strand::Reverse => right
                .anchor
                .target_position
                .cmp(&left.anchor.target_position),
        })
        .then_with(|| left.anchor.hash.cmp(&right.anchor.hash))
        .then_with(|| left.source_index.cmp(&right.source_index))
}

fn chain_sort_key(left: &AnchorChain, right: &AnchorChain) -> std::cmp::Ordering {
    right
        .score
        .cmp(&left.score)
        .then_with(|| right.anchors.len().cmp(&left.anchors.len()))
        .then_with(|| left.contig_id.cmp(&right.contig_id))
        .then_with(|| strand_rank(left.strand).cmp(&strand_rank(right.strand)))
        .then_with(|| left.k.cmp(&right.k))
        .then_with(|| left.linear_query_start.cmp(&right.linear_query_start))
        .then_with(|| left.target_interval.start.cmp(&right.target_interval.start))
}

const fn strand_rank(strand: Strand) -> u8 {
    match strand {
        Strand::Forward => 0,
        Strand::Reverse => 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn anchor(query: u64, target: u64, strand: Strand) -> Anchor {
        Anchor {
            query_position: query,
            target_position: target,
            contig_id: 0,
            strand,
            k: 3,
            hash: query + 1,
            canonical_kmer: query + 9,
            query_reverse: false,
            target_reverse: false,
        }
    }

    fn config() -> ChainConfig {
        ChainConfig {
            max_chains: 4,
            min_anchors: 2,
            max_predecessors: 32,
            max_query_gap: 100,
            max_target_gap: 100,
            gap_penalty: 1,
        }
    }

    #[test]
    fn forward_anchors_form_a_monotone_chain() {
        let anchors = [
            anchor(5, 105, Strand::Forward),
            anchor(15, 115, Strand::Forward),
            anchor(25, 125, Strand::Forward),
        ];
        let chains = chain_anchors(&anchors, 100, config()).unwrap();
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].anchors.len(), 3);
        assert_eq!(
            chains[0].target_interval,
            BaseInterval::new(105, 128).unwrap()
        );
        assert!(!chains[0].origin_crossing);
    }

    #[test]
    fn reverse_anchors_chain_in_decreasing_target_order() {
        let anchors = [
            anchor(5, 125, Strand::Reverse),
            anchor(15, 115, Strand::Reverse),
            anchor(25, 105, Strand::Reverse),
        ];
        let chains = chain_anchors(&anchors, 100, config()).unwrap();
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].anchors.len(), 3);
        assert_eq!(chains[0].strand, Strand::Reverse);
    }

    #[test]
    fn origin_crossing_emits_two_query_segments() {
        let anchors = [
            anchor(90, 190, Strand::Forward),
            anchor(2, 202, Strand::Forward),
            anchor(12, 212, Strand::Forward),
        ];
        let chains = chain_anchors(&anchors, 100, config()).unwrap();
        assert_eq!(chains.len(), 1);
        assert!(chains[0].origin_crossing);
        assert_eq!(chains[0].query_segments.len(), 2);
        assert_eq!(
            chains[0].query_segments[0],
            BaseInterval::new(90, 100).unwrap()
        );
        assert_eq!(
            chains[0].query_segments[1],
            BaseInterval::new(0, 15).unwrap()
        );
    }

    #[test]
    fn duplicate_circular_copy_cannot_reuse_an_anchor() {
        let anchors = [
            anchor(90, 190, Strand::Forward),
            anchor(2, 202, Strand::Forward),
        ];
        let mut chain_config = config();
        chain_config.min_anchors = 3;
        let chains = chain_anchors(&anchors, 100, chain_config).unwrap();
        assert!(chains.is_empty());
    }

    #[test]
    fn limits_and_invalid_coordinates_are_reported() {
        let mut invalid = config();
        invalid.max_predecessors = 0;
        assert!(matches!(
            chain_anchors(&[], 100, invalid),
            Err(ChainError::ZeroLimit {
                field: "max_predecessors"
            })
        ));
        let bad = anchor(99, 1, Strand::Forward);
        assert!(matches!(
            chain_anchors(&[bad], 100, config()),
            Err(ChainError::AnchorOutsideQuery { .. })
        ));
    }

    #[test]
    fn predecessor_limit_counts_rejected_candidates() {
        let anchors = [
            anchor(5, 100, Strand::Forward),
            anchor(15, 150, Strand::Forward),
            anchor(25, 150, Strand::Forward),
            anchor(35, 130, Strand::Forward),
        ];
        let mut limited = config();
        limited.max_predecessors = 2;
        limited.min_anchors = 3;
        assert!(chain_anchors(&anchors, 200, limited).unwrap().is_empty());
    }

    #[test]
    fn global_chain_limit_is_applied_across_contigs() {
        let anchors = [
            anchor(5, 100, Strand::Forward),
            anchor(15, 110, Strand::Forward),
            Anchor {
                contig_id: 1,
                ..anchor(5, 300, Strand::Forward)
            },
            Anchor {
                contig_id: 1,
                ..anchor(15, 310, Strand::Forward)
            },
        ];
        let mut limited = config();
        limited.max_chains = 1;
        let chains = chain_anchors(&anchors, 200, limited).unwrap();
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].contig_id, 0);
    }

    #[test]
    fn predecessor_limit_bounds_repeat_heavy_groups() {
        let anchors = (0..512)
            .map(|index| anchor(5 + index * 3, 100, Strand::Forward))
            .collect::<Vec<_>>();
        let mut repeat_config = config();
        repeat_config.max_predecessors = 4;
        repeat_config.max_query_gap = 10_000;
        assert!(
            chain_anchors(&anchors, 10_000, repeat_config)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn large_target_offsets_are_checked_without_wrapping() {
        let near_limit = u64::MAX - 20;
        let anchors = [
            anchor(5, near_limit, Strand::Forward),
            anchor(15, near_limit + 10, Strand::Forward),
        ];
        let chains = chain_anchors(&anchors, u64::MAX - 100, config()).unwrap();
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].target_interval.start, near_limit);
        assert_eq!(chains[0].target_interval.end, near_limit + 13);

        let overflowing = anchor(5, u64::MAX - 1, Strand::Forward);
        assert!(matches!(
            chain_anchors(&[overflowing], 100, config()),
            Err(ChainError::CoordinateOverflow)
        ));
    }

    #[test]
    fn rearranged_target_order_is_split_into_separate_chains() {
        let anchors = [
            anchor(5, 100, Strand::Forward),
            anchor(15, 115, Strand::Forward),
            anchor(25, 90, Strand::Forward),
            anchor(35, 105, Strand::Forward),
        ];
        let mut chain_config = config();
        chain_config.max_chains = 2;
        let chains = chain_anchors(&anchors, 200, chain_config).unwrap();
        assert_eq!(chains.len(), 2);
        assert!(chains.iter().all(|chain| chain.anchors.len() == 2));
        assert_eq!(chains[0].target_interval.start, 100);
        assert_eq!(chains[1].target_interval.start, 90);
    }
}
