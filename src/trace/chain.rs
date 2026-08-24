//! Bounded, deterministic chaining of positional trace anchors.
//!
//! Chaining is evidence selection only. It groups anchors by contig and
//! strand, scores mixed seed classes with a bounded predecessor window, and
//! optionally exposes wrapped query segments for a path that crosses the
//! stored query origin. It does not perform alignment or make a presence
//! call.

use crate::trace::anchors::Anchor;
use crate::trace::config::SensitivityConfig;
use crate::trace::model::{BaseInterval, CoordinateError, CoordinateModel, Strand};
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
    /// Query-coordinate handling for the chain DP. Only `Wrap` duplicates
    /// anchors one query length apart and permits an origin-crossing path.
    /// `Linear` and `Undetermined` retain the input coordinate line.
    pub coordinate_model: CoordinateModel,
}

/// Evidence class used by the mixed seed-chain scorer.
///
/// The class is deliberately carried beside the stable `Anchor` wire type so
/// callers that still use the legacy k=31/k=21 anchor constructor do not need
/// to rewrite their serialized records.  Experimental schemes can attach a
/// class through [`WeightedAnchor`] without changing the physical anchor
/// layout.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum AnchorClass {
    /// A verified adjacent run produced by exact Gear fragments.
    ExactGearRun,
    /// A low-occurrence canonical k=31 seed.
    SpecificK31,
    /// A low-occurrence spaced exact seed.
    SpecificSpaced,
    /// A low-occurrence paired or strobemer seed.
    SpecificPaired,
    /// A low-occurrence canonical k=21 rescue seed.
    SpecificK21,
    /// An occurrence retained only for bounded repeat-aware diagnostics.
    Repetitive,
}

impl AnchorClass {
    #[must_use]
    pub const fn score(self) -> i64 {
        match self {
            Self::ExactGearRun => 400,
            Self::SpecificK31 => 140,
            Self::SpecificSpaced => 120,
            Self::SpecificPaired => 115,
            Self::SpecificK21 => 80,
            Self::Repetitive => 10,
        }
    }

    #[must_use]
    pub const fn high_specificity(self) -> bool {
        matches!(self, Self::ExactGearRun | Self::SpecificK31)
    }

    #[must_use]
    pub const fn lower_specificity(self) -> bool {
        !matches!(self, Self::ExactGearRun | Self::SpecificK31)
    }
}

/// An anchor with explicit seed evidence class and optional occurrence count.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct WeightedAnchor {
    pub anchor: Anchor,
    pub class: AnchorClass,
    /// Verified span in query and target coordinates. Exact Gear runs may
    /// exceed the one-byte `Anchor.k` seed width.
    pub span: u64,
    /// Number of exact occurrences represented by the seed.  This is used
    /// only for diagnostics and never substitutes for exact key checking.
    pub occurrence_count: u32,
}

impl WeightedAnchor {
    #[must_use]
    pub const fn new(anchor: Anchor, class: AnchorClass) -> Self {
        Self {
            anchor,
            class,
            span: anchor.k as u64,
            occurrence_count: 1,
        }
    }

    /// Attach an exact verified run span, retaining at least the seed width.
    #[must_use]
    pub const fn with_span(mut self, span: u64) -> Self {
        self.span = if span < self.anchor.k as u64 {
            self.anchor.k as u64
        } else {
            span
        };
        self
    }

    #[must_use]
    pub const fn with_occurrence_count(mut self, occurrence_count: u32) -> Self {
        self.occurrence_count = occurrence_count;
        self
    }
}

/// Mixed chain output retaining the classes used by the bounded scorer.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct WeightedAnchorChain {
    pub contig_id: u32,
    pub strand: Strand,
    pub anchors: Vec<WeightedAnchor>,
    pub score: i64,
    pub query_interval: BaseInterval,
    pub query_segments: Vec<BaseInterval>,
    pub target_interval: BaseInterval,
    pub origin_crossing: bool,
    pub linear_query_start: u64,
    pub linear_query_end: u64,
}

impl WeightedAnchorChain {
    #[must_use]
    pub fn anchor_count(&self) -> usize {
        self.anchors.len()
    }

    #[must_use]
    pub fn has_high_specificity_anchor(&self) -> bool {
        self.anchors
            .iter()
            .any(|anchor| anchor.class.high_specificity())
    }
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
            // Preserve the pre-topology chain behavior for callers that do
            // not yet resolve a model. The runner selects Linear or Wrap
            // explicitly when it has a topology model to evaluate.
            coordinate_model: CoordinateModel::Wrap,
        }
    }

    /// Return this configuration with explicit query-coordinate handling.
    #[must_use]
    pub const fn with_coordinate_model(mut self, coordinate_model: CoordinateModel) -> Self {
        self.coordinate_model = coordinate_model;
        self
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
            coordinate_model: CoordinateModel::Wrap,
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
    #[error("anchor span {span} is shorter than k={k}")]
    InvalidAnchorSpan { span: u64, k: u8 },
    #[error(
        "anchor query interval is outside the query: position={position}, k={k}, length={length}"
    )]
    AnchorOutsideQuery { position: u64, k: u8, length: u64 },
    #[error("chain coordinate overflow")]
    CoordinateOverflow,
    #[error("invalid chain interval: {0}")]
    InvalidInterval(#[from] CoordinateError),
}

#[derive(Clone, Copy)]
struct ExpandedAnchor {
    weighted: WeightedAnchor,
    query_position: u64,
    source_index: usize,
}

/// Build up to `config.max_chains` chains per contig/strand group.
///
/// For a wrapped coordinate model, query positions are duplicated one query
/// length apart, allowing one origin crossing while rejecting paths that use
/// the same source anchor twice. Linear and undetermined models keep one copy
/// of each anchor, so they cannot emit origin-crossing chains. The predecessor
/// scan is capped by `max_predecessors`, so runtime and intermediate state
/// remain bounded by the caller's anchor limit.
pub fn chain_anchors(
    anchors: &[Anchor],
    query_length: u64,
    config: ChainConfig,
) -> Result<Vec<AnchorChain>, ChainError> {
    let weighted = anchors
        .iter()
        .copied()
        .map(|anchor| WeightedAnchor::new(anchor, class_for_legacy_anchor(anchor)))
        .collect::<Vec<_>>();
    Ok(chain_weighted_internal(&weighted, query_length, config)?
        .into_iter()
        .map(|chain| chain.base)
        .collect())
}

/// Build mixed weighted chains while retaining the seed classes used by the
/// scorer. At least one high-specificity anchor, or several consistent
/// lower-specificity anchors, is required for an accepted chain.
pub fn chain_weighted_anchors(
    anchors: &[WeightedAnchor],
    query_length: u64,
    config: ChainConfig,
) -> Result<Vec<WeightedAnchorChain>, ChainError> {
    Ok(chain_weighted_internal(anchors, query_length, config)?
        .into_iter()
        .map(|chain| chain.weighted)
        .collect())
}

/// Alias for callers that use the shorter mixed-chain operation name.
pub fn build_weighted_chains(
    anchors: &[WeightedAnchor],
    query_length: u64,
    config: ChainConfig,
) -> Result<Vec<WeightedAnchorChain>, ChainError> {
    chain_weighted_anchors(anchors, query_length, config)
}

fn chain_weighted_internal(
    anchors: &[WeightedAnchor],
    query_length: u64,
    config: ChainConfig,
) -> Result<Vec<WeightedChainResult>, ChainError> {
    validate_config(config, query_length)?;

    let mut validated = Vec::with_capacity(anchors.len());
    for (source_index, &weighted) in anchors.iter().enumerate() {
        let anchor = weighted.anchor;
        if anchor.k == 0 {
            return Err(ChainError::ZeroK);
        }
        if weighted.span < u64::from(anchor.k) {
            return Err(ChainError::InvalidAnchorSpan {
                span: weighted.span,
                k: anchor.k,
            });
        }
        let end = anchor
            .query_position
            .checked_add(weighted.span)
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
            .checked_add(weighted.span)
            .ok_or(ChainError::CoordinateOverflow)?;
        validated.push((source_index, weighted));
    }

    // Primary and rescue coordinates intentionally share a group. Class and
    // k-mer size influence score and tie breaks, but never prevent a mixed
    // monotone path from joining compatible evidence.
    validated.sort_by_key(|(_, weighted)| weighted_anchor_sort_key(weighted));
    let max_chains = usize::try_from(config.max_chains).unwrap_or(usize::MAX);
    let mut all_chains = Vec::with_capacity(validated.len().min(max_chains));
    let mut group_start = 0usize;
    while group_start < validated.len() {
        let current_group_key = group_key(&validated[group_start].1.anchor);
        let mut group_end = group_start + 1;
        while group_end < validated.len()
            && group_key(&validated[group_end].1.anchor) == current_group_key
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
            all_chains.push(chain);
            all_chains.sort_by(weighted_chain_sort_key);
            all_chains.truncate(max_chains);
        }
        group_start = group_end;
    }

    all_chains.sort_by(weighted_chain_sort_key);
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

struct WeightedChainResult {
    base: AnchorChain,
    weighted: WeightedAnchorChain,
    source_indices: Vec<usize>,
}

fn class_for_legacy_anchor(anchor: Anchor) -> AnchorClass {
    match anchor.k {
        crate::trace::seeds::PRIMARY_K => AnchorClass::SpecificK31,
        crate::trace::seeds::RESCUE_K => AnchorClass::SpecificK21,
        // Synthetic and future legacy anchors do not carry a class. Keep
        // their historical chain behavior; callers that know a seed is
        // repetitive should use `chain_weighted_anchors` explicitly.
        _ => AnchorClass::SpecificPaired,
    }
}

fn best_chain(
    group: &[(usize, WeightedAnchor)],
    used: &[bool],
    query_length: u64,
    config: ChainConfig,
) -> Result<Option<WeightedChainResult>, ChainError> {
    if group.is_empty() {
        return Ok(None);
    }

    let expanded_capacity = if permits_wrap(config.coordinate_model) {
        group
            .len()
            .checked_mul(2)
            .ok_or(ChainError::CoordinateOverflow)?
    } else {
        group.len()
    };
    let mut expanded = Vec::with_capacity(expanded_capacity);
    for (group_index, &(source_index, weighted)) in group.iter().enumerate() {
        if used[group_index] {
            continue;
        }
        let anchor = weighted.anchor;
        expanded.push(ExpandedAnchor {
            weighted,
            query_position: anchor.query_position,
            source_index: group_index,
        });
        // Keep `source_index` tied to the group position, not the caller's
        // original slice index, so duplicate circular copies can be blocked.
        let _ = source_index;
        if permits_wrap(config.coordinate_model) {
            let wrapped = anchor
                .query_position
                .checked_add(query_length)
                .ok_or(ChainError::CoordinateOverflow)?;
            expanded.push(ExpandedAnchor {
                weighted,
                query_position: wrapped,
                source_index: group_index,
            });
        }
    }
    expanded.sort_by(|left, right| expanded_cmp(*left, *right));

    let length = expanded.len();
    let mut scores = vec![i64::MIN; length];
    let mut counts = vec![0u32; length];
    let mut predecessors: Vec<Option<usize>> = vec![None; length];

    for i in 0..length {
        let current = expanded[i];
        scores[i] = anchor_base_score(current.weighted);
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
            let target_gap = match current.weighted.anchor.strand {
                Strand::Forward => {
                    if current.weighted.anchor.target_position
                        <= previous.weighted.anchor.target_position
                    {
                        continue;
                    }
                    current.weighted.anchor.target_position
                        - previous.weighted.anchor.target_position
                }
                Strand::Reverse => {
                    if current.weighted.anchor.target_position
                        >= previous.weighted.anchor.target_position
                    {
                        continue;
                    }
                    previous.weighted.anchor.target_position
                        - current.weighted.anchor.target_position
                }
            };
            if target_gap == 0 || target_gap > config.max_target_gap {
                continue;
            }
            if scores[j] == i64::MIN {
                continue;
            }
            let gap = query_gap.abs_diff(target_gap);
            let extension = gap_extension(current.weighted.class, gap, config.gap_penalty);
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
            if !path_has_required_evidence(index, &predecessors, &expanded, config) {
                continue;
            }
        } else if !path_has_required_evidence(index, &predecessors, &expanded, config) {
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
    let linear_start = first.query_position;
    let linear_end = path.iter().try_fold(linear_start, |end, &index| {
        let anchor_end = expanded[index]
            .query_position
            .checked_add(expanded[index].weighted.span)
            .ok_or(ChainError::CoordinateOverflow)?;
        Ok::<_, ChainError>(end.max(anchor_end))
    })?;
    let span = linear_end
        .checked_sub(linear_start)
        .ok_or(ChainError::CoordinateOverflow)?;
    if span > query_length {
        return Ok(None);
    }

    let (query_segments, origin_crossing) = if permits_wrap(config.coordinate_model) {
        let segments = circular_segments(linear_start, span, query_length)?;
        (segments, true)
    } else {
        let interval = BaseInterval::new(linear_start, linear_end)?;
        (vec![interval], false)
    };
    let origin_crossing = origin_crossing && query_segments.len() == 2;
    let query_interval = *query_segments
        .first()
        .expect("non-empty chain has a query segment");

    let (target_start, target_end) =
        path.iter()
            .try_fold((u64::MAX, 0u64), |(start, end), &index| {
                let anchor = expanded[index].weighted.anchor;
                let anchor_end = anchor
                    .target_position
                    .checked_add(expanded[index].weighted.span)
                    .ok_or(ChainError::CoordinateOverflow)?;
                Ok::<_, ChainError>((start.min(anchor.target_position), end.max(anchor_end)))
            })?;
    let target_interval = BaseInterval::new(target_start, target_end)?;
    let weighted_anchors = path
        .iter()
        .map(|&index| expanded[index].weighted)
        .collect::<Vec<_>>();
    let anchors = weighted_anchors
        .iter()
        .map(|weighted| weighted.anchor)
        .collect::<Vec<_>>();
    let source_indices = path
        .iter()
        .map(|&index| expanded[index].source_index)
        .collect();
    let score = scores[endpoint];
    Ok(Some(WeightedChainResult {
        base: AnchorChain {
            contig_id: first.weighted.anchor.contig_id,
            strand: first.weighted.anchor.strand,
            k: first.weighted.anchor.k,
            anchors,
            score,
            query_interval,
            query_segments: query_segments.clone(),
            target_interval,
            origin_crossing,
            linear_query_start: linear_start,
            linear_query_end: linear_end,
        },
        weighted: WeightedAnchorChain {
            contig_id: first.weighted.anchor.contig_id,
            strand: first.weighted.anchor.strand,
            anchors: weighted_anchors,
            score,
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

fn anchor_base_score(anchor: WeightedAnchor) -> i64 {
    anchor.class.score()
}

fn gap_extension(class: AnchorClass, gap: u64, gap_penalty: i64) -> i64 {
    let gap = i64::try_from(gap).unwrap_or(i64::MAX);
    class
        .score()
        .saturating_sub(gap.saturating_mul(gap_penalty))
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
    let candidate_anchor = expanded[candidate].weighted.anchor;
    let current_anchor = expanded[current].weighted.anchor;
    (
        expanded[candidate].query_position,
        candidate_anchor.target_position,
        candidate_anchor.hash,
        expanded[candidate].weighted.class,
    ) < (
        expanded[current].query_position,
        current_anchor.target_position,
        current_anchor.hash,
        expanded[current].weighted.class,
    )
}

fn path_has_required_evidence(
    endpoint: usize,
    predecessors: &[Option<usize>],
    expanded: &[ExpandedAnchor],
    config: ChainConfig,
) -> bool {
    let mut lower_specificity = 0u32;
    let mut cursor = Some(endpoint);
    while let Some(index) = cursor {
        let class = expanded[index].weighted.class;
        if class.high_specificity() {
            return true;
        }
        if class.lower_specificity() && !matches!(class, AnchorClass::Repetitive) {
            lower_specificity = lower_specificity.saturating_add(1);
        }
        cursor = predecessors[index];
    }
    lower_specificity >= config.min_anchors.max(2)
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

#[must_use]
const fn permits_wrap(coordinate_model: CoordinateModel) -> bool {
    matches!(coordinate_model, CoordinateModel::Wrap)
}

fn group_key(anchor: &Anchor) -> (u32, u8) {
    (anchor.contig_id, strand_rank(anchor.strand))
}

fn weighted_anchor_sort_key(
    weighted: &WeightedAnchor,
) -> (u32, u8, u64, u64, u8, AnchorClass, u64, u64) {
    let anchor = weighted.anchor;
    (
        anchor.contig_id,
        strand_rank(anchor.strand),
        anchor.query_position,
        anchor.target_position,
        anchor.k,
        weighted.class,
        anchor.hash,
        anchor.canonical_kmer,
    )
}

fn expanded_cmp(left: ExpandedAnchor, right: ExpandedAnchor) -> std::cmp::Ordering {
    left.query_position
        .cmp(&right.query_position)
        .then_with(|| match left.weighted.anchor.strand {
            Strand::Forward => left
                .weighted
                .anchor
                .target_position
                .cmp(&right.weighted.anchor.target_position),
            Strand::Reverse => right
                .weighted
                .anchor
                .target_position
                .cmp(&left.weighted.anchor.target_position),
        })
        .then_with(|| left.weighted.class.cmp(&right.weighted.class))
        .then_with(|| left.weighted.anchor.hash.cmp(&right.weighted.anchor.hash))
        .then_with(|| left.source_index.cmp(&right.source_index))
}

fn weighted_chain_sort_key(
    left: &WeightedChainResult,
    right: &WeightedChainResult,
) -> std::cmp::Ordering {
    let left_base = &left.base;
    let right_base = &right.base;
    right_base
        .score
        .cmp(&left_base.score)
        .then_with(|| right_base.anchors.len().cmp(&left_base.anchors.len()))
        .then_with(|| left_base.contig_id.cmp(&right_base.contig_id))
        .then_with(|| strand_rank(left_base.strand).cmp(&strand_rank(right_base.strand)))
        .then_with(|| left_base.k.cmp(&right_base.k))
        .then_with(|| {
            left_base
                .linear_query_start
                .cmp(&right_base.linear_query_start)
        })
        .then_with(|| {
            left_base
                .target_interval
                .start
                .cmp(&right_base.target_interval.start)
        })
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
            coordinate_model: CoordinateModel::Wrap,
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

    #[test]
    fn linear_end_covers_an_earlier_long_gear_run() {
        let long_run = WeightedAnchor::new(
            Anchor {
                k: 31,
                ..anchor(10, 110, Strand::Forward)
            },
            AnchorClass::ExactGearRun,
        )
        .with_span(300);
        let later_seed = WeightedAnchor::new(
            Anchor {
                k: 21,
                ..anchor(50, 150, Strand::Forward)
            },
            AnchorClass::SpecificK21,
        );
        let mut chain_config = config();
        chain_config.coordinate_model = CoordinateModel::Linear;
        let chains = chain_weighted_anchors(&[long_run, later_seed], 1_000, chain_config).unwrap();
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].linear_query_start, 10);
        assert_eq!(chains[0].linear_query_end, 310);
        assert_eq!(
            chains[0].query_interval,
            BaseInterval::new(10, 310).unwrap()
        );
    }
}
