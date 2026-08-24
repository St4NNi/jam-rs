//! Bounded refinement helpers around the existing corridor aligner.
//!
//! This module deliberately does not duplicate the affine-gap matrix.  It
//! supplies the state machine needed by a score-only retry integration: probe
//! widths first, choose one deterministic final width, and invoke traceback
//! exactly once.  It also owns conservative chain-bridge and single-specific
//! anchor gates used by the trace runner.

use super::{AlignmentError, AlignmentOptions, AlignmentResult, AlignmentWorkspace};
use crate::trace::alignment::exact_blocks::{ExactBlock, ExactBlockError, extend_anchor};
use crate::trace::anchors::Anchor;
use crate::trace::config::AlignmentScoring;
use crate::trace::model::{BaseInterval, Strand};
use std::cmp::min;
use thiserror::Error;

/// Result of one score-only corridor pass.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ScoreProbe {
    pub score: i32,
    pub band_edge_touched: bool,
    pub under_explained: bool,
    pub predicted_drift: bool,
}

#[derive(Clone, Copy, Debug, Default)]
struct RollingPath {
    query_start: usize,
    target_start: usize,
    edge_touched: bool,
}

#[derive(Clone, Copy, Debug, Default)]
struct RollingCell {
    scores: [i32; 3],
    paths: [RollingPath; 3],
}

/// Reusable rolling-row storage for score-only corridor passes.
///
/// The two vectors are bounded by the widest corridor row, not by the number
/// of query bases or total DP cells.  Their capacity is intentionally exposed
/// so callers can include this scratch in their memory estimate.
#[derive(Debug, Default)]
pub struct ScoreOnlyWorkspace {
    previous: Vec<RollingCell>,
    current: Vec<RollingCell>,
    max_row_width: usize,
    rows_processed: usize,
}

impl ScoreOnlyWorkspace {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of rolling cells reserved across both rows.
    #[must_use]
    pub fn capacity_cells(&self) -> usize {
        self.previous
            .capacity()
            .saturating_add(self.current.capacity())
    }

    /// Widest row observed by the most recent pass.
    #[must_use]
    pub const fn max_row_width(&self) -> usize {
        self.max_row_width
    }

    /// Number of query rows processed by the most recent pass.
    #[must_use]
    pub const fn rows_processed(&self) -> usize {
        self.rows_processed
    }
}

/// Score a local affine alignment through a piecewise corridor without
/// allocating traceback cells.
///
/// This is numerically equivalent to the score recurrence in
/// `AlignmentWorkspace::align_raw_with_layout`: it uses the same local
/// zero-reset, affine open/extend scores, transition tie order, and corridor
/// bounds. Only two variable-width rows and per-state path-start metadata are
/// retained. `ScoreProbe::under_explained` is derived from the winning path's
/// query/target span in the same way as corridor refinement.
pub fn score_only_corridor(
    query: &[u8],
    target: &[u8],
    corridor: &super::AlignmentCorridor,
    options: impl Into<AlignmentOptions>,
) -> Result<ScoreProbe, AlignmentError> {
    let mut scratch = ScoreOnlyWorkspace::new();
    score_only_corridor_with_scratch(&mut scratch, query, target, corridor, options)
}

/// Reuse caller-owned rolling rows for one score-only corridor pass.
pub fn score_only_corridor_with_scratch(
    scratch: &mut ScoreOnlyWorkspace,
    query: &[u8],
    target: &[u8],
    corridor: &super::AlignmentCorridor,
    options: impl Into<AlignmentOptions>,
) -> Result<ScoreProbe, AlignmentError> {
    let options = options.into();
    if query.is_empty() {
        return Err(AlignmentError::EmptyQuery);
    }
    if target.is_empty() {
        return Err(AlignmentError::EmptyTarget);
    }

    let query_len = query.len();
    let target_len = target.len();
    let mut total_cells = 0usize;
    let mut widest_row = 0usize;
    for query_index in 0..=query_len {
        let (start, end) = corridor.bounds(
            query_index,
            target_len,
            options.scoring.band_width,
            options.diagonal_offset,
        )?;
        let width = if start > end { 0 } else { end - start + 1 };
        total_cells = total_cells
            .checked_add(width)
            .ok_or(AlignmentError::LengthOverflow)?;
        widest_row = widest_row.max(width);
    }
    if total_cells == 0 {
        return Err(AlignmentError::BandExcludesInput);
    }
    if options.max_cells == 0 || total_cells > options.max_cells {
        return Err(AlignmentError::MatrixTooLarge {
            cells: total_cells,
            max_cells: options.max_cells,
        });
    }

    scratch.previous.resize(widest_row, RollingCell::default());
    scratch.current.resize(widest_row, RollingCell::default());
    scratch.max_row_width = widest_row;
    scratch.rows_processed = 0;

    let mut previous_start = 0usize;
    let mut previous_width = 0usize;
    let mut best_score = 0i32;
    let mut best_query_index = 0usize;
    let mut best_target_index = 0usize;
    let mut best_state = super::STATE_UNREACHABLE;
    let mut best_path = RollingPath::default();

    for query_index in 0..=query_len {
        let (start, end) = corridor.bounds(
            query_index,
            target_len,
            options.scoring.band_width,
            options.diagonal_offset,
        )?;
        let (row_start, row_width) = if start > end {
            (0, 0)
        } else {
            (start, end - start + 1)
        };
        scratch.current[..row_width].fill(RollingCell::default());
        for offset in 0..row_width {
            let target_index = row_start + offset;
            if query_index == 0 && target_index == 0 {
                continue;
            }

            let mut cell = RollingCell::default();
            let edge_here = row_width > 0
                && ((target_index == row_start && row_start > 0)
                    || (target_index == row_start + row_width - 1 && target_index < target_len));

            if query_index > 0
                && target_index > 0
                && let Some(previous) = rolling_cell(
                    &scratch.previous,
                    previous_start,
                    previous_width,
                    target_index - 1,
                )
            {
                let (previous_score, previous_state) = rolling_best(previous);
                let substitution =
                    if super::equal_base(query[query_index - 1], target[target_index - 1]) {
                        options.scoring.match_score
                    } else {
                        options.scoring.mismatch_score
                    };
                let score = previous_score.saturating_add(substitution);
                if score > 0 {
                    let path = if previous_score > 0 {
                        previous.paths[usize::from(previous_state)]
                    } else {
                        RollingPath {
                            query_start: query_index - 1,
                            target_start: target_index - 1,
                            edge_touched: false,
                        }
                    };
                    cell.scores[usize::from(super::STATE_MATCH)] = score;
                    cell.paths[usize::from(super::STATE_MATCH)] = RollingPath {
                        edge_touched: path.edge_touched || edge_here,
                        ..path
                    };
                }
            }

            if target_index > 0
                && let Some(previous) =
                    rolling_cell(&scratch.current, row_start, row_width, target_index - 1)
            {
                let extension = previous.scores[usize::from(super::STATE_INSERTION)]
                    .saturating_add(options.scoring.gap_extend_score);
                let opened_from_match = previous.scores[usize::from(super::STATE_MATCH)]
                    .saturating_add(super::gap_open_score(options.scoring));
                let opened_from_deletion = previous.scores[usize::from(super::STATE_DELETION)]
                    .saturating_add(super::gap_open_score(options.scoring));
                let (score, previous_state) = super::choose_transition([
                    (extension, super::STATE_INSERTION),
                    (opened_from_match, super::STATE_MATCH),
                    (opened_from_deletion, super::STATE_DELETION),
                ]);
                if score > 0 {
                    let previous_score = previous.scores[usize::from(previous_state)];
                    let path = if previous_score > 0 {
                        previous.paths[usize::from(previous_state)]
                    } else {
                        RollingPath {
                            query_start: query_index,
                            target_start: target_index - 1,
                            edge_touched: false,
                        }
                    };
                    cell.scores[usize::from(super::STATE_INSERTION)] = score;
                    cell.paths[usize::from(super::STATE_INSERTION)] = RollingPath {
                        edge_touched: path.edge_touched || edge_here,
                        ..path
                    };
                }
            }

            if query_index > 0
                && let Some(previous) = rolling_cell(
                    &scratch.previous,
                    previous_start,
                    previous_width,
                    target_index,
                )
            {
                let extension = previous.scores[usize::from(super::STATE_DELETION)]
                    .saturating_add(options.scoring.gap_extend_score);
                let opened_from_match = previous.scores[usize::from(super::STATE_MATCH)]
                    .saturating_add(super::gap_open_score(options.scoring));
                let opened_from_insertion = previous.scores[usize::from(super::STATE_INSERTION)]
                    .saturating_add(super::gap_open_score(options.scoring));
                let (score, previous_state) = super::choose_transition([
                    (extension, super::STATE_DELETION),
                    (opened_from_match, super::STATE_MATCH),
                    (opened_from_insertion, super::STATE_INSERTION),
                ]);
                if score > 0 {
                    let previous_score = previous.scores[usize::from(previous_state)];
                    let path = if previous_score > 0 {
                        previous.paths[usize::from(previous_state)]
                    } else {
                        RollingPath {
                            query_start: query_index - 1,
                            target_start: target_index,
                            edge_touched: false,
                        }
                    };
                    cell.scores[usize::from(super::STATE_DELETION)] = score;
                    cell.paths[usize::from(super::STATE_DELETION)] = RollingPath {
                        edge_touched: path.edge_touched || edge_here,
                        ..path
                    };
                }
            }

            scratch.current[offset] = cell;
            let (score, state) = rolling_best(cell);
            if score > best_score
                || (score == best_score
                    && score > 0
                    && (query_index < best_query_index
                        || (query_index == best_query_index
                            && (target_index < best_target_index
                                || (target_index == best_target_index && state < best_state)))))
            {
                best_score = score;
                best_query_index = query_index;
                best_target_index = target_index;
                best_state = state;
                best_path = cell.paths[usize::from(state)];
            }
        }

        std::mem::swap(&mut scratch.previous, &mut scratch.current);
        previous_start = row_start;
        previous_width = row_width;
        scratch.rows_processed = scratch.rows_processed.saturating_add(1);
    }

    if best_score <= 0 {
        return Err(AlignmentError::NoAlignment);
    }
    let query_span = best_query_index.saturating_sub(best_path.query_start) as u64;
    let target_span = best_target_index.saturating_sub(best_path.target_start) as u64;
    let margin = u64::from(corridor.safety_margin()).saturating_mul(2);
    let under_explained = query_span.saturating_add(margin) < corridor.query_span()
        || target_span.saturating_add(margin) < corridor.target_span();
    let predicted_drift = corridor.required_half_width() > options.scoring.band_width;
    Ok(ScoreProbe {
        score: best_score,
        band_edge_touched: best_path.edge_touched,
        under_explained,
        predicted_drift,
    })
}

fn rolling_cell(
    row: &[RollingCell],
    start: usize,
    width: usize,
    target_index: usize,
) -> Option<RollingCell> {
    if target_index < start || target_index >= start.checked_add(width)? {
        return None;
    }
    row.get(target_index - start).copied()
}

fn rolling_best(cell: RollingCell) -> (i32, u8) {
    super::choose_transition([
        (
            cell.scores[usize::from(super::STATE_MATCH)],
            super::STATE_MATCH,
        ),
        (
            cell.scores[usize::from(super::STATE_INSERTION)],
            super::STATE_INSERTION,
        ),
        (
            cell.scores[usize::from(super::STATE_DELETION)],
            super::STATE_DELETION,
        ),
    ])
}

/// Deterministic choice made after score-only corridor probes.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RetrySelection {
    pub attempted_widths: Vec<u32>,
    pub selected_width: u32,
    pub retries: u8,
    pub band_edge_touched: bool,
    pub chain_span_under_explained: bool,
    pub retry_cap_reached: bool,
}

/// Alignment returned after a score-only selection and one traceback call.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RefinedAlignment {
    pub result: AlignmentResult,
    pub selection: RetrySelection,
}

/// Select a corridor width without constructing traceback output.
///
/// `score_only` is called at most once per unique width, in the order supplied
/// by the caller. The first probe that does not touch an edge, under-explain
/// the chain span, or predict drift is selected. If every probe is suspect,
/// the final probe is retained and explicitly marked as capped.
pub fn select_retry_width<F>(
    widths: &[u32],
    score_only: F,
) -> Result<RetrySelection, AlignmentError>
where
    F: FnMut(u32) -> Result<ScoreProbe, AlignmentError>,
{
    let mut score_only = score_only;
    let mut attempted_widths = Vec::new();
    let mut selected = None;
    let mut aggregate_edge = false;
    let mut aggregate_under_explained = false;
    let mut last_probe = None;

    for &width in widths {
        if width == 0 || attempted_widths.contains(&width) {
            continue;
        }
        let probe = score_only(width)?;
        attempted_widths.push(width);
        aggregate_edge |= probe.band_edge_touched;
        aggregate_under_explained |= probe.under_explained;
        last_probe = Some(probe);
        if !probe.band_edge_touched && !probe.under_explained && !probe.predicted_drift {
            selected = Some(width);
            break;
        }
    }

    let Some(last_probe) = last_probe else {
        return Err(AlignmentError::NoAlignment);
    };
    let selected_width =
        selected.unwrap_or_else(|| *attempted_widths.last().expect("probe exists"));
    let retry_cap_reached = selected.is_none()
        && (last_probe.band_edge_touched
            || last_probe.under_explained
            || last_probe.predicted_drift);
    Ok(RetrySelection {
        retries: u8::try_from(attempted_widths.len().saturating_sub(1)).unwrap_or(u8::MAX),
        attempted_widths,
        selected_width,
        band_edge_touched: aggregate_edge,
        chain_span_under_explained: aggregate_under_explained,
        retry_cap_reached,
    })
}

/// Probe widths first and run `traceback` exactly once at the selected width.
pub fn refine_with_score_only<F, G>(
    widths: &[u32],
    score_only: F,
    mut traceback: G,
) -> Result<RefinedAlignment, AlignmentError>
where
    F: FnMut(u32) -> Result<ScoreProbe, AlignmentError>,
    G: FnMut(u32) -> Result<AlignmentResult, AlignmentError>,
{
    let selection = select_retry_width(widths, score_only)?;
    let mut result = traceback(selection.selected_width)?;
    result.retry_metadata.attempted_widths = selection.attempted_widths.clone();
    result.retry_metadata.selected_width = selection.selected_width;
    result.retry_metadata.retries = selection.retries;
    result.retry_metadata.band_edge_touched |= selection.band_edge_touched;
    result.retry_metadata.chain_span_under_explained |= selection.chain_span_under_explained;
    result.retry_metadata.retry_cap_reached = selection.retry_cap_reached;
    Ok(RefinedAlignment { result, selection })
}

/// Compatibility wrapper for the current public workspace API.  It preserves
/// the existing piecewise-corridor behavior while callers migrate to the
/// score-only callback above.
pub fn align_with_existing_corridor(
    workspace: &mut AlignmentWorkspace,
    query: &[u8],
    target: &[u8],
    corridor: &super::AlignmentCorridor,
    options: impl Into<AlignmentOptions>,
) -> Result<AlignmentResult, AlignmentError> {
    workspace.align_with_retries(query, target, corridor, options)
}

/// Endpoint state used by bounded chain bridging.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ChainEndpoint {
    pub contig_id: u32,
    pub strand: Strand,
    pub query_interval: BaseInterval,
    /// Forward-stored target interval, even for reverse-strand chains.
    pub target_interval: BaseInterval,
    pub support_bases: u64,
}

/// Limits for a chain-bridge attempt.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BridgeConfig {
    pub max_query_gap: u64,
    pub max_target_gap: u64,
    pub max_gap_difference: u64,
    pub min_flank_support: u64,
}

impl Default for BridgeConfig {
    fn default() -> Self {
        Self {
            max_query_gap: 250_000,
            max_target_gap: 250_000,
            max_gap_difference: 2_048,
            min_flank_support: 21,
        }
    }
}

/// A bounded, same-contig/same-strand bridge between two retained chains.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ChainBridge {
    pub contig_id: u32,
    pub strand: Strand,
    pub query_gap: BaseInterval,
    pub target_gap: BaseInterval,
    pub query_gap_bases: u64,
    pub target_gap_bases: u64,
}

/// Validate and construct one bridge.  The target length is required only to
/// compare reverse-oriented chain coordinates in their aligned order.
#[must_use]
pub fn bridge_chains(
    left: ChainEndpoint,
    right: ChainEndpoint,
    target_length: u64,
    config: BridgeConfig,
) -> Option<ChainBridge> {
    if left.contig_id != right.contig_id
        || left.strand != right.strand
        || left.support_bases < config.min_flank_support
        || right.support_bases < config.min_flank_support
        || left.query_interval.end > right.query_interval.start
    {
        return None;
    }
    let (_left_target_start, left_target_end) = oriented_target_interval(left, target_length)?;
    let (right_target_start, _right_target_end) = oriented_target_interval(right, target_length)?;
    if left_target_end > right_target_start {
        return None;
    }
    let query_gap = BaseInterval::new(left.query_interval.end, right.query_interval.start).ok()?;
    let target_gap_oriented = BaseInterval::new(left_target_end, right_target_start).ok()?;
    let query_gap_bases = query_gap.len();
    let target_gap_bases = target_gap_oriented.len();
    if query_gap_bases > config.max_query_gap
        || target_gap_bases > config.max_target_gap
        || query_gap_bases.abs_diff(target_gap_bases) > config.max_gap_difference
    {
        return None;
    }
    let target_gap = match left.strand {
        Strand::Forward => target_gap_oriented,
        Strand::Reverse => BaseInterval::new(
            target_length.checked_sub(target_gap_oriented.end)?,
            target_length.checked_sub(target_gap_oriented.start)?,
        )
        .ok()?,
    };
    Some(ChainBridge {
        contig_id: left.contig_id,
        strand: left.strand,
        query_gap,
        target_gap,
        query_gap_bases,
        target_gap_bases,
    })
}

/// Build deterministic bounded bridges between endpoint pairs.
#[must_use]
pub fn bridge_endpoint_pairs(
    endpoints: &[ChainEndpoint],
    target_length: u64,
    config: BridgeConfig,
    max_bridges: usize,
) -> Vec<ChainBridge> {
    if max_bridges == 0 {
        return Vec::new();
    }
    let mut ordered = endpoints.to_vec();
    ordered.sort_by_key(|endpoint| {
        (
            endpoint.contig_id,
            strand_rank(endpoint.strand),
            endpoint.query_interval.start,
            endpoint.target_interval.start,
            endpoint.target_interval.end,
        )
    });
    let mut bridges = Vec::new();
    for (left_index, left) in ordered.iter().copied().enumerate() {
        for right in ordered.iter().copied().skip(left_index + 1) {
            if right.contig_id != left.contig_id || right.strand != left.strand {
                continue;
            }
            if right.query_interval.start < left.query_interval.end {
                continue;
            }
            if right.query_interval.start - left.query_interval.end > config.max_query_gap {
                break;
            }
            if let Some(bridge) = bridge_chains(left, right, target_length, config) {
                bridges.push(bridge);
                if bridges.len() == max_bridges {
                    return bridges;
                }
            }
        }
    }
    bridges
}

/// Strict admission settings for a one-specific-anchor extension.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SingleAnchorConfig {
    pub max_occurrences: u32,
    pub max_window_bases: u64,
    pub min_identity: f64,
    pub min_anchor_k: u8,
}

impl Default for SingleAnchorConfig {
    fn default() -> Self {
        Self {
            max_occurrences: 2,
            max_window_bases: 512,
            min_identity: 0.99,
            min_anchor_k: 21,
        }
    }
}

/// Accepted single-anchor refinement, retaining the verified exact block and
/// the stricter local alignment result separately.
#[derive(Clone, Debug, PartialEq)]
pub struct SingleAnchorResult {
    pub block: ExactBlock,
    pub alignment: AlignmentResult,
    pub query_window: BaseInterval,
    pub target_window: BaseInterval,
}

/// Attempt one bounded, rare/low-copy anchor extension.
pub fn try_single_anchor(
    query: &[u8],
    target_forward: &[u8],
    anchor: Anchor,
    occurrence_count: u32,
    scoring: AlignmentScoring,
    config: SingleAnchorConfig,
) -> Result<Option<SingleAnchorResult>, SingleAnchorError> {
    validate_single_config(config)?;
    if occurrence_count > config.max_occurrences || anchor.k < config.min_anchor_k {
        return Ok(None);
    }
    let max_window = usize::try_from(config.max_window_bases)
        .map_err(|_| SingleAnchorError::CoordinateOverflow)?;
    let k = usize::from(anchor.k);
    if max_window < k {
        return Ok(None);
    }
    let block = extend_anchor(query, target_forward, anchor)?;
    let query_anchor_start = usize::try_from(anchor.query_position)
        .map_err(|_| SingleAnchorError::CoordinateOverflow)?;
    let query_anchor_end = query_anchor_start
        .checked_add(k)
        .ok_or(SingleAnchorError::CoordinateOverflow)?;
    let query_window = bounded_window(
        query_anchor_start,
        query_anchor_end,
        query.len(),
        max_window,
    )
    .ok_or(SingleAnchorError::CoordinateOverflow)?;

    let (oriented_target_start, oriented_target_end) = oriented_target_interval(
        ChainEndpoint {
            contig_id: anchor.contig_id,
            strand: anchor.strand,
            query_interval: BaseInterval::new(
                anchor.query_position,
                anchor
                    .query_position
                    .checked_add(u64::from(anchor.k))
                    .ok_or(SingleAnchorError::CoordinateOverflow)?,
            )
            .map_err(|_| SingleAnchorError::CoordinateOverflow)?,
            target_interval: BaseInterval::new(
                anchor.target_position,
                anchor
                    .target_position
                    .checked_add(u64::from(anchor.k))
                    .ok_or(SingleAnchorError::CoordinateOverflow)?,
            )
            .map_err(|_| SingleAnchorError::CoordinateOverflow)?,
            support_bases: u64::from(anchor.k),
        },
        target_forward.len() as u64,
    )
    .ok_or(SingleAnchorError::CoordinateOverflow)?;
    let oriented_window = bounded_window(
        usize::try_from(oriented_target_start)
            .map_err(|_| SingleAnchorError::CoordinateOverflow)?,
        usize::try_from(oriented_target_end).map_err(|_| SingleAnchorError::CoordinateOverflow)?,
        target_forward.len(),
        max_window,
    )
    .ok_or(SingleAnchorError::CoordinateOverflow)?;
    let target_window = match anchor.strand {
        Strand::Forward => BaseInterval::new(oriented_window.0 as u64, oriented_window.1 as u64),
        Strand::Reverse => BaseInterval::new(
            (target_forward.len() - oriented_window.1) as u64,
            (target_forward.len() - oriented_window.0) as u64,
        ),
    }
    .map_err(|_| SingleAnchorError::CoordinateOverflow)?;

    let query_slice = &query[query_window.0..query_window.1];
    let target_start =
        usize::try_from(target_window.start).map_err(|_| SingleAnchorError::CoordinateOverflow)?;
    let target_end =
        usize::try_from(target_window.end).map_err(|_| SingleAnchorError::CoordinateOverflow)?;
    let target_slice = target_forward
        .get(target_start..target_end)
        .ok_or(SingleAnchorError::CoordinateOverflow)?;
    let mut workspace = AlignmentWorkspace::new();
    let mut alignment = workspace.align_oriented(
        query_slice,
        target_slice,
        target_window.start,
        anchor.strand,
        scoring,
    )?;
    alignment.query_interval = BaseInterval::new(
        alignment
            .query_interval
            .start
            .checked_add(query_window.0 as u64)
            .ok_or(SingleAnchorError::CoordinateOverflow)?,
        alignment
            .query_interval
            .end
            .checked_add(query_window.0 as u64)
            .ok_or(SingleAnchorError::CoordinateOverflow)?,
    )
    .map_err(|_| SingleAnchorError::CoordinateOverflow)?;
    alignment.query_segments = vec![alignment.query_interval];
    let denominator = alignment
        .matches
        .saturating_add(alignment.substitutions)
        .saturating_add(alignment.insertions)
        .saturating_add(alignment.deletions);
    let identity = if denominator == 0 {
        0.0
    } else {
        alignment.matches as f64 / denominator as f64
    };
    let anchor_query = BaseInterval::new(
        anchor.query_position,
        anchor
            .query_position
            .checked_add(u64::from(anchor.k))
            .ok_or(SingleAnchorError::CoordinateOverflow)?,
    )
    .map_err(|_| SingleAnchorError::CoordinateOverflow)?;
    let anchor_target = BaseInterval::new(
        anchor.target_position,
        anchor
            .target_position
            .checked_add(u64::from(anchor.k))
            .ok_or(SingleAnchorError::CoordinateOverflow)?,
    )
    .map_err(|_| SingleAnchorError::CoordinateOverflow)?;
    if identity < config.min_identity
        || alignment.query_interval.start > anchor_query.start
        || alignment.query_interval.end < anchor_query.end
        || alignment.target_interval.start > anchor_target.start
        || alignment.target_interval.end < anchor_target.end
    {
        return Ok(None);
    }
    Ok(Some(SingleAnchorResult {
        block,
        alignment,
        query_window: BaseInterval::new(query_window.0 as u64, query_window.1 as u64)
            .map_err(|_| SingleAnchorError::CoordinateOverflow)?,
        target_window,
    }))
}

fn validate_single_config(config: SingleAnchorConfig) -> Result<(), SingleAnchorError> {
    if config.max_occurrences == 0 {
        return Err(SingleAnchorError::InvalidConfig(
            "max_occurrences must be > 0",
        ));
    }
    if config.max_window_bases == 0 {
        return Err(SingleAnchorError::InvalidConfig(
            "max_window_bases must be > 0",
        ));
    }
    if !(0.0..=1.0).contains(&config.min_identity) {
        return Err(SingleAnchorError::InvalidConfig(
            "min_identity must be in [0, 1]",
        ));
    }
    Ok(())
}

fn bounded_window(
    start: usize,
    end: usize,
    length: usize,
    max_width: usize,
) -> Option<(usize, usize)> {
    if start > end || end > length || max_width == 0 {
        return None;
    }
    let width = min(max_width, length);
    if end - start >= width {
        return Some((start, min(end, start.saturating_add(width))));
    }
    let flank = width.saturating_sub(end - start) / 2;
    let mut window_start = start.saturating_sub(flank);
    let mut window_end = min(length, window_start.saturating_add(width));
    if window_end < end {
        window_end = end;
        window_start = window_end.saturating_sub(width);
    }
    if window_end - window_start < end - start {
        return None;
    }
    Some((window_start, window_end))
}

fn oriented_target_interval(endpoint: ChainEndpoint, target_length: u64) -> Option<(u64, u64)> {
    if endpoint.target_interval.start > endpoint.target_interval.end
        || endpoint.target_interval.end > target_length
    {
        return None;
    }
    match endpoint.strand {
        Strand::Forward => Some((endpoint.target_interval.start, endpoint.target_interval.end)),
        Strand::Reverse => Some((
            target_length.checked_sub(endpoint.target_interval.end)?,
            target_length.checked_sub(endpoint.target_interval.start)?,
        )),
    }
}

const fn strand_rank(strand: Strand) -> u8 {
    match strand {
        Strand::Forward => 0,
        Strand::Reverse => 1,
    }
}

/// Errors from conservative single-anchor refinement.
#[derive(Debug, Error)]
pub enum SingleAnchorError {
    #[error(transparent)]
    Alignment(#[from] AlignmentError),
    #[error(transparent)]
    Exact(#[from] ExactBlockError),
    #[error("single-anchor coordinate overflow")]
    CoordinateOverflow,
    #[error("invalid single-anchor configuration: {0}")]
    InvalidConfig(&'static str),
}
