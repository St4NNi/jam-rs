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
use crate::trace::model::{BaseInterval, EditOperation, EditRun, Strand};
use std::cmp::min;
use std::fmt::Write as _;
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

/// Limits for stitching exact blocks and bounded affine-gap segments.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ExactBlockRefinementConfig {
    /// Maximum query or oriented-target span between two exact blocks.
    pub max_gap_bases: u64,
    /// Maximum terminal query or oriented-target flank to align. Larger
    /// flanks are left outside the local result rather than expanding DP.
    pub max_terminal_flank_bases: u64,
}

impl Default for ExactBlockRefinementConfig {
    fn default() -> Self {
        Self {
            max_gap_bases: 250_000,
            max_terminal_flank_bases: 2_048,
        }
    }
}

/// Refine one ordered exact-block path while bypassing DP for every verified
/// exact block. DP is invoked only for bounded intervening/terminal segments.
///
/// Blocks must describe one contig and one relative strand. Reverse blocks
/// retain forward target coordinates, while their oriented intervals define
/// the sequence order used by the stitcher. If a gap cannot be represented by
/// the bounded local workspace, callers should fall back to the existing
/// corridor path.
pub fn refine_with_exact_blocks(
    query: &[u8],
    target_forward: &[u8],
    blocks: &[ExactBlock],
    options: impl Into<AlignmentOptions>,
    config: ExactBlockRefinementConfig,
) -> Result<AlignmentResult, ExactBlockRefinementError> {
    let options = options.into();
    if query.is_empty() {
        return Err(ExactBlockRefinementError::Alignment(
            AlignmentError::EmptyQuery,
        ));
    }
    if target_forward.is_empty() {
        return Err(ExactBlockRefinementError::Alignment(
            AlignmentError::EmptyTarget,
        ));
    }
    if blocks.is_empty() {
        return Err(ExactBlockRefinementError::NoBlocks);
    }

    let mut ordered = blocks.to_vec();
    validate_exact_block_path(query.len(), target_forward.len(), &ordered)?;
    // Adjacent/overlapping anchors are common. Normalizing here makes the
    // direct run unambiguous while retaining the input path's contig/strand
    // validation above.
    let normalized = crate::trace::alignment::exact_blocks::merge_exact_blocks(&mut ordered);
    validate_exact_block_path(query.len(), target_forward.len(), &normalized)?;

    let strand = normalized[0].strand;
    let oriented_target = if strand == Strand::Forward {
        target_forward.to_vec()
    } else {
        target_forward
            .iter()
            .rev()
            .copied()
            .map(super::complement)
            .collect::<Vec<_>>()
    };

    let first = normalized[0].clone();
    let last = normalized.last().expect("validated non-empty path").clone();
    let terminal_limit = usize::try_from(config.max_terminal_flank_bases)
        .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
    let first_query_start = usize::try_from(first.query_interval.start)
        .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
    let first_target_start = usize::try_from(first.oriented_target_interval.start)
        .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
    let last_query_end = usize::try_from(last.query_interval.end)
        .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
    let last_target_end = usize::try_from(last.oriented_target_interval.end)
        .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;

    let query_start = if first_query_start <= terminal_limit {
        0
    } else {
        first_query_start
    };
    let target_start = if first_target_start <= terminal_limit {
        0
    } else {
        first_target_start
    };
    let query_end = if query.len().saturating_sub(last_query_end) <= terminal_limit {
        query.len()
    } else {
        last_query_end
    };
    let target_end = if oriented_target.len().saturating_sub(last_target_end) <= terminal_limit {
        oriented_target.len()
    } else {
        last_target_end
    };

    let mut operations = Vec::<EditRun>::new();
    let mut workspace = AlignmentWorkspace::new();
    let mut metadata = super::AlignmentRetryMetadata {
        corridor_used: true,
        ..super::AlignmentRetryMetadata::default()
    };
    append_bounded_segment(
        query,
        &oriented_target,
        query_start,
        first_query_start,
        target_start,
        first_target_start,
        &mut workspace,
        options,
        config.max_terminal_flank_bases,
        &mut operations,
        &mut metadata,
    )?;
    for pair in normalized.windows(2) {
        let left = &pair[0];
        let right = &pair[1];
        let left_query_end = usize::try_from(left.query_interval.end)
            .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
        let right_query_start = usize::try_from(right.query_interval.start)
            .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
        let left_target_end = usize::try_from(left.oriented_target_interval.end)
            .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
        let right_target_start = usize::try_from(right.oriented_target_interval.start)
            .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
        append_exact_run(&mut operations, left.len())?;
        append_bounded_segment(
            query,
            &oriented_target,
            left_query_end,
            right_query_start,
            left_target_end,
            right_target_start,
            &mut workspace,
            options,
            config.max_gap_bases,
            &mut operations,
            &mut metadata,
        )?;
    }

    append_exact_run(&mut operations, last.len())?;
    append_bounded_segment(
        query,
        &oriented_target,
        last_query_end,
        query_end,
        last_target_end,
        target_end,
        &mut workspace,
        options,
        config.max_terminal_flank_bases,
        &mut operations,
        &mut metadata,
    )?;
    let (query_consumed, target_consumed) = run_consumption(&operations);
    let query_expected = query_end.saturating_sub(query_start);
    let target_expected = target_end.saturating_sub(target_start);
    if query_consumed != query_expected || target_consumed != target_expected {
        return Err(ExactBlockRefinementError::StitchConsumption {
            query_consumed,
            query_expected,
            target_consumed,
            target_expected,
        });
    }

    let (matches, substitutions, insertions, deletions) = run_counts(&operations);
    let score = score_runs(&operations, options.scoring);
    let cigar = cigar_runs(&operations)?;
    let query_interval = BaseInterval::new(query_start as u64, query_end as u64)
        .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
    let target_interval = match strand {
        Strand::Forward => BaseInterval::new(target_start as u64, target_end as u64),
        Strand::Reverse => BaseInterval::new(
            (target_forward.len() - target_end) as u64,
            (target_forward.len() - target_start) as u64,
        ),
    }
    .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
    metadata.selected_width = metadata.selected_width.max(options.scoring.band_width);
    Ok(AlignmentResult {
        score,
        strand,
        query_interval,
        query_segments: vec![query_interval],
        target_interval,
        origin_crossing: false,
        query_length: query_end.saturating_sub(query_start) as u64,
        target_length: target_end.saturating_sub(target_start) as u64,
        matches,
        substitutions,
        insertions,
        deletions,
        cigar,
        edit_script: operations,
        chain_score: 0,
        retry_metadata: metadata,
    })
}

fn validate_exact_block_path(
    query_length: usize,
    target_length: usize,
    blocks: &[ExactBlock],
) -> Result<(), ExactBlockRefinementError> {
    let first = blocks.first().ok_or(ExactBlockRefinementError::NoBlocks)?;
    for block in blocks {
        if block.contig_id != first.contig_id {
            return Err(ExactBlockRefinementError::ContigMismatch {
                expected: first.contig_id,
                observed: block.contig_id,
            });
        }
        if block.strand != first.strand {
            return Err(ExactBlockRefinementError::StrandMismatch);
        }
        let query_end = usize::try_from(block.query_interval.end)
            .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
        let target_end = usize::try_from(block.oriented_target_interval.end)
            .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
        if query_end > query_length
            || target_end > target_length
            || block.query_interval.start > block.query_interval.end
            || block.oriented_target_interval.start > block.oriented_target_interval.end
            || block.query_interval.len() != block.oriented_target_interval.len()
        {
            return Err(ExactBlockRefinementError::BlockOutOfBounds);
        }
        let mapped_target = match block.strand {
            Strand::Forward => block.oriented_target_interval,
            Strand::Reverse => BaseInterval::new(
                target_length as u64 - block.oriented_target_interval.end,
                target_length as u64 - block.oriented_target_interval.start,
            )
            .map_err(|_| ExactBlockRefinementError::BlockOutOfBounds)?,
        };
        if mapped_target != block.target_interval {
            return Err(ExactBlockRefinementError::BlockTargetMismatch);
        }
    }
    for pair in blocks.windows(2) {
        if pair[0].query_interval.start > pair[1].query_interval.start
            || pair[0].oriented_target_interval.start > pair[1].oriented_target_interval.start
        {
            return Err(ExactBlockRefinementError::BlocksNotOrdered);
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn append_bounded_segment(
    query: &[u8],
    oriented_target: &[u8],
    query_start: usize,
    query_end: usize,
    target_start: usize,
    target_end: usize,
    workspace: &mut AlignmentWorkspace,
    options: AlignmentOptions,
    max_segment_bases: u64,
    operations: &mut Vec<EditRun>,
    metadata: &mut super::AlignmentRetryMetadata,
) -> Result<(), ExactBlockRefinementError> {
    if query_start > query_end || target_start > target_end {
        return Err(ExactBlockRefinementError::GapOutOfBounds);
    }
    let query_len = query_end - query_start;
    let target_len = target_end - target_start;
    if u64::try_from(query_len)
        .ok()
        .is_some_and(|length| length > max_segment_bases)
        || u64::try_from(target_len)
            .ok()
            .is_some_and(|length| length > max_segment_bases)
    {
        return Err(ExactBlockRefinementError::GapTooLarge {
            query_bases: query_len as u64,
            target_bases: target_len as u64,
        });
    }
    if query_len == 0 {
        append_run(operations, EditOperation::Insertion, target_len)?;
        return Ok(());
    }
    if target_len == 0 {
        append_run(operations, EditOperation::Deletion, query_len)?;
        return Ok(());
    }

    let result = workspace.align(
        &query[query_start..query_end],
        &oriented_target[target_start..target_end],
        options,
    )?;
    metadata.band_edge_touched |= result.retry_metadata.band_edge_touched;
    metadata.memory_bound_reached |= result.retry_metadata.memory_bound_reached;
    metadata.semiglobal_attempted |= result.retry_metadata.semiglobal_attempted;
    metadata.semiglobal_accepted |= result.retry_metadata.semiglobal_accepted;
    metadata
        .attempted_widths
        .extend(result.retry_metadata.attempted_widths);
    metadata.retries = metadata
        .retries
        .saturating_add(result.retry_metadata.retries);
    metadata.selected_width = metadata
        .selected_width
        .max(result.retry_metadata.selected_width);

    let query_core_start = usize::try_from(result.query_interval.start)
        .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
    let query_core_end = usize::try_from(result.query_interval.end)
        .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
    let target_core_start = usize::try_from(result.target_interval.start)
        .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
    let target_core_end = usize::try_from(result.target_interval.end)
        .map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?;
    append_run(operations, EditOperation::Deletion, query_core_start)?;
    append_run(operations, EditOperation::Insertion, target_core_start)?;
    for run in result.edit_script {
        append_run(
            operations,
            run.operation,
            usize::try_from(run.length).unwrap_or(usize::MAX),
        )?;
    }
    append_run(
        operations,
        EditOperation::Deletion,
        query_len.saturating_sub(query_core_end),
    )?;
    append_run(
        operations,
        EditOperation::Insertion,
        target_len.saturating_sub(target_core_end),
    )?;
    Ok(())
}

fn append_exact_run(
    operations: &mut Vec<EditRun>,
    length: u64,
) -> Result<(), ExactBlockRefinementError> {
    append_run(
        operations,
        EditOperation::Equal,
        usize::try_from(length).map_err(|_| ExactBlockRefinementError::CoordinateOverflow)?,
    )
}

fn append_run(
    operations: &mut Vec<EditRun>,
    operation: EditOperation,
    length: usize,
) -> Result<(), ExactBlockRefinementError> {
    if length == 0 {
        return Ok(());
    }
    let length = u32::try_from(length).map_err(|_| ExactBlockRefinementError::RunTooLong)?;
    if let Some(last) = operations.last_mut()
        && last.operation == operation
    {
        last.length = last
            .length
            .checked_add(length)
            .ok_or(ExactBlockRefinementError::RunTooLong)?;
    } else {
        operations.push(EditRun { operation, length });
    }
    Ok(())
}

fn run_consumption(operations: &[EditRun]) -> (usize, usize) {
    operations
        .iter()
        .fold((0usize, 0usize), |(query, target), run| {
            let length = usize::try_from(run.length).unwrap_or(usize::MAX);
            match run.operation {
                EditOperation::Equal | EditOperation::Substitution => {
                    (query.saturating_add(length), target.saturating_add(length))
                }
                EditOperation::Insertion | EditOperation::SoftClip => {
                    (query, target.saturating_add(length))
                }
                EditOperation::Deletion => (query.saturating_add(length), target),
            }
        })
}

fn run_counts(operations: &[EditRun]) -> (u64, u64, u64, u64) {
    let mut matches = 0u64;
    let mut substitutions = 0u64;
    let mut insertions = 0u64;
    let mut deletions = 0u64;
    for run in operations {
        match run.operation {
            EditOperation::Equal => matches = matches.saturating_add(u64::from(run.length)),
            EditOperation::Substitution => {
                substitutions = substitutions.saturating_add(u64::from(run.length))
            }
            EditOperation::Insertion => {
                insertions = insertions.saturating_add(u64::from(run.length))
            }
            EditOperation::Deletion => deletions = deletions.saturating_add(u64::from(run.length)),
            EditOperation::SoftClip => {}
        }
    }
    (matches, substitutions, insertions, deletions)
}

fn score_runs(operations: &[EditRun], scoring: AlignmentScoring) -> i32 {
    operations.iter().fold(0i32, |score, run| {
        let length = i32::try_from(run.length).unwrap_or(i32::MAX);
        let delta = match run.operation {
            EditOperation::Equal => scoring.match_score.saturating_mul(length),
            EditOperation::Substitution => scoring.mismatch_score.saturating_mul(length),
            EditOperation::Insertion | EditOperation::Deletion => scoring
                .gap_open_score
                .saturating_add(scoring.gap_extend_score.saturating_mul(length)),
            EditOperation::SoftClip => 0,
        };
        score.saturating_add(delta)
    })
}

fn cigar_runs(operations: &[EditRun]) -> Result<String, ExactBlockRefinementError> {
    let mut cigar = String::new();
    for run in operations {
        let code = match run.operation {
            EditOperation::Equal => '=',
            EditOperation::Substitution => 'X',
            EditOperation::Insertion => 'I',
            EditOperation::Deletion => 'D',
            EditOperation::SoftClip => 'S',
        };
        write!(&mut cigar, "{}{}", run.length, code)
            .map_err(|_| ExactBlockRefinementError::CigarWrite)?;
    }
    Ok(cigar)
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

/// Errors that tell the runner to use the existing corridor fallback.
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum ExactBlockRefinementError {
    #[error(transparent)]
    Alignment(#[from] AlignmentError),
    #[error(transparent)]
    Exact(#[from] ExactBlockError),
    #[error("exact-block refinement requires at least one block")]
    NoBlocks,
    #[error("exact blocks use multiple contigs: expected {expected}, observed {observed}")]
    ContigMismatch { expected: u32, observed: u32 },
    #[error("exact blocks use multiple strands")]
    StrandMismatch,
    #[error("exact block is outside the supplied query or target")]
    BlockOutOfBounds,
    #[error("exact block target interval does not match its oriented interval")]
    BlockTargetMismatch,
    #[error("exact blocks are not ordered in query and oriented-target coordinates")]
    BlocksNotOrdered,
    #[error("exact-block gap coordinates are invalid")]
    GapOutOfBounds,
    #[error(
        "exact-block gap exceeds the configured bound: query={query_bases}, target={target_bases}"
    )]
    GapTooLarge { query_bases: u64, target_bases: u64 },
    #[error(
        "exact-block stitch consumed query={query_consumed}/{query_expected}, target={target_consumed}/{target_expected}"
    )]
    StitchConsumption {
        query_consumed: usize,
        query_expected: usize,
        target_consumed: usize,
        target_expected: usize,
    },
    #[error("exact-block coordinate overflow")]
    CoordinateOverflow,
    #[error("exact-block edit run exceeds the serialized u32 limit")]
    RunTooLong,
    #[error("failed to write stitched CIGAR")]
    CigarWrite,
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
