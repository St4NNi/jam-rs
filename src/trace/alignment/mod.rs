//! Banded affine-gap alignment for trace candidates.
//!
//! The trace workflow aligns a query sequence against a bounded sequence window
//! from an assembly. This module deliberately exposes alignment evidence
//! (intervals, operations, and counters) rather than a presence call. The
//! dynamic-programming matrix is banded and owned by [`AlignmentWorkspace`] so
//! callers can reuse its allocations for a stream of candidate windows.

use crate::trace::config::AlignmentScoring;
use crate::trace::model::{
    AlignmentRole, BaseAlignment, BaseInterval, EditOperation, EditRun, SeedEvidence, Strand,
};
use serde::{Deserialize, Serialize};
use std::fmt::Write as _;
use thiserror::Error;

pub mod exact_blocks;
pub mod refine;

const STATE_MATCH: u8 = 0;
const STATE_INSERTION: u8 = 1;
const STATE_DELETION: u8 = 2;
const STATE_START: u8 = 3;
const STATE_UNREACHABLE: u8 = 4;

/// Default upper bound on the number of band cells allocated by an alignment.
/// Callers handling larger windows should set an explicit bound and split the
/// windows before aligning them.
pub const DEFAULT_MAX_CELLS: usize = 4_000_000;

/// Bounded band widths used by [`AlignmentWorkspace::align_with_retries`].
///
/// Keeping this list in the alignment module makes retry behavior stable
/// across callers and keeps a pathological candidate from widening without a
/// hard upper bound.
pub const DEFAULT_RETRY_WIDTHS: [u32; 6] = [64, 128, 256, 512, 1024, 2048];

/// One oriented anchor used to construct a piecewise alignment corridor.
///
/// Coordinates are relative to the sequences passed to the alignment call.
/// For a reverse-strand alignment, callers should provide target coordinates
/// in the reverse-complemented target window. The query and target spans are
/// used both for interpolation and for estimating diagonal drift between
/// consecutive anchors.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CorridorAnchor {
    pub query_position: u64,
    pub target_position: u64,
    pub span: u32,
}

impl CorridorAnchor {
    #[must_use]
    pub const fn new(query_position: u64, target_position: u64, span: u32) -> Self {
        Self {
            query_position,
            target_position,
            span,
        }
    }

    fn query_end(self) -> Option<u64> {
        self.query_position.checked_add(u64::from(self.span))
    }

    fn target_end(self) -> Option<u64> {
        self.target_position.checked_add(u64::from(self.span))
    }
}

/// A piecewise-linear target corridor following a monotone anchor path.
///
/// Between anchors the center follows the endpoint-to-endpoint slope. The
/// observed query/target gap discrepancy is added to the local half-width,
/// together with `safety_margin`, so an inferred insertion or deletion does
/// not force the DP through a constant diagonal. Widths supplied to the
/// retry API widen this corridor further while retaining the same centerline.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentCorridor {
    anchors: Vec<CorridorAnchor>,
    safety_margin: u32,
    max_half_width: u32,
}

impl AlignmentCorridor {
    /// Build a corridor from anchors ordered in both oriented coordinates.
    pub fn new(
        anchors: impl Into<Vec<CorridorAnchor>>,
        safety_margin: u32,
    ) -> Result<Self, AlignmentError> {
        let anchors = anchors.into();
        if anchors.is_empty() {
            return Err(AlignmentError::InvalidCorridor("no anchors"));
        }
        for anchor in &anchors {
            if anchor.span == 0 {
                return Err(AlignmentError::InvalidCorridor("zero anchor span"));
            }
            anchor.query_end().ok_or(AlignmentError::LengthOverflow)?;
            anchor.target_end().ok_or(AlignmentError::LengthOverflow)?;
        }
        for pair in anchors.windows(2) {
            let left = pair[0];
            let right = pair[1];
            if right.query_position < left.query_position
                || right.target_position < left.target_position
            {
                return Err(AlignmentError::InvalidCorridor("anchors are not monotone"));
            }
            // Overlapping k-mers are valid. Their starts still define the
            // interpolation order, while only positive gaps contribute drift.
        }
        Ok(Self {
            anchors,
            safety_margin,
            max_half_width: 2048,
        })
    }

    /// Set a hard half-width cap for callers with a stricter memory budget.
    #[must_use]
    pub const fn with_max_half_width(mut self, max_half_width: u32) -> Self {
        self.max_half_width = max_half_width;
        self
    }

    #[must_use]
    pub fn anchors(&self) -> &[CorridorAnchor] {
        &self.anchors
    }

    #[must_use]
    pub const fn safety_margin(&self) -> u32 {
        self.safety_margin
    }

    /// Maximum width required by the observed gaps before a retry margin.
    #[must_use]
    pub fn required_half_width(&self) -> u32 {
        self.anchors
            .windows(2)
            .map(|pair| self.local_half_width(pair[0], pair[1]))
            .max()
            .unwrap_or(self.safety_margin)
    }

    #[must_use]
    fn query_span(&self) -> u64 {
        self.anchors
            .first()
            .and_then(|first| {
                self.anchors
                    .last()
                    .and_then(|last| last.query_end()?.checked_sub(first.query_position))
            })
            .unwrap_or(0)
    }

    #[must_use]
    fn target_span(&self) -> u64 {
        self.anchors
            .first()
            .and_then(|first| {
                self.anchors
                    .last()
                    .and_then(|last| last.target_end()?.checked_sub(first.target_position))
            })
            .unwrap_or(0)
    }

    fn local_half_width(&self, left: CorridorAnchor, right: CorridorAnchor) -> u32 {
        let query_gap = right
            .query_position
            .saturating_sub(left.query_end().unwrap_or(u64::MAX));
        let target_gap = right
            .target_position
            .saturating_sub(left.target_end().unwrap_or(u64::MAX));
        let discrepancy = query_gap.abs_diff(target_gap);
        self.safety_margin.saturating_add(
            u32::try_from(discrepancy)
                .unwrap_or(u32::MAX)
                .min(self.max_half_width),
        )
    }

    fn local_extra(&self, query_index: u64) -> u32 {
        for pair in self.anchors.windows(2) {
            let left = pair[0];
            let right = pair[1];
            let left_end = left.query_end().unwrap_or(u64::MAX);
            if query_index >= left_end && query_index <= right.query_position {
                return self.local_half_width(left, right);
            }
        }
        self.safety_margin
    }

    fn center(&self, query_index: u64) -> Result<i128, AlignmentError> {
        let first = *self
            .anchors
            .first()
            .ok_or(AlignmentError::InvalidCorridor("no anchors"))?;
        let first_end = first.query_end().ok_or(AlignmentError::LengthOverflow)?;
        if query_index <= first.query_position {
            return signed_add(first.target_position, query_index, first.query_position);
        }
        if query_index <= first_end {
            return signed_add(first.target_position, query_index, first.query_position);
        }

        for pair in self.anchors.windows(2) {
            let left = pair[0];
            let right = pair[1];
            let left_query_end = left.query_end().ok_or(AlignmentError::LengthOverflow)?;
            if query_index >= left.query_position && query_index <= left_query_end {
                return signed_add(left.target_position, query_index, left.query_position);
            }
            let right_query = right.query_position;
            if query_index > right_query {
                continue;
            }
            let left_target_end = left.target_end().ok_or(AlignmentError::LengthOverflow)?;
            let target_delta = i128::from(right.target_position)
                .checked_sub(i128::from(left_target_end))
                .ok_or(AlignmentError::LengthOverflow)?;
            let query_delta = right_query
                .checked_sub(left_query_end)
                .ok_or(AlignmentError::LengthOverflow)?;
            if query_delta == 0 {
                return Ok(i128::from(left_target_end));
            }
            let query_offset = query_index
                .checked_sub(left_query_end)
                .ok_or(AlignmentError::LengthOverflow)?;
            let scaled = target_delta
                .checked_mul(i128::from(query_offset))
                .ok_or(AlignmentError::LengthOverflow)?
                / i128::from(query_delta);
            return i128::from(left_target_end)
                .checked_add(scaled)
                .ok_or(AlignmentError::LengthOverflow);
        }

        let last = *self.anchors.last().expect("validated non-empty corridor");
        let last_query_end = last.query_end().ok_or(AlignmentError::LengthOverflow)?;
        let last_target_end = last.target_end().ok_or(AlignmentError::LengthOverflow)?;
        signed_add(last_target_end, query_index, last_query_end)
    }

    fn bounds(
        &self,
        query_index: usize,
        target_len: usize,
        width: u32,
        diagonal_offset: isize,
    ) -> Result<(usize, usize), AlignmentError> {
        let query_index = u64::try_from(query_index).map_err(|_| AlignmentError::LengthOverflow)?;
        let center = self
            .center(query_index)?
            .checked_add(
                i128::try_from(diagonal_offset).map_err(|_| AlignmentError::LengthOverflow)?,
            )
            .ok_or(AlignmentError::LengthOverflow)?;
        let local_width = self.local_extra(query_index);
        let half_width = i128::from(width.max(local_width).min(self.max_half_width));
        band_bounds(center, half_width, target_len)
    }

    fn explains(&self, result: &AlignmentResult) -> bool {
        let query_span = self.query_span();
        let target_span = self.target_span();
        let margin = u64::from(self.safety_margin).saturating_mul(2);
        result.query_interval.len().saturating_add(margin) >= query_span
            && result.target_interval.len().saturating_add(margin) >= target_span
    }
}

fn signed_add(base: u64, value: u64, subtract: u64) -> Result<i128, AlignmentError> {
    i128::from(base)
        .checked_add(
            i128::from(value)
                .checked_sub(i128::from(subtract))
                .ok_or(AlignmentError::LengthOverflow)?,
        )
        .ok_or(AlignmentError::LengthOverflow)
}

fn band_bounds(
    center: i128,
    half_width: i128,
    target_len: usize,
) -> Result<(usize, usize), AlignmentError> {
    let target_len = i128::try_from(target_len).map_err(|_| AlignmentError::LengthOverflow)?;
    let lower = center
        .checked_sub(half_width)
        .ok_or(AlignmentError::LengthOverflow)?;
    let upper = center
        .checked_add(half_width)
        .ok_or(AlignmentError::LengthOverflow)?;
    if upper < 0 || lower > target_len {
        return Ok((1, 0));
    }
    Ok((
        usize::try_from(lower.max(0)).map_err(|_| AlignmentError::LengthOverflow)?,
        usize::try_from(upper.min(target_len)).map_err(|_| AlignmentError::LengthOverflow)?,
    ))
}

/// Details of bounded band retries and anchored refinement.
#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct AlignmentRetryMetadata {
    pub attempted_widths: Vec<u32>,
    pub selected_width: u32,
    pub retries: u8,
    pub band_edge_touched: bool,
    pub predicted_drift: bool,
    pub chain_span_under_explained: bool,
    pub retry_cap_reached: bool,
    pub memory_bound_reached: bool,
    pub corridor_used: bool,
    pub semiglobal_attempted: bool,
    pub semiglobal_accepted: bool,
}

/// Options that affect one local alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AlignmentOptions {
    /// Affine scoring and the half-band width.
    pub scoring: AlignmentScoring,
    /// The band is centered on `target_index - query_index == diagonal_offset`.
    /// This is useful when an anchor has already supplied an approximate offset.
    pub diagonal_offset: isize,
    /// Maximum number of cells (not bytes) retained in the banded matrix.
    pub max_cells: usize,
}

impl AlignmentOptions {
    #[must_use]
    pub const fn new(scoring: AlignmentScoring) -> Self {
        Self {
            scoring,
            diagonal_offset: 0,
            max_cells: DEFAULT_MAX_CELLS,
        }
    }

    #[must_use]
    pub const fn with_diagonal_offset(mut self, diagonal_offset: isize) -> Self {
        self.diagonal_offset = diagonal_offset;
        self
    }

    #[must_use]
    pub const fn with_max_cells(mut self, max_cells: usize) -> Self {
        self.max_cells = max_cells;
        self
    }
}

impl From<AlignmentScoring> for AlignmentOptions {
    fn from(scoring: AlignmentScoring) -> Self {
        Self::new(scoring)
    }
}

/// Compact reusable storage for the banded dynamic-programming matrix.
///
/// A workspace is not thread-safe by design. Give each concurrent alignment
/// worker its own workspace, then reuse it for that worker's candidate windows.
#[derive(Debug, Default)]
pub struct AlignmentWorkspace {
    cells: Vec<Cell>,
    semiglobal_cells: Vec<SemiCell>,
    row_offsets: Vec<usize>,
    row_starts: Vec<usize>,
    row_widths: Vec<usize>,
    operations: Vec<EditOperation>,
    query_buffer: Vec<u8>,
    target_buffer: Vec<u8>,
}

impl AlignmentWorkspace {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of matrix cells currently reserved by this workspace.
    #[must_use]
    pub fn capacity_cells(&self) -> usize {
        self.cells.capacity()
    }

    /// Align a query and an already oriented target window.
    pub fn align(
        &mut self,
        query: &[u8],
        target: &[u8],
        options: impl Into<AlignmentOptions>,
    ) -> Result<AlignmentResult, AlignmentError> {
        let options = options.into();
        let raw = self.align_raw(query, target, &options)?;
        raw.finish_linear(Strand::Forward, 0, target.len() as u64)
    }

    /// Align through a piecewise corridor without automatic widening.
    pub fn align_with_corridor(
        &mut self,
        query: &[u8],
        target: &[u8],
        corridor: &AlignmentCorridor,
        options: impl Into<AlignmentOptions>,
    ) -> Result<AlignmentResult, AlignmentError> {
        let options = options.into();
        let raw = self.align_raw_with_layout(query, target, &options, Some(corridor))?;
        let target_len = u64::try_from(target.len()).map_err(|_| AlignmentError::LengthOverflow)?;
        let mut result = raw.finish_linear(Strand::Forward, 0, target_len)?;
        result.retry_metadata.corridor_used = true;
        result.retry_metadata.predicted_drift =
            corridor.required_half_width() > options.scoring.band_width;
        result.retry_metadata.chain_span_under_explained = !corridor.explains(&result);
        Ok(result)
    }

    /// Corridor alignment for a forward-stored target, including reverse
    /// strand coordinate projection.
    pub fn align_with_corridor_oriented(
        &mut self,
        query: &[u8],
        target_forward: &[u8],
        target_start: u64,
        strand: Strand,
        corridor: &AlignmentCorridor,
        options: impl Into<AlignmentOptions>,
    ) -> Result<AlignmentResult, AlignmentError> {
        let options = options.into();
        let target_buffer = if strand == Strand::Reverse {
            let mut buffer = std::mem::take(&mut self.target_buffer);
            buffer.clear();
            buffer.reserve(target_forward.len());
            buffer.extend(target_forward.iter().rev().map(|base| complement(*base)));
            Some(buffer)
        } else {
            None
        };
        let oriented_target = target_buffer.as_deref().unwrap_or(target_forward);
        let result = self
            .align_with_corridor(query, oriented_target, corridor, options)
            .and_then(|result| {
                map_target_coordinates(result, target_forward.len(), target_start, strand)
            });
        if let Some(buffer) = target_buffer {
            self.target_buffer = buffer;
        }
        result
    }

    /// Align with bounded corridor widening. A successful alignment that still
    /// touches the final band edge is returned with `retry_cap_reached=true`;
    /// callers can decide whether their chain evidence is sufficient.
    pub fn align_with_retries(
        &mut self,
        query: &[u8],
        target: &[u8],
        corridor: &AlignmentCorridor,
        options: impl Into<AlignmentOptions>,
    ) -> Result<AlignmentResult, AlignmentError> {
        let options = options.into();
        let predicted_drift = corridor.required_half_width() > options.scoring.band_width;
        let mut last_error = None;
        let mut selected: Option<AlignmentResult> = None;
        let mut metadata = AlignmentRetryMetadata {
            predicted_drift,
            corridor_used: true,
            ..AlignmentRetryMetadata::default()
        };

        for &width in &DEFAULT_RETRY_WIDTHS {
            let width = width.min(corridor.max_half_width);
            if metadata.attempted_widths.last().copied() == Some(width) {
                continue;
            }
            metadata.attempted_widths.push(width);
            let mut scoring = options.scoring;
            scoring.band_width = width;
            let attempt_options = AlignmentOptions { scoring, ..options };
            let raw =
                match self.align_raw_with_layout(query, target, &attempt_options, Some(corridor)) {
                    Ok(raw) => raw,
                    Err(error @ AlignmentError::MatrixTooLarge { .. }) => {
                        metadata.memory_bound_reached = true;
                        last_error = Some(error);
                        break;
                    }
                    Err(error) => {
                        last_error = Some(error);
                        continue;
                    }
                };
            let target_len =
                u64::try_from(target.len()).map_err(|_| AlignmentError::LengthOverflow)?;
            let result = raw.finish_linear(Strand::Forward, 0, target_len)?;
            let edge_touched = result.retry_metadata.band_edge_touched;
            let under_explained = !corridor.explains(&result);
            let predicted = predicted_drift && width < corridor.required_half_width();
            metadata.band_edge_touched |= edge_touched;
            metadata.chain_span_under_explained |= under_explained;
            metadata.selected_width = width;
            if !(edge_touched || under_explained || predicted) {
                selected = Some(result);
                break;
            }
            selected = Some(result);
        }

        let Some(mut result) = selected else {
            return Err(last_error.unwrap_or(AlignmentError::NoAlignment));
        };
        metadata.retries =
            u8::try_from(metadata.attempted_widths.len().saturating_sub(1)).unwrap_or(u8::MAX);
        metadata.retry_cap_reached = metadata
            .attempted_widths
            .last()
            .is_some_and(|width| *width >= corridor.max_half_width || *width >= 2048)
            && (result.retry_metadata.band_edge_touched
                || !corridor.explains(&result)
                || (predicted_drift
                    && metadata.selected_width.lt(&corridor.required_half_width())));
        result.retry_metadata = metadata;
        Ok(result)
    }

    /// Select a corridor width with rolling score rows, then traceback once.
    pub fn align_with_score_only_retries(
        &mut self,
        query: &[u8],
        target: &[u8],
        corridor: &AlignmentCorridor,
        options: impl Into<AlignmentOptions>,
    ) -> Result<AlignmentResult, AlignmentError> {
        let options = options.into();
        let mut widths = DEFAULT_RETRY_WIDTHS
            .iter()
            .map(|width| (*width).min(corridor.max_half_width))
            .collect::<Vec<_>>();
        widths.dedup();
        let mut scratch = refine::ScoreOnlyWorkspace::new();
        refine::refine_with_score_only(
            &widths,
            |width| {
                let mut scoring = options.scoring;
                scoring.band_width = width;
                refine::score_only_corridor_with_scratch(
                    &mut scratch,
                    query,
                    target,
                    corridor,
                    AlignmentOptions { scoring, ..options },
                )
            },
            |width| {
                let mut scoring = options.scoring;
                scoring.band_width = width;
                self.align_with_corridor(
                    query,
                    target,
                    corridor,
                    AlignmentOptions { scoring, ..options },
                )
            },
        )
        .map(|refined| refined.result)
    }

    /// Bounded corridor retries for a forward-stored target on either strand.
    pub fn align_with_retries_oriented(
        &mut self,
        query: &[u8],
        target_forward: &[u8],
        target_start: u64,
        strand: Strand,
        corridor: &AlignmentCorridor,
        options: impl Into<AlignmentOptions>,
    ) -> Result<AlignmentResult, AlignmentError> {
        let options = options.into();
        let target_buffer = if strand == Strand::Reverse {
            let mut buffer = std::mem::take(&mut self.target_buffer);
            buffer.clear();
            buffer.reserve(target_forward.len());
            buffer.extend(target_forward.iter().rev().map(|base| complement(*base)));
            Some(buffer)
        } else {
            None
        };
        let oriented_target = target_buffer.as_deref().unwrap_or(target_forward);
        let result = self
            .align_with_retries(query, oriented_target, corridor, options)
            .and_then(|result| {
                map_target_coordinates(result, target_forward.len(), target_start, strand)
            });
        if let Some(buffer) = target_buffer {
            self.target_buffer = buffer;
        }
        result
    }

    /// Refine a high-confidence local core with an anchored semiglobal pass.
    ///
    /// The semiglobal pass consumes the complete query while leaving target
    /// prefixes and suffixes free. If forced end-to-end consumption explains
    /// less score than the local core, the core is retained, which preserves
    /// true clipping at contig ends and at genuinely shorter fragments.
    pub fn align_anchored_semiglobal(
        &mut self,
        query: &[u8],
        target: &[u8],
        corridor: &AlignmentCorridor,
        options: impl Into<AlignmentOptions>,
    ) -> Result<AlignmentResult, AlignmentError> {
        let options = options.into();
        let mut core = self.align_with_retries(query, target, corridor, options)?;
        core.retry_metadata.semiglobal_attempted = true;
        let width = core
            .retry_metadata
            .selected_width
            .max(options.scoring.band_width);
        let mut scoring = options.scoring;
        scoring.band_width = width;
        let semiglobal_options = AlignmentOptions { scoring, ..options };
        let refined = self
            .align_semiglobal_raw(query, target, &semiglobal_options, Some(corridor))
            .and_then(|raw| {
                let target_len =
                    u64::try_from(target.len()).map_err(|_| AlignmentError::LengthOverflow)?;
                raw.finish_linear(Strand::Forward, 0, target_len)
            });
        match refined {
            Ok(mut refined) if refined.score >= core.score => {
                refined.retry_metadata = core.retry_metadata;
                refined.retry_metadata.semiglobal_attempted = true;
                refined.retry_metadata.semiglobal_accepted = true;
                Ok(refined)
            }
            Ok(_) | Err(_) => {
                core.retry_metadata.semiglobal_accepted = false;
                Ok(core)
            }
        }
    }

    /// Anchored semiglobal refinement for a forward-stored target on either
    /// strand.
    pub fn align_anchored_semiglobal_oriented(
        &mut self,
        query: &[u8],
        target_forward: &[u8],
        target_start: u64,
        strand: Strand,
        corridor: &AlignmentCorridor,
        options: impl Into<AlignmentOptions>,
    ) -> Result<AlignmentResult, AlignmentError> {
        let options = options.into();
        let target_buffer = if strand == Strand::Reverse {
            let mut buffer = std::mem::take(&mut self.target_buffer);
            buffer.clear();
            buffer.reserve(target_forward.len());
            buffer.extend(target_forward.iter().rev().map(|base| complement(*base)));
            Some(buffer)
        } else {
            None
        };
        let oriented_target = target_buffer.as_deref().unwrap_or(target_forward);
        let result = self
            .align_anchored_semiglobal(query, oriented_target, corridor, options)
            .and_then(|result| {
                map_target_coordinates(result, target_forward.len(), target_start, strand)
            });
        if let Some(buffer) = target_buffer {
            self.target_buffer = buffer;
        }
        result
    }

    /// Align against a forward-stored target window, reporting a forward target
    /// interval even when the selected strand is reverse.
    pub fn align_oriented(
        &mut self,
        query: &[u8],
        target_forward: &[u8],
        target_start: u64,
        strand: Strand,
        options: impl Into<AlignmentOptions>,
    ) -> Result<AlignmentResult, AlignmentError> {
        let options = options.into();
        let target_buffer = if strand == Strand::Reverse {
            let mut buffer = std::mem::take(&mut self.target_buffer);
            buffer.clear();
            buffer.reserve(target_forward.len());
            buffer.extend(target_forward.iter().rev().map(|base| complement(*base)));
            Some(buffer)
        } else {
            None
        };
        let oriented_target = target_buffer.as_deref().unwrap_or(target_forward);
        let result = self
            .align_raw(query, oriented_target, &options)
            .and_then(|raw| {
                raw.finish_oriented(
                    strand,
                    target_start,
                    target_forward.len() as u64,
                    strand == Strand::Reverse,
                )
            });
        if let Some(buffer) = target_buffer {
            self.target_buffer = buffer;
        }
        result
    }

    /// Align a linearized window from a circular plasmid query.
    ///
    /// `query_start` identifies the first base in the linearized query and
    /// `query_span` is limited to one plasmid length so the result can be
    /// represented by at most two non-wrapping query segments.
    #[allow(clippy::too_many_arguments)]
    pub fn align_circular(
        &mut self,
        query: &[u8],
        query_start: u64,
        query_span: u64,
        target_forward: &[u8],
        target_start: u64,
        strand: Strand,
        options: impl Into<AlignmentOptions>,
    ) -> Result<AlignmentResult, AlignmentError> {
        let options = options.into();
        let query_len = u64::try_from(query.len()).map_err(|_| AlignmentError::LengthOverflow)?;
        if query_len == 0 {
            return Err(AlignmentError::EmptyQuery);
        }
        if query_start >= query_len {
            return Err(AlignmentError::CircularStart {
                start: query_start,
                length: query_len,
            });
        }
        if query_span == 0 || query_span > query_len {
            return Err(AlignmentError::CircularSpan {
                span: query_span,
                length: query_len,
            });
        }

        let mut query_buffer = std::mem::take(&mut self.query_buffer);
        query_buffer.clear();
        let span = usize::try_from(query_span).map_err(|_| AlignmentError::LengthOverflow)?;
        query_buffer.reserve(span);
        for index in 0..span {
            let absolute = usize::try_from(query_start)
                .ok()
                .and_then(|start| start.checked_add(index))
                .ok_or(AlignmentError::LengthOverflow)?;
            query_buffer.push(query[absolute % query.len()]);
        }

        let mut target_buffer = if strand == Strand::Reverse {
            let mut buffer = std::mem::take(&mut self.target_buffer);
            buffer.clear();
            buffer.reserve(target_forward.len());
            buffer.extend(target_forward.iter().rev().map(|base| complement(*base)));
            Some(buffer)
        } else {
            None
        };
        let oriented_target = target_buffer.as_deref().unwrap_or(target_forward);
        let result = self
            .align_raw(&query_buffer, oriented_target, &options)
            .and_then(|raw| {
                raw.finish_oriented(
                    strand,
                    target_start,
                    target_forward.len() as u64,
                    strand == Strand::Reverse,
                )
            });
        self.query_buffer = query_buffer;
        if let Some(buffer) = target_buffer.take() {
            self.target_buffer = buffer;
        }
        let mut result = result?;

        let linear_start = result.query_interval.start;
        let linear_end = result.query_interval.end;
        let absolute_start = query_start
            .checked_add(linear_start)
            .ok_or(AlignmentError::LengthOverflow)?;
        let absolute_end = query_start
            .checked_add(linear_end)
            .ok_or(AlignmentError::LengthOverflow)?;
        result.query_segments = circular_segments(absolute_start, absolute_end, query_len)?;
        result.origin_crossing = result.query_segments.len() == 2;
        Ok(result)
    }

    /// Circular-query variant of [`Self::align_with_retries`]. The corridor
    /// coordinates are relative to the linearized query and oriented target
    /// window, matching the other corridor entry points.
    #[allow(clippy::too_many_arguments)]
    pub fn align_circular_with_retries(
        &mut self,
        query: &[u8],
        query_start: u64,
        query_span: u64,
        target_forward: &[u8],
        target_start: u64,
        strand: Strand,
        corridor: &AlignmentCorridor,
        options: impl Into<AlignmentOptions>,
    ) -> Result<AlignmentResult, AlignmentError> {
        let query_len = u64::try_from(query.len()).map_err(|_| AlignmentError::LengthOverflow)?;
        if query_len == 0 {
            return Err(AlignmentError::EmptyQuery);
        }
        if query_start >= query_len {
            return Err(AlignmentError::CircularStart {
                start: query_start,
                length: query_len,
            });
        }
        if query_span == 0 || query_span > query_len {
            return Err(AlignmentError::CircularSpan {
                span: query_span,
                length: query_len,
            });
        }
        let span = usize::try_from(query_span).map_err(|_| AlignmentError::LengthOverflow)?;
        let mut query_buffer = std::mem::take(&mut self.query_buffer);
        query_buffer.clear();
        query_buffer.reserve(span);
        for index in 0..span {
            let absolute = usize::try_from(query_start)
                .ok()
                .and_then(|start| start.checked_add(index))
                .ok_or(AlignmentError::LengthOverflow)?;
            query_buffer.push(query[absolute % query.len()]);
        }
        let mut target_buffer = if strand == Strand::Reverse {
            let mut buffer = std::mem::take(&mut self.target_buffer);
            buffer.clear();
            buffer.reserve(target_forward.len());
            buffer.extend(target_forward.iter().rev().map(|base| complement(*base)));
            Some(buffer)
        } else {
            None
        };
        let oriented_target = target_buffer.as_deref().unwrap_or(target_forward);
        let options = options.into();
        let result = self
            .align_with_retries(&query_buffer, oriented_target, corridor, options)
            .and_then(|result| {
                map_target_coordinates(result, target_forward.len(), target_start, strand)
            });
        self.query_buffer = query_buffer;
        if let Some(buffer) = target_buffer.take() {
            self.target_buffer = buffer;
        }
        let mut result = result?;
        let absolute_start = query_start
            .checked_add(result.query_interval.start)
            .ok_or(AlignmentError::LengthOverflow)?;
        let absolute_end = query_start
            .checked_add(result.query_interval.end)
            .ok_or(AlignmentError::LengthOverflow)?;
        result.query_segments = circular_segments(absolute_start, absolute_end, query_len)?;
        result.origin_crossing = result.query_segments.len() == 2;
        Ok(result)
    }

    /// Circular-query variant of [`Self::align_anchored_semiglobal`].
    #[allow(clippy::too_many_arguments)]
    pub fn align_circular_anchored_semiglobal(
        &mut self,
        query: &[u8],
        query_start: u64,
        query_span: u64,
        target_forward: &[u8],
        target_start: u64,
        strand: Strand,
        corridor: &AlignmentCorridor,
        options: impl Into<AlignmentOptions>,
    ) -> Result<AlignmentResult, AlignmentError> {
        let query_len = u64::try_from(query.len()).map_err(|_| AlignmentError::LengthOverflow)?;
        if query_len == 0 {
            return Err(AlignmentError::EmptyQuery);
        }
        if query_start >= query_len {
            return Err(AlignmentError::CircularStart {
                start: query_start,
                length: query_len,
            });
        }
        if query_span == 0 || query_span > query_len {
            return Err(AlignmentError::CircularSpan {
                span: query_span,
                length: query_len,
            });
        }
        let span = usize::try_from(query_span).map_err(|_| AlignmentError::LengthOverflow)?;
        let mut query_buffer = std::mem::take(&mut self.query_buffer);
        query_buffer.clear();
        query_buffer.reserve(span);
        for index in 0..span {
            let absolute = usize::try_from(query_start)
                .ok()
                .and_then(|start| start.checked_add(index))
                .ok_or(AlignmentError::LengthOverflow)?;
            query_buffer.push(query[absolute % query.len()]);
        }
        let mut target_buffer = if strand == Strand::Reverse {
            let mut buffer = std::mem::take(&mut self.target_buffer);
            buffer.clear();
            buffer.reserve(target_forward.len());
            buffer.extend(target_forward.iter().rev().map(|base| complement(*base)));
            Some(buffer)
        } else {
            None
        };
        let oriented_target = target_buffer.as_deref().unwrap_or(target_forward);
        let options = options.into();
        let result = self
            .align_anchored_semiglobal(&query_buffer, oriented_target, corridor, options)
            .and_then(|result| {
                map_target_coordinates(result, target_forward.len(), target_start, strand)
            });
        self.query_buffer = query_buffer;
        if let Some(buffer) = target_buffer.take() {
            self.target_buffer = buffer;
        }
        let mut result = result?;
        let absolute_start = query_start
            .checked_add(result.query_interval.start)
            .ok_or(AlignmentError::LengthOverflow)?;
        let absolute_end = query_start
            .checked_add(result.query_interval.end)
            .ok_or(AlignmentError::LengthOverflow)?;
        result.query_segments = circular_segments(absolute_start, absolute_end, query_len)?;
        result.origin_crossing = result.query_segments.len() == 2;
        Ok(result)
    }

    fn align_raw(
        &mut self,
        query: &[u8],
        target: &[u8],
        options: &AlignmentOptions,
    ) -> Result<RawAlignment, AlignmentError> {
        self.align_raw_with_layout(query, target, options, None)
    }

    fn align_raw_with_layout(
        &mut self,
        query: &[u8],
        target: &[u8],
        options: &AlignmentOptions,
        corridor: Option<&AlignmentCorridor>,
    ) -> Result<RawAlignment, AlignmentError> {
        if query.is_empty() {
            return Err(AlignmentError::EmptyQuery);
        }
        if target.is_empty() {
            return Err(AlignmentError::EmptyTarget);
        }
        let query_len = query.len();
        let target_len = target.len();

        self.row_offsets.clear();
        self.row_starts.clear();
        self.row_widths.clear();
        self.row_offsets.reserve(query_len.saturating_add(2));
        self.row_starts.reserve(query_len.saturating_add(1));
        self.row_widths.reserve(query_len.saturating_add(1));

        let mut total_cells = 0usize;
        for query_index in 0..=query_len {
            let (start, end) = if let Some(corridor) = corridor {
                corridor.bounds(
                    query_index,
                    target_len,
                    options.scoring.band_width,
                    options.diagonal_offset,
                )?
            } else {
                let center = i128::try_from(query_index)
                    .map_err(|_| AlignmentError::LengthOverflow)?
                    .checked_add(
                        i128::try_from(options.diagonal_offset)
                            .map_err(|_| AlignmentError::LengthOverflow)?,
                    )
                    .ok_or(AlignmentError::LengthOverflow)?;
                let band = i128::from(options.scoring.band_width);
                band_bounds(center, band, target_len)?
            };
            let (start, width) = if start > end {
                (0, 0)
            } else {
                (start, end - start + 1)
            };
            total_cells = total_cells
                .checked_add(width)
                .ok_or(AlignmentError::LengthOverflow)?;
            self.row_offsets.push(total_cells - width);
            self.row_starts.push(start);
            self.row_widths.push(width);
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

        let default_cell = Cell::default();
        if self.cells.len() < total_cells {
            self.cells.resize(total_cells, default_cell);
        } else {
            self.cells[..total_cells].fill(default_cell);
            self.cells.truncate(total_cells);
        }
        self.operations.clear();

        let mut best = BestCell::default();
        for query_index in 0..=query_len {
            let start = self.row_starts[query_index];
            let width = self.row_widths[query_index];
            for target_index in start..start + width {
                if query_index == 0 && target_index == 0 {
                    continue;
                }

                let cell_index = self.cell_index(query_index, target_index);
                let mut cell = Cell::default();

                if query_index > 0
                    && target_index > 0
                    && let Some(previous_index) =
                        self.cell_index_checked(query_index - 1, target_index - 1)
                {
                    let previous = self.cells[previous_index];
                    let (previous_score, previous_state) = previous.best_score();
                    let substitution =
                        if equal_base(query[query_index - 1], target[target_index - 1]) {
                            options.scoring.match_score
                        } else {
                            options.scoring.mismatch_score
                        };
                    let score = previous_score.saturating_add(substitution);
                    if score > 0 {
                        cell.scores[usize::from(STATE_MATCH)] = score;
                        cell.previous[usize::from(STATE_MATCH)] = if previous_score == 0 {
                            STATE_START
                        } else {
                            previous_state
                        };
                    }
                }

                if target_index > 0
                    && let Some(previous_index) =
                        self.cell_index_checked(query_index, target_index - 1)
                {
                    let previous = self.cells[previous_index];
                    let extension = previous.scores[usize::from(STATE_INSERTION)]
                        .saturating_add(options.scoring.gap_extend_score);
                    let opened_from_match = previous.scores[usize::from(STATE_MATCH)]
                        .saturating_add(gap_open_score(options.scoring));
                    let opened_from_deletion = previous.scores[usize::from(STATE_DELETION)]
                        .saturating_add(gap_open_score(options.scoring));
                    let (score, previous_state) = choose_transition([
                        (extension, STATE_INSERTION),
                        (opened_from_match, STATE_MATCH),
                        (opened_from_deletion, STATE_DELETION),
                    ]);
                    if score > 0 {
                        cell.scores[usize::from(STATE_INSERTION)] = score;
                        cell.previous[usize::from(STATE_INSERTION)] = previous_state;
                    }
                }

                if query_index > 0
                    && let Some(previous_index) =
                        self.cell_index_checked(query_index - 1, target_index)
                {
                    let previous = self.cells[previous_index];
                    let extension = previous.scores[usize::from(STATE_DELETION)]
                        .saturating_add(options.scoring.gap_extend_score);
                    let opened_from_match = previous.scores[usize::from(STATE_MATCH)]
                        .saturating_add(gap_open_score(options.scoring));
                    let opened_from_insertion = previous.scores[usize::from(STATE_INSERTION)]
                        .saturating_add(gap_open_score(options.scoring));
                    let (score, previous_state) = choose_transition([
                        (extension, STATE_DELETION),
                        (opened_from_match, STATE_MATCH),
                        (opened_from_insertion, STATE_INSERTION),
                    ]);
                    if score > 0 {
                        cell.scores[usize::from(STATE_DELETION)] = score;
                        cell.previous[usize::from(STATE_DELETION)] = previous_state;
                    }
                }

                self.cells[cell_index] = cell;
                best.consider(query_index, target_index, cell);
            }
        }

        if best.score <= 0 {
            return Err(AlignmentError::NoAlignment);
        }

        let (query_start, target_start, _) = self.traceback(query, target, best)?;
        let query_interval = BaseInterval::new(query_start as u64, best.query_index as u64)
            .map_err(|_| AlignmentError::LengthOverflow)?;
        let target_interval = BaseInterval::new(target_start as u64, best.target_index as u64)
            .map_err(|_| AlignmentError::LengthOverflow)?;
        let operations = self.operations.clone();
        let (edit_script, cigar, matches, substitutions, insertions, deletions) =
            summarize_operations(&operations)?;
        let band_edge_touched =
            self.path_touches_row_edge(query_start, target_start, target.len(), &operations);
        Ok(RawAlignment {
            score: best.score,
            query_interval,
            target_interval,
            matches,
            substitutions,
            insertions,
            deletions,
            cigar,
            edit_script,
            band_edge_touched,
            band_width: options.scoring.band_width,
        })
    }

    fn path_touches_row_edge(
        &self,
        mut query_index: usize,
        mut target_index: usize,
        target_len: usize,
        operations: &[EditOperation],
    ) -> bool {
        let mut touched = false;
        for &operation in operations {
            if let Some(start) = self.row_starts.get(query_index)
                && let Some(width) = self.row_widths.get(query_index)
                && *width > 0
                && ((target_index == *start && *start > 0)
                    || (target_index == start.saturating_add(width.saturating_sub(1))
                        && target_index < target_len))
            {
                touched = true;
            }
            match operation {
                EditOperation::Equal | EditOperation::Substitution => {
                    query_index = query_index.saturating_add(1);
                    target_index = target_index.saturating_add(1);
                }
                EditOperation::Insertion => {
                    target_index = target_index.saturating_add(1);
                }
                EditOperation::Deletion => {
                    query_index = query_index.saturating_add(1);
                }
                EditOperation::SoftClip => {}
            }
        }
        if let Some(start) = self.row_starts.get(query_index)
            && let Some(width) = self.row_widths.get(query_index)
            && *width > 0
            && ((target_index == *start && *start > 0)
                || (target_index == start.saturating_add(width.saturating_sub(1))
                    && target_index < target_len))
        {
            touched = true;
        }
        touched
    }

    fn align_semiglobal_raw(
        &mut self,
        query: &[u8],
        target: &[u8],
        options: &AlignmentOptions,
        corridor: Option<&AlignmentCorridor>,
    ) -> Result<RawAlignment, AlignmentError> {
        if query.is_empty() {
            return Err(AlignmentError::EmptyQuery);
        }
        if target.is_empty() {
            return Err(AlignmentError::EmptyTarget);
        }
        let query_len = query.len();
        let target_len = target.len();
        self.row_offsets.clear();
        self.row_starts.clear();
        self.row_widths.clear();
        self.row_offsets.reserve(query_len.saturating_add(2));
        self.row_starts.reserve(query_len.saturating_add(1));
        self.row_widths.reserve(query_len.saturating_add(1));

        let mut total_cells = 0usize;
        for query_index in 0..=query_len {
            let (start, end) = if let Some(corridor) = corridor {
                corridor.bounds(
                    query_index,
                    target_len,
                    options.scoring.band_width,
                    options.diagonal_offset,
                )?
            } else {
                let center = i128::try_from(query_index)
                    .map_err(|_| AlignmentError::LengthOverflow)?
                    .checked_add(
                        i128::try_from(options.diagonal_offset)
                            .map_err(|_| AlignmentError::LengthOverflow)?,
                    )
                    .ok_or(AlignmentError::LengthOverflow)?;
                let band = i128::from(options.scoring.band_width);
                band_bounds(center, band, target_len)?
            };
            let (start, width) = if start > end {
                (0, 0)
            } else {
                (start, end - start + 1)
            };
            total_cells = total_cells
                .checked_add(width)
                .ok_or(AlignmentError::LengthOverflow)?;
            self.row_offsets.push(total_cells - width);
            self.row_starts.push(start);
            self.row_widths.push(width);
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
        let default_cell = SemiCell::default();
        if self.semiglobal_cells.len() < total_cells {
            self.semiglobal_cells.resize(total_cells, default_cell);
        } else {
            self.semiglobal_cells[..total_cells].fill(default_cell);
            self.semiglobal_cells.truncate(total_cells);
        }

        for query_index in 0..=query_len {
            let start = self.row_starts[query_index];
            let width = self.row_widths[query_index];
            for target_index in start..start + width {
                let cell_index = self.cell_index(query_index, target_index);
                if query_index == 0 {
                    // Any target prefix may be skipped at no cost. The match
                    // state is used as the zero-cost start marker.
                    self.semiglobal_cells[cell_index].scores[usize::from(STATE_MATCH)] = 0;
                    self.semiglobal_cells[cell_index].previous[usize::from(STATE_MATCH)] =
                        STATE_START;
                    continue;
                }
                let mut cell = SemiCell::default();
                if target_index > 0
                    && let Some(previous_index) =
                        self.cell_index_checked(query_index - 1, target_index - 1)
                {
                    let previous = self.semiglobal_cells[previous_index];
                    let (previous_score, previous_state) = previous.best_score();
                    if previous_score > SEMI_NEG_INF {
                        let substitution =
                            if equal_base(query[query_index - 1], target[target_index - 1]) {
                                options.scoring.match_score
                            } else {
                                options.scoring.mismatch_score
                            };
                        cell.scores[usize::from(STATE_MATCH)] =
                            previous_score.saturating_add(substitution);
                        cell.previous[usize::from(STATE_MATCH)] = previous_state;
                    }
                }
                if target_index > 0
                    && let Some(previous_index) =
                        self.cell_index_checked(query_index, target_index - 1)
                {
                    let previous = self.semiglobal_cells[previous_index];
                    let extension = semi_add(
                        previous.scores[usize::from(STATE_INSERTION)],
                        options.scoring.gap_extend_score,
                    );
                    let opened_from_match = semi_add(
                        previous.scores[usize::from(STATE_MATCH)],
                        gap_open_score(options.scoring),
                    );
                    let opened_from_deletion = semi_add(
                        previous.scores[usize::from(STATE_DELETION)],
                        gap_open_score(options.scoring),
                    );
                    let (score, previous_state) = choose_semi_transition([
                        (extension, STATE_INSERTION),
                        (opened_from_match, STATE_MATCH),
                        (opened_from_deletion, STATE_DELETION),
                    ]);
                    cell.scores[usize::from(STATE_INSERTION)] = score;
                    cell.previous[usize::from(STATE_INSERTION)] = previous_state;
                }
                if let Some(previous_index) = self.cell_index_checked(query_index - 1, target_index)
                {
                    let previous = self.semiglobal_cells[previous_index];
                    let extension = semi_add(
                        previous.scores[usize::from(STATE_DELETION)],
                        options.scoring.gap_extend_score,
                    );
                    let opened_from_match = semi_add(
                        previous.scores[usize::from(STATE_MATCH)],
                        gap_open_score(options.scoring),
                    );
                    let opened_from_insertion = semi_add(
                        previous.scores[usize::from(STATE_INSERTION)],
                        gap_open_score(options.scoring),
                    );
                    let (score, previous_state) = choose_semi_transition([
                        (extension, STATE_DELETION),
                        (opened_from_match, STATE_MATCH),
                        (opened_from_insertion, STATE_INSERTION),
                    ]);
                    cell.scores[usize::from(STATE_DELETION)] = score;
                    cell.previous[usize::from(STATE_DELETION)] = previous_state;
                }
                self.semiglobal_cells[cell_index] = cell;
            }
        }

        let final_start = self.row_starts[query_len];
        let final_width = self.row_widths[query_len];
        let mut best = (SEMI_NEG_INF, final_start, STATE_UNREACHABLE);
        for target_index in final_start..final_start + final_width {
            let cell = self.semiglobal_cells[self.cell_index(query_len, target_index)];
            let (score, state) = cell.best_score();
            if score > best.0
                || (score == best.0
                    && (target_index < best.1 || (target_index == best.1 && state < best.2)))
            {
                best = (score, target_index, state);
            }
        }
        if best.2 >= STATE_START || best.0 <= SEMI_NEG_INF {
            return Err(AlignmentError::NoAlignment);
        }

        let (target_start, operations) =
            self.traceback_semiglobal(query, target, query_len, best.1, best.2)?;
        let (edit_script, cigar, matches, substitutions, insertions, deletions) =
            summarize_operations(&operations)?;
        let band_edge_touched =
            self.path_touches_row_edge(0, target_start, target.len(), &operations);
        Ok(RawAlignment {
            score: best.0,
            query_interval: BaseInterval::new(
                0,
                u64::try_from(query_len).map_err(|_| AlignmentError::LengthOverflow)?,
            )
            .map_err(|_| AlignmentError::LengthOverflow)?,
            target_interval: BaseInterval::new(
                u64::try_from(target_start).map_err(|_| AlignmentError::LengthOverflow)?,
                u64::try_from(best.1).map_err(|_| AlignmentError::LengthOverflow)?,
            )
            .map_err(|_| AlignmentError::LengthOverflow)?,
            matches,
            substitutions,
            insertions,
            deletions,
            cigar,
            edit_script,
            band_edge_touched,
            band_width: options.scoring.band_width,
        })
    }

    fn traceback_semiglobal(
        &mut self,
        query: &[u8],
        target: &[u8],
        mut query_index: usize,
        mut target_index: usize,
        mut state: u8,
    ) -> Result<(usize, Vec<EditOperation>), AlignmentError> {
        let mut operations = Vec::new();
        while query_index > 0 {
            let index = self
                .cell_index_checked(query_index, target_index)
                .ok_or(AlignmentError::TracebackOutsideBand)?;
            let cell = self.semiglobal_cells[index];
            if state >= STATE_START || cell.scores[usize::from(state)] <= SEMI_NEG_INF {
                return Err(AlignmentError::InvalidTraceback);
            }
            match state {
                STATE_MATCH => {
                    if target_index == 0 {
                        return Err(AlignmentError::InvalidTraceback);
                    }
                    operations.push(
                        if equal_base(query[query_index - 1], target[target_index - 1]) {
                            EditOperation::Equal
                        } else {
                            EditOperation::Substitution
                        },
                    );
                    let previous = cell.previous[usize::from(STATE_MATCH)];
                    query_index -= 1;
                    target_index -= 1;
                    state = previous;
                }
                STATE_INSERTION => {
                    if target_index == 0 {
                        return Err(AlignmentError::InvalidTraceback);
                    }
                    operations.push(EditOperation::Insertion);
                    let previous = cell.previous[usize::from(STATE_INSERTION)];
                    target_index -= 1;
                    state = previous;
                }
                STATE_DELETION => {
                    operations.push(EditOperation::Deletion);
                    let previous = cell.previous[usize::from(STATE_DELETION)];
                    query_index -= 1;
                    state = previous;
                }
                _ => return Err(AlignmentError::InvalidTraceback),
            }
        }
        operations.reverse();
        let _ = (query, target);
        Ok((target_index, operations))
    }

    fn cell_index(&self, query_index: usize, target_index: usize) -> usize {
        self.row_offsets[query_index] + target_index - self.row_starts[query_index]
    }

    fn cell_index_checked(&self, query_index: usize, target_index: usize) -> Option<usize> {
        let start = self.row_starts.get(query_index).copied()?;
        let width = self.row_widths.get(query_index).copied()?;
        if target_index < start || target_index >= start + width {
            return None;
        }
        Some(self.row_offsets[query_index] + target_index - start)
    }

    fn traceback(
        &mut self,
        query: &[u8],
        target: &[u8],
        best: BestCell,
    ) -> Result<(usize, usize, u8), AlignmentError> {
        let mut query_index = best.query_index;
        let mut target_index = best.target_index;
        let mut state = best.state;
        self.operations.clear();

        while query_index > 0 || target_index > 0 {
            let index = self
                .cell_index_checked(query_index, target_index)
                .ok_or(AlignmentError::TracebackOutsideBand)?;
            let cell = self.cells[index];
            let state_index = usize::from(state);
            if state_index >= 3 || cell.scores[state_index] <= 0 {
                break;
            }
            match state {
                STATE_MATCH => {
                    if query_index == 0 || target_index == 0 {
                        return Err(AlignmentError::InvalidTraceback);
                    }
                    let operation = if equal_base(query[query_index - 1], target[target_index - 1])
                    {
                        EditOperation::Equal
                    } else {
                        EditOperation::Substitution
                    };
                    self.operations.push(operation);
                    let previous = cell.previous[usize::from(STATE_MATCH)];
                    query_index -= 1;
                    target_index -= 1;
                    if previous == STATE_START {
                        break;
                    }
                    state = previous;
                }
                STATE_INSERTION => {
                    if target_index == 0 {
                        return Err(AlignmentError::InvalidTraceback);
                    }
                    self.operations.push(EditOperation::Insertion);
                    let previous = cell.previous[usize::from(STATE_INSERTION)];
                    target_index -= 1;
                    if previous == STATE_START {
                        break;
                    }
                    state = previous;
                }
                STATE_DELETION => {
                    if query_index == 0 {
                        return Err(AlignmentError::InvalidTraceback);
                    }
                    self.operations.push(EditOperation::Deletion);
                    let previous = cell.previous[usize::from(STATE_DELETION)];
                    query_index -= 1;
                    if previous == STATE_START {
                        break;
                    }
                    state = previous;
                }
                _ => return Err(AlignmentError::InvalidTraceback),
            }
        }

        self.operations.reverse();
        Ok((query_index, target_index, state))
    }
}

/// Result of one local alignment before it is attached to workflow IDs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AlignmentResult {
    /// Optimal local alignment score under the supplied scoring scheme.
    pub score: i32,
    /// Strand of the forward-stored target that was aligned.
    pub strand: Strand,
    /// Query interval in the supplied linear query window.
    pub query_interval: BaseInterval,
    /// Query interval(s) in the original circular coordinate system.
    pub query_segments: Vec<BaseInterval>,
    /// Target interval in forward stored-contig coordinates.
    pub target_interval: BaseInterval,
    /// True only when `query_segments` contains the two sides of origin.
    pub origin_crossing: bool,
    pub query_length: u64,
    pub target_length: u64,
    pub matches: u64,
    pub substitutions: u64,
    pub insertions: u64,
    pub deletions: u64,
    pub cigar: String,
    pub edit_script: Vec<EditRun>,
    /// Chain score is populated by the caller that supplied anchors.
    pub chain_score: i64,
    /// Bounded retry, corridor, and refinement evidence for this result.
    pub retry_metadata: AlignmentRetryMetadata,
}

impl AlignmentResult {
    /// Convert this evidence into the frozen workflow exchange model.
    #[must_use]
    pub fn into_base_alignment(
        self,
        plasmid_id: impl Into<String>,
        metagenome_id: impl Into<String>,
        contig_id: impl Into<String>,
        chain_score: i64,
        primary: bool,
    ) -> BaseAlignment {
        let identity_denominator = self
            .matches
            .saturating_add(self.substitutions)
            .saturating_add(self.insertions)
            .saturating_add(self.deletions);
        let identity = if identity_denominator == 0 {
            0.0
        } else {
            self.matches as f64 / identity_denominator as f64
        };
        BaseAlignment {
            alignment_id: String::new(),
            plasmid_id: plasmid_id.into(),
            metagenome_id: metagenome_id.into(),
            contig_id: contig_id.into(),
            strand: self.strand,
            query_segments: self.query_segments,
            target_interval: self.target_interval,
            query_length: self.query_length,
            target_length: self.target_length,
            origin_crossing: self.origin_crossing,
            score: i64::from(self.score),
            matches: self.matches,
            substitutions: self.substitutions,
            insertions: self.insertions,
            deletions: self.deletions,
            cigar: self.cigar,
            edit_script: self.edit_script,
            chain_score,
            identity,
            seed_evidence: SeedEvidence::default(),
            primary_supported_bases: 0,
            secondary_supported_bases: 0,
            newly_supported_bases: 0,
            role: AlignmentRole::AlternativeMapping,
            primary,
        }
    }

    /// Reconstruct the aligned strings from the original linear inputs.
    /// `query` and `target` must be the same sequences passed to `align`.
    pub fn reconstruct(
        &self,
        query: &[u8],
        target: &[u8],
    ) -> Result<(Vec<u8>, Vec<u8>), AlignmentError> {
        let query_start = usize::try_from(self.query_interval.start)
            .map_err(|_| AlignmentError::LengthOverflow)?;
        let query_end =
            usize::try_from(self.query_interval.end).map_err(|_| AlignmentError::LengthOverflow)?;
        let target_start = usize::try_from(self.target_interval.start)
            .map_err(|_| AlignmentError::LengthOverflow)?;
        let target_end = usize::try_from(self.target_interval.end)
            .map_err(|_| AlignmentError::LengthOverflow)?;
        if query_end > query.len() || target_end > target.len() || query_start > query_end {
            return Err(AlignmentError::ReconstructionBounds);
        }
        let mut query_index = query_start;
        let mut target_index = target_start;
        let mut aligned_query = Vec::with_capacity(self.query_length as usize);
        let mut aligned_target = Vec::with_capacity(self.target_length as usize);
        for run in &self.edit_script {
            for _ in 0..run.length {
                match run.operation {
                    EditOperation::Equal | EditOperation::Substitution => {
                        if query_index >= query_end || target_index >= target_end {
                            return Err(AlignmentError::ReconstructionBounds);
                        }
                        aligned_query.push(query[query_index]);
                        aligned_target.push(target[target_index]);
                        query_index += 1;
                        target_index += 1;
                    }
                    EditOperation::Insertion => {
                        if target_index >= target_end {
                            return Err(AlignmentError::ReconstructionBounds);
                        }
                        aligned_query.push(b'-');
                        aligned_target.push(target[target_index]);
                        target_index += 1;
                    }
                    EditOperation::Deletion => {
                        if query_index >= query_end {
                            return Err(AlignmentError::ReconstructionBounds);
                        }
                        aligned_query.push(query[query_index]);
                        aligned_target.push(b'-');
                        query_index += 1;
                    }
                    EditOperation::SoftClip => {
                        return Err(AlignmentError::UnsupportedOperation(
                            EditOperation::SoftClip,
                        ));
                    }
                }
            }
        }
        if query_index != query_end || target_index != target_end {
            return Err(AlignmentError::ReconstructionBounds);
        }
        Ok((aligned_query, aligned_target))
    }
}

/// Align with a fresh reusable workspace.
pub fn align(
    query: &[u8],
    target: &[u8],
    options: impl Into<AlignmentOptions>,
) -> Result<AlignmentResult, AlignmentError> {
    AlignmentWorkspace::new().align(query, target, options)
}

/// Align through a chain-guided corridor with bounded automatic retries.
pub fn align_with_retries(
    query: &[u8],
    target: &[u8],
    corridor: &AlignmentCorridor,
    options: impl Into<AlignmentOptions>,
) -> Result<AlignmentResult, AlignmentError> {
    AlignmentWorkspace::new().align_with_retries(query, target, corridor, options)
}

/// Align a chain-supported window and retain genuine fragment-end clipping.
pub fn align_anchored_semiglobal(
    query: &[u8],
    target: &[u8],
    corridor: &AlignmentCorridor,
    options: impl Into<AlignmentOptions>,
) -> Result<AlignmentResult, AlignmentError> {
    AlignmentWorkspace::new().align_anchored_semiglobal(query, target, corridor, options)
}

/// Align a forward-stored target and return a forward target interval.
pub fn align_oriented(
    query: &[u8],
    target_forward: &[u8],
    target_start: u64,
    strand: Strand,
    options: impl Into<AlignmentOptions>,
) -> Result<AlignmentResult, AlignmentError> {
    AlignmentWorkspace::new().align_oriented(query, target_forward, target_start, strand, options)
}

/// Align a one-wrap circular query window with a forward-stored target.
pub fn align_circular(
    query: &[u8],
    query_start: u64,
    query_span: u64,
    target_forward: &[u8],
    target_start: u64,
    strand: Strand,
    options: impl Into<AlignmentOptions>,
) -> Result<AlignmentResult, AlignmentError> {
    AlignmentWorkspace::new().align_circular(
        query,
        query_start,
        query_span,
        target_forward,
        target_start,
        strand,
        options,
    )
}

#[derive(Clone, Copy, Debug, Default)]
struct Cell {
    scores: [i32; 3],
    previous: [u8; 3],
}

const SEMI_NEG_INF: i32 = i32::MIN / 4;

#[derive(Clone, Copy, Debug)]
struct SemiCell {
    scores: [i32; 3],
    previous: [u8; 3],
}

impl Default for SemiCell {
    fn default() -> Self {
        Self {
            scores: [SEMI_NEG_INF; 3],
            previous: [STATE_UNREACHABLE; 3],
        }
    }
}

impl SemiCell {
    fn best_score(self) -> (i32, u8) {
        choose_semi_transition([
            (self.scores[usize::from(STATE_MATCH)], STATE_MATCH),
            (self.scores[usize::from(STATE_INSERTION)], STATE_INSERTION),
            (self.scores[usize::from(STATE_DELETION)], STATE_DELETION),
        ])
    }
}

impl Cell {
    fn best_score(self) -> (i32, u8) {
        choose_transition([
            (self.scores[usize::from(STATE_MATCH)], STATE_MATCH),
            (self.scores[usize::from(STATE_INSERTION)], STATE_INSERTION),
            (self.scores[usize::from(STATE_DELETION)], STATE_DELETION),
        ])
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct BestCell {
    score: i32,
    query_index: usize,
    target_index: usize,
    state: u8,
}

impl BestCell {
    fn consider(&mut self, query_index: usize, target_index: usize, cell: Cell) {
        let (score, state) = cell.best_score();
        if score > self.score
            || (score == self.score
                && score > 0
                && (query_index < self.query_index
                    || (query_index == self.query_index
                        && (target_index < self.target_index
                            || (target_index == self.target_index && state < self.state)))))
        {
            self.score = score;
            self.query_index = query_index;
            self.target_index = target_index;
            self.state = state;
        }
    }
}

struct RawAlignment {
    score: i32,
    query_interval: BaseInterval,
    target_interval: BaseInterval,
    matches: u64,
    substitutions: u64,
    insertions: u64,
    deletions: u64,
    cigar: String,
    edit_script: Vec<EditRun>,
    band_edge_touched: bool,
    band_width: u32,
}

impl RawAlignment {
    fn finish_linear(
        self,
        strand: Strand,
        target_start: u64,
        target_len: u64,
    ) -> Result<AlignmentResult, AlignmentError> {
        self.finish_oriented(strand, target_start, target_len, false)
    }

    fn finish_oriented(
        self,
        strand: Strand,
        target_start: u64,
        target_len: u64,
        target_was_reversed: bool,
    ) -> Result<AlignmentResult, AlignmentError> {
        let target_interval = if target_was_reversed {
            let start = target_len
                .checked_sub(self.target_interval.end)
                .ok_or(AlignmentError::LengthOverflow)?;
            let end = target_len
                .checked_sub(self.target_interval.start)
                .ok_or(AlignmentError::LengthOverflow)?;
            BaseInterval::new(
                target_start
                    .checked_add(start)
                    .ok_or(AlignmentError::LengthOverflow)?,
                target_start
                    .checked_add(end)
                    .ok_or(AlignmentError::LengthOverflow)?,
            )
            .map_err(|_| AlignmentError::LengthOverflow)?
        } else {
            BaseInterval::new(
                target_start
                    .checked_add(self.target_interval.start)
                    .ok_or(AlignmentError::LengthOverflow)?,
                target_start
                    .checked_add(self.target_interval.end)
                    .ok_or(AlignmentError::LengthOverflow)?,
            )
            .map_err(|_| AlignmentError::LengthOverflow)?
        };
        Ok(AlignmentResult {
            score: self.score,
            strand,
            query_interval: self.query_interval,
            query_segments: vec![self.query_interval],
            target_interval,
            origin_crossing: false,
            query_length: self.query_interval.len(),
            target_length: self.target_interval.len(),
            matches: self.matches,
            substitutions: self.substitutions,
            insertions: self.insertions,
            deletions: self.deletions,
            cigar: self.cigar,
            edit_script: self.edit_script,
            chain_score: 0,
            retry_metadata: AlignmentRetryMetadata {
                attempted_widths: vec![self.band_width],
                selected_width: self.band_width,
                retries: 0,
                band_edge_touched: self.band_edge_touched,
                ..AlignmentRetryMetadata::default()
            },
        })
    }
}

fn choose_transition<const N: usize>(transitions: [(i32, u8); N]) -> (i32, u8) {
    let mut best = (0, STATE_UNREACHABLE);
    for (score, state) in transitions {
        if score > best.0 || (score == best.0 && state < best.1) {
            best = (score, state);
        }
    }
    best
}

fn choose_semi_transition<const N: usize>(transitions: [(i32, u8); N]) -> (i32, u8) {
    let mut best = (SEMI_NEG_INF, STATE_UNREACHABLE);
    for (score, state) in transitions {
        if score > best.0 || (score == best.0 && state < best.1) {
            best = (score, state);
        }
    }
    best
}

fn semi_add(score: i32, delta: i32) -> i32 {
    if score <= SEMI_NEG_INF {
        SEMI_NEG_INF
    } else {
        score.saturating_add(delta)
    }
}

fn gap_open_score(scoring: AlignmentScoring) -> i32 {
    scoring
        .gap_open_score
        .saturating_add(scoring.gap_extend_score)
}

fn equal_base(left: u8, right: u8) -> bool {
    left.eq_ignore_ascii_case(&right)
}

fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        _ => b'N',
    }
}

fn map_target_coordinates(
    mut result: AlignmentResult,
    target_len: usize,
    target_start: u64,
    strand: Strand,
) -> Result<AlignmentResult, AlignmentError> {
    let relative = result.target_interval;
    let target_len = u64::try_from(target_len).map_err(|_| AlignmentError::LengthOverflow)?;
    let mapped = match strand {
        Strand::Forward => BaseInterval::new(
            target_start
                .checked_add(relative.start)
                .ok_or(AlignmentError::LengthOverflow)?,
            target_start
                .checked_add(relative.end)
                .ok_or(AlignmentError::LengthOverflow)?,
        )
        .map_err(|_| AlignmentError::LengthOverflow)?,
        Strand::Reverse => {
            let start = target_len
                .checked_sub(relative.end)
                .ok_or(AlignmentError::LengthOverflow)?;
            let end = target_len
                .checked_sub(relative.start)
                .ok_or(AlignmentError::LengthOverflow)?;
            BaseInterval::new(
                target_start
                    .checked_add(start)
                    .ok_or(AlignmentError::LengthOverflow)?,
                target_start
                    .checked_add(end)
                    .ok_or(AlignmentError::LengthOverflow)?,
            )
            .map_err(|_| AlignmentError::LengthOverflow)?
        }
    };
    result.target_interval = mapped;
    result.strand = strand;
    Ok(result)
}

fn circular_segments(
    absolute_start: u64,
    absolute_end: u64,
    circular_length: u64,
) -> Result<Vec<BaseInterval>, AlignmentError> {
    if absolute_start > absolute_end || circular_length == 0 {
        return Err(AlignmentError::LengthOverflow);
    }
    let start = absolute_start % circular_length;
    if absolute_end <= circular_length {
        return Ok(vec![
            BaseInterval::new(start, absolute_end).map_err(|_| AlignmentError::LengthOverflow)?,
        ]);
    }
    let end = absolute_end % circular_length;
    if start == end {
        return Ok(vec![
            BaseInterval::new(start, circular_length)
                .map_err(|_| AlignmentError::LengthOverflow)?,
            BaseInterval::new(0, start).map_err(|_| AlignmentError::LengthOverflow)?,
        ]);
    }
    if start <= end {
        return Ok(vec![
            BaseInterval::new(start, end).map_err(|_| AlignmentError::LengthOverflow)?,
        ]);
    }
    Ok(vec![
        BaseInterval::new(start, circular_length).map_err(|_| AlignmentError::LengthOverflow)?,
        BaseInterval::new(0, end).map_err(|_| AlignmentError::LengthOverflow)?,
    ])
}

type OperationSummary = (Vec<EditRun>, String, u64, u64, u64, u64);

fn summarize_operations(operations: &[EditOperation]) -> Result<OperationSummary, AlignmentError> {
    let mut runs: Vec<EditRun> = Vec::new();
    for &operation in operations {
        if let Some(last) = runs.last_mut()
            && last.operation == operation
        {
            last.length = last
                .length
                .checked_add(1)
                .ok_or(AlignmentError::RunTooLong)?;
            continue;
        }
        runs.push(EditRun {
            operation,
            length: 1,
        });
    }

    let mut cigar = String::new();
    let mut matches = 0u64;
    let mut substitutions = 0u64;
    let mut insertions = 0u64;
    let mut deletions = 0u64;
    for run in &runs {
        let code = match run.operation {
            EditOperation::Equal => {
                matches = matches.saturating_add(u64::from(run.length));
                '='
            }
            EditOperation::Substitution => {
                substitutions = substitutions.saturating_add(u64::from(run.length));
                'X'
            }
            EditOperation::Insertion => {
                insertions = insertions.saturating_add(u64::from(run.length));
                'I'
            }
            EditOperation::Deletion => {
                deletions = deletions.saturating_add(u64::from(run.length));
                'D'
            }
            EditOperation::SoftClip => 'S',
        };
        write!(&mut cigar, "{}{}", run.length, code).map_err(|_| AlignmentError::CigarWrite)?;
    }
    Ok((runs, cigar, matches, substitutions, insertions, deletions))
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AlignmentError {
    #[error("alignment query is empty")]
    EmptyQuery,
    #[error("alignment target is empty")]
    EmptyTarget,
    #[error("alignment inputs do not fit in the requested matrix")]
    LengthOverflow,
    #[error("alignment band excludes every input cell")]
    BandExcludesInput,
    #[error("invalid alignment corridor: {0}")]
    InvalidCorridor(&'static str),
    #[error("alignment matrix requires {cells} cells, exceeding limit {max_cells}")]
    MatrixTooLarge { cells: usize, max_cells: usize },
    #[error("no positive-scoring local alignment was found")]
    NoAlignment,
    #[error("traceback left the alignment band")]
    TracebackOutsideBand,
    #[error("alignment traceback is inconsistent")]
    InvalidTraceback,
    #[error("circular query start {start} is outside query length {length}")]
    CircularStart { start: u64, length: u64 },
    #[error("circular query span {span} is invalid for query length {length}")]
    CircularSpan { span: u64, length: u64 },
    #[error("alignment edit run exceeds u32 length")]
    RunTooLong,
    #[error("cannot reconstruct an alignment outside the supplied sequences")]
    ReconstructionBounds,
    #[error("alignment operation {0:?} is not reconstructible")]
    UnsupportedOperation(EditOperation),
    #[error("failed to write CIGAR text")]
    CigarWrite,
}

#[cfg(test)]
mod unit_tests {
    use super::*;

    fn scoring() -> AlignmentScoring {
        AlignmentScoring {
            match_score: 2,
            mismatch_score: -3,
            gap_open_score: -5,
            gap_extend_score: -1,
            band_width: 8,
        }
    }

    #[test]
    fn affine_gap_result_reconstructs() {
        let result = align(b"ACGTACGT", b"ACGTGACGT", scoring()).unwrap();
        assert_eq!(result.cigar, "4=1I4=");
        assert_eq!(result.insertions, 1);
        let (query, target) = result.reconstruct(b"ACGTACGT", b"ACGTGACGT").unwrap();
        assert_eq!(query, b"ACGT-ACGT");
        assert_eq!(target, b"ACGTGACGT");
    }

    #[test]
    fn reverse_target_interval_is_reported_in_forward_coordinates() {
        let target = b"AAGTCCT";
        let query = b"AGGACTT";
        let result = align_oriented(query, target, 100, Strand::Reverse, scoring()).unwrap();
        assert_eq!(result.strand, Strand::Reverse);
        assert_eq!(result.target_interval, BaseInterval::new(100, 107).unwrap());
        assert!(result.matches >= 5);
    }

    #[test]
    fn circular_alignment_has_two_query_segments() {
        let result =
            align_circular(b"ACGTAC", 4, 4, b"ACAC", 0, Strand::Forward, scoring()).unwrap();
        assert!(result.origin_crossing);
        assert_eq!(
            result.query_segments,
            vec![
                BaseInterval::new(4, 6).unwrap(),
                BaseInterval::new(0, 2).unwrap()
            ]
        );
    }

    #[test]
    fn workspace_reuses_matrix_capacity() {
        let mut workspace = AlignmentWorkspace::new();
        let _ = workspace.align(b"ACGT", b"ACGT", scoring()).unwrap();
        let first = workspace.capacity_cells();
        let _ = workspace.align(b"ACG", b"ACG", scoring()).unwrap();
        assert_eq!(workspace.capacity_cells(), first);
    }

    #[test]
    fn rows_after_target_end_are_empty_not_out_of_bounds() {
        let query = vec![b'A'; 80];
        let target = vec![b'A'; 64];
        let result = align(&query, &target, scoring()).unwrap();
        assert_eq!(result.matches, 64);
    }
}
