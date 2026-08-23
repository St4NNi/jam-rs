//! Banded affine-gap alignment for trace candidates.
//!
//! The trace workflow aligns a plasmid query against a bounded sequence window
//! from an assembly. This module deliberately exposes alignment evidence
//! (intervals, operations, and counters) rather than a presence call. The
//! dynamic-programming matrix is banded and owned by [`AlignmentWorkspace`] so
//! callers can reuse its allocations for a stream of candidate windows.

use crate::trace::config::AlignmentScoring;
use crate::trace::model::{BaseAlignment, BaseInterval, EditOperation, EditRun, Strand};
use std::fmt::Write as _;
use thiserror::Error;

const STATE_MATCH: u8 = 0;
const STATE_INSERTION: u8 = 1;
const STATE_DELETION: u8 = 2;
const STATE_START: u8 = 3;
const STATE_UNREACHABLE: u8 = 4;

/// Default upper bound on the number of band cells allocated by an alignment.
/// Callers handling larger windows should set an explicit bound and split the
/// windows before aligning them.
pub const DEFAULT_MAX_CELLS: usize = 4_000_000;

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

    fn align_raw(
        &mut self,
        query: &[u8],
        target: &[u8],
        options: &AlignmentOptions,
    ) -> Result<RawAlignment, AlignmentError> {
        if query.is_empty() {
            return Err(AlignmentError::EmptyQuery);
        }
        if target.is_empty() {
            return Err(AlignmentError::EmptyTarget);
        }
        let band = usize::try_from(options.scoring.band_width)
            .map_err(|_| AlignmentError::LengthOverflow)?;
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
            let center = (query_index as isize)
                .checked_add(options.diagonal_offset)
                .ok_or(AlignmentError::LengthOverflow)?;
            let start = center.saturating_sub(band as isize).max(0) as usize;
            let end = center
                .saturating_add(band as isize)
                .max(0)
                .min(target_len as isize) as usize;
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
        })
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
        BaseAlignment {
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
