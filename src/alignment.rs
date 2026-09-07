use serde::ser::SerializeStruct;
use serde::{Deserialize, Serialize};
use std::fmt::Write as _;
use thiserror::Error;

const MATCH: u8 = 0;
const INSERTION: u8 = 1;
const DELETION: u8 = 2;
const START: u8 = 3;
const UNREACHABLE: u8 = 4;

pub const DEFAULT_MAX_CELLS: usize = 4_000_000;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Interval {
    pub start: u64,
    pub end: u64,
}

impl Interval {
    pub fn new(start: u64, end: u64) -> Result<Self, AlignmentError> {
        if start > end {
            return Err(AlignmentError::ReversedInterval { start, end });
        }
        Ok(Self { start, end })
    }

    pub fn len(self) -> u64 {
        self.end - self.start
    }

    pub fn is_empty(self) -> bool {
        self.start == self.end
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum EditOperation {
    #[serde(rename = "=")]
    Equal,
    #[serde(rename = "X")]
    Substitution,
    #[serde(rename = "I")]
    Insertion,
    #[serde(rename = "D")]
    Deletion,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct EditRun {
    pub operation: EditOperation,
    pub length: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AlignmentConfig {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open_score: i32,
    pub gap_extend_score: i32,
    pub band_width: u32,
    pub diagonal_offset: i64,
    pub max_cells: usize,
}

impl Default for AlignmentConfig {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -3,
            gap_open_score: -5,
            gap_extend_score: -1,
            band_width: 128,
            diagonal_offset: 0,
            max_cells: DEFAULT_MAX_CELLS,
        }
    }
}

impl AlignmentConfig {
    fn validate(self) -> Result<(), AlignmentError> {
        if self.match_score <= 0
            || self.mismatch_score > 0
            || self.gap_open_score > 0
            || self.gap_extend_score > 0
            || self.max_cells == 0
        {
            return Err(AlignmentError::InvalidConfig);
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Deserialize)]
pub struct Alignment {
    pub score: i32,
    pub strand: Strand,
    pub query_interval: Interval,
    pub target_interval: Interval,
    pub matches: u64,
    pub substitutions: u64,
    pub insertions: u64,
    pub deletions: u64,
    pub cigar: String,
    pub edit_script: Vec<EditRun>,
}

impl Serialize for Alignment {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut record = serializer.serialize_struct("Alignment", 11)?;
        record.serialize_field("score", &self.score)?;
        record.serialize_field("strand", &self.strand)?;
        record.serialize_field("query_interval", &self.query_interval)?;
        record.serialize_field("target_interval", &self.target_interval)?;
        record.serialize_field("matches", &self.matches)?;
        record.serialize_field("substitutions", &self.substitutions)?;
        record.serialize_field("insertions", &self.insertions)?;
        record.serialize_field("deletions", &self.deletions)?;
        record.serialize_field("cigar", &self.cigar)?;
        record.serialize_field("edit_script", &self.edit_script)?;
        record.serialize_field("identity", &self.identity())?;
        record.end()
    }
}

impl Alignment {
    pub fn identity(&self) -> f64 {
        let total = self
            .matches
            .saturating_add(self.substitutions)
            .saturating_add(self.insertions)
            .saturating_add(self.deletions);
        if total == 0 {
            0.0
        } else {
            self.matches as f64 / total as f64
        }
    }

    pub fn validate_cigar(&self) -> Result<(), AlignmentError> {
        let runs = parse_cigar(&self.cigar)?;
        if runs != self.edit_script {
            return Err(AlignmentError::CigarMismatch);
        }
        let summary = summarize_runs(&runs)?;
        if summary.query_bases != self.query_interval.len()
            || summary.target_bases != self.target_interval.len()
        {
            return Err(AlignmentError::CigarSpanMismatch);
        }
        if (
            summary.matches,
            summary.substitutions,
            summary.insertions,
            summary.deletions,
        ) != (
            self.matches,
            self.substitutions,
            self.insertions,
            self.deletions,
        ) {
            return Err(AlignmentError::CigarCountMismatch);
        }
        Ok(())
    }
}

#[derive(Debug, Default)]
pub struct AlignmentWorkspace {
    cells: Vec<Cell>,
    endpoint_cells: Vec<EndpointCell>,
    row_offsets: Vec<usize>,
    row_starts: Vec<usize>,
    row_widths: Vec<usize>,
    operations: Vec<EditOperation>,
    reverse: Vec<u8>,
}

impl AlignmentWorkspace {
    pub fn capacity_cells(&self) -> usize {
        self.cells.capacity()
    }

    pub fn align(
        &mut self,
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
    ) -> Result<Alignment, AlignmentError> {
        let raw = self.align_raw(query, target, config)?;
        finish(raw, Strand::Forward, 0, target.len())
    }

    pub fn align_oriented(
        &mut self,
        query: &[u8],
        target: &[u8],
        target_offset: u64,
        strand: Strand,
        config: AlignmentConfig,
    ) -> Result<Alignment, AlignmentError> {
        let raw = match strand {
            Strand::Forward => self.align_raw(query, target, config)?,
            Strand::Reverse => {
                let mut reverse = std::mem::take(&mut self.reverse);
                reverse.clear();
                reverse.reserve(target.len());
                reverse.extend(target.iter().rev().map(|base| complement(*base)));
                let result = self.align_raw(query, &reverse, config);
                self.reverse = reverse;
                result?
            }
        };
        finish(raw, strand, target_offset, target.len())
    }

    pub fn complete_endpoints(
        &mut self,
        core: Alignment,
        query: &[u8],
        target: &[u8],
        target_offset: u64,
        max_extension: usize,
        config: AlignmentConfig,
    ) -> Result<EndpointCompletion, AlignmentError> {
        config.validate()?;
        match core.strand {
            Strand::Forward => complete_endpoints(
                &mut self.endpoint_cells,
                core,
                query,
                target,
                target_offset,
                max_extension,
                config,
            ),
            Strand::Reverse => {
                let mut reverse = std::mem::take(&mut self.reverse);
                reverse.clear();
                reverse.reserve(target.len());
                reverse.extend(target.iter().rev().map(|base| complement(*base)));
                let result = complete_endpoints(
                    &mut self.endpoint_cells,
                    core,
                    query,
                    &reverse,
                    target_offset,
                    max_extension,
                    config,
                );
                self.reverse = reverse;
                result
            }
        }
    }

    fn align_raw(
        &mut self,
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
    ) -> Result<RawAlignment, AlignmentError> {
        config.validate()?;
        if query.is_empty() {
            return Err(AlignmentError::EmptyQuery);
        }
        if target.is_empty() {
            return Err(AlignmentError::EmptyTarget);
        }

        self.prepare_rows(query.len(), target.len(), config)?;
        let total_cells = self
            .row_offsets
            .last()
            .copied()
            .unwrap_or(0)
            .checked_add(self.row_widths.last().copied().unwrap_or(0))
            .ok_or(AlignmentError::LengthOverflow)?;
        if total_cells == 0 {
            return Err(AlignmentError::BandExcludesInput);
        }
        if total_cells > config.max_cells {
            return Err(AlignmentError::MatrixTooLarge {
                cells: total_cells,
                max_cells: config.max_cells,
            });
        }
        if self.cells.len() < total_cells {
            self.cells.resize(total_cells, Cell::default());
        } else {
            self.cells[..total_cells].fill(Cell::default());
            self.cells.truncate(total_cells);
        }

        let mut best = BestCell::default();
        for query_index in 0..=query.len() {
            let start = self.row_starts[query_index];
            let width = self.row_widths[query_index];
            for target_index in start..start + width {
                if query_index == 0 && target_index == 0 {
                    continue;
                }
                let mut cell = Cell::default();
                if query_index > 0
                    && target_index > 0
                    && let Some(previous) =
                        self.cell_index_checked(query_index - 1, target_index - 1)
                {
                    let previous = self.cells[previous];
                    let (score, state) = previous.best_score();
                    let substitution =
                        if query[query_index - 1].eq_ignore_ascii_case(&target[target_index - 1]) {
                            config.match_score
                        } else {
                            config.mismatch_score
                        };
                    let score = score.saturating_add(substitution);
                    if score > 0 {
                        cell.scores[MATCH as usize] = score;
                        cell.previous[MATCH as usize] =
                            if score == substitution { START } else { state };
                    }
                }
                if target_index > 0
                    && let Some(previous) = self.cell_index_checked(query_index, target_index - 1)
                {
                    let previous = self.cells[previous];
                    let (score, state) = choose([
                        (
                            previous.scores[INSERTION as usize]
                                .saturating_add(config.gap_extend_score),
                            INSERTION,
                        ),
                        (
                            previous.scores[MATCH as usize].saturating_add(gap_open(config)),
                            MATCH,
                        ),
                        (
                            previous.scores[DELETION as usize].saturating_add(gap_open(config)),
                            DELETION,
                        ),
                    ]);
                    if score > 0 {
                        cell.scores[INSERTION as usize] = score;
                        cell.previous[INSERTION as usize] = state;
                    }
                }
                if query_index > 0
                    && let Some(previous) = self.cell_index_checked(query_index - 1, target_index)
                {
                    let previous = self.cells[previous];
                    let (score, state) = choose([
                        (
                            previous.scores[DELETION as usize]
                                .saturating_add(config.gap_extend_score),
                            DELETION,
                        ),
                        (
                            previous.scores[MATCH as usize].saturating_add(gap_open(config)),
                            MATCH,
                        ),
                        (
                            previous.scores[INSERTION as usize].saturating_add(gap_open(config)),
                            INSERTION,
                        ),
                    ]);
                    if score > 0 {
                        cell.scores[DELETION as usize] = score;
                        cell.previous[DELETION as usize] = state;
                    }
                }
                let index = self.cell_index(query_index, target_index);
                self.cells[index] = cell;
                best.consider(query_index, target_index, cell);
            }
        }
        if best.score <= 0 {
            return Err(AlignmentError::NoAlignment);
        }

        let (query_start, target_start) = self.traceback(query, target, best)?;
        let edit_script = runs_from_operations(&self.operations)?;
        let summary = summarize_runs(&edit_script)?;
        Ok(RawAlignment {
            score: best.score,
            query_interval: Interval::new(query_start as u64, best.query_index as u64)?,
            target_interval: Interval::new(target_start as u64, best.target_index as u64)?,
            cigar: cigar_from_runs(&edit_script)?,
            edit_script,
            summary,
        })
    }

    fn prepare_rows(
        &mut self,
        query_len: usize,
        target_len: usize,
        config: AlignmentConfig,
    ) -> Result<(), AlignmentError> {
        self.row_offsets.clear();
        self.row_starts.clear();
        self.row_widths.clear();
        let target_len = i128::try_from(target_len).map_err(|_| AlignmentError::LengthOverflow)?;
        let band = i128::from(config.band_width);
        let diagonal = i128::from(config.diagonal_offset);
        let mut total = 0usize;
        for query_index in 0..=query_len {
            let center =
                i128::try_from(query_index).map_err(|_| AlignmentError::LengthOverflow)? + diagonal;
            let low = center - band;
            let high = center + band;
            let (start, width) = if high < 0 || low > target_len {
                (0, 0)
            } else {
                let start =
                    usize::try_from(low.max(0)).map_err(|_| AlignmentError::LengthOverflow)?;
                let end = usize::try_from(high.min(target_len))
                    .map_err(|_| AlignmentError::LengthOverflow)?;
                (start, end - start + 1)
            };
            self.row_offsets.push(total);
            self.row_starts.push(start);
            self.row_widths.push(width);
            total = total
                .checked_add(width)
                .ok_or(AlignmentError::LengthOverflow)?;
        }
        Ok(())
    }

    fn cell_index(&self, query_index: usize, target_index: usize) -> usize {
        self.row_offsets[query_index] + target_index - self.row_starts[query_index]
    }

    fn cell_index_checked(&self, query_index: usize, target_index: usize) -> Option<usize> {
        let start = *self.row_starts.get(query_index)?;
        let width = *self.row_widths.get(query_index)?;
        (target_index >= start && target_index < start + width)
            .then(|| self.row_offsets[query_index] + target_index - start)
    }

    fn traceback(
        &mut self,
        query: &[u8],
        target: &[u8],
        best: BestCell,
    ) -> Result<(usize, usize), AlignmentError> {
        let mut query_index = best.query_index;
        let mut target_index = best.target_index;
        let mut state = best.state;
        self.operations.clear();
        while query_index > 0 || target_index > 0 {
            let index = self
                .cell_index_checked(query_index, target_index)
                .ok_or(AlignmentError::TracebackOutsideBand)?;
            let cell = self.cells[index];
            if state > DELETION || cell.scores[state as usize] <= 0 {
                break;
            }
            let previous = cell.previous[state as usize];
            match state {
                MATCH if query_index > 0 && target_index > 0 => {
                    self.operations.push(
                        if query[query_index - 1].eq_ignore_ascii_case(&target[target_index - 1]) {
                            EditOperation::Equal
                        } else {
                            EditOperation::Substitution
                        },
                    );
                    query_index -= 1;
                    target_index -= 1;
                }
                INSERTION if target_index > 0 => {
                    self.operations.push(EditOperation::Insertion);
                    target_index -= 1;
                }
                DELETION if query_index > 0 => {
                    self.operations.push(EditOperation::Deletion);
                    query_index -= 1;
                }
                _ => return Err(AlignmentError::InvalidTraceback),
            }
            if previous == START {
                break;
            }
            state = previous;
        }
        self.operations.reverse();
        Ok((query_index, target_index))
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct Cell {
    scores: [i32; 3],
    previous: [u8; 3],
}

impl Cell {
    fn best_score(self) -> (i32, u8) {
        choose([
            (self.scores[MATCH as usize], MATCH),
            (self.scores[INSERTION as usize], INSERTION),
            (self.scores[DELETION as usize], DELETION),
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
                && (query_index, target_index, state)
                    < (self.query_index, self.target_index, self.state))
        {
            *self = Self {
                score,
                query_index,
                target_index,
                state,
            };
        }
    }
}

struct RawAlignment {
    score: i32,
    query_interval: Interval,
    target_interval: Interval,
    cigar: String,
    edit_script: Vec<EditRun>,
    summary: RunSummary,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct EndpointMetrics {
    pub left_query_bases: usize,
    pub left_target_bases: usize,
    pub right_query_bases: usize,
    pub right_target_bases: usize,
    pub matrix_cells: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EndpointCompletion {
    pub alignment: Alignment,
    pub metrics: EndpointMetrics,
}

const NEGATIVE: i32 = i32::MIN / 4;

#[derive(Clone, Copy, Debug)]
struct EndpointCell {
    scores: [i32; 3],
    previous: [u8; 3],
}

impl Default for EndpointCell {
    fn default() -> Self {
        Self {
            scores: [NEGATIVE; 3],
            previous: [START; 3],
        }
    }
}

struct EndpointResult {
    runs: Vec<EditRun>,
    query_bases: usize,
    target_bases: usize,
    matrix_cells: usize,
}

fn complete_endpoints(
    cells: &mut Vec<EndpointCell>,
    mut core: Alignment,
    query: &[u8],
    target: &[u8],
    target_offset: u64,
    max_extension: usize,
    config: AlignmentConfig,
) -> Result<EndpointCompletion, AlignmentError> {
    let forward_target_len = target.len();
    core.validate_cigar()?;
    let query_start =
        usize::try_from(core.query_interval.start).map_err(|_| AlignmentError::LengthOverflow)?;
    let query_end =
        usize::try_from(core.query_interval.end).map_err(|_| AlignmentError::LengthOverflow)?;
    let local_start = core
        .target_interval
        .start
        .checked_sub(target_offset)
        .and_then(|value| usize::try_from(value).ok())
        .ok_or(AlignmentError::EndpointOutsideWindow)?;
    let local_end = core
        .target_interval
        .end
        .checked_sub(target_offset)
        .and_then(|value| usize::try_from(value).ok())
        .ok_or(AlignmentError::EndpointOutsideWindow)?;
    if query_start > query_end
        || query_end > query.len()
        || local_start > local_end
        || local_end > forward_target_len
    {
        return Err(AlignmentError::EndpointOutsideWindow);
    }
    let (target_start, target_end) = match core.strand {
        Strand::Forward => (local_start, local_end),
        Strand::Reverse => (
            forward_target_len - local_end,
            forward_target_len - local_start,
        ),
    };

    let left_query_len = query_start.min(max_extension);
    let left_target_len = target_start.min(max_extension);
    let left_query: Vec<_> = query[query_start - left_query_len..query_start]
        .iter()
        .rev()
        .copied()
        .collect();
    let left_target: Vec<_> = target[target_start - left_target_len..target_start]
        .iter()
        .rev()
        .copied()
        .collect();
    let mut left = anchored_semiglobal(cells, &left_query, &left_target, config)?;
    left.runs.reverse();

    let query_limit = query_end.saturating_add(max_extension).min(query.len());
    let target_limit = target_end.saturating_add(max_extension).min(target.len());
    let right = anchored_semiglobal(
        cells,
        &query[query_end..query_limit],
        &target[target_end..target_limit],
        config,
    )?;

    let mut runs = Vec::with_capacity(left.runs.len() + core.edit_script.len() + right.runs.len());
    append_runs(&mut runs, left.runs)?;
    append_runs(&mut runs, core.edit_script)?;
    append_runs(&mut runs, right.runs)?;
    let summary = summarize_runs(&runs)?;
    core.score = score_runs(&runs, config);
    core.query_interval = Interval::new(
        u64::try_from(query_start - left.query_bases)
            .map_err(|_| AlignmentError::LengthOverflow)?,
        u64::try_from(query_end + right.query_bases).map_err(|_| AlignmentError::LengthOverflow)?,
    )?;
    let oriented_start = target_start - left.target_bases;
    let oriented_end = target_end + right.target_bases;
    let (forward_start, forward_end) = match core.strand {
        Strand::Forward => (oriented_start, oriented_end),
        Strand::Reverse => (
            forward_target_len - oriented_end,
            forward_target_len - oriented_start,
        ),
    };
    core.target_interval = Interval::new(
        target_offset
            .checked_add(u64::try_from(forward_start).map_err(|_| AlignmentError::LengthOverflow)?)
            .ok_or(AlignmentError::LengthOverflow)?,
        target_offset
            .checked_add(u64::try_from(forward_end).map_err(|_| AlignmentError::LengthOverflow)?)
            .ok_or(AlignmentError::LengthOverflow)?,
    )?;
    core.matches = summary.matches;
    core.substitutions = summary.substitutions;
    core.insertions = summary.insertions;
    core.deletions = summary.deletions;
    core.cigar = cigar_from_runs(&runs)?;
    core.edit_script = runs;
    core.validate_cigar()?;

    Ok(EndpointCompletion {
        alignment: core,
        metrics: EndpointMetrics {
            left_query_bases: left.query_bases,
            left_target_bases: left.target_bases,
            right_query_bases: right.query_bases,
            right_target_bases: right.target_bases,
            matrix_cells: left.matrix_cells.saturating_add(right.matrix_cells),
        },
    })
}

fn anchored_semiglobal(
    cells: &mut Vec<EndpointCell>,
    query: &[u8],
    target: &[u8],
    config: AlignmentConfig,
) -> Result<EndpointResult, AlignmentError> {
    if query.is_empty() || target.is_empty() {
        return Ok(EndpointResult {
            runs: Vec::new(),
            query_bases: 0,
            target_bases: 0,
            matrix_cells: 0,
        });
    }
    let rows = query
        .len()
        .checked_add(1)
        .ok_or(AlignmentError::LengthOverflow)?;
    let columns = target
        .len()
        .checked_add(1)
        .ok_or(AlignmentError::LengthOverflow)?;
    let matrix_cells = rows
        .checked_mul(columns)
        .ok_or(AlignmentError::LengthOverflow)?;
    if matrix_cells > config.max_cells {
        return Err(AlignmentError::MatrixTooLarge {
            cells: matrix_cells,
            max_cells: config.max_cells,
        });
    }
    cells.resize(matrix_cells, EndpointCell::default());
    cells.fill(EndpointCell::default());
    cells[0].scores[MATCH as usize] = 0;
    for (column, cell) in cells.iter_mut().enumerate().take(columns).skip(1) {
        cell.scores[INSERTION as usize] = config.gap_open_score.saturating_add(
            config
                .gap_extend_score
                .saturating_mul(i32::try_from(column).unwrap_or(i32::MAX)),
        );
        cell.previous[INSERTION as usize] = if column == 1 { MATCH } else { INSERTION };
    }
    for row in 1..rows {
        let first = row * columns;
        cells[first].scores[DELETION as usize] = config.gap_open_score.saturating_add(
            config
                .gap_extend_score
                .saturating_mul(i32::try_from(row).unwrap_or(i32::MAX)),
        );
        cells[first].previous[DELETION as usize] = if row == 1 { MATCH } else { DELETION };
        for column in 1..columns {
            let index = first + column;
            let (score, state) = maximum(cells[index - columns - 1].scores);
            cells[index].scores[MATCH as usize] = score.saturating_add(
                if query[row - 1].eq_ignore_ascii_case(&target[column - 1]) {
                    config.match_score
                } else {
                    config.mismatch_score
                },
            );
            cells[index].previous[MATCH as usize] = state;

            let above = cells[index - columns].scores;
            let (score, state) = maximum([
                above[MATCH as usize].saturating_add(gap_open(config)),
                above[INSERTION as usize].saturating_add(gap_open(config)),
                above[DELETION as usize].saturating_add(config.gap_extend_score),
            ]);
            cells[index].scores[DELETION as usize] = score;
            cells[index].previous[DELETION as usize] = state;

            let left = cells[index - 1].scores;
            let (score, state) = maximum([
                left[MATCH as usize].saturating_add(gap_open(config)),
                left[INSERTION as usize].saturating_add(config.gap_extend_score),
                left[DELETION as usize].saturating_add(gap_open(config)),
            ]);
            cells[index].scores[INSERTION as usize] = score;
            cells[index].previous[INSERTION as usize] = state;
        }
    }

    let mut best = (NEGATIVE, 0usize, 0usize, START);
    for row in 0..rows {
        for column in 0..columns {
            if (row + 1 != rows && column + 1 != columns) || (row == 0 && column == 0) {
                continue;
            }
            let (score, state) = maximum(cells[row * columns + column].scores);
            let candidate = (
                row.saturating_add(column),
                row,
                column,
                std::cmp::Reverse(state),
            );
            let current = (
                best.1.saturating_add(best.2),
                best.1,
                best.2,
                std::cmp::Reverse(best.3),
            );
            if score > best.0 || (score == best.0 && candidate > current) {
                best = (score, row, column, state);
            }
        }
    }
    if best.3 > DELETION {
        return Err(AlignmentError::NoAlignment);
    }

    let mut row = best.1;
    let mut column = best.2;
    let mut state = best.3;
    let mut operations = Vec::with_capacity(row.saturating_add(column));
    while row > 0 || column > 0 {
        let index = row * columns + column;
        let previous = cells[index].previous[state as usize];
        match state {
            MATCH if row > 0 && column > 0 => {
                operations.push(
                    if query[row - 1].eq_ignore_ascii_case(&target[column - 1]) {
                        EditOperation::Equal
                    } else {
                        EditOperation::Substitution
                    },
                );
                row -= 1;
                column -= 1;
            }
            INSERTION if column > 0 => {
                operations.push(EditOperation::Insertion);
                column -= 1;
            }
            DELETION if row > 0 => {
                operations.push(EditOperation::Deletion);
                row -= 1;
            }
            _ => return Err(AlignmentError::InvalidTraceback),
        }
        state = previous;
    }
    operations.reverse();
    Ok(EndpointResult {
        runs: runs_from_operations(&operations)?,
        query_bases: best.1,
        target_bases: best.2,
        matrix_cells,
    })
}

fn maximum(values: [i32; 3]) -> (i32, u8) {
    values
        .into_iter()
        .enumerate()
        .skip(1)
        .fold((values[0], MATCH), |best, (state, score)| {
            if score > best.0 {
                (score, state as u8)
            } else {
                best
            }
        })
}

fn append_runs(output: &mut Vec<EditRun>, runs: Vec<EditRun>) -> Result<(), AlignmentError> {
    for run in runs {
        if let Some(last) = output.last_mut()
            && last.operation == run.operation
        {
            last.length = last
                .length
                .checked_add(run.length)
                .ok_or(AlignmentError::RunTooLong)?;
        } else {
            output.push(run);
        }
    }
    Ok(())
}

fn score_runs(runs: &[EditRun], config: AlignmentConfig) -> i32 {
    runs.iter().fold(0, |score, run| {
        let length = i32::try_from(run.length).unwrap_or(i32::MAX);
        let delta = match run.operation {
            EditOperation::Equal => config.match_score.saturating_mul(length),
            EditOperation::Substitution => config.mismatch_score.saturating_mul(length),
            EditOperation::Insertion | EditOperation::Deletion => config
                .gap_open_score
                .saturating_add(config.gap_extend_score.saturating_mul(length)),
        };
        score.saturating_add(delta)
    })
}

fn finish(
    raw: RawAlignment,
    strand: Strand,
    target_offset: u64,
    target_len: usize,
) -> Result<Alignment, AlignmentError> {
    let target_len = u64::try_from(target_len).map_err(|_| AlignmentError::LengthOverflow)?;
    let (start, end) = if strand == Strand::Reverse {
        (
            target_len
                .checked_sub(raw.target_interval.end)
                .ok_or(AlignmentError::LengthOverflow)?,
            target_len
                .checked_sub(raw.target_interval.start)
                .ok_or(AlignmentError::LengthOverflow)?,
        )
    } else {
        (raw.target_interval.start, raw.target_interval.end)
    };
    let alignment = Alignment {
        score: raw.score,
        strand,
        query_interval: raw.query_interval,
        target_interval: Interval::new(
            target_offset
                .checked_add(start)
                .ok_or(AlignmentError::LengthOverflow)?,
            target_offset
                .checked_add(end)
                .ok_or(AlignmentError::LengthOverflow)?,
        )?,
        matches: raw.summary.matches,
        substitutions: raw.summary.substitutions,
        insertions: raw.summary.insertions,
        deletions: raw.summary.deletions,
        cigar: raw.cigar,
        edit_script: raw.edit_script,
    };
    alignment.validate_cigar()?;
    Ok(alignment)
}

fn choose<const N: usize>(values: [(i32, u8); N]) -> (i32, u8) {
    values.into_iter().fold((0, UNREACHABLE), |best, value| {
        if value.0 > best.0 || (value.0 == best.0 && value.1 < best.1) {
            value
        } else {
            best
        }
    })
}

fn gap_open(config: AlignmentConfig) -> i32 {
    config
        .gap_open_score
        .saturating_add(config.gap_extend_score)
}

fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
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

fn runs_from_operations(operations: &[EditOperation]) -> Result<Vec<EditRun>, AlignmentError> {
    let mut runs: Vec<EditRun> = Vec::new();
    for &operation in operations {
        if let Some(run) = runs.last_mut()
            && run.operation == operation
        {
            run.length = run
                .length
                .checked_add(1)
                .ok_or(AlignmentError::RunTooLong)?;
        } else {
            runs.push(EditRun {
                operation,
                length: 1,
            });
        }
    }
    Ok(runs)
}

pub fn parse_cigar(cigar: &str) -> Result<Vec<EditRun>, AlignmentError> {
    let mut runs = Vec::new();
    let mut length = 0u32;
    for byte in cigar.bytes() {
        if byte.is_ascii_digit() {
            length = length
                .checked_mul(10)
                .and_then(|value| value.checked_add(u32::from(byte - b'0')))
                .ok_or(AlignmentError::RunTooLong)?;
            continue;
        }
        if length == 0 {
            return Err(AlignmentError::InvalidCigar);
        }
        let operation = match byte {
            b'=' => EditOperation::Equal,
            b'X' => EditOperation::Substitution,
            b'I' => EditOperation::Insertion,
            b'D' => EditOperation::Deletion,
            _ => return Err(AlignmentError::InvalidCigar),
        };
        if runs
            .last()
            .is_some_and(|run: &EditRun| run.operation == operation)
        {
            return Err(AlignmentError::InvalidCigar);
        }
        runs.push(EditRun { operation, length });
        length = 0;
    }
    if length != 0 || runs.is_empty() {
        return Err(AlignmentError::InvalidCigar);
    }
    Ok(runs)
}

pub(crate) fn cigar_from_runs(runs: &[EditRun]) -> Result<String, AlignmentError> {
    let mut cigar = String::new();
    for run in runs {
        let operation = match run.operation {
            EditOperation::Equal => '=',
            EditOperation::Substitution => 'X',
            EditOperation::Insertion => 'I',
            EditOperation::Deletion => 'D',
        };
        write!(&mut cigar, "{}{operation}", run.length).map_err(|_| AlignmentError::CigarWrite)?;
    }
    Ok(cigar)
}

#[derive(Clone, Copy, Debug, Default)]
struct RunSummary {
    query_bases: u64,
    target_bases: u64,
    matches: u64,
    substitutions: u64,
    insertions: u64,
    deletions: u64,
}

fn summarize_runs(runs: &[EditRun]) -> Result<RunSummary, AlignmentError> {
    let mut summary = RunSummary::default();
    for run in runs {
        if run.length == 0 {
            return Err(AlignmentError::InvalidCigar);
        }
        let length = u64::from(run.length);
        match run.operation {
            EditOperation::Equal => {
                summary.query_bases += length;
                summary.target_bases += length;
                summary.matches += length;
            }
            EditOperation::Substitution => {
                summary.query_bases += length;
                summary.target_bases += length;
                summary.substitutions += length;
            }
            EditOperation::Insertion => {
                summary.target_bases += length;
                summary.insertions += length;
            }
            EditOperation::Deletion => {
                summary.query_bases += length;
                summary.deletions += length;
            }
        }
    }
    Ok(summary)
}

#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum AlignmentError {
    #[error("alignment query is empty")]
    EmptyQuery,
    #[error("alignment target is empty")]
    EmptyTarget,
    #[error("alignment configuration is invalid")]
    InvalidConfig,
    #[error("alignment length overflows")]
    LengthOverflow,
    #[error("alignment band excludes the input")]
    BandExcludesInput,
    #[error("alignment matrix requires {cells} cells, exceeding {max_cells}")]
    MatrixTooLarge { cells: usize, max_cells: usize },
    #[error("no positive-scoring local alignment was found")]
    NoAlignment,
    #[error("alignment traceback left the band")]
    TracebackOutsideBand,
    #[error("alignment traceback is inconsistent")]
    InvalidTraceback,
    #[error("alignment interval is reversed: {start}..{end}")]
    ReversedInterval { start: u64, end: u64 },
    #[error("alignment edit run exceeds u32")]
    RunTooLong,
    #[error("CIGAR text is invalid")]
    InvalidCigar,
    #[error("CIGAR text differs from the edit script")]
    CigarMismatch,
    #[error("CIGAR span differs from alignment coordinates")]
    CigarSpanMismatch,
    #[error("CIGAR counts differ from alignment counts")]
    CigarCountMismatch,
    #[error("failed to write CIGAR text")]
    CigarWrite,
    #[error("alignment core is outside the supplied endpoint window")]
    EndpointOutsideWindow,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> AlignmentConfig {
        AlignmentConfig {
            band_width: 8,
            ..AlignmentConfig::default()
        }
    }

    #[test]
    fn affine_alignment_validates_its_cigar() {
        let mut workspace = AlignmentWorkspace::default();
        let alignment = workspace
            .align(b"ACGTACGT", b"ACGTGACGT", config())
            .unwrap();
        assert_eq!(alignment.cigar, "4=1I4=");
        assert_eq!(alignment.insertions, 1);
        assert_eq!(alignment.identity(), 8.0 / 9.0);
        let json = serde_json::to_value(&alignment).unwrap();
        assert_eq!(json["identity"].as_f64().unwrap(), alignment.identity());
        assert_eq!(
            serde_json::from_value::<Alignment>(json).unwrap(),
            alignment
        );
        alignment.validate_cigar().unwrap();
    }

    #[test]
    fn reverse_interval_uses_forward_coordinates() {
        let mut workspace = AlignmentWorkspace::default();
        let alignment = workspace
            .align_oriented(b"AGGACTT", b"AAGTCCT", 100, Strand::Reverse, config())
            .unwrap();
        assert_eq!(alignment.target_interval, Interval::new(100, 107).unwrap());
        assert_eq!(alignment.strand, Strand::Reverse);
    }

    #[test]
    fn complements_every_accepted_iupac_symbol() {
        for (&base, &expected) in b"ACGTURYSWKMBDHVN".iter().zip(b"TGCAAYRSWMKVHDBN") {
            assert_eq!(complement(base), expected);
            assert_eq!(complement(base.to_ascii_lowercase()), expected);
        }
    }

    #[test]
    fn workspace_reuses_matrix_capacity() {
        let mut workspace = AlignmentWorkspace::default();
        workspace.align(b"ACGT", b"ACGT", config()).unwrap();
        let capacity = workspace.capacity_cells();
        workspace.align(b"ACG", b"ACG", config()).unwrap();
        assert_eq!(workspace.capacity_cells(), capacity);
    }

    #[test]
    fn malformed_cigar_is_rejected() {
        assert_eq!(
            parse_cigar("4=1I4").unwrap_err(),
            AlignmentError::InvalidCigar
        );
        assert_eq!(
            parse_cigar("4=1I1I").unwrap_err(),
            AlignmentError::InvalidCigar
        );
    }

    fn core(strand: Strand, target_interval: Interval) -> Alignment {
        Alignment {
            score: 8,
            strand,
            query_interval: Interval::new(2, 6).unwrap(),
            target_interval,
            matches: 4,
            substitutions: 0,
            insertions: 0,
            deletions: 0,
            cigar: "4=".into(),
            edit_script: vec![EditRun {
                operation: EditOperation::Equal,
                length: 4,
            }],
        }
    }

    #[test]
    fn completes_both_endpoints_within_bound() {
        let mut workspace = AlignmentWorkspace::default();
        let completed = workspace
            .complete_endpoints(
                core(Strand::Forward, Interval::new(102, 106).unwrap()),
                b"TTACGTAA",
                b"TTACGTAA",
                100,
                2,
                config(),
            )
            .unwrap();
        assert_eq!(
            completed.alignment.query_interval,
            Interval::new(0, 8).unwrap()
        );
        assert_eq!(
            completed.alignment.target_interval,
            Interval::new(100, 108).unwrap()
        );
        assert_eq!(completed.alignment.cigar, "8=");
        assert_eq!(completed.metrics.left_query_bases, 2);
        assert_eq!(completed.metrics.right_query_bases, 2);
    }

    #[test]
    fn completion_preserves_reverse_coordinates() {
        let mut workspace = AlignmentWorkspace::default();
        let completed = workspace
            .complete_endpoints(
                core(Strand::Reverse, Interval::new(102, 106).unwrap()),
                b"TTACGCAA",
                b"TTGCGTAA",
                100,
                2,
                config(),
            )
            .unwrap();
        assert_eq!(
            completed.alignment.target_interval,
            Interval::new(100, 108).unwrap()
        );
        assert_eq!(completed.alignment.cigar, "8=");
    }

    #[test]
    fn semiglobal_completion_leaves_outer_overhang_free() {
        let mut cells = Vec::new();
        let result =
            anchored_semiglobal(&mut cells, b"ACGTACGT", b"ACGTACGTCCCC", config()).unwrap();
        assert_eq!((result.query_bases, result.target_bases), (8, 8));
        assert_eq!(cigar_from_runs(&result.runs).unwrap(), "8=");
    }

    #[test]
    fn semiglobal_workspace_growth_matches_fresh_workspace() {
        let mut reused = Vec::new();
        anchored_semiglobal(&mut reused, b"A", b"A", config()).unwrap();
        let reused_result = anchored_semiglobal(&mut reused, b"C", b"AAA", config()).unwrap();
        let fresh_result = anchored_semiglobal(&mut Vec::new(), b"C", b"AAA", config()).unwrap();

        assert_eq!(reused_result.runs, fresh_result.runs);
        assert_eq!(reused_result.query_bases, fresh_result.query_bases);
        assert_eq!(reused_result.target_bases, fresh_result.target_bases);
    }
}
