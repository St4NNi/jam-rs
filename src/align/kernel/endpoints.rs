//! Endpoint completion: extends a local alignment through at most a few bases beyond each end
//! with an anchored semi-global alignment, so a short divergent end is not cut off.

use super::geometry::{DELETION, INSERTION, MATCH, START};
use super::{AlignmentConfig, AlignmentInput, AlignmentWorkspace};
use crate::align::cigar::{append_runs, runs_of};
use crate::align::reverse_complement_into;
use crate::align::{Alignment, AlignmentError, EditOperation, EditRun, Interval, Strand};

/// A score no path reaches; adding scoring terms to it cannot overflow.
const NEGATIVE: i32 = i32::MIN / 4;

/// The previous states of a cell before any is known: two bits per state, the top two zero.
const INITIAL_PREVIOUS: u8 = START | (START << 2) | (START << 4);

/// Reusable buffers of the endpoint matrices.
#[derive(Debug, Default)]
pub(super) struct EndpointWorkspace {
    /// Two rows of the three state scores.
    rows: [Vec<[i32; 3]>; 2],
    /// The previous state of each state of each cell, two bits each.
    previous: Vec<u8>,
    operations: Vec<EditOperation>,
    left_query: Vec<u8>,
    left_contig: Vec<u8>,
}

/// The runs of one end and the bases they take from the query and the contig.
struct EndpointResult {
    runs: Vec<EditRun>,
    query_bases: usize,
    contig_bases: usize,
}

impl AlignmentWorkspace {
    /// Extends `core`, an alignment from `align` over the same `input`, by an anchored
    /// semi-global alignment of at most `max_extension` bases beyond each end. The result is
    /// rescored over all its runs.
    pub fn complete_endpoints(
        &mut self,
        core: Alignment,
        input: AlignmentInput<'_>,
        max_extension: usize,
        config: AlignmentConfig,
    ) -> Result<Alignment, AlignmentError> {
        config.validate()?;
        match core.strand {
            Strand::Forward => complete(
                &mut self.endpoint,
                core,
                input,
                input.contig,
                max_extension,
                config,
            ),
            Strand::Reverse => {
                let mut reverse = std::mem::take(&mut self.reverse);
                reverse_complement_into(input.contig, &mut reverse);
                let result = complete(
                    &mut self.endpoint,
                    core,
                    input,
                    &reverse,
                    max_extension,
                    config,
                );
                self.reverse = reverse;
                result
            }
        }
    }
}

/// `complete_endpoints` with the contig window on the core's strand.
fn complete(
    workspace: &mut EndpointWorkspace,
    mut core: Alignment,
    input: AlignmentInput<'_>,
    contig: &[u8],
    max_extension: usize,
    config: AlignmentConfig,
) -> Result<Alignment, AlignmentError> {
    core.validate()?;
    let query = input.query;
    let local = |position: u64, offset: u64| {
        position
            .checked_sub(offset)
            .and_then(|value| usize::try_from(value).ok())
            .ok_or(AlignmentError::EndpointOutsideWindow)
    };
    let query_start = local(core.query_interval.start, input.query_offset)?;
    let query_end = local(core.query_interval.end, input.query_offset)?;
    let local_start = local(core.contig_interval.start, input.contig_offset)?;
    let local_end = local(core.contig_interval.end, input.contig_offset)?;
    if query_end > query.len() || local_end > contig.len() {
        return Err(AlignmentError::EndpointOutsideWindow);
    }
    let (contig_start, contig_end) = match core.strand {
        Strand::Forward => (local_start, local_end),
        Strand::Reverse => (contig.len() - local_end, contig.len() - local_start),
    };

    let left_query_len = query_start.min(max_extension);
    let left_contig_len = contig_start.min(max_extension);
    let mut left_query = std::mem::take(&mut workspace.left_query);
    let mut left_contig = std::mem::take(&mut workspace.left_contig);
    left_query.clear();
    left_query.extend(
        query[query_start - left_query_len..query_start]
            .iter()
            .rev(),
    );
    left_contig.clear();
    left_contig.extend(
        contig[contig_start - left_contig_len..contig_start]
            .iter()
            .rev(),
    );
    let left = anchored_semiglobal(workspace, &left_query, &left_contig, config);
    (workspace.left_query, workspace.left_contig) = (left_query, left_contig);
    let mut left = left?;
    left.runs.reverse();

    let query_limit = query_end.saturating_add(max_extension).min(query.len());
    let contig_limit = contig_end.saturating_add(max_extension).min(contig.len());
    let right = anchored_semiglobal(
        workspace,
        &query[query_end..query_limit],
        &contig[contig_end..contig_limit],
        config,
    )?;

    let mut runs = Vec::with_capacity(left.runs.len() + core.runs.len() + right.runs.len());
    append_runs(&mut runs, left.runs)?;
    append_runs(&mut runs, std::mem::take(&mut core.runs))?;
    append_runs(&mut runs, right.runs)?;
    core.score = score_runs(&runs, config);
    let position = |offset: u64, at: usize| {
        u64::try_from(at)
            .ok()
            .and_then(|at| offset.checked_add(at))
            .ok_or(AlignmentError::LengthOverflow)
    };
    core.query_interval = Interval::new(
        position(input.query_offset, query_start - left.query_bases)?,
        position(input.query_offset, query_end + right.query_bases)?,
    )?;
    let oriented_start = contig_start - left.contig_bases;
    let oriented_end = contig_end + right.contig_bases;
    let (forward_start, forward_end) = match core.strand {
        Strand::Forward => (oriented_start, oriented_end),
        Strand::Reverse => (contig.len() - oriented_end, contig.len() - oriented_start),
    };
    core.contig_interval = Interval::new(
        position(input.contig_offset, forward_start)?,
        position(input.contig_offset, forward_end)?,
    )?;
    core.runs = runs;
    core.validate()?;
    Ok(core)
}

/// The best end cell of an endpoint matrix so far.
#[derive(Clone, Copy)]
struct EndCell {
    score: i32,
    row: usize,
    column: usize,
    state: u8,
}

impl EndCell {
    /// Keeps the better cell: higher score, then the later (row + column, row, column) and the
    /// lower state. The comparison is a total order, so the visiting order does not matter.
    fn consider(&mut self, row: usize, column: usize, scores: [i32; 3]) {
        let (score, state) = maximum(scores);
        let order = |cell: &EndCell| {
            let (row, column) = (cell.row, cell.column);
            (
                row.saturating_add(column),
                row,
                column,
                std::cmp::Reverse(cell.state),
            )
        };
        let candidate = EndCell {
            score,
            row,
            column,
            state,
        };
        if score > self.score || (score == self.score && order(&candidate) > order(self)) {
            *self = candidate;
        }
    }
}

/// The score of a gap of `bases` bases.
fn gap_score(config: AlignmentConfig, bases: usize) -> i32 {
    let bases = i32::try_from(bases).unwrap_or(i32::MAX);
    let extend = config.gap_extend_score.saturating_mul(bases);
    config.gap_open_score.saturating_add(extend)
}

/// The highest score and its state; ties go to the lower state.
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

/// The best alignment that starts at the first base of both `query` and `contig` and ends at
/// the last base of one of them.
fn anchored_semiglobal(
    workspace: &mut EndpointWorkspace,
    query: &[u8],
    contig: &[u8],
    config: AlignmentConfig,
) -> Result<EndpointResult, AlignmentError> {
    if query.is_empty() || contig.is_empty() {
        return Ok(EndpointResult {
            runs: Vec::new(),
            query_bases: 0,
            contig_bases: 0,
        });
    }
    let (rows, columns) = (query.len() + 1, contig.len() + 1);
    let matrix_cells = rows
        .checked_mul(columns)
        .ok_or(AlignmentError::LengthOverflow)?;
    if matrix_cells > config.max_cells {
        return Err(AlignmentError::MatrixTooLarge {
            cells: matrix_cells,
            max_cells: config.max_cells,
        });
    }
    workspace.previous.clear();
    workspace.previous.resize(matrix_cells, INITIAL_PREVIOUS);
    for row in &mut workspace.rows {
        row.resize(columns, [NEGATIVE; 3]);
    }
    workspace.rows[0][0] = [0, NEGATIVE, NEGATIVE];
    for column in 1..columns {
        workspace.rows[0][column] = [NEGATIVE; 3];
        workspace.rows[0][column][DELETION as usize] = gap_score(config, column);
        let state = if column == 1 { MATCH } else { DELETION };
        workspace.previous[column] = (INITIAL_PREVIOUS & !(3 << 2)) | (state << 2);
    }
    let mut best = EndCell {
        score: NEGATIVE,
        row: 0,
        column: 0,
        state: START,
    };
    best.consider(0, columns - 1, workspace.rows[0][columns - 1]);
    for row in 1..rows {
        fill_endpoint_row(workspace, query, contig, config, row);
        let current = &workspace.rows[1][..columns];
        if row + 1 == rows {
            for (column, &scores) in current.iter().enumerate() {
                best.consider(row, column, scores);
            }
        } else {
            best.consider(row, columns - 1, current[columns - 1]);
        }
        workspace.rows.swap(0, 1);
    }
    if best.state > INSERTION {
        return Err(AlignmentError::NoAlignment);
    }
    trace_endpoint(workspace, query, contig, best)
}

/// Fills row `row` of the endpoint matrix into `rows[1]` from `rows[0]`.
fn fill_endpoint_row(
    workspace: &mut EndpointWorkspace,
    query: &[u8],
    contig: &[u8],
    config: AlignmentConfig,
    row: usize,
) {
    let columns = contig.len() + 1;
    let gap_open = config.gap_open();
    let first = row * columns;
    let [previous, current] = &mut workspace.rows;
    let previous = &previous[..columns];
    let current = &mut current[..columns];
    let traceback = &mut workspace.previous[first..first + columns];
    current[0] = [NEGATIVE; 3];
    current[0][INSERTION as usize] = gap_score(config, row);
    let state = if row == 1 { MATCH } else { INSERTION };
    traceback[0] = (INITIAL_PREVIOUS & !(3 << 4)) | (state << 4);
    for column in 1..columns {
        let (score, state) = maximum(previous[column - 1]);
        let left = current[column - 1];
        let cell = &mut current[column];
        cell[MATCH as usize] = score.saturating_add(
            if query[row - 1].eq_ignore_ascii_case(&contig[column - 1]) {
                config.match_score
            } else {
                config.mismatch_score
            },
        );
        let mut packed = state;

        let above = previous[column];
        let (score, state) = maximum([
            above[MATCH as usize].saturating_add(gap_open),
            above[DELETION as usize].saturating_add(gap_open),
            above[INSERTION as usize].saturating_add(config.gap_extend_score),
        ]);
        cell[INSERTION as usize] = score;
        packed |= state << 4;

        let (score, state) = maximum([
            left[MATCH as usize].saturating_add(gap_open),
            left[DELETION as usize].saturating_add(config.gap_extend_score),
            left[INSERTION as usize].saturating_add(gap_open),
        ]);
        cell[DELETION as usize] = score;
        traceback[column] = packed | (state << 2);
    }
}

/// The runs from the start to the best endpoint cell.
fn trace_endpoint(
    workspace: &mut EndpointWorkspace,
    query: &[u8],
    contig: &[u8],
    best: EndCell,
) -> Result<EndpointResult, AlignmentError> {
    let columns = contig.len() + 1;
    let (mut row, mut column, mut state) = (best.row, best.column, best.state);
    let operations = &mut workspace.operations;
    operations.clear();
    while row > 0 || column > 0 {
        let packed = workspace.previous[row * columns + column];
        let previous = (packed >> (2 * state)) & 3;
        match state {
            MATCH if row > 0 && column > 0 => {
                operations.push(
                    if query[row - 1].eq_ignore_ascii_case(&contig[column - 1]) {
                        EditOperation::Match
                    } else {
                        EditOperation::Substitution
                    },
                );
                row -= 1;
                column -= 1;
            }
            DELETION if column > 0 => {
                operations.push(EditOperation::Deletion);
                column -= 1;
            }
            INSERTION if row > 0 => {
                operations.push(EditOperation::Insertion);
                row -= 1;
            }
            _ => return Err(AlignmentError::InvalidTraceback),
        }
        state = previous;
    }
    operations.reverse();
    Ok(EndpointResult {
        runs: runs_of(operations)?,
        query_bases: best.row,
        contig_bases: best.column,
    })
}

/// The affine score of `runs`.
fn score_runs(runs: &[EditRun], config: AlignmentConfig) -> i32 {
    runs.iter().fold(0, |score, run| {
        let length = i32::try_from(run.length).unwrap_or(i32::MAX);
        let delta = match run.operation {
            EditOperation::Match => config.match_score.saturating_mul(length),
            EditOperation::Substitution => config.mismatch_score.saturating_mul(length),
            EditOperation::Insertion | EditOperation::Deletion => config
                .gap_open_score
                .saturating_add(config.gap_extend_score.saturating_mul(length)),
        };
        score.saturating_add(delta)
    })
}
