//! Joins of short walks across divergent pieces. A local alignment ends at a conserved block
//! shorter than the fragment minimum when the piece to the next block costs more than one block
//! scores, although the joined path has enough columns. A walk is linked to its nearest walk on
//! each side, ordered on query and contig, within the reach of `gap_continuation`. Only the
//! pieces between linked walks are aligned, globally; the best qualifying run of links through
//! the walk is kept.

use super::anchors::QueryTable;
use super::chain_alignment::BandConfig;
use super::dense::Walk;
use super::walk_alignment::{append, oriented_alignment, walk_alignment};
use super::{Alignment, AlignmentError, EditOperation, EditRun};
use std::mem::size_of;

/// The least columns and identity a join must reach to be kept.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct JoinLimits {
    pub min_columns: u64,
    pub min_identity: f64,
}

/// Columns, matches and score of a walk or of the piece between two linked walks.
#[derive(Clone, Copy, Default)]
struct Part {
    score: i64,
    columns: u64,
    matches: u64,
}

impl std::ops::Add for Part {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            score: self.score + other.score,
            columns: self.columns + other.columns,
            matches: self.matches + other.matches,
        }
    }
}

/// Reusable buffers of the joins of one window's walks.
#[derive(Default)]
pub struct JoinScratch {
    limits: JoinLimits,
    starts: Vec<u32>,
    ends: Vec<u32>,
    chain: Vec<u32>,
    parts: Vec<[Part; 2]>,
    traceback: Vec<u8>,
    scores: Vec<[i64; 3]>,
    piece: Vec<u8>,
}

/// The query and contig bases a link may cross, as for `gap_continuation`.
fn reach(table: &QueryTable, config: BandConfig) -> u64 {
    u64::from(config.max_drift) + u64::from(table.anchors().a)
}

impl JoinScratch {
    /// The bytes `prepare` holds for `walks` walks: walk orders, one chain, and one piece
    /// alignment of at most the reach on each side with its operations.
    pub fn bytes(walks: usize, table: &QueryTable, config: BandConfig) -> usize {
        let reach = reach(table, config) as usize;
        walks * (3 * size_of::<u32>() + size_of::<[Part; 2]>())
            + (reach + 1) * (reach + 1)
            + 2 * (reach + 1) * size_of::<[i64; 3]>()
            + reach
            + 2 * reach * size_of::<EditRun>()
    }

    /// Orders `walks` by strand and contig start and end, at the capacities `bytes` priced.
    pub fn prepare(
        &mut self,
        walks: &[Walk],
        table: &QueryTable,
        config: BandConfig,
        limits: JoinLimits,
    ) {
        let reach = reach(table, config) as usize;
        let key = |index: u32, end: bool| {
            let walk = walks[index as usize];
            let contig = walk.contig + if end { walk.length } else { 0 };
            (walk.reverse, contig, index)
        };
        self.limits = limits;
        self.starts = (0..walks.len() as u32).collect();
        self.starts.sort_unstable_by_key(|&index| key(index, false));
        self.ends = (0..walks.len() as u32).collect();
        self.ends.sort_unstable_by_key(|&index| key(index, true));
        self.chain = Vec::with_capacity(walks.len());
        self.parts = Vec::with_capacity(walks.len());
        self.traceback = Vec::with_capacity((reach + 1) * (reach + 1));
        self.scores = Vec::with_capacity(2 * (reach + 1));
        self.piece = Vec::with_capacity(reach);
    }

    /// The nearest walk after (or before) walk `at` on its strand, ordered on query and contig,
    /// with both pieces between them within the reach and drift, and the query piece length.
    fn link(
        &self,
        walks: &[Walk],
        at: usize,
        after: bool,
        table: &QueryTable,
        config: BandConfig,
    ) -> Option<(usize, u64)> {
        let reach = reach(table, config);
        let length = u64::from(table.query_length());
        let walk = walks[at];
        let (order, low, high) = if after {
            let end = u64::from(walk.contig + walk.length);
            (&self.starts, end, end + reach)
        } else {
            let start = u64::from(walk.contig);
            (&self.ends, start.saturating_sub(reach), start)
        };
        let key = |other: &Walk| {
            let contig = u64::from(other.contig) + if after { 0 } else { u64::from(other.length) };
            (other.reverse, contig)
        };
        let first =
            order.partition_point(|&index| key(&walks[index as usize]) < (walk.reverse, low));
        let mut best: Option<((u64, u64, u32), u64)> = None;
        for &index in &order[first..] {
            let other = walks[index as usize];
            if key(&other) > (walk.reverse, high) {
                break;
            }
            let (before, next) = if after { (walk, other) } else { (other, walk) };
            let contig_gap = u64::from(next.contig) - u64::from(before.contig + before.length);
            let end = i64::from(before.query) + i64::from(before.length);
            let query_gap = i64::from(next.query) - end;
            let query_gap = if table.circular() {
                query_gap.rem_euclid(length as i64)
            } else {
                query_gap
            };
            let Ok(query_gap) = u64::try_from(query_gap) else {
                continue;
            };
            if query_gap > reach || query_gap.abs_diff(contig_gap) > u64::from(config.max_drift) {
                continue;
            }
            let rank = (query_gap + contig_gap, contig_gap, index);
            if best.is_none_or(|(saved, _)| rank < saved) {
                best = Some((rank, query_gap));
            }
        }
        best.map(|((_, _, index), gap)| (index as usize, gap))
    }

    /// The best qualifying join through walk `at`, of at least two walks, from the oriented
    /// window `oriented` that starts `offset` bases into the contig. Its score may be zero or
    /// negative: only the join limits decide.
    pub fn join_through(
        &mut self,
        walks: &[Walk],
        at: usize,
        oriented: &[u8],
        offset: u64,
        table: &QueryTable,
        config: BandConfig,
    ) -> Result<Option<Alignment>, AlignmentError> {
        self.link_chain(walks, at, table, config);
        if self.chain.len() < 2 {
            return Ok(None);
        }
        let middle = self
            .chain
            .iter()
            .position(|&index| index as usize == at)
            .unwrap_or(0);
        self.score_parts(walks, oriented, table, config);
        let Some((score, first, last)) = self.best_run(middle) else {
            return Ok(None);
        };
        let mut runs: Vec<EditRun> = Vec::new();
        for position in first..=last {
            let walk = walks[self.chain[position] as usize];
            let alignment = walk_alignment(&walk, oriented, offset, table, config)?;
            append(&mut runs, &alignment.runs);
            if position < last {
                let next = walks[self.chain[position + 1] as usize];
                self.piece(walk, next, oriented, table, config, Some(&mut runs));
            }
        }
        let first = walks[self.chain[first] as usize];
        oriented_alignment(runs, score, &first, oriented.len(), offset).map(Some)
    }

    /// Links walks to walk `at` on both sides; a circular join covers the query at most once.
    fn link_chain(&mut self, walks: &[Walk], at: usize, table: &QueryTable, config: BandConfig) {
        let length = u64::from(table.query_length());
        self.chain.clear();
        self.chain.push(at as u32);
        let mut span = u64::from(walks[at].length);
        for after in [false, true] {
            let mut current = at;
            while let Some((next, gap)) = self.link(walks, current, after, table, config) {
                span += gap + u64::from(walks[next].length);
                if span > length {
                    break;
                }
                self.chain.push(next as u32);
                current = next;
            }
            if !after {
                self.chain.reverse();
            }
        }
    }

    /// The parts of every linked walk and of the piece after it.
    fn score_parts(
        &mut self,
        walks: &[Walk],
        oriented: &[u8],
        table: &QueryTable,
        config: BandConfig,
    ) {
        let query = table.query();
        self.parts.clear();
        for position in 0..self.chain.len() {
            let walk = walks[self.chain[position] as usize];
            let mut part = Part::default();
            for step in 0..walk.length as usize {
                let same = oriented[walk.contig as usize + step]
                    == query[(walk.query as usize + step) % query.len()];
                part = part
                    + Part {
                        score: i64::from(if same {
                            config.alignment.match_score
                        } else {
                            config.alignment.mismatch_score
                        }),
                        columns: 1,
                        matches: u64::from(same),
                    };
            }
            let piece = match self.chain.get(position + 1) {
                Some(&next) => {
                    self.piece(walk, walks[next as usize], oriented, table, config, None)
                }
                None => Part::default(),
            };
            self.parts.push([part, piece]);
        }
    }

    /// The best run of walks `first..=last` around position `middle` that meets the limits, by
    /// score, then the first found.
    fn best_run(&self, middle: usize) -> Option<(i64, usize, usize)> {
        let limits = self.limits;
        let mut best: Option<(i64, usize, usize)> = None;
        let mut left = Part::default();
        for first in (0..=middle).rev() {
            left = left
                + self.parts[first][0]
                + if first < middle {
                    self.parts[first][1]
                } else {
                    Part::default()
                };
            let mut whole = left;
            for last in middle..self.chain.len() {
                if last > middle {
                    whole = whole + self.parts[last - 1][1] + self.parts[last][0];
                }
                let qualifies = whole.columns >= limits.min_columns
                    && whole.matches as f64 >= limits.min_identity * whole.columns as f64;
                if first < last && qualifies && best.is_none_or(|(score, _, _)| whole.score > score)
                {
                    best = Some((whole.score, first, last));
                }
            }
        }
        best
    }

    /// The global alignment of the piece between linked walks: a gap of n bases scores the gap
    /// open score plus n gap extend scores. Its runs are appended to `runs` when given.
    fn piece(
        &mut self,
        before: Walk,
        next: Walk,
        oriented: &[u8],
        table: &QueryTable,
        config: BandConfig,
        runs: Option<&mut Vec<EditRun>>,
    ) -> Part {
        let query = table.query();
        let end = before.query as usize + before.length as usize;
        let query_gap = (next.query as usize + query.len() - end % query.len()) % query.len();
        self.piece.clear();
        self.piece
            .extend((0..query_gap).map(|step| query[(end + step) % query.len()]));
        let contig = &oriented[(before.contig + before.length) as usize..next.contig as usize];
        let score = self.fill_piece(contig, config);
        self.trace_piece(contig, score, runs)
    }

    /// Fills the global piece matrix of `self.piece` against `contig` and returns the end score
    /// and state. States: 0 paired, 1 contig base only (deletion), 2 query base only (insertion).
    /// Each cell keeps the state before each of its states, two bits each.
    fn fill_piece(&mut self, contig: &[u8], config: BandConfig) -> (i64, u8) {
        let scoring = config.alignment;
        let (rows, columns) = (self.piece.len() + 1, contig.len() + 1);
        let (open, extend) = (
            i64::from(scoring.gap_open_score),
            i64::from(scoring.gap_extend_score),
        );
        const NONE: i64 = i64::MIN / 4;
        let best = |scores: [i64; 3]| {
            (0..3).fold((NONE, 0u8), |saved, state| {
                if scores[state] > saved.0 {
                    (scores[state], state as u8)
                } else {
                    saved
                }
            })
        };
        self.traceback.clear();
        self.traceback.resize(rows * columns, 0);
        self.scores.clear();
        self.scores.resize(2 * columns, [NONE; 3]);
        let (previous, current) = self.scores.split_at_mut(columns);
        previous[0] = [0, NONE, NONE];
        for (column, cell) in previous.iter_mut().enumerate().skip(1) {
            *cell = [NONE, open + extend * column as i64, NONE];
            self.traceback[column] = u8::from(column > 1) << 2;
        }
        for row in 1..rows {
            current[0] = [NONE, NONE, open + extend * row as i64];
            self.traceback[row * columns] = (if row > 1 { 2 } else { 0 }) << 4;
            for column in 1..columns {
                let (paired, from) = best(previous[column - 1]);
                let paired = paired
                    + i64::from(if self.piece[row - 1] == contig[column - 1] {
                        scoring.match_score
                    } else {
                        scoring.mismatch_score
                    });
                let left = current[column - 1];
                let (deleted, deleted_from) = best([
                    left[0] + open + extend,
                    left[1] + extend,
                    left[2] + open + extend,
                ]);
                let up = previous[column];
                let (inserted, inserted_from) =
                    best([up[0] + open + extend, up[1] + open + extend, up[2] + extend]);
                current[column] = [paired, deleted, inserted];
                self.traceback[row * columns + column] =
                    from | deleted_from << 2 | inserted_from << 4;
            }
            previous.copy_from_slice(current);
        }
        if rows == 1 && columns == 1 {
            (0, 0)
        } else {
            best(previous[columns - 1])
        }
    }

    /// Traces the filled piece matrix back from its end cell.
    fn trace_piece(
        &self,
        contig: &[u8],
        (score, mut state): (i64, u8),
        runs: Option<&mut Vec<EditRun>>,
    ) -> Part {
        let columns = contig.len() + 1;
        let mut part = Part {
            score,
            ..Part::default()
        };
        let (mut row, mut column) = (self.piece.len(), contig.len());
        let mut operations = Vec::new();
        while row > 0 || column > 0 {
            let before = self.traceback[row * columns + column] >> (2 * state) & 3;
            let operation = match state {
                0 => {
                    let same = self.piece[row - 1] == contig[column - 1];
                    part.matches += u64::from(same);
                    (row, column) = (row - 1, column - 1);
                    if same {
                        EditOperation::Match
                    } else {
                        EditOperation::Substitution
                    }
                }
                1 => {
                    column -= 1;
                    EditOperation::Deletion
                }
                _ => {
                    row -= 1;
                    EditOperation::Insertion
                }
            };
            part.columns += 1;
            if runs.is_some() {
                operations.push(EditRun {
                    operation,
                    length: 1,
                });
            }
            state = before;
        }
        if let Some(runs) = runs {
            operations.reverse();
            append(runs, &operations);
        }
        part
    }
}
