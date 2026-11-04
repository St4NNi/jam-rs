//! Aligns the stretch around a walk again with gaps when its ends ask for it.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::dense::Walk;
use crate::align::BandConfig;
use crate::align::anchors::QueryTable;
use crate::align::kernel::AlignmentInput;
use crate::align::{
    Alignment, AlignmentConfig, AlignmentError, AlignmentWorkspace, EditOperation, Interval, Strand,
};
use rustc_hash::{FxHashMap, FxHasher};
use std::hash::{Hash, Hasher};
use std::ops::Range;

/// A walk alignment proven to be the best for its query range and contig range, so the next
/// walk with the same inputs and its diagonals inside its band need not be aligned again. Only
/// alignments a walk returns are saved, never those of a range it widens.
#[derive(Clone, Copy, Debug)]
struct DenseOptimum {
    strand: Strand,
    /// The query range `start..end` that was aligned.
    query: (u64, u64),
    /// The oriented window range `low..high` that was aligned.
    contig: (u64, u64),
    /// The lowest and highest diagonal of the alignment path.
    diagonals: (i64, i64),
}

/// The kernel results of the dense walks of one window, kept so later walks can use them.
#[derive(Debug, Default)]
pub struct DenseResults {
    optimum: Option<DenseOptimum>,
    /// Aligned cell boxes and their best alignment; `None` when no cell scores positive.
    kernels: Vec<(KernelBox, Option<Alignment>)>,
    /// The boxes of `kernels` by a hash of their kernel input, from the `SAME_FROM`th box on.
    same: FxHashMap<u64, Vec<usize>>,
    /// Kernel results computed ahead on other threads, by their input; one is taken in place of
    /// running the kernel on that input.
    ahead: FxHashMap<KernelInput, Option<Alignment>>,
    /// Whether the kernel results made are kept in `made`, for another loop to take.
    keep_made: bool,
    made: Vec<(KernelInput, Option<Alignment>)>,
}

/// The input of one kernel call in a window: strand, scores and band, and the query and contig
/// positions and lengths, which fix the bases.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct KernelInput {
    reverse: bool,
    scores: AlignmentConfig,
    query: (u64, usize),
    contig: (u64, usize),
}

impl DenseResults {
    /// Results that keep every kernel result they make, to be taken with `take_made`.
    pub fn keeping_made() -> Self {
        Self {
            keep_made: true,
            ..Self::default()
        }
    }

    /// The kernel results made since the last call, by their input.
    pub fn take_made(&mut self) -> Vec<(KernelInput, Option<Alignment>)> {
        std::mem::take(&mut self.made)
    }

    /// Sets the kernel results computed ahead; those of the last call not taken are dropped.
    pub fn set_ahead(&mut self, made: Vec<(KernelInput, Option<Alignment>)>) {
        self.ahead.clear();
        self.ahead.extend(made);
    }

    /// The kernel result of an input: computed ahead, or made now.
    fn kernel(
        &mut self,
        workspace: &mut AlignmentWorkspace,
        aligned: AlignmentInput<'_>,
        strand: Strand,
        scores: AlignmentConfig,
    ) -> Result<Option<Alignment>, AlignmentError> {
        let input = (self.keep_made || !self.ahead.is_empty()).then(|| KernelInput {
            reverse: strand.is_reverse(),
            scores,
            query: (aligned.query_offset, aligned.query.len()),
            contig: (aligned.contig_offset, aligned.contig.len()),
        });
        if let Some(found) = input.and_then(|input| self.ahead.remove(&input)) {
            return Ok(found);
        }
        let found = kernel(workspace, aligned, strand, scores)?;
        if let (true, Some(input)) = (self.keep_made, input) {
            self.made.push((input, found.clone()));
        }
        Ok(found)
    }
}

/// Boxes a window aligns before boxes with the same kernel input are looked up. A dense repeat
/// aligns the same bases at many positions; other windows align few boxes.
const SAME_FROM: usize = 256;

/// Turns the lookup of boxes with the same kernel input and of the fragments that hold a short
/// walk off, so tests can compare.
#[cfg(test)]
pub(crate) static SAME_OFF: std::sync::atomic::AtomicBool =
    std::sync::atomic::AtomicBool::new(false);

/// Kernel results taken from an earlier box with the same input, for tests.
#[cfg(test)]
pub(crate) static SAME_FOUND: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// One walk of an alignment window: the walks after it, the forward window, the contig position
/// where the window starts, and the query.
#[derive(Clone, Copy)]
pub struct DenseWalk<'a> {
    /// The walk to align.
    pub walk: &'a Walk,
    /// The walks of the window after `walk`; the near ones share its kernel call.
    pub later: &'a [Walk],
    /// The forward window bases.
    pub window: &'a [u8],
    /// The contig position of the first window base.
    pub offset: u64,
    /// The query.
    pub table: &'a QueryTable,
}

/// The gapped alignment around a walk, widened on a linear query while it runs into an edge.
/// `None` when no positive alignment exists, or when `results` already holds the best alignment
/// of the same inputs.
pub fn align_dense_walk(
    input: DenseWalk<'_>,
    config: BandConfig,
    workspace: &mut AlignmentWorkspace,
    results: &mut DenseResults,
) -> Result<Option<Alignment>, AlignmentError> {
    let table = input.table;
    let length = u64::from(table.query_length());
    let margin = config.drift_margin();
    let mut extension = 2 * margin;
    loop {
        let range = query_range(input.walk, table, extension);
        let Some((alignment, proof)) =
            align_query_range(input, extension, config, workspace, results)?
        else {
            return Ok(None);
        };
        let (start, end) = (range.start, range.end);
        let runs_into_edge = !table.circular() && {
            let aligned_start = alignment.query_interval.start - start;
            let aligned_end = alignment.query_interval.end - start;
            (start > 0 && aligned_start <= margin)
                || (end < length && end - start - aligned_end <= margin)
        };
        if !runs_into_edge {
            results.optimum = proof.or(results.optimum);
            return Ok(Some(alignment));
        }
        extension = extension.saturating_mul(4);
    }
}

/// The query positions aligned around `walk`: `extension` bases on each side on a linear query,
/// the whole query from the walk on (running on over the origin) on a circular one.
fn query_range(walk: &Walk, table: &QueryTable, extension: u64) -> Range<u64> {
    let (start, length) = (u64::from(walk.query), u64::from(table.query_length()));
    if table.circular() {
        return start..start + length;
    }
    start.saturating_sub(extension)..(start + u64::from(walk.length) + extension).min(length)
}

/// Aligns a dense walk over the query positions of `extension`, with the proof that the result is
/// the best of its inputs when there is one.
fn align_query_range(
    input: DenseWalk<'_>,
    extension: u64,
    config: BandConfig,
    workspace: &mut AlignmentWorkspace,
    results: &mut DenseResults,
) -> Result<Option<(Alignment, Option<DenseOptimum>)>, AlignmentError> {
    let DenseWalk {
        walk,
        window,
        offset,
        table,
        ..
    } = input;
    let range = query_range(walk, table, extension);
    let (query_start, span) = (range.start, range.end - range.start);
    let Some(cells) = KernelBox::around(walk, range.clone(), config, table, window.len() as u64)
    else {
        return Ok(None);
    };
    let diagonal = walk_diagonal(walk) + query_start as i64;
    let (low, high) = (cells.columns.start, cells.columns.end);
    let band = i64::from(config.max_band());
    if results.optimum.as_ref().is_some_and(|saved| {
        saved.strand == walk.strand
            && saved.query == (range.start, range.end)
            && saved.contig == (low, high)
            && saved.diagonals.0 >= diagonal - band
            && saved.diagonals.1 <= diagonal + band
    }) {
        return Ok(None);
    }
    let mut bases = Vec::new();
    table.copy_bases(
        cells.rows.start,
        cells.rows.end - cells.rows.start,
        &mut bases,
    )?;
    let (aligned, scores) = cells.kernel_input(&bases, window, offset, config);
    let found = match shared_result(input, &cells, extension, config, workspace, results) {
        Some(found) => found,
        None => box_kernel(input, &cells, &bases, config, workspace, results)?,
    };
    let Some(alignment) = found else {
        return Ok(None);
    };
    let kernel_score = alignment.score;
    // A completion reaches a rows beyond the alignment, which the window keeps.
    let a = u64::from(table.anchors().a);
    let alignment = if alignment.query_interval.start - query_start <= a
        && span - (alignment.query_interval.end - query_start) <= a
    {
        workspace.complete_endpoints(alignment, aligned, a as usize, scores)?
    } else {
        alignment
    };
    let diagonals = path_diagonals(&alignment, window.len() as u64, offset, query_start);
    let outside = outside_pairs(
        diagonal - low as i64,
        (high - low) as i64,
        span as i64,
        band,
    );
    let outside = outside * i64::from(config.alignment.match_score);
    let proven = alignment.score >= kernel_score && i64::from(alignment.score) > outside;
    let proof = proven.then_some(DenseOptimum {
        strand: walk.strand,
        query: (range.start, range.end),
        contig: (low, high),
        diagonals,
    });
    Ok(Some((alignment, proof)))
}

/// The best local alignment of a kernel input, `None` when no cell scores positive.
fn kernel(
    workspace: &mut AlignmentWorkspace,
    aligned: AlignmentInput<'_>,
    strand: Strand,
    scores: AlignmentConfig,
) -> Result<Option<Alignment>, AlignmentError> {
    match workspace.align(aligned, strand, scores) {
        Ok(alignment) => Ok(Some(alignment)),
        Err(AlignmentError::NoAlignment) => Ok(None),
        Err(error) => Err(error),
    }
}

/// The kernel result of `cells` taken from a larger box, `None` when it must be aligned alone.
/// When no aligned box holds `cells`, one box over the near later walks is aligned for all.
fn shared_result(
    input: DenseWalk<'_>,
    cells: &KernelBox,
    extension: u64,
    config: BandConfig,
    workspace: &mut AlignmentWorkspace,
    results: &mut DenseResults,
) -> Option<Option<Alignment>> {
    let (window_length, offset) = (input.window.len() as u64, input.offset);
    // The kernel's scores only grow with more cells, and it breaks ties by fixed positions and
    // state order. So the best path of a larger box that lies inside `cells` is also the path
    // the kernel finds in `cells`, and a box without a positive cell has none in its parts.
    let fits = |found: &Option<Alignment>| {
        found
            .as_ref()
            .is_none_or(|alignment| cells.holds_path(alignment, window_length, offset))
    };
    let mut held = results
        .kernels
        .iter()
        .filter(|(larger, _)| larger.holds(cells))
        .peekable();
    if held.peek().is_some() {
        return held
            .find(|(_, found)| fits(found))
            .map(|(_, found)| found.clone());
    }
    let union = near_walks_box(input, cells, extension, config);
    if union == *cells {
        return None;
    }
    let mut bases = Vec::new();
    let span = union.rows.end - union.rows.start;
    // A box the kernel cannot align is left out: the walk's own call reports its errors.
    input
        .table
        .copy_bases(union.rows.start, span, &mut bases)
        .ok()?;
    let found = box_kernel(input, &union, &bases, config, workspace, results).ok()?;
    fits(&found).then_some(found)
}

/// The kernel result of `cells`, whose query rows hold `bases`, added to `results`. In a window
/// that aligned many boxes, the result of an earlier box with the same bases and scores is moved
/// to the positions of `cells`: the kernel reads only the bases and scores, and adds the
/// positions where its inputs start to every result.
fn box_kernel(
    input: DenseWalk<'_>,
    cells: &KernelBox,
    bases: &[u8],
    config: BandConfig,
    workspace: &mut AlignmentWorkspace,
    results: &mut DenseResults,
) -> Result<Option<Alignment>, AlignmentError> {
    let (window, offset) = (input.window, input.offset);
    let (aligned, scores) = cells.kernel_input(bases, window, offset, config);
    #[cfg(test)]
    let off = SAME_OFF.load(std::sync::atomic::Ordering::Relaxed);
    #[cfg(not(test))]
    let off = false;
    let key = (results.kernels.len() >= SAME_FROM && !off).then(|| {
        let mut hasher = FxHasher::default();
        (
            cells.strand.is_reverse(),
            scores.band_width,
            scores.diagonal_offset,
        )
            .hash(&mut hasher);
        (aligned.query, aligned.contig).hash(&mut hasher);
        hasher.finish()
    });
    let earlier = key.and_then(|key| {
        let mut other_bases = Vec::new();
        results.same.get(&key)?.iter().find_map(|&at| {
            let (other, found) = &results.kernels[at];
            if other.strand != cells.strand
                || other.rows.end - other.rows.start != cells.rows.end - cells.rows.start
                || other.columns.end - other.columns.start
                    != cells.columns.end - cells.columns.start
            {
                return None;
            }
            other_bases.clear();
            let rows = other.rows.end - other.rows.start;
            input
                .table
                .copy_bases(other.rows.start, rows, &mut other_bases)
                .ok()?;
            let (then, then_scores) = other.kernel_input(&other_bases, window, offset, config);
            if then_scores != scores || then.query != aligned.query || then.contig != aligned.contig
            {
                return None;
            }
            moved(found, then, aligned)
        })
    });
    #[cfg(test)]
    if earlier.is_some() {
        SAME_FOUND.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    }
    let found = match earlier {
        Some(found) => found,
        None => {
            let found = results.kernel(workspace, aligned, cells.strand, scores)?;
            if let Some(key) = key {
                results
                    .same
                    .entry(key)
                    .or_default()
                    .push(results.kernels.len());
            }
            found
        }
    };
    results.kernels.push((cells.clone(), found.clone()));
    Ok(found)
}

/// The kernel result `found` of the input `then`, moved to the positions of the same bases in
/// `now`. `None` when a position would not fit.
fn moved(
    found: &Option<Alignment>,
    then: AlignmentInput<'_>,
    now: AlignmentInput<'_>,
) -> Option<Option<Alignment>> {
    let Some(alignment) = found else {
        return Some(None);
    };
    let shift = |interval: Interval, from: u64, to: u64| {
        let start = (interval.start - from).checked_add(to)?;
        let end = (interval.end - from).checked_add(to)?;
        Interval::new(start, end).ok()
    };
    Some(Some(Alignment {
        query_interval: shift(
            alignment.query_interval,
            then.query_offset,
            now.query_offset,
        )?,
        contig_interval: shift(
            alignment.contig_interval,
            then.contig_offset,
            now.contig_offset,
        )?,
        runs: alignment.runs.clone(),
        ..*alignment
    }))
}

/// The smallest box that holds `cells` and the boxes of the later walks near its walk (on the
/// same strand, the diagonal within one band, the start within one margin), up to four times
/// the cells of `cells`.
fn near_walks_box(
    input: DenseWalk<'_>,
    cells: &KernelBox,
    extension: u64,
    config: BandConfig,
) -> KernelBox {
    let DenseWalk {
        walk, later, table, ..
    } = input;
    let (band, margin) = (i64::from(config.max_band()), config.drift_margin() as i64);
    let (diagonal, query) = (walk_diagonal(walk), i64::from(walk.query));
    let limit = cells.cells().saturating_mul(4);
    let window_length = input.window.len() as u64;
    let mut union = cells.clone();
    let near = later.iter().take_while(|other| {
        other.strand == walk.strand && (table.circular() || walk_diagonal(other) <= diagonal + band)
    });
    for other in near {
        if (walk_diagonal(other) - diagonal).abs() > band
            || (i64::from(other.query) - query).abs() > margin
        {
            continue;
        }
        let range = query_range(other, table, extension);
        let Some(other) = KernelBox::around(other, range, config, table, window_length) else {
            continue;
        };
        let wider = union.union(&other);
        if wider.cells() <= limit {
            union = wider;
        }
    }
    union
}

/// Oriented window position minus query position at the start of a walk.
fn walk_diagonal(walk: &Walk) -> i64 {
    i64::from(walk.contig) - i64::from(walk.query)
}

/// The cells of one kernel call: query rows (counted on over the origin of a circular query),
/// oriented window columns, and the band, the diagonals (column minus row) from `diagonals.0`
/// to `diagonals.1`.
#[derive(Clone, Debug, Eq, PartialEq)]
struct KernelBox {
    strand: Strand,
    rows: Range<u64>,
    columns: Range<u64>,
    diagonals: (i64, i64),
}

impl KernelBox {
    /// The cells aligned around `walk` over the query positions `range`, `None` when the window
    /// holds no column for them.
    fn around(
        walk: &Walk,
        range: Range<u64>,
        config: BandConfig,
        table: &QueryTable,
        window_length: u64,
    ) -> Option<Self> {
        let (query_start, span) = (range.start, range.end - range.start);
        let margin = config.drift_margin() as i64;
        let diagonal = walk_diagonal(walk) + query_start as i64;
        let low = (diagonal - margin).max(0) as u64;
        let high = ((diagonal + span as i64 + margin).max(0) as u64).min(window_length);
        if low >= high {
            return None;
        }
        // Rows before the first and after the last that can meet the contig within the band,
        // and an endpoint completion's reach, hold no cells: the kernel reads only the rows
        // between.
        let band = i64::from(config.max_band());
        let reach = band + i64::from(table.anchors().a);
        let first_row = (low as i64 - diagonal - reach).clamp(0, span as i64) as u64;
        let last_row = (high as i64 - diagonal + reach).clamp(first_row as i64, span as i64);
        Some(Self {
            strand: walk.strand,
            rows: query_start + first_row..query_start + last_row as u64,
            columns: low..high,
            diagonals: (walk_diagonal(walk) - band, walk_diagonal(walk) + band),
        })
    }

    /// The smallest box that holds both boxes, on the strand of `self`.
    fn union(&self, other: &Self) -> Self {
        Self {
            strand: self.strand,
            rows: self.rows.start.min(other.rows.start)..self.rows.end.max(other.rows.end),
            columns: self.columns.start.min(other.columns.start)
                ..self.columns.end.max(other.columns.end),
            diagonals: (
                self.diagonals.0.min(other.diagonals.0),
                self.diagonals.1.max(other.diagonals.1),
            ),
        }
    }

    /// An upper bound of the cells: rows times band diagonals.
    fn cells(&self) -> u64 {
        let width = (self.diagonals.1 - self.diagonals.0 + 1) as u64;
        (self.rows.end - self.rows.start).saturating_mul(width)
    }

    /// Whether every cell of `other` lies in this box.
    fn holds(&self, other: &Self) -> bool {
        self.strand == other.strand
            && self.rows.start <= other.rows.start
            && other.rows.end <= self.rows.end
            && self.columns.start <= other.columns.start
            && other.columns.end <= self.columns.end
            && self.diagonals.0 <= other.diagonals.0
            && other.diagonals.1 <= self.diagonals.1
    }

    /// Whether every cell of the path of `alignment` (of the window of `window_length` bases
    /// from contig position `offset`) lies in this box.
    fn holds_path(&self, alignment: &Alignment, window_length: u64, offset: u64) -> bool {
        let (start, end) = (
            alignment.contig_interval.start - offset,
            alignment.contig_interval.end - offset,
        );
        let columns = if alignment.strand.is_reverse() {
            window_length - end..window_length - start
        } else {
            start..end
        };
        let diagonals = path_diagonals(alignment, window_length, offset, 0);
        self.strand == alignment.strand
            && self.rows.start <= alignment.query_interval.start
            && alignment.query_interval.end <= self.rows.end
            && self.columns.start <= columns.start
            && columns.end <= self.columns.end
            && self.diagonals.0 <= diagonals.0
            && diagonals.1 <= self.diagonals.1
    }

    /// The kernel input and scores of this box, with `bases` the query bases of its rows.
    fn kernel_input<'a>(
        &self,
        bases: &'a [u8],
        window: &'a [u8],
        offset: u64,
        config: BandConfig,
    ) -> (AlignmentInput<'a>, AlignmentConfig) {
        let length = window.len() as u64;
        let (low, high) = if self.strand.is_reverse() {
            (length - self.columns.end, length - self.columns.start)
        } else {
            (self.columns.start, self.columns.end)
        };
        // The narrowest band around a whole diagonal that reaches both outer diagonals.
        let width = (self.diagonals.1 - self.diagonals.0 + 1) / 2;
        let centre = self.diagonals.0 + width;
        let scores = AlignmentConfig {
            band_width: width as u32,
            diagonal_offset: centre - self.columns.start as i64 + self.rows.start as i64,
            ..config.alignment
        };
        let aligned = AlignmentInput {
            query: bases,
            query_offset: self.rows.start,
            contig: &window[low as usize..high as usize],
            contig_offset: offset + low,
        };
        (aligned, scores)
    }
}

/// The most paired bases of any diagonal outside the band of width `band` around `centre`, in a
/// matrix of `rows` query and `columns` contig bases. Every other operation scores at most zero,
/// so an alignment that scores more than these bases of matches is the best for its inputs.
fn outside_pairs(centre: i64, columns: i64, rows: i64, band: i64) -> i64 {
    let left = (rows - (band + 1 - centre).max(0)).max(0).min(columns);
    let right = (columns - (centre + band + 1).max(0)).max(0).min(rows);
    left.max(right)
}

/// The lowest and highest diagonal (oriented window position minus query position from
/// `query_start`) that the path of `alignment` visits.
fn path_diagonals(
    alignment: &Alignment,
    window_length: u64,
    offset: u64,
    query_start: u64,
) -> (i64, i64) {
    let mut contig_at = if alignment.strand.is_reverse() {
        window_length - (alignment.contig_interval.end - offset)
    } else {
        alignment.contig_interval.start - offset
    };
    let mut query_at = alignment.query_interval.start - query_start;
    let first = contig_at as i64 - query_at as i64;
    let mut diagonals = (first, first);
    for run in &alignment.runs {
        let length = u64::from(run.length);
        match run.operation {
            EditOperation::Deletion => contig_at += length,
            EditOperation::Insertion => query_at += length,
            EditOperation::Match | EditOperation::Substitution => {
                query_at += length;
                contig_at += length;
            }
        }
        let diagonal = contig_at as i64 - query_at as i64;
        diagonals = (diagonals.0.min(diagonal), diagonals.1.max(diagonal));
    }
    diagonals
}
