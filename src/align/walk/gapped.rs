//! Aligns the stretch around a walk again with gaps when its ends ask for it.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::dense::Walk;
use crate::align::BandConfig;
use crate::align::anchors::QueryTable;
use crate::align::{
    Alignment, AlignmentConfig, AlignmentError, AlignmentInput, AlignmentWorkspace, EditOperation,
    Strand,
};

/// A walk alignment proven to be the best for its query window and contig range, so the next
/// walk with the same inputs and its diagonals inside its band need not be aligned again.
#[derive(Clone, Copy, Debug)]
pub struct DenseOptimum {
    strand: Strand,
    query_start: u64,
    contig: (u64, u64),
    diagonals: (i64, i64),
}

/// One walk of an alignment window: the forward window, the contig position where it starts,
/// and the query.
#[derive(Clone, Copy)]
struct WalkInput<'a> {
    walk: &'a Walk,
    window: &'a [u8],
    offset: u64,
    table: &'a QueryTable,
}

/// The gapped alignment around `walk` of `window` (the forward contig bases from `offset`),
/// widened on a linear query while it runs into an edge. `None` when no positive alignment
/// exists, or when `optimum` already holds the best alignment of the same inputs.
pub fn align_dense_walk(
    walk: &Walk,
    window: &[u8],
    offset: u64,
    table: &QueryTable,
    config: BandConfig,
    workspace: &mut AlignmentWorkspace,
    optimum: &mut Option<DenseOptimum>,
) -> Result<Option<Alignment>, AlignmentError> {
    let input = WalkInput {
        walk,
        window,
        offset,
        table,
    };
    let length = u64::from(table.query_length());
    if table.circular() {
        let start = u64::from(walk.query);
        return align_query_range(input, start..start + length, config, workspace, optimum);
    }
    let margin = config.drift_margin();
    let walk_end = u64::from(walk.query) + u64::from(walk.length);
    let mut extension = 2 * margin;
    loop {
        let start = u64::from(walk.query).saturating_sub(extension);
        let end = (walk_end + extension).min(length);
        let Some(alignment) = align_query_range(input, start..end, config, workspace, optimum)?
        else {
            return Ok(None);
        };
        let (aligned_start, aligned_end) = (
            alignment.query_interval.start - start,
            alignment.query_interval.end - start,
        );
        let runs_into_edge = (start > 0 && aligned_start <= margin)
            || (end < length && end - start - aligned_end <= margin);
        if !runs_into_edge {
            return Ok(Some(alignment));
        }
        extension = extension.saturating_mul(4);
    }
}

/// Aligns a dense walk over the query positions `range`, which run on over the origin of a
/// circular query.
fn align_query_range(
    input: WalkInput<'_>,
    range: std::ops::Range<u64>,
    config: BandConfig,
    workspace: &mut AlignmentWorkspace,
    optimum: &mut Option<DenseOptimum>,
) -> Result<Option<Alignment>, AlignmentError> {
    let (query_start, span) = (range.start, range.end - range.start);
    let WalkInput {
        walk,
        window,
        offset,
        table,
    } = input;
    let margin = config.drift_margin();
    let diagonal = i64::from(walk.contig) - (u64::from(walk.query) - query_start) as i64;
    let low = (diagonal - margin as i64).max(0) as u64;
    let high = (diagonal + span as i64 + margin as i64).max(0) as u64;
    let high = high.min(window.len() as u64);
    if low >= high {
        return Ok(None);
    }
    let band = i64::from(config.max_band());
    if optimum.as_ref().is_some_and(|saved| {
        saved.strand == walk.strand
            && saved.query_start == query_start
            && saved.contig == (low, high)
            && saved.diagonals.0 >= diagonal - band
            && saved.diagonals.1 <= diagonal + band
    }) {
        return Ok(None);
    }
    let (forward_low, forward_high) = if walk.strand.is_reverse() {
        (window.len() as u64 - high, window.len() as u64 - low)
    } else {
        (low, high)
    };
    // Rows before the first and after the last that can meet the contig within the band, and
    // an endpoint completion's reach, hold no cells: the kernel reads only the rows between.
    let a = u64::from(table.anchors().a);
    let reach = band + a as i64;
    let first_row = (low as i64 - diagonal - reach).clamp(0, span as i64) as u64;
    let last_row = (high as i64 - diagonal + reach).clamp(first_row as i64, span as i64);
    let rows = last_row as u64 - first_row;
    let mut bases = Vec::new();
    table.copy_bases(query_start + first_row, rows, &mut bases)?;
    let scores = AlignmentConfig {
        band_width: config.max_band(),
        diagonal_offset: diagonal - low as i64 + first_row as i64,
        ..config.alignment
    };
    let aligned = AlignmentInput {
        query: &bases,
        query_offset: query_start + first_row,
        contig: &window[forward_low as usize..forward_high as usize],
        contig_offset: offset + forward_low,
    };
    let alignment = match workspace.align(aligned, walk.strand, scores) {
        Ok(alignment) => alignment,
        Err(AlignmentError::NoAlignment) => return Ok(None),
        Err(error) => return Err(error),
    };
    let kernel_score = alignment.score;
    // A completion reaches a rows beyond the alignment, which the window keeps.
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
    if alignment.score >= kernel_score && i64::from(alignment.score) > outside {
        *optimum = Some(DenseOptimum {
            strand: walk.strand,
            query_start,
            contig: (low, high),
            diagonals,
        });
    }
    Ok(Some(alignment))
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
