//! Aligning the dense walks of one window: copies in dense repeats and low-complexity runs.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::contig::{ContigJob, WindowOutput};
use super::paths::{AlignedPath, PathIndex};
use super::place::query_room;
use super::window::PreparedWindow;
use super::{AlignStageError, AlignStageWorkspace};
use crate::align::{
    Alignment, DenseOptimum, EditOperation, JoinContext, JoinLimits, JoinScratch, Walk,
    align_dense_walk, gap_continuation, holds_walk, long_enough, walk_alignment,
};
use crate::decide::Fragment;
use crate::seeds::base_code;
use std::cmp::Reverse;

/// Aligns the walks of `window`: every long walk is an ungapped fragment, walks with changes are
/// also aligned with gaps (kept once per copy), and short walks are joined across their
/// divergent pieces when the join reaches the fragment limits.
pub(super) fn align_walks(
    job: &ContigJob<'_>,
    window: &PreparedWindow,
    output: &mut WindowOutput<'_>,
    workspace: &mut AlignStageWorkspace,
) -> Result<(), AlignStageError> {
    if window.walks.is_empty() {
        return Ok(());
    }
    let gapped = walk_fragments(job, window, output, workspace)?;
    keep_gapped(job, window, output, gapped)?;
    join_short_walks(job, window, output, &mut workspace.joins)
}

/// Keeps the ungapped fragment of every long walk that no earlier alignment covers, and returns
/// the gapped alignments of walks that ask for one, with the walk index.
fn walk_fragments(
    job: &ContigJob<'_>,
    window: &PreparedWindow,
    output: &mut WindowOutput<'_>,
    workspace: &mut AlignStageWorkspace,
) -> Result<Vec<(Alignment, usize)>, AlignStageError> {
    let (query, config) = (job.query, job.config);
    let mut optimum = None;
    let mut gapped = Vec::new();
    for (index, walk) in window.walks.iter().enumerate() {
        job.control.check()?;
        // A walk shorter than a chain only starts a gapped alignment or a join.
        let short = !long_enough(walk, &query.table);
        if output
            .paths
            .covers_walk(walk, query.length(), query.circular())
        {
            output.work.redundant_chains += u64::from(!short);
            continue;
        }
        if !short {
            output.work.chains += 1;
            output.work.alignments += 1;
        }
        let oriented = window.oriented(walk);
        let alignment = walk_alignment(walk, oriented, window.start, &query.table, config.band)?;
        let walked = WalkAlignment {
            walk,
            alignment: &alignment,
        };
        if let Some(joined) =
            gapped_alignment(job, window, output, workspace, walked, &mut optimum)?
        {
            gapped.push((joined, index));
        }
        if short {
            continue;
        }
        let room = walk_room(walk, job);
        if output.keep_or_widen(job, window, alignment, room)? {
            output.paths.add_walk(walk);
        }
    }
    Ok(gapped)
}

/// A walk and its ungapped alignment.
#[derive(Clone, Copy)]
struct WalkAlignment<'a> {
    walk: &'a Walk,
    alignment: &'a Alignment,
}

/// Whether the query has bases before and after a walk.
fn walk_room(walk: &Walk, job: &ContigJob<'_>) -> (bool, bool) {
    let (length, span) = (job.query.length(), u64::from(walk.length));
    if job.query.circular() {
        (span < length, span < length)
    } else {
        let start = u64::from(walk.query);
        (start > 0, start + span < length)
    }
}

/// Whether the ungapped alignment of a walk has a substitution of two real bases, or a run of
/// several substitutions. A single substitution against another byte (such as an IUPAC code)
/// keeps its compact ungapped path.
fn mismatched(walked: WalkAlignment<'_>, oriented: &[u8], query: &[u8]) -> bool {
    let walk = walked.walk;
    let mut at = 0usize;
    walked.alignment.runs.iter().any(|run| {
        let from = at;
        at += run.length as usize;
        run.operation == EditOperation::Substitution
            && (run.length > 1
                || (from..at).any(|at| {
                    base_code(oriented[walk.contig as usize + at]).is_some()
                        && base_code(query[(walk.query as usize + at) % query.len()]).is_some()
                }))
    })
}

/// The gapped alignment around a walk that differs from the query, or whose end can reach a
/// conserved block on another diagonal. `None` when the walk asks for none, kernel alignments
/// already pair it, the alignment has no gap or scores below the walk, or it widens the window.
fn gapped_alignment(
    job: &ContigJob<'_>,
    window: &PreparedWindow,
    output: &mut WindowOutput<'_>,
    workspace: &mut AlignStageWorkspace,
    walked: WalkAlignment<'_>,
    optimum: &mut Option<DenseOptimum>,
) -> Result<Option<Alignment>, AlignStageError> {
    let (query, config) = (job.query, job.config);
    let walk = walked.walk;
    let oriented = window.oriented(walk);
    let asks = mismatched(walked, oriented, query.bases())
        || (u64::from(walk.length) < query.length()
            && gap_continuation(
                walk,
                oriented,
                &query.table,
                config.band,
                &mut workspace.probe,
            )?);
    let band = i64::from(config.band.max_band());
    if !asks
        || output
            .paths
            .aligned_around(walk, band, query.length(), query.circular())
    {
        return Ok(None);
    }
    let aligned = align_dense_walk(
        walk,
        &window.bases,
        window.start,
        &query.table,
        config.band,
        &mut workspace.alignment,
        optimum,
    )?;
    let Some(joined) = aligned else {
        return Ok(None);
    };
    let counts = joined.counts();
    if counts.insertions + counts.deletions == 0 || joined.score < walked.alignment.score {
        return Ok(None);
    }
    let room = query_room(joined.query_interval, query);
    if output.widen_at_edges(job, window, &joined, room) {
        return Ok(None);
    }
    output.paths.add(&AlignedPath::of(
        &joined,
        window.bases.len() as u64,
        window.start,
    ));
    Ok(Some(joined))
}

/// Keeps the gapped walk alignments, best first, except those whose paired bases all lie on a
/// better one. Kernels started from several walks of one copy (as on other rotations of a
/// circular query) find parts of its best alignment: duplicate support, not another copy.
fn keep_gapped(
    job: &ContigJob<'_>,
    window: &PreparedWindow,
    output: &mut WindowOutput<'_>,
    mut gapped: Vec<(Alignment, usize)>,
) -> Result<(), AlignStageError> {
    gapped.sort_unstable_by_key(|(alignment, index)| (Reverse(alignment.score), *index));
    let length = job.query.length();
    let mut kept = PathIndex::default();
    for (alignment, _) in gapped {
        job.control.check()?;
        let path = AlignedPath::of(&alignment, window.bases.len() as u64, window.start);
        let runs: Vec<Walk> = path.walks(length).collect();
        if runs.iter().all(|run| kept.covers_walk(run, length, true)) {
            continue;
        }
        runs.iter().for_each(|run| kept.add_walk(run));
        output.work.alignments += 1;
        output.fragments.push(Fragment::place(
            job.contig,
            alignment,
            job.query.length(),
            job.query.circular(),
        )?)?;
    }
    Ok(())
}

/// Joins walks shorter than a fragment that no counted fragment of the contig holds, longest
/// first, across the divergent pieces between them. A join is kept when it reaches the fragment
/// limits of the rule, whatever its score.
fn join_short_walks(
    job: &ContigJob<'_>,
    window: &PreparedWindow,
    output: &mut WindowOutput<'_>,
    scratch: &mut JoinScratch,
) -> Result<(), AlignStageError> {
    // The fragments of a round that widens are dropped, so no join is needed.
    if *output.grown {
        return Ok(());
    }
    let (query, rule, walks) = (job.query, job.rule, &window.walks);
    let held = |output: &WindowOutput<'_>, walk: &Walk| {
        output.fragments.since(output.first).iter().any(|fragment| {
            rule.counts_fragment(&fragment.alignment)
                && holds_walk(
                    &fragment.alignment,
                    walk,
                    window.bases.len(),
                    window.start,
                    query.length(),
                )
        })
    };
    let mut short: Vec<usize> = (0..walks.len())
        .filter(|&index| {
            u64::from(walks[index].length) < rule.min_fragment_columns
                && !held(output, &walks[index])
        })
        .collect();
    if short.is_empty() {
        return Ok(());
    }
    let context = JoinContext {
        table: &query.table,
        band: job.config.band,
        limits: JoinLimits {
            min_columns: rule.min_fragment_columns,
            min_identity: rule.min_fragment_identity,
        },
    };
    let bytes = JoinScratch::bytes(walks.len(), &context);
    if bytes > job.config.max_prepared_bytes {
        return Err(AlignStageError::TooLarge {
            what: "the joins of one window",
            needed: bytes as u64,
            limit: job.config.max_prepared_bytes as u64,
        });
    }
    scratch.prepare(walks, &context);
    // Longer walks first, so shifted copies inside a kept join are not joined again.
    short.sort_unstable_by_key(|&index| (Reverse(walks[index].length), index));
    for index in short {
        job.control.check()?;
        if held(output, &walks[index]) {
            continue;
        }
        let oriented = window.oriented(&walks[index]);
        let joined = scratch.join_through(walks, index, oriented, window.start, &context)?;
        if let Some(joined) = joined {
            output.work.alignments += 1;
            output.fragments.push(Fragment::place(
                job.contig,
                joined,
                query.length(),
                query.circular(),
            )?)?;
        }
    }
    Ok(())
}
