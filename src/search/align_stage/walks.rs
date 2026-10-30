//! Aligning the dense walks of one window: copies in dense repeats and low-complexity runs.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::contig::{ContigJob, WindowOutput};
use super::paths::{AlignedPath, PathIndex};
use super::place::query_room;
use super::window::PreparedWindow;
use super::{AlignStageError, WindowWorkspace};
use crate::align::{
    Alignment, DenseResults, DenseWalk, EditOperation, JoinContext, JoinLimits, JoinScratch, Walk,
    align_dense_walk, gap_continuation, holds_walk, long_enough, walk_alignment,
};
use crate::decide::Fragment;
use crate::seeds::base_code;
use std::cmp::Reverse;

/// Walks shorter than this, or than the rule's fragment minimum when it is larger, are tried for
/// joins with their neighbours. A join counts from the rule's fragment minimum.
const SHORT_WALK_LENGTH: u64 = 100;

/// Aligns the walks of `window`: every long walk is an ungapped fragment, and long walks with
/// changes are also aligned with gaps (kept once per copy). Short walks are joined later, by
/// `join_short_walks`.
pub(super) fn align_walks(
    job: &ContigJob<'_>,
    window: &PreparedWindow,
    output: &mut WindowOutput,
    workspace: &mut WindowWorkspace,
) -> Result<(), AlignStageError> {
    if window.walks.is_empty() {
        return Ok(());
    }
    let gapped = walk_fragments(job, window, output, workspace)?;
    keep_gapped(job, window, output, gapped)
}

/// Keeps the ungapped fragment of every long walk that no earlier alignment covers, and returns
/// the gapped alignments of walks that ask for one, with the walk index.
fn walk_fragments(
    job: &ContigJob<'_>,
    window: &PreparedWindow,
    output: &mut WindowOutput,
    workspace: &mut WindowWorkspace,
) -> Result<Vec<(Alignment, usize)>, AlignStageError> {
    let (query, config) = (job.query, job.config);
    let mut results = DenseResults::default();
    let mut gapped = Vec::new();
    for (index, walk) in window.walks.iter().enumerate() {
        job.control.check()?;
        // A walk shorter than a chain is only joined. A gapped alignment around it would fill
        // a box of millions of cells, mostly to find a few dozen matching bases.
        if !long_enough(walk, &query.table) {
            continue;
        }
        if output
            .paths
            .covers_walk(walk, query.length(), query.circular())
        {
            output.work.redundant_chains += 1;
            continue;
        }
        output.work.chains += 1;
        output.work.alignments += 1;
        let oriented = window.oriented(walk);
        let alignment = walk_alignment(walk, oriented, window.start, &query.table, config.band)?;
        let walked = WalkAlignment {
            walk,
            later: &window.walks[index + 1..],
            alignment: &alignment,
        };
        if let Some(joined) =
            gapped_alignment(job, window, output, workspace, walked, &mut results)?
        {
            gapped.push((joined, index));
        }
        let room = walk_room(walk, job);
        if output.keep_or_widen(job, window, alignment, room)? {
            output.paths.add_walk(walk);
        }
    }
    Ok(gapped)
}

/// A walk, the walks of its window after it, and its ungapped alignment.
#[derive(Clone, Copy)]
struct WalkAlignment<'a> {
    walk: &'a Walk,
    later: &'a [Walk],
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
    output: &mut WindowOutput,
    workspace: &mut WindowWorkspace,
    walked: WalkAlignment<'_>,
    results: &mut DenseResults,
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
    let input = DenseWalk {
        walk,
        later: walked.later,
        window: &window.bases,
        offset: window.start,
        table: &query.table,
    };
    let aligned = align_dense_walk(input, config.band, &mut workspace.alignment, results)?;
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
    output: &mut WindowOutput,
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
        )?);
    }
    Ok(())
}

/// Joins walks shorter than a fragment that no counted fragment of the window (in `aligned` or a
/// join made before) holds, longest first, across the divergent pieces between them. Returns the
/// joins that reach the fragment limits of the rule, whatever their score.
pub(super) fn join_short_walks(
    job: &ContigJob<'_>,
    window: &PreparedWindow,
    aligned: &[Fragment],
    scratch: &mut JoinScratch,
) -> Result<Vec<Fragment>, AlignStageError> {
    let (query, rule, walks) = (job.query, job.rule, &window.walks);
    let mut joined: Vec<Fragment> = Vec::new();
    let holds = |fragment: &Fragment, walk: &Walk| {
        rule.counts_fragment(&fragment.alignment)
            && holds_walk(
                &fragment.alignment,
                walk,
                window.bases.len(),
                window.start,
                query.length(),
            )
    };
    let held = |fragments: &[Fragment], walk: &Walk| {
        fragments.iter().any(|fragment| holds(fragment, walk))
    };
    // A window with many fragments looks up the fragments that can hold a walk on the contig.
    #[cfg(test)]
    let off = crate::align::SAME_OFF.load(std::sync::atomic::Ordering::Relaxed);
    #[cfg(not(test))]
    let off = false;
    let index = (held_index_from(aligned.len()) && !off).then(|| HeldIndex::new(aligned));
    let held_aligned = |walk: &Walk| match &index {
        Some(index) => index
            .near(walk, window.bases.len() as u64, window.start)
            .any(|at| holds(&aligned[at], walk)),
        None => held(aligned, walk),
    };
    // The deadline is checked for each walk, as a window with many walks takes long here.
    let mut short: Vec<usize> = Vec::new();
    let shorter = SHORT_WALK_LENGTH.max(rule.min_fragment_columns);
    for (index, walk) in walks.iter().enumerate() {
        job.control.check()?;
        if u64::from(walk.length) < shorter && !held_aligned(walk) {
            short.push(index);
        }
    }
    if short.is_empty() {
        return Ok(joined);
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
    // No fragment of `aligned` holds a short walk, so only the joins made before are asked.
    for index in short {
        job.control.check()?;
        if held(&joined, &walks[index]) {
            continue;
        }
        let oriented = window.oriented(&walks[index]);
        if let Some(join) = scratch.join_through(walks, index, oriented, window.start, &context)? {
            joined.push(Fragment::place(
                job.contig,
                join,
                query.length(),
                query.circular(),
            )?);
        }
    }
    Ok(joined)
}

/// Fragments of a window above this many are looked up through `HeldIndex`.
const HELD_INDEX_FROM: usize = 64;

/// Whether the fragments of a window that hold a walk are looked up through `HeldIndex`.
fn held_index_from(fragments: usize) -> bool {
    fragments > HELD_INDEX_FROM
}

/// Fragments in one group of `HeldIndex`, which keeps the largest contig end of each group.
const HELD_GROUP: usize = 32;

/// The fragments of a window by strand and contig start, with the largest contig end of each
/// group of `HELD_GROUP`, to find those that can hold a walk without asking every fragment.
struct HeldIndex {
    /// Per strand (forward, reverse): (contig start, contig end, fragment) by contig start.
    sorted: [Vec<(u64, u64, usize)>; 2],
    ends: [Vec<u64>; 2],
}

impl HeldIndex {
    fn new(fragments: &[Fragment]) -> Self {
        let mut sorted: [Vec<(u64, u64, usize)>; 2] = [Vec::new(), Vec::new()];
        for (at, fragment) in fragments.iter().enumerate() {
            let contig = fragment.alignment.contig_interval;
            let strand = usize::from(fragment.alignment.strand.is_reverse());
            sorted[strand].push((contig.start, contig.end, at));
        }
        let ends = sorted.each_mut().map(|list| {
            list.sort_unstable();
            list.chunks(HELD_GROUP)
                .map(|group| group.iter().map(|&(_, end, _)| end).max().unwrap_or(0))
                .collect()
        });
        Self { sorted, ends }
    }

    /// The fragments on the strand of `walk` whose contig positions hold the walk's, in no
    /// particular order: all that `holds_walk` can accept, and maybe others.
    fn near<'a>(
        &'a self,
        walk: &Walk,
        length: u64,
        offset: u64,
    ) -> impl Iterator<Item = usize> + 'a {
        // As in `holds_walk`: the walk's contig positions on the forward strand.
        let (low, high) = (u64::from(walk.contig), u64::from(walk.contig + walk.length));
        let (low, high) = if walk.strand.is_reverse() {
            (length - high, length - low)
        } else {
            (low, high)
        };
        let (low, high) = (offset + low, offset + high);
        let strand = usize::from(walk.strand.is_reverse());
        let (sorted, ends) = (&self.sorted[strand], &self.ends[strand]);
        let before = sorted.partition_point(|&(start, _, _)| start <= low);
        (0..before.div_ceil(HELD_GROUP))
            .rev()
            .filter(move |&group| ends[group] >= high)
            .flat_map(move |group| {
                let end = ((group + 1) * HELD_GROUP).min(before);
                sorted[group * HELD_GROUP..end].iter()
            })
            .filter(move |&&(_, end, _)| end >= high)
            .map(|&(_, _, at)| at)
    }
}

#[cfg(test)]
mod held_tests {
    use super::*;
    use crate::align::{EditOperation, EditRun, Interval, Strand};
    use crate::test_support::rng;
    use rand::RngExt;

    #[test]
    fn index_from_65_fragments() {
        assert!(!held_index_from(64));
        assert!(held_index_from(65));
    }

    #[test]
    fn index_finds_every_holding_fragment() {
        // Fragments and walks on both strands of a window of 5,000 bases from contig position
        // 1,000 and a circular query of 800 bases: the index offers every fragment that holds.
        let (length, offset, query_length) = (5_000u64, 1_000u64, 800u64);
        let mut random = rng(70);
        let strand = |reverse: bool| {
            if reverse {
                Strand::Reverse
            } else {
                Strand::Forward
            }
        };
        let fragments: Vec<Fragment> = (0..300)
            .map(|_| {
                let span = random.random_range(20..400u64);
                let contig = offset + random.random_range(0..length - span);
                let query = random.random_range(0..query_length);
                let alignment = Alignment {
                    score: 2 * span as i32,
                    strand: strand(random.random_range(0..2) == 1),
                    query_interval: Interval::new(query, query + span).unwrap(),
                    contig_interval: Interval::new(contig, contig + span).unwrap(),
                    runs: vec![EditRun {
                        operation: EditOperation::Match,
                        length: span as u32,
                    }],
                };
                Fragment::place(0, alignment, query_length, true).unwrap()
            })
            .collect();
        let index = HeldIndex::new(&fragments);
        let mut held = 0;
        for _ in 0..2_000 {
            let span = random.random_range(1..60u32);
            let walk = Walk {
                strand: strand(random.random_range(0..2) == 1),
                contig: random.random_range(0..length as u32 - span),
                query: random.random_range(0..query_length as u32),
                length: span,
            };
            let holds = |fragment: &Fragment| {
                holds_walk(
                    &fragment.alignment,
                    &walk,
                    length as usize,
                    offset,
                    query_length,
                )
            };
            let mut offered: Vec<usize> = index.near(&walk, length, offset).collect();
            offered.sort_unstable();
            let holding: Vec<usize> = (0..fragments.len())
                .filter(|&at| holds(&fragments[at]))
                .collect();
            held += holding.len();
            assert!(holding.iter().all(|at| offered.binary_search(at).is_ok()));
        }
        assert!(held > 100, "{held}");
    }
}
