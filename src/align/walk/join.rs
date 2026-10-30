//! Joins short walks across the divergent pieces between them.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

// A walk is linked to its nearest walk on each side; only the pieces between are aligned, and
// the best run of links that reaches the join limits is kept.

use super::dense::Walk;
use super::ungapped::{oriented_alignment, walk_alignment};
use crate::align::BandConfig;
use crate::align::anchors::QueryTable;
use crate::align::cigar::{append_runs, counts};
use crate::align::kernel::{EndpointWorkspace, PathEnd, affine_path};
use crate::align::{Alignment, AlignmentError, EditCounts, EditOperation, EditRun, bases_match};
use std::mem::size_of;

/// The least columns and identity a join must reach to be kept.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct JoinLimits {
    /// The fewest alignment columns of a join.
    pub min_columns: u64,
    /// The lowest identity of a join.
    pub min_identity: f64,
}

/// What the joins of one alignment window read: the query table, the band and the limits.
#[derive(Clone, Copy)]
pub struct JoinContext<'a> {
    /// The query and its anchor settings.
    pub table: &'a QueryTable,
    /// The scores and the drift a link may cross.
    pub band: BandConfig,
    /// The least columns and identity of a kept join.
    pub limits: JoinLimits,
}

/// Reusable buffers of the joins of one window's walks.
#[derive(Default)]
pub struct JoinScratch {
    starts: Vec<u32>,
    ends: Vec<u32>,
    /// The walks linked through one walk, in contig order.
    linked: Vec<u32>,
    parts: Vec<[Part; 2]>,
    matrix: EndpointWorkspace,
    piece: Vec<u8>,
}

impl JoinScratch {
    /// The bytes `prepare` holds for `walks` walks: walk orders, one set of linked walks, and one piece
    /// alignment of at most the reach on each side with its operations.
    pub fn bytes(walks: usize, context: &JoinContext<'_>) -> usize {
        let reach = reach(context) as usize;
        walks * 2 * size_of::<u32>()
            + walks.min(most_linked(context)) * (size_of::<u32>() + size_of::<[Part; 2]>())
            + (reach + 1) * (reach + 1)
            + 2 * (reach + 1) * size_of::<[i32; 3]>()
            + reach
            + 2 * reach * (size_of::<EditOperation>() + size_of::<EditRun>())
    }

    /// Orders `walks` by strand and contig start and end. The piece matrix grows to at most
    /// what `bytes` priced.
    pub fn prepare(&mut self, walks: &[Walk], context: &JoinContext<'_>) {
        let key = |index: u32, end: bool| {
            let walk = walks[index as usize];
            let contig = walk.contig + if end { walk.length } else { 0 };
            (walk.strand.is_reverse(), contig, index)
        };
        self.starts = (0..walks.len() as u32).collect();
        self.starts.sort_unstable_by_key(|&index| key(index, false));
        self.ends = (0..walks.len() as u32).collect();
        self.ends.sort_unstable_by_key(|&index| key(index, true));
        // Exactly the linked walks priced, so neither grows past them nor keeps more from an
        // earlier window.
        let linked = walks.len().min(most_linked(context));
        self.linked = Vec::with_capacity(linked);
        self.parts = Vec::with_capacity(linked);
        self.piece = Vec::with_capacity(reach(context) as usize);
    }

    /// The bytes the walk orders, linked walks and parts hold now.
    #[cfg(test)]
    pub(crate) fn walk_bytes(&self) -> usize {
        (self.starts.capacity() + self.ends.capacity() + self.linked.capacity()) * size_of::<u32>()
            + self.parts.capacity() * size_of::<[Part; 2]>()
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
        context: &JoinContext<'_>,
    ) -> Result<Option<Alignment>, AlignmentError> {
        self.link_walks(walks, at, context);
        if self.linked.len() < 2 {
            return Ok(None);
        }
        let middle = self
            .linked
            .iter()
            .position(|&index| index as usize == at)
            .unwrap_or(0);
        self.score_parts(walks, oriented, context)?;
        let Some(Selection { score, first, last }) = self.best_run(middle, context.limits) else {
            return Ok(None);
        };
        let mut runs: Vec<EditRun> = Vec::new();
        for position in first..=last {
            let walk = walks[self.linked[position] as usize];
            let alignment = walk_alignment(&walk, oriented, offset, context.table, context.band)?;
            append_runs(&mut runs, alignment.runs)?;
            if position < last {
                let next = walks[self.linked[position + 1] as usize];
                self.piece(walk, next, oriented, context, Some(&mut runs))?;
            }
        }
        let first = walks[self.linked[first] as usize];
        oriented_alignment(runs, score, &first, oriented.len(), offset).map(Some)
    }

    /// Links walks to walk `at` on both sides; a circular join covers the query at most once.
    fn link_walks(&mut self, walks: &[Walk], at: usize, context: &JoinContext<'_>) {
        let length = u64::from(context.table.query_length());
        self.linked.clear();
        self.linked.push(at as u32);
        let mut span = u64::from(walks[at].length);
        for after in [false, true] {
            let mut current = at;
            while let Some(Link {
                walk: next,
                query_gap,
            }) = self.link(walks, current, after, context)
            {
                span += query_gap + u64::from(walks[next].length);
                if span > length {
                    break;
                }
                self.linked.push(next as u32);
                current = next;
            }
            if !after {
                self.linked.reverse();
            }
        }
    }

    /// The nearest walk after (or before) walk `at` on its strand, ordered on query and contig,
    /// with both pieces between them within the reach and drift, and the query piece length.
    fn link(
        &self,
        walks: &[Walk],
        at: usize,
        after: bool,
        context: &JoinContext<'_>,
    ) -> Option<Link> {
        let reach = reach(context);
        let length = u64::from(context.table.query_length());
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
            (other.strand.is_reverse(), contig)
        };
        let first = order.partition_point(|&index| {
            key(&walks[index as usize]) < (walk.strand.is_reverse(), low)
        });
        let mut best: Option<((u64, u64, u32), u64)> = None;
        for &index in &order[first..] {
            let other = walks[index as usize];
            if key(&other) > (walk.strand.is_reverse(), high) {
                break;
            }
            let (before, next) = if after { (walk, other) } else { (other, walk) };
            let contig_gap = u64::from(next.contig) - u64::from(before.contig + before.length);
            let end = i64::from(before.query) + i64::from(before.length);
            let query_gap = i64::from(next.query) - end;
            let query_gap = if context.table.circular() {
                query_gap.rem_euclid(length as i64)
            } else {
                query_gap
            };
            let Ok(query_gap) = u64::try_from(query_gap) else {
                continue;
            };
            if query_gap > reach
                || query_gap.abs_diff(contig_gap) > u64::from(context.band.max_drift)
            {
                continue;
            }
            let rank = (query_gap + contig_gap, contig_gap, index);
            if best.is_none_or(|(saved, _)| rank < saved) {
                best = Some((rank, query_gap));
            }
        }
        best.map(|((_, _, index), query_gap)| Link {
            walk: index as usize,
            query_gap,
        })
    }

    /// The parts of every linked walk and of the piece after it.
    fn score_parts(
        &mut self,
        walks: &[Walk],
        oriented: &[u8],
        context: &JoinContext<'_>,
    ) -> Result<(), AlignmentError> {
        let query = context.table.query();
        self.parts.clear();
        for position in 0..self.linked.len() {
            let walk = walks[self.linked[position] as usize];
            let mut part = Part::default();
            for step in 0..walk.length as usize {
                let same = bases_match(
                    query[(walk.query as usize + step) % query.len()],
                    oriented[walk.contig as usize + step],
                );
                part.score += i64::from(if same {
                    context.band.alignment.match_score
                } else {
                    context.band.alignment.mismatch_score
                });
                part.counts.matches += u64::from(same);
                part.counts.substitutions += u64::from(!same);
            }
            let piece = match self.linked.get(position + 1) {
                Some(&next) => self.piece(walk, walks[next as usize], oriented, context, None)?,
                None => Part::default(),
            };
            self.parts.push([part, piece]);
        }
        Ok(())
    }

    /// The best run of walks `first..=last` around position `middle` that meets the limits, by
    /// score, then the first found.
    fn best_run(&self, middle: usize, limits: JoinLimits) -> Option<Selection> {
        let mut best: Option<Selection> = None;
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
            for last in middle..self.linked.len() {
                if last > middle {
                    whole = whole + self.parts[last - 1][1] + self.parts[last][0];
                }
                let qualifies = whole.counts.columns() >= limits.min_columns
                    && whole.counts.reaches_identity(limits.min_identity);
                if first < last
                    && qualifies
                    && best.as_ref().is_none_or(|saved| whole.score > saved.score)
                {
                    best = Some(Selection {
                        score: whole.score,
                        first,
                        last,
                    });
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
        context: &JoinContext<'_>,
        runs: Option<&mut Vec<EditRun>>,
    ) -> Result<Part, AlignmentError> {
        let query = context.table.query();
        let end = before.query as usize + before.length as usize;
        let query_gap = (next.query as usize + query.len() - end % query.len()) % query.len();
        self.piece.clear();
        context
            .table
            .copy_bases(end as u64, query_gap as u64, &mut self.piece)?;
        let contig = &oriented[(before.contig + before.length) as usize..next.contig as usize];
        let alignment = context.band.alignment;
        let path = affine_path(
            &mut self.matrix,
            &self.piece,
            contig,
            alignment,
            PathEnd::Corner,
        )?;
        let part = Part {
            score: i64::from(path.score),
            counts: counts(&path.runs),
        };
        if let Some(runs) = runs {
            append_runs(runs, path.runs)?;
        }
        Ok(part)
    }
}

/// Score and base counts of a walk or of the piece between two linked walks.
#[derive(Clone, Copy, Default)]
struct Part {
    score: i64,
    counts: EditCounts,
}

impl std::ops::Add for Part {
    type Output = Self;
    fn add(mut self, other: Self) -> Self {
        self.score += other.score;
        self.counts += other.counts;
        self
    }
}

/// The walk that a walk links to, and the query bases between them.
struct Link {
    walk: usize,
    query_gap: u64,
}

/// The linked walks `first..=last` of the best join, and its score.
struct Selection {
    score: i64,
    first: usize,
    last: usize,
}

/// The most walks one join links: each walk starts with an exact anchor of a bases, and the
/// linked walks cover the query at most once.
fn most_linked(context: &JoinContext<'_>) -> usize {
    (context.table.query_length() / u32::from(context.table.anchors().a)) as usize + 1
}

/// The query and contig bases a link may cross, as for `gap_continuation`.
fn reach(context: &JoinContext<'_>) -> u64 {
    u64::from(context.band.max_drift) + u64::from(context.table.anchors().a)
}
