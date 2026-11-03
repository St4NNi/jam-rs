//! The windows of one contig, aligned in rounds: when a chain or an alignment reaches a cut window
//! edge with query bases left beyond it, the window grows and the contig is aligned again.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::paths::PathIndex;
use super::place::Collected;
use super::window::{PreparedWindow, merge_spans, prepare_window, widen};
use super::{
    AlignStageConfig, AlignStageError, AlignStageWorkspace, AlignWork, Limit, PreparedQuery,
    SearchControl, WindowWorkspace, chains, walks,
};
use crate::align::{Alignment, Interval};
use crate::decide::DecisionRule;
use crate::decide::Fragment;
use crate::index::PackedSequence;
use rayon::prelude::*;
use std::ops::Range;

/// Bases from a cut window edge within which a chain asks for a wider window: the chaining gap,
/// the extension and one anchor.
fn edge_reach(config: &AlignStageConfig) -> u64 {
    u64::from(config.chains.max_gap) + u64::from(config.extension) + u64::from(config.anchors.a)
}

/// Diagonals beyond a chain's anchors that its alignments can reach: the widest band of a
/// drifting chain reaches 64 past its anchors, a narrow one 128, and endpoint completion adds
/// at most one anchor length.
const CHAIN_SLACK: i64 = 192;

/// Diagonals beyond a walk that its alignments can reach: a gapped box over walks within one
/// band of diagonal (two bands of 1,088), and joins through pieces of at most 1,039 bases.
const WALK_SLACK: i64 = 4_096;

/// A contig's windows after the rounds that grew before aligning, with their chains and walks
/// unless they were let go to save memory.
pub(super) struct ContigRound {
    spans: Vec<Range<u64>>,
    step: u64,
    pub prepared: Option<Vec<PreparedWindow>>,
}

/// One contig of an entry and everything its windows are aligned with.
pub(super) struct ContigJob<'a> {
    pub packed: &'a PackedSequence,
    pub contig: u32,
    pub contig_length: u64,
    pub query: &'a PreparedQuery,
    pub config: &'a AlignStageConfig,
    pub rule: &'a DecisionRule,
    pub control: SearchControl<'a>,
}

/// What aligning one window of a round gives: the window's span (grown when an alignment reaches
/// a cut edge), whether it grew, the paths aligned in the window, its fragments and its counters.
pub(super) struct WindowOutput {
    pub span: Range<u64>,
    pub step: u64,
    pub grown: bool,
    pub paths: PathIndex,
    pub fragments: Vec<Fragment>,
    pub work: AlignWork,
}

impl<'a> ContigJob<'a> {
    pub fn new(
        packed: &'a PackedSequence,
        contig: u32,
        query: &'a PreparedQuery,
        config: &'a AlignStageConfig,
        rule: &'a DecisionRule,
        control: SearchControl<'a>,
    ) -> Result<Self, AlignStageError> {
        let range = packed.contig_range(contig as usize)?;
        Ok(Self {
            packed,
            contig,
            contig_length: range.end - range.start,
            query,
            config,
            rule,
            control,
        })
    }

    /// Reads the windows `ranges` of the contig and finds their chains and walks, widening them
    /// until no chain or walk lies near a cut edge. Nothing is aligned yet.
    pub fn prepare(
        &self,
        ranges: &[Range<u64>],
        workspace: &AlignStageWorkspace,
        work: &mut AlignWork,
    ) -> Result<ContigRound, AlignStageError> {
        let spans = self.checked_spans(ranges)?;
        let step = self.config.widening.max(edge_reach(self.config));
        self.prepare_until_still(spans, step, workspace, work)
    }

    /// Prepare rounds from `spans` until no window grows before aligning.
    fn prepare_until_still(
        &self,
        mut spans: Vec<Range<u64>>,
        mut step: u64,
        workspace: &AlignStageWorkspace,
        work: &mut AlignWork,
    ) -> Result<ContigRound, AlignStageError> {
        loop {
            self.control.check()?;
            let (prepared, grown) = self.prepare_round(&mut spans, step, workspace, work)?;
            if !grown {
                return Ok(ContigRound {
                    spans,
                    step,
                    prepared: Some(prepared),
                });
            }
            work.widenings += 1;
            merge_spans(&mut spans);
            step = step.saturating_mul(2);
        }
    }

    /// Aligns a prepared round, widening and aligning again until no alignment is cut by a
    /// window edge, and adds the fragments. A round without its windows is read again.
    pub fn finish(
        &self,
        round: ContigRound,
        workspace: &AlignStageWorkspace,
        fragments: &mut Collected,
        work: &mut AlignWork,
    ) -> Result<(), AlignStageError> {
        let ContigRound {
            mut spans,
            mut step,
            prepared,
        } = round;
        let mut prepared = match prepared {
            Some(prepared) => prepared,
            None => {
                let round = self.prepare_until_still(spans, step, workspace, work)?;
                (spans, step) = (round.spans, round.step);
                round.prepared.unwrap_or_default()
            }
        };
        let first = fragments.mark();
        loop {
            if !self.align_round(&prepared, &mut spans, step, workspace, fragments, work)? {
                return Ok(());
            }
            fragments.reset(first);
            work.widenings += 1;
            merge_spans(&mut spans);
            step = step.saturating_mul(2);
            let round = self.prepare_until_still(spans, step, workspace, work)?;
            (spans, step) = (round.spans, round.step);
            prepared = round.prepared.unwrap_or_default();
        }
    }

    /// The query bases that the fragments of a prepared round can cover, if the round is
    /// aligned without growing: each chain and walk paired along its band, over the whole
    /// window. `None` when a window with a chain or walk has a cut edge, as it may still grow.
    pub fn reach(&self, prepared: &[PreparedWindow]) -> Option<Vec<Interval>> {
        let mut reach = Vec::new();
        let length = self.query.length() as i64;
        for window in prepared {
            let width = window.bases.len() as i64;
            let chained = window.chains.iter().map(|chain| {
                let diagonals = chain
                    .anchors
                    .iter()
                    .map(|anchor| i64::from(anchor.contig) - i64::from(anchor.query));
                let low = diagonals.clone().min().unwrap_or(0);
                (
                    low - CHAIN_SLACK,
                    diagonals.max().unwrap_or(0) + CHAIN_SLACK,
                )
            });
            let walked = window.walks.iter().map(|walk| {
                let diagonal = i64::from(walk.contig) - i64::from(walk.query);
                (diagonal - WALK_SLACK, diagonal + WALK_SLACK)
            });
            let mut bands = chained.chain(walked).peekable();
            let range = window.range();
            if bands.peek().is_some() && (range.start > 0 || range.end < self.contig_length) {
                return None;
            }
            // Query position = window position - diagonal, for window positions 0..width.
            for (low, high) in bands {
                let start = (-high).clamp(0, length) as u64;
                let end = (width - low).clamp(0, length) as u64;
                if start < end {
                    reach.push(Interval { start, end });
                }
            }
        }
        Some(Interval::union(reach))
    }

    /// The windows, checked against the contig and merged.
    fn checked_spans(&self, ranges: &[Range<u64>]) -> Result<Vec<Range<u64>>, AlignStageError> {
        let mut spans = Vec::with_capacity(ranges.len());
        for range in ranges {
            if range.start >= range.end || range.end > self.contig_length {
                return Err(AlignStageError::WindowOutsideContig {
                    contig: self.contig,
                    start: range.start,
                    end: range.end,
                });
            }
            spans.push(range.clone());
        }
        merge_spans(&mut spans);
        Ok(spans)
    }

    /// Reads every window of a round and finds its chains and walks, on several threads. Windows
    /// with a chain or a walk near a cut edge grow before anything is aligned.
    fn prepare_round(
        &self,
        spans: &mut [Range<u64>],
        step: u64,
        workspace: &AlignStageWorkspace,
        work: &mut AlignWork,
    ) -> Result<(Vec<PreparedWindow>, bool), AlignStageError> {
        let config = self.config;
        let bases: u64 = spans.iter().map(|span| span.end - span.start).sum();
        if bases > config.max_window_bases {
            return Err(AlignStageError::TooLarge {
                what: "the alignment windows of one contig",
                needed: bases,
                limit: config.max_window_bases,
                kind: Limit::WindowBases,
            });
        }
        let read = |span: &Range<u64>| {
            self.control.check()?;
            let mut read_work = AlignWork::default();
            let window = workspace.lend(|buffers| {
                let range = span.clone();
                let (packed, query) = (self.packed, self.query);
                prepare_window(
                    packed,
                    self.contig,
                    range,
                    query,
                    config,
                    &mut buffers.windows,
                    &mut read_work,
                    self.control,
                )
            })?;
            Ok((window, read_work))
        };
        let read: Vec<Result<_, AlignStageError>> = spans.par_iter().map(read).collect();
        let (k, reach) = (u64::from(config.anchors.a), edge_reach(config));
        let mut prepared = Vec::with_capacity(spans.len());
        let (mut bytes, mut grown) = (0usize, false);
        // In window order, so the first failure and the byte limit act as one thread would.
        let mut read = read.into_iter();
        for (span, result) in spans.iter_mut().zip(read.by_ref()) {
            let (window, read_work) = result?;
            work.add(&read_work);
            bytes += window.bytes();
            if bytes > config.max_prepared_bytes {
                // The windows after it count too, so the error names what the round needs.
                let rest: usize = read.flatten().map(|(window, _)| window.bytes()).sum();
                return Err(AlignStageError::TooLarge {
                    what: "the prepared windows of one contig",
                    needed: (bytes + rest) as u64,
                    limit: config.max_prepared_bytes as u64,
                    kind: Limit::PreparedBytes,
                });
            }
            let edges = window.near_edges(k, reach);
            grown |= widen(span, &window.range(), edges, step, self.contig_length);
            prepared.push(window);
        }
        Ok((prepared, grown))
    }

    /// Aligns the chains and dense walks of every window on several threads, then joins the short
    /// walks of the windows before the first one that grows, and adds the fragments in window
    /// order. Returns whether a window grew.
    fn align_round(
        &self,
        prepared: &[PreparedWindow],
        spans: &mut [Range<u64>],
        step: u64,
        workspace: &AlignStageWorkspace,
        fragments: &mut Collected,
        work: &mut AlignWork,
    ) -> Result<bool, AlignStageError> {
        let aligned: Vec<Result<WindowOutput, AlignStageError>> = prepared
            .par_iter()
            .zip(spans.par_iter())
            .map(|(window, span)| {
                let output = WindowOutput {
                    span: span.clone(),
                    step,
                    grown: false,
                    paths: PathIndex::default(),
                    fragments: Vec::new(),
                    work: AlignWork::default(),
                };
                workspace.lend(|buffers| self.align_window(window, output, buffers))
            })
            .collect();
        // Joins are made only while no window has grown, as their fragments would be dropped.
        let joined = aligned
            .iter()
            .position(|output| output.as_ref().map_or(true, |output| output.grown))
            .unwrap_or(aligned.len());
        let joins: Vec<Result<Vec<Fragment>, AlignStageError>> = prepared[..joined]
            .par_iter()
            .zip(aligned[..joined].par_iter())
            .map(|(window, output)| {
                let held = output.as_ref().map_or(&[][..], |output| &output.fragments);
                workspace
                    .lend(|buffers| walks::join_short_walks(self, window, held, &mut buffers.joins))
            })
            .collect();
        let mut joins = joins.into_iter();
        let mut grown = false;
        for (span, output) in spans.iter_mut().zip(aligned) {
            let output = output?;
            *span = output.span;
            grown |= output.grown;
            work.add(&output.work);
            for fragment in output.fragments {
                fragments.push(fragment);
            }
            if let Some(join) = joins.next() {
                let joined = join?;
                work.alignments += joined.len() as u64;
                for fragment in joined {
                    fragments.push(fragment);
                }
            }
        }
        Ok(grown)
    }

    /// Aligns the chains and then the dense walks of one window.
    fn align_window(
        &self,
        window: &PreparedWindow,
        mut output: WindowOutput,
        workspace: &mut WindowWorkspace,
    ) -> Result<WindowOutput, AlignStageError> {
        self.control.check()?;
        chains::align_chains(self, window, &mut output, &mut workspace.alignment)?;
        walks::align_walks(self, window, &mut output, workspace)?;
        Ok(output)
    }
}

impl WindowOutput {
    /// Grows the window when `alignment` reaches a cut window edge with query bases left beyond
    /// it (`room`: before and after it on the query). Returns whether it did.
    pub fn widen_at_edges(
        &mut self,
        job: &ContigJob<'_>,
        window: &PreparedWindow,
        alignment: &Alignment,
        room: (bool, bool),
    ) -> bool {
        let read = window.range();
        let k = u64::from(job.config.anchors.a);
        let contig = alignment.contig_interval;
        // On the reverse strand the query runs from the window end to its start.
        let (before, after) = if alignment.strand == crate::align::Strand::Reverse {
            (room.1, room.0)
        } else {
            room
        };
        let left = read.start > 0 && before && contig.start <= read.start + k;
        let right = read.end < job.contig_length && after && contig.end + k >= read.end;
        if left || right {
            widen(
                &mut self.span,
                &read,
                (left, right),
                self.step,
                job.contig_length,
            );
            self.grown = true;
        }
        left || right
    }

    /// Keeps the fragment of `alignment` unless it widens the window. Returns whether it was
    /// kept.
    pub fn keep_or_widen(
        &mut self,
        job: &ContigJob<'_>,
        window: &PreparedWindow,
        alignment: Alignment,
        room: (bool, bool),
    ) -> Result<bool, AlignStageError> {
        if self.widen_at_edges(job, window, &alignment, room) {
            return Ok(false);
        }
        self.fragments.push(Fragment::place(
            job.contig,
            alignment,
            job.query.length(),
            job.query.circular(),
        )?);
        Ok(true)
    }
}
