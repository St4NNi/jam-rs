//! The windows of one contig, aligned in rounds: when a chain or an alignment reaches a cut window
//! edge with query bases left beyond it, the window grows and the contig is aligned again.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::paths::PathIndex;
use super::place::Collected;
use super::window::{PreparedWindow, merge_spans, prepare_window, widen};
use super::{
    AlignStageConfig, AlignStageError, AlignStageWorkspace, AlignWork, PreparedQuery,
    SearchControl, chains, walks,
};
use crate::align::Alignment;
use crate::decide::DecisionRule;
use crate::decide::Fragment;
use crate::index::PackedSequence;
use std::ops::Range;

/// Bases from a cut window edge within which a chain asks for a wider window: the chaining gap,
/// the extension and one anchor.
fn edge_reach(config: &AlignStageConfig) -> u64 {
    u64::from(config.chains.max_gap) + u64::from(config.extension) + u64::from(config.anchors.a)
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

/// What aligning one window of a round changes: the window's span, whether the round grew a
/// window, the paths aligned in the window, the entry's fragments and the counters.
pub(super) struct WindowOutput<'a> {
    pub span: &'a mut Range<u64>,
    pub step: u64,
    pub grown: &'a mut bool,
    pub paths: PathIndex,
    pub fragments: &'a mut Collected,
    /// Where the fragments of this contig start.
    pub first: super::place::Mark,
    pub work: &'a mut AlignWork,
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

    /// Aligns the query in the windows `ranges` of the contig, widening them until no alignment
    /// is cut by a window edge, and adds the fragments.
    pub fn align(
        &self,
        ranges: &[Range<u64>],
        workspace: &mut AlignStageWorkspace,
        fragments: &mut Collected,
        work: &mut AlignWork,
    ) -> Result<(), AlignStageError> {
        let mut spans = self.checked_spans(ranges)?;
        let mut step = self.config.widening.max(edge_reach(self.config));
        let first = fragments.mark();
        loop {
            self.control.check()?;
            let (prepared, mut grown) = self.prepare_round(&mut spans, step, workspace, work)?;
            if !grown {
                for (window, span) in prepared.iter().zip(spans.iter_mut()) {
                    self.control.check()?;
                    let output = WindowOutput {
                        span,
                        step,
                        grown: &mut grown,
                        paths: PathIndex::default(),
                        fragments,
                        first,
                        work,
                    };
                    self.align_window(window, output, workspace)?;
                }
            }
            for window in prepared {
                workspace.windows.recycle(window);
            }
            if !grown {
                return Ok(());
            }
            fragments.reset(first);
            work.widenings += 1;
            merge_spans(&mut spans);
            step = step.saturating_mul(2);
        }
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

    /// Reads every window of a round and finds its chains and walks. Windows with a chain or a
    /// walk near a cut edge grow before anything is aligned.
    fn prepare_round(
        &self,
        spans: &mut [Range<u64>],
        step: u64,
        workspace: &mut AlignStageWorkspace,
        work: &mut AlignWork,
    ) -> Result<(Vec<PreparedWindow>, bool), AlignStageError> {
        let config = self.config;
        let bases: u64 = spans.iter().map(|span| span.end - span.start).sum();
        if bases > config.max_window_bases {
            return Err(AlignStageError::TooLarge {
                what: "the alignment windows of one contig",
                needed: bases,
                limit: config.max_window_bases,
            });
        }
        let (k, reach) = (u64::from(config.anchors.a), edge_reach(config));
        let mut prepared = Vec::with_capacity(spans.len());
        let (mut bytes, mut grown) = (0usize, false);
        for span in spans.iter_mut() {
            self.control.check()?;
            let window = prepare_window(
                self.packed,
                self.contig,
                span.clone(),
                self.query,
                config,
                &mut workspace.windows,
                work,
            )?;
            bytes += window.bytes();
            if bytes > config.max_prepared_bytes {
                return Err(AlignStageError::TooLarge {
                    what: "the prepared windows of one contig",
                    needed: bytes as u64,
                    limit: config.max_prepared_bytes as u64,
                });
            }
            let edges = window.near_edges(k, reach);
            grown |= widen(span, &window.range(), edges, step, self.contig_length);
            prepared.push(window);
        }
        Ok((prepared, grown))
    }

    /// Aligns the chains and then the dense walks of one window.
    fn align_window(
        &self,
        window: &PreparedWindow,
        mut output: WindowOutput<'_>,
        workspace: &mut AlignStageWorkspace,
    ) -> Result<(), AlignStageError> {
        chains::align_chains(self, window, &mut output, &mut workspace.alignment)?;
        walks::align_walks(self, window, &mut output, workspace)
    }
}

impl WindowOutput<'_> {
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
                self.span,
                &read,
                (left, right),
                self.step,
                job.contig_length,
            );
            *self.grown = true;
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
        )?)?;
        Ok(true)
    }
}
