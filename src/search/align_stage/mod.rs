//! Stage 3: reads the alignment windows of a located entry, aligns the query there and decides
//! the entry.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

mod chains;
mod config;
mod contig;
#[cfg(test)]
mod copy_tests;
#[cfg(test)]
mod dense_tests;
mod paths;
mod place;
mod query;
#[cfg(test)]
mod reach_tests;
#[cfg(test)]
mod tests;
mod walks;
mod window;

pub use config::{AlignStageConfig, Limit};
pub use query::PreparedQuery;

use super::{LocatedEntry, SearchControl, Stopped};
use crate::align::{AlignmentError, AlignmentWorkspace, Budgeted, Interval, JoinScratch};
use crate::decide::{DecideError, Decision, DecisionRule, Fragment};
use crate::index::{IndexError, PackedSequence};
use place::Collected;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::sync::{Mutex, PoisonError};
use thiserror::Error;
use window::PreparedWindow;

/// Why stage 3 could not decide an entry.
#[derive(Debug, Error)]
pub enum AlignStageError {
    #[error(transparent)]
    Stopped(#[from] Stopped),
    #[error(transparent)]
    Index(#[from] IndexError),
    #[error(transparent)]
    Alignment(#[from] AlignmentError),
    #[error(transparent)]
    Decide(#[from] DecideError),
    #[error("the query is empty")]
    EmptyQuery,
    #[error("the stage 3 settings are invalid: {0}")]
    InvalidConfig(&'static str),
    #[error("alignment window {start}..{end} lies outside contig {contig}")]
    WindowOutsideContig { contig: u32, start: u64, end: u64 },
    #[error("{what} need {needed} bytes or bases, more than the limit of {limit}")]
    TooLarge {
        what: &'static str,
        needed: u64,
        limit: u64,
        kind: Limit,
    },
}

impl AlignStageError {
    /// The limit that an entry passed, with the bytes, bases or cells it needed when known;
    /// `None` for an error that is not a limit.
    pub fn passed_limit(&self) -> Option<(Limit, Option<u64>)> {
        match self {
            Self::TooLarge { needed, kind, .. } => Some((*kind, Some(*needed))),
            Self::Alignment(AlignmentError::MatrixTooLarge { cells, .. }) => {
                Some((Limit::Cells, Some(*cells as u64)))
            }
            Self::Alignment(AlignmentError::ChainingBudget(Budgeted::ChainingWork)) => {
                Some((Limit::Examined, None))
            }
            Self::Alignment(AlignmentError::ChainingBudget(_)) => Some((Limit::DenseBytes, None)),
            _ => None,
        }
    }
}

/// Reusable buffers of stage 3, kept across entries. The windows of a contig are aligned on
/// several threads, so the buffers are lent to one thread per window; all workers of a search
/// can share one workspace.
#[derive(Default)]
pub struct AlignStageWorkspace {
    idle: Mutex<Vec<WindowWorkspace>>,
}

impl AlignStageWorkspace {
    /// Runs `work` with buffers that no other thread uses meanwhile.
    fn lend<T>(&self, work: impl FnOnce(&mut WindowWorkspace) -> T) -> T {
        let idle = || self.idle.lock().unwrap_or_else(PoisonError::into_inner);
        let mut buffers = idle().pop().unwrap_or_default();
        let result = work(&mut buffers);
        idle().push(buffers);
        result
    }
}

/// The buffers that aligning one window uses.
#[derive(Default)]
struct WindowWorkspace {
    alignment: AlignmentWorkspace,
    windows: window::WindowBuffers,
    /// Query k-mers near the end of a dense walk, for `gap_continuation`.
    probe: rustc_hash::FxHashMap<u64, (u32, u32)>,
    joins: JoinScratch,
}

/// The work stage 3 did for one entry, for the search summary.
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
pub struct AlignWork {
    /// Alignment windows read, counting every widening round.
    pub windows: u64,
    /// Bases read from the packed sequence (about four per byte read).
    pub bases_read: u64,
    /// Anchors found in the windows.
    pub anchors: u64,
    /// Chains and long dense walks aligned.
    pub chains: u64,
    /// Chains and walks not aligned because an earlier alignment already covers them.
    pub redundant_chains: u64,
    /// Alignments made, including those of rounds that widened a window.
    pub alignments: u64,
    /// Widening rounds: a contig aligned again because a window grew.
    pub widenings: u64,
}

impl AlignWork {
    /// Adds the work of another entry.
    pub fn add(&mut self, other: &Self) {
        self.windows += other.windows;
        self.bases_read += other.bases_read;
        self.anchors += other.anchors;
        self.chains += other.chains;
        self.redundant_chains += other.redundant_chains;
        self.alignments += other.alignments;
        self.widenings += other.widenings;
    }
}

/// A located entry after stage 3: its decision, with the mosaic of its fragments.
#[derive(Clone, Debug, PartialEq)]
pub struct AlignedEntry {
    /// The entry's index in the shard.
    pub entry: u32,
    /// Whether the entry is reported, with its mosaic, counts and identity.
    pub decision: Decision,
    /// The work stage 3 did for the entry.
    pub work: AlignWork,
}

/// Aligns the query in the windows of a located entry of `shard` and decides the entry.
#[cfg(test)]
pub fn align_entry(
    shard: &crate::index::Shard,
    query: &PreparedQuery,
    located: &LocatedEntry,
    config: &AlignStageConfig,
    rule: &DecisionRule,
    workspace: &AlignStageWorkspace,
    control: SearchControl<'_>,
) -> Result<AlignedEntry, AlignStageError> {
    control.check()?;
    let packed = shard.packed(located.entry)?;
    align_packed_entry(&packed, query, located, config, rule, workspace, control)
}

/// Aligns the query in the windows of a located entry, read from an open packed sequence (any
/// source, such as S3 or a block cache), and decides the entry from the fragments. Chains and
/// dense walks are aligned per window; a window grows while an alignment reaches its edge.
pub fn align_packed_entry(
    packed: &PackedSequence,
    query: &PreparedQuery,
    located: &LocatedEntry,
    config: &AlignStageConfig,
    rule: &DecisionRule,
    workspace: &AlignStageWorkspace,
    control: SearchControl<'_>,
) -> Result<AlignedEntry, AlignStageError> {
    let request = EntryRequest {
        packed,
        query,
        located,
        config,
        rule,
        workspace,
        control,
    };
    let required = rule.min_query_coverage.required_bases(query.length());
    // An entry whose fragments cannot reach the coverage is decided without them.
    let (fragments, work) = request.fragments(Some(required))?;
    let fragments = fragments.unwrap_or_default();
    control.check()?;
    let contig_length = |contig: u32| {
        packed
            .contig_range(contig as usize)
            .map_or(0, |range| range.end - range.start)
    };
    let stopped = || control.check().is_err();
    let decided = rule.decide(
        query.length(),
        fragments,
        &query.low_complexity,
        stopped,
        contig_length,
    );
    // A stop inside the decision shows here as `Stopped`.
    control.check()?;
    let decision = decided?;
    Ok(AlignedEntry {
        entry: located.entry,
        decision,
        work,
    })
}

/// Aligns the query in the windows of a located entry and returns the fragments without a
/// decision, for callers that decide over other units, such as the screen. The fragment limits
/// of `rule` still apply to dense walks.
pub fn entry_fragments(
    packed: &PackedSequence,
    query: &PreparedQuery,
    located: &LocatedEntry,
    config: &AlignStageConfig,
    rule: &DecisionRule,
    workspace: &AlignStageWorkspace,
    control: SearchControl<'_>,
) -> Result<(Vec<Fragment>, AlignWork), AlignStageError> {
    let request = EntryRequest {
        packed,
        query,
        located,
        config,
        rule,
        workspace,
        control,
    };
    let (fragments, work) = request.fragments(None)?;
    Ok((fragments.unwrap_or_default(), work))
}

/// Prepared windows an entry keeps at once while its reach is computed; past them, a contig's
/// windows are read again when it is aligned.
const KEPT_PREPARED_BYTES: usize = 256 << 20;

/// Batches in which the contigs of an entry are aligned between checks of its reach.
const REACH_BATCHES: usize = 64;

/// The fewest contigs of one batch, so that the threads stay busy.
const LEAST_BATCH: usize = 16;

/// One located entry with everything its contigs are aligned with.
#[derive(Clone, Copy)]
struct EntryRequest<'a> {
    packed: &'a PackedSequence,
    query: &'a PreparedQuery,
    located: &'a LocatedEntry,
    config: &'a AlignStageConfig,
    rule: &'a DecisionRule,
    workspace: &'a AlignStageWorkspace,
    control: SearchControl<'a>,
}

/// A contig of an entry while its entry is aligned: its prepared round, or its result.
enum ContigState {
    Prepared(contig::ContigRound, Option<Vec<Interval>>),
    Done(Result<Vec<Fragment>, AlignStageError>),
}

impl EntryRequest<'_> {
    fn job(&self, contig: u32) -> Result<contig::ContigJob<'_>, AlignStageError> {
        let Self {
            packed,
            query,
            config,
            rule,
            control,
            ..
        } = *self;
        contig::ContigJob::new(packed, contig, query, config, rule, control)
    }

    /// The fragments of the entry in contig order. With `required`, `None` when the fragments
    /// cannot cover that many query bases: the contigs not aligned yet are counted by their
    /// reach, and the entry stops once aligned and unaligned contigs together fall short.
    fn fragments(
        &self,
        required: Option<u64>,
    ) -> Result<(Option<Vec<Fragment>>, AlignWork), AlignStageError> {
        let (config, workspace, control) = (self.config, self.workspace, self.control);
        config.validate()?;
        self.rule.validate()?;
        self.query.check(config)?;
        let required = required.filter(|_| !self.query.circular());
        let windows = &self.located.alignment_windows;
        let kept = std::sync::atomic::AtomicUsize::new(0);
        // Contigs are prepared and aligned on all threads, as an entry may hold thousands of
        // contigs with one window each.
        let prepared: Vec<(ContigState, AlignWork)> = windows
            .par_iter()
            .map(|windows| {
                let mut work = AlignWork::default();
                let round = control
                    .check()
                    .map_err(AlignStageError::from)
                    .and_then(|()| {
                        let job = self.job(windows.contig)?;
                        let mut round = job.prepare(&windows.ranges, workspace, &mut work)?;
                        let reach =
                            required.and(round.prepared.as_deref().and_then(|p| job.reach(p)));
                        let bytes = round
                            .prepared
                            .iter()
                            .flatten()
                            .map(PreparedWindow::bytes)
                            .sum();
                        let order = std::sync::atomic::Ordering::Relaxed;
                        if kept.fetch_add(bytes, order) + bytes > KEPT_PREPARED_BYTES {
                            round.prepared = None;
                            kept.fetch_sub(bytes, order);
                        }
                        Ok((round, reach))
                    });
                let state = match round {
                    Ok((round, reach)) => ContigState::Prepared(round, reach),
                    Err(error) => ContigState::Done(Err(error)),
                };
                (state, work)
            })
            .collect();
        let (mut states, works): (Vec<ContigState>, Vec<AlignWork>) = prepared.into_iter().unzip();
        let mut work = AlignWork::default();
        works.iter().for_each(|contig_work| work.add(contig_work));
        // Contigs that may still grow first, then those that reach the most query bases.
        let mut order: Vec<usize> = (0..states.len()).collect();
        order.sort_by_key(|&at| match &states[at] {
            ContigState::Prepared(_, reach) => (
                reach.is_some(),
                std::cmp::Reverse(reach.as_deref().map_or(0, Interval::total_length)),
            ),
            ContigState::Done(_) => (false, std::cmp::Reverse(0)),
        });
        let batch = order.len().div_ceil(REACH_BATCHES).max(LEAST_BATCH);
        for part in order.chunks(batch) {
            if required.is_some_and(|required| self.short_of(&states, required)) {
                return Ok((None, work));
            }
            let taken: Vec<(usize, contig::ContigRound)> = part
                .iter()
                .filter_map(
                    |&at| match std::mem::replace(&mut states[at], placeholder()) {
                        ContigState::Prepared(round, _) => Some((at, round)),
                        done => {
                            states[at] = done;
                            None
                        }
                    },
                )
                .collect();
            let aligned: Vec<(usize, Result<Vec<Fragment>, AlignStageError>, AlignWork)> = taken
                .into_par_iter()
                .map(|(at, contig_round)| {
                    let mut contig_work = AlignWork::default();
                    let mut fragments = Collected::new(config.max_fragment_bytes);
                    let result = control
                        .check()
                        .map_err(AlignStageError::from)
                        .and_then(|()| self.job(windows[at].contig))
                        .and_then(|job| {
                            job.finish(contig_round, workspace, &mut fragments, &mut contig_work)
                        });
                    (at, result.map(|()| fragments.into_list()), contig_work)
                })
                .collect();
            for (at, result, contig_work) in aligned {
                work.add(&contig_work);
                states[at] = ContigState::Done(result);
            }
        }
        // Their fragments are added in contig order, so the limit acts as on one thread.
        let mut fragments = Collected::new(config.max_fragment_bytes);
        for state in states {
            let ContigState::Done(result) = state else {
                unreachable!("every contig is aligned");
            };
            for fragment in result? {
                fragments.push(fragment)?;
            }
        }
        Ok((Some(fragments.into_list()), work))
    }

    /// Whether the aligned fragments and the reach of the contigs not aligned yet cover fewer
    /// than `required` query bases. Never after a contig failed, so failures stay as they are.
    fn short_of(&self, states: &[ContigState], required: u64) -> bool {
        let mut covered = Vec::new();
        for state in states {
            match state {
                ContigState::Prepared(_, Some(reach)) => covered.extend_from_slice(reach),
                ContigState::Prepared(_, None) | ContigState::Done(Err(_)) => return false,
                ContigState::Done(Ok(fragments)) => covered.extend(
                    fragments
                        .iter()
                        .flat_map(|fragment| fragment.query_segments.iter().copied()),
                ),
            }
        }
        Interval::total_length(&Interval::union(covered)) < required
    }
}

/// A state to hold a contig's place while its round is taken out.
fn placeholder() -> ContigState {
    ContigState::Done(Ok(Vec::new()))
}
