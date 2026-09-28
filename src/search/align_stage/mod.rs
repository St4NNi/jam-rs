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
mod tests;
mod walks;
mod window;

pub use config::AlignStageConfig;
pub use query::PreparedQuery;

use super::{LocatedEntry, SearchControl, Stopped};
use crate::align::{AlignmentError, AlignmentWorkspace, JoinScratch};
use crate::decide::{DecideError, Decision, DecisionRule, Fragment};
use crate::index::{IndexError, PackedSequence};
use place::Collected;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::sync::{Mutex, PoisonError};
use thiserror::Error;

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
    },
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
    let (fragments, work) =
        entry_fragments(packed, query, located, config, rule, workspace, control)?;
    control.check()?;
    let contig_length = |contig: u32| {
        packed
            .contig_range(contig as usize)
            .map_or(0, |range| range.end - range.start)
    };
    let decision = rule.decide(query.length(), fragments, contig_length)?;
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
    config.validate()?;
    rule.validate()?;
    query.check(config)?;
    // Contigs are aligned on all threads, as an entry may hold thousands of contigs with one
    // window each; their fragments are added in contig order, so the limit acts as on one thread.
    let aligned: Vec<Result<(Vec<Fragment>, AlignWork), AlignStageError>> = located
        .alignment_windows
        .par_iter()
        .map(|windows| {
            control.check()?;
            let job = contig::ContigJob::new(packed, windows.contig, query, config, rule, control)?;
            let mut fragments = Collected::new(config.max_fragment_bytes);
            let mut work = AlignWork::default();
            job.align(&windows.ranges, workspace, &mut fragments, &mut work)?;
            Ok((fragments.into_list(), work))
        })
        .collect();
    let mut fragments = Collected::new(config.max_fragment_bytes);
    let mut work = AlignWork::default();
    for result in aligned {
        let (contig_fragments, contig_work) = result?;
        work.add(&contig_work);
        for fragment in contig_fragments {
            fragments.push(fragment)?;
        }
    }
    Ok((fragments.into_list(), work))
}
