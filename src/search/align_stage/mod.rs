//! Stage 3: reads the alignment windows of a located entry, aligns the query there and decides
//! the entry.
//!
//! For every contig the windows are read, anchors are chained and the chains aligned (`chains`),
//! and copies in dense repeats are walked and aligned (`walks`). When an alignment or a chain
//! reaches a cut window edge with query bases left beyond it, the window is widened and the
//! contig is aligned again (`contig`), so an alignment is never cut by a window edge. The
//! alignments become fragments, and `DecisionRule` decides the entry from them.
//!
//! Memory of one entry is bounded by the limits of `AlignStageConfig`: the window bases of one
//! contig, the bytes of its prepared windows, and the bytes of the entry's fragments. The kernel
//! keeps at most `max_cells` traceback cells resident per pass.

mod chains;
mod config;
mod contig;
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
use crate::index::{IndexError, PackedSequence, Shard};
use place::Collected;
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

/// Reusable buffers of one worker thread, kept across entries.
#[derive(Default)]
pub struct AlignStageWorkspace {
    alignment: AlignmentWorkspace,
    windows: window::WindowBuffers,
    /// Query k-mers near the end of a dense walk, for `gap_continuation`.
    probe: Vec<(u64, u32)>,
    joins: JoinScratch,
}

/// The work stage 3 did for one entry, for the search summary.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
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
pub fn align_entry(
    shard: &Shard,
    query: &PreparedQuery,
    located: &LocatedEntry,
    config: &AlignStageConfig,
    rule: &DecisionRule,
    workspace: &mut AlignStageWorkspace,
    control: SearchControl<'_>,
) -> Result<AlignedEntry, AlignStageError> {
    control.check()?;
    let packed = shard.packed(located.entry)?;
    align_packed_entry(&packed, query, located, config, rule, workspace, control)
}

/// `align_entry` on an open packed sequence, which may read from any source, such as S3 or a
/// block cache.
pub fn align_packed_entry(
    packed: &PackedSequence,
    query: &PreparedQuery,
    located: &LocatedEntry,
    config: &AlignStageConfig,
    rule: &DecisionRule,
    workspace: &mut AlignStageWorkspace,
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
    workspace: &mut AlignStageWorkspace,
    control: SearchControl<'_>,
) -> Result<(Vec<Fragment>, AlignWork), AlignStageError> {
    config.validate()?;
    rule.validate()?;
    query.check(config)?;
    let mut fragments = Collected::new(config.max_fragment_bytes);
    let mut work = AlignWork::default();
    for windows in &located.windows {
        control.check()?;
        let job = contig::ContigJob::new(packed, windows.contig, query, config, rule, control)?;
        job.align(&windows.ranges, workspace, &mut fragments, &mut work)?;
    }
    Ok((fragments.into_list(), work))
}
