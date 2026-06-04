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

mod config;
#[cfg(test)]
mod paths;
#[cfg(test)]
mod place;
#[cfg(test)]
mod query;
#[cfg(test)]
mod window;

pub use config::AlignStageConfig;
#[cfg(test)]
use query::PreparedQuery;

use super::Stopped;
use crate::align::AlignmentError;
use crate::decide::DecideError;
use crate::index::IndexError;
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
