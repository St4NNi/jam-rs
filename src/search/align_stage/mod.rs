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
