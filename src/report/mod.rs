//! The output records of a search: JSON lines for every entry event and the summary, and PAF
//! lines with the CIGAR in the `cg:Z` tag. The command line and the network use the same records.
//!
//! A record holds names and numbers only, never a reference into a shard, so it can be sent to
//! another process and read back with serde.

mod decision;
mod event;
mod paf;

pub use decision::{ContigInfo, DecisionRecord, FragmentRecord, Span, Strand};
pub use event::{
    AlignWorkRecord, Ending, LocateWorkRecord, Reason, Record, SketchWorkRecord, Stage,
};
pub use paf::paf_lines;

use thiserror::Error;

/// What can go wrong while building a record.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum ReportError {
    /// A fragment names a contig that the contig table does not hold.
    #[error("contig {0} is not in the contig table of the entry")]
    UnknownContig(u32),
}
