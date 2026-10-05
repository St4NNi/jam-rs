//! The output records of a search: JSON lines for every entry event and the summary, and PAF
//! lines with the CIGAR in the `cg:Z` tag. The command line and the network use the same records.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

mod compact;
mod decision;
mod event;
mod paf;

pub use compact::{QueryTotals, result_line};
pub use decision::{DecisionRecord, FragmentRecord};
pub use event::{
    Confirmed, Dropped, Located, Possible, Reason, Record, ShardRecord, Status, Summary,
    Unfinished, Window,
};
pub use paf::paf_lines;

use thiserror::Error;

/// What can go wrong while building a record.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum ReportError {
    /// A fragment names a contig that the contig table does not hold.
    #[error("contig {0} is not in the contig table of the entry")]
    UnknownContig(u32),
    /// An event names an entry that the shard does not hold.
    #[error("entry {0} is not in the shard")]
    UnknownEntry(u32),
}
