//! The three stages of one search against one shard: sketch, locate and align.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

mod align_stage;
mod bound;
mod control;
mod driver;
mod error;
mod events;
mod in_order;
mod locate_stage;
mod located;
mod query;
mod seed_mismatches;
mod sketch_stage;
mod times;

#[cfg(test)]
pub use align_stage::align_entry;
pub use align_stage::{
    AlignStageConfig, AlignStageError, AlignStageWorkspace, AlignWork, AlignedEntry, Limited,
    PreparedQuery, entry_fragments,
};
pub use bound::CoverageBound;
pub use control::{SearchControl, Stopped};
pub use driver::{AlignStageWork, PreparedKmers, SearchConfig, SearchSummary, Searcher};
pub use error::SearchError;
pub use events::{Completion, DropReason, EntryEvent, SearchStage};
pub use in_order::in_order;
pub use locate_stage::{
    EntryLocations, LocateOutcome, LocateRequest, LocateWork, LocateWorkspace, locate_entry,
    read_entry_locations,
};
pub use located::{ContigWindows, LocatedEntry};
pub use query::QueryKmers;
pub use sketch_stage::{
    FoundSeeds, PossibleEntry, SketchLimits, SketchWork, find_query_seeds, sketch_stage,
};
pub use times::{PartTime, SearchTimes};

#[cfg(test)]
mod locate_stage_tests;
#[cfg(test)]
mod miss_simulation;
#[cfg(test)]
mod s3_reads_tests;
#[cfg(test)]
mod search_tests;
#[cfg(test)]
mod sketch_stage_tests;
#[cfg(test)]
mod warm_latency_tests;
