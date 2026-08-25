//! The three stages of one search against one shard: sketch, locate and align.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

pub mod align_stage;
pub mod bound;
pub mod control;
pub mod driver;
pub mod error;
pub mod events;
pub mod locate_stage;
pub mod located;
pub mod query;
pub mod sketch_stage;

pub use control::{SearchControl, Stopped};
pub use driver::{AlignStageWork, SearchConfig, SearchSummary, Searcher};
pub use error::SearchError;
pub use events::{Completion, DropReason, EntryEvent, SearchStage};
pub use locate_stage::LocateWork;
pub use located::{ContigWindows, LocatedEntry};
pub use sketch_stage::SketchWork;

#[cfg(test)]
mod locate_stage_tests;
#[cfg(test)]
mod search_tests;
#[cfg(test)]
mod sketch_stage_tests;
