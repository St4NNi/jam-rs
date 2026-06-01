//! What a search reports while it runs: one event each time an entry changes its state.

use super::LocatedEntry;
use super::control::Stopped;

/// A stage of the search. As a search level it is the last stage that runs.
#[derive(Clone, Copy, Debug, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum SearchStage {
    Sketch,
    Locate,
    Align,
}

/// Why an entry left the search before it was confirmed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DropReason {
    /// The coverage bound fell below the requested query coverage.
    CoverageBound,
    /// The positions hold none of the entry's seed hits, so there is nothing to align.
    NoSeedHits,
}

/// A change of one entry's state. Every possible entry ends with exactly one of `Located` (when
/// the search stops after stage 2), `Dropped` or `Unfinished`, and its events arrive in order.
#[derive(Clone, Debug, PartialEq)]
pub enum EntryEvent {
    /// Stage 1 kept the entry.
    Possible {
        entry: u32,
        coverage_bound: u64,
        containment_estimate: f64,
    },
    /// Stage 2 found where to align the entry.
    Located(LocatedEntry),
    /// The entry left the search at `stage`.
    Dropped {
        entry: u32,
        stage: SearchStage,
        reason: DropReason,
    },
    /// The search ended before `stage` finished for the entry, so the entry is not decided.
    Unfinished { entry: u32, stage: SearchStage },
}

impl EntryEvent {
    /// The entry the event is about.
    pub fn entry(&self) -> u32 {
        match self {
            Self::Possible { entry, .. }
            | Self::Dropped { entry, .. }
            | Self::Unfinished { entry, .. } => *entry,
            Self::Located(located) => located.entry,
        }
    }
}

/// How a search ended. Only `Complete` means every possible entry was decided.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Completion {
    Complete,
    Stopped(Stopped),
    /// Reading the index failed; the message names the cause.
    Failed(String),
}
