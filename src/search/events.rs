//! What a search reports while it runs: one event each time an entry changes its state.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::LocatedEntry;
use super::align_stage::{AlignStageError, AlignedEntry};
use super::control::Stopped;
use super::error::SearchError;
use crate::index::Contig;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// A stage of the search; as the last stage of a search it is its level.
#[derive(
    Clone,
    Copy,
    Debug,
    Deserialize,
    Eq,
    Hash,
    Ord,
    PartialEq,
    PartialOrd,
    Serialize,
    clap::ValueEnum,
)]
#[serde(rename_all = "snake_case")]
pub enum SearchStage {
    /// Stage 1: the entries that may hold the query
    Sketch,
    /// Stage 2: where in those entries the query may lie
    Locate,
    /// Stage 3: aligned and decided entries
    Align,
}

/// Why an entry left the search before it was confirmed.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum DropReason {
    /// The coverage bound fell below the requested query coverage.
    CoverageBound,
    /// The locator part holds none of the entry's seed hits, so there is nothing to align.
    NoSeedHits,
    /// The alignments of the entry do not reach the decision rule.
    Decision,
    /// Aligning the entry needs more memory or work than a limit of stage 3 allows, as when a
    /// periodic query meets many periodic regions. The entry is not decided, so the search ends
    /// `ResourceIncomplete`.
    ResourceLimit,
}

/// A change of one entry's state; the events of an entry arrive in order. Past stage 1, every
/// possible entry ends with one of `Located` (at the locate level), `Confirmed`, `Dropped` or
/// `Unfinished`.
#[derive(Clone, Debug, PartialEq)]
pub enum EntryEvent {
    /// Stage 1 kept the entry.
    Possible {
        entry: u32,
        coverage_bound: u64,
        containment_estimate: f64,
    },
    /// Stage 2 found where to align the entry. At the locate level the map holds the contigs of
    /// the alignment windows, by contig index; else it is empty.
    Located(LocatedEntry, BTreeMap<u32, Contig>),
    /// Stage 3 aligned the entry and it qualifies; the decision holds its mosaic, and `contigs`
    /// the contigs that its fragments refer to, by contig index.
    Confirmed {
        aligned: Box<AlignedEntry>,
        contigs: BTreeMap<u32, Contig>,
    },
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
            Self::Located(located, _) => located.entry,
            Self::Confirmed { aligned, .. } => aligned.entry,
        }
    }
}

/// How a search ended. Only `Complete` means every possible entry was decided.
#[derive(Debug)]
pub enum Completion {
    Complete,
    /// Every other entry was decided, but at least one was dropped for a resource limit, so a
    /// qualifying entry may be missing.
    ResourceIncomplete,
    Stopped(Stopped),
    /// A stage failed for one entry, for example reading the index or aligning; the entries not
    /// yet decided are unfinished.
    Failed(SearchError),
}

impl From<SearchError> for Completion {
    fn from(error: SearchError) -> Self {
        match error {
            SearchError::Stopped(stopped)
            | SearchError::Align(AlignStageError::Stopped(stopped)) => Self::Stopped(stopped),
            error => Self::Failed(error),
        }
    }
}
