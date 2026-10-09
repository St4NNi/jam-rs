//! The JSON line of an entry event and the closing summary line of a search.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{DecisionRecord, ReportError};
use crate::index::Contig;
use crate::search::{
    AlignStageWork, Completion, ContigWindows, DropReason, EntryEvent, LocateWork, SearchStage,
    SearchSummary, SearchTimes, SketchWork, Stopped,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// How a search ended.
#[derive(Clone, Copy, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Status {
    /// Every possible entry was decided.
    Complete,
    /// Every other entry was decided, but at least one was dropped for a resource limit
    /// (`align.limited_entries`), so a qualifying entry may be missing.
    Incomplete,
    /// The search stopped early, at its deadline or when asked.
    Stopped,
    /// A stage, or the server, failed; the summary names the error.
    Failed,
}

/// Why a search did not complete.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Reason {
    /// The search passed its deadline.
    Deadline,
    /// The search was asked to stop.
    Asked,
    /// An entry needed more memory or work than a limit allows.
    ResourceLimit,
    /// A stage, or the server, failed.
    Error,
}

/// Stage 1 kept the entry.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct Possible {
    pub query: String,
    pub entry: String,
    pub coverage_bound: u64,
    /// The share of the query held by the entry, 0 to 1, from the seed hits of the sketch (frequent
    /// seeds count as hit in their color sets); mutations lower it.
    pub containment_estimate: f64,
}

/// Stage 2 found where to align the entry.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct Located {
    pub query: String,
    pub entry: String,
    pub coverage_bound: u64,
    /// The share of the query held by the entry, 0 to 1, from the seed hits with seed locations,
    /// so false seed hits no longer count; mutations lower it.
    pub containment_estimate: f64,
    /// Alignment windows and their bases together.
    pub alignment_windows: u64,
    pub window_bases: u64,
    /// The windows on named contigs; only at the locate level, else empty.
    pub windows: Vec<Window>,
}

/// An alignment window on a named contig, in contig coordinates.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct Window {
    pub contig: String,
    pub start: u64,
    pub end: u64,
}

/// Stage 3 aligned the entry and it qualifies.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct Confirmed {
    pub query: String,
    pub entry: String,
    pub decision: DecisionRecord,
}

/// The entry left the search at `stage`.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct Dropped {
    pub query: String,
    pub entry: String,
    pub stage: SearchStage,
    pub reason: DropReason,
}

/// The search ended before `stage` finished for the entry.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct Unfinished {
    pub query: String,
    pub entry: String,
    pub stage: SearchStage,
}

/// The last line of a search.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct Summary {
    pub query: String,
    pub query_length: u64,
    /// The lookup level: 0 exact, 1 one mismatch, 2 two mismatches, 3 also one indel.
    pub seed_mismatches: u8,
    pub status: Status,
    /// Why the search did not complete; `None` when it did.
    pub reason: Option<Reason>,
    /// The cause when the status is `failed`.
    pub error: Option<String>,
    pub sketch: SketchWork,
    pub locate: LocateWork,
    pub align: AlignStageWork,
    pub times: SearchTimes,
}

/// One line of output. `kind` names the record, `query` the query record and `entry` the entry.
/// It holds names and numbers only, so another process can read it back with serde.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Record {
    Possible(Possible),
    Located(Located),
    Confirmed(Confirmed),
    Dropped(Dropped),
    Unfinished(Unfinished),
    /// Boxed, as it is much larger than the other records.
    Summary(Box<Summary>),
}

/// A record with the shard it comes from, as `jam search` prints it. Entry names may repeat
/// across shards, so every line names its shard.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct ShardRecord {
    /// The name of the shard directory.
    pub shard: String,
    #[serde(flatten)]
    pub record: Record,
}

impl Record {
    /// The record of an entry event of query `query` in entry `entry`.
    pub fn from_event(query: &str, entry: &str, event: &EntryEvent) -> Result<Self, ReportError> {
        let (query, entry) = (query.to_string(), entry.to_string());
        Ok(match event {
            EntryEvent::Possible {
                coverage_bound,
                containment_estimate,
                ..
            } => Self::Possible(Possible {
                query,
                entry,
                coverage_bound: *coverage_bound,
                containment_estimate: *containment_estimate,
            }),
            EntryEvent::Located(located, contigs) => Self::Located(Located {
                query,
                entry,
                coverage_bound: located.coverage_bound,
                containment_estimate: located.containment_estimate,
                alignment_windows: located
                    .alignment_windows
                    .iter()
                    .map(|c| c.ranges.len() as u64)
                    .sum(),
                window_bases: located
                    .alignment_windows
                    .iter()
                    .flat_map(|c| &c.ranges)
                    .map(|range| range.end - range.start)
                    .sum(),
                windows: named_windows(&located.alignment_windows, contigs),
            }),
            EntryEvent::Confirmed { aligned, contigs } => Self::Confirmed(Confirmed {
                query,
                entry,
                decision: DecisionRecord::new(&aligned.decision, contigs)?,
            }),
            EntryEvent::Dropped { stage, reason, .. } => Self::Dropped(Dropped {
                query,
                entry,
                stage: *stage,
                reason: *reason,
            }),
            EntryEvent::Unfinished { stage, .. } => Self::Unfinished(Unfinished {
                query,
                entry,
                stage: *stage,
            }),
        })
    }

    /// The summary line of a search that failed before it ran, or whose server failed.
    pub fn failed(query: &str, query_length: u64, error: &dyn std::error::Error) -> Self {
        Self::Summary(Box::new(Summary {
            query: query.to_string(),
            query_length,
            seed_mismatches: 0,
            status: Status::Failed,
            reason: Some(Reason::Error),
            error: Some(error.to_string()),
            sketch: SketchWork::default(),
            locate: LocateWork::default(),
            align: AlignStageWork::default(),
            times: SearchTimes::default(),
        }))
    }

    /// The summary line of a finished search.
    pub fn summary(query: &str, summary: &SearchSummary) -> Self {
        let (status, reason, error) = match &summary.completion {
            Completion::Complete => (Status::Complete, None, None),
            Completion::ResourceIncomplete => {
                (Status::Incomplete, Some(Reason::ResourceLimit), None)
            }
            Completion::Stopped(Stopped::Deadline) => {
                (Status::Stopped, Some(Reason::Deadline), None)
            }
            Completion::Stopped(Stopped::Asked) => (Status::Stopped, Some(Reason::Asked), None),
            Completion::Failed(error) => {
                (Status::Failed, Some(Reason::Error), Some(error.to_string()))
            }
        };
        Self::Summary(Box::new(Summary {
            query: query.to_string(),
            query_length: summary.query_length,
            seed_mismatches: summary.seed_mismatches,
            status,
            reason,
            error,
            sketch: summary.sketch,
            locate: summary.locate,
            align: summary.align,
            times: summary.times,
        }))
    }
}

/// The windows of the contigs that `contigs` names.
fn named_windows(windows: &[ContigWindows], contigs: &BTreeMap<u32, Contig>) -> Vec<Window> {
    let named = windows
        .iter()
        .filter_map(|windows| Some((contigs.get(&windows.contig)?, &windows.ranges)));
    named
        .flat_map(|(contig, ranges)| {
            ranges.iter().map(|range| Window {
                contig: contig.name.clone(),
                start: range.start,
                end: range.end,
            })
        })
        .collect()
}

#[cfg(test)]
#[path = "event_tests.rs"]
mod tests;
