//! The records that result frames carry: one per entry event and one per shard summary.
//!
//! They stand in for the records of `report` until that module lands; the network and the command
//! line will then share one record type. Every record names its entry, so a client needs no list
//! of entry names.

use crate::align::{EditCounts, Strand, format_cigar};
use crate::decide::{Decision, Fragment};
use crate::index::Shard;
use crate::search::{
    Completion, DropReason, EntryEvent, LocatedEntry, SearchStage, SearchSummary, Stopped,
};
use serde::{Deserialize, Serialize};

/// A stage of the search, by name.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StageName {
    Sketch,
    Locate,
    Align,
}

impl From<SearchStage> for StageName {
    fn from(stage: SearchStage) -> Self {
        match stage {
            SearchStage::Sketch => Self::Sketch,
            SearchStage::Locate => Self::Locate,
            SearchStage::Align => Self::Align,
        }
    }
}

impl From<StageName> for SearchStage {
    fn from(stage: StageName) -> Self {
        match stage {
            StageName::Sketch => Self::Sketch,
            StageName::Locate => Self::Locate,
            StageName::Align => Self::Align,
        }
    }
}

/// One entry event of one shard.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EventRecord {
    /// The entry's index in its shard.
    pub entry: u32,
    /// The entry's name from the manifest.
    pub name: String,
    /// The new state of the entry.
    pub state: EntryState,
}

/// The state an entry event moves an entry to.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "state", rename_all = "snake_case")]
pub enum EntryState {
    /// Stage 1 kept the entry.
    Possible {
        coverage_bound: u64,
        containment_estimate: f64,
    },
    /// Stage 2 found where to align the entry. It holds all of the located entry, so stage 3 can
    /// later run on another shard.
    Located(LocatedRecord),
    /// Stage 3 aligned the entry and it qualifies.
    Confirmed(ConfirmedRecord),
    /// The entry left the search at `stage`.
    Dropped { stage: StageName, reason: DropName },
    /// The search ended before `stage` finished for the entry.
    Unfinished { stage: StageName },
}

/// Why an entry was dropped, by name.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DropName {
    CoverageBound,
    NoSeedHits,
    Decision,
}

/// A located entry: its bounds and its alignment windows.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LocatedRecord {
    pub coverage_bound: u64,
    pub containment_estimate: f64,
    /// Per contig, its alignment windows as `[start, end)` pairs.
    pub windows: Vec<(u32, Vec<(u64, u64)>)>,
}

/// A confirmed entry: its query coverage, counts and the fragments of its mosaic.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ConfirmedRecord {
    pub query_length: u64,
    pub covered_bases: u64,
    pub matches: u64,
    pub substitutions: u64,
    pub insertions: u64,
    pub deletions: u64,
    pub identity: f64,
    pub dropped_fragments: u64,
    /// Primary fragments first, best first, then the alternatives.
    pub fragments: Vec<FragmentRecord>,
}

/// One fragment of a mosaic.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FragmentRecord {
    /// Whether the fragment adds query bases no better fragment covers.
    pub primary: bool,
    /// The contig's index within its entry.
    pub contig: u32,
    /// The query bases, as `[start, end)` pairs in alignment order.
    pub query_segments: Vec<(u64, u64)>,
    pub contig_start: u64,
    pub contig_end: u64,
    /// Whether the query aligns to the reverse strand of the contig.
    pub reverse: bool,
    pub score: i32,
    pub cigar: String,
}

/// How a shard's search ended and the work it did.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SummaryRecord {
    pub query_length: u64,
    pub completion: CompletionRecord,
    pub possible_entries: u64,
    pub located_entries: u64,
    pub aligned_entries: u64,
    pub confirmed_entries: u64,
    pub seed_hits: u64,
    pub seed_locations: u64,
    pub candidate_bins: u64,
    pub bases_read: u64,
    pub alignments: u64,
    /// Wall time of stage 1, and the times of stages 2 and 3 summed over threads, in seconds.
    pub sketch_seconds: f64,
    pub locate_seconds: f64,
    pub align_seconds: f64,
}

/// How a search ended.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(tag = "completion", rename_all = "snake_case")]
pub enum CompletionRecord {
    Complete,
    /// The deadline passed (`deadline`) or the search was asked to stop (`asked`).
    Stopped {
        deadline: bool,
    },
    Failed {
        message: String,
    },
}

impl SummaryRecord {
    /// The record of a search summary.
    pub fn new(summary: &SearchSummary) -> Self {
        let completion = match &summary.completion {
            Completion::Complete => CompletionRecord::Complete,
            Completion::Stopped(stopped) => CompletionRecord::Stopped {
                deadline: *stopped == Stopped::Deadline,
            },
            Completion::Failed(message) => CompletionRecord::Failed {
                message: message.clone(),
            },
        };
        let (sketch, locate, align) = (&summary.sketch, &summary.locate, &summary.align);
        Self {
            query_length: summary.query_length,
            completion,
            possible_entries: sketch.possible_entries,
            located_entries: locate.entries.saturating_sub(locate.dropped_entries),
            aligned_entries: align.entries,
            confirmed_entries: align.confirmed_entries,
            seed_hits: sketch.seed_hits,
            seed_locations: locate.seed_locations,
            candidate_bins: locate.candidate_bins,
            bases_read: align.counts.bases_read,
            alignments: align.counts.alignments,
            sketch_seconds: sketch.time.as_secs_f64(),
            locate_seconds: locate.time.as_secs_f64(),
            align_seconds: align.time.as_secs_f64(),
        }
    }

    /// Whether every possible entry was decided.
    pub fn is_complete(&self) -> bool {
        self.completion == CompletionRecord::Complete
    }
}

impl EventRecord {
    /// The record of an event of `shard`; the entry name comes from the shard's manifest.
    pub fn new(shard: &Shard, event: EntryEvent) -> Self {
        let entry = event.entry();
        let name = shard
            .manifest()
            .entries
            .get(entry as usize)
            .map(|record| record.name.clone())
            .unwrap_or_default();
        let state = match event {
            EntryEvent::Possible {
                coverage_bound,
                containment_estimate,
                ..
            } => EntryState::Possible {
                coverage_bound,
                containment_estimate,
            },
            EntryEvent::Located(located) => EntryState::Located(LocatedRecord::new(&located)),
            EntryEvent::Confirmed(aligned) => {
                EntryState::Confirmed(ConfirmedRecord::new(&aligned.decision))
            }
            EntryEvent::Dropped { stage, reason, .. } => EntryState::Dropped {
                stage: stage.into(),
                reason: match reason {
                    DropReason::CoverageBound => DropName::CoverageBound,
                    DropReason::NoSeedHits => DropName::NoSeedHits,
                    DropReason::Decision => DropName::Decision,
                },
            },
            EntryEvent::Unfinished { stage, .. } => EntryState::Unfinished {
                stage: stage.into(),
            },
        };
        Self { entry, name, state }
    }

    /// Whether the event ends the entry: every entry gets at most one such event.
    pub fn is_final(&self) -> bool {
        !matches!(self.state, EntryState::Possible { .. })
    }
}

impl LocatedRecord {
    /// The record of a located entry.
    pub fn new(located: &LocatedEntry) -> Self {
        Self {
            coverage_bound: located.coverage_bound,
            containment_estimate: located.containment_estimate,
            windows: located
                .windows
                .iter()
                .map(|contig| {
                    let ranges = contig.ranges.iter().map(|r| (r.start, r.end)).collect();
                    (contig.contig, ranges)
                })
                .collect(),
        }
    }
}

impl ConfirmedRecord {
    /// The record of the decision of a confirmed entry.
    pub fn new(decision: &Decision) -> Self {
        let mosaic = &decision.mosaic;
        let EditCounts {
            matches,
            substitutions,
            insertions,
            deletions,
        } = decision.counts;
        let primary = mosaic.primary.iter().map(|s| (true, &s.fragment));
        let alternatives = mosaic.alternatives.iter().map(|fragment| (false, fragment));
        Self {
            query_length: mosaic.query_length,
            covered_bases: mosaic.covered_bases,
            matches,
            substitutions,
            insertions,
            deletions,
            identity: decision.identity(),
            dropped_fragments: decision.dropped_fragments as u64,
            fragments: primary
                .chain(alternatives)
                .map(|(primary, fragment)| FragmentRecord::new(primary, fragment))
                .collect(),
        }
    }
}

impl FragmentRecord {
    fn new(primary: bool, fragment: &Fragment) -> Self {
        let alignment = &fragment.alignment;
        Self {
            primary,
            contig: fragment.contig_index,
            query_segments: fragment
                .query_segments
                .iter()
                .map(|segment| (segment.start, segment.end))
                .collect(),
            contig_start: alignment.contig_interval.start,
            contig_end: alignment.contig_interval.end,
            reverse: alignment.strand == Strand::Reverse,
            score: alignment.score,
            cigar: format_cigar(&alignment.runs),
        }
    }
}
