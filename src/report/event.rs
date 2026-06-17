//! The JSON line of an entry event and the closing summary line of a search.

use super::DecisionRecord;
use crate::report::ContigInfo;
use crate::report::ReportError;
use crate::search::align_stage::AlignedEntry;
use crate::search::{Completion, DropReason, EntryEvent, SearchStage, SearchSummary, Stopped};
use serde::{Deserialize, Serialize};

/// A stage of the search.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Stage {
    Sketch,
    Locate,
    Align,
}

impl From<SearchStage> for Stage {
    fn from(stage: SearchStage) -> Self {
        match stage {
            SearchStage::Sketch => Self::Sketch,
            SearchStage::Locate => Self::Locate,
            SearchStage::Align => Self::Align,
        }
    }
}

/// Why an entry left the search before it was confirmed.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Reason {
    CoverageBound,
    NoSeedHits,
    Decision,
}

impl From<DropReason> for Reason {
    fn from(reason: DropReason) -> Self {
        match reason {
            DropReason::CoverageBound => Self::CoverageBound,
            DropReason::NoSeedHits => Self::NoSeedHits,
            DropReason::Decision => Self::Decision,
        }
    }
}

/// How a search ended.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Ending {
    /// Every possible entry was decided.
    Complete,
    /// The search passed its deadline.
    Deadline,
    /// The search was asked to stop.
    Asked,
    /// Reading the index failed; the summary names the error.
    Failed,
}

/// The work of stage 1.
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
pub struct SketchWorkRecord {
    pub query_kmers: u64,
    pub lookups: u64,
    pub seeds_found: u64,
    pub frequent_seeds: u64,
    pub expanded_seeds: u64,
    pub seed_hits: u64,
    pub entries_with_seed_hits: u64,
    pub entries_without_seed_hits: u64,
    pub possible_entries: u64,
    pub seconds: f64,
}

/// The work of stage 2.
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
pub struct LocateWorkRecord {
    pub entries: u64,
    pub lookups: u64,
    pub seed_locations: u64,
    pub false_seed_hits: u64,
    pub dropped_entries: u64,
    pub candidate_bins: u64,
    pub windows: u64,
    /// Summed over threads.
    pub seconds: f64,
}

/// The work of stage 3.
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
pub struct AlignWorkRecord {
    pub entries: u64,
    pub confirmed_entries: u64,
    pub windows: u64,
    pub bases_read: u64,
    pub anchors: u64,
    pub chains: u64,
    pub alignments: u64,
    /// Summed over threads.
    pub seconds: f64,
}

/// One line of output. `kind` names the record; `query` is the name of the query record and
/// `entry` the name of the entry.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Record {
    /// Stage 1 kept the entry.
    Possible {
        query: String,
        entry: String,
        coverage_bound: u64,
        containment_estimate: f64,
    },
    /// Stage 2 found where to align the entry.
    Located {
        query: String,
        entry: String,
        coverage_bound: u64,
        containment_estimate: f64,
        /// Alignment windows and their bases together.
        windows: u64,
        window_bases: u64,
    },
    /// Stage 3 aligned the entry and it qualifies.
    Confirmed {
        query: String,
        entry: String,
        decision: DecisionRecord,
    },
    /// The entry left the search at `stage`.
    Dropped {
        query: String,
        entry: String,
        stage: Stage,
        reason: Reason,
    },
    /// The search ended before `stage` finished for the entry.
    Unfinished {
        query: String,
        entry: String,
        stage: Stage,
    },
    /// The last line of a search.
    Summary {
        query: String,
        query_length: u64,
        ending: Ending,
        /// The cause when the ending is `failed`.
        error: Option<String>,
        sketch: SketchWorkRecord,
        locate: LocateWorkRecord,
        align: AlignWorkRecord,
    },
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
    /// The record of an entry event. `contigs` gives the contig table of the entry and is only
    /// called for a confirmed entry.
    pub fn from_event<E: From<ReportError>>(
        query: &str,
        entry: &str,
        event: &EntryEvent,
        contigs: impl FnOnce() -> Result<Vec<ContigInfo>, E>,
    ) -> Result<Self, E> {
        let (query, entry) = (query.to_string(), entry.to_string());
        Ok(match event {
            EntryEvent::Possible {
                coverage_bound,
                containment_estimate,
                ..
            } => Self::Possible {
                query,
                entry,
                coverage_bound: *coverage_bound,
                containment_estimate: *containment_estimate,
            },
            EntryEvent::Located(located) => Self::Located {
                query,
                entry,
                coverage_bound: located.coverage_bound,
                containment_estimate: located.containment_estimate,
                windows: located.windows.iter().map(|c| c.ranges.len() as u64).sum(),
                window_bases: located
                    .windows
                    .iter()
                    .flat_map(|c| &c.ranges)
                    .map(|range| range.end - range.start)
                    .sum(),
            },
            EntryEvent::Confirmed(AlignedEntry { decision, .. }) => Self::Confirmed {
                query,
                entry,
                decision: DecisionRecord::new(decision.mosaic.query_length, decision, &contigs()?)?,
            },
            EntryEvent::Dropped { stage, reason, .. } => Self::Dropped {
                query,
                entry,
                stage: (*stage).into(),
                reason: (*reason).into(),
            },
            EntryEvent::Unfinished { stage, .. } => Self::Unfinished {
                query,
                entry,
                stage: (*stage).into(),
            },
        })
    }

    /// The summary line of a finished search.
    pub fn summary(query: &str, summary: &SearchSummary) -> Self {
        let (ending, error) = match &summary.completion {
            Completion::Complete => (Ending::Complete, None),
            Completion::Stopped(Stopped::Deadline) => (Ending::Deadline, None),
            Completion::Stopped(Stopped::Asked) => (Ending::Asked, None),
            Completion::Failed(message) => (Ending::Failed, Some(message.clone())),
        };
        let (sketch, locate, align) = (&summary.sketch, &summary.locate, &summary.align);
        Self::Summary {
            query: query.to_string(),
            query_length: summary.query_length,
            ending,
            error,
            sketch: SketchWorkRecord {
                query_kmers: sketch.query_kmers,
                lookups: sketch.lookups,
                seeds_found: sketch.seeds_found,
                frequent_seeds: sketch.frequent_seeds,
                expanded_seeds: sketch.expanded_seeds,
                seed_hits: sketch.seed_hits,
                entries_with_seed_hits: sketch.entries_with_seed_hits,
                entries_without_seed_hits: sketch.entries_without_seed_hits,
                possible_entries: sketch.possible_entries,
                seconds: sketch.time.as_secs_f64(),
            },
            locate: LocateWorkRecord {
                entries: locate.entries,
                lookups: locate.lookups,
                seed_locations: locate.seed_locations,
                false_seed_hits: locate.false_seed_hits,
                dropped_entries: locate.dropped_entries,
                candidate_bins: locate.candidate_bins,
                windows: locate.windows,
                seconds: locate.time.as_secs_f64(),
            },
            align: AlignWorkRecord {
                entries: align.entries,
                confirmed_entries: align.confirmed_entries,
                windows: align.counts.windows,
                bases_read: align.counts.bases_read,
                anchors: align.counts.anchors,
                chains: align.counts.chains,
                alignments: align.counts.alignments,
                seconds: align.time.as_secs_f64(),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::search::LocatedEntry;

    fn round_trip(record: &Record) {
        let text = serde_json::to_string(record).unwrap();
        assert!(!text.contains('\n'));
        assert_eq!(&serde_json::from_str::<Record>(&text).unwrap(), record);
    }

    #[test]
    fn event_records_round_trip_through_json() {
        let events = [
            EntryEvent::Possible {
                entry: 3,
                coverage_bound: 900,
                containment_estimate: 0.75,
            },
            EntryEvent::Located(LocatedEntry {
                entry: 3,
                coverage_bound: 900,
                containment_estimate: 0.75,
                windows: vec![crate::search::ContigWindows {
                    contig: 0,
                    ranges: vec![10..20, 40..45],
                }],
            }),
            EntryEvent::Dropped {
                entry: 3,
                stage: SearchStage::Align,
                reason: DropReason::Decision,
            },
            EntryEvent::Unfinished {
                entry: 3,
                stage: SearchStage::Locate,
            },
        ];
        for event in &events {
            let record =
                Record::from_event::<ReportError>("q1", "plasmid_a", event, || Ok(vec![])).unwrap();
            round_trip(&record);
        }
        let located = Record::from_event::<ReportError>("q1", "e", &events[1], || Ok(vec![]));
        assert!(matches!(
            located,
            Ok(Record::Located {
                windows: 2,
                window_bases: 15,
                ..
            })
        ));
    }

    #[test]
    fn a_record_with_its_shard_keeps_the_flat_json_line() {
        let record = Record::Dropped {
            query: "q1".into(),
            entry: "e".into(),
            stage: Stage::Align,
            reason: Reason::Decision,
        };
        let shard_record = ShardRecord {
            shard: "s1".into(),
            record,
        };
        let text = serde_json::to_string(&shard_record).unwrap();
        let json: serde_json::Value = serde_json::from_str(&text).unwrap();
        assert_eq!(
            (json["shard"].as_str(), json["kind"].as_str()),
            (Some("s1"), Some("dropped"))
        );
        assert_eq!(
            serde_json::from_str::<ShardRecord>(&text).unwrap(),
            shard_record
        );
    }

    #[test]
    fn a_summary_round_trips_and_names_the_ending() {
        let summary = SearchSummary {
            query_length: 5_000,
            completion: Completion::Stopped(Stopped::Deadline),
            sketch: Default::default(),
            locate: Default::default(),
            align: Default::default(),
        };
        let record = Record::summary("q1", &summary);
        round_trip(&record);
        let json = serde_json::to_value(&record).unwrap();
        assert_eq!(json["kind"], "summary");
        assert_eq!(json["ending"], "deadline");
    }
}
