//! The JSON line of an entry event and the closing summary line of a search.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{ContigInfo, DecisionRecord, ReportError};
use crate::search::align_stage::AlignedEntry;
use crate::search::{Completion, DropReason, EntryEvent, SearchStage, SearchSummary, Stopped};
use serde::{Deserialize, Serialize};

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
    pub seeded_entries: u64,
    pub unseeded_entries: u64,
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

/// Stage 1 kept the entry.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct Possible {
    pub query: String,
    pub entry: String,
    pub coverage_bound: u64,
    pub containment_estimate: f64,
}

/// Stage 2 found where to align the entry.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct Located {
    pub query: String,
    pub entry: String,
    pub coverage_bound: u64,
    pub containment_estimate: f64,
    /// Alignment windows and their bases together.
    pub windows: u64,
    pub window_bases: u64,
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
    pub ending: Ending,
    /// The cause when the ending is `failed`.
    pub error: Option<String>,
    pub sketch: SketchWorkRecord,
    pub locate: LocateWorkRecord,
    pub align: AlignWorkRecord,
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
    Summary(Summary),
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
            } => Self::Possible(Possible {
                query,
                entry,
                coverage_bound: *coverage_bound,
                containment_estimate: *containment_estimate,
            }),
            EntryEvent::Located(located) => Self::Located(Located {
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
            }),
            EntryEvent::Confirmed(AlignedEntry { decision, .. }) => Self::Confirmed(Confirmed {
                query,
                entry,
                decision: DecisionRecord::new(decision.mosaic.query_length, decision, &contigs()?)?,
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

    /// The summary line of a finished search.
    pub fn summary(query: &str, summary: &SearchSummary) -> Self {
        let (ending, error) = match &summary.completion {
            Completion::Complete => (Ending::Complete, None),
            Completion::Stopped(Stopped::Deadline) => (Ending::Deadline, None),
            Completion::Stopped(Stopped::Asked) => (Ending::Asked, None),
            Completion::Failed(message) => (Ending::Failed, Some(message.clone())),
        };
        let (sketch, locate, align) = (&summary.sketch, &summary.locate, &summary.align);
        Self::Summary(Summary {
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
                seeded_entries: sketch.seeded_entries,
                unseeded_entries: sketch.unseeded_entries,
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
        })
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
    fn events_round_trip() {
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
            Ok(Record::Located(Located {
                windows: 2,
                window_bases: 15,
                ..
            }))
        ));
    }

    #[test]
    fn shard_record_flat() {
        let record = Record::Dropped(Dropped {
            query: "q1".into(),
            entry: "e".into(),
            stage: SearchStage::Align,
            reason: DropReason::Decision,
        });
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
    fn summary_names_ending() {
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
