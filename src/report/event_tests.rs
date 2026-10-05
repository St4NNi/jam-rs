//! JSON round trips of the event and summary records.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::search::LocatedEntry;

fn round_trip(record: &Record) {
    let text = serde_json::to_string(record).unwrap();
    assert!(!text.contains('\n'));
    assert_eq!(&serde_json::from_str::<Record>(&text).unwrap(), record);
}

fn contig_a() -> Contig {
    Contig {
        name: "contig_a".into(),
        length: 100,
    }
}

#[test]
fn events_round_trip() {
    let events = [
        EntryEvent::Possible {
            entry: 3,
            coverage_bound: 900,
            containment_estimate: 0.75,
        },
        EntryEvent::Located(
            LocatedEntry {
                entry: 3,
                coverage_bound: 900,
                containment_estimate: 0.75,
                alignment_windows: vec![crate::search::ContigWindows {
                    contig: 0,
                    ranges: vec![10..20, 40..45],
                }],
            },
            BTreeMap::from([(0, contig_a())]),
        ),
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
        let record = Record::from_event("q1", "plasmid_a", event).unwrap();
        round_trip(&record);
    }
    let located = Record::from_event("q1", "e", &events[1]);
    assert!(matches!(
        located,
        Ok(Record::Located(Located {
            alignment_windows: 2,
            window_bases: 15,
            ..
        }))
    ));
    let Ok(Record::Located(located)) = located else {
        panic!("not located");
    };
    let window = |start, end| Window {
        contig: "contig_a".into(),
        start,
        end,
    };
    assert_eq!(located.windows, [window(10, 20), window(40, 45)]);
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
fn summary_names_status() {
    let summary = SearchSummary {
        query_length: 5_000,
        seed_mismatches: 1,
        completion: Completion::Stopped(Stopped::Deadline),
        sketch: Default::default(),
        locate: Default::default(),
        align: Default::default(),
        times: Default::default(),
    };
    let record = Record::summary("q1", &summary);
    round_trip(&record);
    let json = serde_json::to_value(&record).unwrap();
    assert_eq!(json["kind"], "summary");
    assert_eq!(json["status"], "stopped");
    assert_eq!(json["reason"], "deadline");
    assert_eq!(json["seed_mismatches"], 1);
}

/// A resource drop ends the search as resource-incomplete, never complete.
#[test]
fn resource_drop_named() {
    let summary = SearchSummary {
        query_length: 1_000,
        seed_mismatches: 0,
        completion: Completion::ResourceIncomplete,
        sketch: Default::default(),
        locate: Default::default(),
        align: Default::default(),
        times: Default::default(),
    };
    let record = Record::summary("q1", &summary);
    round_trip(&record);
    let json = serde_json::to_value(&record).unwrap();
    assert_eq!(json["status"], "incomplete");
    assert_eq!(json["reason"], "resource_limit");
}
