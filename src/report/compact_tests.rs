//! Tests of the compact result and summary lines.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::align::Interval;
use crate::report::{DecisionRecord, Located, Possible, Window};
use serde_json::{Value, json};

fn line(record: &Record, level: SearchStage) -> Option<Value> {
    let text = result_line(record, level).unwrap()?;
    Some(serde_json::from_str(&text).unwrap())
}

#[test]
fn sketch_line() {
    let possible = Record::Possible(Possible {
        query: "q1".into(),
        entry: "a1".into(),
        coverage_bound: 900,
        containment_estimate: 0.834_56,
    });
    let expected = json!({"q": "q1", "a": "a1", "s": 1, "c": 0.8346});
    assert_eq!(line(&possible, SearchStage::Sketch), Some(expected));
    assert_eq!(line(&possible, SearchStage::Align), None);
}

#[test]
fn locate_line() {
    let window = |start, end| Window {
        contig: "c1".into(),
        start,
        end,
    };
    let located = Record::Located(Located {
        query: "q1".into(),
        entry: "a1".into(),
        coverage_bound: 900,
        containment_estimate: 0.81,
        alignment_windows: 2,
        window_bases: 15,
        windows: vec![window(10, 20), window(40, 45)],
    });
    let expected = json!({"q": "q1", "a": "a1", "s": 2, "c": 0.81,
        "w": [["c1", 10, 20], ["c1", 40, 45]]});
    assert_eq!(line(&located, SearchStage::Locate), Some(expected));
    assert_eq!(line(&located, SearchStage::Align), None);
}

fn fragment(strand: Strand, start: u64) -> FragmentRecord {
    FragmentRecord {
        contig: "c1".into(),
        contig_length: 5_000,
        query_segments: vec![Interval::new(start, start + 8).unwrap()],
        contig_interval: Interval::new(100, 109).unwrap(),
        strand,
        score: 20,
        new_bases: 8,
        matches: 7,
        substitutions: 1,
        insertions: 0,
        deletions: 1,
        identity: 7.0 / 9.0,
        cigar: "3=1X3=1D1=".into(),
    }
}

#[test]
fn align_line() {
    let decision = DecisionRecord {
        query_length: 10,
        covered_bases: 8,
        query_coverage: 0.8,
        identity: 7.0 / 9.0,
        matches: 7,
        substitutions: 1,
        insertions: 0,
        deletions: 1,
        dropped_fragments: 0,
        fragments: vec![fragment(Strand::Forward, 2)],
        alternatives: vec![fragment(Strand::Reverse, 0)],
    };
    let confirmed = Record::Confirmed(Confirmed {
        query: "q1".into(),
        entry: "a1".into(),
        decision,
    });
    let expected = json!({"q": "q1", "a": "a1", "s": 3, "cov": 0.8, "id": 0.7778, "bases": 8,
        "f": [["c1", 2, 10, 100, 109, "+", "3=1X3=1D1="]],
        "alt": [["c1", 0, 8, 100, 109, "-", "3=1X3=1D1="]]});
    assert_eq!(line(&confirmed, SearchStage::Align), Some(expected));
    let Record::Confirmed(confirmed) = confirmed else {
        unreachable!()
    };
    assert_eq!(tsv_line(&confirmed), "q1\ta1\t0.8\t0.7778\t8\t1");
}

fn summary(status: Status, reason: Option<Reason>) -> Summary {
    let mut summary = Summary {
        query: "q1".into(),
        query_length: 1_000,
        seed_mismatches: 0,
        status,
        reason,
        error: None,
        sketch: Default::default(),
        locate: Default::default(),
        align: Default::default(),
        times: Default::default(),
    };
    summary.sketch.seeded_entries = 5;
    summary.sketch.possible_entries = 3;
    summary.sketch.total_miss = 0.000_01;
    summary.locate.entries = 3;
    summary.locate.dropped_entries = 1;
    summary.locate.seconds = 0.25;
    summary.align.confirmed_entries = 1;
    summary.times.sketch.wall_seconds = 0.5;
    summary
}

#[test]
fn summary_over_shards() {
    let mut totals = QueryTotals::new("q1");
    totals.add(&summary(Status::Complete, None));
    let complete: Value = serde_json::from_str(&totals.line().unwrap()).unwrap();
    let expected = json!({"q": "q1", "status": "complete", "n": [5, 3, 2, 1],
        "miss": 0.000_01, "t": [0.5, 0.25, 0.0]});
    assert_eq!(complete, expected);
    let mut unpruned = summary(Status::Stopped, Some(Reason::Deadline));
    unpruned.sketch.miss_over_budget = true;
    unpruned.sketch.total_miss = 0.02;
    totals.add(&unpruned);
    totals.add(&summary(Status::Incomplete, Some(Reason::ResourceLimit)));
    assert_eq!(totals.status(), Status::Stopped);
    let stopped: Value = serde_json::from_str(&totals.line().unwrap()).unwrap();
    assert_eq!(stopped["status"], "stopped");
    assert_eq!(stopped["reason"], "deadline");
    assert_eq!(stopped["n"], json!([15, 9, 6, 3]));
    assert_eq!(
        (stopped["miss"].as_f64(), stopped["prune"].as_bool()),
        (Some(0.02), Some(false))
    );
    assert_eq!(
        totals.progress(),
        "q1: stopped (deadline); 15 seeded, 9 possible, 6 located, 3 confirmed; \
         sketch 1.500 s, locate 0.750 s, align 0.000 s; pruning off"
    );
}
