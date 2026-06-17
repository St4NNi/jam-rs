//! End-to-end tests of all three stages on small generated shards.

use super::control::SearchControl;
use super::driver::{SearchConfig, SearchSummary, Searcher};
use super::events::{Completion, DropReason, EntryEvent, SearchStage};
use super::test_shard::{Random, build, model, reverse_complement};
use crate::align::{Strand, format_cigar};
use crate::decide::{Fragment, QueryCoverage};
use crate::index::Shard;
use std::collections::BTreeMap;
use std::sync::atomic::AtomicBool;
use std::sync::mpsc::sync_channel;

const QUERY_LENGTH: usize = 5_000;

fn config() -> SearchConfig {
    let mut config = SearchConfig {
        threads: 3,
        ..SearchConfig::default()
    };
    config.rule.min_identity = 0.95;
    config.rule.min_query_coverage = QueryCoverage::Share {
        numerator: 9,
        decimal_places: 1,
    };
    config
}

fn search(
    shard: &Shard,
    query: &[u8],
    circular: bool,
    config: &SearchConfig,
) -> (SearchSummary, BTreeMap<u32, Vec<EntryEvent>>) {
    let searcher = Searcher::new(shard, model()).unwrap();
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let (sender, receiver) = sync_channel(16);
    let mut events = BTreeMap::<u32, Vec<EntryEvent>>::new();
    let summary = std::thread::scope(|scope| {
        let search = scope.spawn(|| {
            let summary = searcher.search(query, circular, config, control, &sender);
            drop(sender);
            summary
        });
        for event in receiver {
            events.entry(event.entry()).or_default().push(event);
        }
        search.join().unwrap().unwrap()
    });
    assert_eq!(summary.completion, Completion::Complete);
    (summary, events)
}

/// Checks that the events of an entry run possible, located, then the last one, and returns it.
fn last_event(events: &BTreeMap<u32, Vec<EntryEvent>>, entry: u32) -> &EntryEvent {
    let list = &events[&entry];
    assert_eq!(list.len(), 3, "entry {entry}: {list:?}");
    assert!(matches!(list[0], EntryEvent::Possible { .. }));
    assert!(matches!(list[1], EntryEvent::Located(_)));
    &list[2]
}

/// The fragments of a confirmed entry's mosaic.
fn fragments(event: &EntryEvent) -> Vec<&Fragment> {
    let EntryEvent::Confirmed(aligned) = event else {
        panic!("not confirmed: {event:?}");
    };
    assert!(aligned.decision.qualifies);
    aligned
        .decision
        .mosaic
        .primary
        .iter()
        .map(|selected| &selected.fragment)
        .collect()
}

/// Entry 0: the query with three substitutions in its second contig, at 10,000. Entry 1: the
/// reverse complement of the query at 8,000. Entry 2: unrelated.
fn planted_shard(directory: &std::path::Path) -> (Shard, Vec<u8>) {
    let mut random = Random(0xa119);
    let query = random.bases(QUERY_LENGTH);
    let mut changed = query.clone();
    for at in [1_000, 2_000, 3_000] {
        changed[at] = if changed[at] == b'A' { b'C' } else { b'A' };
    }
    let entries = vec![
        vec![
            random.bases(3_000),
            [random.bases(10_000), changed, random.bases(10_000)].concat(),
        ],
        vec![
            [
                random.bases(8_000),
                reverse_complement(&query),
                random.bases(4_000),
            ]
            .concat(),
        ],
        vec![random.bases(20_000)],
    ];
    (build(directory, &entries, None), query)
}

#[test]
fn a_planted_query_is_confirmed_with_its_place_strand_and_cigar() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = planted_shard(directory.path());
    let (summary, events) = search(&shard, &query, false, &config());
    // The unrelated entry never reaches stage 3.
    assert_eq!(summary.align.entries, 2);
    assert_eq!(summary.align.confirmed_entries, 2);
    assert!(summary.align.counts.alignments >= 2);

    let forward = fragments(last_event(&events, 0));
    assert_eq!(forward.len(), 1);
    let fragment = forward[0];
    assert_eq!(fragment.contig_index, 1);
    assert_eq!(fragment.alignment.strand, Strand::Forward);
    let contig = fragment.alignment.contig_interval;
    assert_eq!((contig.start, contig.end), (10_000, 15_000));
    let query_interval = fragment.alignment.query_interval;
    assert_eq!((query_interval.start, query_interval.end), (0, 5_000));
    assert_eq!(
        format_cigar(&fragment.alignment.runs),
        "1000=1X999=1X999=1X1999="
    );
    let EntryEvent::Confirmed(aligned) = last_event(&events, 0) else {
        unreachable!()
    };
    assert!((aligned.decision.identity() - 4_997.0 / 5_000.0).abs() < 1e-12);

    let reverse = fragments(last_event(&events, 1));
    assert_eq!(reverse.len(), 1);
    assert_eq!(reverse[0].contig_index, 0);
    assert_eq!(reverse[0].alignment.strand, Strand::Reverse);
    let contig = reverse[0].alignment.contig_interval;
    assert_eq!((contig.start, contig.end), (8_000, 13_000));
    assert_eq!(format_cigar(&reverse[0].alignment.runs), "5000=");
}

#[test]
fn an_entry_below_the_decision_rule_is_dropped_in_stage_3() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = planted_shard(directory.path());
    // The query covers a fifth of the contig of entry 0 and a third of that of entry 1.
    let mut config = config();
    config.rule.min_contig_coverage = 0.5;
    let (summary, events) = search(&shard, &query, false, &config);
    for entry in [0, 1] {
        assert_eq!(
            *last_event(&events, entry),
            EntryEvent::Dropped {
                entry,
                stage: SearchStage::Align,
                reason: DropReason::Decision,
            }
        );
    }
    assert_eq!(summary.align.entries, 2);
    assert_eq!(summary.align.confirmed_entries, 0);
}

#[test]
fn a_circular_query_is_confirmed_over_its_end() {
    let directory = tempfile::tempdir().unwrap();
    let mut random = Random(0xc1c);
    let query = random.bases(QUERY_LENGTH);
    let rotated = [&query[1_700..], &query[..1_700]].concat();
    let entries = vec![
        vec![[random.bases(6_000), rotated, random.bases(6_000)].concat()],
        vec![random.bases(15_000)],
    ];
    let shard = build(directory.path(), &entries, None);
    let (summary, events) = search(&shard, &query, true, &config());
    assert_eq!(summary.align.confirmed_entries, 1);
    let EntryEvent::Confirmed(aligned) = last_event(&events, 0) else {
        panic!("{:?}", events[&0]);
    };
    assert_eq!(aligned.decision.mosaic.covered_bases, QUERY_LENGTH as u64);
    assert!((aligned.decision.identity() - 1.0).abs() < 1e-12);
    let over_end = aligned
        .decision
        .mosaic
        .primary
        .iter()
        .any(|selected| selected.fragment.query_segments.len() == 2);
    assert!(over_end, "{:?}", aligned.decision.mosaic.primary);
    // The unrelated entry never reaches stage 3.
    assert_eq!(summary.align.entries, 1);
}
