//! End-to-end tests of the search on small generated shards.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::control::{SearchControl, Stopped};
use super::driver::{SearchConfig, SearchSummary, Searcher};
use super::events::{Completion, DropReason, EntryEvent, SearchStage};
use crate::align::{Strand, format_cigar};
use crate::decide::{Fragment, QueryCoverage};
use crate::index::{BuildOptions, EntryInput, Shard, build_shard};
use crate::seeds::reverse_complement;
use crate::test_support::{
    build_options, diverged, miss_model, no_stop, one_mismatch_model, random_bases, rng, test_shard,
};
use rand::RngExt;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::sync_channel;

const QUERY_LENGTH: usize = 5_000;

/// 90% query coverage at 95% identity, on three threads.
fn config(last_stage: SearchStage) -> SearchConfig {
    let mut config = SearchConfig {
        last_stage,
        threads: 3,
        flank: 1_000,
        ..SearchConfig::default()
    };
    config.rule.min_identity = 0.95;
    config.rule.min_fragment_identity = 0.0;
    config.rule.min_query_coverage = QueryCoverage::Share {
        numerator: 9,
        decimal_places: 1,
    };
    config
}

/// Runs a search and returns the summary and the events of each entry, in order. With
/// `stop_early` the stop flag is set at the first event.
fn search(
    shard: &Shard,
    query: &[u8],
    circular: bool,
    config: &SearchConfig,
    stop_early: bool,
) -> (SearchSummary, BTreeMap<u32, Vec<EntryEvent>>) {
    let searcher = Searcher::new(shard, miss_model(), one_mismatch_model()).unwrap();
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    // Without a buffer each event is handed over only when it is read, so a stop set at the first
    // event is seen before stage 2 starts.
    let (sender, receiver) = sync_channel(if stop_early { 0 } else { 16 });
    let mut events = BTreeMap::<u32, Vec<EntryEvent>>::new();
    let summary = std::thread::scope(|scope| {
        let search = scope.spawn(|| {
            let summary = searcher.search(query, circular, config, control, &sender);
            drop(sender);
            summary
        });
        for event in receiver {
            if stop_early {
                stop.store(true, Ordering::Relaxed);
            }
            events.entry(event.entry()).or_default().push(event);
        }
        search.join().unwrap().unwrap()
    });
    (summary, events)
}

/// Entries 0 to 3 hold the query at 97% identity, entries 4 to 7 are unrelated.
fn small_shard(directory: &std::path::Path) -> (Shard, Vec<u8>) {
    let mut random = rng(0xd1e);
    let query = random_bases(&mut random, QUERY_LENGTH);
    let mut entries = Vec::new();
    for _ in 0..4 {
        let copy = diverged(&mut random, &query, 0.97);
        entries.push(vec![
            [
                random_bases(&mut random, 2_000),
                copy,
                random_bases(&mut random, 2_000),
            ]
            .concat(),
        ]);
    }
    for _ in 0..4 {
        entries.push(vec![random_bases(&mut random, 9_000)]);
    }
    (test_shard(directory, &entries), query)
}

#[test]
fn level_ends_search() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = small_shard(directory.path());
    let (summary, events) = search(&shard, &query, false, &config(SearchStage::Sketch), false);
    assert!(matches!(summary.completion, Completion::Complete));
    assert_eq!(summary.locate.entries, 0);
    assert_eq!(events.keys().copied().collect::<Vec<_>>(), [0, 1, 2, 3]);
    assert!(
        events
            .values()
            .all(|list| list.len() == 1 && matches!(list[0], EntryEvent::Possible { .. }))
    );
    let (summary, events) = search(&shard, &query, false, &config(SearchStage::Locate), false);
    assert!(matches!(summary.completion, Completion::Complete));
    assert_eq!(summary.locate.entries, summary.sketch.possible_entries);
    for (entry, list) in &events {
        assert_eq!(list.len(), 2, "{list:?}");
        assert!(matches!(list[0], EntryEvent::Possible { .. }));
        let EntryEvent::Located(located) = &list[1] else {
            panic!("entry {entry}: {list:?}");
        };
        assert_eq!(located.alignment_windows.len(), 1);
        assert!(located.coverage_bound >= 4_500);
    }
}

#[test]
fn sketch_only_answers() {
    let directory = tempfile::tempdir().unwrap();
    let (_, query) = small_shard(directory.path());
    let output = directory.path().join("sketch");
    let options = BuildOptions {
        sketch_only: true,
        ..build_options(directory.path(), 2)
    };
    let inputs: Vec<EntryInput> = (0..8)
        .map(|entry| EntryInput::from_path(directory.path().join(format!("entry{entry}.fa"))))
        .collect();
    build_shard(&inputs, &output, &options).unwrap();
    let shard = Shard::open(&output).unwrap();
    let (summary, events) = search(&shard, &query, false, &config(SearchStage::Sketch), false);
    assert!(matches!(summary.completion, Completion::Complete));
    assert_eq!(events.keys().copied().collect::<Vec<_>>(), [0, 1, 2, 3]);
    let searcher = Searcher::new(&shard, miss_model(), one_mismatch_model()).unwrap();
    let (sender, _receiver) = sync_channel(1);
    let config = config(SearchStage::Locate);
    assert!(
        searcher
            .search(&query, false, &config, no_stop(), &sender)
            .is_err()
    );
}

#[test]
fn stop_marks_unfinished() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = small_shard(directory.path());
    let (summary, events) = search(&shard, &query, false, &config(SearchStage::Locate), true);
    assert!(matches!(
        summary.completion,
        Completion::Stopped(Stopped::Asked)
    ));
    assert_eq!(events.len(), 4);
    for list in events.values() {
        assert_eq!(list.len(), 2, "{list:?}");
        assert!(matches!(list[0], EntryEvent::Possible { .. }));
        assert!(matches!(
            list[1],
            EntryEvent::Unfinished {
                stage: SearchStage::Locate,
                ..
            }
        ));
    }
}

#[test]
fn bad_input_fails() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, _) = small_shard(directory.path());
    let searcher = Searcher::new(&shard, miss_model(), one_mismatch_model()).unwrap();
    let control = no_stop();
    let config = config(SearchStage::Locate);
    let (sender, _receiver) = sync_channel(1);
    assert!(
        searcher
            .search(b"ACGT", false, &config, control, &sender)
            .is_err()
    );
    let no_flank = SearchConfig { flank: 0, ..config };
    let query = vec![b'A'; 100];
    assert!(
        searcher
            .search(&query, false, &no_flank, control, &sender)
            .is_err()
    );
}

/// A caller that drops the receiver ends the search as if it asked to stop.
#[test]
fn dropped_receiver_stops() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = small_shard(directory.path());
    let searcher = Searcher::new(&shard, miss_model(), one_mismatch_model()).unwrap();
    let control = no_stop();
    let (sender, receiver) = sync_channel(1);
    drop(receiver);
    let summary = searcher
        .search(&query, false, &config(SearchStage::Align), control, &sender)
        .unwrap();
    assert!(matches!(
        summary.completion,
        Completion::Stopped(Stopped::Asked)
    ));
}

/// Every entry holds the query at exactly the accepted identity, linear or rotated; none may be
/// dropped by stages 1 and 2, while unrelated entries are.
#[test]
fn accepted_identity_kept() {
    let mut random = rng(0x1055);
    for (trial, circular) in [(0, false), (1, true)] {
        let directory = tempfile::tempdir().unwrap();
        let query = random_bases(&mut random, QUERY_LENGTH);
        let mut entries = Vec::new();
        for _ in 0..150 {
            let copy = diverged(&mut random, &query, 0.95);
            let copy = if circular {
                let cut = random.random_range(0..QUERY_LENGTH);
                [&copy[cut..], &copy[..cut]].concat()
            } else {
                copy
            };
            entries.push(vec![
                [
                    random_bases(&mut random, 500),
                    copy,
                    random_bases(&mut random, 500),
                ]
                .concat(),
            ]);
        }
        for _ in 0..20 {
            entries.push(vec![random_bases(&mut random, 6_000)]);
        }
        let shard = test_shard(directory.path(), &entries);
        for frequent_seed_entries in [None, Some(u64::MAX)] {
            let config = SearchConfig {
                frequent_seed_entries,
                ..config(SearchStage::Locate)
            };
            let (summary, events) = search(&shard, &query, circular, &config, false);
            assert!(matches!(summary.completion, Completion::Complete));
            for entry in 0..170u32 {
                let last = events.get(&entry).and_then(|list| list.last());
                if entry < 150 {
                    assert!(
                        matches!(last, Some(EntryEvent::Located(_))),
                        "trial {trial}, entry {entry}: {last:?}"
                    );
                } else {
                    assert!(
                        last.is_none()
                            || matches!(
                                last,
                                Some(EntryEvent::Dropped {
                                    reason: DropReason::CoverageBound,
                                    ..
                                })
                            ),
                        "trial {trial}, entry {entry}: {last:?}"
                    );
                }
            }
        }
    }
}

/// Each entry holds one piece of a long query at exactly the accepted identity, with long
/// stretches without seed hits around it; none may be dropped by stages 1 and 2. Entries with a
/// few short scattered pieces of the query must be dropped.
#[test]
fn long_stretch_kept() {
    const PIECE: usize = 4_000;
    let mut random = rng(0x9ece);
    let directory = tempfile::tempdir().unwrap();
    let query = random_bases(&mut random, 20_000);
    let mut entries = Vec::new();
    for _ in 0..150 {
        let start = random.random_range(0..query.len() - PIECE);
        let piece = diverged(&mut random, &query[start..start + PIECE], 0.95);
        entries.push(vec![
            [
                random_bases(&mut random, 500),
                piece,
                random_bases(&mut random, 500),
            ]
            .concat(),
        ]);
    }
    for _ in 0..20 {
        let mut sequence = random_bases(&mut random, 500);
        for _ in 0..3 {
            let start = random.random_range(0..query.len() - 50);
            sequence.extend_from_slice(&query[start..start + 50]);
            sequence.extend(random_bases(&mut random, 500));
        }
        entries.push(vec![sequence]);
    }
    let shard = test_shard(directory.path(), &entries);
    for frequent_seed_entries in [None, Some(u64::MAX)] {
        let mut config = SearchConfig {
            frequent_seed_entries,
            ..config(SearchStage::Locate)
        };
        config.rule.min_query_coverage = QueryCoverage::Bases(PIECE as u64);
        let (summary, events) = search(&shard, &query, false, &config, false);
        assert!(matches!(summary.completion, Completion::Complete));
        for entry in 0..170u32 {
            let last = events.get(&entry).and_then(|list| list.last());
            if entry < 150 {
                assert!(
                    matches!(last, Some(EntryEvent::Located(_))),
                    "entry {entry}: {last:?}"
                );
            } else {
                assert!(
                    last.is_none() || matches!(last, Some(EntryEvent::Dropped { .. })),
                    "entry {entry}: {last:?}"
                );
            }
        }
        // The scattered pieces have seed hits, so the bound, not their absence, drops them.
        assert!(summary.sketch.seeded_entries >= 160);
    }
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
    let EntryEvent::Confirmed { aligned, .. } = event else {
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
    let mut random = rng(0xa119);
    let query = random_bases(&mut random, QUERY_LENGTH);
    let mut changed = query.clone();
    for at in [1_000, 2_000, 3_000] {
        changed[at] = if changed[at] == b'A' { b'C' } else { b'A' };
    }
    let entries = vec![
        vec![
            random_bases(&mut random, 3_000),
            [
                random_bases(&mut random, 10_000),
                changed,
                random_bases(&mut random, 10_000),
            ]
            .concat(),
        ],
        vec![
            [
                random_bases(&mut random, 8_000),
                reverse_complement(&query),
                random_bases(&mut random, 4_000),
            ]
            .concat(),
        ],
        vec![random_bases(&mut random, 20_000)],
    ];
    (test_shard(directory, &entries), query)
}

#[test]
fn planted_query_confirmed() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = planted_shard(directory.path());
    let (summary, events) = search(&shard, &query, false, &config(SearchStage::Align), false);
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
    let EntryEvent::Confirmed { aligned, .. } = last_event(&events, 0) else {
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
fn decision_drops_entry() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = planted_shard(directory.path());
    // The query covers a fifth of the contig of entry 0 and a third of that of entry 1.
    let mut config = config(SearchStage::Align);
    config.rule.min_contig_coverage = 0.5;
    let (summary, events) = search(&shard, &query, false, &config, false);
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

/// An entry whose fragments pass the memory limit is dropped; the search still completes.
#[test]
fn oversized_entry_dropped() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = planted_shard(directory.path());
    let mut config = config(SearchStage::Align);
    config.align.max_fragment_bytes = 1;
    let (summary, events) = search(&shard, &query, false, &config, false);
    assert!(matches!(summary.completion, Completion::Complete));
    for entry in [0, 1] {
        assert_eq!(
            *last_event(&events, entry),
            EntryEvent::Dropped {
                entry,
                stage: SearchStage::Align,
                reason: DropReason::MemoryLimit,
            }
        );
    }
}

/// A query made of a repeat cannot reach the coverage with its other bases, so it has no
/// answer; with a unique part long enough it is searched as usual.
#[test]
fn repeat_query_unanswered() {
    let directory = tempfile::tempdir().unwrap();
    let mut random = rng(0x1c);
    let mut repeat: Vec<u8> = b"ACG".iter().copied().cycle().take(6_000).collect();
    for at in (40..6_000).step_by(70) {
        repeat[at] = b'T';
    }
    let unique = random_bases(&mut random, 5_000);
    let entries = vec![vec![[unique.clone(), repeat.clone()].concat()]];
    let shard = test_shard(directory.path(), &entries);
    let config = config(SearchStage::Align);
    let (summary, events) = search(&shard, &repeat[1_000..5_000], false, &config, false);
    assert!(matches!(summary.completion, Completion::Complete));
    assert_eq!(summary.sketch.low_complexity_bases, 4_000);
    assert!(events.is_empty(), "{events:?}");
    let mixed = [&unique[1_000..], &repeat[..300]].concat();
    let (summary, events) = search(&shard, &mixed, false, &config, false);
    // A DUST window that reaches over the border is masked whole.
    assert!((300..364).contains(&summary.sketch.low_complexity_bases));
    assert!(matches!(
        last_event(&events, 0),
        EntryEvent::Confirmed { .. }
    ));
}

#[test]
fn circular_query_confirmed() {
    let directory = tempfile::tempdir().unwrap();
    let mut random = rng(0xc1c);
    let query = random_bases(&mut random, QUERY_LENGTH);
    let rotated = [&query[1_700..], &query[..1_700]].concat();
    let entries = vec![
        vec![
            [
                random_bases(&mut random, 6_000),
                rotated,
                random_bases(&mut random, 6_000),
            ]
            .concat(),
        ],
        vec![random_bases(&mut random, 15_000)],
    ];
    let shard = test_shard(directory.path(), &entries);
    let (summary, events) = search(&shard, &query, true, &config(SearchStage::Align), false);
    assert_eq!(summary.align.confirmed_entries, 1);
    let EntryEvent::Confirmed { aligned, .. } = last_event(&events, 0) else {
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

/// A packed sequence that disappears after the shard is open fails stage 3 for its entry.
#[test]
fn failure_marks_unfinished() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = planted_shard(directory.path());
    for entry in 0..2 {
        std::fs::remove_file(shard.packed_path(entry).unwrap()).unwrap();
    }
    let (summary, events) = search(&shard, &query, false, &config(SearchStage::Align), false);
    assert!(
        matches!(summary.completion, Completion::Failed(_)),
        "{:?}",
        summary.completion
    );
    for entry in [0, 1] {
        let last = events[&entry].last();
        assert!(
            matches!(last, Some(EntryEvent::Unfinished { .. })),
            "entry {entry}: {last:?}"
        );
    }
}
