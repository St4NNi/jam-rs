use super::control::{SearchControl, Stopped};
use super::driver::{SearchConfig, SearchSummary, Searcher};
use super::events::{Completion, DropReason, EntryEvent, SearchStage};
use super::test_shard::{Random, build, model};
use crate::decide::QueryCoverage;
use crate::index::Shard;
use std::collections::BTreeMap;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, Ordering};

const QUERY_LENGTH: usize = 5_000;

/// A copy of `sequence` with exactly `1 - identity` of its bases changed, at random places.
fn diverged(random: &mut Random, sequence: &[u8], identity: f64) -> Vec<u8> {
    let mut copy = sequence.to_vec();
    let changes = ((1.0 - identity) * sequence.len() as f64).round() as usize;
    let mut changed = vec![false; sequence.len()];
    let mut done = 0;
    while done < changes {
        let at = (random.next() % sequence.len() as u64) as usize;
        if !changed[at] {
            changed[at] = true;
            let others: Vec<u8> = b"ACGT".iter().copied().filter(|&b| b != copy[at]).collect();
            copy[at] = others[(random.next() % 3) as usize];
            done += 1;
        }
    }
    copy
}

fn config(level: SearchStage) -> SearchConfig {
    let mut config = SearchConfig {
        level,
        threads: 3,
        flank: 1_000,
        ..SearchConfig::default()
    };
    config.rule.min_identity = 0.95;
    config.rule.min_query_coverage = QueryCoverage::Share {
        numerator: 9,
        decimal_places: 1,
    };
    config
}

/// Runs a search and returns the summary and the events of each entry, in order.
fn search(
    shard: &Shard,
    query: &[u8],
    circular: bool,
    config: &SearchConfig,
    stop_at_first_event: bool,
) -> (SearchSummary, BTreeMap<u32, Vec<EntryEvent>>) {
    let searcher = Searcher::new(shard, model()).unwrap();
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let events = Mutex::new(BTreeMap::<u32, Vec<EntryEvent>>::new());
    let on_event = |event: EntryEvent| {
        if stop_at_first_event {
            stop.store(true, Ordering::Relaxed);
        }
        let mut events = events.lock().unwrap();
        events.entry(event.entry()).or_default().push(event);
    };
    let summary = searcher
        .search(query, circular, config, control, &on_event)
        .unwrap();
    (summary, events.into_inner().unwrap())
}

/// Entries 0 to 3 hold the query at 97% identity, entries 4 to 7 are unrelated.
fn small_shard(directory: &std::path::Path) -> (Shard, Vec<u8>) {
    let mut random = Random(0xd1e);
    let query = random.bases(QUERY_LENGTH);
    let mut entries = Vec::new();
    for _ in 0..4 {
        let copy = diverged(&mut random, &query, 0.97);
        entries.push(vec![
            [random.bases(2_000), copy, random.bases(2_000)].concat(),
        ]);
    }
    for _ in 0..4 {
        entries.push(vec![random.bases(9_000)]);
    }
    (build(directory, &entries, None), query)
}

#[test]
fn the_search_stops_after_the_requested_stage() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = small_shard(directory.path());
    let (summary, events) = search(&shard, &query, false, &config(SearchStage::Sketch), false);
    assert_eq!(summary.completion, Completion::Complete);
    assert_eq!(summary.locate.entries, 0);
    assert_eq!(events.keys().copied().collect::<Vec<_>>(), [0, 1, 2, 3]);
    assert!(
        events
            .values()
            .all(|list| list.len() == 1 && matches!(list[0], EntryEvent::Possible { .. }))
    );
    let (summary, events) = search(&shard, &query, false, &config(SearchStage::Locate), false);
    assert_eq!(summary.completion, Completion::Complete);
    assert_eq!(summary.locate.entries, summary.sketch.possible_entries);
    for (entry, list) in &events {
        assert_eq!(list.len(), 2, "{list:?}");
        assert!(matches!(list[0], EntryEvent::Possible { .. }));
        let EntryEvent::Located(located) = &list[1] else {
            panic!("entry {entry}: {list:?}");
        };
        assert_eq!(located.windows.len(), 1);
        assert!(located.coverage_bound >= 4_500);
    }
}

#[test]
fn a_stopped_search_reports_the_rest_as_unfinished() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = small_shard(directory.path());
    let (summary, events) = search(&shard, &query, false, &config(SearchStage::Locate), true);
    assert_eq!(summary.completion, Completion::Stopped(Stopped::Asked));
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
fn bad_queries_and_settings_fail_at_once() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, _) = small_shard(directory.path());
    let searcher = Searcher::new(&shard, model()).unwrap();
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let config = config(SearchStage::Locate);
    assert!(
        searcher
            .search(b"ACGT", false, &config, control, &|_| {})
            .is_err()
    );
    let no_flank = SearchConfig { flank: 0, ..config };
    let query = vec![b'A'; 100];
    assert!(
        searcher
            .search(&query, false, &no_flank, control, &|_| {})
            .is_err()
    );
}

/// Every entry holds the query at exactly the accepted identity, linear or rotated; none may be
/// dropped by stages 1 and 2, while unrelated entries are.
#[test]
fn entries_at_the_accepted_identity_are_never_dropped() {
    let mut random = Random(0x1055);
    for (trial, circular) in [(0, false), (1, true)] {
        let directory = tempfile::tempdir().unwrap();
        let query = random.bases(QUERY_LENGTH);
        let mut entries = Vec::new();
        for _ in 0..150 {
            let copy = diverged(&mut random, &query, 0.95);
            let copy = if circular {
                let cut = (random.next() % QUERY_LENGTH as u64) as usize;
                [&copy[cut..], &copy[..cut]].concat()
            } else {
                copy
            };
            entries.push(vec![[random.bases(500), copy, random.bases(500)].concat()]);
        }
        for _ in 0..20 {
            entries.push(vec![random.bases(6_000)]);
        }
        let shard = build(directory.path(), &entries, None);
        for frequent_seed_entries in [None, Some(u64::MAX)] {
            let config = SearchConfig {
                frequent_seed_entries,
                ..config(SearchStage::Locate)
            };
            let (summary, events) = search(&shard, &query, circular, &config, false);
            assert_eq!(summary.completion, Completion::Complete);
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
fn a_present_piece_inside_a_long_stretch_is_never_dropped() {
    const PIECE: usize = 4_000;
    let mut random = Random(0x9ece);
    let directory = tempfile::tempdir().unwrap();
    let query = random.bases(20_000);
    let mut entries = Vec::new();
    for _ in 0..150 {
        let start = (random.next() % (query.len() - PIECE) as u64) as usize;
        let piece = diverged(&mut random, &query[start..start + PIECE], 0.95);
        entries.push(vec![[random.bases(500), piece, random.bases(500)].concat()]);
    }
    for _ in 0..20 {
        let mut sequence = random.bases(500);
        for _ in 0..3 {
            let start = (random.next() % (query.len() - 50) as u64) as usize;
            sequence.extend_from_slice(&query[start..start + 50]);
            sequence.extend(random.bases(500));
        }
        entries.push(vec![sequence]);
    }
    let shard = build(directory.path(), &entries, None);
    for frequent_seed_entries in [None, Some(u64::MAX)] {
        let mut config = SearchConfig {
            frequent_seed_entries,
            ..config(SearchStage::Locate)
        };
        config.rule.min_query_coverage = QueryCoverage::Bases(PIECE as u64);
        let (summary, events) = search(&shard, &query, false, &config, false);
        assert_eq!(summary.completion, Completion::Complete);
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
        assert!(summary.sketch.entries_with_seed_hits >= 160);
    }
}
