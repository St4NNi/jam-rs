use super::bound::CoverageBound;
use super::control::SearchControl;
use super::events::DropReason;
use super::locate_stage::{
    LocateOutcome, LocateRequest, LocateWork, LocateWorkspace, locate_entry,
};
use super::query::QueryKmers;
use super::sketch_stage::{PossibleEntry, SketchLimits, sketch_stage};
use super::test_shard::{Random, build, model, parameters};
use crate::index::Shard;
use crate::loss::LossBudget;
use crate::seeds::Seeder;
use std::ops::Range;
use std::sync::atomic::AtomicBool;

const FLANK: u64 = 1_000;

fn kmers(query: &[u8]) -> QueryKmers {
    QueryKmers::new(&Seeder::new(parameters().scheme()).unwrap(), query, false)
}

fn bound(length: u64, required_bases: u64) -> CoverageBound<'static> {
    CoverageBound::new(
        model(),
        LossBudget::default(),
        0.95,
        length,
        false,
        required_bases,
    )
}

/// Runs stage 2 for `possible` and returns the outcome and the work.
fn locate(
    shard: &Shard,
    query: &QueryKmers,
    bound: &CoverageBound,
    possible: &PossibleEntry,
) -> (LocateOutcome, LocateWork) {
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let request = LocateRequest {
        query,
        bound,
        frequent_seeds: &[],
        flank: FLANK,
    };
    let mut work = LocateWork::default();
    let outcome = locate_entry(
        shard,
        &request,
        possible,
        &mut LocateWorkspace::default(),
        &mut work,
        control,
    )
    .unwrap();
    (outcome, work)
}

/// Stages 1 and 2 for every candidate entry.
fn search(
    shard: &Shard,
    query: &[u8],
    required: u64,
) -> Vec<(PossibleEntry, LocateOutcome, LocateWork)> {
    let query_kmers = kmers(query);
    let bound = bound(query.len() as u64, required);
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let limits = SketchLimits {
        frequent_seed_entries: u64::MAX,
        max_seed_hits: u64::MAX,
    };
    let outcome = sketch_stage(shard.sketch(), &query_kmers, &bound, limits, control).unwrap();
    outcome
        .candidates
        .into_iter()
        .map(|possible| {
            let (located, work) = locate(shard, &query_kmers, &bound, &possible);
            (possible, located, work)
        })
        .collect()
}

fn windows(outcome: &LocateOutcome) -> Vec<(u32, Vec<Range<u64>>)> {
    match outcome {
        LocateOutcome::Located(located) => located
            .windows
            .iter()
            .map(|contig| (contig.contig, contig.ranges.clone()))
            .collect(),
        LocateOutcome::Dropped(reason) => panic!("dropped: {reason:?}"),
    }
}

#[test]
fn false_seed_hits_of_the_sketch_have_no_seed_locations() {
    let directory = tempfile::tempdir().unwrap();
    let mut random = Random(0xfa15e);
    let query = random.bases(20_000);
    // The fewest buckets put many seeds of the large unrelated entry into each bucket, so some
    // share a fingerprint with a query k-mer.
    let entries = vec![vec![random.bases(400_000)], vec![query.clone()]];
    let shard = build(directory.path(), &entries, Some(8));
    let results = search(&shard, &query, 0);
    let (possible, outcome, work) = results.iter().find(|(p, ..)| p.entry == 0).unwrap();
    assert!(
        !possible.seed_hits.is_empty(),
        "no false seed hit was forced"
    );
    assert_eq!(work.false_seed_hits, possible.seed_hits.len() as u64);
    assert_eq!(work.seed_locations, 0);
    assert_eq!(*outcome, LocateOutcome::Dropped(DropReason::NoSeedHits));
    let (_, outcome, work) = results.iter().find(|(p, ..)| p.entry == 1).unwrap();
    assert_eq!(work.false_seed_hits, 0);
    assert_eq!(
        windows(outcome),
        [(0, std::iter::once(0..20_000).collect())]
    );
}

#[test]
fn the_bound_from_the_seed_locations_drops_entries_again() {
    let directory = tempfile::tempdir().unwrap();
    let mut random = Random(0xb0b);
    let query = random.bases(5_000);
    let entries = vec![vec![random.bases(10_000)]];
    let shard = build(directory.path(), &entries, None);
    let query_kmers = kmers(&query);
    // As if the sketch had placed every query k-mer in the unrelated entry.
    let possible = PossibleEntry {
        entry: 0,
        coverage_bound: 5_000,
        containment_estimate: 1.0,
        seed_hits: (0..query_kmers.seeds().len() as u32).collect(),
    };
    let (outcome, work) = locate(&shard, &query_kmers, &bound(5_000, 2_500), &possible);
    assert_eq!(outcome, LocateOutcome::Dropped(DropReason::CoverageBound));
    assert_eq!(work.false_seed_hits, possible.seed_hits.len() as u64);
    assert_eq!(work.dropped_entries, 1);
}

#[test]
fn windows_widen_the_bins_by_the_flank_and_stop_at_contig_ends() {
    let directory = tempfile::tempdir().unwrap();
    let mut random = Random(0x3d);
    let query = random.bases(4_000);
    let entries = vec![
        // The query in the middle of the second contig.
        vec![
            random.bases(3_000),
            [random.bases(20_000), query.clone(), random.bases(20_000)].concat(),
        ],
        // The query split over the end of one contig and the start of the next.
        vec![
            [random.bases(10_000), query[..2_000].to_vec()].concat(),
            [query[2_000..].to_vec(), random.bases(10_000)].concat(),
        ],
    ];
    let shard = build(directory.path(), &entries, None);
    let results = search(&shard, &query, 2_000);
    assert_eq!(results.len(), 2);
    for (possible, outcome, _) in &results {
        let windows = windows(outcome);
        if possible.entry == 0 {
            // Bins hold 256 bases and start at multiples of 256 in entry coordinates.
            assert_eq!(windows.len(), 1, "{windows:?}");
            let (contig, ranges) = &windows[0];
            assert_eq!((*contig, ranges.len()), (1, 1));
            let range = &ranges[0];
            assert!(range.start <= 20_000 - FLANK && range.start + 256 + FLANK >= 20_000);
            assert!(range.end >= 24_000 + FLANK - 21 && range.end <= 24_000 + FLANK + 256);
        } else {
            assert_eq!(windows.len(), 2, "{windows:?}");
            assert_eq!(windows[0].0, 0);
            assert_eq!(windows[0].1.len(), 1);
            assert!(windows[0].1[0].start <= 10_000 - FLANK);
            assert_eq!(windows[0].1[0].end, 12_000);
            assert_eq!(windows[1].0, 1);
            assert_eq!(windows[1].1[0].start, 0);
            assert!(windows[1].1[0].end >= 2_000 + FLANK - 21);
        }
    }
}
