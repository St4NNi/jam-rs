//! Tests of stage 2: false seed hits, the second pruning and the alignment windows.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::bound::CoverageBound;
use super::driver::ENTRIES_AHEAD;
use super::events::DropReason;
use super::locate_stage::{
    LocateOutcome, LocateRequest, LocateWork, LocateWorkspace, ahead_bytes, locate_entry,
    read_entry_locations,
};
use super::query::QueryKmers;
use super::sketch_stage::{PossibleEntry, SketchLimits, find_query_seeds, sketch_stage};
use crate::index::{Shard, build_shard};
use crate::loss::MaxMiss;
use crate::seeds::Seeder;
use crate::test_support::{
    build_options, miss_model, no_stop, parameters, random_bases, rng, test_shard, write_entries,
};
use std::ops::Range;

const FLANK: u64 = 1_000;

fn kmers(query: &[u8]) -> QueryKmers {
    QueryKmers::new(&Seeder::new(parameters().scheme).unwrap(), query, false)
}

fn bound(length: u64, required_bases: u64) -> CoverageBound<'static> {
    CoverageBound::new(
        miss_model(),
        MaxMiss::default(),
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
    let control = no_stop();
    let request = LocateRequest {
        query,
        bound,
        frequent_seeds: &[],
        flank: FLANK,
        verify: None,
    };
    let mut work = LocateWork::default();
    let found = find_query_seeds(shard.seed_file(), query, control).unwrap();
    let candidates = std::slice::from_ref(possible);
    let file = shard.seed_file();
    let locations = read_entry_locations(file, &request, &found, candidates, control, &mut work);
    let outcome = locate_entry(
        shard,
        &request,
        possible,
        &locations.unwrap()[0],
        &mut LocateWorkspace::default(),
        &mut work,
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
    let control = no_stop();
    let limits = SketchLimits {
        frequent_seed_entries: u64::MAX,
        max_seed_hits: u64::MAX,
    };
    let outcome = sketch_stage(shard.seed_file(), &query_kmers, &bound, limits, control).unwrap();
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
            .alignment_windows
            .iter()
            .map(|contig| (contig.contig, contig.ranges.clone()))
            .collect(),
        LocateOutcome::Dropped(reason) => panic!("dropped: {reason:?}"),
    }
}

#[test]
fn false_hits_dropped() {
    let directory = tempfile::tempdir().unwrap();
    let mut random = rng(0xfa15e5);
    let query = random_bases(&mut random, 20_000);
    // The fewest buckets put many seeds of the large unrelated entry into each bucket, so a short
    // fingerprint would match some query k-mers; the full check value matches none.
    let entries = vec![
        vec![random_bases(&mut random, 400_000)],
        vec![query.clone()],
    ];
    let mut options = build_options(directory.path(), 2);
    options.bucket_bits = Some(8);
    let output = directory.path().join("shard.jam");
    build_shard(
        &write_entries(directory.path(), &entries),
        &output,
        &options,
    )
    .unwrap();
    let shard = Shard::open(&output).unwrap();
    let results = search(&shard, &query, 0);
    assert!(results.iter().all(|(possible, ..)| possible.entry != 0));
    let (_, outcome, _) = results.iter().find(|(p, ..)| p.entry == 1).unwrap();
    assert_eq!(
        windows(outcome),
        [(0, std::iter::once(0..20_000).collect())]
    );
}

#[test]
fn second_pruning_drops() {
    let directory = tempfile::tempdir().unwrap();
    let mut random = rng(0xb0b);
    let query = random_bases(&mut random, 5_000);
    let entries = vec![vec![random_bases(&mut random, 10_000)]];
    let shard = test_shard(directory.path(), &entries);
    let query_kmers = kmers(&query);
    // As if the sketch had placed every query k-mer in the unrelated entry.
    let possible = PossibleEntry {
        entry: 0,
        coverage_bound: 5_000,
        containment_estimate: 1.0,
        seed_hits: (0..query_kmers.distinct().len() as u32).collect(),
    };
    let (outcome, work) = locate(&shard, &query_kmers, &bound(5_000, 2_500), &possible);
    assert_eq!(outcome, LocateOutcome::Dropped(DropReason::CoverageBound));
    assert_eq!(work.seed_locations, 0);
    assert_eq!(work.dropped_entries, 1);
}

#[test]
fn flank_widens_bins() {
    let directory = tempfile::tempdir().unwrap();
    let mut random = rng(0x3d);
    let query = random_bases(&mut random, 4_000);
    let entries = vec![
        // The query in the middle of the second contig.
        vec![
            random_bases(&mut random, 3_000),
            [
                random_bases(&mut random, 20_000),
                query.clone(),
                random_bases(&mut random, 20_000),
            ]
            .concat(),
        ],
        // The query split over the end of one contig and the start of the next.
        vec![
            [random_bases(&mut random, 10_000), query[..2_000].to_vec()].concat(),
            [query[2_000..].to_vec(), random_bases(&mut random, 10_000)].concat(),
        ],
    ];
    let shard = test_shard(directory.path(), &entries);
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

/// Stage 2 of entry 0 with sequence checks, after `edit` changes its seed locations.
fn checked(
    shard: &Shard,
    query: &QueryKmers,
    edit: impl Fn(&mut Vec<(u32, crate::index::SeedLocation)>),
) -> (LocateOutcome, LocateWork) {
    let (control, seeder) = (no_stop(), Seeder::new(parameters().scheme).unwrap());
    let bound = bound(query.length(), 0);
    let request = LocateRequest {
        query,
        bound: &bound,
        frequent_seeds: &[],
        flank: FLANK,
        verify: Some(&seeder),
    };
    let found = find_query_seeds(shard.seed_file(), query, control).unwrap();
    let limits = SketchLimits {
        frequent_seed_entries: u64::MAX,
        max_seed_hits: u64::MAX,
    };
    let sketch = sketch_stage(shard.seed_file(), query, &bound, limits, control).unwrap();
    let possible = sketch.candidates.iter().find(|p| p.entry == 0).unwrap();
    let mut work = LocateWork::default();
    let candidates = std::slice::from_ref(possible);
    let file = shard.seed_file();
    let read = read_entry_locations(file, &request, &found, candidates, control, &mut work);
    let mut locations = read.unwrap().remove(0);
    edit(&mut locations);
    let mut workspace = LocateWorkspace::default();
    let outcome = locate_entry(
        shard,
        &request,
        possible,
        &locations,
        &mut workspace,
        &mut work,
    );
    (outcome.unwrap(), work)
}

/// A weak entry whose seed hits lie in its sequence is kept. Moved to a bin that does not hold
/// the k-mers, as a chance match of a seed key would be, it is dropped.
#[test]
fn chance_keys_dropped() {
    let directory = tempfile::tempdir().unwrap();
    let mut random = rng(0xc011);
    let query = random_bases(&mut random, 2_000);
    let piece = query[700..760].to_vec();
    let entries = vec![vec![
        [
            random_bases(&mut random, 5_000),
            piece,
            random_bases(&mut random, 15_000),
        ]
        .concat(),
    ]];
    let shard = test_shard(directory.path(), &entries);
    let query = kmers(&query);
    let (kept, work) = checked(&shard, &query, |_| {});
    assert!(matches!(kept, LocateOutcome::Located(_)), "{kept:?}");
    assert_eq!((work.unverified_entries, work.verified_bins), (0, 1));
    // Bin 60 of 256 bases starts at 15,360, far from the piece.
    let moved = |locations: &mut Vec<(u32, crate::index::SeedLocation)>| {
        locations
            .iter_mut()
            .for_each(|(_, location)| location.bin = 60);
    };
    let (dropped, work) = checked(&shard, &query, moved);
    assert_eq!(dropped, LocateOutcome::Dropped(DropReason::NoSeedHits));
    assert_eq!((work.unverified_entries, work.verified_bins), (1, 1));
}

/// The entries ahead of the workers ask for as many bytes in all as 128 entries of 1/1024 of the
/// block cache, at least 1 MiB each, did: an eighth of the cache, at least 128 MiB.
#[test]
fn entries_ahead_bytes() {
    for ram in [64 << 20, 1 << 30, 4 << 30, 16 << 30, 100 << 30] {
        let before = 128 * (ram / 1024).max(1 << 20);
        assert_eq!(
            ENTRIES_AHEAD as u64 * ahead_bytes(Some(ram)),
            before,
            "{ram}"
        );
    }
    assert_eq!(ahead_bytes(None), u64::MAX);
}
