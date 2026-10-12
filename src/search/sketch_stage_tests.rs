//! Tests of stage 1: seed hits, frequent seeds, the memory limit and the candidate order.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::bound::CoverageBound;
use super::control::{SearchControl, Stopped};
use super::error::SearchError;
use super::query::QueryKmers;
use super::sketch_stage::{SketchLimits, keep_found, sketch_stage};
use crate::index::{Shard, build_shard};
use crate::loss::MaxMiss;
use crate::seeds::Seeder;
use crate::test_support::{
    build_options, miss_model, no_stop, parameters, random_bases, rng, test_shard,
    with_substitutions, write_entries,
};
use std::sync::atomic::AtomicBool;

const UNLIMITED: SketchLimits = SketchLimits {
    frequent_seed_entries: u64::MAX,
    max_seed_hits: u64::MAX,
};

/// Entry 0 holds the query, entry 1 a copy at 97% identity, entry 2 the first 30% of it, and
/// entry 3 is unrelated.
fn shard(directory: &std::path::Path) -> (Shard, Vec<u8>) {
    let mut random = rng(0x5eed);
    let query = random_bases(&mut random, 5_000);
    let close = with_substitutions(&mut random, &query, 0.97);
    let entries = vec![
        vec![
            [
                random_bases(&mut random, 3_000),
                query.clone(),
                random_bases(&mut random, 2_000),
            ]
            .concat(),
        ],
        vec![
            random_bases(&mut random, 1_000),
            [close, random_bases(&mut random, 500)].concat(),
        ],
        vec![[random_bases(&mut random, 4_000), query[..1_500].to_vec()].concat()],
        vec![random_bases(&mut random, 8_000)],
    ];
    (test_shard(directory, &entries), query)
}

fn bound(length: u64, required_bases: u64) -> CoverageBound<'static> {
    bound_at(length, required_bases, 0.95)
}

fn bound_at(length: u64, required_bases: u64, identity: f64) -> CoverageBound<'static> {
    CoverageBound::new(
        miss_model(),
        MaxMiss::default(),
        identity,
        length,
        false,
        required_bases,
    )
}

fn run(shard: &Shard, query: &[u8], required: u64, limits: SketchLimits) -> Vec<(u32, u64)> {
    run_with(shard, query, &bound(query.len() as u64, required), limits).0
}

/// Runs stage 1 and returns the candidate entries with their bounds, and the entries left out
/// for lack of any seed hit.
fn run_with(
    shard: &Shard,
    query: &[u8],
    bound: &CoverageBound,
    limits: SketchLimits,
) -> (Vec<(u32, u64)>, u64) {
    let seeder = Seeder::new(parameters().scheme).unwrap();
    let kmers = QueryKmers::new(&seeder, query, false);
    let control = no_stop();
    let outcome = sketch_stage(shard.seed_file(), &kmers, bound, limits, control).unwrap();
    assert_eq!(
        outcome.work.possible_entries,
        outcome.candidates.len() as u64
    );
    let kept = outcome
        .candidates
        .iter()
        .map(|candidate| (candidate.entry, candidate.coverage_bound))
        .collect();
    (kept, outcome.work.unseeded_entries)
}

#[test]
fn low_bounds_dropped() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = shard(directory.path());
    let kept = run(&shard, &query, 3_500, UNLIMITED);
    assert_eq!(kept.len(), 2, "{kept:?}");
    assert_eq!(kept[0], (0, 5_000));
    assert_eq!(kept[1].0, 1);
    // With a low limit the partial entry stays, ordered after the full ones.
    let kept = run(&shard, &query, 1_000, UNLIMITED);
    assert_eq!(
        kept.iter().map(|(entry, _)| *entry).collect::<Vec<_>>(),
        [0, 1, 2]
    );
    assert!(kept[2].1 >= 1_500 && kept[2].1 < 3_500, "{kept:?}");
}

/// The containment estimate of each entry of the 4 entry shard; 0 when it is not a candidate.
fn estimates(shard: &Shard, query: &[u8], limits: SketchLimits) -> Vec<f64> {
    let kmers = QueryKmers::new(&Seeder::new(parameters().scheme).unwrap(), query, false);
    let bound = bound(query.len() as u64, 0);
    let outcome = sketch_stage(shard.seed_file(), &kmers, &bound, limits, no_stop()).unwrap();
    let mut estimates = vec![0.0; 4];
    for candidate in &outcome.candidates {
        estimates[candidate.entry as usize] = candidate.containment_estimate;
    }
    estimates
}

#[test]
fn containment_follows_share() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = shard(directory.path());
    let estimates = estimates(&shard, &query, UNLIMITED);
    assert!(estimates[0] > 0.8, "{estimates:?}");
    assert!(estimates[2] > 0.15 && estimates[2] < 0.45, "{estimates:?}");
    // The unrelated entry has no seed hit, so it is not a candidate at all.
    assert_eq!(estimates[3], 0.0);
}

#[test]
fn frequent_seeds_count() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = shard(directory.path());
    // Seeds held by more than one entry are frequent; they keep the entries of their color sets
    // possible over the stretches they cover.
    let limits = SketchLimits {
        frequent_seed_entries: 1,
        max_seed_hits: u64::MAX,
    };
    let kept = run(&shard, &query, 2_500, limits);
    assert!(kept.iter().any(|(entry, _)| *entry == 0));
    // With no room for seed hits every seed is frequent. The entries in their color sets keep
    // the full bound; the unrelated entry is in none of them and is left out.
    let limits = SketchLimits {
        frequent_seed_entries: u64::MAX,
        max_seed_hits: 0,
    };
    let kept = run(&shard, &query, 2_500, limits);
    assert_eq!(
        kept.iter().map(|(entry, _)| *entry).collect::<Vec<_>>(),
        [0, 1, 2]
    );
    assert!(kept.iter().all(|(_, bound)| *bound == 5_000));
    // Every seed is frequent, and each entry still counts only the seeds of its color sets.
    let estimates = estimates(&shard, &query, limits);
    assert!(estimates[0] > 0.8, "{estimates:?}");
    assert!(estimates[2] > 0.15 && estimates[2] < 0.45, "{estimates:?}");
}

#[test]
fn stop_ends_sketch() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = shard(directory.path());
    let seeder = Seeder::new(parameters().scheme).unwrap();
    let kmers = QueryKmers::new(&seeder, &query, false);
    let stop = AtomicBool::new(true);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let outcome = sketch_stage(
        shard.seed_file(),
        &kmers,
        &bound(5_000, 0),
        UNLIMITED,
        control,
    );
    assert!(matches!(outcome, Err(SearchError::Stopped(Stopped::Asked))));
}

#[test]
fn seedless_entries_dropped() {
    let directory = tempfile::tempdir().unwrap();
    let mut random = rng(0x10e);
    let query = random_bases(&mut random, 5_000);
    let entries = vec![
        vec![[random_bases(&mut random, 1_000), query.clone()].concat()],
        // One short piece of the query: a single seed hit or very few.
        vec![
            [
                random_bases(&mut random, 3_000),
                query[2_000..2_060].to_vec(),
                random_bases(&mut random, 3_000),
            ]
            .concat(),
        ],
        vec![random_bases(&mut random, 6_000)],
    ];
    let shard = test_shard(directory.path(), &entries);
    // At 80% a stretch without seed hits may hold the query, so every bound is the whole query.
    let low = bound_at(5_000, 4_500, 0.8);
    let (kept, without) = run_with(&shard, &query, &low, UNLIMITED);
    assert_eq!(kept, [(0, 5_000), (1, 5_000)]);
    assert_eq!(without, 1);
}

#[test]
fn frequent_union_only() {
    let directory = tempfile::tempdir().unwrap();
    let mut random = rng(0xf2e);
    let (first, second) = (
        random_bases(&mut random, 3_000),
        random_bases(&mut random, 3_000),
    );
    let query = [first.clone(), second.clone()].concat();
    // The first half is in entries 0, 1, 2 and 4, so its seeds are frequent at a limit of 2. The
    // second half is in entries 3 and 4 only, so its seeds are expanded.
    let mut entries: Vec<Vec<Vec<u8>>> = (0..3)
        .map(|_| vec![[random_bases(&mut random, 500), first.clone()].concat()])
        .collect();
    entries.push(vec![
        [random_bases(&mut random, 500), second.clone()].concat(),
    ]);
    entries.push(vec![query.clone()]);
    // Many buckets, so no seed of entry 3 shares a stored key with a seed of the first half.
    let mut options = build_options(directory.path(), 2);
    options.bucket_bits = Some(20);
    let output = directory.path().join("shard.jam");
    build_shard(
        &write_entries(directory.path(), &entries),
        &output,
        &options,
    )
    .unwrap();
    let shard = Shard::open(&output).unwrap();
    let limits = SketchLimits {
        frequent_seed_entries: 2,
        max_seed_hits: u64::MAX,
    };
    let (kept, without) = run_with(&shard, &query, &bound(6_000, 3_000), limits);
    let bound_of = |entry: u32| kept.iter().find(|(e, _)| *e == entry).unwrap().1;
    assert_eq!(without, 0);
    assert_eq!(bound_of(4), 6_000);
    for entry in 0..3 {
        assert!(bound_of(entry) >= 3_000 && bound_of(entry) < 5_000);
    }
    // Entry 3 is in no frequent color set: only its own seed hits on the second half count.
    assert!(
        bound_of(3) >= 3_000 && bound_of(3) < 5_000,
        "{}",
        bound_of(3)
    );
}

/// Keeping only the lookups the index holds leaves the same candidates as looking up all of
/// them, with the one-substitution neighbours as well.
#[test]
fn found_lookups_kept() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = shard(directory.path());
    let seeder = Seeder::new(parameters().scheme).unwrap();
    let all = QueryKmers::with_neighbours(&seeder, &query, false);
    let plan = super::seed_mismatches::LookupPlan::full(1);
    let keep = |lookups: &mut Vec<(u64, u32)>| keep_found(shard.seed_file(), lookups);
    let kept = QueryKmers::filtered(&seeder, &query, false, &plan, keep).unwrap();
    assert!(kept.distinct().len() < all.distinct().len() / 10);
    let bound = bound(query.len() as u64, 1_000);
    let entries = |kmers: &QueryKmers| {
        let outcome = sketch_stage(shard.seed_file(), kmers, &bound, UNLIMITED, no_stop());
        let outcome = outcome.unwrap();
        let entries: Vec<(u32, u64)> = (outcome.candidates.iter())
            .map(|possible| (possible.entry, possible.coverage_bound))
            .collect();
        (entries, outcome.work.seed_hits, outcome.work.lookups)
    };
    let (kept_entries, kept_hits, lookups) = entries(&kept);
    let (all_entries, all_hits, _) = entries(&all);
    assert_eq!((kept_entries, kept_hits), (all_entries, all_hits));
    assert_eq!(lookups, all.lookups());
}
