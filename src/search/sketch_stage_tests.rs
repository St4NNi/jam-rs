use super::bound::CoverageBound;
use super::control::{SearchControl, Stopped};
use super::query::QueryKmers;
use super::sketch_stage::{SketchLimits, sketch_stage};
use super::test_shard::{Random, build, model, parameters};
use crate::index::Shard;
use crate::loss::LossBudget;
use crate::seeds::Seeder;
use std::sync::atomic::AtomicBool;

const UNLIMITED: SketchLimits = SketchLimits {
    frequent_seed_entries: u64::MAX,
    max_seed_hits: u64::MAX,
};

/// Entry 0 holds the query, entry 1 a copy at 97% identity, entry 2 the first 30% of it, and
/// entry 3 is unrelated.
fn shard(directory: &std::path::Path) -> (Shard, Vec<u8>) {
    let mut random = Random(0x5eed);
    let query = random.bases(5_000);
    let close = random.mutate(&query, 0.97);
    let entries = vec![
        vec![[random.bases(3_000), query.clone(), random.bases(2_000)].concat()],
        vec![random.bases(1_000), [close, random.bases(500)].concat()],
        vec![[random.bases(4_000), query[..1_500].to_vec()].concat()],
        vec![random.bases(8_000)],
    ];
    (build(directory, &entries, None), query)
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

fn run(shard: &Shard, query: &[u8], required: u64, limits: SketchLimits) -> Vec<(u32, u64)> {
    let seeder = Seeder::new(parameters().scheme()).unwrap();
    let kmers = QueryKmers::new(&seeder, query, false);
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let bound = bound(query.len() as u64, required);
    let outcome = sketch_stage(shard.sketch(), &kmers, &bound, limits, control).unwrap();
    assert_eq!(
        outcome.work.possible_entries,
        outcome.candidates.len() as u64
    );
    outcome
        .candidates
        .iter()
        .map(|candidate| (candidate.entry, candidate.coverage_bound))
        .collect()
}

#[test]
fn entries_below_the_coverage_limit_are_left_out() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = shard(directory.path());
    let kept = run(&shard, &query, 2_500, UNLIMITED);
    assert_eq!(kept.len(), 2, "{kept:?}");
    assert_eq!(kept[0], (0, 5_000));
    assert_eq!(kept[1].0, 1);
    // With a low limit the partial entry stays, ordered after the full ones.
    let kept = run(&shard, &query, 1_000, UNLIMITED);
    assert_eq!(
        kept.iter().map(|(entry, _)| *entry).collect::<Vec<_>>(),
        [0, 1, 2]
    );
    assert!(kept[2].1 >= 1_500 && kept[2].1 < 2_500, "{kept:?}");
}

#[test]
fn the_containment_estimate_follows_the_shared_share() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = shard(directory.path());
    let seeder = Seeder::new(parameters().scheme()).unwrap();
    let kmers = QueryKmers::new(&seeder, &query, false);
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let outcome =
        sketch_stage(shard.sketch(), &kmers, &bound(5_000, 0), UNLIMITED, control).unwrap();
    let estimate = |entry: u32| {
        outcome
            .candidates
            .iter()
            .find(|candidate| candidate.entry == entry)
            .map(|candidate| candidate.containment_estimate)
            .unwrap()
    };
    assert!(estimate(0) > 0.8, "{}", estimate(0));
    assert!(estimate(2) > 0.15 && estimate(2) < 0.45, "{}", estimate(2));
    assert!(estimate(3) < 0.05, "{}", estimate(3));
}

#[test]
fn frequent_seeds_count_as_possibly_hit_in_every_entry() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = shard(directory.path());
    // Seeds held by more than one entry are frequent; they keep every entry possible over the
    // stretches they cover, so the partial entry and the unrelated one cannot be left out there.
    let limits = SketchLimits {
        frequent_seed_entries: 1,
        max_seed_hits: u64::MAX,
    };
    let kept = run(&shard, &query, 2_500, limits);
    assert!(kept.iter().any(|(entry, _)| *entry == 0));
    // With no room for seed hits every seed is frequent, and every entry keeps the full bound.
    let limits = SketchLimits {
        frequent_seed_entries: u64::MAX,
        max_seed_hits: 0,
    };
    let kept = run(&shard, &query, 2_500, limits);
    assert_eq!(kept.len(), 4);
    assert!(kept.iter().all(|(_, bound)| *bound == 5_000));
}

#[test]
fn a_stopped_search_ends_stage_1() {
    let directory = tempfile::tempdir().unwrap();
    let (shard, query) = shard(directory.path());
    let seeder = Seeder::new(parameters().scheme()).unwrap();
    let kmers = QueryKmers::new(&seeder, &query, false);
    let stop = AtomicBool::new(true);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let outcome = sketch_stage(shard.sketch(), &kmers, &bound(5_000, 0), UNLIMITED, control);
    assert_eq!(outcome.unwrap_err(), Stopped::Asked);
}
