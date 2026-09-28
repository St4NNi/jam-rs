//! Tests of the exact layout chance, the averages and the table of the miss model.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use std::sync::OnceLock;

/// The model of these tests, built once from 100 layouts per length.
fn model() -> &'static MissModel {
    static MODEL: OnceLock<MissModel> = OnceLock::new();
    MODEL.get_or_init(|| MissModel::with_layouts(SeedScheme::default(), 100, 0).unwrap())
}

/// The chance that all seeds are mutated, by trying every pattern of mutated bases.
fn brute_force(length: usize, k: usize, seed_ends: &[usize], identity: f64) -> f64 {
    let mut total = 0.0;
    for mask in 0u32..(1 << length) {
        let mutated = mask.count_ones() as i32;
        let weight = identity.powi(length as i32 - mutated) * (1.0 - identity).powi(mutated);
        let every_seed_mutated = seed_ends.iter().all(|&end| {
            let seed_mask = ((1u32 << k) - 1) << (end + 1 - k);
            mask & seed_mask != 0
        });
        if every_seed_mutated {
            total += weight;
        }
    }
    total
}

#[test]
fn layout_chance_exact() {
    for (length, k, seeds) in [
        (12, 3, vec![2, 5, 6, 11]),
        (10, 4, vec![3, 9]),
        (14, 5, vec![4, 7, 8, 13]),
        (9, 3, vec![8]),
    ] {
        for identity in [0.5, 0.8, 0.95] {
            let exact = layout_miss_chance(&seeds, k, identity);
            let brute = brute_force(length, k, &seeds, identity);
            assert!((exact - brute).abs() < 1e-12, "{seeds:?} {identity}");
        }
    }
}

/// As `brute_force`, for seeds hit with at most one mismatch.
fn brute_force_one(length: usize, k: usize, seed_ends: &[usize], identity: f64) -> f64 {
    let mut total = 0.0;
    for mask in 0u32..(1 << length) {
        let mutated = mask.count_ones() as i32;
        let weight = identity.powi(length as i32 - mutated) * (1.0 - identity).powi(mutated);
        let every_seed_missed = seed_ends.iter().all(|&end| {
            let seed_mask = ((1u32 << k) - 1) << (end + 1 - k);
            (mask & seed_mask).count_ones() >= 2
        });
        if every_seed_missed {
            total += weight;
        }
    }
    total
}

#[test]
fn one_mismatch_exact() {
    for (length, k, seeds) in [
        (12, 3, vec![2, 5, 6, 11]),
        (10, 4, vec![3, 9]),
        (14, 5, vec![4, 7, 8, 13]),
        (16, 6, vec![5, 6, 12, 15]),
        (9, 3, vec![8]),
    ] {
        for identity in [0.5, 0.8, 0.95] {
            let exact = one_mismatch_chance(&seeds, k, identity);
            let brute = brute_force_one(length, k, &seeds, identity);
            assert!(
                (exact - brute).abs() < 1e-12,
                "{seeds:?} {identity} {exact} {brute}"
            );
            assert!(exact <= layout_miss_chance(&seeds, k, identity) + 1e-15);
        }
    }
}

#[test]
fn layouts_follow_windows() {
    let mut random = Pcg64Mcg::seed_from_u64(7);
    // A stretch shorter than one window may lie between two seeds; one window long, never.
    let empty = (0..1000)
        .filter(|_| random_layout(60, 25, 100, &mut random).is_empty())
        .count();
    assert!((100..900).contains(&empty), "{empty}");
    assert!((0..1000).all(|_| !random_layout(124, 25, 100, &mut random).is_empty()));
    let long = random_layout(20_000, 25, 100, &mut random);
    assert!(long.windows(2).all(|pair| pair[0] < pair[1]));
    assert!(long.iter().all(|&end| (24..20_000).contains(&end)));
    let density = long.len() as f64 / (20_000 - 24) as f64;
    assert!((density - 2.0 / 101.0).abs() < 0.004, "{density}");
}

/// Values of an independent simulation that mutates random stretches inside longer random
/// contigs and seeds them with the seeder (k 25, w 100, 400,000 trials each,
/// `search::miss_simulation::stretch_miss_shares`): length, identity, chance, standard error.
const SIMULATION: [(usize, f64, f64, f64); 5] = [
    (1000, 0.93, 4.489e-2, 3.3e-4),
    (1000, 0.95, 4.000e-3, 1.0e-4),
    (1000, 0.96, 5.975e-4, 3.9e-5),
    (2000, 0.93, 1.925e-3, 6.9e-5),
    (2000, 0.95, 2.250e-5, 7.5e-6),
];

#[test]
fn simulation_agrees() {
    const LAYOUTS: usize = 5000;
    [1000, 2000].into_par_iter().for_each(|length| {
        let points: Vec<_> = SIMULATION
            .iter()
            .filter(|point| point.0 == length)
            .collect();
        let identities: Vec<f64> = points.iter().map(|point| point.1).collect();
        let row = average_row(length, 25, 100, LAYOUTS, 0, &identities);
        for (&&(_, identity, chance, error), (mine, _)) in points.iter().zip(row) {
            // Both sides are sample means: allow three times their combined standard error.
            let model_error = layout_spread(length, identity) / (LAYOUTS as f64).sqrt();
            let tolerance = 3.0 * (error.powi(2) + model_error.powi(2)).sqrt();
            assert!(
                (mine - chance).abs() < tolerance,
                "{length} {identity}: {mine:e} against {chance:e}, tolerance {tolerance:e}"
            );
        }
    });
}

/// The standard deviation of the miss chance between layouts, from an independent sample.
fn layout_spread(length: usize, identity: f64) -> f64 {
    let mut random = Pcg64Mcg::seed_from_u64(99);
    let values: Vec<f64> = (0..1000)
        .map(|_| layout_miss_chance(&random_layout(length, 25, 100, &mut random), 25, identity))
        .collect();
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / (values.len() - 1) as f64;
    variance.sqrt()
}

#[test]
fn lookups_round_down() {
    let model = model();
    assert_eq!((model.k(), model.w()), (25, 100));
    // A point of the table holds its upper value; between points the lower length and the lower
    // identity hold, so a lookup is never below the value of the table at the point itself.
    let (mean, error) = average_row(1029, 25, 100, 100, 0, &[0.95])[0];
    let ratio = model.miss_chance(1029, 0.95) / upper_value(mean, error);
    assert!((0.5..=1.0 + 1e-9).contains(&ratio), "{ratio}");
    assert_eq!(
        model.miss_chance(1130, 0.9549),
        model.miss_chance(1029, 0.95)
    );
    for length in [200, 1000, 3000] {
        assert!(model.miss_chance(length, 0.95) >= model.miss_chance(length + 1, 0.95));
        assert!(model.miss_chance(length, 0.95) >= model.miss_chance(length, 0.96));
    }
}

/// A single seed at identity 0.9525 misses with 1 - 0.9525^25 = 0.703773; interpolating the
/// logarithm between the columns 0.95 and 0.955 gave 0.702892, below it.
#[test]
fn single_seed_covered() {
    let identities = identity_grid();
    let ln_miss = (identities.iter())
        .map(|&identity| layout_miss_chance(&[24], 25, identity).ln())
        .collect();
    let single = MissModel {
        k: 25,
        w: 100,
        lengths: vec![25.0],
        identities,
        mismatches: 0,
        ln_miss,
    };
    let exact = 1.0 - 0.9525f64.powi(25);
    assert!((layout_miss_chance(&[24], 25, 0.9525) - exact).abs() < 1e-12);
    assert!(single.miss_chance(25, 0.9525) >= exact);
}

#[test]
fn edges_are_safe() {
    let model = model();
    assert!(model.miss_chance(1000, 0.95) > model.miss_chance(1500, 0.95));
    assert!(model.miss_chance(1000, 0.95) > model.miss_chance(1000, 0.97));
    assert_eq!(model.miss_chance(24, 0.99), 1.0);
    assert_eq!(model.miss_chance(1000, 0.5), 1.0);
    // At identity 1 only a stretch of a whole window (w + k - 1 = 124 bases) surely has a seed.
    assert_eq!(model.miss_chance(124, 1.0), 0.0);
    assert_eq!(model.miss_chance(123, 1.0), 1.0);
    // A shorter stretch inside a contig may lie between two seeds.
    assert!(model.miss_chance(100, 0.999) > 0.1);
    assert_eq!(model.miss_chance(1000, f64::NAN), 1.0);
    // Beyond the table the chance of its last point is kept.
    assert_eq!(
        model.miss_chance(100_000, 0.95),
        model.miss_chance(7612, 0.95)
    );
    assert_eq!(
        model.miss_chance(1000, 0.9995),
        model.miss_chance(1000, 0.999)
    );
}

/// The proven bound: q^floor(L / 124), q the chance that one seed has more mismatches than the
/// lookups allow.
#[test]
fn block_upper_bound() {
    let scheme = SeedScheme::default();
    let (exact, one) = (model(), MissModel::one_mismatch(scheme).unwrap());
    let q0 = 1.0 - 0.95f64.powi(25);
    let q1 = q0 - 25.0 * 0.05 * 0.95f64.powi(24);
    assert_eq!(exact.miss_upper_bound(123, 0.95), 1.0);
    assert!((exact.miss_upper_bound(750, 0.95) - q0.powi(6)).abs() < 1e-12);
    assert!((one.miss_upper_bound(750, 0.95) - q1.powi(6)).abs() < 1e-12);
    // The values of the review at L = 750.
    assert!((exact.miss_upper_bound(750, 0.95) - 0.142_372_25).abs() < 1e-8);
    assert!((one.miss_upper_bound(750, 0.95) - 0.002_092_00).abs() < 1e-8);
    assert_eq!(exact.miss_upper_bound(124, 1.0), 0.0);
    // The bound is never below the table's estimate where both apply.
    for length in [500, 1000, 2000, 5000] {
        assert!(exact.miss_upper_bound(length, 0.95) >= exact.miss_chance(length, 0.95));
    }
}

#[test]
fn bad_seeding_rejected() {
    for (k, w) in [(0, 100), (33, 100), (25, 0)] {
        let scheme = SeedScheme {
            k,
            w,
            ..SeedScheme::default()
        };
        let error = MissModel::new(scheme).unwrap_err();
        assert!(matches!(error, LossError::Seeds(_)), "{k} {w}");
    }
    let error = MissModel::with_layouts(SeedScheme::default(), 0, 0);
    assert_eq!(error.unwrap_err(), LossError::NoLayouts);
}

#[test]
fn default_budget() {
    assert_eq!(LossBudget::default().chance(), 0.0001);
}

#[test]
fn budget_is_chance() {
    for bad in [0.0, 1.0, -0.1, f64::NAN, 2.0] {
        assert_eq!(LossBudget::new(bad).unwrap_err(), LossError::InvalidBudget);
    }
}

#[test]
fn default_table_stored() {
    let scheme = SeedScheme::default();
    let (stored, built) = (
        MissModel::new(scheme).unwrap(),
        MissModel::with_layouts(scheme, LAYOUTS, 0).unwrap(),
    );
    assert_eq!((stored.k, stored.w), (built.k, built.w));
    assert_eq!(stored.lengths, built.lengths);
    assert_eq!(stored.identities, built.identities);
    assert_eq!(stored.ln_miss, built.ln_miss);
}

/// The stored rows of the one-mismatch table for some lengths, against a new average.
#[test]
fn mismatch_rows_stored() {
    let stored = MissModel::one_mismatch(SeedScheme::default()).unwrap();
    assert_eq!(stored.mismatches(), 1);
    let width = stored.identities.len();
    for length in [25, 271, 1029] {
        let row = stored
            .lengths
            .iter()
            .position(|&l| l == length as f64)
            .unwrap();
        let built = average_row(length, 25, 100, LAYOUTS, 1, &stored.identities);
        let built: Vec<f64> = (built.iter())
            .map(|&(mean, error)| upper_value(mean, error).max(1e-300).ln())
            .collect();
        // A stored value is the computed one, or lower when a shorter length or a lower
        // identity had a lower one.
        let row = &stored.ln_miss[row * width..(row + 1) * width];
        assert!(
            row.iter()
                .zip(&built)
                .all(|(stored, built)| stored <= built)
        );
        let equal = row
            .iter()
            .zip(&built)
            .filter(|(stored, built)| stored == built);
        assert!(equal.count() > width / 2, "{length}");
    }
}

/// The whole stored one-mismatch table against a new average. It takes minutes, so it runs
/// only on request (`--ignored`).
#[test]
#[ignore]
fn mismatch_table_stored() {
    let scheme = SeedScheme::default();
    let stored = MissModel::one_mismatch(scheme).unwrap();
    let built = MissModel::with_layouts(scheme, LAYOUTS, 1).unwrap();
    assert_eq!(stored.lengths, built.lengths);
    assert_eq!(stored.ln_miss, built.ln_miss);
}

/// Writes both stored tables of the default k and w anew (`--ignored`, takes minutes).
#[test]
#[ignore]
fn write_stored_tables() {
    let scheme = SeedScheme::default();
    let folder = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/loss");
    for (mismatches, name) in [(0, "default_table.bin"), (1, "one_mismatch_table.bin")] {
        let built = MissModel::with_layouts(scheme, LAYOUTS, mismatches).unwrap();
        let bytes: Vec<u8> = built.ln_miss.iter().flat_map(|v| v.to_le_bytes()).collect();
        std::fs::write(folder.join(name), bytes).unwrap();
    }
}
