//! Tests of the exact layout chance, the averages and the table of the miss model.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use std::sync::OnceLock;

/// The model of these tests, built once from 100 layouts per length.
fn model() -> &'static MissModel {
    static MODEL: OnceLock<MissModel> = OnceLock::new();
    MODEL.get_or_init(|| MissModel::with_layouts(SeedScheme::default(), 100).unwrap())
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
    let short = random_layout(60, 25, 100, &mut random);
    assert_eq!(short.len(), 1);
    assert!((24..60).contains(&short[0]));
    let long = random_layout(20_000, 25, 100, &mut random);
    assert!(long.windows(2).all(|pair| pair[0] < pair[1]));
    assert!(long.iter().all(|&end| (24..20_000).contains(&end)));
    let density = long.len() as f64 / (20_000 - 24) as f64;
    assert!((density - 2.0 / 101.0).abs() < 0.004, "{density}");
}

/// Values of an independent simulation that mutates random sequences (k 25, w 100): length,
/// identity, chance, standard error.
const SIMULATION: [(usize, f64, f64, f64); 7] = [
    (1000, 0.93, 5.19e-2, 1.1e-4),
    (1000, 0.95, 5.08e-3, 3.6e-5),
    (1000, 0.96, 7.80e-4, 1.4e-5),
    (1000, 0.97, 5.62e-5, 3.7e-6),
    (1000, 0.98, 1.25e-6, 5.6e-7),
    (2000, 0.93, 2.16e-3, 3.3e-5),
    (2000, 0.95, 2.20e-5, 3.3e-6),
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
        for (&&(_, identity, chance, error), mine) in points.iter().zip(row) {
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
fn table_interpolation_holds() {
    let model = model();
    assert_eq!((model.k(), model.w()), (25, 100));
    // A point of the table is returned as stored.
    let on_grid = average_row(1000, 25, 100, 100, 0, &[0.95])[0];
    assert!((model.miss_chance(1000, 0.95) / on_grid - 1.0).abs() < 1e-9);
    // Between points the interpolation stays near a direct average; 100 layouts leave some noise.
    for (length, identity) in [(1100, 0.953), (700, 0.912), (2200, 0.9)] {
        let direct = average_row(length, 25, 100, 100, 0, &[identity])[0];
        let ratio = model.miss_chance(length as u64, identity) / direct;
        assert!(
            (0.75..1.33).contains(&ratio),
            "{length} {identity}: {ratio}"
        );
    }
}

#[test]
fn edges_are_safe() {
    let model = model();
    assert!(model.miss_chance(1000, 0.95) > model.miss_chance(1500, 0.95));
    assert!(model.miss_chance(1000, 0.95) > model.miss_chance(1000, 0.97));
    assert_eq!(model.miss_chance(24, 0.99), 1.0);
    assert_eq!(model.miss_chance(1000, 0.5), 1.0);
    assert_eq!(model.miss_chance(1000, 1.0), 0.0);
    assert_eq!(model.miss_chance(1000, f64::NAN), 1.0);
    // Beyond the table the chance of its last point is kept.
    assert_eq!(
        model.miss_chance(100_000, 0.95),
        model.miss_chance(8000, 0.95)
    );
    assert_eq!(
        model.miss_chance(1000, 0.9995),
        model.miss_chance(1000, 0.999)
    );
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
    let error = MissModel::with_layouts(SeedScheme::default(), 0);
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
        MissModel::with_layouts(scheme, LAYOUTS).unwrap(),
    );
    assert_eq!((stored.k, stored.w), (built.k, built.w));
    assert_eq!(stored.lengths, built.lengths);
    assert_eq!(stored.identities, built.identities);
    assert_eq!(stored.ln_miss, built.ln_miss);
}
