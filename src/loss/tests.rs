//! Tests of the exact layout chance, the averages and the table of the miss model.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::test_support::miss_model;

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

#[test]
fn layouts_follow_windows() {
    let mut random = Random(7);
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
    [1000, 2000].into_par_iter().for_each(|length| {
        let points: Vec<_> = SIMULATION
            .iter()
            .filter(|point| point.0 == length)
            .collect();
        let identities: Vec<f64> = points.iter().map(|point| point.1).collect();
        let row = average_row(length, 25, 100, 5000, &identities);
        for (&&(_, identity, chance, error), mine) in points.iter().zip(row) {
            assert!(
                (mine - chance).abs() < 2.0 * error,
                "{length} {identity}: {mine:e} against {chance:e}"
            );
        }
    });
}

#[test]
fn table_interpolation_holds() {
    let model = miss_model();
    assert_eq!((model.k(), model.window()), (21, 20));
    // A point of the table is returned as stored.
    let on_grid = average_row(1000, 21, 20, 200, &[0.95])[0];
    assert!((model.miss_chance(1000, 0.95) / on_grid - 1.0).abs() < 1e-9);
    // Between points the interpolation stays near the direct average of the same layouts.
    for (length, identity) in [(1100, 0.953), (700, 0.912), (2200, 0.9)] {
        let direct = average_row(length, 21, 20, 200, &[identity])[0];
        let ratio = model.miss_chance(length as u64, identity) / direct;
        assert!(
            (0.85..1.15).contains(&ratio),
            "{length} {identity}: {ratio}"
        );
    }
}

#[test]
fn edges_are_safe() {
    let model = miss_model();
    assert!(model.miss_chance(1000, 0.95) > model.miss_chance(1500, 0.95));
    assert!(model.miss_chance(1000, 0.95) > model.miss_chance(1000, 0.97));
    assert_eq!(model.miss_chance(20, 0.99), 1.0);
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
    assert_eq!(
        MissModel::new(0, 100).unwrap_err(),
        LossError::InvalidSeeding
    );
    assert_eq!(
        MissModel::new(33, 100).unwrap_err(),
        LossError::InvalidSeeding
    );
    assert_eq!(
        MissModel::new(25, 0).unwrap_err(),
        LossError::InvalidSeeding
    );
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
