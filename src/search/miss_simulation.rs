//! Simulated losses of seeding and the coverage bound for matches outside the miss model.
//! They take minutes, so they run only on request (`--ignored`).
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::bound::CoverageBound;
use super::query::QueryKmers;
use crate::loss::{MaxMiss, MissModel};
use crate::seeds::{SeedScheme, Seeder};
use crate::test_support::{random_bases, rng, with_substitutions};
use rand::RngExt;
use rand::rngs::StdRng;
use rayon::prelude::*;
use std::collections::HashSet;

/// How the entry's copy of the query differs from it.
#[derive(Clone, Copy, Debug)]
enum Change {
    /// Independent substitutions; the matching piece covers the required half of the query.
    Substitutions,
    /// As `Substitutions`, with the whole query present.
    WholeQuery,
    /// Half of the changes are insertions or deletions of 1 to 3 bases.
    Indels,
    /// Changes come in bursts of 5 within 10 bases.
    Clustered,
    /// Changes are evenly spaced, 1 / (1 - identity) bases apart, give or take 2.
    Spread,
    /// The piece is a motif of this many bases, repeated; every copy carries the same changes.
    Repeat(usize),
    /// Present pieces of this many bases alternate with as many unrelated bases.
    Fragments(usize),
}

impl Change {
    fn name(self) -> String {
        match self {
            Self::Substitutions => "substitutions".into(),
            Self::WholeQuery => "whole-query".into(),
            Self::Indels => "indels".into(),
            Self::Clustered => "clustered".into(),
            Self::Spread => "spread".into(),
            Self::Repeat(motif) => format!("repeat-{motif}"),
            Self::Fragments(piece) => format!("fragments-{piece}"),
        }
    }
}

/// One query and its entry: the entry holds the changed piece between unrelated flanks.
fn trial(change: Change, length: usize, identity: f64, random: &mut StdRng) -> (Vec<u8>, Vec<u8>) {
    let q = 1.0 - identity;
    let half = length / 2;
    let flank = |random: &mut StdRng| random_bases(random, 500);
    let (query, piece) = match change {
        Change::WholeQuery => {
            let query = random_bases(random, length);
            let piece = with_substitutions(random, &query, identity);
            (query, piece)
        }
        Change::Fragments(size) => {
            let query = random_bases(random, length);
            let mut piece = Vec::new();
            for (index, part) in query.chunks(size).enumerate() {
                if index % 2 == 0 {
                    piece.extend(with_substitutions(random, part, identity));
                } else {
                    piece.extend(random_bases(random, part.len()));
                }
            }
            (query, piece)
        }
        Change::Repeat(motif) => {
            let unit = random_bases(random, motif);
            let changed = with_substitutions(random, &unit, identity);
            let copies = half.div_ceil(motif);
            let start = random.random_range(0..=length - half);
            let mut query = random_bases(random, length);
            let repeat: Vec<u8> = unit.repeat(copies)[..half].to_vec();
            query[start..start + half].copy_from_slice(&repeat);
            (query, changed.repeat(copies)[..half].to_vec())
        }
        _ => {
            let query = random_bases(random, length);
            let start = random.random_range(0..=length - half);
            let original = &query[start..start + half];
            let piece = match change {
                Change::Indels => with_indels(random, original, q),
                Change::Clustered => with_bursts(random, original, q),
                Change::Spread => spread(random, original, q),
                _ => with_substitutions(random, original, identity),
            };
            (query, piece)
        }
    };
    let entry = [flank(random), piece, flank(random)].concat();
    (query, entry)
}

/// A base other than `base`.
fn other(random: &mut StdRng, base: u8) -> u8 {
    let others: Vec<u8> = b"ACGT".iter().copied().filter(|&b| b != base).collect();
    others[random.random_range(0..3)]
}

/// Changes at rate `q` per base: half substitutions, a quarter insertions and a quarter
/// deletions of 1 to 3 bases.
fn with_indels(random: &mut StdRng, bases: &[u8], q: f64) -> Vec<u8> {
    let mut copy = Vec::with_capacity(bases.len());
    let mut at = 0;
    while at < bases.len() {
        if !random.random_bool(q) {
            copy.push(bases[at]);
            at += 1;
            continue;
        }
        let size = random.random_range(1..=3);
        match random.random_range(0..4) {
            0 => copy.extend(random_bases(random, size)),
            1 => at += size,
            _ => {
                copy.push(other(random, bases[at]));
                at += 1;
            }
        }
    }
    copy
}

/// Substitutions at rate `q`, in bursts of 5 at random places within 10 bases.
fn with_bursts(random: &mut StdRng, bases: &[u8], q: f64) -> Vec<u8> {
    let mut copy = bases.to_vec();
    let bursts = (q * bases.len() as f64 / 5.0).round() as usize;
    for _ in 0..bursts {
        let start = random.random_range(0..bases.len() - 10);
        for _ in 0..5 {
            let at = start + random.random_range(0..10);
            copy[at] = other(random, bases[at]);
        }
    }
    copy
}

/// Substitutions evenly spaced, `1 / q` bases apart, each moved by up to 2 bases.
fn spread(random: &mut StdRng, bases: &[u8], q: f64) -> Vec<u8> {
    let mut copy = bases.to_vec();
    let step = 1.0 / q;
    let mut place = random.random_range(0.0..step);
    while (place as usize) < bases.len() {
        let at = (place as i64 + random.random_range(-2..=2)).clamp(0, bases.len() as i64 - 1);
        copy[at as usize] = other(random, bases[at as usize]);
        place += step;
    }
    copy
}

/// The sorted query starts with a seed hit in `entry`.
fn hit_starts(seeder: &Seeder, query: &[u8], entry: &[u8], one_mismatch: bool) -> Vec<u32> {
    let mut seeds = HashSet::new();
    seeder.seeds(entry, |kmer| {
        seeds.insert(kmer.hash);
    });
    let kmers = if one_mismatch {
        QueryKmers::with_neighbours(seeder, query, false)
    } else {
        QueryKmers::new(seeder, query, false)
    };
    let hits: Vec<u32> = (0u32..)
        .zip(kmers.distinct())
        .filter(|(_, kmer)| seeds.contains(&kmer.hash))
        .map(|(seed, _)| seed)
        .collect();
    kmers.sorted_positions(&hits)
}

/// Losses of one setting: entries without a seed hit, and entries that the bound drops.
#[derive(Clone, Copy, Debug, Default)]
struct Losses {
    unseeded: u64,
    dropped: u64,
}

fn measure(change: Change, length: usize, identity: f64, one_mismatch: bool, trials: u64) {
    let scheme = SeedScheme::default();
    let seeder = Seeder::new(scheme).unwrap();
    let model = MissModel::new(scheme, usize::from(one_mismatch)).unwrap();
    let required = length as u64 / 2;
    let budget = MaxMiss::default();
    let bound = CoverageBound::new(&model, budget, identity, length as u64, false, required);
    let losses = (0..trials)
        .into_par_iter()
        .map(|index| {
            let mut random = rng(index ^ (length as u64) << 32);
            let (query, entry) = trial(change, length, identity, &mut random);
            let starts = hit_starts(&seeder, &query, &entry, one_mismatch);
            Losses {
                unseeded: u64::from(starts.is_empty()),
                dropped: u64::from(!bound.keeps(bound.bound(&starts))),
            }
        })
        .reduce(Losses::default, |a, b| Losses {
            unseeded: a.unseeded + b.unseeded,
            dropped: a.dropped + b.dropped,
        });
    println!(
        "{}\t{identity}\t{length}\t{}\t{trials}\t{}\t{}\t{:.3e}",
        change.name(),
        u8::from(one_mismatch),
        losses.unseeded,
        losses.dropped,
        bound.seed_miss_bound(),
    );
}

/// Prints, per change, identity, query length and lookup mode: trials, entries without a seed
/// hit, entries the bound drops (counted against the required half of the query), and the
/// proven chance of no seed hit. TRIALS sets the trials per setting (default 2,000).
#[test]
#[ignore]
fn adversarial_losses() {
    let trials = std::env::var("TRIALS").map_or(2_000, |value| value.parse().unwrap());
    println!("change\tidentity\tlength\tone_mismatch\ttrials\tunseeded\tdropped\tmiss_bound");
    let changes = [
        Change::Substitutions,
        Change::WholeQuery,
        Change::Indels,
        Change::Clustered,
        Change::Spread,
        Change::Repeat(100),
        Change::Repeat(300),
        Change::Fragments(50),
        Change::Fragments(200),
    ];
    for identity in [0.95, 0.8] {
        for length in [1_000, 2_000, 5_000, 10_000] {
            for change in changes {
                for one_mismatch in [false, true] {
                    measure(change, length, identity, one_mismatch, trials);
                }
            }
        }
    }
}

/// The repeated 100-base motif of the review: a 98% copy whose every stored seed carries the
/// change, so neither lookup mode finds a hit.
#[test]
fn repeated_motif_unseeded() {
    let motif = b"TGCCTAAACCTATAAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCTTGCGACAGTGACGCTTTCGCCGT";
    let reference = motif.repeat(20);
    let query = [&b"AT"[..], &motif[2..]].concat().repeat(20);
    let seeder = Seeder::new(SeedScheme::default()).unwrap();
    for one_mismatch in [false, true] {
        assert!(hit_starts(&seeder, &query, &reference, one_mismatch).is_empty());
    }
}

/// Prints the share of random query stretches, inside longer random contigs at independent
/// substitutions, whose entry shows no exact seed hit, with its standard error: the values that
/// `loss::tests::simulation_agrees` compares the table with.
#[test]
#[ignore]
fn stretch_miss_shares() {
    let trials: u64 = std::env::var("TRIALS").map_or(400_000, |value| value.parse().unwrap());
    let seeder = Seeder::new(SeedScheme::default()).unwrap();
    for (length, identity) in [
        (1000, 0.93),
        (1000, 0.95),
        (1000, 0.96),
        (2000, 0.93),
        (2000, 0.95),
    ] {
        let missed: u64 = (0..trials)
            .into_par_iter()
            .map(|index| {
                let mut random = rng(index ^ 0x5157 << 40 ^ (length as u64) << 20);
                let (query, entry) = trial(Change::WholeQuery, length, identity, &mut random);
                u64::from(hit_starts(&seeder, &query, &entry, false).is_empty())
            })
            .sum();
        let share = missed as f64 / trials as f64;
        let error = (share * (1.0 - share) / trials as f64).sqrt();
        println!("({length}, {identity}, {share:.3e}, {error:.1e}),");
    }
}
