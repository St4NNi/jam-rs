//! Stage 1 of the screen: which entries share enough seeds with which contigs of a batch.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::ScreenError;
use super::batch::Batch;
use super::key_filter::KeyFilter;
use crate::index::{ColorSetId, ColorSets, SeedFile, seed_key};
use crate::search::QueryKmers;
use crate::search::SearchControl;
use crate::seeds::Seeder;
use rayon::prelude::*;
use std::collections::HashMap;

/// The part of the expected seeds that a contig or a set of contigs must share with an entry to
/// go on: seeds spread unevenly, some differences are insertions, and a contig may end inside a
/// bin.
pub(super) const SEEDS_PART: f64 = 0.6;

/// With a closest margin, a piece is kept only for the entries that share at least the seeds that
/// an entry at the margin from its best one would, less this part and three standard deviations
/// of the count.
const CLOSE_SLACK: f64 = 0.075;

/// The seeds a stretch of `bases` bases shares with an entry at `identity`: the seed density
/// 2 / (w + 1) times the chance that a k-mer has no difference.
pub(super) fn expected_seeds(seeder: &Seeder, bases: f64, identity: f64) -> f64 {
    let scheme = seeder.scheme();
    bases * 2.0 / (f64::from(scheme.w) + 1.0) * identity.powi(i32::from(scheme.k))
}

/// A distinct k-mer of one piece of a batch that the sketch part holds.
#[derive(Clone, Copy, Debug)]
pub(super) struct PieceSeed {
    pub hash: u64,
    pub piece: u32,
    pub color_set: ColorSetId,
}

/// An entry that a piece of a batch counts for, with the seeds they share.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct Surveyed {
    pub entry: u32,
    /// The piece of the batch.
    pub piece: u32,
    /// The contig of the assembly.
    pub contig: u32,
    pub seeds: u64,
}

/// The limits of the survey.
#[derive(Clone, Copy, Debug)]
pub(super) struct SurveyLimits {
    /// The identity at which seeds are counted: the least that is asked for.
    pub identity: f64,
    /// The share of a contig that must align to an entry for the contig to count.
    pub min_contig_coverage: f64,
    pub closest_margin: Option<f64>,
}

/// Hashes the k-mers of each piece of a batch and looks up in the sketch part those that `filter`
/// may hold. Returns the distinct k-mers that it holds, sorted by piece and then color set, and
/// their color sets.
pub(super) fn found_seeds(
    seed_file: &SeedFile,
    filter: &KeyFilter,
    seeder: &Seeder,
    batch: &Batch,
    control: SearchControl<'_>,
) -> Result<(Vec<PieceSeed>, ColorSets), ScreenError> {
    let pieces = batch
        .pieces
        .par_iter()
        .enumerate()
        .map(|(index, piece)| {
            control.check()?;
            let mut hashes = Vec::new();
            let bases = &batch.bases[piece.start as usize..piece.end() as usize];
            seeder.query_kmers(bases, |kmer| {
                let key = seed_key(kmer.hash);
                if filter.may_hold(key) {
                    hashes.push((key, kmer.hash));
                }
            });
            hashes.sort_unstable();
            hashes.dedup();
            let keys: Vec<u64> = hashes.iter().map(|&(key, _)| key).collect();
            let stored = seed_file.find_seeds(&keys)?;
            let mut found: Vec<PieceSeed> = (hashes.iter().zip(stored))
                .filter_map(|(&(_, hash), stored)| {
                    Some(PieceSeed {
                        hash,
                        piece: index as u32,
                        color_set: stored?.color_set,
                    })
                })
                .collect();
            found.sort_unstable_by_key(|seed| (seed.color_set.0, seed.hash));
            Ok(found)
        })
        .collect::<Result<Vec<_>, ScreenError>>()?;
    let found = pieces.concat();
    let color_sets = seed_file.read_color_sets(found.iter().map(|seed| seed.color_set))?;
    Ok((found, color_sets))
}

/// The entries that each piece of the batch counts for. A piece counts for an entry only when it
/// shares enough seeds for the share that the contig rule requires, so a contig that shares a
/// gene with many entries counts for none of them.
pub(super) fn survey(
    seeder: &Seeder,
    batch: &Batch,
    (found, color_sets): &(Vec<PieceSeed>, ColorSets),
    limits: SurveyLimits,
    control: SearchControl<'_>,
) -> Result<Vec<Surveyed>, ScreenError> {
    let pieces: Vec<&[PieceSeed]> = found.chunk_by(|a, b| a.piece == b.piece).collect();
    let surveyed = pieces
        .into_par_iter()
        .map(|found| {
            control.check()?;
            let piece = found[0].piece;
            let shared = shared_seeds(color_sets, found);
            let length = batch.pieces[piece as usize].length as f64;
            let counted = counted_entries(seeder, shared, length, limits);
            let contig = batch.pieces[piece as usize].contig;
            Ok(counted
                .into_iter()
                .map(|(entry, seeds)| Surveyed {
                    entry,
                    piece,
                    contig,
                    seeds,
                })
                .collect::<Vec<_>>())
        })
        .collect::<Result<Vec<_>, ScreenError>>()?;
    Ok(surveyed.concat())
}

/// The seeds of one piece that each entry holds, by the color sets of the stored seeds. Their
/// keys keep enough hash bits that false seed hits are rare.
fn shared_seeds(color_sets: &ColorSets, found: &[PieceSeed]) -> HashMap<u32, u64> {
    let mut shared: HashMap<u32, u64> = HashMap::new();
    for same in found.chunk_by(|a, b| a.color_set == b.color_set) {
        for &entry in color_sets.entries(same[0].color_set).unwrap_or_default() {
            *shared.entry(entry).or_default() += same.len() as u64;
        }
    }
    shared
}

/// The entries that a piece of `length` bases counts for, from the seeds it shares with each.
fn counted_entries(
    seeder: &Seeder,
    shared: HashMap<u32, u64>,
    length: f64,
    limits: SurveyLimits,
) -> Vec<(u32, u64)> {
    let needed =
        SEEDS_PART * limits.min_contig_coverage * expected_seeds(seeder, length, limits.identity);
    let mut counted: Vec<(u32, u64)> = shared
        .into_iter()
        .filter(|&(_, seeds)| seeds as f64 >= needed)
        .collect();
    if let Some(margin) = limits.closest_margin {
        // An entry that holds the piece several times (a multimer) shares more seeds with it, but
        // is not closer: count at most the seeds that the piece has.
        let cap = expected_seeds(seeder, length, 1.0);
        let capped = |seeds: u64| (seeds as f64).min(cap);
        let best = counted
            .iter()
            .map(|&(_, seeds)| capped(seeds))
            .fold(0.0, f64::max);
        // A difference of `margin` in identity spoils this share of the seeds.
        let ratio = expected_seeds(seeder, 1.0, 1.0 - margin) / expected_seeds(seeder, 1.0, 1.0);
        let floor = (ratio - CLOSE_SLACK).max(0.0) * best - 3.0 * best.sqrt();
        counted.retain(|&(_, seeds)| capped(seeds) >= floor);
    }
    counted.sort_unstable();
    counted
}

/// The seed hits (indices into `QueryKmers::distinct` of the batch) of each kept entry, from the
/// pieces it is kept for. `kept` maps an entry to its kept pieces, ascending.
pub(super) fn kept_seed_hits(
    kmers: &QueryKmers,
    (found, color_sets): &(Vec<PieceSeed>, ColorSets),
    kept: &HashMap<u32, Vec<u32>>,
) -> Vec<(u32, Vec<u32>)> {
    let seed_of = |hash: u64| {
        kmers
            .distinct()
            .binary_search_by_key(&hash, |seed| seed.hash)
            .ok()
            .map(|seed| seed as u32)
    };
    let mut by_color_set: Vec<&PieceSeed> = found.iter().collect();
    by_color_set.sort_unstable_by_key(|seed| seed.color_set.0);
    let mut seed_hits: HashMap<u32, Vec<u32>> = HashMap::new();
    for same in by_color_set.chunk_by(|a, b| a.color_set == b.color_set) {
        let entries: Vec<u32> = color_sets
            .entries(same[0].color_set)
            .unwrap_or_default()
            .iter()
            .copied()
            .filter(|entry| kept.contains_key(entry))
            .collect();
        for entry in entries {
            let pieces = &kept[&entry];
            seed_hits.entry(entry).or_default().extend(
                same.iter()
                    .filter(|seed| pieces.binary_search(&seed.piece).is_ok())
                    .filter_map(|seed| seed_of(seed.hash)),
            );
        }
    }
    let mut seed_hits: Vec<(u32, Vec<u32>)> = seed_hits
        .into_iter()
        .filter(|(_, seeds)| !seeds.is_empty())
        .collect();
    for (_, seeds) in &mut seed_hits {
        seeds.sort_unstable();
        seeds.dedup();
    }
    seed_hits.sort_unstable_by_key(|(entry, _)| *entry);
    seed_hits
}
