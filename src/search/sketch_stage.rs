//! Stage 1: seed hits per entry from the sketch part, the containment estimate and the coverage
//! bound.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::bound::CoverageBound;
use super::control::SearchControl;
use super::error::SearchError;
use super::query::QueryKmers;
use crate::index::{ColorSetId, ColorSets, SeedFile, StoredSeed, seed_key};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::time::Instant;

/// Distinct query k-mers that one parallel step looks up or expands.
const CHUNK_SEEDS: usize = 4096;

/// The limits of stage 1.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SketchLimits {
    /// A seed held by more entries than this is not expanded.
    pub frequent_seed_entries: u64,
    /// The most seed hits (pairs of entry and distinct query k-mer) held at once.
    pub max_seed_hits: u64,
}

/// A candidate entry after stage 1.
#[derive(Clone, Debug, PartialEq)]
pub struct PossibleEntry {
    pub entry: u32,
    /// The highest query coverage, in bases, that the entry can still reach.
    pub coverage_bound: u64,
    /// The containment estimate from the seed hits of the sketch part; for an entry in a frequent
    /// seed's color set, that seed counts as hit.
    pub containment_estimate: f64,
    /// The distinct query k-mers (indices into `QueryKmers::distinct`) that the sketch part places
    /// in the entry, ascending. Frequent seeds are not listed.
    pub seed_hits: Vec<u32>,
}

/// The work of stage 1.
#[derive(Clone, Copy, Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct SketchWork {
    /// Query bases in low-complexity stretches, which count toward no query coverage.
    pub low_complexity_bases: u64,
    /// Query k-mers, counting repeats.
    pub query_kmers: u64,
    /// Distinct query k-mers looked up in the sketch part.
    pub lookups: u64,
    /// Distinct query k-mers that the sketch part holds.
    pub seeds_found: u64,
    /// Seeds not expanded, as they are frequent or over the memory limit.
    pub frequent_seeds: u64,
    /// Seeds whose color sets were expanded.
    pub expanded_seeds: u64,
    /// Seed hits from the expanded seeds (pairs of entry and distinct query k-mer).
    pub seed_hits: u64,
    /// Entries with at least one expanded seed hit.
    pub seeded_entries: u64,
    /// Entries left out as they have no seed hit and are in no frequent seed's color set.
    pub unseeded_entries: u64,
    /// Candidate entries kept by the coverage bound.
    pub possible_entries: u64,
    /// Wall seconds of the stage.
    pub seconds: f64,
    /// A proven upper bound on the chance that a present piece of the required coverage shows no
    /// seed hit, under independent substitutions, for one continuous piece without indels.
    pub miss_upper_bound: f64,
    /// The modeled total chance of missing a qualifying entry: `miss_upper_bound` plus the chance
    /// the coverage pruning spends; at most --max-miss while the pruning is on.
    pub total_miss: f64,
    /// Whether `miss_upper_bound` leaves no room under --max-miss, so the pruning was off.
    pub miss_over_budget: bool,
}

/// What stage 1 hands on: the candidate entries, highest coverage bound first, the seeds that
/// count as possibly hit in every entry, and where the found seeds lie.
#[derive(Clone, Debug)]
pub struct SketchOutcome {
    pub candidates: Vec<PossibleEntry>,
    /// The frequent seeds (indices into `QueryKmers::distinct`), ascending.
    pub frequent_seeds: Vec<u32>,
    pub found: FoundSeeds,
    pub work: SketchWork,
}

/// The distinct query k-mers that the sketch part holds, with their stored seeds and color sets.
#[derive(Clone, Debug, Default)]
pub struct FoundSeeds {
    /// Pairs of distinct query k-mer and stored seed, ascending by k-mer.
    seeds: Vec<(u32, StoredSeed)>,
    /// The color sets of the stored seeds.
    pub color_sets: ColorSets,
}

impl FoundSeeds {
    /// The stored seed of a distinct query k-mer, if the sketch part holds it.
    pub fn stored(&self, seed: u32) -> Option<StoredSeed> {
        let index = self.seeds.binary_search_by_key(&seed, |(seed, _)| *seed);
        index.ok().map(|index| self.seeds[index].1)
    }

    /// Distinct query k-mers that the sketch part holds.
    pub fn len(&self) -> usize {
        self.seeds.len()
    }

    /// Whether no query k-mer is held.
    pub fn is_empty(&self) -> bool {
        self.seeds.is_empty()
    }
}

/// A distinct query k-mer that the sketch part holds.
#[derive(Clone, Copy, Debug)]
struct FoundSeed {
    seed: u32,
    color_set: ColorSetId,
    entries: u64,
}

/// Runs stage 1 for a query against the sketch part of a shard.
pub fn sketch_stage(
    seed_file: &SeedFile,
    query: &QueryKmers,
    bound: &CoverageBound,
    limits: SketchLimits,
    control: SearchControl<'_>,
) -> Result<SketchOutcome, SearchError> {
    let started = Instant::now();
    let found = find_query_seeds(seed_file, query, control)?;
    let mut work = SketchWork {
        query_kmers: query.kmers(),
        lookups: query.distinct().len() as u64,
        seeds_found: found.len() as u64,
        ..SketchWork::default()
    };
    let color_sets = &found.color_sets;
    let with_sizes = (found.seeds.iter())
        .map(|&(seed, stored)| {
            let entries = color_sets.entries(stored.color_set).map_or(0, <[u32]>::len);
            FoundSeed {
                seed,
                color_set: stored.color_set,
                entries: entries as u64,
            }
        })
        .collect();
    let (expanded, frequent) = split_frequent(with_sizes, limits);
    let mut frequent_seeds: Vec<u32> = frequent.iter().map(|found| found.seed).collect();
    frequent_seeds.sort_unstable();
    work.frequent_seeds = frequent_seeds.len() as u64;
    work.expanded_seeds = expanded.len() as u64;
    let seed_hits = expand(color_sets, &expanded, control)?;
    work.seed_hits = seed_hits.len() as u64;
    control.check()?;
    let entries = seed_file.entry_count() as u32;
    let (held, held_starts) = frequent_starts(color_sets, &frequent, query, entries);
    let frequent = Frequent {
        positions: query.sorted_positions(&frequent_seeds),
        held,
        held_starts,
    };
    let mut candidates: Vec<PossibleEntry> = seed_hits
        .par_chunk_by(|a, b| a.0 == b.0)
        .filter_map(|hits| {
            let (entry, seeds) = (hits[0].0, hits.iter().map(|&(_, seed)| seed).collect());
            possible_entry(entry, seeds, query, bound, &frequent)
        })
        .collect();
    let with_hits: Vec<u32> = seed_hits
        .chunk_by(|a, b| a.0 == b.0)
        .map(|h| h[0].0)
        .collect();
    work.seeded_entries = with_hits.len() as u64;
    // An entry in a frequent color set without an expanded seed hit has only the frequent seeds.
    // Entries in neither are left out.
    let only_frequent: Vec<u32> = (0..entries)
        .filter(|&entry| frequent.holds(entry) && with_hits.binary_search(&entry).is_err())
        .collect();
    work.unseeded_entries =
        u64::from(entries) - with_hits.len() as u64 - only_frequent.len() as u64;
    candidates.extend(
        only_frequent
            .into_iter()
            .filter_map(|entry| possible_entry(entry, Vec::new(), query, bound, &frequent)),
    );
    candidates.sort_by(|a, b| (b.coverage_bound, a.entry).cmp(&(a.coverage_bound, b.entry)));
    work.possible_entries = candidates.len() as u64;
    work.seconds = started.elapsed().as_secs_f64();
    Ok(SketchOutcome {
        candidates,
        frequent_seeds,
        found,
        work,
    })
}

/// Looks up every distinct query k-mer in the sketch part, in parallel chunks of neighbouring
/// seed keys, and reads the color sets of the stored seeds found.
pub fn find_query_seeds(
    seed_file: &SeedFile,
    query: &QueryKmers,
    control: SearchControl<'_>,
) -> Result<FoundSeeds, SearchError> {
    let mut keys: Vec<(u64, u32)> = (query.distinct().iter().zip(0..))
        .map(|(query_seed, seed)| (seed_key(query_seed.hash), seed))
        .collect();
    keys.par_sort_unstable();
    let chunks = keys
        .par_chunks(CHUNK_SEEDS)
        .map(|chunk| {
            control.check()?;
            let mut ascending: Vec<u64> = chunk.iter().map(|&(key, _)| key).collect();
            ascending.dedup();
            let stored = seed_file.find_seeds(&ascending)?;
            let found = chunk.iter().filter_map(|&(key, seed)| {
                let index = ascending.binary_search(&key).ok()?;
                Some((seed, stored[index]?))
            });
            Ok(found.collect::<Vec<_>>())
        })
        .collect::<Result<Vec<_>, SearchError>>()?;
    let mut seeds = chunks.concat();
    seeds.par_sort_unstable_by_key(|(seed, _)| *seed);
    let mut ids: Vec<ColorSetId> = seeds.iter().map(|(_, stored)| stored.color_set).collect();
    ids.par_sort_unstable();
    ids.dedup();
    let parts = ids
        .par_chunks(CHUNK_SEEDS)
        .map(|ids| {
            control.check()?;
            Ok(seed_file.read_color_sets(ids.iter().copied())?)
        })
        .collect::<Result<Vec<_>, SearchError>>()?;
    Ok(FoundSeeds {
        seeds,
        color_sets: ColorSets::concat(parts),
    })
}

/// Splits the found seeds into those to expand, rarest first, and the frequent ones: those held
/// by more entries than the frequent-seed limit, and those over the memory limit of seed hits.
fn split_frequent(
    mut found: Vec<FoundSeed>,
    limits: SketchLimits,
) -> (Vec<FoundSeed>, Vec<FoundSeed>) {
    found.sort_unstable_by_key(|found| (found.entries, found.seed));
    let mut seed_hits = 0u64;
    let fits = found
        .iter()
        .take_while(|found| {
            seed_hits += found.entries;
            found.entries <= limits.frequent_seed_entries && seed_hits <= limits.max_seed_hits
        })
        .count();
    let frequent = found.split_off(fits);
    (found, frequent)
}

/// What stage 1 knows of the frequent seeds for every entry.
struct Frequent {
    /// The query starts of all frequent seeds, ascending; an entry in any frequent color set
    /// counts them all as possibly hit, which can only raise its bound.
    positions: Vec<u32>,
    /// For each entry, the index in `held_starts` of the query starts of the frequent seeds
    /// whose color sets hold it, or `NOT_HELD`.
    held: Vec<u32>,
    /// The distinct lists of query starts that entries hold, each ascending without repeats.
    held_starts: Vec<Vec<u32>>,
}

/// The `held` value of an entry that no frequent color set holds.
const NOT_HELD: u32 = u32::MAX;

impl Frequent {
    fn holds(&self, entry: u32) -> bool {
        self.held
            .get(entry as usize)
            .is_some_and(|&held| held != NOT_HELD)
    }

    /// The query starts of the frequent seeds that hold `entry`.
    fn starts(&self, entry: u32) -> &[u32] {
        self.held
            .get(entry as usize)
            .and_then(|&held| self.held_starts.get(held as usize))
            .map_or(&[], Vec::as_slice)
    }
}

/// For each entry, the query starts of the frequent seeds whose color sets hold it. Each distinct
/// color set is read once, and a query start that several frequent seeds share is listed once.
fn frequent_starts(
    color_sets: &ColorSets,
    frequent: &[FoundSeed],
    query: &QueryKmers,
    entries: u32,
) -> (Vec<u32>, Vec<Vec<u32>>) {
    let mut by_color_set: Vec<(ColorSetId, u32)> = frequent
        .iter()
        .map(|found| (found.color_set, found.seed))
        .collect();
    by_color_set.sort_unstable();
    let mut sets: Vec<Vec<u32>> = vec![Vec::new(); entries as usize];
    let mut starts_of_set: Vec<Vec<u32>> = Vec::new();
    for (index, group) in by_color_set.chunk_by(|a, b| a.0 == b.0).enumerate() {
        let seeds: Vec<u32> = group.iter().map(|&(_, seed)| seed).collect();
        starts_of_set.push(query.sorted_positions(&seeds));
        for &entry in color_sets.entries(group[0].0).into_iter().flatten() {
            if let Some(held) = sets.get_mut(entry as usize) {
                held.push(index as u32);
            }
        }
    }
    let mut ids: HashMap<Vec<u32>, u32> = HashMap::new();
    let mut held_starts: Vec<Vec<u32>> = Vec::new();
    let held = sets
        .into_iter()
        .map(|members| {
            if members.is_empty() {
                return NOT_HELD;
            }
            *ids.entry(members).or_insert_with_key(|members| {
                let mut starts: Vec<u32> = (members.iter())
                    .flat_map(|&member| starts_of_set[member as usize].iter().copied())
                    .collect();
                starts.sort_unstable();
                starts.dedup();
                held_starts.push(starts);
                held_starts.len() as u32 - 1
            })
        })
        .collect();
    (held, held_starts)
}

/// Expands the color sets of the seeds into seed hits, sorted by entry and then seed.
fn expand(
    color_sets: &ColorSets,
    seeds: &[FoundSeed],
    control: SearchControl<'_>,
) -> Result<Vec<(u32, u32)>, SearchError> {
    let chunks = seeds
        .par_chunks(CHUNK_SEEDS)
        .map(|seeds| {
            control.check()?;
            let mut hits = Vec::new();
            for found in seeds {
                if let Some(entries) = color_sets.entries(found.color_set) {
                    hits.extend(entries.iter().map(|&entry| (entry, found.seed)));
                }
            }
            Ok(hits)
        })
        .collect::<Result<Vec<_>, SearchError>>()?;
    let mut hits = chunks.concat();
    hits.par_sort_unstable();
    Ok(hits)
}

/// The candidate entry for `entry` with seed hits `seeds`, or `None` when its bound is too low or
/// it has no seed hit and is in no frequent color set (then stage 2 could build no alignment
/// window, though at low identities the bound is the whole query).
fn possible_entry(
    entry: u32,
    seeds: Vec<u32>,
    query: &QueryKmers,
    bound: &CoverageBound,
    frequent: &Frequent,
) -> Option<PossibleEntry> {
    let in_frequent = frequent.holds(entry);
    if seeds.is_empty() && !in_frequent {
        return None;
    }
    let mut positions = query.sorted_positions(&seeds);
    let mut hit_starts = positions.len() as u64;
    if in_frequent {
        positions.extend_from_slice(frequent.starts(entry));
        positions.sort_unstable();
        positions.dedup();
        hit_starts = positions.len() as u64;
        positions.extend_from_slice(&frequent.positions);
        positions.sort_unstable();
    }
    let coverage_bound = bound.bound(&positions);
    bound.keeps(coverage_bound).then(|| PossibleEntry {
        entry,
        coverage_bound,
        containment_estimate: query.containment_estimate(hit_starts),
        seed_hits: seeds,
    })
}
