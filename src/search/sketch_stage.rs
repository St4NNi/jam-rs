//! Stage 1: seed hits per entry from the sketch, the containment estimate and the coverage bound.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::bound::CoverageBound;
use super::control::{SearchControl, Stopped};
use super::query::QueryKmers;
use crate::index::{ColorSetId, Sketch};
use rayon::prelude::*;
use std::time::{Duration, Instant};

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
    /// The share of the query the entry holds, estimated from its expanded seed hits.
    pub containment_estimate: f64,
    /// The distinct query k-mers (indices into `QueryKmers::seeds`) that the sketch places in the
    /// entry, ascending. Frequent seeds are not listed.
    pub seed_hits: Vec<u32>,
}

/// The work of stage 1.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SketchWork {
    /// Query k-mers, counting repeats.
    pub query_kmers: u64,
    /// Distinct query k-mers looked up in the sketch.
    pub lookups: u64,
    /// Distinct query k-mers that the sketch holds.
    pub seeds_found: u64,
    /// Seeds not expanded, as they are frequent or over the memory limit.
    pub frequent_seeds: u64,
    /// Seeds whose color sets were expanded.
    pub expanded_seeds: u64,
    /// Seed hits from the expanded seeds (pairs of entry and distinct query k-mer).
    pub seed_hits: u64,
    /// Entries with at least one expanded seed hit.
    pub entries_with_seed_hits: u64,
    /// Entries left out as they have no seed hit and are in no frequent seed's color set.
    pub entries_without_seed_hits: u64,
    /// Candidate entries kept by the coverage bound.
    pub possible_entries: u64,
    /// Wall time of the stage.
    pub time: Duration,
}

/// What stage 1 hands on: the candidate entries, highest coverage bound first, and the seeds that
/// count as possibly hit in every entry.
#[derive(Clone, Debug)]
pub struct SketchOutcome {
    pub candidates: Vec<PossibleEntry>,
    /// The frequent seeds (indices into `QueryKmers::seeds`), ascending.
    pub frequent_seeds: Vec<u32>,
    pub work: SketchWork,
}

/// A distinct query k-mer that the sketch holds.
#[derive(Clone, Copy, Debug)]
struct FoundSeed {
    seed: u32,
    color_set: ColorSetId,
    entries: u64,
}

/// Runs stage 1 for a query against the sketch of a shard.
pub fn sketch_stage(
    sketch: &Sketch,
    query: &QueryKmers,
    bound: &CoverageBound,
    limits: SketchLimits,
    control: SearchControl<'_>,
) -> Result<SketchOutcome, Stopped> {
    let started = Instant::now();
    let found = look_up(sketch, query, control)?;
    let mut work = SketchWork {
        query_kmers: query.kmers(),
        lookups: query.seeds().len() as u64,
        seeds_found: found.len() as u64,
        ..SketchWork::default()
    };
    let (expanded, frequent) = split_frequent(found, limits);
    let mut frequent_seeds: Vec<u32> = frequent.iter().map(|found| found.seed).collect();
    frequent_seeds.sort_unstable();
    work.frequent_seeds = frequent_seeds.len() as u64;
    work.expanded_seeds = expanded.len() as u64;
    let seed_hits = expand(sketch, &expanded, control)?;
    work.seed_hits = seed_hits.len() as u64;
    control.check()?;
    let entries = sketch.entry_count() as u32;
    let in_frequent = FrequentUnion::new(sketch, &frequent, entries);
    let frequent_positions = query.sorted_positions(&frequent_seeds);
    let possible_hits = |entry: u32| match in_frequent.holds(entry) {
        true => frequent_positions.as_slice(),
        false => &[],
    };
    let mut candidates: Vec<PossibleEntry> = seed_hits
        .par_chunk_by(|a, b| a.0 == b.0)
        .filter_map(|hits| {
            let (entry, seeds) = (hits[0].0, hits.iter().map(|&(_, seed)| seed).collect());
            possible_entry(entry, seeds, query, bound, possible_hits(entry))
        })
        .collect();
    let with_hits: Vec<u32> = seed_hits
        .chunk_by(|a, b| a.0 == b.0)
        .map(|h| h[0].0)
        .collect();
    work.entries_with_seed_hits = with_hits.len() as u64;
    // An entry in the frequent union without an expanded seed hit has only the frequent seeds;
    // their bound is shared. Entries in neither are left out.
    let only_frequent: Vec<u32> = (0..entries)
        .filter(|&entry| in_frequent.holds(entry) && with_hits.binary_search(&entry).is_err())
        .collect();
    work.entries_without_seed_hits =
        u64::from(entries) - with_hits.len() as u64 - only_frequent.len() as u64;
    if let Some(shared) = possible_entry(0, Vec::new(), query, bound, &frequent_positions) {
        candidates.extend(only_frequent.into_iter().map(|entry| PossibleEntry {
            entry,
            ..shared.clone()
        }));
    }
    candidates.sort_by(|a, b| (b.coverage_bound, a.entry).cmp(&(a.coverage_bound, b.entry)));
    work.possible_entries = candidates.len() as u64;
    work.time = started.elapsed();
    Ok(SketchOutcome {
        candidates,
        frequent_seeds,
        work,
    })
}

/// Looks up every distinct query k-mer in the sketch, in parallel chunks.
fn look_up(
    sketch: &Sketch,
    query: &QueryKmers,
    control: SearchControl<'_>,
) -> Result<Vec<FoundSeed>, Stopped> {
    let color_sets = sketch.color_sets();
    let chunks = query
        .seeds()
        .par_chunks(CHUNK_SEEDS)
        .enumerate()
        .map(|(chunk, seeds)| {
            control.check()?;
            let first = (chunk * CHUNK_SEEDS) as u32;
            Ok(seeds
                .iter()
                .zip(first..)
                .filter_map(|(query_seed, seed)| {
                    let color_set = sketch.color_set(query_seed.hash)?;
                    let entries = color_sets.entries(color_set)?.len();
                    Some(FoundSeed {
                        seed,
                        color_set,
                        entries,
                    })
                })
                .collect::<Vec<_>>())
        })
        .collect::<Result<Vec<_>, Stopped>>()?;
    Ok(chunks.concat())
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

/// The union of the frequent seeds' color sets, as a bitmap over the shard's entries. An entry in
/// it counts every frequent seed as possibly hit, which can only raise its bound.
struct FrequentUnion {
    words: Vec<u64>,
}

impl FrequentUnion {
    /// Marks the entries of every frequent seed's color set, and stops once all are marked, so the
    /// cost is at most the sum of the frequent color sets' sizes.
    fn new(sketch: &Sketch, frequent: &[FoundSeed], entries: u32) -> Self {
        let mut words = vec![0u64; (entries as usize).div_ceil(64)];
        let mut marked = 0u32;
        let color_sets = sketch.color_sets();
        for found in frequent {
            if marked == entries {
                break;
            }
            for entry in color_sets.entries(found.color_set).into_iter().flatten() {
                let (word, bit) = (entry as usize / 64, 1u64 << (entry % 64));
                if let Some(word) = words.get_mut(word).filter(|word| **word & bit == 0) {
                    *word |= bit;
                    marked += 1;
                }
            }
        }
        Self { words }
    }

    fn holds(&self, entry: u32) -> bool {
        self.words
            .get(entry as usize / 64)
            .is_some_and(|word| word & (1u64 << (entry % 64)) != 0)
    }
}

/// Expands the color sets of the seeds into seed hits, sorted by entry and then seed.
fn expand(
    sketch: &Sketch,
    seeds: &[FoundSeed],
    control: SearchControl<'_>,
) -> Result<Vec<(u32, u32)>, Stopped> {
    let color_sets = sketch.color_sets();
    let chunks = seeds
        .par_chunks(CHUNK_SEEDS)
        .map(|seeds| {
            control.check()?;
            let mut hits = Vec::new();
            for found in seeds {
                if let Some(entries) = color_sets.entries(found.color_set) {
                    hits.extend(entries.map(|entry| (entry, found.seed)));
                }
            }
            Ok(hits)
        })
        .collect::<Result<Vec<_>, Stopped>>()?;
    let mut hits = chunks.concat();
    hits.par_sort_unstable();
    Ok(hits)
}

/// The candidate entry for `entry` with seed hits `seeds` and possible seed hits at the query
/// starts `possible_hits`, or `None` when its bound is too low or it has neither (then stage 2
/// could build no alignment window, though at low identities the bound is the whole query).
fn possible_entry(
    entry: u32,
    seeds: Vec<u32>,
    query: &QueryKmers,
    bound: &CoverageBound,
    possible_hits: &[u32],
) -> Option<PossibleEntry> {
    if seeds.is_empty() && possible_hits.is_empty() {
        return None;
    }
    let hit_positions = query.sorted_positions(&seeds);
    let containment_estimate = query.containment_estimate(hit_positions.len() as u64);
    let mut positions = hit_positions;
    if !possible_hits.is_empty() {
        positions.extend_from_slice(possible_hits);
        positions.sort_unstable();
    }
    let coverage_bound = bound.bound(&positions);
    bound.keeps(coverage_bound).then_some(PossibleEntry {
        entry,
        coverage_bound,
        containment_estimate,
        seed_hits: seeds,
    })
}
