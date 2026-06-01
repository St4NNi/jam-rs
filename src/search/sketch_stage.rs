//! Stage 1: seed hits per entry from the sketch, the containment estimate and the coverage bound.
//!
//! Every distinct query k-mer is looked up in the sketch. A k-mer whose color set holds more than
//! the frequent-seed limit is not expanded: its query starts count as possibly hit in every entry,
//! which can only raise bounds. The other color sets are expanded rare seeds first, until the
//! memory limit of seed hits is reached; the seeds left over count as frequent too. Entries whose
//! coverage bound stays below the requested query coverage are left out; the others become
//! candidate entries, ordered by coverage bound.

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
    let (expanded, mut frequent_seeds) = split_frequent(found, limits);
    frequent_seeds.sort_unstable();
    work.frequent_seeds = frequent_seeds.len() as u64;
    work.expanded_seeds = expanded.len() as u64;
    let seed_hits = expand(sketch, &expanded, control)?;
    work.seed_hits = seed_hits.len() as u64;
    control.check()?;
    let frequent_positions = query.sorted_positions(&frequent_seeds);
    let mut candidates: Vec<PossibleEntry> = seed_hits
        .par_chunk_by(|a, b| a.0 == b.0)
        .filter_map(|hits| {
            let seeds: Vec<u32> = hits.iter().map(|&(_, seed)| seed).collect();
            possible_entry(hits[0].0, seeds, query, bound, &frequent_positions)
        })
        .collect();
    let with_hits: Vec<u32> = seed_hits
        .chunk_by(|a, b| a.0 == b.0)
        .map(|h| h[0].0)
        .collect();
    work.entries_with_seed_hits = with_hits.len() as u64;
    // An entry without an expanded seed hit has only the frequent seeds; their bound is shared.
    if let Some(shared) = possible_entry(0, Vec::new(), query, bound, &frequent_positions) {
        let entries = sketch.entries() as u32;
        candidates.extend(
            (0..entries)
                .filter(|entry| with_hits.binary_search(entry).is_err())
                .map(|entry| PossibleEntry {
                    entry,
                    ..shared.clone()
                }),
        );
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

/// Splits the found seeds into those to expand, rarest first, and the frequent ones.
fn split_frequent(mut found: Vec<FoundSeed>, limits: SketchLimits) -> (Vec<FoundSeed>, Vec<u32>) {
    found.sort_unstable_by_key(|found| (found.entries, found.seed));
    let mut seed_hits = 0u64;
    let fits = found
        .iter()
        .take_while(|found| {
            seed_hits += found.entries;
            found.entries <= limits.frequent_seed_entries && seed_hits <= limits.max_seed_hits
        })
        .count();
    let frequent = found[fits..].iter().map(|found| found.seed).collect();
    found.truncate(fits);
    (found, frequent)
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

/// The candidate entry for `entry` with seed hits `seeds`, or `None` when its bound is too low.
fn possible_entry(
    entry: u32,
    seeds: Vec<u32>,
    query: &QueryKmers,
    bound: &CoverageBound,
    frequent_positions: &[u32],
) -> Option<PossibleEntry> {
    let hit_positions = query.sorted_positions(&seeds);
    let containment_estimate = query.containment_estimate(hit_positions.len() as u64);
    let mut positions = hit_positions;
    if !frequent_positions.is_empty() {
        positions.extend_from_slice(frequent_positions);
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
