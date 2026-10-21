//! Stage 2: the seed locations of each candidate entry, the second pruning and the alignment
//! windows.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::bound::CoverageBound;
use super::control::SearchControl;
use super::error::SearchError;
use super::events::DropReason;
use super::located::{ContigWindows, LocatedEntry};
use super::query::QueryKmers;
use super::sketch_stage::{FoundSeeds, PossibleEntry};
use crate::index::{SeedFile, SeedLocation, Shard, StoredSeed};
use crate::seeds::Seeder;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::ops::Range;
use std::time::Instant;

/// Seeds whose seed locations one parallel step reads.
const CHUNK_SEEDS: usize = 4096;

/// The most query starts with a located seed hit at which an entry counts as weak, so its seed
/// hits are checked on its sequence before it is aligned: the seed hits that a plan expects of a
/// qualifying copy.
const WEAK_HIT_STARTS: usize = 8;
/// Bytes of packed bases asked for ahead for one candidate entry, at most: a share of the block
/// cache, so that 128 entries ahead in each of 4 shards searched at once fill half of it, but at
/// least 1 MiB, the windows of about 256 seed bins apart. Without a cache, hints only reach
/// local files and have no bound.
const AHEAD_SHARE: u64 = 1024;
const MIN_AHEAD: u64 = 1 << 20;

/// The seed locations of one candidate entry: pairs of seed (an index into
/// `QueryKmers::distinct`) and seed location, ascending.
pub type EntryLocations = Vec<(u32, SeedLocation)>;

/// Seed locations read for several candidate entries: triples of the candidate's place among the
/// candidates, the seed and the seed location.
type CandidateLocations = Vec<(u32, u32, SeedLocation)>;

/// The work of stage 2.
#[derive(Clone, Copy, Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct LocateWork {
    /// Candidate entries that stage 2 looked at.
    pub entries: u64,
    /// Pairs of seed and candidate entry looked up in the locator part.
    pub lookups: u64,
    /// Seed locations read.
    pub seed_locations: u64,
    /// Entries dropped by the second pruning or for lack of seed hits.
    pub dropped_entries: u64,
    /// Weak entries dropped as none of their seed hits is on their sequence: chance matches of
    /// seed keys.
    #[serde(default)]
    pub unverified_entries: u64,
    /// Bins read from the sequences to check the seed hits of weak entries.
    #[serde(default)]
    pub verified_bins: u64,
    /// Bins with seed hits, over all located entries.
    pub candidate_bins: u64,
    /// Alignment windows, over all located entries.
    pub alignment_windows: u64,
    /// Seconds spent in the stage, summed over threads.
    pub seconds: f64,
}

impl LocateWork {
    /// Adds the work of another part of the search.
    pub fn add(&mut self, other: &Self) {
        self.entries += other.entries;
        self.lookups += other.lookups;
        self.seed_locations += other.seed_locations;
        self.dropped_entries += other.dropped_entries;
        self.unverified_entries += other.unverified_entries;
        self.verified_bins += other.verified_bins;
        self.candidate_bins += other.candidate_bins;
        self.alignment_windows += other.alignment_windows;
        self.seconds += other.seconds;
    }
}

/// The result of stage 2 for one entry.
#[derive(Clone, Debug, PartialEq)]
pub enum LocateOutcome {
    Located(LocatedEntry),
    Dropped(DropReason),
}

/// Buffers of stage 2, reused across the entries of one thread.
#[derive(Debug, Default)]
pub struct LocateWorkspace {
    positions: Vec<u32>,
    bins: Vec<u32>,
    bases: Vec<u8>,
    /// Bins read to check seed hits, for the work of the entry.
    verified_bins: u64,
    /// Whether the entry was dropped as none of its seed hits checked out.
    unverified: bool,
}

/// The query, bound and settings that every entry of one search shares in stage 2.
#[derive(Clone, Copy)]
pub struct LocateRequest<'a> {
    pub query: &'a QueryKmers,
    pub bound: &'a CoverageBound<'a>,
    /// The frequent seeds of stage 1, looked up in every candidate entry.
    pub frequent_seeds: &'a [u32],
    /// Bases added on each side of a candidate bin.
    pub flank: u64,
    /// The seeder of the shard, to check the seed hits of weak entries on their sequence; `None`
    /// checks none.
    pub verify: Option<&'a Seeder>,
}

/// Runs stage 2 for one candidate entry with its seed locations from `read_entry_locations`. The
/// coverage bound is computed again from the seed hits with seed locations; the bins of those
/// seed hits, widened by the flank, give the windows.
pub fn locate_entry(
    shard: &Shard,
    request: &LocateRequest<'_>,
    possible: &PossibleEntry,
    locations: &[(u32, SeedLocation)],
    workspace: &mut LocateWorkspace,
    work: &mut LocateWork,
) -> Result<LocateOutcome, SearchError> {
    let started = Instant::now();
    let outcome = locate(shard, request, possible, locations, workspace);
    work.entries += 1;
    work.seconds += started.elapsed().as_secs_f64();
    match &outcome {
        Ok(LocateOutcome::Located(located)) => {
            work.candidate_bins += workspace.bins.len() as u64;
            work.alignment_windows += located
                .alignment_windows
                .iter()
                .map(|contig| contig.ranges.len() as u64)
                .sum::<u64>();
        }
        Ok(LocateOutcome::Dropped(_)) => work.dropped_entries += 1,
        Err(_) => {}
    }
    work.verified_bins += workspace.verified_bins;
    work.unverified_entries += u64::from(workspace.unverified);
    outcome
}

fn locate(
    shard: &Shard,
    request: &LocateRequest<'_>,
    possible: &PossibleEntry,
    locations: &[(u32, SeedLocation)],
    workspace: &mut LocateWorkspace,
) -> Result<LocateOutcome, SearchError> {
    workspace.verified_bins = 0;
    workspace.unverified = false;
    let LocateWorkspace {
        positions, bins, ..
    } = workspace;
    positions.clear();
    bins.clear();
    for same in locations.chunk_by(|a, b| a.0 == b.0) {
        positions.extend_from_slice(request.query.positions(same[0].0));
        bins.extend(same.iter().map(|(_, location)| location.bin));
    }
    positions.sort_unstable();
    positions.dedup();
    let coverage_bound = request.bound.bound(positions);
    if !request.bound.keeps(coverage_bound) {
        return Ok(LocateOutcome::Dropped(DropReason::CoverageBound));
    }
    if bins.is_empty() {
        return Ok(LocateOutcome::Dropped(DropReason::NoSeedHits));
    }
    bins.sort_unstable();
    bins.dedup();
    let bin_bits = u32::from(shard.manifest().parameters.bin_bits);
    let hit_starts = positions.len() as u64;
    if let Some(seeder) = request
        .verify
        .filter(|_| hit_starts as usize <= WEAK_HIT_STARTS)
    {
        let checked = Check {
            shard,
            seeder,
            query: request.query,
            entry: possible.entry,
            bin_bits,
        };
        if !checked.any_on_sequence(locations, workspace)? {
            workspace.unverified = true;
            return Ok(LocateOutcome::Dropped(DropReason::NoSeedHits));
        }
    }
    let LocateWorkspace { bins, .. } = workspace;
    let contig_table = shard.contig_table(possible.entry)?;
    Ok(LocateOutcome::Located(LocatedEntry {
        entry: possible.entry,
        coverage_bound,
        containment_estimate: request.query.containment_estimate(hit_starts),
        alignment_windows: alignment_windows(contig_table.starts(), bins, bin_bits, request.flank),
    }))
}

/// Asks for what the later stages will read of a candidate entry with seed locations
/// `locations`: its contig table, its header page, and the bases of its seed bins with `reach`
/// bases (the flank, at least k) on each side, from the first, as far as the share of the block
/// cache lasts. These hold each bin with k bases on each side, which the seed check reads, and
/// the alignment window around it.
pub fn will_need_entry(shard: &Shard, entry: u32, locations: &[(u32, SeedLocation)], reach: u64) {
    shard.seed_file().hint_contig_table(entry);
    let bin_bits = u32::from(shard.manifest().parameters.bin_bits);
    let mut bins: Vec<u32> = locations.iter().map(|(_, location)| location.bin).collect();
    bins.sort_unstable();
    bins.dedup();
    let most = (shard.cache_ram()).map_or(u64::MAX, |ram| (ram / AHEAD_SHARE).max(MIN_AHEAD));
    let (mut ranges, mut bases, mut end) = (Vec::new(), 0u64, 0u64);
    for &bin in &bins {
        let start = u64::from(bin) << bin_bits;
        let range = start.saturating_sub(reach)..start + (1 << bin_bits) + reach;
        // Windows of close bins overlap; their bases count once. Four bases take a byte.
        bases += range.end.saturating_sub(range.start.max(end));
        end = end.max(range.end);
        if bases / 4 > most {
            break;
        }
        ranges.push(range);
    }
    shard.will_need_bases(entry, ranges);
}

/// What checking the seed hits of one weak entry on its sequence needs.
struct Check<'a> {
    shard: &'a Shard,
    seeder: &'a Seeder,
    query: &'a QueryKmers,
    entry: u32,
    bin_bits: u32,
}

impl Check<'_> {
    /// Whether any seed hit is on the sequence: whether the bin of a seed location, with k bases
    /// more on each side, holds the looked-up k-mer (with its substitutions or indel) on either
    /// strand. A seed hit whose 44-bit key only matches by chance fails, as the full hash differs.
    /// Bins are read one by one until a seed hit checks out.
    fn any_on_sequence(
        &self,
        locations: &[(u32, SeedLocation)],
        workspace: &mut LocateWorkspace,
    ) -> Result<bool, SearchError> {
        let k = u64::from(self.seeder.scheme().k);
        let packed = self.shard.packed(self.entry)?;
        let mut by_bin: Vec<(u32, u64)> = (locations.iter())
            .map(|(seed, location)| (location.bin, self.query.distinct()[*seed as usize].hash))
            .collect();
        by_bin.sort_unstable();
        by_bin.dedup();
        let mut hashes = Vec::new();
        for same in by_bin.chunk_by(|a, b| a.0 == b.0) {
            let start = u64::from(same[0].0) << self.bin_bits;
            let end = (start + (1 << self.bin_bits) + k).min(packed.bases());
            let start = start.saturating_sub(k);
            if start >= end {
                continue;
            }
            packed.read_bases(start..end, &mut workspace.bases)?;
            workspace.verified_bins += 1;
            hashes.clear();
            (self.seeder).query_kmers(&workspace.bases, |kmer| hashes.push(kmer.hash));
            hashes.sort_unstable();
            if same
                .iter()
                .any(|(_, hash)| hashes.binary_search(hash).is_ok())
            {
                return Ok(true);
            }
        }
        Ok(false)
    }
}

/// Reads the seed locations of the seed hits of all candidate entries, and of the frequent seeds
/// in every candidate entry, at once: each locator bucket is read once. Returns them by candidate.
pub fn read_entry_locations(
    seed_file: &SeedFile,
    request: &LocateRequest<'_>,
    found: &FoundSeeds,
    candidates: &[PossibleEntry],
    control: SearchControl<'_>,
    work: &mut LocateWork,
) -> Result<Vec<EntryLocations>, SearchError> {
    let started = Instant::now();
    let mut slots = vec![u32::MAX; seed_file.entry_count() as usize];
    for (slot, possible) in candidates.iter().enumerate() {
        if let Some(entry_slot) = slots.get_mut(possible.entry as usize) {
            *entry_slot = slot as u32;
        }
    }
    let mut seeds: Vec<u32> = candidates
        .iter()
        .flat_map(|possible| possible.seed_hits.iter().copied())
        .chain(request.frequent_seeds.iter().copied())
        .collect();
    seeds.par_sort_unstable();
    seeds.dedup();
    let mut stored: Vec<(StoredSeed, u32)> = seeds
        .into_iter()
        .filter_map(|seed| Some((found.stored(seed)?, seed)))
        .collect();
    stored.par_sort_unstable_by_key(|(seed, _)| (seed.bucket, seed.place));
    let chunks = stored
        .par_chunks(CHUNK_SEEDS)
        .map(|chunk| {
            control.check()?;
            read_chunk(seed_file, found, chunk, &slots)
        })
        .collect::<Result<Vec<_>, SearchError>>()?;
    let mut read = CandidateLocations::new();
    for (lookups, part) in chunks {
        work.lookups += lookups;
        read.extend(part);
    }
    work.seed_locations += read.len() as u64;
    read.par_sort_unstable();
    let mut by_candidate = vec![EntryLocations::new(); candidates.len()];
    for same in read.chunk_by(|a, b| a.0 == b.0) {
        let pairs = same.iter().map(|&(_, seed, location)| (seed, location));
        by_candidate[same[0].0 as usize] = pairs.collect();
    }
    work.seconds += started.elapsed().as_secs_f64();
    Ok(by_candidate)
}

/// Reads the seed locations of a chunk of stored seeds (sorted by bucket) in the candidate
/// entries, which `slots` maps to their place among the candidates. Returns the pairs of seed
/// and entry looked up, and (candidate, seed, seed location) triples.
fn read_chunk(
    seed_file: &SeedFile,
    found: &FoundSeeds,
    chunk: &[(StoredSeed, u32)],
    slots: &[u32],
) -> Result<(u64, CandidateLocations), SearchError> {
    let mut buckets: Vec<u32> = chunk.iter().map(|(seed, _)| seed.bucket).collect();
    buckets.dedup();
    let mut read = Vec::new();
    let mut values = Vec::new();
    let mut lookups = 0u64;
    let mut next = 0;
    seed_file.read_locator_buckets(&buckets, |bucket, part| {
        while let Some(&(stored, seed)) = chunk.get(next).filter(|(s, _)| s.bucket == bucket) {
            next += 1;
            let entries = found
                .color_sets
                .entries(stored.color_set)
                .unwrap_or_default();
            for (rank, &entry) in entries.iter().enumerate() {
                let Some(&slot) = slots.get(entry as usize).filter(|&&slot| slot != u32::MAX)
                else {
                    continue;
                };
                values.clear();
                part.locations(stored.place, rank as u32, &mut values)?;
                lookups += 1;
                let locations = values.iter().map(|&value| SeedLocation::from_value(value));
                read.extend(locations.map(|location| (slot, seed, location)));
            }
        }
        Ok(())
    })?;
    Ok((lookups, read))
}

/// The alignment windows of an entry: each candidate bin widened by `flank` on both sides, cut
/// at the ends of the contigs it touches, and merged per contig. `bins` must be ascending, and
/// contig `i` spans `contig_starts[i]..contig_starts[i + 1]` (the last start is the entry end).
fn alignment_windows(
    contig_starts: &[u64],
    bins: &[u32],
    bin_bits: u32,
    flank: u64,
) -> Vec<ContigWindows> {
    let mut windows: Vec<ContigWindows> = Vec::new();
    for &bin in bins {
        let start = u64::from(bin) << bin_bits;
        let end = start + (1u64 << bin_bits);
        // The first contig that ends after the bin start; a bin can span several contigs.
        let first = contig_starts.partition_point(|&contig_start| contig_start <= start);
        for contig in first.saturating_sub(1)..contig_starts.len().saturating_sub(1) {
            let range = contig_starts[contig]..contig_starts[contig + 1];
            if range.start >= end {
                break;
            }
            let window = start.saturating_sub(flank).max(range.start)..(end + flank).min(range.end);
            if !window.is_empty() {
                let local = window.start - range.start..window.end - range.start;
                push_window(&mut windows, contig as u32, local);
            }
        }
    }
    windows
}

/// Appends a window, merging it with the last window of the same contig when they touch.
/// Windows arrive by ascending contig, and by ascending start within a contig.
fn push_window(windows: &mut Vec<ContigWindows>, contig: u32, window: Range<u64>) {
    match windows.last_mut() {
        Some(last) if last.contig == contig => match last.ranges.last_mut() {
            Some(previous) if window.start <= previous.end => {
                previous.end = previous.end.max(window.end);
            }
            _ => last.ranges.push(window),
        },
        _ => windows.push(ContigWindows {
            contig,
            ranges: vec![window],
        }),
    }
}
