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
use super::sketch_stage::PossibleEntry;
use crate::index::{SeedLocation, Shard};
use serde::{Deserialize, Serialize};
use std::ops::Range;
use std::time::Instant;

/// Seeds looked up between two checks of the search control.
const LOOKUPS_PER_CHECK: usize = 1024;

/// The work of stage 2.
#[derive(Clone, Copy, Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct LocateWork {
    /// Candidate entries that stage 2 looked at.
    pub entries: u64,
    /// Seeds looked up in the positions, over all entries.
    pub lookups: u64,
    /// Seed locations read.
    pub seed_locations: u64,
    /// Seed hits of the sketch without seed locations.
    pub false_seed_hits: u64,
    /// Entries dropped by the second pruning or for lack of seed hits.
    pub dropped_entries: u64,
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
        self.false_seed_hits += other.false_seed_hits;
        self.dropped_entries += other.dropped_entries;
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
    locations: Vec<SeedLocation>,
    positions: Vec<u32>,
    bins: Vec<u32>,
}

/// The query, bound and settings that every entry of one search shares in stage 2.
#[derive(Clone, Copy, Debug)]
pub struct LocateRequest<'a> {
    pub query: &'a QueryKmers,
    pub bound: &'a CoverageBound<'a>,
    /// The frequent seeds of stage 1, looked up in every entry.
    pub frequent_seeds: &'a [u32],
    /// Bases added on each side of a candidate bin.
    pub flank: u64,
}

/// Runs stage 2 for one candidate entry. The coverage bound is computed again from the seed hits
/// with seed locations; the bins of those seed hits, widened by the flank, give the windows.
pub fn locate_entry(
    shard: &Shard,
    request: &LocateRequest<'_>,
    possible: &PossibleEntry,
    workspace: &mut LocateWorkspace,
    work: &mut LocateWork,
    control: SearchControl<'_>,
) -> Result<LocateOutcome, SearchError> {
    let started = Instant::now();
    let outcome = locate(shard, request, possible, workspace, work, control);
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
    outcome
}

fn locate(
    shard: &Shard,
    request: &LocateRequest<'_>,
    possible: &PossibleEntry,
    workspace: &mut LocateWorkspace,
    work: &mut LocateWork,
    control: SearchControl<'_>,
) -> Result<LocateOutcome, SearchError> {
    read_seed_hits(shard, request, possible, workspace, work, control)?;
    let LocateWorkspace {
        positions, bins, ..
    } = workspace;
    positions.sort_unstable();
    let coverage_bound = request.bound.bound(positions);
    if !request.bound.keeps(coverage_bound) {
        return Ok(LocateOutcome::Dropped(DropReason::CoverageBound));
    }
    if bins.is_empty() {
        return Ok(LocateOutcome::Dropped(DropReason::NoSeedHits));
    }
    bins.sort_unstable();
    bins.dedup();
    let contig_starts = shard.positions().contig_starts(possible.entry)?;
    let bin_bits = shard.positions().bin_bits();
    Ok(LocateOutcome::Located(LocatedEntry {
        entry: possible.entry,
        coverage_bound,
        containment_estimate: request.query.containment_estimate(positions.len() as u64),
        alignment_windows: alignment_windows(&contig_starts, bins, bin_bits, request.flank),
    }))
}

/// Reads the seed locations of the entry's seed hits and of the frequent seeds. The positions keep
/// more hash bits than the fingerprint, so a seed hit without seed locations is a false seed hit.
/// Fills the workspace with the query starts and the bins of the seed hits that have locations.
fn read_seed_hits(
    shard: &Shard,
    request: &LocateRequest<'_>,
    possible: &PossibleEntry,
    workspace: &mut LocateWorkspace,
    work: &mut LocateWork,
    control: SearchControl<'_>,
) -> Result<(), SearchError> {
    workspace.positions.clear();
    workspace.bins.clear();
    let expanded = possible.seed_hits.iter().map(|&seed| (seed, true));
    let frequent = request.frequent_seeds.iter().map(|&seed| (seed, false));
    let hashes = expanded.clone().chain(frequent.clone());
    let distinct = request.query.distinct();
    shard.positions().prefetch_seeds(
        possible.entry,
        hashes.map(|(seed, _)| distinct[seed as usize].hash),
    )?;
    for (index, (seed, from_sketch)) in expanded.chain(frequent).enumerate() {
        if index % LOOKUPS_PER_CHECK == 0 {
            control.check()?;
        }
        let hash = request.query.distinct()[seed as usize].hash;
        workspace.locations.clear();
        shard
            .positions()
            .locations(possible.entry, hash, &mut workspace.locations)?;
        work.lookups += 1;
        work.seed_locations += workspace.locations.len() as u64;
        if workspace.locations.is_empty() {
            work.false_seed_hits += u64::from(from_sketch);
            continue;
        }
        workspace
            .positions
            .extend_from_slice(request.query.positions(seed));
        workspace
            .bins
            .extend(workspace.locations.iter().map(|location| location.bin));
    }
    Ok(())
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
