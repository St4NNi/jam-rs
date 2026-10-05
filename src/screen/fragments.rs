//! Stages 2 and 3 of the screen for one batch: the kept entries are located and aligned with the
//! batch as the query, and their fragments are cut at the contig borders and turned around.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::ScreenError;
use super::batch::{Batch, Piece};
use super::turn::{slice_alignment, swap_roles};
use crate::align::{Alignment, AlignmentConfig, Interval};
use crate::decide::{DecisionRule, Fragment};
use crate::index::{SeedLocation, Shard};
use crate::loss::{MaxMiss, MissModel};
use crate::search::PossibleEntry;
use crate::search::QueryKmers;
use crate::search::SearchControl;
use crate::search::{
    AlignStageConfig, AlignStageWorkspace, AlignWork, PreparedQuery, entry_fragments,
};
use crate::search::CoverageBound;
use crate::search::{LocateOutcome, LocateRequest, LocateWork, LocateWorkspace, locate_entry};
use crate::search::{find_query_seeds, read_entry_locations};
use rayon::prelude::*;

/// An entry contig, the query of the screen: the entry and the contig within it.
pub(super) type EntryContig = (u32, u32);

/// A fragment of an entry contig (the query of the screen) on an assembly contig.
#[derive(Clone, Debug)]
pub(super) struct ScreenFragment {
    pub entry: u32,
    /// The contig of the entry.
    pub entry_contig: u32,
    /// Its `contig_index` is the assembly contig; its query is the entry contig.
    pub fragment: Fragment,
}

/// What the alignment of a batch needs besides the batch.
pub(super) struct BatchAlignment<'a> {
    pub shard: &'a Shard,
    pub model: &'a MissModel,
    pub rule: &'a DecisionRule,
    pub align: &'a AlignStageConfig,
    pub flank: u64,
    pub control: SearchControl<'a>,
}

/// The work of stages 2 and 3 over the batches of an assembly.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct FragmentWork {
    pub locate: LocateWork,
    pub align: AlignWork,
    /// Entries aligned, counted once per batch.
    pub aligned_entries: u64,
}

impl FragmentWork {
    pub fn add(&mut self, other: &Self) {
        self.locate.add(&other.locate);
        self.align.add(&other.align);
        self.aligned_entries += other.aligned_entries;
    }
}

impl BatchAlignment<'_> {
    /// Locates and aligns each entry with its seed hits in the batch, and turns its fragments.
    pub fn fragments(
        &self,
        batch: &Batch,
        kmers: &QueryKmers,
        seed_hits: Vec<(u32, Vec<u32>)>,
    ) -> Result<(Vec<ScreenFragment>, FragmentWork), ScreenError> {
        // The screen decides an entry by its own coverage, not by the coverage of the batch, so
        // the bound of stage 2 keeps every entry.
        let bound = CoverageBound::new(
            self.model,
            MaxMiss::default(),
            self.rule.min_identity,
            kmers.length(),
            false,
            0,
        );
        let request = LocateRequest {
            query: kmers,
            bound: &bound,
            frequent_seeds: &[],
            flank: self.flank,
        };
        let query = PreparedQuery::new(&batch.bases, false, self.align)?;
        let candidates: Vec<PossibleEntry> = (seed_hits.into_iter())
            .map(|(entry, seed_hits)| PossibleEntry {
                entry,
                coverage_bound: 0,
                containment_estimate: 0.0,
                seed_hits,
            })
            .collect();
        let seed_file = self.shard.seed_file();
        let found = find_query_seeds(seed_file, kmers, self.control)?;
        let mut work = FragmentWork::default();
        let locations = read_entry_locations(
            seed_file,
            &request,
            &found,
            &candidates,
            self.control,
            &mut work.locate,
        )?;
        let results = (candidates.par_iter().zip(&locations))
            .map_init(
                || (LocateWorkspace::default(), AlignStageWorkspace::default()),
                |(locate, align), (possible, locations)| {
                    let workspaces = (locate, align);
                    self.entry(batch, &request, &query, possible, locations, workspaces)
                },
            )
            .collect::<Result<Vec<_>, ScreenError>>()?;
        let mut turned = Vec::new();
        for (fragments, part) in results {
            turned.extend(fragments);
            work.add(&part);
        }
        Ok((turned, work))
    }

    fn entry(
        &self,
        batch: &Batch,
        request: &LocateRequest<'_>,
        query: &PreparedQuery,
        possible: &PossibleEntry,
        locations: &[(u32, SeedLocation)],
        (locate, align): (&mut LocateWorkspace, &mut AlignStageWorkspace),
    ) -> Result<(Vec<ScreenFragment>, FragmentWork), ScreenError> {
        let mut work = FragmentWork::default();
        let outcome = locate_entry(
            self.shard,
            request,
            possible,
            locations,
            locate,
            &mut work.locate,
        )?;
        let LocateOutcome::Located(located) = outcome else {
            return Ok((Vec::new(), work));
        };
        let packed = self.shard.packed(possible.entry)?;
        let (fragments, align_work) = entry_fragments(
            &packed,
            query,
            &located,
            self.align,
            self.rule,
            align,
            self.control,
        )?;
        work.align = align_work;
        work.aligned_entries = 1;
        let scores = &self.align.band.alignment;
        let mut turned = Vec::new();
        for fragment in fragments {
            for (contig, alignment) in cut_at_contigs(batch, &fragment.alignment, scores)? {
                turned.push(ScreenFragment {
                    entry: possible.entry,
                    entry_contig: fragment.contig_index,
                    fragment: Fragment {
                        contig_index: contig,
                        query_segments: vec![alignment.query_interval],
                        alignment,
                    },
                });
            }
        }
        Ok((turned, work))
    }
}

/// The parts of a batch alignment on each piece of the batch, in the coordinates of the assembly
/// contig, turned around so that the entry contig is the query. Returns the assembly contig of
/// each part. Of an alignment that overlapping pieces both hold, only one copy is kept.
pub(super) fn cut_at_contigs(
    batch: &Batch,
    alignment: &Alignment,
    scores: &AlignmentConfig,
) -> Result<Vec<(u32, Alignment)>, ScreenError> {
    let interval = alignment.query_interval;
    let first = batch
        .pieces
        .partition_point(|piece| piece.end() <= interval.start);
    let mut parts = Vec::new();
    for piece in batch.pieces[first..]
        .iter()
        .take_while(|piece| piece.start < interval.end)
    {
        let (start, end) = (
            interval.start.max(piece.start),
            interval.end.min(piece.end()),
        );
        let Some((from, to)) = owned_part(piece, start, end) else {
            continue;
        };
        let part = if (from, to) == (interval.start, interval.end) {
            Some(alignment.clone())
        } else {
            slice_alignment(alignment, from, to, scores)?
        };
        let Some(mut part) = part else {
            continue;
        };
        let shift = |position: u64| position - piece.start + piece.offset;
        part.query_interval = Interval::new(
            shift(part.query_interval.start),
            shift(part.query_interval.end),
        )?;
        parts.push((piece.contig, swap_roles(&part)));
    }
    Ok(parts)
}

/// The part of an alignment at batch positions `start..end` of `piece` that the piece keeps. A
/// part starting in the owned positions is kept whole, unless the piece cuts it before the
/// contig end; one starting before them is kept from there, unless the previous piece holds it.
fn owned_part(piece: &Piece, start: u64, end: u64) -> Option<(u64, u64)> {
    let to_contig = |position: u64| position - piece.start + piece.offset;
    let to_batch = |position: u64| position - piece.offset + piece.start;
    let (start, end) = (to_contig(start), to_contig(end));
    let (owned_start, owned_end) = piece.owned;
    let piece_end = piece.offset + piece.length;
    let kept_end = match end == piece_end && owned_end < piece_end {
        true => owned_end,
        false => end,
    };
    let kept_start = match start {
        _ if start >= owned_end => return None,
        _ if start >= owned_start => start,
        _ if end < piece.previous_end => return None,
        _ => owned_start,
    };
    (kept_start < kept_end).then(|| (to_batch(kept_start), to_batch(kept_end)))
}
