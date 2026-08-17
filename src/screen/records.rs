//! The output records of a screen: one entry record for each entry contig found in an assembly,
//! and one summary record for each assembly. Both are written as JSON lines.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::align::Interval;
use crate::decide::{Fragment, SelectedFragment};
use serde::Serialize;

/// An entry contig (the query of the screen) found in an assembly.
#[derive(Clone, Debug, Serialize)]
pub struct ScreenEntry {
    /// Always `entry`.
    #[serde(rename = "type")]
    pub kind: &'static str,
    /// The entry name, followed by `/` and the contig name when the entry has more than one
    /// contig.
    pub entry: String,
    /// Bases of the entry contig.
    pub entry_length: u64,
    pub assembly: String,
    /// Bases of the entry contig covered by the mosaic, and their share of it.
    pub covered_bases: u64,
    pub query_coverage: f64,
    /// Matches over alignment columns of the counted fragments together.
    pub identity: f64,
    pub matches: u64,
    pub substitutions: u64,
    pub insertions: u64,
    pub deletions: u64,
    /// The assembly contigs that the fragments of the mosaic lie on, by id.
    pub contigs: Vec<ScreenContig>,
    /// The fragments of the mosaic: query intervals on the entry contig, contig intervals on the
    /// assembly contig named by `contig_index` (its id), CIGAR with the assembly contig as the
    /// reference.
    pub primary: Vec<SelectedFragment>,
    pub alternatives: Vec<Fragment>,
    /// The stretches of the entry contig that no fragment covers.
    pub gaps: Vec<Interval>,
}

/// An assembly contig that a reported entry contig lies on.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ScreenContig {
    /// The id of the contig: its rank among the searched contigs of the assembly.
    pub id: u32,
    /// The first word of the FASTA header.
    pub name: String,
    pub length: u64,
    /// Bases of the contig in the fragments of the mosaic, counted once.
    pub aligned: u64,
    /// Matches over columns of all fragments of the contig with the entry contig.
    pub identity: f64,
    /// How close the entry contig is to the contig: the aligned share of the contig times the
    /// identity. `best_score` is the highest score any entry contig reaches for this contig in the
    /// assembly, so a loose screen can be reduced to the closest entries later.
    pub score: f64,
    pub best_score: f64,
    /// The entry contigs that this contig counts for in the assembly, this one included. A contig
    /// that counts for many (an insertion sequence, a common backbone) says little about any one.
    pub entries: u32,
}

/// The summary of one screened assembly.
#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct AssemblySummary {
    /// Always `assembly`.
    #[serde(rename = "type")]
    pub kind: &'static str,
    pub assembly: String,
    /// Searched contigs and their bases, and contigs skipped as too short or unknown.
    pub contigs: u64,
    pub bases: u64,
    pub skipped_contigs: u64,
    pub batches: u64,
    /// Entries that enough contigs share enough seeds with to be aligned.
    pub surveyed_entries: u64,
    /// Entries aligned in a batch, counted once per batch.
    pub aligned_entries: u64,
    /// Fragments after the cut at contig borders.
    pub fragments: u64,
    pub alignments: u64,
    pub reported_entries: u64,
    /// Contigs that share enough seeds with a kept entry, and their bases.
    pub matched_contigs: u64,
    pub matched_bases: u64,
    /// Contigs that reported entry contigs lie on, and their bases.
    pub reported_contigs: u64,
    pub reported_bases: u64,
    /// Wall time of the survey and of the alignments.
    pub survey_seconds: f64,
    pub align_seconds: f64,
}
