//! The output records of a screen, written as JSON lines.

use serde::Serialize;

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
