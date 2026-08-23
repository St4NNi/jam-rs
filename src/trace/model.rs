//! Shared coordinate, alignment, coverage, and JSON record models.

use super::config::{SensitivityConfig, TraceAlgorithmMetadata};
use crate::resource::ResourceMetrics;
use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct BaseInterval {
    pub start: u64,
    pub end: u64,
}

impl BaseInterval {
    pub fn new(start: u64, end: u64) -> Result<Self, CoordinateError> {
        if start > end {
            return Err(CoordinateError::Reversed { start, end });
        }
        Ok(Self { start, end })
    }

    #[must_use]
    pub fn len(self) -> u64 {
        self.end - self.start
    }

    #[must_use]
    pub fn is_empty(self) -> bool {
        self.start == self.end
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EditOperation {
    Equal,
    Substitution,
    Insertion,
    Deletion,
    SoftClip,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct EditRun {
    pub operation: EditOperation,
    pub length: u32,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BaseAlignment {
    pub plasmid_id: String,
    pub metagenome_id: String,
    pub contig_id: String,
    pub strand: Strand,
    /// One segment normally, two ordered segments when the circular query crosses origin.
    pub query_segments: Vec<BaseInterval>,
    pub target_interval: BaseInterval,
    pub query_length: u64,
    pub target_length: u64,
    pub origin_crossing: bool,
    pub score: i64,
    pub matches: u64,
    pub substitutions: u64,
    pub insertions: u64,
    pub deletions: u64,
    pub cigar: String,
    pub edit_script: Vec<EditRun>,
    pub chain_score: i64,
    pub primary: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CandidateResult {
    pub metagenome_id: String,
    pub shared_hashes: u64,
    pub plasmid_hashes: u64,
    pub metagenome_hashes: u64,
    pub plasmid_containment: f64,
    pub metagenome_containment: f64,
    pub rank: u32,
    pub score_mode: String,
    pub bias_weighted_plasmid_containment: Option<f64>,
    pub uniform_hash_e_value: Option<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GapRecord {
    pub interval: BaseInterval,
    pub length: u64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoverageSummary {
    pub plasmid_length: u64,
    pub supported_bases: u64,
    pub supported_fraction: f64,
    pub primary_intervals: Vec<BaseInterval>,
    pub secondary_intervals: Vec<BaseInterval>,
    pub gaps: Vec<GapRecord>,
    pub largest_gap: u64,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct InputResource {
    pub role: String,
    pub redacted_locator: String,
    pub sha256: Option<String>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct TraceFailure {
    pub stage: String,
    pub code: String,
    pub message: String,
    pub resource: Option<String>,
    pub retryable: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TraceStatus {
    Complete,
    Partial,
    NoCandidate,
    Failed,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TraceRunHeader {
    pub schema_version: String,
    pub run_id: String,
    pub jam_rs_version: String,
    pub source_commit: Option<String>,
    pub started_at_utc: String,
    pub command: Vec<String>,
    pub plasmid_id: String,
    pub plasmid_length: u64,
    pub sensitivity: SensitivityConfig,
    pub algorithm: TraceAlgorithmMetadata,
    pub inputs: Vec<InputResource>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TraceMetagenomeResult {
    pub schema_version: String,
    pub run_id: String,
    pub plasmid_id: String,
    pub metagenome_id: String,
    pub algorithm: TraceAlgorithmMetadata,
    pub status: TraceStatus,
    pub candidate: Option<CandidateResult>,
    pub alignments: Vec<BaseAlignment>,
    pub coverage: Option<CoverageSummary>,
    pub failures: Vec<TraceFailure>,
    pub resource_metrics: ResourceMetrics,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct TraceRunFooter {
    pub schema_version: String,
    pub run_id: String,
    pub completed_at_utc: String,
    pub metagenomes_total: u64,
    pub metagenomes_with_candidates: u64,
    pub metagenomes_aligned: u64,
    pub metagenomes_failed: u64,
    pub alignments_total: u64,
    pub resource_metrics: ResourceMetrics,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "record_type", rename_all = "snake_case")]
pub enum TraceRecord {
    RunHeader(TraceRunHeader),
    MetagenomeResult(TraceMetagenomeResult),
    RunFooter(TraceRunFooter),
}

#[derive(Debug, Error)]
pub enum CoordinateError {
    #[error("base interval is reversed: start={start}, end={end}")]
    Reversed { start: u64, end: u64 },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intervals_are_zero_based_half_open() {
        let interval = BaseInterval::new(3, 8).unwrap();
        assert_eq!(interval.len(), 5);
        assert!(BaseInterval::new(8, 3).is_err());
    }

    #[test]
    fn trace_records_have_stable_record_type() {
        let record = TraceRecord::RunFooter(TraceRunFooter {
            schema_version: super::super::TRACE_JSON_SCHEMA_VERSION.to_string(),
            run_id: "run".to_string(),
            completed_at_utc: "1970-01-01T00:00:00Z".to_string(),
            metagenomes_total: 0,
            metagenomes_with_candidates: 0,
            metagenomes_aligned: 0,
            metagenomes_failed: 0,
            alignments_total: 0,
            resource_metrics: ResourceMetrics::default(),
        });
        let value = serde_json::to_value(record).unwrap();
        assert_eq!(value["record_type"], "run_footer");
    }
}
