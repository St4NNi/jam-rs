//! Shared coordinate, alignment, coverage, and JSON record models.

use super::config::{SensitivityConfig, TraceAlgorithmMetadata};
use crate::resource::ResourceMetrics;
use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum QueryKind {
    Plasmid,
    Phage,
    Other,
    #[default]
    Unknown,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TopologyRequested {
    Linear,
    Circular,
    #[default]
    Auto,
    Unknown,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CoordinateModel {
    Linear,
    Wrap,
    Undetermined,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TopologyEvidence {
    LinearSupported,
    WrapSupported,
    BothCompatible,
    Insufficient,
    Undetermined,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AlignmentRole {
    PrimaryMosaic,
    OverlappingSupport,
    AlternativeMapping,
    CommonSequence,
    RepeatOnly,
    OriginCrossing,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct BaseInterval {
    pub start: u64,
    pub end: u64,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct QueryDescriptor {
    pub query_id: String,
    pub query_length: u64,
    pub query_kind: QueryKind,
    pub topology_requested: TopologyRequested,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct TraceAlgorithmIdentifiers {
    pub screen_algorithm: String,
    pub local_alignment_algorithm: String,
    pub mosaic_algorithm: String,
    pub trace_workflow: String,
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct SeedEvidence {
    pub primary_anchor_count: u64,
    pub rescue_anchor_count: u64,
    pub nonrepetitive_anchor_count: u64,
    pub common_anchor_count: u64,
    pub repetitive_seed_count: u64,
    pub query_occurrence_max: u32,
    pub candidate_occurrence_max: u32,
    pub collection_document_frequency_max: Option<u32>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct AlignmentMosaicEvidence {
    pub alignment_id: String,
    pub primary_supported_bases: u64,
    pub secondary_supported_bases: u64,
    pub newly_supported_bases: u64,
    pub role: AlignmentRole,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct MosaicAtomicInterval {
    pub interval: BaseInterval,
    pub primary_alignment_id: Option<String>,
    pub alternative_alignment_ids: Vec<String>,
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct MosaicSelectionComponents {
    pub local_alignment_score_sum: i64,
    pub newly_supported_query_bases: u64,
    pub nonrepetitive_anchor_evidence: u64,
    pub redundant_overlap_bases: u64,
    pub coordinate_contradictions: u64,
    pub sequence_contradictions: u64,
    pub fragment_count: u64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FragmentMosaicSummary {
    pub mosaic_algorithm: String,
    pub coordinate_model: CoordinateModel,
    pub base_covered_bases: u64,
    pub base_coverage_fraction: f64,
    pub aligned_span_bases: u64,
    pub aligned_span_fraction: f64,
    pub covered_intervals: Vec<BaseInterval>,
    pub unsupported_gaps: Vec<GapRecord>,
    pub alignment_deletions: Vec<BaseInterval>,
    pub common_sequence_intervals: Vec<BaseInterval>,
    pub repeat_only_intervals: Vec<BaseInterval>,
    pub supporting_contigs: Vec<String>,
    pub accepted_alignment_count: u64,
    pub alternative_alignment_count: u64,
    pub nonrepetitive_supported_bases: u64,
    pub common_sequence_supported_bases: u64,
    pub repeat_only_supported_bases: u64,
    pub atomic_intervals: Vec<MosaicAtomicInterval>,
    pub alignment_evidence: Vec<AlignmentMosaicEvidence>,
    pub selection_components: MosaicSelectionComponents,
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct TopologyModelEvidence {
    pub newly_supported_query_bases: u64,
    pub alignment_quality_sum: i64,
    pub nonrepetitive_anchor_support: u64,
    pub origin_crossing_alignment_count: u64,
    pub unexplained_terminal_gap_bases: u64,
    pub contradictory_alignment_count: u64,
    pub fragment_count: u64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TopologyModelSummary {
    pub coordinate_model: CoordinateModel,
    pub mosaic: FragmentMosaicSummary,
    pub evidence: TopologyModelEvidence,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TopologyAssessment {
    pub topology_requested: TopologyRequested,
    pub coordinate_model: CoordinateModel,
    pub topology_evidence: TopologyEvidence,
    pub selection_margin_bases: u64,
    pub linear_model: TopologyModelSummary,
    pub wrap_model: Option<TopologyModelSummary>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct RescueRoundMetrics {
    pub round: u8,
    pub seed_k: u8,
    pub seed_scale: u64,
    pub target_gaps: Vec<BaseInterval>,
    pub seed_buckets_requested: u64,
    pub seed_keys_tested: u64,
    pub anchors_created: u64,
    pub chains_accepted: u64,
    pub sequence_blocks_fetched: u64,
    pub alignment_windows_attempted: u64,
    pub new_query_bases_supported: u64,
    pub elapsed_millis: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TraceFailureCategory {
    CandidateMiss,
    NoRetainedSeed,
    NoMatchingSeed,
    RepetitiveSeedSuppression,
    AnchorCap,
    NoValidChain,
    ChainLimit,
    SequenceBudget,
    AlignmentBandFailure,
    AlignmentRejection,
    AlignmentCap,
    MosaicSelection,
    NoRetainedRouterWitness,
    NoMatchingPositionalWitness,
    OccurrenceCap,
    NoAcceptedChain,
    ChainBridgeAbsent,
    SequenceWindowAbsent,
    CorridorEdge,
    LocalClipping,
    RefinementFailure,
    OutputCap,
    Other,
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct SeedSchemeDiagnostics {
    pub scheme_id: u32,
    pub query_seeds_emitted: u64,
    pub matching_target_keys: u64,
    pub matching_occurrences_before_limits: u64,
    pub matching_occurrences_after_limits: u64,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct TruthIntervalDiagnostics {
    pub truth_interval: BaseInterval,
    pub candidate_selected: bool,
    pub seed_schemes: Vec<SeedSchemeDiagnostics>,
    pub anchors_before_limits: u64,
    pub anchors_after_limits: u64,
    pub best_chain_anchor_count: u64,
    pub best_chain_query_span: u64,
    pub best_chain_target_span: u64,
    pub sequence_blocks_requested: u64,
    pub sequence_window_decoded: Option<BaseInterval>,
    pub alignment_attempted: bool,
    pub alignment_band_widths: Vec<u32>,
    pub alignment_band_edge_touched: bool,
    pub alignment_accepted: bool,
    pub alignment_query_span: u64,
    pub mosaic_supported_bases: u64,
    pub alignment_or_output_cap_reached: bool,
    pub failure_category: Option<TraceFailureCategory>,
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct TraceStageMetrics {
    pub stage: u8,
    pub name: String,
    pub seed_pages_read: u64,
    pub keys_tested: u64,
    pub occurrences_decoded: u64,
    pub anchors_retained: u64,
    pub chains_retained: u64,
    pub sequence_blocks_read: u64,
    pub alignment_attempts: u64,
    pub new_query_bases_supported: u64,
    pub wall_micros: u64,
    pub cpu_micros: u64,
    pub bytes_read: u64,
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct TraceDiagnosticReport {
    pub metagenome_id: String,
    pub truth_intervals: Vec<TruthIntervalDiagnostics>,
    pub stages: Vec<TraceStageMetrics>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct CandidatePerformanceCounters {
    pub candidates_processed: u64,
    pub contigs_considered: u64,
    pub windows_retrieved: u64,
    pub alignments_attempted: u64,
    pub alignments_succeeded: u64,
    pub failures: u64,
    pub elapsed_millis: u64,
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
    pub alignment_id: String,
    #[serde(rename = "query_id", alias = "plasmid_id")]
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
    pub identity: f64,
    pub seed_evidence: SeedEvidence,
    pub primary_supported_bases: u64,
    pub secondary_supported_bases: u64,
    pub newly_supported_bases: u64,
    pub role: AlignmentRole,
    pub primary: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CandidateResult {
    pub metagenome_id: String,
    pub shared_hashes: u64,
    #[serde(rename = "query_hashes", alias = "plasmid_hashes")]
    pub plasmid_hashes: u64,
    pub metagenome_hashes: u64,
    #[serde(default)]
    pub shared_spatial_signatures: u64,
    #[serde(default)]
    pub rare_shared_signatures: u64,
    #[serde(default)]
    pub shared_whole_sample_signatures: u64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub screen_policy: Option<String>,
    #[serde(
        default,
        rename = "candidate_entry_reason",
        skip_serializing_if = "Option::is_none"
    )]
    pub admission_source: Option<CandidateAdmissionSource>,
    #[serde(rename = "query_containment", alias = "plasmid_containment")]
    pub plasmid_containment: f64,
    pub metagenome_containment: f64,
    pub rank: u32,
    pub score_mode: String,
    pub bias_weighted_plasmid_containment: Option<f64>,
    pub uniform_hash_e_value: Option<f64>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CandidateAdmissionSource {
    Standard,
    WindowSpread,
    RareRescue,
    WholeSampleFallback,
    Hamming1Rescue,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GapRecord {
    /// First ordinary segment retained for one-release compatibility.
    pub interval: BaseInterval,
    /// One segment for a linear gap; two terminal segments for a gap wrapping origin.
    pub segments: Vec<BaseInterval>,
    pub wraps_origin: bool,
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
    #[serde(rename = "query_id", alias = "plasmid_id")]
    pub plasmid_id: String,
    #[serde(rename = "query_length", alias = "plasmid_length")]
    pub plasmid_length: u64,
    pub query_kind: QueryKind,
    pub topology_requested: TopologyRequested,
    pub threads: usize,
    pub io_concurrency: usize,
    pub sensitivity: SensitivityConfig,
    pub algorithms: TraceAlgorithmIdentifiers,
    pub algorithm: TraceAlgorithmMetadata,
    pub inputs: Vec<InputResource>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TraceMetagenomeResult {
    pub schema_version: String,
    pub run_id: String,
    #[serde(rename = "query_id", alias = "plasmid_id")]
    pub plasmid_id: String,
    pub metagenome_id: String,
    pub query_kind: QueryKind,
    pub topology_requested: TopologyRequested,
    pub coordinate_model: CoordinateModel,
    pub topology_evidence: TopologyEvidence,
    pub algorithms: TraceAlgorithmIdentifiers,
    pub algorithm: TraceAlgorithmMetadata,
    pub status: TraceStatus,
    pub candidate: Option<CandidateResult>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub router_candidate: Option<crate::router::RoutedCandidate>,
    pub alignments: Vec<BaseAlignment>,
    pub primary_fragment_mosaic: Option<FragmentMosaicSummary>,
    pub topology: Option<TopologyAssessment>,
    pub rescue_rounds: Vec<RescueRoundMetrics>,
    pub stages: Vec<TraceStageMetrics>,
    pub alignment_retries: Vec<super::alignment::AlignmentRetryMetadata>,
    pub performance_counters: CandidatePerformanceCounters,
    pub coverage: Option<CoverageSummary>,
    pub warnings: Vec<String>,
    pub failures: Vec<TraceFailure>,
    pub archive_metrics: Option<crate::archive::ArchiveMetrics>,
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
#[allow(clippy::large_enum_variant)]
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

    #[test]
    fn generic_query_and_topology_values_have_stable_names() {
        assert_eq!(
            serde_json::to_string(&QueryKind::Phage).unwrap(),
            "\"phage\""
        );
        assert_eq!(
            serde_json::to_string(&TopologyRequested::Unknown).unwrap(),
            "\"unknown\""
        );
        assert_eq!(
            serde_json::to_string(&CoordinateModel::Wrap).unwrap(),
            "\"wrap\""
        );
        assert_eq!(
            serde_json::to_string(&TopologyEvidence::BothCompatible).unwrap(),
            "\"both_compatible\""
        );
    }
}
