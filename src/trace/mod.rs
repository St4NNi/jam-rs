//! One-query-to-many-metagenomes fragment-search contracts.

pub mod alignment;
pub mod anchors;
pub mod batch;
pub mod catalog;
pub mod chain;
pub mod config;
pub mod contracts;
pub mod coverage;
pub mod intervals;
pub mod memory;
pub mod model;
pub mod mosaic;
pub mod output;
pub mod query_plan;
pub mod raw;
pub mod runner;
pub mod scheduler;
pub mod screen;
pub mod seeds;

pub use config::{
    AlignmentScoring, GapRescueConfig, SeedSensitivity, SensitivityConfig, SensitivityProfile,
    TraceAlgorithmMetadata, TraceAlgorithmParameters, algorithm_identifiers,
};
pub use contracts::{
    BatchCandidateWork, BatchResultKey, CandidateMemoryBudget, MemoryBudgetError,
    PreparedBatchQuery, TraceMemoryEstimate, WitnessPlanningRequest, WitnessTierPlan,
};
pub use model::{
    AlignmentMosaicEvidence, AlignmentRole, BaseAlignment, BaseInterval,
    CandidatePerformanceCounters, CandidateResult, CoordinateModel, CoverageSummary, EditOperation,
    EditRun, FragmentMosaicSummary, MosaicAtomicInterval, MosaicSelectionComponents,
    QueryDescriptor, QueryKind, RescueRoundMetrics, SeedEvidence, SeedSchemeDiagnostics, Strand,
    TopologyAssessment, TopologyEvidence, TopologyModelEvidence, TopologyModelSummary,
    TopologyRequested, TraceAlgorithmIdentifiers, TraceDiagnosticReport, TraceFailureCategory,
    TraceMetagenomeResult, TraceRecord, TraceRunFooter, TraceRunHeader, TraceStageMetrics,
    TruthIntervalDiagnostics,
};

pub const TRACE_ALGORITHM_ID: &str = "jam-seed-chain-align-v1";
pub const TRACE_ALGORITHM_VERSION: u16 = 1;
pub const TRACE_JSON_SCHEMA_VERSION: &str = "2.0.0";
