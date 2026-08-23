//! One-plasmid-to-many-metagenomes trace-search contracts.

pub mod alignment;
pub mod anchors;
pub mod catalog;
pub mod chain;
pub mod config;
pub mod coverage;
pub mod intervals;
pub mod model;
pub mod mosaic;
pub mod output;
pub mod raw;
pub mod runner;
pub mod screen;
pub mod seeds;

pub use config::{AlignmentScoring, SeedSensitivity, SensitivityConfig, SensitivityProfile};
pub use model::{
    BaseAlignment, BaseInterval, CandidateResult, CoverageSummary, EditOperation, EditRun, Strand,
    TraceMetagenomeResult, TraceRecord, TraceRunFooter, TraceRunHeader,
};

pub const TRACE_JSON_SCHEMA_VERSION: &str = "1.0.0";
