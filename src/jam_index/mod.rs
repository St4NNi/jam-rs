//! Local append-only Jam Index datasets.
//!
//! A Jam Index is a manifest plus independently buildable/searchable parts.
//! Every part owns a pure `.jam` screening shard, compact contig postings,
//! metadata, and checked references to external assembly files. It stores no
//! nucleotide sequence or nucleotide seed positions.

pub mod archive;
pub mod batch;
pub mod builder;
pub mod contig_search;
pub mod distributed;
pub(crate) mod external;
pub mod manifest;
pub mod part;
pub mod screen;
pub mod signature;
pub mod trace;

pub use archive::JamIndexArchive;
pub use batch::{
    JamIndexBatchConfig, JamIndexBatchError, JamIndexBatchExecution, JamIndexBatchMetrics,
    JamIndexBatchPlan, JamIndexBatchQuery, JamIndexBatchQueryStatus, JamIndexBatchStatusKind,
    JamIndexBatchWork, JamIndexBatchWorkGroup, JamIndexPlannedQuery, execute_batch,
    execute_batch_serialized, plan_batch,
};
pub use builder::{
    JamIndexBuildConfig, JamIndexBuildError, JamIndexBuildStats, append_jam_index, build_jam_index,
    load_manifest,
};
pub use contig_search::{
    JamIndexBatchContigSearchResult, JamIndexContigPlan, JamIndexContigSearchConfig,
    JamIndexContigSearchError, JamIndexContigSearchMetrics, JamIndexContigSearchResult,
    RankedJamIndexContig, select_candidate_contigs, select_candidate_contigs_batch,
};
pub use distributed::{
    DISTRIBUTED_PLAN_SCHEMA, DistributedIndexError, FRAGMENT_MANIFEST_SCHEMA,
    IndexBuildFragmentPlan, IndexBuildPartPlan, IndexBuildPlan, IndexFragmentManifest,
    IndexFragmentSourceResult, IndexPlanSource, MERGED_PART_MANIFEST_SCHEMA, MergedPartManifest,
    build_fragment, finalize_index, load_plan, merge_part, plan_index, write_plan_atomic,
};
pub use manifest::{
    ContigSignatureBudget, JamIndexManifest, JamIndexManifestError, JamIndexPart,
    ScreenSelectionPolicy,
};
pub use part::{
    JamIndexPartError, JamIndexPartReader, LoadedPartContig, MetagenomeSource, PartMetagenome,
    PartReadResult, PartScreenSample, PartWriteResult, PublishedMetagenomeSource,
    StagedMetagenomeSource, merge_part_fragments, write_part, write_part_staged,
};
pub use screen::{
    JamIndexBatchScreenResult, JamIndexCandidate, JamIndexScreenConfig, JamIndexScreenError,
    JamIndexScreenMetrics, JamIndexScreenResult, PreparedJamIndexQuery, QueryHashOccurrence,
    SharedScreenHash, prepare_screen_query, search_jam_index, search_jam_index_batch,
};
pub use signature::{
    ContigSignature, MetagenomeSignatureBuilder, SignatureSelectionError, SpatialSignature,
};
pub use trace::{
    JamIndexTraceConfig, JamIndexTraceError, JamIndexTraceMetrics, JamIndexTraceOutput, trace_index,
};
