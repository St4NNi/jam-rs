//! Local append-only Jam Index datasets.
//!
//! A Jam Index is a manifest plus independently buildable/searchable parts.
//! Every part owns a pure `.jam` screening shard, compact contig postings,
//! metadata, and checked references to external assembly files. It stores no
//! nucleotide sequence or nucleotide seed positions.

pub mod archive;
pub mod builder;
pub mod contig_search;
pub(crate) mod external;
pub mod manifest;
pub mod part;
pub mod screen;
pub mod signature;
pub mod trace;

pub use archive::JamIndexArchive;
pub use builder::{
    JamIndexBuildConfig, JamIndexBuildError, JamIndexBuildStats, append_jam_index, build_jam_index,
    load_manifest,
};
pub use contig_search::{
    JamIndexContigPlan, JamIndexContigSearchConfig, JamIndexContigSearchError,
    JamIndexContigSearchMetrics, JamIndexContigSearchResult, RankedJamIndexContig,
    select_candidate_contigs,
};
pub use manifest::{
    ContigSignatureBudget, JamIndexManifest, JamIndexManifestError, JamIndexPart,
    ScreenSelectionPolicy,
};
pub use part::{
    JamIndexPartReader, MetagenomeSource, PartContig, PartMetagenome, PartScreenSample,
    PartWriteResult, write_part,
};
pub use screen::{
    JamIndexCandidate, JamIndexScreenConfig, JamIndexScreenError, JamIndexScreenMetrics,
    JamIndexScreenResult, PreparedJamIndexQuery, QueryHashOccurrence, SharedScreenHash,
    prepare_screen_query, search_jam_index,
};
pub use signature::{ContigSignature, MetagenomeSignatureBuilder, SignatureSelectionError};
pub use trace::{
    JamIndexTraceConfig, JamIndexTraceError, JamIndexTraceMetrics, JamIndexTraceOutput, trace_index,
};
