//! Local append-only Jam Index datasets.
//!
//! A Jam Index is a manifest plus independently buildable/searchable parts.
//! Every part owns a pure `.jam` screening shard, compact contig signatures,
//! and complete two-bit contig sequence. It never stores nucleotide seed
//! positions.

pub mod manifest;
pub mod signature;

pub use manifest::{
    ContigSignatureBudget, JamIndexManifest, JamIndexManifestError, JamIndexPart,
    ScreenSelectionPolicy,
};
pub use signature::{ContigSignature, MetagenomeSignatureBuilder, SignatureSelectionError};
