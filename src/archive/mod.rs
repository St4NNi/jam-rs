//! Backend-neutral random-access contracts for trace archives.

use crate::resource::ResourceMetrics;
use serde::{Deserialize, Serialize};
use thiserror::Error;

pub mod native;
pub mod pithos;
pub use native::NativeJmaArchive;
pub use pithos::{
    PithosArchiveOptions, PithosBiosequenceArchive, PithosBiosequenceSource, PithosProfile,
    PithosSequenceOrganization,
};

pub const JAMHASH_ALGORITHM_ID: &str = "jamhash_u64_v1";

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct SeedSchemeId(pub u32);

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SeedSchemeDescriptor {
    pub scheme_id: u32,
    pub algorithm_id: u32,
    pub span: u16,
    pub informative_bases: u16,
    pub density_parameter: u32,
    pub bucket_bits: u8,
    pub key_encoding: u8,
    pub occurrence_encoding: u8,
    pub flags: u32,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct ArchiveContig {
    pub id: u32,
    pub name: String,
    pub length: u64,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct ArchiveMetadata {
    pub format_identifier: String,
    pub format_version: u16,
    pub layout_identifier: u32,
    pub source_assembly_sha256: [u8; 32],
    pub archive_sha256: Option<[u8; 32]>,
    pub builder_version: String,
    pub source_commit: String,
    pub hash_algorithm: String,
    pub total_bases: u64,
    pub contigs: Vec<ArchiveContig>,
}

/// A digest lookup is only a candidate. `verification` contains the exact
/// packed selected bases required before an occurrence can become an anchor.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct SeedKey {
    pub digest: u64,
    pub verification: Vec<u8>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedOccurrence {
    pub contig_id: u32,
    pub position: u64,
    pub span: u16,
    pub reverse: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SeedMatch {
    pub key_index: usize,
    pub occurrences: Vec<SeedOccurrence>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct SeedLookupMetrics {
    pub pages_read: u64,
    pub keys_tested: u64,
    pub occurrences_before_limits: u64,
    pub occurrences_after_limits: u64,
    pub bytes_read: u64,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct SeedLookupResult {
    pub matches: Vec<SeedMatch>,
    pub metrics: SeedLookupMetrics,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SequenceRequest {
    pub contig_id: u32,
    pub start: u64,
    pub end: u64,
    pub reverse_complement: bool,
}

impl SequenceRequest {
    pub fn new(
        contig_id: u32,
        start: u64,
        end: u64,
        reverse_complement: bool,
    ) -> ArchiveResult<Self> {
        if start > end {
            return Err(ArchiveError::InvalidSequenceRange { start, end });
        }
        Ok(Self {
            contig_id,
            start,
            end,
            reverse_complement,
        })
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

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SequenceSlice {
    pub request: SequenceRequest,
    pub bases: Vec<u8>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct ArchiveMetrics {
    pub resource: ResourceMetrics,
    pub mapped_bytes: u64,
    pub resident_bytes: u64,
    pub metadata_bytes_read: u64,
    pub seed_bytes_read: u64,
    pub sequence_bytes_read: u64,
    pub decoded_sequence_bases: u64,
    pub coalesced_range_requests: u64,
}

/// Synchronous batch operations match the current bounded range-reader and
/// CPU-pool execution model. Implementations must not eagerly read the object.
pub trait TraceArchive: Send + Sync {
    fn metadata(&self) -> ArchiveResult<ArchiveMetadata>;

    fn available_seed_schemes(&self) -> ArchiveResult<Vec<SeedSchemeDescriptor>>;

    fn lookup_seeds(
        &self,
        scheme: SeedSchemeId,
        keys: &[SeedKey],
    ) -> ArchiveResult<SeedLookupResult>;

    fn read_sequences(&self, requests: &[SequenceRequest]) -> ArchiveResult<Vec<SequenceSlice>>;

    fn metrics(&self) -> ArchiveMetrics;
}

pub type ArchiveResult<T> = Result<T, ArchiveError>;

#[derive(Debug, Error)]
pub enum ArchiveError {
    #[error("invalid archive sequence range: start={start}, end={end}")]
    InvalidSequenceRange { start: u64, end: u64 },
    #[error("unknown archive contig id {0}")]
    UnknownContig(u32),
    #[error("unknown archive seed scheme {0}")]
    UnknownSeedScheme(u32),
    #[error("archive metadata is corrupt: {0}")]
    CorruptMetadata(String),
    #[error("archive checksum mismatch in {0}")]
    ChecksumMismatch(String),
    #[error("archive backend failed: {0}")]
    Backend(String),
    #[error(transparent)]
    Resource(#[from] crate::resource::ResourceError),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sequence_requests_are_zero_based_half_open() {
        let request = SequenceRequest::new(2, 4, 11, false).unwrap();
        assert_eq!(request.len(), 7);
        assert!(SequenceRequest::new(2, 11, 4, false).is_err());
    }

    #[test]
    fn seed_keys_include_exact_verification_material() {
        let left = SeedKey {
            digest: 7,
            verification: vec![1, 2],
        };
        let right = SeedKey {
            digest: 7,
            verification: vec![1, 3],
        };
        assert_ne!(left, right);
    }
}
