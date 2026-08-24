//! Versioned random-access metagenome assembly archives.
//!
//! JMA is separate from the existing `.jam` sketch database. `.jam` remains v3.

use crate::resource::RangeReader;
use serde::{Deserialize, Serialize};
use thiserror::Error;

pub mod builder;
pub mod contigs;
pub mod format;
pub mod header;
pub mod index;
pub mod reader;
pub mod seed_builder;
pub mod sequence;
pub mod sequence_builder;
pub mod writer;

pub const JMA_MAGIC: [u8; 8] = *b"JMAF1\0\0\0";
pub const JMA_FORMAT_VERSION: u16 = 1;
pub const JMA_LAYOUT_IDENTIFIER: u32 = 0x3141_4d4a;
pub type ContigId = u32;

pub use crate::sequence::SequenceBlockRecord;
pub use format::{ArchiveMetadataFields, SectionDescriptor};
pub use index::SeedIndexDirectory;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SequenceRange {
    pub start: u64,
    pub end: u64,
}

impl SequenceRange {
    pub fn new(start: u64, end: u64) -> JmaResult<Self> {
        if start > end {
            return Err(JmaError::InvalidSequenceRange { start, end });
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
pub struct SeedLevel {
    pub k: u8,
    pub scale: u64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ArchiveHeader {
    pub format_version: u16,
    pub flags: u32,
    pub contig_count: u32,
    pub total_bases: u64,
    pub source_sha256: [u8; 32],
    pub seed_levels: Vec<SeedLevel>,
    pub algorithm_id: Option<String>,
    pub algorithm_version: Option<u16>,
    pub min_entropy: Option<f64>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct ContigMetadata {
    pub id: ContigId,
    pub name: String,
    pub length: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SeedQuery {
    pub k: u8,
    pub hash: u64,
    pub canonical_kmer: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SeedOccurrence {
    pub contig_id: ContigId,
    pub position: u64,
    pub reverse: bool,
}

/// Random-access contract shared by local and remote JMA readers.
pub trait ArchiveReader: Send + Sync {
    fn resource(&self) -> &dyn RangeReader;
    fn header(&self) -> &ArchiveHeader;
    fn contigs(&self) -> &[ContigMetadata];
    fn read_sequence(&self, contig_id: ContigId, range: SequenceRange) -> JmaResult<Vec<u8>>;
    fn seed_occurrences(&self, query: SeedQuery) -> JmaResult<Vec<SeedOccurrence>>;
    fn metrics(&self) -> crate::resource::ResourceMetrics {
        self.resource().metrics()
    }
}

pub type JmaResult<T> = Result<T, JmaError>;

#[derive(Debug, Error)]
pub enum JmaError {
    #[error("invalid JMA magic")]
    InvalidMagic,
    #[error("unsupported JMA format version {0}")]
    UnsupportedVersion(u16),
    #[error("invalid JMA sequence range: start={start}, end={end}")]
    InvalidSequenceRange { start: u64, end: u64 },
    #[error("unknown JMA contig id {0}")]
    UnknownContig(ContigId),
    #[error("JMA offset or length overflow")]
    OffsetOverflow,
    #[error("JMA section is corrupt: {0}")]
    CorruptSection(String),
    #[error("JMA checksum mismatch in {0}")]
    ChecksumMismatch(String),
    #[error("JMA resource error: {0}")]
    Resource(#[from] crate::resource::ResourceError),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jma_and_jam_versions_are_independent() {
        assert_eq!(JMA_FORMAT_VERSION, 1);
        assert_eq!(crate::format::VERSION, 3);
    }

    #[test]
    fn sequence_ranges_are_half_open() {
        let range = SequenceRange::new(4, 11).unwrap();
        assert_eq!(range.len(), 7);
        assert!(SequenceRange::new(2, 1).is_err());
    }
}
