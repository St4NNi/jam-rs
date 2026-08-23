use jam_rs::jma::reader::JmaReader;
use jam_rs::jma::sequence::{SequenceBlock, pack_bases};
use jam_rs::jma::writer::{ArchiveParts, SeedLevelData, SeedRecord, SeedSection, encode_archive};
use jam_rs::jma::{
    ArchiveReader, ContigMetadata, JmaResult, SeedLevel, SeedOccurrence, SeedQuery, SequenceRange,
};
use jam_rs::resource::{
    ByteRange, RangeReader, ResourceError, ResourceLocator, ResourceMetadata, ResourceMetrics,
};
use std::io::{Cursor, Read};
use std::sync::Arc;

struct MemoryResource {
    locator: ResourceLocator,
    bytes: Arc<Vec<u8>>,
}

impl MemoryResource {
    fn new(bytes: Vec<u8>) -> Self {
        Self {
            locator: ResourceLocator::parse("memory://jma-fixture").unwrap_or_else(|_| {
                ResourceLocator::parse("file://jma-fixture").expect("file locator")
            }),
            bytes: Arc::new(bytes),
        }
    }
}

impl RangeReader for MemoryResource {
    fn locator(&self) -> &ResourceLocator {
        &self.locator
    }

    fn metadata(&self) -> Result<ResourceMetadata, ResourceError> {
        Ok(ResourceMetadata {
            size: self.bytes.len() as u64,
            etag: Some("fixture".to_string()),
            last_modified: None,
            accepts_ranges: true,
        })
    }

    fn read_range(&self, range: ByteRange) -> Result<Vec<u8>, ResourceError> {
        let end = range.end()?;
        if end > self.bytes.len() as u64 {
            return Err(ResourceError::RangeOutOfBounds {
                offset: range.offset,
                length: range.length,
                size: self.bytes.len() as u64,
            });
        }
        Ok(self.bytes[range.offset as usize..end as usize].to_vec())
    }

    fn stream(&self) -> Result<Box<dyn Read + Send>, ResourceError> {
        Ok(Box::new(Cursor::new(self.bytes.as_ref().clone())))
    }

    fn metrics(&self) -> ResourceMetrics {
        ResourceMetrics::default()
    }
}

fn fixture() -> ArchiveParts {
    let first = b"ACGTNNACGTAC";
    let second = b"TTTTGCA";
    let (packed_first, mask_first) = pack_bases(first);
    let (packed_second, mask_second) = pack_bases(second);
    ArchiveParts {
        flags: 7,
        source_sha256: [9u8; 32],
        contigs: vec![
            ContigMetadata {
                id: 0,
                name: "contig-a".to_string(),
                length: first.len() as u64,
            },
            ContigMetadata {
                id: 1,
                name: "contig-b".to_string(),
                length: second.len() as u64,
            },
        ],
        sequence_blocks: vec![
            SequenceBlock {
                contig_id: 0,
                start: 0,
                base_length: first.len() as u64,
                packed: packed_first,
                unknown_mask: mask_first,
            },
            SequenceBlock {
                contig_id: 1,
                start: 0,
                base_length: second.len() as u64,
                packed: packed_second,
                unknown_mask: mask_second,
            },
        ],
        seed_sections: vec![
            SeedSection {
                k: 21,
                levels: vec![SeedLevelData {
                    level: SeedLevel { k: 21, scale: 500 },
                    records: vec![SeedRecord {
                        query: SeedQuery {
                            k: 21,
                            hash: 11,
                            canonical_kmer: 3,
                        },
                        occurrence: SeedOccurrence {
                            contig_id: 1,
                            position: 2,
                            reverse: true,
                        },
                    }],
                }],
            },
            SeedSection {
                k: 31,
                levels: vec![SeedLevelData {
                    level: SeedLevel { k: 31, scale: 200 },
                    records: vec![SeedRecord {
                        query: SeedQuery {
                            k: 31,
                            hash: 17,
                            canonical_kmer: 5,
                        },
                        occurrence: SeedOccurrence {
                            contig_id: 0,
                            position: 4,
                            reverse: false,
                        },
                    }],
                }],
            },
        ],
    }
}

#[test]
fn roundtrip_preserves_header_contigs_sequences_and_exact_seeds() -> JmaResult<()> {
    let bytes = encode_archive(&fixture())?;
    assert_eq!(&bytes[0..4], b"JMA\0");
    let reader = JmaReader::from_resource(MemoryResource::new(bytes))?;
    assert_eq!(reader.header().format_version, 1);
    assert_eq!(reader.header().flags, 7);
    assert_eq!(reader.contigs()[0].name, "contig-a");
    assert_eq!(
        reader.read_sequence(0, SequenceRange::new(2, 10)?)?,
        b"GTNNACGT".to_vec()
    );
    assert_eq!(
        reader.seed_occurrences(SeedQuery {
            k: 31,
            hash: 17,
            canonical_kmer: 5,
        })?,
        vec![SeedOccurrence {
            contig_id: 0,
            position: 4,
            reverse: false,
        }]
    );
    Ok(())
}

#[test]
fn packed_sequence_roundtrip_handles_ambiguity_and_lowercase() {
    let input = b"acgtnuX";
    let (packed, mask) = pack_bases(input);
    let decoded = jam_rs::jma::sequence::unpack_bases(&packed, &mask, input.len() as u64).unwrap();
    assert_eq!(decoded, b"ACGTNTN");
}
