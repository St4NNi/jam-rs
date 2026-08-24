use jam_rs::jma::ArchiveReader;
use jam_rs::jma::format::JMA_FORMAT_MAGIC;
use jam_rs::jma::index::{decode_seed_index_directory, encode_seed_index};
use jam_rs::jma::reader::JmaReader;
use jam_rs::jma::writer::{ArchiveParts, SeedLevelData, SeedRecord, SeedSection, encode_archive};
use jam_rs::jma::{ContigMetadata, SeedLevel, SeedOccurrence, SeedQuery, SequenceRange};
use jam_rs::resource::{
    ByteRange, RangeReader, ResourceError, ResourceLocator, ResourceMetadata, ResourceMetrics,
};
use jam_rs::sequence::{BlockCodec, encode_sequence_block};
use jam_rs::sequence::{decode_reverse_complement_range, encode_contig};
use std::io::{Cursor, Read};
use std::sync::Arc;

struct MemoryResource {
    locator: ResourceLocator,
    bytes: Arc<Vec<u8>>,
}

impl MemoryResource {
    fn new(bytes: Vec<u8>) -> Self {
        Self {
            locator: ResourceLocator::parse("file:///jma-format-roundtrip").unwrap(),
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
            etag: None,
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
    let sequence = b"ACGTNNACGTACGTTTGCAACGTACGTACGT";
    let sequence_block = encode_sequence_block(0, 0, sequence, BlockCodec::Raw2Bit).unwrap();
    ArchiveParts {
        flags: 7,
        source_sha256: [9; 32],
        contigs: vec![ContigMetadata {
            id: 0,
            name: "contig-a".into(),
            length: sequence.len() as u64,
        }],
        sequence_blocks: vec![sequence_block],
        seed_sections: vec![
            SeedSection {
                k: 21,
                levels: vec![SeedLevelData {
                    level: SeedLevel { k: 21, scale: 1 },
                    records: vec![SeedRecord {
                        query: SeedQuery {
                            k: 21,
                            hash: 0x1000,
                            canonical_kmer: 3,
                        },
                        occurrence: SeedOccurrence {
                            contig_id: 0,
                            position: 2,
                            reverse: true,
                        },
                    }],
                }],
            },
            SeedSection {
                k: 22,
                levels: vec![SeedLevelData {
                    level: SeedLevel { k: 22, scale: 2 },
                    records: vec![SeedRecord {
                        query: SeedQuery {
                            k: 22,
                            hash: 0x2000,
                            canonical_kmer: 7,
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
fn format_one_is_deterministic_and_self_contained() {
    let first = encode_archive(&fixture()).unwrap();
    let second = encode_archive(&fixture()).unwrap();
    assert_eq!(first, second);
    assert_eq!(&first[..JMA_FORMAT_MAGIC.len()], &JMA_FORMAT_MAGIC);
    let reader = JmaReader::open(MemoryResource::new(first)).unwrap();
    assert_eq!(reader.metadata_fields().source_assembly_sha256, [9; 32]);
    assert_eq!(reader.metadata_fields().hash_algorithm, "jamhash_u64_v1");
    assert!(reader.metadata_fields().archive_sha256.is_some());
    assert_eq!(reader.seed_schemes().count(), 2);
    assert_eq!(
        reader
            .read_sequence(0, SequenceRange::new(4, 12).unwrap())
            .unwrap(),
        b"NNACGTAC"
    );
}

#[test]
fn seed_lookup_checks_exact_key_after_digest() {
    let mut parts = fixture();
    parts.seed_sections[0].levels[0].records.push(SeedRecord {
        query: SeedQuery {
            k: 21,
            hash: 0x1000,
            canonical_kmer: 4,
        },
        occurrence: SeedOccurrence {
            contig_id: 0,
            position: 7,
            reverse: false,
        },
    });
    let reader = JmaReader::open(MemoryResource::new(encode_archive(&parts).unwrap())).unwrap();
    let matches = reader
        .seed_occurrences(SeedQuery {
            k: 21,
            hash: 0x1000,
            canonical_kmer: 3,
        })
        .unwrap();
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0].position, 2);
}

#[test]
fn seed_pages_never_cross_hash_prefix_boundaries() {
    let records = (0..4097u64)
        .map(|index| {
            let prefix = if index < 4096 { 0x100u64 } else { 0x101u64 };
            SeedRecord {
                query: SeedQuery {
                    k: 21,
                    hash: (prefix << 52) | index,
                    canonical_kmer: index & ((1u64 << 42) - 1),
                },
                occurrence: SeedOccurrence {
                    contig_id: 0,
                    position: index,
                    reverse: false,
                },
            }
        })
        .collect();
    let bytes = encode_seed_index(&[SeedSection {
        k: 21,
        levels: vec![SeedLevelData {
            level: SeedLevel { k: 21, scale: 1 },
            records,
        }],
    }])
    .unwrap();
    let index = decode_seed_index_directory(&bytes).unwrap();
    let pages = &index.schemes[0].pages;
    assert_eq!(pages.len(), 2);
    assert!(
        pages
            .iter()
            .all(|page| page.first_hash >> 52 == page.last_hash >> 52)
    );
}

#[test]
fn mixed_raw_zstd_blocks_preserve_iupac_and_reverse_complement() {
    let left = b"ACGTRYK";
    let right = b"MBDHVNAC";
    let mut complete = left.to_vec();
    complete.extend_from_slice(right);
    let raw = encode_sequence_block(0, 0, left, BlockCodec::Raw2Bit).unwrap();
    let zstd = encode_sequence_block(0, left.len() as u64, right, BlockCodec::Zstd2Bit).unwrap();
    let reader = JmaReader::open(MemoryResource::new(
        encode_archive(&ArchiveParts {
            flags: 0,
            source_sha256: [4; 32],
            contigs: vec![ContigMetadata {
                id: 0,
                name: "mixed".into(),
                length: complete.len() as u64,
            }],
            sequence_blocks: vec![raw, zstd],
            seed_sections: Vec::new(),
        })
        .unwrap(),
    ))
    .unwrap();
    assert_eq!(
        reader
            .read_sequence(0, SequenceRange::new(0, complete.len() as u64).unwrap())
            .unwrap(),
        complete
    );
    let encoded = encode_contig(&complete).unwrap();
    assert_eq!(
        reader
            .read_sequence_reverse_complement(
                0,
                SequenceRange::new(0, complete.len() as u64).unwrap(),
            )
            .unwrap(),
        decode_reverse_complement_range(&encoded, 0..complete.len() as u64).unwrap()
    );
}
