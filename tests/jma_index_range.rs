use jam_rs::jma::builder::{ArchiveBuildConfig, write_archive_from_fasta};
use jam_rs::jma::index::{build_index, encode_index};
use jam_rs::jma::reader::JmaReader;
use jam_rs::jma::sequence::{SequenceBlock, pack_bases};
use jam_rs::jma::writer::{ArchiveParts, SeedLevelData, SeedRecord, SeedSection, encode_archive};
use jam_rs::jma::{
    ArchiveReader, ContigMetadata, SeedLevel, SeedOccurrence, SeedQuery, SequenceRange,
};
use jam_rs::resource::{
    ByteRange, RangeReader, ResourceError, ResourceLocator, ResourceMetadata, ResourceMetrics,
};
use sha2::{Digest, Sha256};
use std::io::{Cursor, Read};
use std::sync::{Arc, Mutex};
use tempfile::tempdir;

#[derive(Clone)]
struct TrackingResource {
    locator: ResourceLocator,
    bytes: Arc<Vec<u8>>,
    ranges: Arc<Mutex<Vec<ByteRange>>>,
}

impl TrackingResource {
    fn new(name: &str, bytes: Vec<u8>) -> Self {
        Self {
            locator: ResourceLocator::parse(&format!("file:///{name}")).unwrap(),
            bytes: Arc::new(bytes),
            ranges: Arc::new(Mutex::new(Vec::new())),
        }
    }

    fn ranges(&self) -> Vec<ByteRange> {
        self.ranges.lock().unwrap().clone()
    }
}

impl RangeReader for TrackingResource {
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
        self.ranges.lock().unwrap().push(range);
        Ok(self.bytes[range.offset as usize..end as usize].to_vec())
    }

    fn stream(&self) -> Result<Box<dyn Read + Send>, ResourceError> {
        Ok(Box::new(Cursor::new(self.bytes.as_ref().clone())))
    }

    fn metrics(&self) -> ResourceMetrics {
        ResourceMetrics::default()
    }
}

fn fixture() -> Vec<u8> {
    let first = b"ACGTACGT";
    let second = b"TTGCAAAA";
    let (first_packed, first_mask) = pack_bases(first);
    let (second_packed, second_mask) = pack_bases(second);
    encode_archive(&ArchiveParts {
        flags: 0,
        source_sha256: [7; 32],
        contigs: vec![ContigMetadata {
            id: 0,
            name: "c".to_string(),
            length: (first.len() + second.len()) as u64,
        }],
        sequence_blocks: vec![
            SequenceBlock {
                contig_id: 0,
                start: 0,
                base_length: first.len() as u64,
                packed: first_packed,
                unknown_mask: first_mask,
            },
            SequenceBlock {
                contig_id: 0,
                start: first.len() as u64,
                base_length: second.len() as u64,
                packed: second_packed,
                unknown_mask: second_mask,
            },
        ],
        seed_sections: vec![SeedSection {
            k: 21,
            levels: vec![SeedLevelData {
                level: SeedLevel { k: 21, scale: 1 },
                records: vec![
                    SeedRecord {
                        query: SeedQuery {
                            k: 21,
                            hash: 0x1000_0000_0000_0000,
                            canonical_kmer: 1,
                        },
                        occurrence: SeedOccurrence {
                            contig_id: 0,
                            position: 2,
                            reverse: false,
                        },
                    },
                    SeedRecord {
                        query: SeedQuery {
                            k: 21,
                            hash: 0x9000_0000_0000_0000,
                            canonical_kmer: 2,
                        },
                        occurrence: SeedOccurrence {
                            contig_id: 0,
                            position: 10,
                            reverse: true,
                        },
                    },
                ],
            }],
        }],
    })
    .unwrap()
}

#[test]
fn indexed_open_reads_only_metadata_until_ranges_are_requested() {
    let archive_bytes = fixture();
    let index = build_index(&archive_bytes).unwrap();
    let mut digest = Sha256::new();
    digest.update(&archive_bytes);
    let mut archive_sha256 = [0u8; 32];
    archive_sha256.copy_from_slice(&digest.finalize());
    assert_eq!(index.archive_sha256, archive_sha256);
    let sidecar_bytes = encode_index(&index).unwrap();
    let archive = TrackingResource::new("archive.jma", archive_bytes);
    let index = TrackingResource::new("archive.jma.idx.json", sidecar_bytes);
    let reader = JmaReader::open_indexed(archive.clone(), index.clone()).unwrap();

    let open_ranges = archive.ranges();
    assert!(open_ranges.iter().all(|range| range.length < 1024));
    assert!(index.ranges().iter().any(|range| range.length > 0));

    let sequence = reader
        .read_sequence(0, SequenceRange::new(1, 5).unwrap())
        .unwrap();
    assert_eq!(sequence, b"CGTA");
    let metrics = reader.metrics();
    assert_eq!(metrics.sequence_blocks_read, 1);
    assert!(metrics.decoded_bytes > 0);

    let occurrences = reader
        .seed_occurrences(SeedQuery {
            k: 21,
            hash: 0x1000_0000_0000_0000,
            canonical_kmer: 1,
        })
        .unwrap();
    assert_eq!(occurrences.len(), 1);
    assert_eq!(occurrences[0].position, 2);
    assert_eq!(reader.metrics().seed_buckets_read, 1);

    let fetched = archive.ranges();
    assert_eq!(fetched.len(), open_ranges.len() + 2);
    assert!(fetched.iter().any(|range| range.length > 36));
}

#[test]
fn indexed_payload_checksum_is_checked_before_decode() {
    let mut archive_bytes = fixture();
    let index = build_index(&archive_bytes).unwrap();
    let sequence = index.sequence_blocks[0].clone();
    archive_bytes[sequence.offset as usize] ^= 1;
    let archive = TrackingResource::new("archive.jma", archive_bytes);
    let sidecar = TrackingResource::new("archive.jma.idx.json", encode_index(&index).unwrap());
    let reader = JmaReader::open_indexed(archive, sidecar).unwrap();
    let error = reader
        .read_sequence(0, SequenceRange::new(0, 4).unwrap())
        .unwrap_err();
    assert!(error.to_string().contains("checksum"));
}

#[test]
fn builder_writes_only_the_derived_sidecar() {
    let directory = tempdir().unwrap();
    let input = directory.path().join("assembly.fa");
    let archive = directory.path().join("assembly.jma");
    let unrelated = directory.path().join("notes.txt");
    std::fs::write(
        &input,
        b">contig\nACGTTGCATGTCAGTAGGCATCAGTACCGATGCTAGCTAGGCTAACGTTACGATCGATGCA\n",
    )
    .unwrap();
    std::fs::write(&unrelated, b"keep me").unwrap();
    write_archive_from_fasta(&input, &archive, ArchiveBuildConfig::default()).unwrap();
    let sidecar = archive.with_file_name("assembly.jma.idx.json");
    assert!(sidecar.is_file());
    assert_eq!(std::fs::read(&unrelated).unwrap(), b"keep me");
    let parsed: serde_json::Value =
        serde_json::from_slice(&std::fs::read(sidecar).unwrap()).unwrap();
    assert_eq!(parsed["index_version"], 1);
}
