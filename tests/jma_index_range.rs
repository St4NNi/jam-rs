use jam_rs::jma::ArchiveReader;
use jam_rs::jma::JmaError;
use jam_rs::jma::format::{SECTION_SEEDS, decode_header, decode_section_directory};
use jam_rs::jma::reader::JmaReader;
use jam_rs::jma::writer::{ArchiveParts, SeedLevelData, SeedRecord, SeedSection, encode_archive};
use jam_rs::jma::{ContigMetadata, SeedLevel, SeedOccurrence, SeedQuery, SequenceRange};
use jam_rs::resource::{
    ByteRange, RangeReader, ResourceError, ResourceLocator, ResourceMetadata, ResourceMetrics,
};
use jam_rs::sequence::{BlockCodec, encode_sequence_block};
use std::io::{Cursor, Read};
use std::sync::{Arc, Mutex};

#[derive(Clone)]
struct TrackingResource {
    locator: ResourceLocator,
    bytes: Arc<Vec<u8>>,
    ranges: Arc<Mutex<Vec<ByteRange>>>,
}
impl TrackingResource {
    fn new(bytes: Vec<u8>) -> Self {
        Self {
            locator: ResourceLocator::parse("file:///jma-range").unwrap(),
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
            etag: Some("fixture".into()),
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

fn fixture() -> ArchiveParts {
    let first = b"ACGTACGTACGTACGT";
    let second = b"TTGCAAAATTGCAAAA";
    let first_block = encode_sequence_block(0, 0, first, BlockCodec::Raw2Bit).unwrap();
    let second_block = encode_sequence_block(0, 16, second, BlockCodec::Raw2Bit).unwrap();
    ArchiveParts {
        flags: 0,
        source_sha256: [7; 32],
        contigs: vec![ContigMetadata {
            id: 0,
            name: "c".into(),
            length: 32,
        }],
        sequence_blocks: vec![first_block, second_block],
        seed_sections: vec![SeedSection {
            k: 21,
            levels: vec![SeedLevelData {
                level: SeedLevel { k: 21, scale: 1 },
                records: vec![SeedRecord {
                    query: SeedQuery {
                        k: 21,
                        hash: jam_rs::jamhash_u64_v1(1),
                        canonical_kmer: 1,
                    },
                    occurrence: SeedOccurrence {
                        contig_id: 0,
                        position: 2,
                        reverse: false,
                    },
                }],
            }],
        }],
    }
}

#[test]
fn open_reads_directories_but_not_sequence_or_seed_payloads() {
    let resource = TrackingResource::new(encode_archive(&fixture()).unwrap());
    let reader = JmaReader::open(resource.clone()).unwrap();
    let opened = resource.ranges();
    assert!(opened.iter().all(|range| range.length < 64 * 1024));
    let sequence = reader
        .read_sequence(0, SequenceRange::new(1, 5).unwrap())
        .unwrap();
    assert_eq!(sequence, b"CGTA");
    let occurrences = reader
        .seed_occurrences(SeedQuery {
            k: 21,
            hash: jam_rs::jamhash_u64_v1(1),
            canonical_kmer: 1,
        })
        .unwrap();
    assert_eq!(occurrences.len(), 1);
    assert!(reader.metrics().seed_buckets_read >= 1);
    assert!(reader.metrics().sequence_blocks_read >= 1);
    assert!(resource.ranges().len() > opened.len());
}

#[test]
fn compact_seed_directory_corruption_is_rejected_without_reading_pages() {
    let mut bytes = encode_archive(&fixture()).unwrap();
    let header = decode_header(&bytes[..256]).unwrap();
    let directory_start = usize::try_from(header.section_directory_offset).unwrap();
    let directory_end = directory_start + usize::try_from(header.section_directory_length).unwrap();
    let entries =
        decode_section_directory(&bytes[directory_start..directory_end], header.section_count)
            .unwrap();
    let seed = entries
        .iter()
        .find(|entry| entry.kind == SECTION_SEEDS)
        .unwrap();
    bytes[usize::try_from(seed.offset).unwrap()] ^= 1;
    assert!(matches!(
        JmaReader::open(TrackingResource::new(bytes)),
        Err(JmaError::InvalidMagic)
    ));
}
