pub use jam_rs::jamhash_u64_v1;

pub mod jma {
    pub use jam_rs::jma::{JmaError, JmaResult, SeedOccurrence, SeedQuery};

    pub mod seeds {
        pub use crate::compact_seeds::*;
    }

    pub mod format {
        pub use jam_rs::jma::format::*;
    }

    pub mod writer {
        pub use jam_rs::jma::writer::*;
    }
}

pub mod resource {
    pub use jam_rs::resource::*;
}

#[path = "../src/jma/seeds/mod.rs"]
mod compact_seeds;

#[path = "../src/jma/reader/seed_compact.rs"]
mod remote_seed;

#[path = "../src/jma/writer/seed_compact.rs"]
mod compact_writer;

use compact_seeds::{
    KeyEncoding, LookupOptions, OccurrenceEncoding, PostingClass, SeedBinding, SeedBuildConfig,
    SeedCollection, SeedIndex, encode_seed_collection, encode_seed_section,
};
use jam_rs::jma::writer::{SeedLevelData, SeedSection};
use jam_rs::jma::{SeedLevel, SeedOccurrence, SeedQuery};
use jam_rs::resource::{
    ByteRange, RangeReader, ResourceError, ResourceLocator, ResourceMetadata, ResourceMetrics,
};
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
            locator: ResourceLocator::parse("file:///compact-seed").unwrap(),
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
            etag: Some("compact".to_string()),
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

fn key(value: u64) -> (u64, u64) {
    let packed = value.max(1);
    let hash = jam_hash(packed);
    (hash, packed)
}

fn jam_hash(packed: u64) -> u64 {
    let hash = jam_rs::jamhash_u64_v1(packed);
    assert_ne!(hash, 0);
    hash
}

fn records() -> Vec<jam_rs::jma::writer::SeedRecord> {
    let (hash_a, packed_a) = key(0x1234);
    let (hash_b, packed_b) = key(0x2345);
    let (hash_c, packed_c) = key(0x3456);
    vec![
        jam_rs::jma::writer::SeedRecord {
            query: SeedQuery {
                k: 21,
                hash: hash_b,
                canonical_kmer: packed_b,
            },
            occurrence: SeedOccurrence {
                contig_id: 1,
                position: 40,
                reverse: true,
            },
        },
        jam_rs::jma::writer::SeedRecord {
            query: SeedQuery {
                k: 21,
                hash: hash_a,
                canonical_kmer: packed_a,
            },
            occurrence: SeedOccurrence {
                contig_id: 0,
                position: 10,
                reverse: false,
            },
        },
        jam_rs::jma::writer::SeedRecord {
            query: SeedQuery {
                k: 21,
                hash: hash_c,
                canonical_kmer: packed_c,
            },
            occurrence: SeedOccurrence {
                contig_id: 2,
                position: 80,
                reverse: false,
            },
        },
        jam_rs::jma::writer::SeedRecord {
            query: SeedQuery {
                k: 21,
                hash: hash_a,
                canonical_kmer: packed_a,
            },
            occurrence: SeedOccurrence {
                contig_id: 0,
                position: 10,
                reverse: false,
            },
        },
    ]
}

fn binding() -> SeedBinding {
    SeedBinding {
        scheme_id: 7,
        k: 21,
        scale: 50,
        catalog_checksum: [1; 32],
        archive_checksum: [2; 32],
        scheme_checksum: [3; 32],
    }
}

fn binding31() -> SeedBinding {
    SeedBinding {
        scheme_id: 11,
        k: 31,
        scale: 100,
        catalog_checksum: [1; 32],
        archive_checksum: [2; 32],
        scheme_checksum: [4; 32],
    }
}

fn query(value: u64) -> SeedQuery {
    let (hash, packed) = key(value);
    SeedQuery {
        k: 21,
        hash,
        canonical_kmer: packed,
    }
}

fn records31() -> Vec<jam_rs::jma::writer::SeedRecord> {
    let (hash, packed) = key(0x789a);
    vec![jam_rs::jma::writer::SeedRecord {
        query: SeedQuery {
            k: 31,
            hash,
            canonical_kmer: packed,
        },
        occurrence: SeedOccurrence {
            contig_id: 4,
            position: 12,
            reverse: true,
        },
    }]
}

#[test]
fn compact_roundtrip_is_deterministic() {
    let first = encode_seed_section(&records(), binding(), SeedBuildConfig::default()).unwrap();
    let mut reversed = records();
    reversed.reverse();
    let second = encode_seed_section(&reversed, binding(), SeedBuildConfig::default()).unwrap();
    assert_eq!(first, second);
    let index = SeedIndex::open(first, Some(binding())).unwrap();
    let result = index
        .lookup(query(0x1234), LookupOptions::default())
        .unwrap();
    assert_eq!(result.matches.len(), 1);
    assert_eq!(result.matches[0].occurrence_count, 1);
    assert_eq!(result.matches[0].occurrences[0].position, 10);
    assert_eq!(result.matches[0].class, PostingClass::Rare);
}

#[test]
fn collection_separates_scheme_ranges() {
    let blob21 = encode_seed_section(&records(), binding(), SeedBuildConfig::default()).unwrap();
    let blob31 =
        encode_seed_section(&records31(), binding31(), SeedBuildConfig::default()).unwrap();
    let bytes = encode_seed_collection(&[(binding(), blob21), (binding31(), blob31)]).unwrap();
    let collection = SeedCollection::open(bytes.clone()).unwrap();
    assert_eq!(collection.entries().len(), 2);
    let local = collection.open_scheme(7, Some(binding())).unwrap();
    assert_eq!(
        local
            .lookup(query(0x1234), LookupOptions::default())
            .unwrap()
            .matches
            .len(),
        1
    );

    let resource = TrackingResource::new(bytes.clone());
    let remote = remote_seed::RemoteSeedCollection::open(
        resource.clone(),
        0,
        u64::try_from(bytes.len()).unwrap(),
    )
    .unwrap();
    assert_eq!(remote.header().scheme_count, 2);
    assert_eq!(remote.entries().len(), 2);
    let other = remote.entry(11).unwrap();
    let selected = remote.open_scheme(7, Some(binding())).unwrap();
    assert_eq!(
        selected
            .lookup(query(0x1234), LookupOptions::default())
            .unwrap()
            .matches
            .len(),
        1
    );
    let other_end = other.offset.checked_add(other.length).unwrap();
    assert!(
        resource
            .ranges()
            .iter()
            .all(|range| { range.end().unwrap() <= other.offset || range.offset >= other_end })
    );
}

#[test]
fn remote_skip_does_not_read_position_pages() {
    let (hash, packed) = key(0x4567);
    let records = (0..8)
        .map(|position| jam_rs::jma::writer::SeedRecord {
            query: SeedQuery {
                k: 21,
                hash,
                canonical_kmer: packed,
            },
            occurrence: SeedOccurrence {
                contig_id: 1,
                position,
                reverse: false,
            },
        })
        .collect::<Vec<_>>();
    let bytes = encode_seed_section(
        &records,
        binding(),
        SeedBuildConfig {
            rare_occurrence_limit: 1,
            common_occurrence_limit: 2,
            suppressed_occurrence_limit: 4,
            ..SeedBuildConfig::default()
        },
    )
    .unwrap();
    let section_header = compact_seeds::SeedIndex::open(bytes.clone(), None)
        .unwrap()
        .header()
        .to_owned();
    let resource = TrackingResource::new(bytes);
    let remote = remote_seed::RemoteSeedIndex::open(
        resource.clone(),
        0,
        section_header.section_length,
        Some(binding()),
    )
    .unwrap();
    let result = remote
        .lookup(
            SeedQuery {
                k: 21,
                hash,
                canonical_kmer: packed,
            },
            LookupOptions::default(),
        )
        .unwrap();
    assert_eq!(remote.header().scheme_id, binding().scheme_id);
    assert_eq!(result.matches[0].occurrence_count, 8);
    assert_eq!(result.metrics.position_pages_read, 0);
    assert!(
        resource
            .ranges()
            .iter()
            .all(|range| range.offset + range.length <= section_header.position_data_offset)
    );
}

#[test]
fn writer_adapter_preserves_binding() {
    let level = SeedLevelData {
        level: SeedLevel { k: 21, scale: 50 },
        records: records().into_iter().take(1).collect(),
    };
    let section = SeedSection {
        k: 21,
        levels: vec![level.clone()],
    };
    let one =
        compact_writer::encode_seed_level(&section, &level, binding(), SeedBuildConfig::default())
            .unwrap();
    assert_eq!(
        SeedIndex::open(one, Some(binding()))
            .unwrap()
            .header()
            .scheme_id,
        7
    );
    let many =
        compact_writer::encode_seed_levels(&section, &[(binding(), SeedBuildConfig::default())])
            .unwrap();
    assert_eq!(many.len(), 1);
    let collection = compact_writer::encode_seed_collection_levels(
        &section,
        &[(binding(), SeedBuildConfig::default())],
    )
    .unwrap();
    assert_eq!(SeedCollection::open(collection).unwrap().entries().len(), 1);
    let section31 = SeedSection {
        k: 31,
        levels: vec![SeedLevelData {
            level: SeedLevel { k: 31, scale: 100 },
            records: records31(),
        }],
    };
    let collection = compact_writer::encode_seed_collection_sections(
        &[section, section31],
        &[
            (binding(), SeedBuildConfig::default()),
            (binding31(), SeedBuildConfig::default()),
        ],
    )
    .unwrap();
    assert_eq!(SeedCollection::open(collection).unwrap().entries().len(), 2);
}

#[test]
fn grouped_pages_are_read_once() {
    let bytes = encode_seed_section(
        &records(),
        binding(),
        SeedBuildConfig {
            key_page_bytes: 16,
            position_page_bytes: 1024,
            ..SeedBuildConfig::default()
        },
    )
    .unwrap();
    let index = SeedIndex::open(bytes, None).unwrap();
    let result = index
        .lookup_batch(
            &[query(0x1234), query(0x2345), query(0x3456)],
            LookupOptions::default(),
        )
        .unwrap();
    assert_eq!(result.matches.len(), 3);
    assert_eq!(result.metrics.key_pages_read, 2);
    assert_eq!(result.metrics.position_pages_read, 1);
    assert_eq!(result.metrics.position_payload_reads, 3);
}

#[test]
fn repetitive_headers_skip_positions() {
    let (hash, packed) = key(0x4567);
    let records = (0..8)
        .map(|position| jam_rs::jma::writer::SeedRecord {
            query: SeedQuery {
                k: 21,
                hash,
                canonical_kmer: packed,
            },
            occurrence: SeedOccurrence {
                contig_id: 1,
                position,
                reverse: false,
            },
        })
        .collect::<Vec<_>>();
    let bytes = encode_seed_section(
        &records,
        binding(),
        SeedBuildConfig {
            rare_occurrence_limit: 1,
            common_occurrence_limit: 2,
            suppressed_occurrence_limit: 4,
            ..SeedBuildConfig::default()
        },
    )
    .unwrap();
    let index = SeedIndex::open(bytes, None).unwrap();
    let result = index
        .lookup(
            SeedQuery {
                k: 21,
                hash,
                canonical_kmer: packed,
            },
            LookupOptions {
                include_positions: true,
                max_occurrences: Some(2),
            },
        )
        .unwrap();
    assert_eq!(result.matches[0].class, PostingClass::Suppressed);
    assert_eq!(result.matches[0].occurrence_count, 8);
    assert!(result.matches[0].occurrences.is_empty());
    assert_eq!(result.metrics.position_pages_read, 0);
    assert_eq!(result.metrics.position_payload_reads, 0);
}

#[test]
fn exact_key_and_encoding_variants_roundtrip() {
    for key_encoding in [KeyEncoding::Packed8, KeyEncoding::Packed6] {
        for occurrence_encoding in [OccurrenceEncoding::Fixed16, OccurrenceEncoding::DeltaVarint] {
            let bytes = encode_seed_section(
                &records(),
                binding(),
                SeedBuildConfig {
                    key_encoding,
                    occurrence_encoding,
                    ..SeedBuildConfig::default()
                },
            )
            .unwrap();
            let index = SeedIndex::open(bytes, None).unwrap();
            let result = index
                .lookup(query(0x2345), LookupOptions::default())
                .unwrap();
            assert_eq!(result.matches[0].occurrences[0].position, 40);
        }
    }
}

#[test]
fn corruption_collision_and_binding_fail_closed() {
    let bytes = encode_seed_section(&records(), binding(), SeedBuildConfig::default()).unwrap();
    let mut corrupt = bytes.clone();
    corrupt[compact_seeds::COMPACT_HEADER_SIZE] ^= 0x01;
    assert!(SeedIndex::open(corrupt, None).is_err());

    let mut wrong_binding = binding();
    wrong_binding.archive_checksum = [9; 32];
    assert!(SeedIndex::open(bytes.clone(), Some(wrong_binding)).is_err());

    let (hash, packed) = key(0x1234);
    let wrong = SeedQuery {
        k: 21,
        hash,
        canonical_kmer: packed ^ 1,
    };
    assert!(
        SeedIndex::open(bytes, None)
            .unwrap()
            .lookup(wrong, LookupOptions::default())
            .is_err()
    );
}
