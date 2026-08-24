use jam_rs::resource::{
    ByteRange, RangeReader, ResourceError, ResourceLocator, ResourceMetadata, ResourceMetrics,
};
use jam_rs::router::format::{
    KeyBlockConfig, MetagenomeEntry, PackedKeyWidth, PostingHeader, RouterFormatError,
    range_in_object,
};
use jam_rs::router::reader::RouterReader;
use jam_rs::router::writer::{RouterBuildInput, RouterKeyInput, RouterWriterOptions, build_router};
use jam_rs::router::{HashAlgorithmId, WITNESS_K, WitnessKey, WitnessScheme};
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
    fn new(bytes: Vec<u8>) -> Self {
        Self {
            locator: ResourceLocator::parse("file:///router-format-roundtrip").unwrap(),
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
            etag: Some("router-fixture".into()),
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

fn fixture() -> RouterBuildInput {
    let scheme = WitnessScheme {
        scheme_id: 7,
        k: WITNESS_K,
        base_scale: 20,
        available_scales: vec![20, 50, 100, 200],
        hash_id: HashAlgorithmId::JamhashU64V1,
        zero_excluded: true,
    };
    let mut keys = Vec::new();
    for (index, packed) in [7, 11, 17, 23, 29].into_iter().enumerate() {
        let key = WitnessKey::from_packed(packed).unwrap();
        let mut input = RouterKeyInput::new(key, vec![index as u8, 0xaa, 0x55]);
        input.document_frequency = (index + 1) as u32;
        if index == 0 {
            input.flags = PostingHeader::FLAG_POSITION_BEARING;
            input.positions = Some(vec![0x10, 0x20, 0x30]);
        }
        keys.push(input);
    }
    RouterBuildInput {
        metagenomes: vec![MetagenomeEntry {
            id: "mg-1".into(),
            object_uri: "s3://bucket/mg-1.jma".into(),
            checksum: [3; 32],
        }],
        schemes: vec![scheme],
        keys,
    }
}

fn options(width: PackedKeyWidth) -> RouterWriterOptions {
    RouterWriterOptions {
        key_blocks: KeyBlockConfig {
            target_block_bytes: 4 * 1024,
            packed_width: width,
            store_jamhash: false,
        },
        ..RouterWriterOptions::default()
    }
}

#[test]
fn deterministic_roundtrip_supports_six_and_eight_byte_keys() {
    for width in [PackedKeyWidth::Six, PackedKeyWidth::Eight] {
        let first = build_router(&fixture(), options(width)).unwrap();
        let second = build_router(&fixture(), options(width)).unwrap();
        assert_eq!(first, second);
        let directory = tempdir().unwrap();
        let path = directory.path().join("collection.jwr");
        std::fs::write(&path, first).unwrap();
        let reader = RouterReader::open_mmap(&path).unwrap();
        assert_eq!(
            reader.header().object_size,
            std::fs::metadata(&path).unwrap().len()
        );
        assert_eq!(reader.metagenomes().len(), 1);
        assert_eq!(reader.schemes()[0].available_scales, vec![20, 50, 100, 200]);
        assert_eq!(reader.key_count(), 5);
        let key = WitnessKey::from_packed(7).unwrap();
        let found = reader.lookup(key).unwrap().unwrap();
        assert_eq!(found.key, key);
        let posting = reader.read_posting(found).unwrap();
        assert_eq!(posting.posting, vec![0, 0xaa, 0x55]);
        assert_eq!(posting.positions, Some(vec![0x10, 0x20, 0x30]));
    }
}

#[test]
fn remote_open_reads_directories_but_filter_miss_skips_key_page() {
    let bytes = build_router(&fixture(), options(PackedKeyWidth::Six)).unwrap();
    let resource = TrackingResource::new(bytes);
    let reader = RouterReader::open(resource.clone()).unwrap();
    let opened = resource.ranges();
    let absent = (1..u64::MAX)
        .map(|hash| WitnessKey {
            packed: 7,
            jamhash: hash,
        })
        .find(|key| !reader.filter_maybe_contains(*key))
        .expect("fixture must have a filter miss");
    assert_eq!(reader.lookup(absent).unwrap(), None);
    assert_eq!(
        resource.ranges(),
        opened,
        "filter miss must not read an exact key page"
    );
    let found = reader.lookup(WitnessKey::from_packed(7).unwrap()).unwrap();
    assert!(found.is_some());
    assert!(resource.ranges().len() > opened.len());
}

#[test]
fn digest_match_without_packed_equality_is_not_a_witness() {
    let bytes = build_router(&fixture(), options(PackedKeyWidth::Six)).unwrap();
    let reader = RouterReader::open(TrackingResource::new(bytes)).unwrap();
    let stored = WitnessKey::from_packed(7).unwrap();
    let collision = WitnessKey {
        packed: stored.packed + 1,
        jamhash: stored.jamhash,
    };
    assert_eq!(reader.lookup(collision).unwrap(), None);
}

#[test]
fn mmap_rejects_superblock_corruption() {
    let mut bytes = build_router(&fixture(), options(PackedKeyWidth::Six)).unwrap();
    bytes[16] ^= 1;
    let directory = tempdir().unwrap();
    let path = directory.path().join("corrupt.jwr");
    std::fs::write(&path, bytes).unwrap();
    let error = match RouterReader::open_mmap(path) {
        Ok(_) => panic!("corrupt superblock unexpectedly opened"),
        Err(error) => error,
    };
    assert!(error.to_string().contains("checksum"));
}

#[test]
fn invalid_key_width_is_rejected_before_lookup() {
    let mut bytes = build_router(&fixture(), options(PackedKeyWidth::Six)).unwrap();
    let reader = RouterReader::open(TrackingResource::new(bytes.clone())).unwrap();
    let key_section = reader
        .sections()
        .get(jam_rs::router::format::SECTION_KEYS)
        .unwrap();
    bytes[key_section.offset as usize + 4] = 5;
    let error = match RouterReader::open(TrackingResource::new(bytes)) {
        Ok(_) => panic!("invalid key width unexpectedly opened"),
        Err(error) => error,
    };
    assert!(matches!(
        error,
        jam_rs::router::reader::RouterReaderError::Format(RouterFormatError::InvalidKeyBlock)
            | jam_rs::router::reader::RouterReaderError::Format(
                RouterFormatError::UnsupportedKeyWidth(5)
            )
    ));
}

#[test]
fn section_and_posting_checksums_fail_closed() {
    let bytes = build_router(&fixture(), options(PackedKeyWidth::Six)).unwrap();
    let reader = RouterReader::open(TrackingResource::new(bytes.clone())).unwrap();
    let filter = reader
        .sections()
        .get(jam_rs::router::format::SECTION_FILTER)
        .unwrap();
    let mut filter_corrupt = bytes.clone();
    filter_corrupt[filter.offset as usize + 32] ^= 1;
    let error = match RouterReader::open(TrackingResource::new(filter_corrupt)) {
        Ok(_) => panic!("corrupt filter unexpectedly opened"),
        Err(error) => error,
    };
    assert!(error.to_string().contains("checksum"));

    let posting_section = reader
        .sections()
        .get(jam_rs::router::format::SECTION_POSTING_PAYLOAD)
        .unwrap();
    let mut posting_corrupt = bytes;
    posting_corrupt[posting_section.offset as usize] ^= 1;
    let posting_reader = RouterReader::open(TrackingResource::new(posting_corrupt)).unwrap();
    let match_ = [7, 11, 17, 23, 29]
        .into_iter()
        .filter_map(|packed| {
            posting_reader
                .lookup(WitnessKey::from_packed(packed).unwrap())
                .unwrap()
        })
        .find(|match_| match_.key_index == 0)
        .expect("fixture has a first posting");
    let error = match posting_reader.read_posting(match_) {
        Ok(_) => panic!("corrupt posting unexpectedly decoded"),
        Err(error) => error,
    };
    assert!(error.to_string().contains("checksum"));
}

#[test]
fn range_bounds_are_checked_before_object_access() {
    assert!(range_in_object(u64::MAX, 1, u64::MAX, "fixture").is_err());
    assert!(range_in_object(9, 2, 10, "fixture").is_err());
    assert!(range_in_object(8, 2, 10, "fixture").is_ok());
}
