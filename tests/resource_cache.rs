use jam_rs::resource::cache::BlockCache;
use jam_rs::resource::{ByteRange, CacheIdentity};

fn identity(version: &str) -> CacheIdentity {
    CacheIdentity {
        redacted_locator: "file:///fixture".to_string(),
        version: version.to_string(),
        size: 12,
    }
}

#[test]
fn cache_is_bounded_and_identity_aware() {
    let cache = BlockCache::new(4, 8).unwrap();
    let first = identity("one");
    cache.prepare(&first).unwrap();
    cache.insert(&first, 0, b"abcd".to_vec()).unwrap();
    cache.insert(&first, 4, b"efgh".to_vec()).unwrap();
    assert_eq!(cache.bytes().unwrap(), 8);
    cache.insert(&first, 8, b"ijkl".to_vec()).unwrap();
    assert_eq!(cache.bytes().unwrap(), 8);
    assert!(cache.get(&first, 0).unwrap().is_none());

    let second = identity("two");
    assert!(cache.get(&second, 4).unwrap().is_none());
    assert!(cache.is_empty());
}

#[test]
fn block_offsets_cover_checked_half_open_ranges() {
    let cache = BlockCache::new(4, 32).unwrap();
    assert_eq!(
        cache.block_offsets(ByteRange::new(3, 7).unwrap()).unwrap(),
        vec![0, 4, 8]
    );
    assert_eq!(
        cache.block_offsets(ByteRange::new(11, 2).unwrap()).unwrap(),
        vec![8, 12]
    );
    assert!(
        cache
            .block_offsets(ByteRange {
                offset: u64::MAX,
                length: 1,
            })
            .is_err()
    );
}
