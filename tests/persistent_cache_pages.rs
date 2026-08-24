use jam_rs::resource::ByteRange;
use jam_rs::resource::persistent_cache::{
    CacheLookup, CacheMissReason, CacheStoreOutcome, PageKey, PersistentCacheOptions,
    PersistentPageCache,
};
use sha2::{Digest, Sha256};
use std::fs;
use std::io::Write;
use tempfile::{Builder, TempDir};

fn cache_directory() -> TempDir {
    Builder::new()
        .prefix("persistent-page-cache-")
        .tempdir_in("target")
        .expect("target must be available for cache fixtures")
}

fn key(object: u8, page: u64, offset: u64, payload: &[u8]) -> PageKey {
    let checksum: [u8; 32] = Sha256::digest(payload).into();
    PageKey::new(
        [object; 32],
        7,
        page,
        ByteRange::new(offset, payload.len() as u64).unwrap(),
        checksum,
    )
    .unwrap()
}

#[test]
fn cold_write_warm_lookup_survives_cache_reconstruction() {
    let directory = cache_directory();
    let payload = b"cold-to-warm";
    let page_key = key(1, 4, 128, payload);
    let first = PersistentPageCache::open(directory.path(), 1024).unwrap();
    assert_eq!(
        first.put(&page_key, payload).unwrap(),
        CacheStoreOutcome::Stored { evicted: 0 }
    );
    assert!(matches!(
        first.get(&page_key).unwrap(),
        CacheLookup::Hit(bytes) if bytes == payload
    ));
    drop(first);

    let second = PersistentPageCache::open(directory.path(), 1024).unwrap();
    assert!(matches!(
        second.get(&page_key).unwrap(),
        CacheLookup::Hit(bytes) if bytes == payload
    ));
    assert_eq!(second.stats().unwrap().entries, 1);
    let on_disk_bytes = fs::metadata(second.entry_path(&page_key)).unwrap().len();
    assert_eq!(second.stats().unwrap().bytes, on_disk_bytes);
}

#[test]
fn changed_object_identity_discards_old_page_and_signals_refetch() {
    let directory = cache_directory();
    let old_payload = b"old-page";
    let new_payload = b"new-page";
    let old_key = key(2, 9, 32, old_payload);
    let new_key = key(3, 9, 32, new_payload);
    let cache = PersistentPageCache::open(directory.path(), 1024).unwrap();
    cache.put(&old_key, old_payload).unwrap();

    assert_eq!(
        cache.get(&new_key).unwrap(),
        CacheLookup::Miss(CacheMissReason::StaleIdentity)
    );
    assert!(!cache.entry_path(&old_key).exists());
    assert_eq!(cache.stats().unwrap().stale_identity_rejections, 1);
    cache.put(&new_key, new_payload).unwrap();
    assert!(matches!(
        cache.get(&new_key).unwrap(),
        CacheLookup::Hit(bytes) if bytes == new_payload
    ));
}

#[test]
fn corrupt_page_is_removed_and_requires_refetch() {
    let directory = cache_directory();
    let payload = b"checksum-protected";
    let page_key = key(4, 2, 0, payload);
    let cache = PersistentPageCache::open(directory.path(), 1024).unwrap();
    cache.put(&page_key, payload).unwrap();
    let path = cache.entry_path(&page_key);
    let mut file = fs::OpenOptions::new().write(true).open(&path).unwrap();
    file.write_all(b"corrupt").unwrap();
    file.flush().unwrap();

    assert_eq!(
        cache.get(&page_key).unwrap(),
        CacheLookup::Miss(CacheMissReason::Corrupt)
    );
    assert!(!path.exists());
    assert_eq!(cache.stats().unwrap().corrupt_rejections, 1);
}

#[test]
fn eviction_is_bounded_and_oldest_sequence_is_deterministic() {
    let directory = cache_directory();
    let first_payload = b"1111";
    let second_payload = b"2222";
    let third_payload = b"3333";
    let first_key = key(5, 1, 0, first_payload);
    let second_key = key(5, 2, 4, second_payload);
    let third_key = key(5, 3, 8, third_payload);
    let cache = PersistentPageCache::open(directory.path(), 250).unwrap();
    cache.put(&first_key, first_payload).unwrap();
    cache.put(&second_key, second_payload).unwrap();
    assert_eq!(
        cache.put(&third_key, third_payload).unwrap(),
        CacheStoreOutcome::Stored { evicted: 1 }
    );
    assert!(cache.stats().unwrap().bytes <= 250);
    assert_eq!(cache.stats().unwrap().entries, 2);
    assert_eq!(
        cache.get(&first_key).unwrap(),
        CacheLookup::Miss(CacheMissReason::NotFound)
    );
    assert!(matches!(
        cache.get(&second_key).unwrap(),
        CacheLookup::Hit(bytes) if bytes == second_payload
    ));
    assert!(matches!(
        cache.get(&third_key).unwrap(),
        CacheLookup::Hit(bytes) if bytes == third_payload
    ));
}

#[test]
fn physical_cache_files_stay_within_the_full_entry_bound() {
    let directory = cache_directory();
    let payloads = [b"aaaa".as_slice(), b"bbbb".as_slice(), b"cccc".as_slice()];
    let keys = [
        key(8, 1, 0, payloads[0]),
        key(8, 2, 4, payloads[1]),
        key(8, 3, 8, payloads[2]),
    ];
    let max_bytes = 250;
    let cache = PersistentPageCache::open(directory.path(), max_bytes).unwrap();
    for (page_key, payload) in keys.iter().zip(payloads) {
        assert!(matches!(
            cache.put(page_key, payload).unwrap(),
            CacheStoreOutcome::Stored { .. }
        ));
    }

    let physical_bytes = fs::read_dir(directory.path())
        .unwrap()
        .map(|entry| {
            let entry = entry.unwrap();
            if entry.path().extension().and_then(|value| value.to_str()) == Some("jpc") {
                entry.metadata().unwrap().len()
            } else {
                0
            }
        })
        .sum::<u64>();
    assert!(physical_bytes <= max_bytes);
    assert_eq!(physical_bytes, cache.stats().unwrap().bytes);
}

#[test]
fn disabled_cache_performs_no_filesystem_operations() {
    let directory = cache_directory();
    let cache_path = directory.path().join("pages");
    let page_key = key(6, 1, 0, b"disabled");
    let cache =
        PersistentPageCache::open_with_options(&cache_path, PersistentCacheOptions::disabled())
            .unwrap();
    assert_eq!(
        cache.put(&page_key, b"disabled").unwrap(),
        CacheStoreOutcome::Disabled
    );
    assert_eq!(
        cache.get(&page_key).unwrap(),
        CacheLookup::Miss(CacheMissReason::Disabled)
    );
    assert!(!cache_path.exists());
}
