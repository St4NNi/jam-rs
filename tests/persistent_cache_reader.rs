use jam_rs::resource::persistent_cache::PersistentPageCache;
use jam_rs::resource::persistent_cache::cached_reader::{
    CachedRangeReader, PageDescriptor, exact_page,
};
use jam_rs::resource::{
    ByteRange, RangeReader, ResourceError, ResourceLocator, ResourceMetadata, ResourceMetrics,
    ResourceResult,
};
use sha2::{Digest, Sha256};
use std::io::{Cursor, Read};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use tempfile::{Builder, TempDir};

#[derive(Clone, Debug)]
struct MemoryReader {
    bytes: Arc<Vec<u8>>,
    reads: Arc<AtomicUsize>,
    locator: ResourceLocator,
}

impl MemoryReader {
    fn new(bytes: &[u8]) -> Self {
        Self {
            bytes: Arc::new(bytes.to_vec()),
            reads: Arc::new(AtomicUsize::new(0)),
            locator: ResourceLocator::parse("https://example.org/object").unwrap(),
        }
    }

    fn read_count(&self) -> usize {
        self.reads.load(Ordering::Relaxed)
    }
}

impl RangeReader for MemoryReader {
    fn locator(&self) -> &ResourceLocator {
        &self.locator
    }

    fn metadata(&self) -> ResourceResult<ResourceMetadata> {
        Ok(ResourceMetadata {
            size: self.bytes.len() as u64,
            etag: Some("fixture".to_string()),
            last_modified: None,
            accepts_ranges: true,
        })
    }

    fn read_range(&self, range: ByteRange) -> ResourceResult<Vec<u8>> {
        let end = range.end()?;
        if end > self.bytes.len() as u64 {
            return Err(ResourceError::RangeOutOfBounds {
                offset: range.offset,
                length: range.length,
                size: self.bytes.len() as u64,
            });
        }
        self.reads.fetch_add(1, Ordering::Relaxed);
        Ok(self.bytes[range.offset as usize..end as usize].to_vec())
    }

    fn stream(&self) -> ResourceResult<Box<dyn Read + Send>> {
        Ok(Box::new(Cursor::new((*self.bytes).clone())))
    }

    fn metrics(&self) -> ResourceMetrics {
        ResourceMetrics::default()
    }
}

fn cache_directory() -> TempDir {
    Builder::new()
        .prefix("persistent-cache-reader-")
        .tempdir_in("target")
        .expect("target must be available for cache fixtures")
}

fn descriptor(payload: &[u8]) -> PageDescriptor {
    PageDescriptor {
        section_id: 4,
        page_id: 9,
        page_checksum: Sha256::digest(payload).into(),
    }
}

fn request() -> ByteRange {
    ByteRange::new(2, 4).unwrap()
}

#[test]
fn warm_reconstructed_reader_avoids_inner_range_read() {
    let directory = cache_directory();
    let cache = Arc::new(PersistentPageCache::open(directory.path(), 1024).unwrap());
    let source = b"abcdefgh";
    let payload = &source[2..6];
    let page = descriptor(payload);
    let first_inner = MemoryReader::new(source);
    let first = CachedRangeReader::new(
        first_inner.clone(),
        Arc::clone(&cache),
        [1; 32],
        exact_page(page),
    );
    assert_eq!(first.read_range(request()).unwrap(), payload);
    assert_eq!(first_inner.read_count(), 1);

    drop(first);
    let second_inner = MemoryReader::new(source);
    let second = CachedRangeReader::new(
        second_inner.clone(),
        Arc::clone(&cache),
        [1; 32],
        exact_page(page),
    );
    assert_eq!(second.read_range(request()).unwrap(), payload);
    assert_eq!(second_inner.read_count(), 0);
    let metrics = second.metrics();
    assert_eq!(metrics.cache_hits, 1);
    assert_eq!(metrics.cache_bytes, payload.len() as u64);
    assert_eq!(metrics.remote_bytes, 0);
}

#[test]
fn stale_identity_and_corrupt_pages_are_refetched() {
    let directory = cache_directory();
    let cache = Arc::new(PersistentPageCache::open(directory.path(), 1024).unwrap());
    let source = b"abcdefgh";
    let payload = &source[2..6];
    let page = descriptor(payload);
    let first_inner = MemoryReader::new(source);
    let first = CachedRangeReader::new(first_inner, Arc::clone(&cache), [2; 32], exact_page(page));
    first.read_range(request()).unwrap();

    let stale_inner = MemoryReader::new(source);
    let stale = CachedRangeReader::new(
        stale_inner.clone(),
        Arc::clone(&cache),
        [3; 32],
        exact_page(page),
    );
    assert_eq!(stale.read_range(request()).unwrap(), payload);
    assert_eq!(stale_inner.read_count(), 1);
    assert_eq!(stale.metrics().stale_cache_rejections, 1);

    let corrupt_inner = MemoryReader::new(source);
    let corrupt_reader = CachedRangeReader::new(
        corrupt_inner.clone(),
        Arc::clone(&cache),
        [3; 32],
        exact_page(page),
    );
    let key_path = cache.entry_path(
        &jam_rs::resource::persistent_cache::PageKey::new(
            [3; 32],
            page.section_id,
            page.page_id,
            request(),
            page.page_checksum,
        )
        .unwrap(),
    );
    std::fs::write(key_path, b"corrupt").unwrap();
    assert_eq!(corrupt_reader.read_range(request()).unwrap(), payload);
    assert_eq!(corrupt_inner.read_count(), 1);
    assert_eq!(corrupt_reader.metrics().stale_cache_rejections, 1);
}

#[test]
fn checksum_mismatch_from_inner_reader_fails_closed_without_caching() {
    let directory = cache_directory();
    let cache = Arc::new(PersistentPageCache::open(directory.path(), 1024).unwrap());
    let expected = b"cdef";
    let page = descriptor(expected);
    let reader = CachedRangeReader::new(
        MemoryReader::new(b"abXXXXgh"),
        Arc::clone(&cache),
        [4; 32],
        exact_page(page),
    );
    assert!(matches!(
        reader.read_range(request()),
        Err(ResourceError::CacheIdentityChanged(_))
    ));
    assert_eq!(cache.stats().unwrap().entries, 0);
}

#[test]
fn unresolved_ranges_bypass_cache_without_claiming_a_hit() {
    let directory = cache_directory();
    let cache = Arc::new(PersistentPageCache::open(directory.path(), 1024).unwrap());
    let inner = MemoryReader::new(b"abcdefgh");
    let reader = CachedRangeReader::new(inner.clone(), cache, [5; 32], |_range: ByteRange| None);
    assert_eq!(reader.read_range(request()).unwrap(), b"cdef");
    assert_eq!(inner.read_count(), 1);
    let metrics = reader.metrics();
    assert_eq!(metrics.cache_hits, 0);
    assert_eq!(metrics.cache_misses, 0);
}
