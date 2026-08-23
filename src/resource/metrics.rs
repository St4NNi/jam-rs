//! Atomic counters used by resource readers.

use super::ResourceMetrics;
use std::sync::atomic::{AtomicU64, Ordering};

/// Thread-safe counters shared by a resource and its cache-backed reads.
#[derive(Debug, Default)]
pub struct MetricsCounter {
    metadata_requests: AtomicU64,
    head_requests: AtomicU64,
    get_requests: AtomicU64,
    range_requests: AtomicU64,
    stream_requests: AtomicU64,
    requested_bytes: AtomicU64,
    returned_bytes: AtomicU64,
    decoded_bytes: AtomicU64,
    remote_bytes: AtomicU64,
    cache_bytes: AtomicU64,
    cache_hits: AtomicU64,
    cache_misses: AtomicU64,
    cache_evictions: AtomicU64,
    stale_cache_rejections: AtomicU64,
    retries: AtomicU64,
    full_object_fallbacks: AtomicU64,
    seed_buckets_read: AtomicU64,
    sequence_blocks_read: AtomicU64,
}

impl MetricsCounter {
    pub fn metadata_request(&self) {
        self.metadata_requests.fetch_add(1, Ordering::Relaxed);
    }

    pub fn head_request(&self) {
        self.head_requests.fetch_add(1, Ordering::Relaxed);
    }

    pub fn get_request(&self) {
        self.get_requests.fetch_add(1, Ordering::Relaxed);
    }

    pub fn range_request(&self) {
        self.range_requests.fetch_add(1, Ordering::Relaxed);
    }

    pub fn stream_request(&self) {
        self.stream_requests.fetch_add(1, Ordering::Relaxed);
    }

    pub fn requested_bytes(&self, bytes: u64) {
        self.requested_bytes.fetch_add(bytes, Ordering::Relaxed);
    }

    pub fn returned_bytes(&self, bytes: u64) {
        self.returned_bytes.fetch_add(bytes, Ordering::Relaxed);
    }

    pub fn decoded_bytes(&self, bytes: u64) {
        self.decoded_bytes.fetch_add(bytes, Ordering::Relaxed);
    }

    pub fn remote_bytes(&self, bytes: u64) {
        self.remote_bytes.fetch_add(bytes, Ordering::Relaxed);
    }

    pub fn cache_bytes(&self, bytes: u64) {
        self.cache_bytes.fetch_add(bytes, Ordering::Relaxed);
    }

    pub fn cache_hit(&self) {
        self.cache_hits.fetch_add(1, Ordering::Relaxed);
    }

    pub fn cache_miss(&self) {
        self.cache_misses.fetch_add(1, Ordering::Relaxed);
    }

    pub fn cache_eviction(&self) {
        self.cache_evictions.fetch_add(1, Ordering::Relaxed);
    }

    pub fn cache_evictions(&self, count: usize) {
        self.cache_evictions
            .fetch_add(u64::try_from(count).unwrap_or(u64::MAX), Ordering::Relaxed);
    }

    pub fn stale_cache_rejection(&self) {
        self.stale_cache_rejections.fetch_add(1, Ordering::Relaxed);
    }

    pub fn retry(&self) {
        self.retries.fetch_add(1, Ordering::Relaxed);
    }

    pub fn full_object_fallback(&self) {
        self.full_object_fallbacks.fetch_add(1, Ordering::Relaxed);
    }

    pub fn seed_buckets_read(&self, count: u64) {
        self.seed_buckets_read.fetch_add(count, Ordering::Relaxed);
    }

    pub fn sequence_blocks_read(&self, count: u64) {
        self.sequence_blocks_read
            .fetch_add(count, Ordering::Relaxed);
    }

    #[must_use]
    pub fn snapshot(&self) -> ResourceMetrics {
        ResourceMetrics {
            metadata_requests: self.metadata_requests.load(Ordering::Relaxed),
            head_requests: self.head_requests.load(Ordering::Relaxed),
            get_requests: self.get_requests.load(Ordering::Relaxed),
            range_requests: self.range_requests.load(Ordering::Relaxed),
            stream_requests: self.stream_requests.load(Ordering::Relaxed),
            requested_bytes: self.requested_bytes.load(Ordering::Relaxed),
            returned_bytes: self.returned_bytes.load(Ordering::Relaxed),
            decoded_bytes: self.decoded_bytes.load(Ordering::Relaxed),
            remote_bytes: self.remote_bytes.load(Ordering::Relaxed),
            cache_bytes: self.cache_bytes.load(Ordering::Relaxed),
            cache_hits: self.cache_hits.load(Ordering::Relaxed),
            cache_misses: self.cache_misses.load(Ordering::Relaxed),
            cache_evictions: self.cache_evictions.load(Ordering::Relaxed),
            stale_cache_rejections: self.stale_cache_rejections.load(Ordering::Relaxed),
            retries: self.retries.load(Ordering::Relaxed),
            full_object_fallbacks: self.full_object_fallbacks.load(Ordering::Relaxed),
            seed_buckets_read: self.seed_buckets_read.load(Ordering::Relaxed),
            sequence_blocks_read: self.sequence_blocks_read.load(Ordering::Relaxed),
        }
    }
}
