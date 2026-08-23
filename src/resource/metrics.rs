//! Atomic counters used by resource readers.

use super::ResourceMetrics;
use std::sync::atomic::{AtomicU64, Ordering};

/// Thread-safe counters shared by a resource and its cache-backed reads.
#[derive(Debug, Default)]
pub struct MetricsCounter {
    metadata_requests: AtomicU64,
    range_requests: AtomicU64,
    stream_requests: AtomicU64,
    remote_bytes: AtomicU64,
    cache_bytes: AtomicU64,
    cache_hits: AtomicU64,
    retries: AtomicU64,
}

impl MetricsCounter {
    pub fn metadata_request(&self) {
        self.metadata_requests.fetch_add(1, Ordering::Relaxed);
    }

    pub fn range_request(&self) {
        self.range_requests.fetch_add(1, Ordering::Relaxed);
    }

    pub fn stream_request(&self) {
        self.stream_requests.fetch_add(1, Ordering::Relaxed);
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

    pub fn retry(&self) {
        self.retries.fetch_add(1, Ordering::Relaxed);
    }

    #[must_use]
    pub fn snapshot(&self) -> ResourceMetrics {
        ResourceMetrics {
            metadata_requests: self.metadata_requests.load(Ordering::Relaxed),
            range_requests: self.range_requests.load(Ordering::Relaxed),
            stream_requests: self.stream_requests.load(Ordering::Relaxed),
            remote_bytes: self.remote_bytes.load(Ordering::Relaxed),
            cache_bytes: self.cache_bytes.load(Ordering::Relaxed),
            cache_hits: self.cache_hits.load(Ordering::Relaxed),
            retries: self.retries.load(Ordering::Relaxed),
        }
    }
}
