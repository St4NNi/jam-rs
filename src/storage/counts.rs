//! Counts of the storage work of the whole process, so a search can report what each of its
//! parts asked of S3, the block cache and the pins.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use serde::{Deserialize, Serialize};
use std::sync::atomic::{AtomicU64, Ordering};

static S3_REQUESTS: AtomicU64 = AtomicU64::new(0);
static S3_BYTES: AtomicU64 = AtomicU64::new(0);
static CACHE_HITS: AtomicU64 = AtomicU64::new(0);
static CACHE_MISSES: AtomicU64 = AtomicU64::new(0);
static PINNED_READS: AtomicU64 = AtomicU64::new(0);

/// Storage work of the process: all threads, hinted fetches included.
#[derive(Clone, Copy, Debug, Default, Deserialize, PartialEq, Eq, Serialize)]
pub struct StorageCounts {
    /// GET requests sent to S3.
    pub s3_requests: u64,
    /// Bytes received from S3.
    pub s3_bytes: u64,
    /// Blocks that reads found in the block cache.
    pub cache_hits: u64,
    /// Blocks that reads did not find in the block cache and fetched or waited for.
    pub cache_misses: u64,
    /// Reads served from pinned ranges.
    pub pinned_reads: u64,
}

impl StorageCounts {
    /// The counts so far.
    pub fn now() -> Self {
        let load = |count: &AtomicU64| count.load(Ordering::Relaxed);
        Self {
            s3_requests: load(&S3_REQUESTS),
            s3_bytes: load(&S3_BYTES),
            cache_hits: load(&CACHE_HITS),
            cache_misses: load(&CACHE_MISSES),
            pinned_reads: load(&PINNED_READS),
        }
    }

    /// The counts since `earlier`.
    pub fn since(&self, earlier: &Self) -> Self {
        Self {
            s3_requests: self.s3_requests.saturating_sub(earlier.s3_requests),
            s3_bytes: self.s3_bytes.saturating_sub(earlier.s3_bytes),
            cache_hits: self.cache_hits.saturating_sub(earlier.cache_hits),
            cache_misses: self.cache_misses.saturating_sub(earlier.cache_misses),
            pinned_reads: self.pinned_reads.saturating_sub(earlier.pinned_reads),
        }
    }
}

/// Counts one GET of `bytes` bytes.
pub(crate) fn count_s3(bytes: u64) {
    S3_REQUESTS.fetch_add(1, Ordering::Relaxed);
    S3_BYTES.fetch_add(bytes, Ordering::Relaxed);
}

/// Counts the blocks of one read: found in the cache, and not found.
pub(crate) fn count_blocks(hits: u64, misses: u64) {
    CACHE_HITS.fetch_add(hits, Ordering::Relaxed);
    CACHE_MISSES.fetch_add(misses, Ordering::Relaxed);
}

/// Counts one read served from a pinned range.
pub(crate) fn count_pinned() {
    PINNED_READS.fetch_add(1, Ordering::Relaxed);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_grow_since() {
        let before = StorageCounts::now();
        count_s3(100);
        count_blocks(2, 1);
        count_pinned();
        let after = StorageCounts::now().since(&before);
        // Other tests of this process may count at the same time.
        assert!(after.s3_requests >= 1 && after.s3_bytes >= 100);
        assert!(after.cache_hits >= 2 && after.cache_misses >= 1 && after.pinned_reads >= 1);
        assert_eq!(
            before.since(&StorageCounts::now()),
            StorageCounts::default()
        );
    }
}
