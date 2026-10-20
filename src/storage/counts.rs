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
/// Successful GETs below 64 KiB, below 1 MiB, and from 1 MiB.
static S3_BY_SIZE: [AtomicU64; 3] = [AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)];
static S3_IN_FLIGHT: AtomicU64 = AtomicU64::new(0);
static S3_MOST_IN_FLIGHT: AtomicU64 = AtomicU64::new(0);

/// Storage work of the process: all threads, hinted fetches included.
#[derive(Clone, Copy, Debug, Default, Deserialize, PartialEq, Eq, Serialize)]
pub struct StorageCounts {
    /// GET attempts sent to S3: retries and failed attempts included.
    pub s3_requests: u64,
    /// Bytes of the successful GET answers; what failed attempts received is not known.
    pub s3_bytes: u64,
    /// Blocks that reads found in the block cache.
    pub cache_hits: u64,
    /// Blocks that reads did not find in the block cache and fetched or waited for.
    pub cache_misses: u64,
    /// Reads served from pinned ranges.
    pub pinned_reads: u64,
    /// Successful GETs of fewer than 64 KiB.
    pub s3_requests_below_64_kib: u64,
    /// Successful GETs of 64 KiB to less than 1 MiB.
    pub s3_requests_below_1_mib: u64,
    /// Successful GETs of 1 MiB or more.
    pub s3_requests_from_1_mib: u64,
    /// The most GETs in flight at once in the process so far, not only in this part.
    pub s3_most_in_flight: u64,
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
            s3_requests_below_64_kib: load(&S3_BY_SIZE[0]),
            s3_requests_below_1_mib: load(&S3_BY_SIZE[1]),
            s3_requests_from_1_mib: load(&S3_BY_SIZE[2]),
            s3_most_in_flight: load(&S3_MOST_IN_FLIGHT),
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
            s3_requests_below_64_kib: (self.s3_requests_below_64_kib)
                .saturating_sub(earlier.s3_requests_below_64_kib),
            s3_requests_below_1_mib: (self.s3_requests_below_1_mib)
                .saturating_sub(earlier.s3_requests_below_1_mib),
            s3_requests_from_1_mib: (self.s3_requests_from_1_mib)
                .saturating_sub(earlier.s3_requests_from_1_mib),
            s3_most_in_flight: self.s3_most_in_flight,
        }
    }

    /// The share `part` (0 to 1) of these counts, rounded; the most in flight stays.
    pub fn share(&self, part: f64) -> Self {
        let count = |value: u64| (value as f64 * part).round() as u64;
        Self {
            s3_requests: count(self.s3_requests),
            s3_bytes: count(self.s3_bytes),
            cache_hits: count(self.cache_hits),
            cache_misses: count(self.cache_misses),
            pinned_reads: count(self.pinned_reads),
            s3_requests_below_64_kib: count(self.s3_requests_below_64_kib),
            s3_requests_below_1_mib: count(self.s3_requests_below_1_mib),
            s3_requests_from_1_mib: count(self.s3_requests_from_1_mib),
            s3_most_in_flight: self.s3_most_in_flight,
        }
    }

    /// These counts and `other` together; the most in flight is the larger.
    pub fn plus(&self, other: &Self) -> Self {
        Self {
            s3_requests: self.s3_requests + other.s3_requests,
            s3_bytes: self.s3_bytes + other.s3_bytes,
            cache_hits: self.cache_hits + other.cache_hits,
            cache_misses: self.cache_misses + other.cache_misses,
            pinned_reads: self.pinned_reads + other.pinned_reads,
            s3_requests_below_64_kib: self.s3_requests_below_64_kib
                + other.s3_requests_below_64_kib,
            s3_requests_below_1_mib: self.s3_requests_below_1_mib + other.s3_requests_below_1_mib,
            s3_requests_from_1_mib: self.s3_requests_from_1_mib + other.s3_requests_from_1_mib,
            s3_most_in_flight: self.s3_most_in_flight.max(other.s3_most_in_flight),
        }
    }
}

/// Counts one successful GET of `bytes` bytes.
pub(crate) fn count_s3(bytes: u64) {
    S3_REQUESTS.fetch_add(1, Ordering::Relaxed);
    S3_BYTES.fetch_add(bytes, Ordering::Relaxed);
    let size = match bytes {
        0..65_536 => 0,
        65_536..1_048_576 => 1,
        _ => 2,
    };
    S3_BY_SIZE[size].fetch_add(1, Ordering::Relaxed);
}

/// Counts a GET from its start until dropped. At most `limit` are sent at once, so the most in
/// flight counts no more.
pub(crate) fn count_in_flight(limit: usize) -> InFlight {
    let now = S3_IN_FLIGHT.fetch_add(1, Ordering::Relaxed) + 1;
    S3_MOST_IN_FLIGHT.fetch_max(now.min(limit as u64), Ordering::Relaxed);
    InFlight
}

/// One GET in flight, counted down when dropped.
pub(crate) struct InFlight;

impl Drop for InFlight {
    fn drop(&mut self) {
        S3_IN_FLIGHT.fetch_sub(1, Ordering::Relaxed);
    }
}

/// Counts one failed GET attempt, retried or not.
pub(crate) fn count_failed() {
    S3_REQUESTS.fetch_add(1, Ordering::Relaxed);
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
        let flight = count_in_flight(8);
        count_blocks(2, 1);
        count_pinned();
        let after = StorageCounts::now().since(&before);
        // Other tests of this process may count at the same time.
        assert!(after.s3_requests >= 1 && after.s3_bytes >= 100);
        assert!(after.cache_hits >= 2 && after.cache_misses >= 1 && after.pinned_reads >= 1);
        assert!(after.s3_requests_below_64_kib >= 1 && after.s3_most_in_flight >= 1);
        drop(flight);
        let both = after.plus(&after);
        assert_eq!(
            both.s3_requests_below_64_kib,
            2 * after.s3_requests_below_64_kib
        );
        assert_eq!(after.share(1.0), after);
        // The most in flight is of the whole process.
        let none = StorageCounts {
            s3_most_in_flight: before.s3_most_in_flight,
            ..StorageCounts::default()
        };
        assert_eq!(before.since(&StorageCounts::now()), none);
    }
}
