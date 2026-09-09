use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicU64};

pub const OWNER_COUNTER_COUNT: usize = 32;

pub const MANIFEST_BYTES: usize = 0;
pub const OWNERS_OPENED: usize = 1;
pub const REQUESTED_BYTES: usize = 2;
pub const CHECKSUM_ROOT_PAGES_HASHED: usize = 3;
pub const FILE_IDENTITY_CHECKS: usize = 4;
pub const READ_CALLS: usize = 5;
pub const LOOKUP_KEYS: usize = 6;
pub const LOOKUP_PRESENT: usize = 7;
pub const LOOKUP_ABSENT: usize = 8;
pub const OWNER_PARTITIONS_TOUCHED: usize = 9;
pub const BLOCK_SEARCH_COMPARISONS: usize = 10;
pub const DIRECTORY_REQUESTS: usize = 11;
pub const DIRECTORY_BYTES: usize = 12;
pub const HOT_BLOCK_REQUESTS: usize = 13;
pub const HOT_CACHE_HITS: usize = 14;
pub const HOT_CACHE_MISSES: usize = 15;
pub const HOT_ENCODED_BYTES: usize = 16;
pub const HOT_DECODE_CHARGE_BYTES: usize = 17;
pub const MEMBER_RECORDS_DECODED: usize = 18;
pub const COLD_MEMBER_REQUESTS: usize = 19;
pub const COLD_ENCODED_BYTES: usize = 20;
pub const POSITIONS_DECODED: usize = 21;
pub const PAGE_VERIFY_ATTEMPTS: usize = 22;
pub const PAGE_HASHES: usize = 23;
pub const NEWLY_VERIFIED_PAGES: usize = 24;
pub const PAGE_HASH_BYTES: usize = 25;
pub const CHECKSUM_TREE_PAGES_HASHED: usize = 26;
pub const DOCUMENT_RECORD_REQUESTS: usize = 27;
pub const CONTIG_RECORD_REQUESTS: usize = 28;
pub const STRING_REQUESTED_BYTES: usize = 29;
pub const GZI_REQUESTED_BYTES: usize = 30;
pub const RETRIEVAL_BLOCK_REUSES: usize = 31;

const MAX_RETRIEVAL_IDS: usize = 4096;

pub struct OwnerReadObserver {
    enabled: bool,
    counters: [AtomicU64; OWNER_COUNTER_COUNT],
    retrievals: Mutex<Vec<u128>>,
    retrieval_cap: usize,
    capped: AtomicBool,
}

#[derive(Clone, Debug, serde::Serialize)]
pub struct OwnerReadSnapshot {
    pub totals: [u64; OWNER_COUNTER_COUNT],
    pub distinct_retrieval_blocks: usize,
    pub retrieval_ids_capped: bool,
    pub observer_bytes: usize,
}

impl OwnerReadObserver {
    pub fn disabled() -> Self {
        Self {
            enabled: false,
            counters: std::array::from_fn(|_| AtomicU64::new(0)),
            retrievals: Mutex::new(Vec::new()),
            retrieval_cap: 0,
            capped: AtomicBool::new(false),
        }
    }

    pub fn enabled(retrieval_cap: usize) -> Result<Self, &'static str> {
        if !(1..=MAX_RETRIEVAL_IDS).contains(&retrieval_cap) {
            return Err("owner retrieval observer capacity");
        }
        let mut retrievals = Vec::new();
        retrievals
            .try_reserve_exact(retrieval_cap)
            .map_err(|_| "owner retrieval observer allocation")?;
        Ok(Self {
            enabled: true,
            counters: std::array::from_fn(|_| AtomicU64::new(0)),
            retrievals: Mutex::new(retrievals),
            retrieval_cap,
            capped: AtomicBool::new(false),
        })
    }

    pub fn is_enabled(&self) -> bool {
        self.enabled
    }

    pub fn add(&self, counter: usize, value: u64) {
        if self.enabled && value != 0 {
            self.counters[counter].fetch_add(value, std::sync::atomic::Ordering::Relaxed);
        }
    }

    pub fn record_open(
        &self,
        manifest_bytes: u64,
        owners: u64,
        header_bytes: u64,
        root_pages_hashed: u64,
        identity_checks: u64,
    ) {
        for (counter, value) in [
            (MANIFEST_BYTES, manifest_bytes),
            (OWNERS_OPENED, owners),
            (REQUESTED_BYTES, manifest_bytes.saturating_add(header_bytes)),
            (CHECKSUM_ROOT_PAGES_HASHED, root_pages_hashed),
            (FILE_IDENTITY_CHECKS, identity_checks),
            (
                READ_CALLS,
                u64::from(manifest_bytes != 0).saturating_add(owners),
            ),
        ] {
            self.add(counter, value);
        }
    }

    pub fn record_route(&self, keys: u64, present: u64, absent: u64, owners_touched: u64) {
        for (counter, value) in [
            (LOOKUP_KEYS, keys),
            (LOOKUP_PRESENT, present),
            (LOOKUP_ABSENT, absent),
            (OWNER_PARTITIONS_TOUCHED, owners_touched),
        ] {
            self.add(counter, value);
        }
    }

    pub fn record_directory(&self, comparisons: u64, requests: u64, bytes: u64) {
        self.add(BLOCK_SEARCH_COMPARISONS, comparisons);
        self.add(DIRECTORY_REQUESTS, requests);
        self.add(DIRECTORY_BYTES, bytes);
    }

    pub fn record_hot_request(&self, owner: u32, block: u64, encoded_bytes: u64, cache_hit: bool) {
        self.add(HOT_BLOCK_REQUESTS, 1);
        self.add(
            if cache_hit {
                HOT_CACHE_HITS
            } else {
                HOT_CACHE_MISSES
            },
            1,
        );
        self.add(HOT_ENCODED_BYTES, encoded_bytes);
        self.record_retrieval(owner, block);
    }

    pub fn record_hot_decode(&self, charge_bytes: u64, members: u64) {
        self.add(HOT_DECODE_CHARGE_BYTES, charge_bytes);
        self.add(MEMBER_RECORDS_DECODED, members);
    }

    pub fn record_cold(&self, encoded_bytes: u64, positions: u64) {
        self.add(COLD_MEMBER_REQUESTS, 1);
        self.add(COLD_ENCODED_BYTES, encoded_bytes);
        self.add(POSITIONS_DECODED, positions);
    }

    pub fn record_integrity(
        &self,
        verify_attempts: u64,
        page_hashes: u64,
        newly_verified: u64,
        page_hash_bytes: u64,
        tree_pages_hashed: u64,
    ) {
        for (counter, value) in [
            (PAGE_VERIFY_ATTEMPTS, verify_attempts),
            (PAGE_HASHES, page_hashes),
            (NEWLY_VERIFIED_PAGES, newly_verified),
            (PAGE_HASH_BYTES, page_hash_bytes),
            (CHECKSUM_TREE_PAGES_HASHED, tree_pages_hashed),
        ] {
            self.add(counter, value);
        }
    }

    pub fn record_metadata(
        &self,
        document_requests: u64,
        contig_requests: u64,
        string_bytes: u64,
        gzi_bytes: u64,
    ) {
        for (counter, value) in [
            (DOCUMENT_RECORD_REQUESTS, document_requests),
            (CONTIG_RECORD_REQUESTS, contig_requests),
            (STRING_REQUESTED_BYTES, string_bytes),
            (GZI_REQUESTED_BYTES, gzi_bytes),
        ] {
            self.add(counter, value);
        }
    }

    pub fn snapshot(&self) -> OwnerReadSnapshot {
        let retrievals = self
            .retrievals
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        OwnerReadSnapshot {
            totals: std::array::from_fn(|index| {
                self.counters[index].load(std::sync::atomic::Ordering::Relaxed)
            }),
            distinct_retrieval_blocks: retrievals.len(),
            retrieval_ids_capped: self.capped.load(std::sync::atomic::Ordering::Relaxed),
            observer_bytes: std::mem::size_of::<Self>()
                + retrievals.capacity() * std::mem::size_of::<u128>(),
        }
    }

    fn record_retrieval(&self, owner: u32, block: u64) {
        if !self.enabled {
            return;
        }
        let id = (u128::from(owner) << 64) | u128::from(block);
        let mut retrievals = self
            .retrievals
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        match retrievals.binary_search(&id) {
            Ok(_) => self.add(RETRIEVAL_BLOCK_REUSES, 1),
            Err(index) if retrievals.len() < self.retrieval_cap => retrievals.insert(index, id),
            Err(_) => self
                .capped
                .store(true, std::sync::atomic::Ordering::Relaxed),
        }
    }
}

impl Default for OwnerReadObserver {
    fn default() -> Self {
        Self::disabled()
    }
}

impl OwnerReadSnapshot {
    pub fn checked_sub(&self, before: &Self) -> Option<Self> {
        let mut totals = [0; OWNER_COUNTER_COUNT];
        for (index, value) in totals.iter_mut().enumerate() {
            *value = self.totals[index].checked_sub(before.totals[index])?;
        }
        Some(Self {
            totals,
            distinct_retrieval_blocks: self
                .distinct_retrieval_blocks
                .checked_sub(before.distinct_retrieval_blocks)?,
            retrieval_ids_capped: self.retrieval_ids_capped,
            observer_bytes: self.observer_bytes,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn disabled_observer_is_zero_and_fixed_part_is_bounded() {
        let observer = OwnerReadObserver::disabled();
        observer.record_hot_request(1, 2, 100, false);
        let snapshot = observer.snapshot();
        assert_eq!(snapshot.totals, [0; OWNER_COUNTER_COUNT]);
        assert_eq!(snapshot.distinct_retrieval_blocks, 0);
        assert!(std::mem::size_of::<OwnerReadObserver>() <= 4096);
    }

    #[test]
    fn retrieval_table_is_bounded_and_delta_is_checked() {
        let observer = OwnerReadObserver::enabled(2).unwrap();
        let before = observer.snapshot();
        observer.record_hot_request(1, 3, 10, false);
        observer.record_hot_request(1, 3, 10, true);
        observer.record_hot_request(2, 4, 20, false);
        observer.record_hot_request(3, 5, 30, false);
        let after = observer.snapshot();
        let delta = after.checked_sub(&before).unwrap();
        assert_eq!(delta.distinct_retrieval_blocks, 2);
        assert_eq!(delta.totals[RETRIEVAL_BLOCK_REUSES], 1);
        assert!(delta.retrieval_ids_capped);
        assert_eq!(
            after.observer_bytes,
            std::mem::size_of::<OwnerReadObserver>() + 32
        );
    }

    #[test]
    fn retrieval_capacity_rejects_unbounded_values() {
        assert!(OwnerReadObserver::enabled(0).is_err());
        assert!(OwnerReadObserver::enabled(MAX_RETRIEVAL_IDS + 1).is_err());
    }
}
