use super::*;
use crate::test_support::{random_bytes, rng};
use rand::RngExt;
use std::io;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Duration;

fn cache(ram_bytes: u64, block_bytes: u64) -> Arc<BlockCache> {
    BlockCache::open(&CacheConfig::new(ram_bytes, block_bytes)).unwrap()
}

/// Waits until `done` holds, for at most five seconds.
fn wait_for(done: impl Fn() -> bool) {
    for _ in 0..500 {
        if done() {
            return;
        }
        std::thread::sleep(Duration::from_millis(10));
    }
}

/// Bytes in memory that remember their largest read and how many reads ran at once.
struct Watched {
    data: Vec<u8>,
    largest: Arc<AtomicUsize>,
    active: Arc<AtomicUsize>,
    peak: Arc<AtomicUsize>,
    delay: Duration,
}

impl Watched {
    fn new(data: Vec<u8>, delay: Duration) -> Self {
        let counter = || Arc::new(AtomicUsize::new(0));
        Self {
            data,
            largest: counter(),
            active: counter(),
            peak: counter(),
            delay,
        }
    }
}

impl ReadAt for Watched {
    fn size(&self) -> u64 {
        self.data.len() as u64
    }

    fn read_exact_at(&self, at: u64, out: &mut [u8]) -> io::Result<()> {
        self.largest.fetch_max(out.len(), Ordering::SeqCst);
        let now = self.active.fetch_add(1, Ordering::SeqCst) + 1;
        self.peak.fetch_max(now, Ordering::SeqCst);
        std::thread::sleep(self.delay);
        self.active.fetch_sub(1, Ordering::SeqCst);
        self.data.read_exact_at(at, out)
    }
}

#[test]
fn large_range_is_fetched_in_bounded_batches() {
    let data = random_bytes(&mut rng(7), 2 << 20);
    let source = Watched::new(data.clone(), Duration::ZERO);
    let largest = Arc::clone(&source.largest);
    let cache = cache(4096, 4096);
    let cached = cache.wrap(source);
    let mut out = vec![0; data.len()];
    cached.read_exact_at(0, &mut out).unwrap();
    assert_eq!(out, data);
    assert!(largest.load(Ordering::SeqCst) as u64 <= FILL_BLOCKS * 4096);
    assert!(cache.held_bytes() <= 4096);
}

#[test]
fn reads_equal_the_source_and_count_requests() {
    let data = random_bytes(&mut rng(1), 100_000);
    let cache = cache(1 << 20, 4096);
    let source = cache.wrap(data.clone());
    let mut buffer = vec![0; 9000];
    source.read_exact_at(5000, &mut buffer).unwrap();
    assert_eq!(buffer, data[5000..14_000]);
    // Blocks 1, 2 and 3 were missing: one request.
    assert_eq!((cache.stats().requests, cache.stats().misses), (1, 3));
    source.read_exact_at(4096, &mut buffer[..100]).unwrap();
    assert_eq!(buffer[..100], data[4096..4196]);
    assert_eq!((cache.stats().requests, cache.stats().hits), (1, 1));
    let mut draws = rng(7);
    for _ in 0..300 {
        let start = draws.random_range(0..100_000);
        let length = draws.random_range(0..10_000usize).min(100_000 - start);
        let mut out = vec![0; length];
        source.read_exact_at(start as u64, &mut out).unwrap();
        assert_eq!(out, data[start..start + length]);
    }
    assert!(source.read_exact_at(99_999, &mut [0; 2]).is_err());
    assert_eq!(source.size(), 100_000);
    let other = cache.wrap(random_bytes(&mut rng(2), 100_000));
    let mut out = vec![0; 10];
    other.read_exact_at(0, &mut out).unwrap();
    assert_eq!(out, random_bytes(&mut rng(2), 100_000)[..10]);
}

#[test]
fn a_repeated_read_is_a_hit() {
    let cache = cache(1 << 20, 4096);
    let source = cache.wrap(random_bytes(&mut rng(3), 50_000));
    let mut out = vec![0; 3000];
    for _ in 0..2 {
        source.read_exact_at(20_000, &mut out).unwrap();
        source.read_exact_at(40_000, &mut out[..100]).unwrap();
    }
    assert_eq!(cache.stats().requests, 2);
}

#[test]
fn short_final_blocks_stay_within_the_budget() {
    // A budget of two full blocks retains the 10-byte final blocks of a few hundred objects at
    // most: each one is charged 256 bytes.
    let cache = cache(2 * 4096, 4096);
    assert_eq!(cache.shards.len(), 1);
    let mut out = [0u8; 1];
    for object in 0..600 {
        let source = cache.wrap(random_bytes(&mut rng(object), 10));
        source.read_exact_at(0, &mut out).unwrap();
        assert!(cache.held_bytes() <= 2 * 4096);
    }
    assert_eq!(cache.shards[0].lock().unwrap().entries.len(), 32);
}

#[test]
fn budget_evicts_the_least_recently_used_block() {
    let data = random_bytes(&mut rng(3), 40_960);
    let cache = cache(2 * 4096, 4096);
    let source = cache.wrap(data);
    let mut out = vec![0; 1];
    for block in [0u64, 1, 2, 1, 3] {
        source.read_exact_at(block * 4096, &mut out).unwrap();
        assert!(cache.held_bytes() <= 2 * 4096);
    }
    let before = cache.stats();
    // Block 1 was used last but one: still held. Block 0 was evicted first.
    source.read_exact_at(4096, &mut out).unwrap();
    assert_eq!(cache.stats().hits, before.hits + 1);
    source.read_exact_at(0, &mut out).unwrap();
    assert_eq!(cache.stats().misses, before.misses + 1);
    assert!(cache.stats().evictions >= 2);
}

#[test]
fn parallel_misses_reach_the_fetch_limit_and_not_beyond() {
    let readers = 32;
    let data = random_bytes(&mut rng(8), readers * 4096);
    let source = Watched::new(data.clone(), Duration::from_millis(50));
    let peak = Arc::clone(&source.peak);
    let mut config = CacheConfig::new(1 << 20, 4096);
    config.max_fetches = 16;
    let cache = BlockCache::open(&config).unwrap();
    let cached = cache.wrap(source);
    std::thread::scope(|scope| {
        for reader in 0..readers {
            let (cached, data) = (&cached, &data);
            scope.spawn(move || {
                let mut out = vec![0; 4096];
                cached
                    .read_exact_at(reader as u64 * 4096, &mut out)
                    .unwrap();
                assert_eq!(out, data[reader * 4096..(reader + 1) * 4096]);
            });
        }
    });
    // Distinct misses of many threads are fetched together, up to the limit.
    let peak = peak.load(Ordering::SeqCst);
    assert!(peak > 4 && peak <= 16, "peak {peak}");
    assert_eq!(cache.stats().requests, readers as u64);
}

#[test]
fn read_ahead_fills_before_reads() {
    let data = random_bytes(&mut rng(9), 1 << 20);
    let mut config = CacheConfig::new(1 << 20, 4096);
    config.read_ahead_threads = 2;
    config.read_ahead_queue = 64;
    let cache = BlockCache::open(&config).unwrap();
    let source = cache.wrap(data.clone());
    source.will_need(100_000, 50_000);
    let wanted = (150_000u64 - 1) / 4096 - 100_000 / 4096 + 1;
    wait_for(|| cache.stats().misses >= wanted);
    let before = cache.stats();
    assert_eq!(before.misses, wanted);
    let mut out = vec![0; 50_000];
    source.read_exact_at(100_000, &mut out).unwrap();
    assert_eq!(out, data[100_000..150_000]);
    assert_eq!(cache.stats().requests, before.requests);
}

#[test]
fn disk_level_serves_blocks_evicted_from_memory() {
    let data = random_bytes(&mut rng(4), 40_960);
    let directory = tempfile::tempdir().unwrap();
    let mut config = CacheConfig::new(2 * 4096, 4096);
    config.disk = Some(DiskConfig {
        directory: directory.path().into(),
        bytes: 8 * 4096,
    });
    let cache = BlockCache::open(&config).unwrap();
    let source = cache.wrap(data.clone());
    let mut out = vec![0; 4096];
    for block in 0..6u64 {
        source.read_exact_at(block * 4096, &mut out).unwrap();
    }
    wait_for(|| cache.stats().disk_writes == 6);
    let before = cache.stats();
    assert_eq!(
        (before.requests, before.disk_hits, before.disk_writes),
        (6, 0, 6)
    );
    assert_eq!((before.disk_skipped, before.disk_errors), (0, 0));
    // Block 0 left memory long ago but is still on disk: no request to the source.
    source.read_exact_at(100, &mut out).unwrap();
    assert_eq!(out, data[100..4196]);
    let after = cache.stats();
    assert_eq!((after.requests, after.disk_hits), (6, 2));
    assert!(cache.held_bytes() <= 2 * 4096);
}

#[test]
fn shards_share_the_budget() {
    let data = random_bytes(&mut rng(6), 1 << 22);
    let cache = cache(1 << 20, 1024);
    assert_eq!(cache.shards.len(), 16);
    let source = cache.wrap(data.clone());
    read_in_parallel(&source, &data, 4000, 3000);
    assert!(cache.held_bytes() <= 1 << 20);
    assert!(cache.stats().evictions > 0);
}

#[test]
fn parallel_readers_agree() {
    let data = random_bytes(&mut rng(5), 1 << 20);
    let cache = cache(256 << 10, 4096);
    let source = cache.wrap(data.clone());
    read_in_parallel(&source, &data, 2000, 5000);
    assert!(cache.held_bytes() <= 256 << 10);
}

/// Reads `reads` ranges of `length` bytes from `source` on 8 threads and compares them to `data`.
fn read_in_parallel(source: &CachedSource, data: &[u8], reads: u64, length: usize) {
    std::thread::scope(|scope| {
        for thread in 0..8u64 {
            scope.spawn(move || {
                for index in (thread..reads).step_by(8) {
                    let start = (index * 7919) % (data.len() as u64 - length as u64);
                    let mut out = vec![0; length];
                    source.read_exact_at(start, &mut out).unwrap();
                    assert_eq!(out, data[start as usize..start as usize + length]);
                }
            });
        }
    });
}

#[test]
fn a_bad_block_size_is_rejected() {
    assert!(BlockCache::open(&CacheConfig::new(1 << 20, 0)).is_err());
    assert!(BlockCache::open(&CacheConfig::new(1 << 30, 2 << 20)).is_err());
}
