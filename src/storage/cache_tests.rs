//! Tests of the block cache: exact bytes, merged runs, bounded batches, read-ahead, disk level.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::storage::{PinnedPages, PinnedSource, SliceSource};
use crate::test_support::{random_bytes, rng, runtime};
use rand::RngExt;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Duration;

/// The digests the test sources are cached under.
const DIGEST: &str = "1111111111111111111111111111111111111111111111111111111111111111";
const OTHER: &str = "2222222222222222222222222222222222222222222222222222222222222222";

fn cache(ram_bytes: u64, block_bytes: u64) -> Arc<BlockCache> {
    BlockCache::open(&CacheConfig::new(ram_bytes, block_bytes), runtime()).unwrap()
}

/// Bytes in memory that count their reads, remember the largest one and the most at once, and
/// take `delay` per read.
struct Watched {
    data: Vec<u8>,
    reads: Arc<AtomicUsize>,
    largest: Arc<AtomicUsize>,
    delay: Duration,
    active: AtomicUsize,
    peak: Arc<AtomicUsize>,
}

impl Watched {
    fn new(data: Vec<u8>) -> Self {
        Self {
            data,
            reads: Arc::default(),
            largest: Arc::default(),
            delay: Duration::ZERO,
            active: AtomicUsize::new(0),
            peak: Arc::default(),
        }
    }
}

impl ReadAt for Watched {
    fn size(&self) -> u64 {
        self.data.len() as u64
    }

    fn read_exact_at(&self, at: u64, out: &mut [u8]) -> io::Result<()> {
        self.reads.fetch_add(1, Ordering::SeqCst);
        self.largest.fetch_max(out.len(), Ordering::SeqCst);
        let now = self.active.fetch_add(1, Ordering::SeqCst) + 1;
        self.peak.fetch_max(now, Ordering::SeqCst);
        std::thread::sleep(self.delay);
        self.active.fetch_sub(1, Ordering::SeqCst);
        self.data.read_exact_at(at, out)
    }
}

/// A source whose bytes and version a test replaces.
struct Replaced {
    data: Mutex<Vec<u8>>,
    version: AtomicU64,
}

impl ReadAt for Replaced {
    fn size(&self) -> u64 {
        lock(&self.data).len() as u64
    }

    fn read_exact_at(&self, at: u64, out: &mut [u8]) -> io::Result<()> {
        lock(&self.data).read_exact_at(at, out)
    }

    fn version(&self) -> io::Result<u64> {
        Ok(self.version.load(Ordering::SeqCst))
    }
}

/// Blocks of an old version are not served for a new one, and `forget` drops them.
#[test]
fn versions_stay_apart() {
    let old = random_bytes(&mut rng(3), 20_000);
    let new = random_bytes(&mut rng(4), 20_000);
    let source = Arc::new(Replaced {
        data: Mutex::new(old.clone()),
        version: AtomicU64::new(1),
    });
    let cache = cache(1 << 20, 4096);
    let cached = cache.wrap(source.clone(), DIGEST, None).unwrap();
    let mut out = vec![0; 20_000];
    cached.read_exact_at(0, &mut out).unwrap();
    assert_eq!(out, old);
    *lock(&source.data) = new.clone();
    source.version.store(2, Ordering::SeqCst);
    cached.read_exact_at(0, &mut out).unwrap();
    assert_eq!(out, new);
    // Back at version 1 the old blocks are still there, until they are forgotten.
    source.version.store(1, Ordering::SeqCst);
    cached.read_exact_at(0, &mut out).unwrap();
    assert_eq!(out, old);
    cache.forget(DIGEST, None, 20_000, 1);
    cached.read_exact_at(0, &mut out).unwrap();
    assert_eq!(out, new);
}

#[test]
fn batches_stay_bounded() {
    let data = random_bytes(&mut rng(7), 2 << 20);
    let source = Watched::new(data.clone());
    let largest = Arc::clone(&source.largest);
    let cached = cache(4096, 4096)
        .wrap(Arc::new(source), DIGEST, None)
        .unwrap();
    let mut out = vec![0; data.len()];
    cached.read_exact_at(0, &mut out).unwrap();
    assert_eq!(out, data);
    assert!(largest.load(Ordering::SeqCst) as u64 <= FILL_BYTES);
}

/// Bytes that cost much per read, as on S3.
struct Costly(Vec<u8>);

impl ReadAt for Costly {
    fn size(&self) -> u64 {
        self.0.size()
    }

    fn read_exact_at(&self, at: u64, out: &mut [u8]) -> io::Result<()> {
        self.0.read_exact_at(at, out)
    }

    fn costly_reads(&self) -> bool {
        true
    }
}

/// The cache, a slice and pinned pages tell that reads of a costly source cost much, so close
/// ranges are joined wider; a pinned range costs what its local copy costs.
#[test]
fn costly_reads_forwarded() {
    let source = Arc::new(Costly(vec![0; 10_000]));
    let cached: Arc<dyn ReadAt> =
        Arc::new(cache(1 << 20, 4096).wrap(source, DIGEST, None).unwrap());
    assert!(cached.costly_reads());
    assert!(
        SliceSource::new(Arc::clone(&cached), 0..100)
            .unwrap()
            .costly_reads()
    );
    let local: Arc<dyn ReadAt> = Arc::new(vec![0u8; 100]);
    let range = PinnedSource::new(Arc::clone(&local), 0..100, Arc::clone(&cached));
    assert!(!range.costly_reads());
    assert!(PinnedPages::new(local, vec![0], 100, cached).costly_reads());
}

#[test]
fn reads_match_source() {
    let data = random_bytes(&mut rng(1), 100_000);
    let cache = cache(1 << 20, 4096);
    let source = Watched::new(data.clone());
    let reads = Arc::clone(&source.reads);
    let cached = cache.wrap(Arc::new(source), DIGEST, None).unwrap();
    let mut buffer = vec![0; 9000];
    cached.read_exact_at(5000, &mut buffer).unwrap();
    assert_eq!(buffer, data[5000..14_000]);
    // Blocks 1, 2 and 3 were missing: one request.
    assert_eq!(reads.load(Ordering::SeqCst), 1);
    cached.read_exact_at(4096, &mut buffer[..100]).unwrap();
    assert_eq!(buffer[..100], data[4096..4196]);
    assert_eq!(reads.load(Ordering::SeqCst), 1);
    let mut draws = rng(7);
    for _ in 0..300 {
        let start = draws.random_range(0..100_000);
        let length = draws.random_range(0..10_000usize).min(100_000 - start);
        let mut out = vec![0; length];
        cached.read_exact_at(start as u64, &mut out).unwrap();
        assert_eq!(out, data[start..start + length]);
    }
    assert!(cached.read_exact_at(99_999, &mut [0; 2]).is_err());
    assert_eq!(cached.size(), 100_000);
    let other_data = random_bytes(&mut rng(2), 100_000);
    let other = cache
        .wrap(Arc::new(other_data.clone()), OTHER, None)
        .unwrap();
    let mut out = vec![0; 10];
    other.read_exact_at(0, &mut out).unwrap();
    assert_eq!(out, other_data[..10]);
}

#[test]
fn repeated_read_hits() {
    let source = Watched::new(random_bytes(&mut rng(3), 50_000));
    let reads = Arc::clone(&source.reads);
    let cached = cache(1 << 20, 4096)
        .wrap(Arc::new(source), DIGEST, None)
        .unwrap();
    let mut out = vec![0; 3000];
    for _ in 0..2 {
        cached.read_exact_at(20_000, &mut out).unwrap();
        cached.read_exact_at(40_000, &mut out[..100]).unwrap();
    }
    assert_eq!(reads.load(Ordering::SeqCst), 2);
}

/// Waits until `count` reaches `wanted`, for at most ten seconds.
fn wait_until(count: &AtomicUsize, wanted: usize) {
    let started = std::time::Instant::now();
    while count.load(Ordering::SeqCst) < wanted && started.elapsed() < Duration::from_secs(10) {
        std::thread::sleep(Duration::from_millis(1));
    }
}

#[test]
fn read_ahead_fills() {
    let data = random_bytes(&mut rng(9), 1 << 20);
    let config = CacheConfig {
        read_ahead: 2,
        ..CacheConfig::new(1 << 20, 4096)
    };
    let source = Watched::new(data.clone());
    let reads = Arc::clone(&source.reads);
    let cache = BlockCache::open(&config, runtime()).unwrap();
    let cached = cache.wrap(Arc::new(source), DIGEST, None).unwrap();
    cached.will_need(100_000, 50_000);
    wait_until(&reads, 1);
    let mut out = vec![0; 50_000];
    cached.read_exact_at(100_000, &mut out).unwrap();
    assert_eq!(out, data[100_000..150_000]);
    assert_eq!(reads.load(Ordering::SeqCst), 1);
}

/// A long read in pieces with the next pieces asked for has many fetches in flight, where a
/// plain read has one, and returns the same bytes.
#[test]
fn long_read_ahead() {
    let data = random_bytes(&mut rng(11), 40 << 20);
    let config = CacheConfig {
        read_ahead: 64,
        ..CacheConfig::new(256 << 20, 64 << 10)
    };
    for ahead in [false, true] {
        let mut source = Watched::new(data.clone());
        source.delay = Duration::from_millis(2);
        let peak = Arc::clone(&source.peak);
        let cache = BlockCache::open(&config, runtime()).unwrap();
        let cached = cache.wrap(Arc::new(source), DIGEST, None).unwrap();
        let mut out = vec![0; data.len()];
        match ahead {
            true => crate::storage::read_exact_ahead(&cached, 0, &mut out).unwrap(),
            false => cached.read_exact_at(0, &mut out).unwrap(),
        }
        assert!(out == data);
        let peak = peak.load(Ordering::SeqCst);
        assert!(if ahead { peak > 4 } else { peak == 1 }, "{ahead} {peak}");
    }
}

#[test]
fn hints_stay_bounded() {
    let data = random_bytes(&mut rng(11), 5 << 20);
    let config = CacheConfig {
        read_ahead: 1,
        ..CacheConfig::new(8 << 20, 4096)
    };
    let mut source = Watched::new(data.clone());
    source.delay = Duration::from_millis(20);
    let reads = Arc::clone(&source.reads);
    let cache = BlockCache::open(&config, runtime()).unwrap();
    let cached = cache.wrap(Arc::new(source), DIGEST, None).unwrap();
    // Each hint is one batch; the fifth waits until the first has been fetched.
    for batch in 0..5u64 {
        cached.will_need(batch * FILL_BYTES, 4096);
    }
    assert!(reads.load(Ordering::SeqCst) >= 1);
    wait_until(&reads, 5);
    cache.close().unwrap();
    assert_eq!(cache.hints.available_permits(), HINTS_PER_PLACE);
}

#[test]
fn reads_join_fetches() {
    let data = random_bytes(&mut rng(10), 16 << 20);
    let config = CacheConfig {
        read_ahead: 64,
        // RAM is split over shards of the cache, so it holds far more than the data.
        ..CacheConfig::new(256 << 20, 4096)
    };
    let source = Watched {
        delay: Duration::from_millis(50),
        ..Watched::new(data.clone())
    };
    let (reads, peak) = (Arc::clone(&source.reads), Arc::clone(&source.peak));
    let cache = BlockCache::open(&config, runtime()).unwrap();
    let cached = cache.wrap(Arc::new(source), DIGEST, None).unwrap();
    // A long hint is fetched in batches at once; reads of hinted blocks wait for them.
    cached.will_need(0, 16 << 20);
    wait_until(&reads, 16);
    let mut out = vec![0; 16 << 20];
    cached.read_exact_at(0, &mut out).unwrap();
    assert_eq!(out, data);
    assert_eq!(reads.load(Ordering::SeqCst), 16);
    assert!(peak.load(Ordering::SeqCst) > 1);
}

fn disk_config(directory: &std::path::Path) -> CacheConfig {
    let disk = DiskConfig {
        directory: directory.into(),
        bytes: 64 << 20,
    };
    CacheConfig {
        disk: Some(disk),
        ..CacheConfig::new(2 * 4096, 4096)
    }
}

#[test]
fn disk_serves_evicted() {
    let data = random_bytes(&mut rng(4), 40_960);
    let directory = tempfile::tempdir().unwrap();
    let cache = BlockCache::open(&disk_config(directory.path()), runtime()).unwrap();
    let source = Watched::new(data.clone());
    let reads = Arc::clone(&source.reads);
    let cached = cache.wrap(Arc::new(source), DIGEST, None).unwrap();
    let mut out = vec![0; 4096];
    for block in 0..6u64 {
        cached.read_exact_at(block * 4096, &mut out).unwrap();
    }
    runtime().block_on(cache.blocks.storage().wait());
    // Block 0 left RAM long ago but is on disk: no request to the source.
    cached.read_exact_at(100, &mut out).unwrap();
    assert_eq!(out, data[100..4196]);
    assert_eq!(reads.load(Ordering::SeqCst), 6);
}

#[test]
fn disk_outlives_process() {
    let data = random_bytes(&mut rng(6), 40_960);
    let directory = tempfile::tempdir().unwrap();
    // RAM holds every block, so only writing on insertion puts them on disk.
    let config = CacheConfig {
        ram_bytes: 1 << 20,
        ..disk_config(directory.path())
    };
    let cache = BlockCache::open(&config, runtime()).unwrap();
    assert!(matches!(
        BlockCache::open(&config, runtime()),
        Err(StorageError::Busy(_))
    ));
    let mut out = vec![0; data.len()];
    let cached = cache.wrap(Arc::new(data.clone()), DIGEST, None).unwrap();
    cached.read_exact_at(0, &mut out).unwrap();
    cache.close().unwrap();
    drop((cached, cache));
    // The next open finds every block on disk, under the same digest only.
    let cache = BlockCache::open(&config, runtime()).unwrap();
    let source = Watched::new(data.clone());
    let reads = Arc::clone(&source.reads);
    let cached = cache.wrap(Arc::new(source), DIGEST, None).unwrap();
    cached.read_exact_at(0, &mut out).unwrap();
    assert_eq!(out, data);
    assert_eq!(reads.load(Ordering::SeqCst), 0);
    let other = Watched::new(data.clone());
    let reads = Arc::clone(&other.reads);
    let other = cache.wrap(Arc::new(other), OTHER, None).unwrap();
    other.read_exact_at(0, &mut out[..10]).unwrap();
    assert_eq!(reads.load(Ordering::SeqCst), 1);
    assert!(cache.wrap(Arc::new(data.clone()), "1234", None).is_err());
    assert!(cache.wrap(Arc::new(data.clone()), DIGEST, Some(0)).is_err());
    // Another block size keeps its own blocks.
    let small_blocks = cache
        .wrap(Arc::new(data.clone()), DIGEST, Some(512))
        .unwrap();
    let mut part = vec![0; 1000];
    small_blocks.read_exact_at(700, &mut part).unwrap();
    assert_eq!(part, data[700..1700]);
    cache.close().unwrap();
    drop((cached, other, small_blocks, cache));
    // A smaller disk budget in a later run still reads the right bytes.
    let mut smaller = config.clone();
    smaller.disk.as_mut().unwrap().bytes = 16 << 20;
    let cache = BlockCache::open(&smaller, runtime()).unwrap();
    let cached = cache.wrap(Arc::new(data.clone()), DIGEST, None).unwrap();
    cached.read_exact_at(0, &mut out).unwrap();
    assert_eq!(out, data);
}

#[test]
fn burst_reaches_disk() {
    let data = random_bytes(&mut rng(7), 64 << 20);
    let directory = tempfile::tempdir().unwrap();
    // Far more blocks than the write buffer holds arrive at once from many threads.
    let config = CacheConfig {
        disk: Some(DiskConfig {
            directory: directory.path().into(),
            bytes: 256 << 20,
        }),
        ..CacheConfig::new(1 << 20, 64 << 10)
    };
    let cache = BlockCache::open(&config, runtime()).unwrap();
    let cached = cache.wrap(Arc::new(data.clone()), DIGEST, None).unwrap();
    let part = data.len() / 8;
    std::thread::scope(|scope| {
        for thread in 0..8 {
            let cached = &cached;
            scope.spawn(move || {
                let mut out = vec![0; part];
                cached
                    .read_exact_at((thread * part) as u64, &mut out)
                    .unwrap();
            });
        }
    });
    cache.close().unwrap();
    drop((cached, cache));
    // The next open finds every block on disk.
    let cache = BlockCache::open(&config, runtime()).unwrap();
    let source = Watched::new(data.clone());
    let reads = Arc::clone(&source.reads);
    let cached = cache.wrap(Arc::new(source), DIGEST, None).unwrap();
    let mut out = vec![0; data.len()];
    cached.read_exact_at(0, &mut out).unwrap();
    assert!(out == data);
    assert_eq!(reads.load(Ordering::SeqCst), 0);
}

#[test]
fn parallel_readers_agree() {
    let data = random_bytes(&mut rng(5), 1 << 20);
    let cached = cache(256 << 10, 4096)
        .wrap(Arc::new(data.clone()), DIGEST, None)
        .unwrap();
    std::thread::scope(|scope| {
        for thread in 0..8u64 {
            let (cached, data) = (&cached, &data);
            scope.spawn(move || {
                for index in (thread..2000).step_by(8) {
                    let start = (index * 7919) % (data.len() as u64 - 5000);
                    let mut out = vec![0; 5000];
                    cached.read_exact_at(start, &mut out).unwrap();
                    assert_eq!(out, data[start as usize..start as usize + 5000]);
                }
            });
        }
    });
}

#[test]
fn bad_block_size() {
    let open = |block_bytes| BlockCache::open(&CacheConfig::new(1 << 20, block_bytes), runtime());
    assert!(open(0).is_err());
    assert!(open(2 << 20).is_err());
}
