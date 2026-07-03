//! Tests of the block cache: exact bytes, merged runs, bounded batches, read-ahead, disk level.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::test_support::{random_bytes, rng, runtime};
use rand::RngExt;
use std::sync::atomic::AtomicUsize;
use std::time::Duration;

fn cache(ram_bytes: u64, block_bytes: u64) -> Arc<BlockCache> {
    BlockCache::open(&CacheConfig::new(ram_bytes, block_bytes), runtime()).unwrap()
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

/// Bytes in memory that count their reads and remember the largest one.
struct Watched {
    data: Vec<u8>,
    reads: Arc<AtomicUsize>,
    largest: Arc<AtomicUsize>,
}

impl Watched {
    fn new(data: Vec<u8>) -> Self {
        let reads = Arc::new(AtomicUsize::new(0));
        let largest = Arc::new(AtomicUsize::new(0));
        Self {
            data,
            reads,
            largest,
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
        self.data.read_exact_at(at, out)
    }
}

#[test]
fn batches_stay_bounded() {
    let data = random_bytes(&mut rng(7), 2 << 20);
    let source = Watched::new(data.clone());
    let largest = Arc::clone(&source.largest);
    let cached = cache(4096, 4096).wrap(source);
    let mut out = vec![0; data.len()];
    cached.read_exact_at(0, &mut out).unwrap();
    assert_eq!(out, data);
    assert!(largest.load(Ordering::SeqCst) as u64 <= FILL_BLOCKS * 4096);
}

#[test]
fn reads_match_source() {
    let data = random_bytes(&mut rng(1), 100_000);
    let cache = cache(1 << 20, 4096);
    let source = Watched::new(data.clone());
    let reads = Arc::clone(&source.reads);
    let cached = cache.wrap(source);
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
    let other = cache.wrap(other_data.clone());
    let mut out = vec![0; 10];
    other.read_exact_at(0, &mut out).unwrap();
    assert_eq!(out, other_data[..10]);
}

#[test]
fn repeated_read_hits() {
    let source = Watched::new(random_bytes(&mut rng(3), 50_000));
    let reads = Arc::clone(&source.reads);
    let cached = cache(1 << 20, 4096).wrap(source);
    let mut out = vec![0; 3000];
    for _ in 0..2 {
        cached.read_exact_at(20_000, &mut out).unwrap();
        cached.read_exact_at(40_000, &mut out[..100]).unwrap();
    }
    assert_eq!(reads.load(Ordering::SeqCst), 2);
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
    let cached = BlockCache::open(&config, runtime()).unwrap().wrap(source);
    cached.will_need(100_000, 50_000);
    wait_for(|| reads.load(Ordering::SeqCst) == 1);
    let mut out = vec![0; 50_000];
    cached.read_exact_at(100_000, &mut out).unwrap();
    assert_eq!(out, data[100_000..150_000]);
    assert_eq!(reads.load(Ordering::SeqCst), 1);
}

#[test]
fn disk_serves_evicted() {
    let data = random_bytes(&mut rng(4), 40_960);
    let directory = tempfile::tempdir().unwrap();
    let disk = DiskConfig {
        directory: directory.path().into(),
        bytes: 64 << 20,
    };
    let config = CacheConfig {
        disk: Some(disk),
        ..CacheConfig::new(2 * 4096, 4096)
    };
    let cache = BlockCache::open(&config, runtime()).unwrap();
    let source = Watched::new(data.clone());
    let reads = Arc::clone(&source.reads);
    let cached = cache.wrap(source);
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
fn parallel_readers_agree() {
    let data = random_bytes(&mut rng(5), 1 << 20);
    let cached = cache(256 << 10, 4096).wrap(data.clone());
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
