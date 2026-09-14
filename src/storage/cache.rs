//! The block cache: fixed-size blocks of slow sources in RAM and on local disk, through foyer.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::counts::count_blocks;
use super::{ReadAt, StorageError, checked_end, lock};
use foyer::{
    BlockEngineConfig, Code, DeviceBuilder, FsDeviceBuilder, HybridCache, HybridCacheBuilder,
    HybridCacheEntry, HybridCachePolicy, PsyncIoEngineConfig, RecoverMode, Spawner,
};
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::fs::File;
use std::io::{self, Read, Write};
use std::ops::RangeInclusive;
use std::path::PathBuf;
use std::sync::{Arc, Condvar, Mutex, PoisonError};
use tokio::runtime::Handle;
use tokio::sync::Semaphore;

/// The largest block size, so that a fetched run of blocks stays small.
const MAX_BLOCK_BYTES: u64 = 1 << 20;
/// Most blocks fetched with one request.
const FILL_BLOCKS: u64 = 16;
/// Every block counts at least this many bytes against the RAM budget, so that many tiny final
/// blocks cannot make the cache index grow without bound.
const MIN_CHARGE: usize = 256;

/// Threads that write blocks to the disk level, one per GiB of it. Each holds a region of the
/// disk level, so a small one gets one writer.
const DISK_WRITERS: u64 = 4;
/// RAM of the write buffers of the disk level, at most a sixteenth of the disk level.
const DISK_BUFFER: u64 = 128 << 20;
/// Bytes of blocks that may wait to be written to the disk level, at most a quarter of it;
/// further blocks are not written. Fetches come in bursts of hundreds of MiB.
const DISK_QUEUE: u64 = 1 << 30;
/// The smallest write buffers and queue, the defaults of foyer.
const DISK_MIN_QUEUE: u64 = 16 << 20;
/// The file in the disk level directory that one process locks while it uses the directory.
const LOCK_FILE: &str = "lock";

/// The local disk level: a directory and the most bytes to keep there. Its blocks stay there
/// when the process ends and are found again by the next one.
#[derive(Clone, Debug)]
pub struct DiskConfig {
    /// The directory of the cache files.
    pub directory: PathBuf,
    /// The most bytes of blocks to keep on disk.
    pub bytes: u64,
}

/// The settings of a block cache.
#[derive(Clone, Debug)]
pub struct CacheConfig {
    /// The RAM budget for blocks in bytes.
    pub ram_bytes: u64,
    /// The size of one block in bytes, at most 1 MiB.
    pub block_bytes: u64,
    /// The optional disk level.
    pub disk: Option<DiskConfig>,
    /// Hinted ranges fetched ahead at once; further hints wait for a place (0: hints only go
    /// to the source).
    pub read_ahead: usize,
}

impl CacheConfig {
    /// A RAM-only cache of `ram_bytes` with blocks of `block_bytes`, without read-ahead.
    pub fn new(ram_bytes: u64, block_bytes: u64) -> Self {
        Self {
            ram_bytes,
            block_bytes,
            disk: None,
            read_ahead: 0,
        }
    }
}

/// A cache of fixed-size blocks shared by every wrapped source of one process. Wrap each source
/// with [`BlockCache::wrap`]. A block is known by the SHA-256 digest of its file and its number,
/// so the same file finds its blocks again in every shard and every run.
pub struct BlockCache {
    blocks: HybridCache<BlockKey, Vec<u8>>,
    block_bytes: u64,
    /// The lock of the disk level directory, held while the cache is open; `None` without a
    /// disk level.
    disk_lock: Option<File>,
    /// Places for hinted ranges fetched ahead, if hints fetch at all.
    read_ahead: Option<Arc<Semaphore>>,
    /// The number of places for hinted ranges.
    places: usize,
    /// Blocks being fetched now, so another read waits for them instead of asking again.
    in_flight: Mutex<HashMap<BlockKey, Arc<Fetch>>>,
    runtime: Handle,
}

impl BlockCache {
    /// Builds the cache; disk reads and read-ahead run on `runtime`. Fails if the block size is
    /// not between 1 byte and 1 MiB or the disk level cannot be created.
    pub fn open(config: &CacheConfig, runtime: Handle) -> Result<Arc<Self>, StorageError> {
        let block_bytes = config.block_bytes;
        if block_bytes == 0 || block_bytes > MAX_BLOCK_BYTES {
            return Err(StorageError::Config(format!(
                "cache block size {block_bytes} is not between 1 and {MAX_BLOCK_BYTES} bytes"
            )));
        }
        let mut storage = HybridCacheBuilder::new()
            .with_policy(HybridCachePolicy::WriteOnInsertion)
            .memory(config.ram_bytes as usize)
            .with_weighter(|_: &BlockKey, block: &Vec<u8>| block.len().max(MIN_CHARGE))
            .storage()
            .with_spawner(Spawner::from(runtime.clone()));
        let mut lock = None;
        if let Some(disk) = &config.disk {
            lock = Some(lock_directory(disk)?);
            let device = FsDeviceBuilder::new(&disk.directory)
                .with_capacity(disk.bytes as usize)
                .build()?;
            storage = storage
                .with_io_engine_config(PsyncIoEngineConfig::new())
                .with_engine_config(
                    BlockEngineConfig::new(device)
                        .with_flushers((disk.bytes >> 30).clamp(1, DISK_WRITERS) as usize)
                        .with_buffer_pool_size(disk_share(disk.bytes / 16, DISK_BUFFER))
                        .with_submit_queue_size_threshold(disk_share(disk.bytes / 4, DISK_QUEUE)),
                )
                .with_recover_mode(RecoverMode::Quiet);
        }
        let blocks = runtime.block_on(storage.build())?;
        Ok(Arc::new(Self {
            blocks,
            block_bytes,
            disk_lock: lock,
            read_ahead: (config.read_ahead > 0)
                .then(|| Arc::new(Semaphore::new(config.read_ahead))),
            places: config.read_ahead,
            in_flight: Mutex::default(),
            runtime,
        }))
    }

    /// Wraps a source whose bytes have the SHA-256 digest `digest` (lower case hex).
    pub fn wrap(
        self: &Arc<Self>,
        inner: Arc<dyn ReadAt>,
        digest: &str,
    ) -> Result<CachedSource, StorageError> {
        let mut file = [0u8; 32];
        hex::decode_to_slice(digest, &mut file)
            .map_err(|_| StorageError::Config(format!("{digest:?} is not a SHA-256 digest")))?;
        Ok(CachedSource {
            inner,
            file,
            cache: Arc::clone(self),
        })
    }

    /// Writes the blocks still on their way to the disk level, so the next run finds them.
    /// Hinted ranges not yet fetched are dropped, and the running fetches end first, so none
    /// runs while the runtime stops.
    pub fn close(&self) -> Result<(), StorageError> {
        if let Some(places) = &self.read_ahead {
            places.close();
            while places.available_permits() < self.places {
                std::thread::sleep(std::time::Duration::from_millis(1));
            }
        }
        Ok(self.runtime.block_on(self.blocks.close())?)
    }

    /// Calls `consume(first, batch)` for the `blocks` of a source, at most
    /// `FILL_BLOCKS` at a time, each batch after the fetch of its missing blocks.
    fn fill_batches(
        &self,
        inner: &dyn ReadAt,
        file: [u8; 32],
        blocks: RangeInclusive<u64>,
        mut consume: impl FnMut(u64, &[Block]),
    ) -> io::Result<()> {
        let last = *blocks.end();
        for start in blocks.step_by(FILL_BLOCKS as usize) {
            let end = last.min(start + FILL_BLOCKS - 1);
            consume(start, &self.fill(inner, file, start, end, true)?);
        }
        Ok(())
    }

    /// The blocks `first..=last` of a source. Each run of missing blocks that no other read
    /// fetches is fetched with one request. Blocks that another read fetches are waited for if
    /// `wait` is set, else left out.
    fn fill(
        &self,
        inner: &dyn ReadAt,
        file: [u8; 32],
        first: u64,
        last: u64,
        wait: bool,
    ) -> io::Result<Vec<Block>> {
        let keys: Vec<BlockKey> = (first..=last).map(|block| self.key(file, block)).collect();
        let mut found: Vec<Option<Block>> = keys.iter().map(|key| self.find(*key)).collect();
        if wait {
            let hits = found.iter().flatten().count() as u64;
            count_blocks(hits, keys.len() as u64 - hits);
        }
        let mut claim = Claim {
            cache: self,
            keys: Vec::new(),
            fetch: Arc::default(),
        };
        let mut others: Vec<Option<Arc<Fetch>>> = vec![None; keys.len()];
        {
            let mut in_flight = lock(&self.in_flight);
            for (at, key) in keys.iter().enumerate() {
                if found[at].is_some() {
                    continue;
                }
                match in_flight.entry(*key) {
                    Entry::Occupied(other) => others[at] = Some(Arc::clone(other.get())),
                    Entry::Vacant(place) => {
                        place.insert(Arc::clone(&claim.fetch));
                        claim.keys.push(*key);
                    }
                }
            }
        }
        let mut index = 0;
        while index < keys.len() {
            if found[index].is_some() || others[index].is_some() {
                index += 1;
                continue;
            }
            let run_end = (index..keys.len())
                .find(|at| found[*at].is_some() || others[*at].is_some())
                .unwrap_or(keys.len());
            self.fetch_run(
                inner,
                file,
                first + index as u64,
                &mut found[index..run_end],
            )?;
            index = run_end;
        }
        drop(claim);
        for (at, other) in others.iter().enumerate() {
            let Some(other) = other.as_ref().filter(|_| wait) else {
                continue;
            };
            other.wait();
            found[at] = self.find(keys[at]);
            if found[at].is_none() {
                // The other fetch failed, or its block already left the cache.
                self.fetch_run(inner, file, first + at as u64, &mut found[at..=at])?;
            }
        }
        Ok(found.into_iter().flatten().collect())
    }

    fn key(&self, file: [u8; 32], block: u64) -> BlockKey {
        BlockKey {
            file,
            block_bytes: self.block_bytes,
            block,
        }
    }

    /// A block held in RAM, or else on the disk level. A failed disk read counts as a miss: the
    /// disk level is best effort. RAM is looked at alone first, as a lookup of the disk level
    /// can join an older one that misses a block inserted since.
    fn find(&self, key: BlockKey) -> Option<Block> {
        if let Some(block) = self.blocks.memory().get(&key) {
            return Some(block);
        }
        self.disk_lock.as_ref()?;
        match self.blocks.get(&key).try_unwrap() {
            Ok(block) => Some(block),
            Err(pending) => self.runtime.block_on(pending.into_future()).ok().flatten(),
        }
    }

    /// Reads the run of missing blocks that starts at block `first` with one request, and fills
    /// `run` with them.
    fn fetch_run(
        &self,
        inner: &dyn ReadAt,
        file: [u8; 32],
        first: u64,
        run: &mut [Option<Block>],
    ) -> io::Result<()> {
        let start = first * self.block_bytes;
        let stop = ((first + run.len() as u64) * self.block_bytes).min(inner.size());
        let mut bytes = vec![0; (stop - start) as usize];
        inner.read_exact_at(start, &mut bytes)?;
        for (at, chunk) in bytes.chunks(self.block_bytes as usize).enumerate() {
            let key = self.key(file, first + at as u64);
            run[at] = Some(self.blocks.insert(key, chunk.to_vec()));
        }
        Ok(())
    }
}

/// A source read through a [`BlockCache`].
pub struct CachedSource {
    inner: Arc<dyn ReadAt>,
    /// The SHA-256 digest of the source's bytes.
    file: [u8; 32],
    cache: Arc<BlockCache>,
}

impl ReadAt for CachedSource {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        let end = checked_end(offset, buffer.len(), self.inner.size())?;
        if buffer.is_empty() {
            return Ok(());
        }
        let block_bytes = self.cache.block_bytes;
        let blocks = offset / block_bytes..=(end - 1) / block_bytes;
        let mut written = 0;
        let inner = self.inner.as_ref();
        self.cache
            .fill_batches(inner, self.file, blocks, |first, batch| {
                for (at, block) in batch.iter().enumerate() {
                    let block_start = (first + at as u64) * block_bytes;
                    let from = (offset.max(block_start) - block_start) as usize;
                    let to = (end.min(block_start + block.len() as u64) - block_start) as usize;
                    buffer[written..written + to - from].copy_from_slice(&block[from..to]);
                    written += to - from;
                }
            })
    }

    fn size(&self) -> u64 {
        self.inner.size()
    }

    fn will_need(&self, offset: u64, length: u64) {
        self.inner.will_need(offset, length);
        let cache = &self.cache;
        let end = offset.saturating_add(length).min(self.inner.size());
        let Some(places) = &cache.read_ahead else {
            return;
        };
        if offset >= end {
            return;
        }
        let last = (end - 1) / cache.block_bytes;
        // Each batch of blocks is fetched on its own, so long ranges use many requests at once.
        for start in (offset / cache.block_bytes..=last).step_by(FILL_BLOCKS as usize) {
            let stop = last.min(start + FILL_BLOCKS - 1);
            let (inner, file, cache) = (Arc::clone(&self.inner), self.file, Arc::clone(cache));
            let places = Arc::clone(places);
            cache.runtime.clone().spawn(async move {
                let Ok(place) = places.acquire_owned().await else {
                    return;
                };
                let fetch = move || {
                    let _place = place;
                    // A failed hint only means the read fetches the blocks and reports the error.
                    let _ = cache.fill(inner.as_ref(), file, start, stop, false);
                };
                let _ = tokio::task::spawn_blocking(fetch).await;
            });
        }
    }
}

/// One read's fetch of some blocks, which other reads of those blocks wait for.
#[derive(Default)]
struct Fetch {
    done: Mutex<bool>,
    finished: Condvar,
}

impl Fetch {
    /// Waits until the fetch ends, with or without its blocks.
    fn wait(&self) {
        let mut done = lock(&self.done);
        while !*done {
            done = (self.finished.wait(done)).unwrap_or_else(PoisonError::into_inner);
        }
    }
}

/// Blocks that one read claimed to fetch. Dropping it, also after an error, releases them and
/// wakes the reads that wait for them.
struct Claim<'a> {
    cache: &'a BlockCache,
    keys: Vec<BlockKey>,
    fetch: Arc<Fetch>,
}

impl Drop for Claim<'_> {
    fn drop(&mut self) {
        let mut in_flight = lock(&self.cache.in_flight);
        for key in &self.keys {
            in_flight.remove(key);
        }
        drop(in_flight);
        *lock(&self.fetch.done) = true;
        self.fetch.finished.notify_all();
    }
}

/// `share` bytes of the disk level, at least the foyer default and at most `most`.
fn disk_share(share: u64, most: u64) -> usize {
    share.clamp(DISK_MIN_QUEUE, most) as usize
}

/// Locks the disk level directory for this process, creating it if needed.
fn lock_directory(disk: &DiskConfig) -> Result<File, StorageError> {
    let failed = |source| StorageError::File {
        path: disk.directory.clone(),
        source,
    };
    std::fs::create_dir_all(&disk.directory).map_err(failed)?;
    let lock = File::create(disk.directory.join(LOCK_FILE)).map_err(failed)?;
    lock.try_lock()
        .map_err(|_| StorageError::Busy(disk.directory.clone()))?;
    Ok(lock)
}

/// A block of a file, known by the file's SHA-256 digest, the block size and its number.
#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
struct BlockKey {
    file: [u8; 32],
    block_bytes: u64,
    block: u64,
}

impl Code for BlockKey {
    fn encode(&self, writer: &mut impl Write) -> foyer::Result<()> {
        writer
            .write_all(&self.file)
            .map_err(foyer::Error::io_error)?;
        writer
            .write_all(&self.block_bytes.to_le_bytes())
            .map_err(foyer::Error::io_error)?;
        writer
            .write_all(&self.block.to_le_bytes())
            .map_err(foyer::Error::io_error)
    }

    fn decode(reader: &mut impl Read) -> foyer::Result<Self> {
        let mut file = [0u8; 32];
        reader
            .read_exact(&mut file)
            .map_err(foyer::Error::io_error)?;
        let mut number = || -> foyer::Result<u64> {
            let mut bytes = [0u8; 8];
            reader
                .read_exact(&mut bytes)
                .map_err(foyer::Error::io_error)?;
            Ok(u64::from_le_bytes(bytes))
        };
        Ok(Self {
            file,
            block_bytes: number()?,
            block: number()?,
        })
    }

    fn estimated_size(&self) -> usize {
        48
    }
}

type Block = HybridCacheEntry<BlockKey, Vec<u8>>;

#[cfg(test)]
#[path = "cache_tests.rs"]
mod tests;
