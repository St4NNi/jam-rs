//! The block cache: fixed-size blocks of slow sources in RAM and on local disk, through foyer.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{ReadAt, StorageError, checked_end};
use foyer::{
    BlockEngineConfig, Code, DeviceBuilder, FsDeviceBuilder, HybridCache, HybridCacheBuilder,
    HybridCacheEntry, PsyncIoEngineConfig, RecoverMode, Spawner,
};
use std::io::{self, Read, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use tokio::runtime::Handle;
use tokio::sync::Semaphore;

/// The largest block size, so that a fetched run of blocks stays small.
const MAX_BLOCK_BYTES: u64 = 1 << 20;
/// Most blocks fetched with one request.
const FILL_BLOCKS: u64 = 16;
/// Every block counts at least this many bytes against the RAM budget, so that many tiny final
/// blocks cannot make the cache index grow without bound.
const MIN_CHARGE: usize = 256;

/// The local disk level: a directory and the most bytes to keep there.
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
    /// Hinted ranges fetched ahead at once; further hints are dropped (0: hints only go to the
    /// source).
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

/// A block of one wrapped source.
#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
struct BlockKey {
    /// The number of the wrapped source.
    source: u64,
    /// The number of the block in that source.
    block: u64,
}

impl Code for BlockKey {
    fn encode(&self, writer: &mut impl Write) -> foyer::Result<()> {
        writer
            .write_all(&self.source.to_le_bytes())
            .map_err(foyer::Error::io_error)?;
        writer
            .write_all(&self.block.to_le_bytes())
            .map_err(foyer::Error::io_error)
    }

    fn decode(reader: &mut impl Read) -> foyer::Result<Self> {
        let mut number = || -> foyer::Result<u64> {
            let mut bytes = [0u8; 8];
            reader
                .read_exact(&mut bytes)
                .map_err(foyer::Error::io_error)?;
            Ok(u64::from_le_bytes(bytes))
        };
        Ok(Self {
            source: number()?,
            block: number()?,
        })
    }

    fn estimated_size(&self) -> usize {
        16
    }
}

type Block = HybridCacheEntry<BlockKey, Vec<u8>>;

/// A cache of fixed-size blocks shared by every wrapped source of one process. Wrap each source
/// with [`BlockCache::wrap`].
pub struct BlockCache {
    blocks: HybridCache<BlockKey, Vec<u8>>,
    block_bytes: u64,
    next_source: AtomicU64,
    /// Places for hinted ranges fetched ahead.
    read_ahead: Arc<Semaphore>,
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
            .memory(config.ram_bytes as usize)
            .with_weighter(|_: &BlockKey, block: &Vec<u8>| block.len().max(MIN_CHARGE))
            .storage()
            .with_spawner(Spawner::from(runtime.clone()));
        if let Some(disk) = &config.disk {
            let device = FsDeviceBuilder::new(&disk.directory)
                .with_capacity(disk.bytes as usize)
                .build()?;
            // The disk level starts empty and dies with the process.
            storage = storage
                .with_io_engine_config(PsyncIoEngineConfig::new())
                .with_engine_config(BlockEngineConfig::new(device))
                .with_recover_mode(RecoverMode::None);
        }
        let blocks = runtime.block_on(storage.build())?;
        Ok(Arc::new(Self {
            blocks,
            block_bytes,
            next_source: AtomicU64::new(0),
            read_ahead: Arc::new(Semaphore::new(config.read_ahead)),
            runtime,
        }))
    }

    /// Wraps a source; each wrapped source gets its own key space in the cache.
    pub fn wrap(self: &Arc<Self>, inner: impl ReadAt + 'static) -> CachedSource {
        CachedSource {
            inner: Arc::new(inner),
            source: self.next_source.fetch_add(1, Ordering::Relaxed),
            cache: Arc::clone(self),
        }
    }

    /// Calls `consume(first, blocks)` for the blocks `first..=last` of a source, at most
    /// `FILL_BLOCKS` at a time, each batch after the fetch of its missing blocks.
    fn fill_batches(
        &self,
        inner: &dyn ReadAt,
        source: u64,
        (first, last): (u64, u64),
        mut consume: impl FnMut(u64, &[Block]),
    ) -> io::Result<()> {
        for start in (first..=last).step_by(FILL_BLOCKS as usize) {
            let end = last.min(start + FILL_BLOCKS - 1);
            consume(start, &self.fill(inner, source, start, end)?);
        }
        Ok(())
    }

    /// The blocks `first..=last` of a source, fetching each run of missing ones with one request.
    fn fill(
        &self,
        inner: &dyn ReadAt,
        source: u64,
        first: u64,
        last: u64,
    ) -> io::Result<Vec<Block>> {
        let mut found: Vec<Option<Block>> = (first..=last)
            .map(|block| self.find(BlockKey { source, block }))
            .collect();
        let mut index = 0;
        while index < found.len() {
            if found[index].is_some() {
                index += 1;
                continue;
            }
            let run_end = (index..found.len())
                .find(|at| found[*at].is_some())
                .unwrap_or(found.len());
            let run_first = first + index as u64;
            self.fetch_run(inner, source, run_first, &mut found[index..run_end])?;
            index = run_end;
        }
        Ok(found.into_iter().flatten().collect())
    }

    /// A block held in RAM, or else on the disk level. A failed disk read counts as a miss: the
    /// disk level is best effort.
    fn find(&self, key: BlockKey) -> Option<Block> {
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
        source: u64,
        first: u64,
        run: &mut [Option<Block>],
    ) -> io::Result<()> {
        let start = first * self.block_bytes;
        let stop = ((first + run.len() as u64) * self.block_bytes).min(inner.size());
        let mut bytes = vec![0; (stop - start) as usize];
        inner.read_exact_at(start, &mut bytes)?;
        for (at, chunk) in bytes.chunks(self.block_bytes as usize).enumerate() {
            let key = BlockKey {
                source,
                block: first + at as u64,
            };
            run[at] = Some(self.blocks.insert(key, chunk.to_vec()));
        }
        Ok(())
    }
}

/// A source read through a [`BlockCache`].
pub struct CachedSource {
    inner: Arc<dyn ReadAt>,
    source: u64,
    cache: Arc<BlockCache>,
}

impl ReadAt for CachedSource {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        let end = checked_end(offset, buffer.len(), self.inner.size())?;
        if buffer.is_empty() {
            return Ok(());
        }
        let block_bytes = self.cache.block_bytes;
        let blocks = (offset / block_bytes, (end - 1) / block_bytes);
        let mut written = 0;
        let inner = self.inner.as_ref();
        self.cache
            .fill_batches(inner, self.source, blocks, |first, batch| {
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
        let cache = &self.cache;
        let length = length.min(FILL_BLOCKS * cache.block_bytes);
        self.inner.will_need(offset, length);
        let end = offset.saturating_add(length).min(self.inner.size());
        // Without a free place the hint is dropped: it never changes what is read.
        let Ok(permit) = Arc::clone(&cache.read_ahead).try_acquire_owned() else {
            return;
        };
        if offset >= end {
            return;
        }
        let blocks = (offset / cache.block_bytes, (end - 1) / cache.block_bytes);
        let (inner, source, cache) = (Arc::clone(&self.inner), self.source, Arc::clone(cache));
        cache.runtime.clone().spawn_blocking(move || {
            let _permit = permit;
            // A failed hint only means the read fetches the blocks later and reports its error.
            let _ = cache.fill_batches(inner.as_ref(), source, blocks, |_, _| {});
        });
    }
}

#[cfg(test)]
#[path = "blocks_tests.rs"]
mod tests;
