//! The block cache: fixed-size blocks of slow sources in RAM and on local disk, through foyer.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{ReadAt, StorageError, checked_end};
use foyer::{
    BlockEngineConfig, Code, DeviceBuilder, FsDeviceBuilder, HybridCache, HybridCacheBuilder,
    HybridCacheEntry, HybridCachePolicy, PsyncIoEngineConfig, RecoverMode, Spawner,
};
use std::fs::File;
use std::io::{self, Read, Write};
use std::ops::RangeInclusive;
use std::path::PathBuf;
use std::sync::Arc;
use tokio::runtime::Handle;
use tokio::sync::Semaphore;

/// The largest block size, so that a fetched run of blocks stays small.
const MAX_BLOCK_BYTES: u64 = 1 << 20;
/// Most blocks fetched with one request.
const FILL_BLOCKS: u64 = 16;
/// Every block counts at least this many bytes against the RAM budget, so that many tiny final
/// blocks cannot make the cache index grow without bound.
const MIN_CHARGE: usize = 256;

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

/// A cache of fixed-size blocks shared by every wrapped source of one process. Wrap each source
/// with [`BlockCache::wrap`]. A block is known by the SHA-256 digest of its file and its number,
/// so the same file finds its blocks again in every shard and every run.
pub struct BlockCache {
    blocks: HybridCache<BlockKey, Vec<u8>>,
    block_bytes: u64,
    /// The lock of the disk level directory, held while the cache is open.
    _lock: Option<File>,
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
                .with_engine_config(BlockEngineConfig::new(device))
                .with_recover_mode(RecoverMode::Quiet);
        }
        let blocks = runtime.block_on(storage.build())?;
        Ok(Arc::new(Self {
            blocks,
            block_bytes,
            _lock: lock,
            read_ahead: Arc::new(Semaphore::new(config.read_ahead)),
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
    pub fn close(&self) -> Result<(), StorageError> {
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
            consume(start, &self.fill(inner, file, start, end)?);
        }
        Ok(())
    }

    /// The blocks `first..=last` of a source, fetching each run of missing ones with one request.
    fn fill(
        &self,
        inner: &dyn ReadAt,
        file: [u8; 32],
        first: u64,
        last: u64,
    ) -> io::Result<Vec<Block>> {
        let mut found: Vec<Option<Block>> = (first..=last)
            .map(|block| self.find(self.key(file, block)))
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
            self.fetch_run(inner, file, run_first, &mut found[index..run_end])?;
            index = run_end;
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
        let cache = &self.cache;
        let length = length.min(FILL_BLOCKS * cache.block_bytes);
        self.inner.will_need(offset, length);
        let end = offset.saturating_add(length).min(self.inner.size());
        if offset >= end {
            return;
        }
        // Without a free place the hint is dropped: it never changes what is read.
        let Ok(permit) = Arc::clone(&cache.read_ahead).try_acquire_owned() else {
            return;
        };
        let blocks = offset / cache.block_bytes..=(end - 1) / cache.block_bytes;
        let (inner, file, cache) = (Arc::clone(&self.inner), self.file, Arc::clone(cache));
        cache.runtime.clone().spawn_blocking(move || {
            let _permit = permit;
            // A failed hint only means the read fetches the blocks later and reports its error.
            let _ = cache.fill_batches(inner.as_ref(), file, blocks, |_, _| {});
        });
    }
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
