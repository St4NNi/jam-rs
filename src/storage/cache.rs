//! The block cache: fixed-size blocks of slow sources in RAM and on local disk, through foyer.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::counts::count_blocks;
use super::{ReadAt, ReadCost, StorageError, checked_end, lock};
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
use tokio::sync::{OwnedSemaphorePermit, Semaphore};

/// The largest block size, so that a fetched run of blocks stays small.
const MAX_BLOCK_BYTES: u64 = 1 << 20;
/// Most bytes of blocks fetched with one request while most read-ahead places are free, so that
/// more requests fly; else up to the part size of the source's read cost.
const FILL_BYTES: u64 = 1 << 20;
/// Bytes of fetches in flight at once: an eighth of the RAM budget, within these bounds. Each
/// fetch takes room for its bytes before it reads, so many large parts at once stay bounded.
const FETCH_ROOM: std::ops::RangeInclusive<u64> = (32 << 20)..=(1 << 30);
/// Parts of a long read asked for ahead of the first, so its blocks stay in the cache until
/// read.
const PARTS_AHEAD: u64 = 16;
/// Every block counts at least this many bytes against the RAM budget, so that many tiny final
/// blocks cannot make the cache index grow without bound.
const MIN_CHARGE: usize = 256;

/// Threads that write blocks to the disk level, one per GiB of it. Each holds a region of the
/// disk level, so a small one gets one writer.
const DISK_WRITERS: u64 = 4;
/// RAM of the write buffers of the disk level, at most a sixteenth of the disk level.
const DISK_BUFFER: u64 = 128 << 20;
/// The smallest write buffers, the default of foyer.
const DISK_MIN_BUFFER: u64 = 16 << 20;
/// Bytes that one block takes in a write buffer of the disk level at most, besides its data:
/// a header and the key, padded to a page.
const DISK_OVERHEAD: u64 = 4096;
/// Hinted ranges per place that may be fetched or wait for a place; a further hint waits until
/// one of them ends, so hints never pile up without bound.
const HINTS_PER_PLACE: usize = 4;
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
    /// Room in the write buffer of the disk level, in bytes. A block takes room before it is
    /// inserted and gives it back once written, so foyer never drops a write for a full buffer.
    disk_room: Option<(Arc<Semaphore>, u32)>,
    /// Places for hinted ranges fetched ahead, if hints fetch at all.
    read_ahead: Option<Arc<Semaphore>>,
    /// Hinted ranges fetched or waiting for a place; each holds one until it ends.
    hints: Arc<Semaphore>,
    /// Places for hinted ranges fetched ahead.
    places: usize,
    /// Room for the bytes of fetches in flight, and all of it.
    fetch_room: (Arc<Semaphore>, u32),
    /// The most hinted ranges fetched or waiting at once.
    most_hints: usize,
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
        let (mut lock, mut disk_room) = (None, None);
        if let Some(disk) = &config.disk {
            lock = Some(lock_directory(disk)?);
            let device = FsDeviceBuilder::new(&disk.directory)
                .with_capacity(disk.bytes as usize)
                .build()?;
            let writers = (disk.bytes >> 30).clamp(1, DISK_WRITERS);
            let buffers = (disk.bytes / 16).clamp(DISK_MIN_BUFFER, DISK_BUFFER);
            // Each writer has a buffer of its share; all room fits the buffer of one writer.
            let room = u32::try_from(buffers / writers).unwrap_or(u32::MAX);
            disk_room = Some((Arc::new(Semaphore::new(room as usize)), room));
            storage = storage
                .with_io_engine_config(PsyncIoEngineConfig::new())
                .with_engine_config(
                    BlockEngineConfig::new(device)
                        .with_flushers(writers as usize)
                        .with_buffer_pool_size(buffers as usize)
                        .with_submit_queue_size_threshold(buffers as usize),
                )
                .with_recover_mode(RecoverMode::Quiet);
        }
        let blocks = runtime.block_on(storage.build())?;
        Ok(Arc::new(Self {
            blocks,
            block_bytes,
            disk_lock: lock,
            disk_room,
            read_ahead: (config.read_ahead > 0)
                .then(|| Arc::new(Semaphore::new(config.read_ahead))),
            hints: Arc::new(Semaphore::new(config.read_ahead * HINTS_PER_PLACE)),
            places: config.read_ahead,
            fetch_room: {
                let room = Self::fetch_room(config.ram_bytes) as u32;
                (Arc::new(Semaphore::new(room as usize)), room)
            },
            most_hints: config.read_ahead * HINTS_PER_PLACE,
            in_flight: Mutex::default(),
            runtime,
        }))
    }

    /// Wraps a source whose bytes have the SHA-256 digest `digest` (lower case hex), read in
    /// blocks of `block_bytes` (at most 1 MiB), or of the cache's block size.
    pub fn wrap(
        self: &Arc<Self>,
        inner: Arc<dyn ReadAt>,
        digest: &str,
        block_bytes: Option<u64>,
    ) -> Result<CachedSource, StorageError> {
        let mut file = [0u8; 32];
        hex::decode_to_slice(digest, &mut file)
            .map_err(|_| StorageError::Config(format!("{digest:?} is not a SHA-256 digest")))?;
        let block_bytes = block_bytes.unwrap_or(self.block_bytes);
        if block_bytes == 0 || block_bytes > MAX_BLOCK_BYTES {
            return Err(StorageError::Config(format!(
                "cache block size {block_bytes} is not between 1 and {MAX_BLOCK_BYTES} bytes"
            )));
        }
        let part_bytes = inner.read_cost().part_bytes();
        Ok(CachedSource {
            inner,
            file: FileKey {
                digest: file,
                version: 0,
                block_bytes,
                part_bytes,
            },
            cache: Arc::clone(self),
        })
    }

    /// Most blocks of `file` fetched with one request: fewer while most read-ahead places are
    /// free, so that the parts of a long read fly at the same time.
    fn batch(&self, file: FileKey) -> u64 {
        let small = (FILL_BYTES / file.block_bytes).max(1);
        match &self.read_ahead {
            Some(places) if places.available_permits() * 2 > self.places => file.batch().min(small),
            _ => file.batch(),
        }
    }

    /// Bytes of fetches in flight at once for a RAM budget of `ram_bytes`, besides the blocks.
    pub fn fetch_room(ram_bytes: u64) -> u64 {
        (ram_bytes / 8).clamp(*FETCH_ROOM.start(), *FETCH_ROOM.end())
    }

    /// Writes the blocks still on their way to the disk level, so the next run finds them.
    /// Hinted ranges not yet fetched are dropped, and the running fetches end first, so none
    /// runs while the runtime stops.
    pub fn close(&self) -> Result<(), StorageError> {
        if let Some(places) = &self.read_ahead {
            places.close();
            self.hints.close();
            while self.hints.available_permits() < self.most_hints {
                std::thread::sleep(std::time::Duration::from_millis(1));
            }
        }
        Ok(self.runtime.block_on(self.blocks.close())?)
    }

    /// Calls `consume(first, batch)` for the `blocks` of a source, at most one part at a time,
    /// each batch after the fetch of its missing blocks.
    fn fill_batches(
        &self,
        inner: &dyn ReadAt,
        file: FileKey,
        blocks: RangeInclusive<u64>,
        mut consume: impl FnMut(u64, &[Block]),
    ) -> io::Result<()> {
        let (last, batch) = (*blocks.end(), self.batch(file));
        for start in blocks.step_by(batch as usize) {
            let end = last.min(start + batch - 1);
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
        file: FileKey,
        first: u64,
        last: u64,
        wait: bool,
    ) -> io::Result<Vec<Block>> {
        let keys: Vec<BlockKey> = (first..=last).map(|block| file.block(block)).collect();
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
        file: FileKey,
        first: u64,
        run: &mut [Option<Block>],
    ) -> io::Result<()> {
        let start = first * file.block_bytes;
        let stop = ((first + run.len() as u64) * file.block_bytes).min(inner.size());
        let (room, most) = &self.fetch_room;
        let wanted = u32::try_from(stop - start).unwrap_or(u32::MAX).min(*most);
        let _room = self
            .runtime
            .block_on(Arc::clone(room).acquire_many_owned(wanted));
        let mut bytes = vec![0; (stop - start) as usize];
        inner.read_exact_at(start, &mut bytes)?;
        let room = self.disk_room(bytes.len() as u64 + run.len() as u64 * DISK_OVERHEAD);
        for (at, chunk) in bytes.chunks(file.block_bytes as usize).enumerate() {
            let key = file.block(first + at as u64);
            run[at] = Some(self.blocks.insert(key, chunk.to_vec()));
        }
        if let Some(room) = room {
            // The room comes back once foyer has written everything inserted so far.
            let blocks = self.blocks.clone();
            self.runtime.spawn(async move {
                blocks.storage().wait().await;
                drop(room);
            });
        }
        Ok(())
    }

    /// Takes `bytes` of room in the write buffer of the disk level, waiting while it is full.
    /// `None` without a disk level.
    fn disk_room(&self, bytes: u64) -> Option<OwnedSemaphorePermit> {
        let (room, most) = self.disk_room.as_ref()?;
        let bytes = u32::try_from(bytes).unwrap_or(u32::MAX).min(*most);
        let room = Arc::clone(room);
        self.runtime.block_on(room.acquire_many_owned(bytes)).ok()
    }
}

impl BlockCache {
    /// Drops the blocks of the file with SHA-256 digest `digest` and ETag version `version`,
    /// for a file of `size` bytes in blocks of `block_bytes` (the cache's size if `None`).
    pub fn forget(&self, digest: &str, block_bytes: Option<u64>, size: u64, version: u64) {
        let mut file = [0u8; 32];
        if hex::decode_to_slice(digest, &mut file).is_err() {
            return;
        }
        let block_bytes = block_bytes.unwrap_or(self.block_bytes).max(1);
        let file = FileKey {
            digest: file,
            version,
            block_bytes,
            part_bytes: FILL_BYTES,
        };
        for block in 0..size.div_ceil(block_bytes) {
            self.blocks.remove(&file.block(block));
        }
    }
}

/// A source read through a [`BlockCache`].
pub struct CachedSource {
    inner: Arc<dyn ReadAt>,
    file: FileKey,
    cache: Arc<BlockCache>,
}

impl ReadAt for CachedSource {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        let end = checked_end(offset, buffer.len(), self.inner.size())?;
        if buffer.is_empty() {
            return Ok(());
        }
        let file = self.versioned()?;
        let block_bytes = file.block_bytes;
        // The parts after the first are asked for, so they are fetched at the same time.
        let part = self.cache.batch(file) * block_bytes;
        let second = (offset / part + 1) * part;
        if self.cache.read_ahead.is_some() && second < end {
            self.will_need(second, (end - second).min(PARTS_AHEAD * part));
        }
        let blocks = offset / block_bytes..=(end - 1) / block_bytes;
        let mut written = 0;
        let inner = self.inner.as_ref();
        self.cache
            .fill_batches(inner, file, blocks, |first, batch| {
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

    fn version(&self) -> io::Result<u64> {
        self.inner.version()
    }

    /// A cached read that misses still reads the source, so it costs as much as the source's.
    fn read_cost(&self) -> ReadCost {
        self.inner.read_cost()
    }

    fn will_need(&self, offset: u64, length: u64) {
        self.inner.will_need(offset, length);
        let cache = &self.cache;
        let end = offset.saturating_add(length).min(self.inner.size());
        let Some(places) = &cache.read_ahead else {
            return;
        };
        let Ok(file) = self.versioned() else {
            return;
        };
        if offset >= end {
            return;
        }
        let (block_bytes, batch) = (file.block_bytes, cache.batch(file));
        let last = (end - 1) / block_bytes;
        // Each batch of blocks is fetched on its own, so long ranges use many requests at once.
        for start in (offset / block_bytes..=last).step_by(batch as usize) {
            let stop = last.min(start + batch - 1);
            let (inner, cache) = (Arc::clone(&self.inner), Arc::clone(cache));
            let places = Arc::clone(places);
            let hints = Arc::clone(&cache.hints);
            let Ok(hint) = cache.runtime.block_on(hints.acquire_owned()) else {
                return;
            };
            cache.runtime.clone().spawn(async move {
                let Ok(place) = places.acquire_owned().await else {
                    return;
                };
                let fetch = move || {
                    let _hold = (place, hint);
                    // A failed hint only means the read fetches the blocks and reports the error.
                    let _ = cache.fill(inner.as_ref(), file, start, stop, false);
                };
                let _ = tokio::task::spawn_blocking(fetch).await;
            });
        }
    }
}

impl CachedSource {
    /// The key of the file with the version that the source has now.
    fn versioned(&self) -> io::Result<FileKey> {
        let version = self.inner.version()?;
        Ok(FileKey {
            version,
            ..self.file
        })
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

/// A file read through the cache: its SHA-256 digest, the version of its bytes (a number for
/// their ETag), the size of its blocks, and the most bytes fetched with one request.
#[derive(Clone, Copy)]
struct FileKey {
    digest: [u8; 32],
    version: u64,
    block_bytes: u64,
    part_bytes: u64,
}

impl FileKey {
    /// The key of block `block` of the file.
    fn block(self, block: u64) -> BlockKey {
        BlockKey {
            file: self.digest,
            version: self.version,
            block_bytes: self.block_bytes,
            block,
        }
    }

    /// Most blocks fetched with one request.
    fn batch(self) -> u64 {
        (self.part_bytes / self.block_bytes).max(1)
    }
}

/// A block of a file, known by the file's SHA-256 digest, its version, the block size and its
/// number. Blocks of an older version are never found for a newer one.
#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
struct BlockKey {
    file: [u8; 32],
    version: u64,
    block_bytes: u64,
    block: u64,
}

impl Code for BlockKey {
    fn encode(&self, writer: &mut impl Write) -> foyer::Result<()> {
        writer
            .write_all(&self.file)
            .map_err(foyer::Error::io_error)?;
        writer
            .write_all(&self.version.to_le_bytes())
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
            version: number()?,
            block_bytes: number()?,
            block: number()?,
        })
    }

    fn estimated_size(&self) -> usize {
        56
    }
}

type Block = HybridCacheEntry<BlockKey, Vec<u8>>;

#[cfg(test)]
#[path = "cache_tests.rs"]
mod tests;
