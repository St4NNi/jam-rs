//! A bounded cache of fixed-size blocks in front of slow sources, shared by every object of one
//! process. Blocks stay in RAM until the byte budget needs their room, least recently used
//! first. A large cache spreads its blocks over up to 64 shards, each with its own lock, use
//! order and equal part of the budget, so parallel readers rarely wait for each other. Adjacent
//! missing blocks of one read are fetched with one request. With an optional disk level, fetched
//! blocks are also kept on local disk, and blocks missing in RAM are looked up there before the
//! slow source is asked.

mod background;
mod disk;
mod fill;

pub use disk::{BlockKey, DiskLevel};
pub use fill::CachedSource;

use super::limit::InFlight;
use super::{ReadAt, StorageError};
use background::{DiskWrite, Hint};
use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::SyncSender;
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};

/// The default block size, 64 KiB.
pub const DEFAULT_BLOCK_BYTES: u64 = 64 << 10;
/// The largest block size, so that a fetched run of blocks stays small.
const MAX_BLOCK_BYTES: u64 = 1 << 20;
/// Most shards of one cache.
const SHARDS: u64 = 64;
/// Most blocks fetched with one request.
const FILL_BLOCKS: u64 = 16;
/// Fetched blocks waiting for the disk level; with 64 KiB blocks at most 64 MiB.
const DISK_QUEUE_BLOCKS: usize = 1024;
/// Most threads that read hinted ranges ahead.
const MAX_READ_AHEAD_THREADS: usize = 8;
/// Most hints that wait for a read-ahead thread; further hints are dropped.
const MAX_HINT_QUEUE: usize = 256;
/// Every block counts at least this many bytes against the budget, so that many tiny final
/// blocks cannot make the index of blocks grow without bound.
const MIN_CHARGE: u64 = 256;
/// A guess of the bytes the index spends on one block, for `memory_bound`.
const INDEX_BYTES_PER_BLOCK: u64 = 200;

/// The local disk level: a directory and the most bytes to keep there.
#[derive(Clone, Debug)]
pub struct DiskConfig {
    /// The directory of the unnamed slot file.
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
    /// Threads that fetch hinted ranges before they are read (0: hints only go to the source).
    pub read_ahead_threads: usize,
    /// Hints that may wait for those threads.
    pub read_ahead_queue: usize,
    /// Fetches of all threads at once. A fetch from S3 takes about 25 ms, so the workers of a
    /// node (28 or more, most of them waiting for the network) need many fetches at once so that
    /// misses do not queue; 32 gave 1,000 reads a second.
    pub max_fetches: usize,
}

impl CacheConfig {
    /// A RAM-only cache of `ram_bytes` with blocks of `block_bytes`, without read-ahead threads.
    pub fn new(ram_bytes: u64, block_bytes: u64) -> Self {
        Self {
            ram_bytes,
            block_bytes,
            disk: None,
            read_ahead_threads: 0,
            read_ahead_queue: 0,
            max_fetches: 256,
        }
    }
}

/// What a cache did so far.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct CacheStats {
    /// Blocks found in RAM or on disk.
    pub hits: u64,
    /// Blocks fetched from the slow source.
    pub misses: u64,
    /// Reads sent to the slow source; a run of missing blocks is one request.
    pub requests: u64,
    /// Bytes fetched from the slow source.
    pub fetched_bytes: u64,
    /// Blocks dropped from RAM to make room.
    pub evictions: u64,
    /// Blocks found on the local disk level, also counted in `hits`.
    pub disk_hits: u64,
    /// Fetched blocks written to the local disk level.
    pub disk_writes: u64,
    /// Fetched blocks not written to the local disk level because its queue was full.
    pub disk_skipped: u64,
    /// Blocks that could not be written to the local disk level.
    pub disk_errors: u64,
}

#[derive(Default)]
struct Counters {
    hits: AtomicU64,
    misses: AtomicU64,
    requests: AtomicU64,
    fetched_bytes: AtomicU64,
    evictions: AtomicU64,
    disk_hits: AtomicU64,
    disk_writes: AtomicU64,
    disk_skipped: AtomicU64,
    disk_errors: AtomicU64,
}

/// The blocks of one shard with their use order.
#[derive(Default)]
struct Blocks {
    /// Each block with the tick of its last use.
    entries: HashMap<BlockKey, (Arc<[u8]>, u64)>,
    /// The ticks of last use, oldest first, with the block of each.
    order: BTreeMap<u64, BlockKey>,
    charged: u64,
    tick: u64,
}

fn charge(data: &[u8]) -> u64 {
    (data.len() as u64).max(MIN_CHARGE)
}

/// A block cache with an optional disk level. Wrap each object with [`BlockCache::wrap`].
pub struct BlockCache {
    /// Budget of one shard.
    budget: u64,
    block_bytes: u64,
    shards: Vec<Mutex<Blocks>>,
    fetches: InFlight,
    max_fetches: usize,
    next_object: AtomicU64,
    counters: Counters,
    hints: Option<SyncSender<Hint>>,
    disk: Option<DiskLevel>,
    disk_writes: Option<SyncSender<DiskWrite>>,
}

impl BlockCache {
    /// Builds the cache and starts its background threads. Fails if the block size is not
    /// between 1 byte and 1 MiB or the disk level cannot be created.
    pub fn open(config: &CacheConfig) -> Result<Arc<Self>, StorageError> {
        let block_bytes = config.block_bytes;
        if block_bytes == 0 || block_bytes > MAX_BLOCK_BYTES {
            return Err(StorageError::Config(format!(
                "cache block size {block_bytes} is not between 1 and {MAX_BLOCK_BYTES} bytes"
            )));
        }
        let disk = config
            .disk
            .as_ref()
            .map(|disk| DiskLevel::create(&disk.directory, disk.bytes, block_bytes))
            .transpose()?;
        let threads = config.read_ahead_threads.min(MAX_READ_AHEAD_THREADS);
        let queue = config.read_ahead_queue.min(MAX_HINT_QUEUE);
        let (hint_sender, hint_receiver) = std::sync::mpsc::sync_channel(queue);
        let (write_sender, write_receiver) = std::sync::mpsc::sync_channel(DISK_QUEUE_BLOCKS);
        let shards = (config.ram_bytes / block_bytes / SHARDS).clamp(1, SHARDS);
        let cache = Arc::new(Self {
            budget: config.ram_bytes / shards,
            block_bytes,
            shards: (0..shards).map(|_| Mutex::default()).collect(),
            fetches: InFlight::new(config.max_fetches),
            max_fetches: config.max_fetches.max(1),
            next_object: AtomicU64::new(0),
            counters: Counters::default(),
            hints: (threads > 0).then_some(hint_sender),
            disk_writes: disk.is_some().then_some(write_sender),
            disk,
        });
        background::start(&cache, threads, hint_receiver, write_receiver);
        Ok(cache)
    }

    /// Wraps a source; each wrapped source gets its own key space in the cache.
    pub fn wrap(self: &Arc<Self>, inner: impl ReadAt + 'static) -> CachedSource {
        CachedSource::new(
            Arc::new(inner),
            self.next_object.fetch_add(1, Ordering::Relaxed),
            Arc::clone(self),
        )
    }

    /// What the cache did so far.
    pub fn stats(&self) -> CacheStats {
        let counters = &self.counters;
        let load = |counter: &AtomicU64| counter.load(Ordering::Relaxed);
        CacheStats {
            hits: load(&counters.hits),
            misses: load(&counters.misses),
            requests: load(&counters.requests),
            fetched_bytes: load(&counters.fetched_bytes),
            evictions: load(&counters.evictions),
            disk_hits: load(&counters.disk_hits),
            disk_writes: load(&counters.disk_writes),
            disk_skipped: load(&counters.disk_skipped),
            disk_errors: load(&counters.disk_errors),
        }
    }

    /// A bound for the RAM the cache uses: the block budget, the index of the blocks, the blocks
    /// of fetches in progress (a fetched run and the blocks cut from it), the blocks queued for
    /// the disk level and the index of the disk level.
    pub fn memory_bound(&self) -> u64 {
        let blocks = self.budget * self.shards.len() as u64;
        let index = blocks / MIN_CHARGE * INDEX_BYTES_PER_BLOCK;
        let fetching = 2 * self.max_fetches as u64 * FILL_BLOCKS * self.block_bytes;
        let queued = if self.disk.is_some() {
            DISK_QUEUE_BLOCKS as u64 * self.block_bytes
        } else {
            0
        };
        blocks + index + fetching + queued + self.disk.as_ref().map_or(0, DiskLevel::memory_bound)
    }

    /// The bytes charged for the blocks held in RAM now; never above the budget.
    pub fn held_bytes(&self) -> u64 {
        self.shards.iter().map(|shard| lock(shard).charged).sum()
    }

    fn shard(&self, key: BlockKey) -> MutexGuard<'_, Blocks> {
        let mixed =
            (key.0.wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ key.1).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        lock(&self.shards[((mixed >> 32) % self.shards.len() as u64) as usize])
    }

    fn get(&self, key: BlockKey) -> Option<Arc<[u8]>> {
        let mut guard = self.shard(key);
        let blocks = &mut *guard;
        blocks.tick += 1;
        let tick = blocks.tick;
        let (data, last) = blocks.entries.get_mut(&key)?;
        let (data, previous) = (Arc::clone(data), std::mem::replace(last, tick));
        blocks.order.remove(&previous);
        blocks.order.insert(tick, key);
        Some(data)
    }

    fn insert(&self, key: BlockKey, data: Arc<[u8]>) {
        let size = charge(&data);
        if size > self.budget {
            return;
        }
        let mut guard = self.shard(key);
        let blocks = &mut *guard;
        if blocks.entries.contains_key(&key) {
            return;
        }
        while blocks.charged + size > self.budget {
            let Some((_, oldest)) = blocks.order.pop_first() else {
                break;
            };
            if let Some((evicted, _)) = blocks.entries.remove(&oldest) {
                blocks.charged -= charge(&evicted);
                self.counters.evictions.fetch_add(1, Ordering::Relaxed);
            }
        }
        blocks.tick += 1;
        let tick = blocks.tick;
        blocks.charged += size;
        blocks.order.insert(tick, key);
        blocks.entries.insert(key, (data, tick));
    }
}

fn lock<T>(mutex: &Mutex<T>) -> MutexGuard<'_, T> {
    mutex.lock().unwrap_or_else(PoisonError::into_inner)
}
