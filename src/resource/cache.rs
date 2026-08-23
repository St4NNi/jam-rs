//! Bounded, identity-aware block caching for range resources.

use super::{ByteRange, CacheIdentity, ResourceError, ResourceResult};
use std::collections::{HashMap, VecDeque};
use std::sync::Mutex;

#[derive(Debug)]
struct CacheState {
    identity: Option<CacheIdentity>,
    blocks: HashMap<u64, Vec<u8>>,
    order: VecDeque<u64>,
    bytes: u64,
}

/// A small in-memory block cache.
///
/// Blocks are keyed by their absolute byte offset.  The identity is replaced
/// atomically with the first access to a new resource version, which prevents
/// bytes from two versions of the same URL from being mixed.  Eviction is
/// deterministic FIFO; callers do not rely on cache eviction for correctness.
pub struct BlockCache {
    block_bytes: u64,
    max_bytes: u64,
    state: Mutex<CacheState>,
}

impl std::fmt::Debug for BlockCache {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let state = self.state.lock().map_err(|_| std::fmt::Error)?;
        formatter
            .debug_struct("BlockCache")
            .field("block_bytes", &self.block_bytes)
            .field("max_bytes", &self.max_bytes)
            .field("identity", &state.identity)
            .field("blocks", &state.blocks.len())
            .field("bytes", &state.bytes)
            .finish()
    }
}

impl BlockCache {
    /// Create a cache with a fixed block size and an overall byte limit.
    pub fn new(block_bytes: u64, max_bytes: u64) -> ResourceResult<Self> {
        if block_bytes == 0 {
            return Err(ResourceError::Io {
                locator: "cache".to_string(),
                message: "cache block size must be greater than zero".to_string(),
            });
        }
        Ok(Self {
            block_bytes,
            max_bytes,
            state: Mutex::new(CacheState {
                identity: None,
                blocks: HashMap::new(),
                order: VecDeque::new(),
                bytes: 0,
            }),
        })
    }

    #[must_use]
    pub fn block_bytes(&self) -> u64 {
        self.block_bytes
    }

    #[must_use]
    pub fn max_bytes(&self) -> u64 {
        self.max_bytes
    }

    /// Select an identity and invalidate blocks from a previous version.
    pub fn prepare(&self, identity: &CacheIdentity) -> ResourceResult<()> {
        let mut state = self.state.lock().map_err(|_| ResourceError::Io {
            locator: identity.redacted_locator.clone(),
            message: "cache lock poisoned".to_string(),
        })?;
        if state.identity.as_ref() != Some(identity) {
            state.identity = Some(identity.clone());
            state.blocks.clear();
            state.order.clear();
            state.bytes = 0;
        }
        Ok(())
    }

    /// Return a full block when present for the selected identity.
    pub fn get(&self, identity: &CacheIdentity, offset: u64) -> ResourceResult<Option<Vec<u8>>> {
        let mut state = self.state.lock().map_err(|_| ResourceError::Io {
            locator: identity.redacted_locator.clone(),
            message: "cache lock poisoned".to_string(),
        })?;
        if state.identity.as_ref() != Some(identity) {
            state.identity = Some(identity.clone());
            state.blocks.clear();
            state.order.clear();
            state.bytes = 0;
            return Ok(None);
        }
        Ok(state.blocks.get(&offset).cloned())
    }

    /// Insert a full block, evicting oldest blocks until the byte bound holds.
    pub fn insert(
        &self,
        identity: &CacheIdentity,
        offset: u64,
        bytes: Vec<u8>,
    ) -> ResourceResult<()> {
        let mut state = self.state.lock().map_err(|_| ResourceError::Io {
            locator: identity.redacted_locator.clone(),
            message: "cache lock poisoned".to_string(),
        })?;
        if state.identity.as_ref() != Some(identity) {
            state.identity = Some(identity.clone());
            state.blocks.clear();
            state.order.clear();
            state.bytes = 0;
        }
        let len = u64::try_from(bytes.len()).map_err(|_| ResourceError::Io {
            locator: identity.redacted_locator.clone(),
            message: "cache block length does not fit in u64".to_string(),
        })?;
        if len == 0 || len > self.max_bytes {
            return Ok(());
        }

        if let Some(previous) = state.blocks.remove(&offset) {
            state.bytes = state
                .bytes
                .saturating_sub(u64::try_from(previous.len()).unwrap_or(u64::MAX));
            state.order.retain(|key| *key != offset);
        }
        while state.bytes.saturating_add(len) > self.max_bytes {
            let Some(oldest) = state.order.pop_front() else {
                break;
            };
            if let Some(previous) = state.blocks.remove(&oldest) {
                state.bytes = state
                    .bytes
                    .saturating_sub(u64::try_from(previous.len()).unwrap_or(u64::MAX));
            }
        }
        if state.bytes.saturating_add(len) <= self.max_bytes {
            state.bytes = state.bytes.saturating_add(len);
            state.order.push_back(offset);
            state.blocks.insert(offset, bytes);
        }
        Ok(())
    }

    /// Number of bytes currently held by the cache.
    pub fn bytes(&self) -> ResourceResult<u64> {
        self.state
            .lock()
            .map(|state| state.bytes)
            .map_err(|_| ResourceError::Io {
                locator: "cache".to_string(),
                message: "cache lock poisoned".to_string(),
            })
    }

    /// Number of cached blocks currently held.
    pub fn len(&self) -> ResourceResult<usize> {
        self.state
            .lock()
            .map(|state| state.blocks.len())
            .map_err(|_| ResourceError::Io {
                locator: "cache".to_string(),
                message: "cache lock poisoned".to_string(),
            })
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.state
            .lock()
            .map(|state| state.blocks.is_empty())
            .unwrap_or(true)
    }

    /// Return the blocks needed to cover a range, in offset order.
    pub fn block_offsets(&self, range: ByteRange) -> ResourceResult<Vec<u64>> {
        let end = range.end()?;
        if range.is_empty() {
            return Ok(Vec::new());
        }
        let first = range.offset / self.block_bytes;
        let last = (end - 1) / self.block_bytes;
        let count = last
            .checked_sub(first)
            .and_then(|n| n.checked_add(1))
            .ok_or(ResourceError::RangeOverflow {
                offset: range.offset,
                length: range.length,
            })?;
        let count = usize::try_from(count).map_err(|_| ResourceError::Io {
            locator: "cache".to_string(),
            message: "range covers too many cache blocks".to_string(),
        })?;
        let mut offsets = Vec::with_capacity(count);
        for index in 0..count {
            let index = u64::try_from(index).map_err(|_| ResourceError::Io {
                locator: "cache".to_string(),
                message: "cache block index does not fit in u64".to_string(),
            })?;
            offsets.push((first + index) * self.block_bytes);
        }
        Ok(offsets)
    }
}

impl Default for BlockCache {
    fn default() -> Self {
        Self::new(1024 * 1024, 4 * 1024 * 1024 * 1024)
            .expect("default cache configuration is valid")
    }
}
