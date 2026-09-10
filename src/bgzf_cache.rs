use std::sync::{Condvar, Mutex};

pub const MAX_BGZF_BLOCK_BYTES: usize = 64 * 1024;
pub const DEFAULT_BATCH_BGZF_CACHE_BYTES: usize = 32 * 1024 * 1024;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct BgzfBlockIdentity {
    pub source_sha256: [u8; 32],
    pub source_bytes: u64,
    pub locator_sha256: [u8; 32],
    pub local_file: Option<[u64; 7]>,
    pub compressed_offset: u64,
    pub uncompressed_offset: u64,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct BgzfCacheStats {
    pub capacity_bytes: usize,
    pub accounted_bytes: usize,
    pub payload_bytes: usize,
    pub reserved_bytes: usize,
    pub retained_entry_bytes: usize,
    pub resident_blocks: usize,
    pub loading_blocks: usize,
    pub hits: u64,
    pub blocks_decoded: u64,
    pub evictions: u64,
    pub waits: u64,
}

pub struct BgzfBlockCache {
    capacity_bytes: usize,
    changed: Condvar,
    inner: Mutex<CacheInner>,
}

#[derive(Default)]
struct CacheInner {
    entries: Vec<CacheEntry>,
    payload_bytes: usize,
    reserved_bytes: usize,
    clock: u64,
    hits: u64,
    blocks_decoded: u64,
    evictions: u64,
    waits: u64,
}

struct CacheEntry {
    identity: BgzfBlockIdentity,
    state: EntryState,
}

enum EntryState {
    Loading,
    Ready { data: Box<[u8]>, last_used: u64 },
}

impl BgzfBlockCache {
    pub fn new(capacity_bytes: usize) -> Option<Self> {
        let entry_capacity = capacity_bytes / (MAX_BGZF_BLOCK_BYTES + size_of::<CacheEntry>());
        if entry_capacity == 0 {
            return None;
        }
        let entries = Vec::with_capacity(entry_capacity);
        if entries.capacity() * size_of::<CacheEntry>() + MAX_BGZF_BLOCK_BYTES > capacity_bytes {
            return None;
        }
        Some(Self {
            capacity_bytes,
            changed: Condvar::new(),
            inner: Mutex::new(CacheInner {
                entries,
                ..CacheInner::default()
            }),
        })
    }

    pub fn stats(&self) -> BgzfCacheStats {
        let inner = self.inner.lock().unwrap_or_else(|error| error.into_inner());
        let retained_entry_bytes = inner.entries.capacity() * size_of::<CacheEntry>();
        BgzfCacheStats {
            capacity_bytes: self.capacity_bytes,
            accounted_bytes: inner.payload_bytes + inner.reserved_bytes + retained_entry_bytes,
            payload_bytes: inner.payload_bytes,
            reserved_bytes: inner.reserved_bytes,
            retained_entry_bytes,
            resident_blocks: inner
                .entries
                .iter()
                .filter(|entry| matches!(entry.state, EntryState::Ready { .. }))
                .count(),
            loading_blocks: inner
                .entries
                .iter()
                .filter(|entry| matches!(entry.state, EntryState::Loading))
                .count(),
            hits: inner.hits,
            blocks_decoded: inner.blocks_decoded,
            evictions: inner.evictions,
            waits: inner.waits,
        }
    }

    pub(crate) fn with_block<E, T>(
        &self,
        identity: BgzfBlockIdentity,
        decode: impl FnOnce() -> Result<Vec<u8>, E>,
        consume: impl FnOnce(&[u8]) -> Result<T, E>,
    ) -> Result<T, E> {
        let mut decode = Some(decode);
        'lookup: loop {
            let mut inner = self.inner.lock().unwrap_or_else(|error| error.into_inner());
            if let Some(index) = inner
                .entries
                .iter()
                .position(|entry| entry.identity == identity)
            {
                if matches!(inner.entries[index].state, EntryState::Loading) {
                    inner.waits = inner.waits.saturating_add(1);
                    inner = self
                        .changed
                        .wait(inner)
                        .unwrap_or_else(|error| error.into_inner());
                    drop(inner);
                    continue 'lookup;
                }
                inner.clock = inner.clock.wrapping_add(1);
                let clock = inner.clock;
                inner.hits = inner.hits.saturating_add(1);
                let EntryState::Ready { data, last_used } = &mut inner.entries[index].state else {
                    unreachable!();
                };
                *last_used = clock;
                return consume(data);
            }

            while inner.entries.len() == inner.entries.capacity()
                || accounted_bytes(&inner)
                    .checked_add(MAX_BGZF_BLOCK_BYTES)
                    .is_none_or(|bytes| bytes > self.capacity_bytes)
            {
                if !evict_oldest(&mut inner) {
                    inner.waits = inner.waits.saturating_add(1);
                    inner = self
                        .changed
                        .wait(inner)
                        .unwrap_or_else(|error| error.into_inner());
                    drop(inner);
                    continue 'lookup;
                }
            }
            inner.reserved_bytes += MAX_BGZF_BLOCK_BYTES;
            inner.entries.push(CacheEntry {
                identity,
                state: EntryState::Loading,
            });
            drop(inner);

            let decoded = decode.take().expect("BGZF decoder called once")();
            let mut inner = self.inner.lock().unwrap_or_else(|error| error.into_inner());
            let index = inner
                .entries
                .iter()
                .position(|entry| entry.identity == identity)
                .expect("loading BGZF block remains present");
            inner.reserved_bytes -= MAX_BGZF_BLOCK_BYTES;
            match decoded {
                Ok(data) => {
                    assert!(data.len() <= MAX_BGZF_BLOCK_BYTES);
                    inner.clock = inner.clock.wrapping_add(1);
                    let clock = inner.clock;
                    inner.payload_bytes += data.len();
                    inner.blocks_decoded = inner.blocks_decoded.saturating_add(1);
                    inner.entries[index].state = EntryState::Ready {
                        data: data.into_boxed_slice(),
                        last_used: clock,
                    };
                    self.changed.notify_all();
                    let EntryState::Ready { data, .. } = &inner.entries[index].state else {
                        unreachable!();
                    };
                    return consume(data);
                }
                Err(error) => {
                    inner.entries.swap_remove(index);
                    self.changed.notify_all();
                    return Err(error);
                }
            }
        }
    }
}

fn accounted_bytes(inner: &CacheInner) -> usize {
    inner.payload_bytes + inner.reserved_bytes + inner.entries.capacity() * size_of::<CacheEntry>()
}

fn evict_oldest(inner: &mut CacheInner) -> bool {
    let Some((index, _)) = inner
        .entries
        .iter()
        .enumerate()
        .filter_map(|(index, entry)| match &entry.state {
            EntryState::Ready { last_used, .. } => Some((index, *last_used)),
            EntryState::Loading => None,
        })
        .min_by_key(|(_, last_used)| *last_used)
    else {
        return false;
    };
    let entry = inner.entries.swap_remove(index);
    let EntryState::Ready { data, .. } = entry.state else {
        unreachable!();
    };
    inner.payload_bytes -= data.len();
    inner.evictions = inner.evictions.saturating_add(1);
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{Arc, mpsc};
    use std::thread;

    fn identity(block: u64) -> BgzfBlockIdentity {
        BgzfBlockIdentity {
            source_sha256: [1; 32],
            source_bytes: 100,
            locator_sha256: [2; 32],
            local_file: Some([3; 7]),
            compressed_offset: block,
            uncompressed_offset: block * 10,
        }
    }

    #[test]
    fn concurrent_readers_decode_a_resident_identity_once() {
        let cache = Arc::new(BgzfBlockCache::new(2 * MAX_BGZF_BLOCK_BYTES).unwrap());
        let decodes = Arc::new(AtomicUsize::new(0));
        let (entered_tx, entered_rx) = mpsc::channel();
        let (release_tx, release_rx) = mpsc::channel();

        let first_cache = Arc::clone(&cache);
        let first_decodes = Arc::clone(&decodes);
        let first = thread::spawn(move || {
            first_cache
                .with_block(
                    identity(7),
                    || {
                        first_decodes.fetch_add(1, Ordering::SeqCst);
                        entered_tx.send(()).unwrap();
                        release_rx.recv().unwrap();
                        Ok::<_, ()>(b"ACGT".to_vec())
                    },
                    |data| Ok::<_, ()>(data.to_vec()),
                )
                .unwrap()
        });
        entered_rx.recv().unwrap();

        let second_cache = Arc::clone(&cache);
        let second_decodes = Arc::clone(&decodes);
        let second = thread::spawn(move || {
            second_cache
                .with_block(
                    identity(7),
                    || {
                        second_decodes.fetch_add(1, Ordering::SeqCst);
                        Ok::<_, ()>(b"wrong".to_vec())
                    },
                    |data| Ok::<_, ()>(data.to_vec()),
                )
                .unwrap()
        });
        for _ in 0..10_000 {
            if cache.stats().waits >= 1 {
                break;
            }
            thread::yield_now();
        }
        let loading = cache.stats();
        assert!(loading.waits >= 1);
        assert_eq!(loading.loading_blocks, 1);
        assert_eq!(loading.reserved_bytes, MAX_BGZF_BLOCK_BYTES);
        assert!(loading.accounted_bytes <= loading.capacity_bytes);
        release_tx.send(()).unwrap();

        assert_eq!(first.join().unwrap(), b"ACGT");
        assert_eq!(second.join().unwrap(), b"ACGT");
        assert_eq!(decodes.load(Ordering::SeqCst), 1);
        let stats = cache.stats();
        assert_eq!(stats.blocks_decoded, 1);
        assert_eq!(stats.hits, 1);
        assert!(stats.waits >= 1);
        assert_eq!(stats.payload_bytes, 4);
        assert_eq!(stats.reserved_bytes, 0);
        assert!(stats.accounted_bytes <= stats.capacity_bytes);
    }

    #[test]
    fn capacity_counts_retained_entries_and_evicts_payload() {
        let cache = BgzfBlockCache::new(MAX_BGZF_BLOCK_BYTES + 1024).unwrap();
        for block in 0..3 {
            cache
                .with_block(
                    identity(block),
                    || Ok::<_, ()>(vec![block as u8; MAX_BGZF_BLOCK_BYTES]),
                    |_| Ok::<_, ()>(()),
                )
                .unwrap();
            let stats = cache.stats();
            assert!(stats.accounted_bytes <= stats.capacity_bytes);
            assert_eq!(stats.reserved_bytes, 0);
            assert_eq!(stats.resident_blocks, 1);
        }
        let stats = cache.stats();
        assert_eq!(stats.blocks_decoded, 3);
        assert_eq!(stats.evictions, 2);
        assert!(stats.retained_entry_bytes >= size_of::<CacheEntry>());
    }

    #[test]
    fn failed_decode_releases_reservation_for_retry() {
        let cache = BgzfBlockCache::new(2 * MAX_BGZF_BLOCK_BYTES).unwrap();
        assert_eq!(
            cache.with_block(identity(4), || Err::<Vec<u8>, _>("decode"), |_| Ok(())),
            Err("decode")
        );
        let failed = cache.stats();
        assert_eq!(failed.reserved_bytes, 0);
        assert_eq!(failed.loading_blocks, 0);
        assert_eq!(failed.blocks_decoded, 0);
        assert_eq!(
            cache
                .with_block(
                    identity(4),
                    || Ok::<_, &str>(b"retry".to_vec()),
                    |data| Ok::<_, &str>(data.to_vec()),
                )
                .unwrap(),
            b"retry"
        );
        assert_eq!(cache.stats().blocks_decoded, 1);
    }

    #[test]
    fn every_source_identity_field_separates_blocks() {
        let cache = BgzfBlockCache::new(2 * MAX_BGZF_BLOCK_BYTES).unwrap();
        let original = identity(0);
        let mut identities = vec![original];
        let mut changed = original;
        changed.source_sha256[0] ^= 1;
        identities.push(changed);
        changed = original;
        changed.source_bytes += 1;
        identities.push(changed);
        changed = original;
        changed.locator_sha256[0] ^= 1;
        identities.push(changed);
        changed = original;
        changed.local_file.as_mut().unwrap()[0] += 1;
        identities.push(changed);
        changed = original;
        changed.compressed_offset += 1;
        identities.push(changed);
        changed = original;
        changed.uncompressed_offset += 1;
        identities.push(changed);

        for (value, identity) in identities.into_iter().enumerate() {
            cache
                .with_block(
                    identity,
                    || Ok::<_, ()>(vec![value as u8]),
                    |data| Ok::<_, ()>(data[0]),
                )
                .unwrap();
        }
        assert_eq!(cache.stats().blocks_decoded, 7);
    }
}
