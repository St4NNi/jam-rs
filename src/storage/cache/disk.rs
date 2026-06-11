//! The local disk level of the block cache. One unnamed file of fixed-size slots holds blocks
//! fetched from slow storage; it disappears with the process. The index, the use order and a
//! CRC-32 of every slot stay in memory (about 100 bytes per slot). A read from a slot counts only
//! if the slot still holds the same block afterwards and the bytes match their CRC-32, so reuse
//! of a slot during the read, or a damaged file, reads as a miss and never as other bytes. Written
//! and read slots leave the page cache at once, so this level holds no memory beyond its index.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io;
use std::os::unix::fs::FileExt;
use std::path::Path;
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};

/// A block is named by its object number and its own number in that object.
pub type BlockKey = (u64, u64);

#[derive(Clone, Copy)]
struct Slot {
    slot: u64,
    length: u32,
    crc: u32,
    generation: u64,
    tick: u64,
}

#[derive(Default)]
struct State {
    index: HashMap<BlockKey, Slot>,
    order: BTreeMap<u64, BlockKey>,
    free: Vec<u64>,
    next: u64,
    tick: u64,
    generation: u64,
}

/// Blocks kept in a file on local disk, least recently used first out.
pub struct DiskLevel {
    file: File,
    block_bytes: u64,
    slots: u64,
    state: Mutex<State>,
}

/// Drops a slot from the page cache; a written slot is first written out, as only clean pages go.
fn release(file: &File, offset: u64, length: usize, written: bool) {
    #[cfg(target_os = "linux")]
    {
        use std::os::fd::AsRawFd;
        let (offset, length) = (offset as libc::off_t, length as libc::off_t);
        // SAFETY: the descriptor stays open for the calls; failures only leave pages cached.
        unsafe {
            if written {
                libc::sync_file_range(
                    file.as_raw_fd(),
                    offset,
                    length,
                    libc::SYNC_FILE_RANGE_WAIT_BEFORE
                        | libc::SYNC_FILE_RANGE_WRITE
                        | libc::SYNC_FILE_RANGE_WAIT_AFTER,
                );
            }
            libc::posix_fadvise(file.as_raw_fd(), offset, length, libc::POSIX_FADV_DONTNEED);
        }
    }
    #[cfg(not(target_os = "linux"))]
    let _ = (file, offset, length, written);
}

impl DiskLevel {
    /// Creates an unnamed slot file of at most `budget` bytes in `directory`.
    pub fn create(directory: &Path, budget: u64, block_bytes: u64) -> io::Result<Self> {
        let block_bytes = block_bytes.max(1);
        let file = tempfile::tempfile_in(directory)?;
        let slots = budget / block_bytes;
        file.set_len(slots * block_bytes)?;
        Ok(Self {
            file,
            block_bytes,
            slots,
            state: Mutex::default(),
        })
    }

    /// The size of every slot, which is the block size of the cache that uses this level.
    pub fn block_bytes(&self) -> u64 {
        self.block_bytes
    }

    fn lock(&self) -> MutexGuard<'_, State> {
        self.state.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// The block stored under `key`, or `None` if it is absent or its slot cannot be trusted.
    pub fn get(&self, key: BlockKey) -> Option<Arc<[u8]>> {
        let entry = {
            let mut guard = self.lock();
            let state = &mut *guard;
            state.tick += 1;
            let tick = state.tick;
            let entry = state.index.get_mut(&key)?;
            state.order.remove(&entry.tick);
            entry.tick = tick;
            state.order.insert(tick, key);
            *entry
        };
        let mut data = vec![0; entry.length as usize];
        let offset = entry.slot * self.block_bytes;
        let read = self.file.read_exact_at(&mut data, offset);
        release(&self.file, offset, data.len(), false);
        read.ok()?;
        let unchanged = self
            .lock()
            .index
            .get(&key)
            .is_some_and(|now| now.generation == entry.generation);
        (unchanged && crc32fast::hash(&data) == entry.crc).then(|| Arc::from(data))
    }

    /// Stores a block unless it is held already; the least recently used block makes room.
    /// Returns whether the block was stored.
    pub fn put(&self, key: BlockKey, data: &[u8]) -> io::Result<bool> {
        if self.slots == 0 || data.len() as u64 > self.block_bytes {
            return Ok(false);
        }
        let (slot, generation) = {
            let mut guard = self.lock();
            let state = &mut *guard;
            if state.index.contains_key(&key) {
                return Ok(false);
            }
            let slot = if let Some(slot) = state.free.pop() {
                slot
            } else if state.next < self.slots {
                state.next += 1;
                state.next - 1
            } else {
                let Some((_, oldest)) = state.order.pop_first() else {
                    return Ok(false);
                };
                let evicted = state.index.remove(&oldest);
                let Some(evicted) = evicted else {
                    return Ok(false);
                };
                evicted.slot
            };
            state.generation += 1;
            (slot, state.generation)
        };
        // The slot belongs to nobody while it is written: it is neither indexed nor free.
        let offset = slot * self.block_bytes;
        let written = self.file.write_all_at(data, offset);
        release(&self.file, offset, data.len(), true);
        if let Err(error) = written {
            self.lock().free.push(slot);
            return Err(error);
        }
        let mut guard = self.lock();
        let state = &mut *guard;
        state.tick += 1;
        let tick = state.tick;
        let entry = Slot {
            slot,
            length: data.len() as u32,
            crc: crc32fast::hash(data),
            generation,
            tick,
        };
        if let Some(previous) = state.index.insert(key, entry) {
            state.order.remove(&previous.tick);
            state.free.push(previous.slot);
        }
        state.order.insert(tick, key);
        Ok(true)
    }

    /// A bound for the memory of the slot index: a hash entry, two use-order entries and a free
    /// slot number per slot, doubled for the growth of the tables.
    pub fn memory_bound(&self) -> u64 {
        use std::mem::size_of;
        let per_slot = size_of::<(BlockKey, Slot)>() + 2 * size_of::<(u64, BlockKey)>() + 8;
        2 * (self.slots + 4) * per_slot as u64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn holds_blocks_and_evicts_least_recent() {
        let directory = tempfile::tempdir().unwrap();
        let level = DiskLevel::create(directory.path(), 3 * 64, 64).unwrap();
        let block = |value: u8| vec![value; 64];
        for key in 0..3u64 {
            assert!(level.put((1, key), &block(key as u8)).unwrap());
        }
        assert_eq!(&*level.get((1, 0)).unwrap(), block(0).as_slice());
        level.put((1, 3), &block(3)).unwrap();
        // Block 1 was the least recently used one.
        assert!(level.get((1, 1)).is_none());
        assert_eq!(&*level.get((1, 0)).unwrap(), block(0).as_slice());
        assert_eq!(&*level.get((1, 3)).unwrap(), block(3).as_slice());
        assert!(level.get((2, 0)).is_none());
        assert!(!level.put((1, 3), &block(9)).unwrap());
        assert_eq!(&*level.get((1, 3)).unwrap(), block(3).as_slice());
        assert!(!level.put((1, 4), &[5; 100]).unwrap());
        assert!(level.get((1, 4)).is_none());
        assert!(DiskLevel::create(&directory.path().join("missing"), 64, 64).is_err());
        // The slot file has no name: the directory stays empty.
        assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
    }

    #[test]
    fn slot_indexes_stay_within_the_slot_count() {
        let directory = tempfile::tempdir().unwrap();
        for slots in [1u64, 7, 8, 16, 65] {
            let level = DiskLevel::create(directory.path(), slots * 64, 64).unwrap();
            for key in 0..3 * slots {
                level.put((1, key), &[key as u8; 64]).unwrap();
                let state = level.lock();
                assert!(state.index.len() as u64 <= slots);
                assert_eq!(state.order.len(), state.index.len());
            }
            assert_eq!(
                &*level.get((1, 3 * slots - 1)).unwrap(),
                &[((3 * slots - 1) as u8); 64]
            );
        }
    }

    #[test]
    fn damaged_slot_reads_as_miss() {
        let directory = tempfile::tempdir().unwrap();
        let level = DiskLevel::create(directory.path(), 4 * 64, 64).unwrap();
        level.put((7, 7), &[3u8; 64]).unwrap();
        let slot = level.lock().index[&(7, 7)].slot;
        level.file.write_all_at(&[4u8; 1], slot * 64 + 10).unwrap();
        assert!(level.get((7, 7)).is_none());
    }
}
