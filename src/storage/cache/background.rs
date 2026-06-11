//! The threads of a block cache: one writes fetched blocks to the disk level, others fetch
//! hinted ranges before they are read. Both end when the cache is gone.

use super::{BlockCache, BlockKey};
use crate::storage::ReadAt;
use std::sync::atomic::Ordering;
use std::sync::mpsc::Receiver;
use std::sync::{Arc, Mutex, Weak};

/// A fetched block waiting for the disk level.
pub(super) type DiskWrite = (BlockKey, Arc<[u8]>);

/// A hinted range of one wrapped source.
pub(super) struct Hint {
    pub source: Arc<dyn ReadAt>,
    pub object: u64,
    pub offset: u64,
    pub length: u64,
}

/// Starts the disk writer (if the cache has a disk level) and `threads` read-ahead threads.
pub(super) fn start(
    cache: &Arc<BlockCache>,
    threads: usize,
    hints: Receiver<Hint>,
    writes: Receiver<DiskWrite>,
) {
    let hints = Arc::new(Mutex::new(hints));
    for _ in 0..threads {
        let (cache, hints) = (Arc::downgrade(cache), Arc::clone(&hints));
        std::thread::spawn(move || read_ahead(&cache, &hints));
    }
    if cache.disk.is_some() {
        let cache = Arc::downgrade(cache);
        std::thread::spawn(move || write_behind(&cache, &writes));
    }
}

/// Writes queued blocks to the disk level until the cache is gone.
fn write_behind(cache: &Weak<BlockCache>, receiver: &Receiver<DiskWrite>) {
    while let Ok((key, data)) = receiver.recv() {
        let Some(cache) = cache.upgrade() else {
            return;
        };
        let Some(disk) = &cache.disk else {
            return;
        };
        let counter = match disk.put(key, &data) {
            Ok(true) => &cache.counters.disk_writes,
            Ok(false) => continue,
            Err(_) => &cache.counters.disk_errors,
        };
        counter.fetch_add(1, Ordering::Relaxed);
    }
}

/// Fetches hinted ranges until the cache is gone.
fn read_ahead(cache: &Weak<BlockCache>, receiver: &Mutex<Receiver<Hint>>) {
    loop {
        let hint = match receiver.lock() {
            Ok(receiver) => receiver.recv(),
            Err(_) => return,
        };
        let (Ok(hint), Some(cache)) = (hint, cache.upgrade()) else {
            return;
        };
        let end = hint
            .offset
            .saturating_add(hint.length)
            .min(hint.source.size());
        if hint.length > 0 && hint.offset < end {
            let block_bytes = cache.block_bytes;
            // A failed hint only means the read fetches the blocks later, and reports its error.
            let _ = cache.fill_batches(
                hint.source.as_ref(),
                hint.object,
                hint.offset / block_bytes,
                (end - 1) / block_bytes,
                |_, _| {},
            );
        }
    }
}
