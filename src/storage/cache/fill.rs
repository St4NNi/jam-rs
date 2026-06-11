//! Reading through the block cache: finding blocks in RAM or on disk and fetching the missing
//! ones from the slow source.

use super::background::Hint;
use super::{BlockCache, BlockKey, FILL_BLOCKS};
use crate::storage::{ReadAt, checked_end};
use std::io;
use std::sync::Arc;
use std::sync::atomic::Ordering;

/// A source read through a [`BlockCache`].
pub struct CachedSource {
    inner: Arc<dyn ReadAt>,
    object: u64,
    cache: Arc<BlockCache>,
}

impl CachedSource {
    pub(super) fn new(inner: Arc<dyn ReadAt>, object: u64, cache: Arc<BlockCache>) -> Self {
        Self {
            inner,
            object,
            cache,
        }
    }
}

impl BlockCache {
    /// Calls `consume(first, blocks)` for the blocks `first..=last` of an object, at most
    /// `FILL_BLOCKS` at a time, each batch after the fetch of its missing blocks.
    pub(super) fn fill_batches(
        &self,
        source: &dyn ReadAt,
        object: u64,
        first: u64,
        last: u64,
        mut consume: impl FnMut(u64, &[Arc<[u8]>]),
    ) -> io::Result<()> {
        for start in (first..=last).step_by(FILL_BLOCKS as usize) {
            let _permit = self.fetches.acquire();
            let end = last.min(start + FILL_BLOCKS - 1);
            consume(start, &self.fill(source, object, start, end)?);
        }
        Ok(())
    }

    /// The blocks `first..=last` of an object, fetching each run of missing ones with one request.
    fn fill(
        &self,
        source: &dyn ReadAt,
        object: u64,
        first: u64,
        last: u64,
    ) -> io::Result<Vec<Arc<[u8]>>> {
        let mut found: Vec<Option<Arc<[u8]>>> = (first..=last)
            .map(|block| self.find((object, block)))
            .collect();
        let mut index = 0;
        while index < found.len() {
            if found[index].is_some() {
                self.counters.hits.fetch_add(1, Ordering::Relaxed);
                index += 1;
                continue;
            }
            let run_end = (index..found.len())
                .find(|at| found[*at].is_some())
                .unwrap_or(found.len());
            self.fetch_run(
                source,
                object,
                first + index as u64,
                &mut found[index..run_end],
            )?;
            index = run_end;
        }
        Ok(found.into_iter().flatten().collect())
    }

    /// A block held in RAM, or else on the disk level (which then also puts it in RAM).
    fn find(&self, key: BlockKey) -> Option<Arc<[u8]>> {
        self.get(key).or_else(|| {
            let data = self.disk.as_ref()?.get(key)?;
            self.counters.disk_hits.fetch_add(1, Ordering::Relaxed);
            self.insert(key, Arc::clone(&data));
            Some(data)
        })
    }

    /// Reads the run of missing blocks that starts at block `first` with one request, and fills
    /// `run` with them. Each block is also queued for the disk level.
    fn fetch_run(
        &self,
        source: &dyn ReadAt,
        object: u64,
        first: u64,
        run: &mut [Option<Arc<[u8]>>],
    ) -> io::Result<()> {
        let start = first * self.block_bytes;
        let stop = ((first + run.len() as u64) * self.block_bytes).min(source.size());
        let mut bytes = vec![0; (stop - start) as usize];
        source.read_exact_at(start, &mut bytes)?;
        self.counters.requests.fetch_add(1, Ordering::Relaxed);
        self.counters
            .fetched_bytes
            .fetch_add(bytes.len() as u64, Ordering::Relaxed);
        for (at, chunk) in bytes.chunks(self.block_bytes as usize).enumerate() {
            let data: Arc<[u8]> = Arc::from(chunk);
            let key = (object, first + at as u64);
            self.counters.misses.fetch_add(1, Ordering::Relaxed);
            if let Some(writes) = &self.disk_writes
                && writes.try_send((key, Arc::clone(&data))).is_err()
            {
                self.counters.disk_skipped.fetch_add(1, Ordering::Relaxed);
            }
            self.insert(key, Arc::clone(&data));
            run[at] = Some(data);
        }
        Ok(())
    }
}

impl ReadAt for CachedSource {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        let end = checked_end(offset, buffer.len(), self.inner.size())?;
        if buffer.is_empty() {
            return Ok(());
        }
        let block_bytes = self.cache.block_bytes;
        let (first, last) = (offset / block_bytes, (end - 1) / block_bytes);
        if first == last
            && let Some(data) = self.cache.get((self.object, first))
        {
            let from = (offset - first * block_bytes) as usize;
            buffer.copy_from_slice(&data[from..from + buffer.len()]);
            self.cache.counters.hits.fetch_add(1, Ordering::Relaxed);
            return Ok(());
        }
        let mut written = 0;
        self.cache.fill_batches(
            self.inner.as_ref(),
            self.object,
            first,
            last,
            |first, blocks| {
                for (at, data) in blocks.iter().enumerate() {
                    let block_start = (first + at as u64) * block_bytes;
                    let from = (offset.max(block_start) - block_start) as usize;
                    let to = (end.min(block_start + data.len() as u64) - block_start) as usize;
                    buffer[written..written + to - from].copy_from_slice(&data[from..to]);
                    written += to - from;
                }
            },
        )
    }

    fn size(&self) -> u64 {
        self.inner.size()
    }

    fn will_need(&self, offset: u64, length: u64) {
        let length = length.min(FILL_BLOCKS * self.cache.block_bytes);
        self.inner.will_need(offset, length);
        if let Some(hints) = &self.cache.hints {
            // A full queue drops the hint: it never changes what is read.
            let _ = hints.try_send(Hint {
                source: Arc::clone(&self.inner),
                object: self.object,
                offset,
                length,
            });
        }
    }
}
