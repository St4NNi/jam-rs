//! Reading bytes by offset from local files, S3 byte ranges, pinned ranges and the block cache.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

mod cache;
mod checked;
mod counts;
mod error;
mod etags;
mod local;
mod mapped;
mod pinned;
mod readahead;
mod s3;
mod slice;

pub use cache::{BlockCache, CacheConfig, CachedSource, DiskConfig};
pub use checked::Checked;
pub use counts::StorageCounts;
pub use error::StorageError;
pub use etags::{ChangedError, is_changed};
pub use local::LocalFile;
pub use mapped::MappedFile;
pub use pinned::{PinnedPages, PinnedSource};
pub use s3::{ChangeHook, S3Config, S3Limits, S3Object};
pub use slice::SliceSource;

use std::io;
use std::sync::{Mutex, MutexGuard, PoisonError};

/// Positional reads from an immutable file or object, shared by many threads.
pub trait ReadAt: Send + Sync {
    /// Fills `buffer` with the bytes that start at `offset`, or fails.
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()>;

    /// The size of the file or object in bytes.
    fn size(&self) -> u64;

    /// Tells the source that a range will be read soon. It never changes what is read.
    fn will_need(&self, _offset: u64, _length: u64) {}

    /// A number that names the version of the bytes; it differs when the object was replaced.
    /// A source that never changes gives 0.
    fn version(&self) -> io::Result<u64> {
        Ok(0)
    }

    /// What each read costs beyond its bytes, which decides how close ranges are read together.
    fn read_cost(&self) -> ReadCost {
        ReadCost::Local
    }
}

/// What one read of a source costs beyond its bytes. Two ranges are read as one when the bytes
/// between them cost less than one more read; the joined read stops growing at a limit, and the
/// block cache fetches it in parts, several at once.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReadCost {
    /// A local disk: one more read costs about the system time of copying 5 KiB more (measured
    /// on worm).
    Local,
    /// A network file system such as CephFS: one more read costs about as much as 64 KiB more.
    Network,
    /// A store where each read is one request, as S3: about 60 ms per request. With 256
    /// requests in flight on a link of 1 GiB/s, one stream moves about 4 MiB/s, so 60 ms of it
    /// are about 256 KiB.
    Request,
}

impl ReadCost {
    /// Ranges with fewer bytes than this between them are read as one.
    pub fn join_gap(self) -> u64 {
        match self {
            Self::Local => 8 << 10,
            Self::Network => 64 << 10,
            Self::Request => 256 << 10,
        }
    }

    /// Joined reads stop growing at this size.
    pub fn join_limit(self) -> u64 {
        match self {
            Self::Local => 1 << 20,
            Self::Network => 4 << 20,
            Self::Request => 8 << 20,
        }
    }
}

impl<T: ReadAt + ?Sized> ReadAt for Box<T> {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        (**self).read_exact_at(offset, buffer)
    }

    fn size(&self) -> u64 {
        (**self).size()
    }

    fn will_need(&self, offset: u64, length: u64) {
        (**self).will_need(offset, length)
    }
    fn version(&self) -> io::Result<u64> {
        (**self).version()
    }

    fn read_cost(&self) -> ReadCost {
        (**self).read_cost()
    }
}

impl ReadAt for Vec<u8> {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        read_slice(self, offset, buffer)
    }

    fn size(&self) -> u64 {
        self.len() as u64
    }
}

/// Ranges closer than this many bytes are hinted as one at least, so reads become fewer and
/// larger.
const HINT_GAP: u64 = 32 << 10;

/// Tells `source` that the byte ranges `start..end` will be read soon, joining ranges that lie
/// close together by its read cost. It never changes what is read.
pub fn will_need_ranges(source: &dyn ReadAt, ranges: &mut [(u64, u64)]) {
    ranges.sort_unstable();
    let cost = source.read_cost();
    let (gap, limit) = (cost.join_gap().max(HINT_GAP), cost.join_limit());
    let mut joined: Option<(u64, u64)> = None;
    for &(start, end) in ranges.iter() {
        match &mut joined {
            Some(last)
                if start <= last.1.saturating_add(gap) && end.max(last.1) - last.0 <= limit =>
            {
                last.1 = last.1.max(end)
            }
            _ => {
                if let Some((first, last)) = joined.replace((start, end)) {
                    source.will_need(first, last - first);
                }
            }
        }
    }
    if let Some((first, last)) = joined {
        source.will_need(first, last - first);
    }
}

/// Bytes of one piece of a long read, and the pieces asked for ahead of the one being read.
const AHEAD_PIECE: usize = 16 << 20;
const AHEAD_PIECES: u64 = 4;

/// Reads like `read_exact_at`, in pieces, and asks for the next pieces before each is read, so a
/// long read of a slow source has many requests in flight but only a few pieces ahead.
pub fn read_exact_ahead(source: &dyn ReadAt, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
    let end = checked_end(offset, buffer.len(), source.size())?;
    let mut hinted = offset;
    for (index, piece) in buffer.chunks_mut(AHEAD_PIECE).enumerate() {
        let start = offset + (index * AHEAD_PIECE) as u64;
        let until = end.min(start + AHEAD_PIECES * AHEAD_PIECE as u64);
        if hinted < until {
            source.will_need(hinted, until - hinted);
            hinted = until;
        }
        source.read_exact_at(start, piece)?;
    }
    Ok(())
}

/// Locks `mutex`. A thread that panicked while holding it leaves data that is still whole here,
/// so the lock is taken anyway.
pub(crate) fn lock<T>(mutex: &Mutex<T>) -> MutexGuard<'_, T> {
    mutex.lock().unwrap_or_else(PoisonError::into_inner)
}

/// Fills `buffer` with the bytes of `bytes` that start at `offset`, or fails past the end.
pub(crate) fn read_slice(bytes: &[u8], offset: u64, buffer: &mut [u8]) -> io::Result<()> {
    let part = usize::try_from(offset)
        .ok()
        .and_then(|start| bytes.get(start..start.checked_add(buffer.len())?))
        .ok_or_else(past_the_end)?;
    buffer.copy_from_slice(part);
    Ok(())
}

/// The error of a read that starts or ends after the last byte.
pub(crate) fn past_the_end() -> io::Error {
    io::Error::new(io::ErrorKind::UnexpectedEof, "read past the end")
}

/// The end of the range `offset..offset + length`, if it lies within `size` bytes.
pub(crate) fn checked_end(offset: u64, length: usize, size: u64) -> io::Result<u64> {
    offset
        .checked_add(length as u64)
        .filter(|end| *end <= size)
        .ok_or_else(past_the_end)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Records the hints it gets; its reads cost `.1`.
    struct Hinted(Mutex<Vec<(u64, u64)>>, ReadCost);

    impl Default for Hinted {
        fn default() -> Self {
            Self(Mutex::default(), ReadCost::Local)
        }
    }

    impl ReadAt for Hinted {
        fn read_exact_at(&self, _offset: u64, _buffer: &mut [u8]) -> io::Result<()> {
            Ok(())
        }

        fn size(&self) -> u64 {
            u64::MAX
        }

        fn will_need(&self, offset: u64, length: u64) {
            lock(&self.0).push((offset, length));
        }

        fn read_cost(&self) -> ReadCost {
            self.1
        }
    }

    #[test]
    fn close_ranges_joined() {
        let source = Hinted::default();
        let mut ranges = [(HINT_GAP * 10, HINT_GAP * 10 + 8), (0, 8), (100, 4_000)];
        will_need_ranges(&source, &mut ranges);
        assert_eq!(*lock(&source.0), vec![(0, 4_000), (HINT_GAP * 10, 8)]);
        will_need_ranges(&source, &mut []);
        assert_eq!(lock(&source.0).len(), 2);
    }

    /// A source of requests joins ranges 200 KiB apart, up to the join limit.
    #[test]
    fn requests_joined_wider() {
        let source = Hinted(Mutex::default(), ReadCost::Request);
        let step = 200 << 10;
        let mut ranges: Vec<(u64, u64)> = (0..50).map(|at| (at * step, at * step + 8)).collect();
        will_need_ranges(&source, &mut ranges);
        let hints = lock(&source.0).clone();
        let limit = ReadCost::Request.join_limit();
        assert!(hints.len() > 1 && hints.len() < 5, "{hints:?}");
        assert!(hints.iter().all(|(_, length)| *length <= limit));
        let local = Hinted::default();
        will_need_ranges(&local, &mut ranges);
        assert_eq!(lock(&local.0).len(), 50);
    }

    #[test]
    fn vec_reads_offsets() {
        let data: Vec<u8> = (0..100).collect();
        let mut buffer = [0; 10];
        data.read_exact_at(90, &mut buffer).unwrap();
        assert_eq!(buffer, [90, 91, 92, 93, 94, 95, 96, 97, 98, 99]);
        assert_eq!(
            data.read_exact_at(91, &mut buffer).unwrap_err().kind(),
            io::ErrorKind::UnexpectedEof
        );
        assert!(data.read_exact_at(u64::MAX, &mut buffer).is_err());
        assert_eq!(data.size(), 100);
    }
}
