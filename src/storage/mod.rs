//! Reading bytes by offset from local files, S3 byte ranges and the block cache.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

mod cache;
mod error;
mod local;
mod mapped;
mod readahead;
mod s3;

pub use cache::{BlockCache, CacheConfig, CachedSource, DiskConfig};
pub use error::StorageError;
pub use local::LocalFile;
pub use mapped::MappedFile;
pub use s3::{S3Config, S3Limits, S3Object};

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
}

impl ReadAt for Vec<u8> {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        read_slice(self, offset, buffer)
    }

    fn size(&self) -> u64 {
        self.len() as u64
    }
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

    #[test]
    fn bytes_in_memory_are_read_by_offset() {
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
