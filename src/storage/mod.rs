//! Reading bytes by offset from local files, S3 byte ranges and the block cache.

use std::io;

/// Positional reads from an immutable file or object, shared by many threads.
pub trait ReadAt: Send + Sync {
    /// Fills `buffer` with the bytes that start at `offset`, or fails.
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()>;

    /// The size of the file or object in bytes.
    fn size(&self) -> u64;

    /// Tells the source that a range will be read soon. It never changes what is read.
    fn will_need(&self, _offset: u64, _length: u64) {}
}
