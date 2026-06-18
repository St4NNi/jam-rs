//! The on-disk index: per shard its manifest, sketch, positions and packed sequences; building
//! shards, and resharding them into an index of several shards.

mod bits;
pub mod build;
pub mod index_file;
pub mod manifest;
pub mod packed;
pub mod positions;
pub mod reshard;
pub mod shard;
pub mod sketch;

use crate::seeds::SeedError;
use std::io;
use std::path::Path;
use thiserror::Error;

pub use build::{BuildOptions, EntryInput, build_shard};
pub use index_file::{IndexFile, ShardRecord};
pub use manifest::{EntryRecord, FileRecord, Manifest, Parameters};
pub use packed::{PackedSequence, PackedWriter};
pub use positions::{Positions, SeedLocation};
pub use reshard::{PackedOutput, ReshardOptions, Split, reshard};
pub use shard::{Shard, ShardStats};
pub use sketch::{ColorSetEntries, ColorSetId, ColorSets, Sketch};

/// Errors of the index module.
#[derive(Debug, Error)]
pub enum IndexError {
    /// Reading or writing a file failed.
    #[error("{0}")]
    Io(#[from] io::Error),
    /// A file is not in the expected format, or is damaged.
    #[error("bad index file: {0}")]
    Format(String),
    /// A checksum or digest does not match the data.
    #[error("checksum mismatch: {0}")]
    Checksum(String),
    /// The caller gave input that cannot be indexed.
    #[error("{0}")]
    Input(String),
    /// The seed scheme is not valid.
    #[error(transparent)]
    Seeds(#[from] SeedError),
}

impl IndexError {
    fn format(message: impl Into<String>) -> Self {
        Self::Format(message.into())
    }
}

/// Bits of a seed's hash that a sketch bucket stores to tell its seeds apart.
pub const FINGERPRINT_BITS: u8 = 16;

/// The most bucket bits a shard may use, so a stored seed key fits the build's spill records.
pub const MAX_BUCKET_BITS: u32 = 28;

/// The bits of a hash that the index keeps: the bucket bits followed by the fingerprint.
pub(crate) fn seed_key(hash: u64, key_bits: u32) -> u64 {
    hash >> (64 - key_bits)
}

/// A local file mapped into memory and read by offset. It stands in for the local sources of the
/// storage module.
pub struct MappedFile {
    map: Option<memmap2::Mmap>,
}

impl MappedFile {
    /// Maps the file at `path`; `populate` reads it into RAM at once.
    pub fn open(path: &Path, populate: bool) -> Result<Self, IndexError> {
        let file = std::fs::File::open(path)?;
        if file.metadata()?.len() == 0 {
            return Ok(Self { map: None });
        }
        let mut options = memmap2::MmapOptions::new();
        if populate {
            options.populate();
        }
        // SAFETY: index files are written once and never changed while a shard is open.
        let map = unsafe { options.map(&file)? };
        Ok(Self { map: Some(map) })
    }

    /// All bytes of the file.
    pub fn bytes(&self) -> &[u8] {
        self.map.as_deref().unwrap_or_default()
    }
}

impl crate::storage::ReadAt for MappedFile {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        let bytes = usize::try_from(offset)
            .ok()
            .and_then(|start| self.bytes().get(start..start.checked_add(buffer.len())?))
            .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "read past the end"))?;
        buffer.copy_from_slice(bytes);
        Ok(())
    }

    fn size(&self) -> u64 {
        self.bytes().len() as u64
    }
}
#[cfg(test)]
mod tests;
