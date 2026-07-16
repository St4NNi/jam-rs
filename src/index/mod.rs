//! The on-disk index: per shard its manifest, sketch, positions and packed sequences; building
//! shards, and resharding them into an index of several shards.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

mod bits;
mod build;
mod color_sets;
mod index_file;
mod manifest;
mod packed;
mod positions;
mod reshard;
mod shard;
mod sketch;
mod sources;

use crate::seeds::SeedError;
use crate::storage::StorageError;
use std::io;
use std::path::{Path, PathBuf};
use thiserror::Error;

pub use build::{BuildOptions, EntryInput, build_shard};
pub use index_file::{IndexFile, ShardRecord};
pub use manifest::{EntryRecord, FileRecord, Manifest, Parameters};
pub use packed::PackedSequence;
pub use positions::{Positions, SeedLocation};
pub use reshard::{PackedOutput, ReshardOptions, Split, reshard};
pub use shard::{Shard, ShardStats};
pub use sketch::{ColorSetEntries, ColorSetId, ColorSets, Sketch};
pub use sources::{PackedLocation, PackedOptions};

/// Errors of the index module.
#[derive(Debug, Error)]
pub enum IndexError {
    /// Reading a shard file through its source failed.
    #[error("reading the index failed: {0}")]
    Io(#[from] io::Error),
    /// A file or directory of the index could not be read, written or created.
    #[error("{}: {source}", path.display())]
    File {
        /// The file or directory.
        path: PathBuf,
        /// What went wrong.
        #[source]
        source: io::Error,
    },
    /// The thread pool of a build could not start.
    #[error("cannot start the build threads: {0}")]
    Threads(#[from] rayon::ThreadPoolBuildError),
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
    /// A packed sequence source could not be opened.
    #[error(transparent)]
    Storage(#[from] StorageError),
}

impl IndexError {
    fn format(message: impl Into<String>) -> Self {
        Self::Format(message.into())
    }
}

/// Names the file of a failed file operation.
pub(crate) trait AtPath<T> {
    /// The result, with an error that names `path`.
    fn at(self, path: &Path) -> Result<T, IndexError>;
}

impl<T> AtPath<T> for io::Result<T> {
    fn at(self, path: &Path) -> Result<T, IndexError> {
        self.map_err(|source| IndexError::File {
            path: path.to_path_buf(),
            source,
        })
    }
}

/// Bits of a seed's hash that a sketch bucket stores to tell its seeds apart.
pub const FINGERPRINT_BITS: u8 = 16;

/// The most bucket bits a shard may use, so a stored seed key fits the build's spill records.
pub const MAX_BUCKET_BITS: u8 = 28;

/// The fewest bucket bits a shard uses.
pub const MIN_BUCKET_BITS: u8 = 8;

/// The bits of a hash that the index keeps: the bucket bits followed by the fingerprint.
pub(crate) fn seed_key(hash: u64, key_bits: u32) -> u64 {
    hash >> (64 - key_bits)
}

#[cfg(test)]
mod tests;
