//! The on-disk index: per shard its manifest, seeds file and packed sequences; building shards,
//! and resharding them into an index of several shards.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

mod bits;
#[allow(dead_code)]
mod buckets;
mod build;
mod contig_table;
mod index_file;
mod manifest;
mod packed;
mod pages;
mod plan;
mod reshard;
mod seed_file;
mod seed_writer;
#[allow(dead_code)]
mod sequence_file;
mod shard;
mod sources;
mod spill;
mod warm;

use crate::seeds::SeedError;
use crate::storage::StorageError;
use std::io;
use std::path::{Path, PathBuf};
use thiserror::Error;

pub use build::{BuildOptions, EntryInput, build_shard};
pub use contig_table::ContigTable;
pub use index_file::{IndexFile, ShardRecord};
pub use manifest::{EntryRecord, FileRecord, Manifest, Parameters};
pub use packed::{Contig, PackedSequence};
pub use plan::{Part, Pin, Place, Plan, part_pins};
pub use reshard::{PackedOutput, ReshardOptions, Split, reshard};
pub use seed_file::{
    ColorSetId, ColorSets, SEEDS_FILE, SeedFile, SeedLocation, StoredSeed, seed_key,
};
pub use shard::{Shard, ShardStats};
pub use sources::{Location, ShardOptions};
pub use warm::warm;

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

/// The most bucket bits a shard may use, so a check value keeps at least 16 bits.
pub const MAX_BUCKET_BITS: u8 = 28;

/// The fewest bucket bits a shard uses.
pub const MIN_BUCKET_BITS: u8 = 8;

#[cfg(test)]
mod tests;

#[cfg(test)]
mod etag_tests;
