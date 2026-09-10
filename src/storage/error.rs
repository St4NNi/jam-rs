//! Why a storage source could not be opened or configured.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use std::io;
use std::path::PathBuf;
use thiserror::Error;

/// Why a source could not be opened or configured. Reads themselves fail with `io::Error`.
#[derive(Debug, Error)]
pub enum StorageError {
    /// A local file could not be opened or mapped.
    #[error("{}: {source}", path.display())]
    File {
        /// The file.
        path: PathBuf,
        /// What went wrong.
        #[source]
        source: io::Error,
    },
    /// Another process uses the disk level directory of the block cache.
    #[error("the cache directory {} is in use by another process", .0.display())]
    Busy(PathBuf),
    /// The location is not `s3://bucket/key`, or the bucket or key is empty.
    #[error("not a usable S3 location: {0}")]
    Location(String),
    /// The block cache could not be built.
    #[error("block cache: {0}")]
    Cache(#[from] foyer::Error),
    /// The S3 settings are not valid, for example an empty endpoint.
    #[error("invalid S3 setting: {0}")]
    Config(String),
    /// The object cannot be read.
    #[error("S3 object {location} cannot be read: {source}")]
    Object {
        /// The `s3://` location that was opened.
        location: String,
        /// What went wrong.
        #[source]
        source: io::Error,
    },
}
