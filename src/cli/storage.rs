//! How the commands open shards: local directories, or S3 prefixes read through a block cache.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::CliError;
use crate::index::{Location, Shard, ShardOptions};
use crate::storage::{BlockCache, CacheConfig, S3Config, S3Limits};
use std::path::Path;
use std::sync::{Arc, OnceLock};
use tokio::runtime::Runtime;

/// RAM of the block cache for shards on S3.
const CACHE_RAM: u64 = 1 << 30;
/// Bytes per block of the block cache, the size of most S3 reads.
const BLOCK_BYTES: u64 = 64 << 10;

/// The storage of the shards that one command opens. Reads of S3 shards wait on its runtime,
/// so it must outlive them.
pub(super) struct ShardStorage {
    runtime: Runtime,
    s3: OnceLock<S3Config>,
    cache: OnceLock<Arc<BlockCache>>,
}

impl ShardStorage {
    /// Storage with its own runtime; the S3 store and the cache are set up at the first S3 shard.
    pub(super) fn new() -> Result<Self, CliError> {
        Ok(Self {
            runtime: Runtime::new().map_err(CliError::Runtime)?,
            s3: OnceLock::new(),
            cache: OnceLock::new(),
        })
    }

    /// Opens the shard at `path`, a local directory or `s3://bucket/prefix`.
    pub(super) fn open(&self, path: &Path) -> Result<Shard, CliError> {
        let failed = |source| CliError::Shard {
            path: path.to_path_buf(),
            source,
        };
        let text = path.to_string_lossy();
        if !text.starts_with("s3://") {
            return Shard::open(path).map_err(failed);
        }
        let location = Location::s3(&text, self.s3()?).map_err(failed)?;
        let options = ShardOptions {
            packed: None,
            cache: Some(self.cache()?),
        };
        Shard::open_with(&location, &options).map_err(failed)
    }

    fn s3(&self) -> Result<&S3Config, CliError> {
        if let Some(config) = self.s3.get() {
            return Ok(config);
        }
        let handle = self.runtime.handle().clone();
        let config = S3Config::from_environment(S3Limits::default(), handle)?;
        Ok(self.s3.get_or_init(|| config))
    }

    fn cache(&self) -> Result<Arc<BlockCache>, CliError> {
        if let Some(cache) = self.cache.get() {
            return Ok(Arc::clone(cache));
        }
        let config = CacheConfig::new(CACHE_RAM, BLOCK_BYTES);
        let cache = BlockCache::open(&config, self.runtime.handle().clone())?;
        Ok(Arc::clone(self.cache.get_or_init(|| cache)))
    }
}
