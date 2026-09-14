//! How the commands open shards: local directories, or S3 prefixes read through pinned parts of
//! a warmed cache directory and a block cache.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::CliError;
use crate::index::{Location, Place, Plan, Shard, ShardOptions};
use crate::storage::{BlockCache, CacheConfig, DiskConfig, S3Config, S3Limits};
use bytesize::ByteSize;
use clap::Args;
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};
use tokio::runtime::Runtime;

/// RAM for pins and the block cache when neither the command line nor a plan gives it.
const DEFAULT_MEMORY: u64 = 1 << 30;
/// The least RAM of the block cache, whatever the pins take.
const MIN_CACHE_RAM: u64 = 64 << 20;
/// Bytes per block of the block cache, the size of most S3 reads.
const BLOCK_BYTES: u64 = 64 << 10;
/// Hinted ranges that the block cache fetches at once, as many as S3 requests may be in flight.
const READ_AHEAD: usize = 256;
/// The directory of the block cache's disk level in the cache directory.
const BLOCKS_DIRECTORY: &str = "blocks";

/// Where shards on S3 keep parts locally.
#[derive(Debug, Args)]
pub(super) struct CacheArgs {
    /// A cache directory, warmed by `jam index warm`: shards on S3 read their pinned parts from
    /// there, and the block cache keeps blocks there for later runs
    #[arg(long)]
    cache: Option<PathBuf>,
    /// RAM for pins and the block cache, for example 16GiB [default: the plan's, else 1GiB]
    #[arg(long)]
    memory: Option<ByteSize>,
    /// Disk in the cache directory for pins and the block cache [default: the plan's, else 0]
    #[arg(long, requires = "cache")]
    disk: Option<ByteSize>,
}

/// The storage of the shards that one command opens. Reads of S3 shards wait on its runtime,
/// so it must outlive them.
pub(super) struct ShardStorage {
    runtime: Runtime,
    s3: OnceLock<S3Config>,
    cache: Option<Arc<BlockCache>>,
}

impl ShardStorage {
    /// Storage with its own runtime; the S3 store is set up at the first S3 shard.
    pub(super) fn new() -> Result<Self, CliError> {
        Ok(Self {
            runtime: Runtime::new().map_err(CliError::Runtime)?,
            s3: OnceLock::new(),
            cache: None,
        })
    }

    /// The location of the shard at `path`, a local directory or `s3://bucket/prefix`.
    pub(super) fn location(&self, path: &Path) -> Result<Location, CliError> {
        let text = path.to_string_lossy();
        if !text.starts_with("s3://") {
            return Ok(Location::Directory(path.to_path_buf()));
        }
        Location::s3(&text, self.s3()?).map_err(|source| CliError::Shard {
            path: path.to_path_buf(),
            source,
        })
    }

    /// Opens the shards at `paths`. Shards on S3 read the pins of their plans in the cache
    /// directory and everything else through one block cache, which gets the budgets that the
    /// pins of all plans leave.
    pub(super) fn open_all(
        &mut self,
        paths: &[PathBuf],
        args: &CacheArgs,
    ) -> Result<Vec<Shard>, CliError> {
        let mut plans = Vec::new();
        let mut locations = Vec::new();
        for path in paths {
            let location = self.location(path)?;
            if let (Location::S3 { .. }, Some(cache)) = (&location, &args.cache) {
                let failed = |source| CliError::Shard {
                    path: path.clone(),
                    source,
                };
                let (_, digest) = location.manifest().map_err(failed)?;
                plans.extend(Plan::read(cache, &digest).map_err(failed)?);
            }
            locations.push(location);
        }
        if locations.iter().any(|at| matches!(at, Location::S3 { .. })) {
            self.cache = Some(self.block_cache(&plans, args)?);
        }
        let options = ShardOptions {
            cache: self.cache.clone(),
            pins: args.cache.clone(),
            ..ShardOptions::default()
        };
        let open = |(path, location): (&PathBuf, Location)| {
            let options = match location {
                Location::Directory(_) => ShardOptions::default(),
                Location::S3 { .. } => options.clone(),
            };
            Shard::open_with(&location, &options).map_err(|source| CliError::Shard {
                path: path.clone(),
                source,
            })
        };
        paths.iter().zip(locations).map(open).collect()
    }

    /// The block cache with the budgets of the command line, else of the plans, less their pins.
    fn block_cache(&self, plans: &[Plan], args: &CacheArgs) -> Result<Arc<BlockCache>, CliError> {
        let largest = |budget: fn(&Plan) -> u64| plans.iter().map(budget).max();
        let memory = (args.memory.map(|size| size.as_u64()))
            .or(largest(|plan| plan.memory))
            .unwrap_or(DEFAULT_MEMORY);
        let disk = (args.disk.map(|size| size.as_u64()))
            .or(largest(|plan| plan.disk))
            .unwrap_or(0);
        let pins = plans.iter().flat_map(|plan| &plan.pins);
        let ram_pins: u64 = (pins.clone())
            .filter(|pin| pin.place == Place::Ram)
            .map(|pin| pin.bytes())
            .sum();
        let disk_pins: u64 = pins.map(|pin| pin.bytes()).sum();
        let mut config = CacheConfig {
            read_ahead: READ_AHEAD,
            ..CacheConfig::new(
                memory.saturating_sub(ram_pins).max(MIN_CACHE_RAM),
                BLOCK_BYTES,
            )
        };
        let disk = disk.saturating_sub(disk_pins);
        if let Some(directory) = &args.cache
            && disk > 0
        {
            let directory = directory.join(BLOCKS_DIRECTORY);
            config.disk = Some(DiskConfig {
                directory,
                bytes: disk,
            });
        }
        Ok(BlockCache::open(&config, self.runtime.handle().clone())?)
    }

    fn s3(&self) -> Result<&S3Config, CliError> {
        if let Some(config) = self.s3.get() {
            return Ok(config);
        }
        let handle = self.runtime.handle().clone();
        let config = S3Config::from_environment(S3Limits::default(), handle)?;
        Ok(self.s3.get_or_init(|| config))
    }
}

impl Drop for ShardStorage {
    fn drop(&mut self) {
        // The disk level is best effort: blocks that fail to reach it are fetched again later.
        if let Some(cache) = &self.cache {
            let _ = cache.close();
        }
    }
}
