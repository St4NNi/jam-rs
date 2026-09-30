//! How the commands open shards: local directories, or S3 prefixes read through pinned parts of
//! a warmed cache directory and a block cache.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::CliError;
use crate::index::{Location, Place, Plan, Shard, ShardOptions, warm};
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
#[derive(Clone, Debug, Default, Args)]
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

    /// Opens the shards at `paths`: pins of their plans from the cache directory, the rest
    /// through one block cache with the budgets the pins leave. With `warm_unplanned`, a shard
    /// without a plan first gets the default plan if `--disk` is given (for a server).
    pub(super) fn open_all(
        &mut self,
        paths: &[PathBuf],
        args: &CacheArgs,
        warm_unplanned: bool,
    ) -> Result<Vec<Shard>, CliError> {
        let mut plans = Vec::new();
        let mut opened = Vec::new();
        for path in paths {
            let failed = |source| CliError::Shard {
                path: path.clone(),
                source,
            };
            let location = self.location(path)?;
            if let Some(cache) = &args.cache {
                location.keep_etags(cache).map_err(failed)?;
            }
            let manifest = match location {
                Location::S3 { .. } => Some(location.manifest().map_err(failed)?),
                Location::Directory(_) => None,
            };
            if let (Some((manifest, digest)), Some(cache)) = (&manifest, &args.cache) {
                let mut plan = Plan::read(cache, digest).map_err(failed)?;
                if let (None, true, Some(disk)) = (&plan, warm_unplanned, args.disk) {
                    let memory = args.memory.map_or(DEFAULT_MEMORY, |size| size.as_u64());
                    let chosen = Plan::choose(manifest, digest, memory, disk.as_u64(), &[])?;
                    warm(&location, manifest, &chosen, cache).map_err(failed)?;
                    plan = Some(chosen);
                }
                plans.extend(plan);
            }
            opened.push((path, location, manifest));
        }
        // A shard opened again, for example after its upload changed, keeps the block cache.
        if self.cache.is_none() && opened.iter().any(|(_, _, manifest)| manifest.is_some()) {
            self.cache = Some(self.block_cache(&plans, args)?);
        }
        let options = ShardOptions {
            cache: self.cache.clone(),
            pins: args.cache.clone(),
            ..ShardOptions::default()
        };
        let open = |(path, location, manifest): (&PathBuf, Location, _)| {
            let shard = match manifest {
                None => Shard::open_with(&location, &ShardOptions::default()),
                Some((manifest, digest)) => {
                    Shard::open_with_manifest(&location, &options, manifest, digest)
                }
            };
            shard.map_err(|source| CliError::Shard {
                path: path.clone(),
                source,
            })
        };
        opened.into_iter().map(open).collect()
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{random_bases, rng, runtime, serve_directory, test_shard};

    #[test]
    fn serve_warms_unplanned() {
        let directory = tempfile::tempdir().unwrap();
        let mut rng = rng(9);
        let entries: Vec<Vec<Vec<u8>>> = (0..2)
            .map(|_| vec![random_bases(&mut rng, 20_000)])
            .collect();
        let local = test_shard(directory.path(), &entries);
        let (server, asked) = serve_directory(directory.path());
        let endpoint = Some(server.address.as_str());
        let config = S3Config::new("local", endpoint, true, S3Limits::default(), runtime());
        let config = config.unwrap().unsigned();
        let cache = directory.path().join("cache");
        let args = CacheArgs {
            cache: Some(cache.clone()),
            memory: Some(ByteSize::gib(1)),
            disk: Some(ByteSize::gib(1)),
        };
        let paths = [PathBuf::from("s3://bucket/shard")];
        let open = |warm_unplanned| {
            let mut storage = ShardStorage::new().unwrap();
            assert!(storage.s3.set(config.clone()).is_ok());
            storage.open_all(&paths, &args, warm_unplanned).unwrap();
        };
        // A search reads through the block cache alone; a server warms the default plan.
        open(false);
        assert_eq!(Plan::read(&cache, local.digest()).unwrap(), None);
        open(true);
        let plan = Plan::read(&cache, local.digest()).unwrap().unwrap();
        let parts: Vec<(&str, Place)> = (plan.pins.iter())
            .map(|pin| (pin.part.as_str(), pin.place))
            .collect();
        assert_eq!(parts, [("sketch", Place::Ram), ("locators", Place::Ram)]);
        // The next start finds the plan and reads the seeds file from the pins alone.
        let before = asked.lock().unwrap().len();
        open(true);
        let asked = asked.lock().unwrap();
        assert!(asked[before..].iter().all(|key| !key.ends_with("seeds")));
    }

    /// A server opens its shard again with the same storage, so with the same block cache.
    #[test]
    fn shard_opened_twice() {
        let directory = tempfile::tempdir().unwrap();
        let mut rng = rng(10);
        let entries = vec![vec![random_bases(&mut rng, 20_000)]];
        test_shard(directory.path(), &entries);
        let (server, _) = serve_directory(directory.path());
        let endpoint = Some(server.address.as_str());
        let config = S3Config::new("local", endpoint, true, S3Limits::default(), runtime());
        let args = CacheArgs {
            cache: Some(directory.path().join("cache")),
            memory: Some(ByteSize::gib(1)),
            disk: Some(ByteSize::gib(1)),
        };
        let mut storage = ShardStorage::new().unwrap();
        assert!(storage.s3.set(config.unwrap().unsigned()).is_ok());
        let paths = [PathBuf::from("s3://bucket/shard")];
        assert_eq!(storage.open_all(&paths, &args, true).unwrap().len(), 1);
        assert_eq!(storage.open_all(&paths, &args, true).unwrap().len(), 1);
    }
}
