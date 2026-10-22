//! How the commands open shards: local directories, or S3 prefixes read through pinned parts of
//! a warmed cache directory and a block cache.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::CliError;
use super::memory::{available, default_budget};
use crate::index::{Location, Place, Plan, Shard, ShardOptions, warm};
use crate::storage::{BlockCache, CacheConfig, DiskConfig, S3Config, S3Limits};
use bytesize::ByteSize;
use clap::Args;
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};
use tokio::runtime::Runtime;

/// The least RAM of the block cache, whatever the pins take.
const MIN_CACHE_RAM: u64 = 64 << 20;
/// Bytes per block of the block cache, the size of most S3 reads.
const BLOCK_BYTES: u64 = 64 << 10;
/// S3 requests in flight at once when the command line does not say, and the most it may say.
/// The block cache fetches as many hinted ranges at once, each on a thread of the runtime.
const REQUESTS: usize = 512;
const MAX_REQUESTS: usize = 4096;
/// Threads of the runtime for other blocking work, such as the disk level of the block cache.
const OTHER_BLOCKING: usize = 256;
/// Open files kept for shards, pins and the disk level; each request in flight needs one more,
/// as each has its own connection.
const OTHER_FILES: u64 = 1024;
/// The directory of the block cache's disk level in the cache directory.
const BLOCKS_DIRECTORY: &str = "blocks";

/// Where shards on S3 keep parts locally.
#[derive(Clone, Debug, Default, Args)]
pub(super) struct CacheArgs {
    /// A cache directory, warmed by `jam index warm`: shards on S3 read their pinned parts from
    /// there, and the block cache keeps blocks there for later runs
    #[arg(long)]
    cache: Option<PathBuf>,
    /// RAM for pins, the block cache and the buffers of S3 reads, for example 16GiB [default:
    /// the plan's, else half the RAM available, less what the searches take, at least 1GiB]
    #[arg(long)]
    memory: Option<ByteSize>,
    /// Disk in the cache directory for pins and the block cache [default: the plan's, else 0]
    #[arg(long, requires = "cache")]
    disk: Option<ByteSize>,
    /// S3 requests in flight at once, up to 4096 [default: 512]
    #[arg(long, value_parser = clap::value_parser!(u16).range(1..=MAX_REQUESTS as i64))]
    s3_requests: Option<u16>,
}

impl CacheArgs {
    /// S3 requests in flight at once.
    fn requests(&self) -> usize {
        self.s3_requests.map_or(REQUESTS, usize::from)
    }
}

/// The storage of the shards that one command opens. Reads of S3 shards wait on its runtime,
/// so it must outlive them.
pub(super) struct ShardStorage {
    runtime: Runtime,
    s3: OnceLock<S3Config>,
    cache: Option<Arc<BlockCache>>,
    /// S3 requests in flight at once, set by the first `open_all`.
    requests: usize,
    /// The most S3 requests in flight that the open files limit leaves room for.
    fit: usize,
    /// Searches that run at the same time, which the default RAM budget leaves room for.
    searches: usize,
    /// The RAM budget, the RAM pins and the RAM of the block cache, once it is opened.
    budget: Option<(u64, u64, u64)>,
}

impl ShardStorage {
    /// Storage with its own runtime; the S3 store is set up at the first S3 shard. The runtime
    /// may start a blocking thread for each request in flight; it starts them only when needed.
    pub(super) fn new() -> Result<Self, CliError> {
        let runtime = tokio::runtime::Builder::new_multi_thread()
            .enable_all()
            .max_blocking_threads(MAX_REQUESTS + OTHER_BLOCKING)
            .build();
        // Each request in flight holds a connection, so the open files limit bounds them.
        let files = rlimit::increase_nofile_limit(u64::MAX).unwrap_or(OTHER_FILES * 2);
        let fit = files.saturating_sub(OTHER_FILES).max(16) as usize;
        Ok(Self {
            runtime: runtime.map_err(CliError::Runtime)?,
            s3: OnceLock::new(),
            cache: None,
            requests: REQUESTS.min(fit),
            fit,
            searches: 1,
            budget: None,
        })
    }

    /// Leaves room in the default RAM budget for `searches` searches at the same time.
    pub(super) fn set_searches(&mut self, searches: usize) {
        self.searches = searches;
    }

    /// A line on the RAM budget and how it is shared, once the block cache is open.
    pub(super) fn budget_line(&self) -> Option<String> {
        let (memory, pins, cache) = self.budget?;
        let size = |bytes| ByteSize(bytes).display().iec().to_string();
        Some(format!(
            "jam: memory budget {}: {} pinned in RAM, {} block cache, {} for reads in flight",
            size(memory),
            size(pins),
            size(cache),
            size(BlockCache::fetch_room(cache))
        ))
    }

    /// The location of the index file at `path`: `NAME.jam` or `s3://bucket/prefix/NAME.jam`.
    pub(super) fn location(&self, path: &Path) -> Result<Location, CliError> {
        let text = path.to_string_lossy();
        if !text.starts_with("s3://") {
            return Ok(Location::File(path.to_path_buf()));
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
        if self.s3.get().is_none() {
            self.requests = args.requests().min(self.fit);
            if args.requests() > self.fit {
                eprintln!(
                    "jam: {} S3 requests at once, as open files are limited",
                    self.fit
                );
            }
        }
        let mut located = Vec::new();
        for path in paths {
            let location = self.location(path)?;
            if let Some(cache) = &args.cache {
                (location.keep_etags(cache)).map_err(|source| CliError::Shard {
                    path: path.clone(),
                    source,
                })?;
            }
            located.push((path, location));
        }
        // The heads of shards on S3 are read at the same time, each with a few requests.
        let heads: Vec<_> = (located.par_iter())
            .map(|(_, location)| match location {
                Location::S3 { .. } => location.head().map(Some),
                Location::File(_) => Ok(None),
            })
            .collect();
        for ((path, location), head) in located.into_iter().zip(heads) {
            let failed = |source| CliError::Shard {
                path: path.clone(),
                source,
            };
            let head = head.map_err(failed)?;
            if let (Some(head), Some(cache)) = (&head, &args.cache) {
                let mut plan = Plan::read(cache, &head.digest).map_err(failed)?;
                if let (None, true, Some(disk)) = (&plan, warm_unplanned, args.disk) {
                    let memory = (args.memory.map(|size| size.as_u64()))
                        .unwrap_or_else(|| default_budget(self.searches));
                    let chosen = Plan::choose(head, memory, disk.as_u64(), &[])?;
                    warm(&location, head, &chosen, cache).map_err(failed)?;
                    plan = Some(chosen);
                }
                plans.extend(plan);
            }
            opened.push((path, location, head));
        }
        // A shard opened again, for example after its upload changed, keeps the block cache.
        if self.cache.is_none() && opened.iter().any(|(_, _, head)| head.is_some()) {
            self.cache = Some(self.block_cache(&plans, args)?);
        }
        let options = ShardOptions {
            cache: self.cache.clone(),
            pins: args.cache.clone(),
        };
        let open = |(path, location, head): (&PathBuf, Location, _)| {
            let shard = match head {
                None => Shard::open_with(&location, &ShardOptions::default()),
                Some(head) => Shard::open_with_head(&location, &options, head),
            };
            shard.map_err(|source| CliError::Shard {
                path: path.clone(),
                source,
            })
        };
        // Shards open at the same time, as opening one on S3 waits for several reads in turn.
        let shards: Vec<_> = opened.into_par_iter().map(open).collect();
        shards.into_iter().collect()
    }

    /// The block cache with the budgets of the command line, else of the plans, less their pins.
    /// The RAM pins, the block cache and the room of the reads in flight stay within the RAM
    /// budget, but the block cache gets at least `MIN_CACHE_RAM`.
    fn block_cache(
        &mut self,
        plans: &[Plan],
        args: &CacheArgs,
    ) -> Result<Arc<BlockCache>, CliError> {
        let largest = |budget: fn(&Plan) -> u64| plans.iter().map(budget).max();
        let given = (args.memory.map(|size| size.as_u64())).or(largest(|plan| plan.memory));
        let memory = given.unwrap_or_else(|| default_budget(self.searches));
        if let (Some(given), Some(available)) = (given, available())
            && given > available
        {
            let size = |bytes| ByteSize(bytes).display().iec().to_string();
            eprintln!(
                "jam: the memory budget {} is more than the {} available",
                size(given),
                size(available)
            );
        }
        let disk = (args.disk.map(|size| size.as_u64()))
            .or(largest(|plan| plan.disk))
            .unwrap_or(0);
        let pins = plans.iter().flat_map(|plan| &plan.pins);
        let ram_pins: u64 = (pins.clone())
            .filter(|pin| pin.place == Place::Ram)
            .map(|pin| pin.bytes())
            .sum();
        let disk_pins: u64 = pins.map(|pin| pin.bytes()).sum();
        let left = memory.saturating_sub(ram_pins);
        let cache_ram = (left - BlockCache::fetch_room(left).min(left)).max(MIN_CACHE_RAM);
        self.budget = Some((memory, ram_pins, cache_ram));
        let mut config = CacheConfig {
            read_ahead: self.requests,
            ..CacheConfig::new(cache_ram, BLOCK_BYTES)
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
        let limits = S3Limits {
            max_in_flight: self.requests,
            ..S3Limits::default()
        };
        let config = S3Config::from_environment(limits, handle)?;
        Ok(self.s3.get_or_init(|| config))
    }
}

/// The name of the index file at `path` in records: its file name without `.jam`.
pub(super) fn shard_name(path: &Path) -> String {
    let name = path
        .file_name()
        .map(|name| name.to_string_lossy().into_owned());
    let name = name.unwrap_or_else(|| path.display().to_string());
    name.strip_suffix(".jam")
        .map_or(name.clone(), str::to_string)
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
            s3_requests: None,
        };
        let paths = [PathBuf::from("s3://bucket/shard.jam")];
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
        let all = [
            ("headers", Place::Ram),
            ("sketch", Place::Ram),
            ("locators", Place::Ram),
        ];
        assert_eq!(parts, all);
        // The next start finds the plan and reads the seeds from the pins alone: the seeds file
        // only gets the reads of its header and manifest, which also check it.
        let before = asked.lock().unwrap().len();
        open(true);
        let asked = asked.lock().unwrap();
        let seeds = asked[before..].iter().filter(|k| *k == "shard.jam").count();
        assert_eq!(seeds, 2);
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
            s3_requests: None,
        };
        let mut storage = ShardStorage::new().unwrap();
        assert!(storage.s3.set(config.unwrap().unsigned()).is_ok());
        let paths = [PathBuf::from("s3://bucket/shard.jam")];
        assert_eq!(storage.open_all(&paths, &args, true).unwrap().len(), 1);
        assert_eq!(storage.open_all(&paths, &args, true).unwrap().len(), 1);
    }
}
