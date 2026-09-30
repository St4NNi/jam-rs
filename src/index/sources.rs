//! Where an open shard reads its files: a local directory or an S3 prefix, optionally through a
//! block cache shared by all files.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::manifest::{FileRecord, MANIFEST_FILE, Manifest, PACKED_DIRECTORY, sha256_hex};
use super::packed::BLOCK_BYTES;
use super::plan::{Pin, Place, Plan, etag_file};
use super::{AtPath, IndexError};
use crate::storage::{
    BlockCache, ChangeHook, Checked, LocalFile, MappedFile, PinnedSource, ReadAt, S3Config,
    S3Object,
};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

/// Where files lie: a local directory or an S3 prefix.
#[derive(Clone)]
pub enum Location {
    /// A local directory.
    Directory(PathBuf),
    /// An S3 prefix that ends with `/`; each file name is appended to it.
    S3 {
        /// The prefix, `s3://bucket/` followed by an optional key prefix.
        prefix: String,
        /// The store, the credentials and the limits.
        config: Box<S3Config>,
    },
}

impl Location {
    /// The S3 prefix `s3://bucket/prefix` on the store of `config`.
    pub fn s3(prefix: &str, config: &S3Config) -> Result<Self, IndexError> {
        let bucket = prefix
            .strip_prefix("s3://")
            .and_then(|rest| rest.split('/').next());
        if bucket.is_none_or(str::is_empty) {
            return Err(IndexError::Input(format!(
                "{prefix} is not an S3 prefix like s3://bucket/prefix/"
            )));
        }
        let mut prefix = prefix.to_string();
        if !prefix.ends_with('/') {
            prefix.push('/');
        }
        let config = Box::new(config.clone());
        Ok(Self::S3 { prefix, config })
    }

    /// The location of the directory `name` in this one.
    pub(crate) fn join(&self, name: &str) -> Self {
        match self {
            Self::Directory(path) => Self::Directory(path.join(name)),
            Self::S3 { prefix, config } => Self::S3 {
                prefix: format!("{prefix}{name}/"),
                config: config.clone(),
            },
        }
    }

    /// The local path of the file `name`, or `None` when it lies on S3.
    pub(crate) fn local_path(&self, name: &str) -> Option<PathBuf> {
        match self {
            Self::Directory(path) => Some(path.join(name)),
            Self::S3 { .. } => None,
        }
    }

    /// All bytes of the small file `name`.
    pub(crate) fn read_whole(&self, name: &str) -> Result<Vec<u8>, IndexError> {
        match self {
            Self::Directory(path) => {
                let path = path.join(name);
                std::fs::read(&path).at(&path)
            }
            Self::S3 { prefix, config } => {
                Ok(S3Object::read_whole(config, &format!("{prefix}{name}"))?)
            }
        }
    }

    /// The manifest of the shard at this location and the SHA-256 digest of its file.
    pub fn manifest(&self) -> Result<(Manifest, String), IndexError> {
        let bytes = self.read_whole(MANIFEST_FILE)?;
        Ok((Manifest::parse(&bytes)?, sha256_hex(&[&bytes])))
    }

    /// The S3 object of the file `name` of `bytes` bytes, or `None` in a directory. Nothing is
    /// sent, as the size is known.
    fn object(&self, name: &str, bytes: u64) -> Result<Option<S3Object>, IndexError> {
        Ok(match self {
            Self::Directory(_) => None,
            Self::S3 { prefix, config } => {
                Some(S3Object::open(config, &format!("{prefix}{name}"), bytes)?)
            }
        })
    }

    /// Keeps the ETags of S3 objects in `directory` too, so a later run finds them.
    pub(crate) fn keep_etags(&self, directory: &Path) -> Result<(), IndexError> {
        if let Self::S3 { config, .. } = self {
            config.keep_etags(directory)?;
        }
        Ok(())
    }

    /// The ETag that the file `name` has on S3 now, learned with one byte; `None` in a
    /// directory.
    pub(crate) fn refresh_etag(
        &self,
        name: &str,
        bytes: u64,
    ) -> Result<Option<String>, IndexError> {
        match self.object(name, bytes)? {
            Some(object) => Ok(object.refresh().at(Path::new(name))?),
            None => Ok(None),
        }
    }

    /// Fails if the file `name` no longer has `etag` on S3, with one conditional byte; a missing
    /// ETag fails too. Nothing is checked in a directory.
    pub(crate) fn verify_etag(
        &self,
        name: &str,
        bytes: u64,
        etag: Option<&str>,
    ) -> Result<(), IndexError> {
        match (self.object(name, bytes)?, etag) {
            (Some(object), Some(etag)) => Ok(object.verify(etag).at(Path::new(name))?),
            (Some(_), None) => Err(IndexError::Input(format!("{name} has no ETag"))),
            (None, _) => Ok(()),
        }
    }

    /// Opens the file `name` of `bytes` bytes. On S3 nothing is sent, as the size is known.
    pub(crate) fn open(&self, name: &str, bytes: u64) -> Result<Arc<dyn ReadAt>, IndexError> {
        Ok(match self {
            Self::Directory(path) => Arc::new(LocalFile::open(path.join(name))?),
            Self::S3 { prefix, config } => {
                Arc::new(S3Object::open(config, &format!("{prefix}{name}"), bytes)?)
            }
        })
    }
}

/// Removes the pin file at `path` and its ETag file. Open sources keep reading what they mapped.
pub(crate) fn remove_pin(path: &Path) {
    for file in [path.to_path_buf(), etag_file(path)] {
        let _ = std::fs::remove_file(file);
    }
}

/// How a shard opens its files.
#[derive(Clone, Default)]
pub struct ShardOptions {
    /// Where the packed sequences lie when not in the shard's own `packed/`, for example a
    /// directory that several shards share.
    pub packed: Option<Location>,
    /// A block cache that every file except the sketch is read through. With it, each entry's
    /// packed sequence stays open after its first read, so its blocks are found again.
    pub cache: Option<Arc<BlockCache>>,
    /// A cache directory warmed by `jam index warm`; the ranges its plan pins are read from
    /// there. Without a plan for the shard nothing is pinned.
    pub pins: Option<PathBuf>,
}

/// The file sources of one open shard.
pub(crate) struct ShardFiles {
    location: Location,
    packed: Location,
    cache: Option<Arc<BlockCache>>,
    /// The pins of the plan with their local files, by file name relative to the shard.
    pins: HashMap<String, Vec<(Pin, PathBuf)>>,
    /// Packed sequences kept open per entry when they are cached or remote; empty otherwise.
    opened: Vec<OnceLock<Arc<dyn ReadAt>>>,
}

impl ShardFiles {
    /// The files of the shard at `location` with `entries` entries and manifest digest `digest`,
    /// opened as `options` say.
    pub(crate) fn new(
        location: &Location,
        options: &ShardOptions,
        entries: usize,
        digest: &str,
    ) -> Result<Self, IndexError> {
        let mut pins: HashMap<String, Vec<(Pin, PathBuf)>> = HashMap::new();
        if let Some(cache) = &options.pins
            && let Some(plan) = Plan::read(cache, digest)?
        {
            for pin in &plan.pins {
                let path = plan.pin_path(cache, pin);
                pins.entry(pin.file.clone())
                    .or_default()
                    .push((pin.clone(), path));
            }
        }
        let packed = (options.packed.clone()).unwrap_or_else(|| location.join(PACKED_DIRECTORY));
        let keep = options.cache.is_some() || matches!(packed, Location::S3 { .. });
        let opened = if keep {
            (0..entries).map(|_| OnceLock::new()).collect()
        } else {
            Vec::new()
        };
        Ok(Self {
            location: location.clone(),
            packed,
            cache: options.cache.clone(),
            pins,
            opened,
        })
    }

    /// The source of a shard file: its pinned range locally, the rest through the cache if
    /// there is one.
    pub(crate) fn open(&self, record: &FileRecord) -> Result<Arc<dyn ReadAt>, IndexError> {
        self.open_at(&self.location, &record.path, record, None)
    }

    /// The file `name` (relative to the shard) of `record` at `location`, with its pins; cached
    /// in blocks of `block_bytes`, or of the cache's block size.
    fn open_at(
        &self,
        location: &Location,
        name: &str,
        record: &FileRecord,
        block_bytes: Option<u64>,
    ) -> Result<Arc<dyn ReadAt>, IndexError> {
        let object = location.object(&record.path, record.bytes)?;
        let hook = self.forget(name, record, block_bytes);
        let object = object.map(|object| Arc::new(object.on_change(hook)));
        let source: Arc<dyn ReadAt> = match &object {
            Some(object) => object.clone(),
            None => location.open(&record.path, record.bytes)?,
        };
        let mut source = self.cached(source, record, block_bytes)?;
        // Each pin takes its range from the source before it, so several parts of one file can
        // be pinned. A pin of the whole file is wrapped too, so its reads are counted.
        for (pin, path) in self.pins.get(name).into_iter().flatten() {
            if !path.exists() || !self.pin_current(object.as_deref(), path) {
                continue;
            }
            let file = LocalFile::open(path)?;
            if file.size() != pin.bytes() || pin.end > record.bytes {
                return Err(IndexError::format(format!(
                    "the pinned file {} differs from its plan",
                    path.display()
                )));
            }
            let local: Arc<dyn ReadAt> = match pin.place {
                Place::Ram => Arc::new(MappedFile::read_from(&file)?),
                Place::Ssd => Arc::new(file),
            };
            source = Arc::new(PinnedSource::new(local, pin.start..pin.end, source));
        }
        Ok(match object {
            Some(object) => Arc::new(Checked::new(source, object)),
            None => source,
        })
    }

    /// Whether the pin file at `path` holds bytes of the version that `object` has. Its ETag is
    /// taken as the known one if none is known. A pin without an ETag, or of another version, is
    /// not used and removed. Pins of a directory have no ETag.
    fn pin_current(&self, object: Option<&S3Object>, path: &Path) -> bool {
        let Some(object) = object else {
            return true;
        };
        let etag = std::fs::read_to_string(etag_file(path)).ok();
        if etag.is_some_and(|etag| object.adopt(&etag)) {
            return true;
        }
        remove_pin(path);
        false
    }

    /// What to do when the file `name` of `record` is found replaced: drop its pins by their
    /// exact paths, and its blocks of the old version.
    fn forget(&self, name: &str, record: &FileRecord, block_bytes: Option<u64>) -> ChangeHook {
        let cache = self.cache.clone();
        let (digest, size) = (record.sha256.clone(), record.bytes);
        let pins: Vec<PathBuf> = (self.pins.get(name).into_iter().flatten())
            .map(|(_, path)| path.clone())
            .collect();
        Arc::new(move |version| {
            pins.iter().for_each(|path| remove_pin(path));
            if let Some(cache) = &cache {
                cache.forget(&digest, block_bytes, size, version);
            }
        })
    }

    /// The source of a shard file without the cache, for example to verify it.
    pub(crate) fn open_uncached(&self, record: &FileRecord) -> Result<Arc<dyn ReadAt>, IndexError> {
        self.location.open(&record.path, record.bytes)
    }

    /// The local path of a packed sequence file, or `None` when it lies on S3.
    pub(crate) fn packed_path(&self, record: &FileRecord) -> Option<PathBuf> {
        self.packed.local_path(&record.path)
    }

    /// The source of entry `entry`'s packed sequence file, through the cache if there is one.
    pub(crate) fn open_packed(
        &self,
        entry: usize,
        record: &FileRecord,
    ) -> Result<Arc<dyn ReadAt>, IndexError> {
        let name = format!("{PACKED_DIRECTORY}/{}", record.path);
        // Packed files are cached in their checked blocks, so a window reads few bytes.
        let blocks = Some(BLOCK_BYTES);
        let Some(slot) = self.opened.get(entry) else {
            return self.open_at(&self.packed, &name, record, blocks);
        };
        if let Some(source) = slot.get() {
            return Ok(Arc::clone(source));
        }
        let source = self.open_at(&self.packed, &name, record, blocks)?;
        // Two threads may open the same entry at once; both get the source that was kept.
        Ok(Arc::clone(slot.get_or_init(|| source)))
    }

    /// The source of a packed sequence file without the cache.
    pub(crate) fn open_packed_uncached(
        &self,
        record: &FileRecord,
    ) -> Result<Arc<dyn ReadAt>, IndexError> {
        self.packed.open(&record.path, record.bytes)
    }

    /// `source` read through the cache, if there is one, under the digest of its record.
    fn cached(
        &self,
        source: Arc<dyn ReadAt>,
        record: &FileRecord,
        block_bytes: Option<u64>,
    ) -> Result<Arc<dyn ReadAt>, IndexError> {
        Ok(match &self.cache {
            Some(cache) => Arc::new(cache.wrap(source, &record.sha256, block_bytes)?),
            None => source,
        })
    }
}
