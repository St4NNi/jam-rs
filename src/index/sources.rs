//! Where an open shard reads its files: a local directory or an S3 prefix, optionally through a
//! block cache shared by all files.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::manifest::{FileRecord, MANIFEST_FILE, Manifest, PACKED_DIRECTORY, sha256_hex};
use super::packed::BLOCK_BYTES;
use super::plan::{Pin, Place, Plan};
use super::{AtPath, IndexError};
use crate::storage::{BlockCache, LocalFile, MappedFile, PinnedSource, ReadAt, S3Config, S3Object};
use rayon::prelude::*;
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

    /// Opens the file of `record` under `name`. On S3 nothing is sent, as size and ETag come
    /// from the record; a record without an ETag cannot be read from S3, since its reads could
    /// not be bound to the object that the manifest describes.
    pub(crate) fn open(
        &self,
        name: &str,
        record: &FileRecord,
    ) -> Result<Arc<dyn ReadAt>, IndexError> {
        Ok(match self {
            Self::Directory(path) => Arc::new(LocalFile::open(path.join(name))?),
            Self::S3 { prefix, config } => {
                let etag = record.etag.as_deref().ok_or_else(|| {
                    IndexError::Input(format!(
                        "the manifest has no ETag for {name}: run jam index etags after the upload"
                    ))
                })?;
                let location = format!("{prefix}{name}");
                Arc::new(S3Object::open(config, &location, record.bytes, etag)?)
            }
        })
    }
}

/// Records in the manifest of the local shard `directory` the ETag that each of its files has
/// now under the S3 prefix `uploaded`, and rewrites the manifest. Run it right after the upload:
/// it trusts that the objects are the uploaded files. Reads one byte of each object.
pub fn record_etags(directory: &Path, uploaded: &Location) -> Result<Manifest, IndexError> {
    let Location::S3 { prefix, config } = uploaded else {
        return Err(IndexError::Input("ETags come from an s3:// prefix".into()));
    };
    let mut manifest = Manifest::read(directory)?;
    let etag = |name: &str| S3Object::etag(config, &format!("{prefix}{name}"));
    manifest.seeds.etag = Some(etag(&manifest.seeds.path)?);
    (manifest.entries.par_iter_mut())
        .filter_map(|entry| entry.packed.as_mut())
        .try_for_each(|packed| {
            packed.etag = Some(etag(&format!("{PACKED_DIRECTORY}/{}", packed.path))?);
            Ok::<_, IndexError>(())
        })?;
    let bytes = serde_json::to_vec_pretty(&manifest)
        .map_err(|error| IndexError::format(format!("manifest: {error}")))?;
    let path = directory.join(MANIFEST_FILE);
    let partial = directory.join(format!("{MANIFEST_FILE}.partial"));
    std::fs::write(&partial, bytes).at(&partial)?;
    std::fs::rename(&partial, &path).at(&path)?;
    Ok(manifest)
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
        let source = location.open(&record.path, record)?;
        let mut source = self.cached(source, record, block_bytes)?;
        // Each pin takes its range from the source before it, so several parts of one file can
        // be pinned. A pin of the whole file is wrapped too, so its reads are counted.
        for (pin, path) in self.pins.get(name).into_iter().flatten() {
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
        Ok(source)
    }

    /// The source of a shard file without the cache, for example to verify it.
    pub(crate) fn open_uncached(&self, record: &FileRecord) -> Result<Arc<dyn ReadAt>, IndexError> {
        self.location.open(&record.path, record)
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
        self.packed.open(&record.path, record)
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
