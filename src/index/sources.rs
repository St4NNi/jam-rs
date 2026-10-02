//! Where an open index reads its files: the index file and the sequence files next to it, each a
//! local file or an S3 object, optionally through a block cache shared by all files.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::manifest::{FileRecord, Manifest};
use super::packed::BLOCK_BYTES;
use super::plan::{Pin, Place, Plan, etag_file};
use super::sequence_file::INDEX_EXTENSION;
use super::{AtPath, IndexError};
use crate::storage::{
    BlockCache, ChangeHook, Checked, LocalFile, MappedFile, PinnedSource, ReadAt, S3Config,
    S3Object,
};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

/// Where an index file lies: a local path or an S3 object. Its sequence files are named
/// relative to its directory.
#[derive(Clone)]
pub enum Location {
    /// A local file.
    File(PathBuf),
    /// An S3 object.
    S3 {
        /// The object, `s3://bucket/key`.
        location: String,
        /// The store, the credentials and the limits.
        config: Box<S3Config>,
    },
}

impl Location {
    /// The S3 object `s3://bucket/key` on the store of `config`.
    pub fn s3(location: &str, config: &S3Config) -> Result<Self, IndexError> {
        let parts = location
            .strip_prefix("s3://")
            .and_then(|rest| rest.split_once('/'));
        if !parts.is_some_and(|(bucket, key)| !bucket.is_empty() && !key.ends_with('/')) {
            return Err(IndexError::Input(format!(
                "{location} is not an S3 object like s3://bucket/prefix/NAME.jam"
            )));
        }
        let config = Box::new(config.clone());
        let location = location.to_string();
        Ok(Self::S3 { location, config })
    }

    /// The file name without `.jam`, which names the index in records.
    pub fn name(&self) -> String {
        let file = self.file_name();
        let suffix = format!(".{INDEX_EXTENSION}");
        file.strip_suffix(&suffix).unwrap_or(&file).to_string()
    }

    /// The file name of the index file.
    fn file_name(&self) -> String {
        let text = self.to_string();
        text.rsplit('/').next().unwrap_or_default().to_string()
    }

    /// The S3 location of `name`, a path relative to the directory of the index file; `..`
    /// leaves a directory of the key.
    fn s3_location(location: &str, name: &str) -> String {
        let rest = location.strip_prefix("s3://").unwrap_or(location);
        let (bucket, key) = rest.split_once('/').unwrap_or((rest, ""));
        let mut parts: Vec<&str> = key.split('/').collect();
        parts.pop();
        for part in name.split('/') {
            match part {
                "" | "." => {}
                ".." => drop(parts.pop()),
                part => parts.push(part),
            }
        }
        format!("s3://{bucket}/{}", parts.join("/"))
    }

    /// The local path of the file `name`, relative to the directory of the index file, or
    /// `None` when it lies on S3.
    pub(crate) fn local_path(&self, name: &str) -> Option<PathBuf> {
        match self {
            Self::File(path) => Some(path.with_file_name("").join(name)),
            Self::S3 { .. } => None,
        }
    }

    /// The head of the index file at this location: its header and manifest. On S3 the first
    /// read learns the ETag of the object, or sends the known one.
    pub fn head(&self) -> Result<IndexHead, IndexError> {
        // The size is not known before the header is read. The reads of an S3 object are bound
        // to one ETag, so they also check the object for this process.
        let file = self.file_name();
        let (manifest, digest, bytes) = match self.object(&file, u64::MAX)? {
            Some(object) => {
                let head = Manifest::read_from(&object)?;
                object.pass_check();
                head
            }
            None => Manifest::read_from(self.open(&file, 0)?.as_ref())?,
        };
        Ok(IndexHead {
            manifest,
            digest,
            bytes,
            file,
        })
    }

    /// The S3 object of the file `name` of `bytes` bytes, or `None` for a local file. Nothing
    /// is sent, as the size is known.
    fn object(&self, name: &str, bytes: u64) -> Result<Option<S3Object>, IndexError> {
        Ok(match self {
            Self::File(_) => None,
            Self::S3 { location, config } => Some(S3Object::open(
                config,
                &Self::s3_location(location, name),
                bytes,
            )?),
        })
    }

    /// How often an S3 object of the store was found replaced; 0 for a local file.
    pub(crate) fn changes(&self) -> u64 {
        match self {
            Self::File(_) => 0,
            Self::S3 { config, .. } => config.changes(),
        }
    }

    /// Forgets the ETags of the S3 objects of `names`.
    fn forget_etags(&self, names: &[String]) {
        if let Self::S3 { location, config } = self {
            for name in names {
                config.forget_etag(&Self::s3_location(location, name));
            }
        }
    }

    /// Keeps the ETags of S3 objects in `directory` too, so a later run finds them.
    pub(crate) fn keep_etags(&self, directory: &Path) -> Result<(), IndexError> {
        if let Self::S3 { config, .. } = self {
            config.keep_etags(directory)?;
        }
        Ok(())
    }

    /// The ETag that the file `name` has on S3 now, learned with one byte; `None` for a local
    /// file.
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
    /// ETag fails too. Nothing is checked for a local file.
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
        Ok(match self.object(name, bytes)? {
            Some(object) => Arc::new(object),
            None => Arc::new(LocalFile::open(self.local_path(name).unwrap_or_default())?),
        })
    }
}

impl std::fmt::Display for Location {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::File(path) => write!(formatter, "{}", path.display()),
            Self::S3 { location, .. } => formatter.write_str(location),
        }
    }
}

/// The head of an index file: its manifest, the SHA-256 digest of the manifest, which names this
/// version of the index, the size of the file and its name.
#[derive(Clone, Debug)]
pub struct IndexHead {
    pub manifest: Manifest,
    pub digest: String,
    pub bytes: u64,
    pub file: String,
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
    /// A block cache that every file except the sketch is read through. With it, each sequence
    /// file stays open after its first read, so its blocks are found again.
    pub cache: Option<Arc<BlockCache>>,
    /// A cache directory warmed by `jam index warm`; the ranges its plan pins are read from
    /// there. Without a plan for the shard nothing is pinned.
    pub pins: Option<PathBuf>,
}

/// The file sources of one open shard.
pub(crate) struct ShardFiles {
    location: Location,
    /// The sequence files.
    files: Vec<FileRecord>,
    cache: Option<Arc<BlockCache>>,
    /// The pins of the plan with their local files, by file name relative to the shard.
    pins: HashMap<String, Vec<(Pin, PathBuf)>>,
    /// The seeds file as a record: its size, and the manifest digest as its key.
    seeds: FileRecord,
    /// The sequence files, each opened at its first use.
    opened: Vec<OnceLock<Arc<dyn ReadAt>>>,
    /// The replaced objects counted when the files were opened.
    changes: u64,
}

impl ShardFiles {
    /// The files of the shard at `location` with head `head`, opened as `options` say.
    pub(crate) fn new(
        location: &Location,
        options: &ShardOptions,
        head: &IndexHead,
    ) -> Result<Self, IndexError> {
        let mut pins: HashMap<String, Vec<(Pin, PathBuf)>> = HashMap::new();
        if let Some(cache) = &options.pins
            && let Some(plan) = Plan::read(cache, &head.digest)?
        {
            for pin in &plan.pins {
                let path = plan.pin_path(cache, pin);
                pins.entry(pin.file.clone())
                    .or_default()
                    .push((pin.clone(), path));
            }
        }
        let files = head.manifest.sequence_files.clone();
        Ok(Self {
            location: location.clone(),
            opened: files.iter().map(|_| OnceLock::new()).collect(),
            files,
            cache: options.cache.clone(),
            pins,
            seeds: FileRecord {
                path: head.file.clone(),
                bytes: head.bytes,
                sha256: head.digest.clone(),
            },
            changes: location.changes(),
        })
    }

    /// Whether an object was found replaced since the files were opened. Their reads then fail.
    pub(crate) fn is_stale(&self) -> bool {
        self.location.changes() != self.changes
    }

    /// The source of the index file: its pinned ranges locally, the rest through the cache if
    /// there is one, cached under the manifest digest.
    pub(crate) fn open_seeds(&self) -> Result<Arc<dyn ReadAt>, IndexError> {
        self.open_at(&self.location, &self.seeds.path, &self.seeds, None)
    }

    /// The file `name` (relative to the index file) of `record` at `location`, with its pins;
    /// cached in blocks of `block_bytes`, or of the cache's block size.
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

    /// What to do when the file `name` of `record` is found replaced: forget the ETags of all
    /// files of the index, since an upload replaces all of them, drop the pins of the file by
    /// their exact paths, and drop its blocks of the old version.
    fn forget(&self, name: &str, record: &FileRecord, block_bytes: Option<u64>) -> ChangeHook {
        let cache = self.cache.clone();
        let (digest, size) = (record.sha256.clone(), record.bytes);
        let pins: Vec<PathBuf> = (self.pins.get(name).into_iter().flatten())
            .map(|(_, path)| path.clone())
            .collect();
        let location = self.location.clone();
        let mut names = vec![self.seeds.path.clone()];
        names.extend(self.files.iter().map(|file| file.path.clone()));
        Arc::new(move |version| {
            location.forget_etags(&names);
            pins.iter().for_each(|path| remove_pin(path));
            if let Some(cache) = &cache {
                cache.forget(&digest, block_bytes, size, version);
            }
        })
    }

    /// The index file without the cache, for example to verify it.
    pub(crate) fn open_seeds_uncached(&self) -> Result<Arc<dyn ReadAt>, IndexError> {
        self.location.open(&self.seeds.path, self.seeds.bytes)
    }

    /// The local path of sequence file `file`, or `None` when it lies on S3.
    pub(crate) fn sequence_path(&self, file: usize) -> Option<PathBuf> {
        self.location.local_path(&self.files.get(file)?.path)
    }

    /// The source of sequence file `file`, opened once, through the cache if there is one. It
    /// is cached in the checked blocks of the packed sequences, so a window reads few bytes.
    pub(crate) fn sequence(&self, file: usize) -> Result<Arc<dyn ReadAt>, IndexError> {
        let (slot, record) = (self.opened.get(file).zip(self.files.get(file)))
            .ok_or_else(|| IndexError::format("the manifest names no such sequence file"))?;
        if let Some(source) = slot.get() {
            return Ok(Arc::clone(source));
        }
        let source = self.open_at(&self.location, &record.path, record, Some(BLOCK_BYTES))?;
        if source.size() != record.bytes {
            let path = &record.path;
            return Err(IndexError::format(format!(
                "{path} differs from the manifest"
            )));
        }
        // Two threads may open the same file at once; both get the source that was kept.
        Ok(Arc::clone(slot.get_or_init(|| source)))
    }

    /// Sequence file `file` without the cache, with its record.
    pub(crate) fn open_uncached(
        &self,
        file: usize,
    ) -> Result<(Arc<dyn ReadAt>, &FileRecord), IndexError> {
        let record = (self.files.get(file))
            .ok_or_else(|| IndexError::format("the manifest names no such sequence file"))?;
        Ok((self.location.open(&record.path, record.bytes)?, record))
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
