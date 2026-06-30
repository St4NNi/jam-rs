//! Where an open shard reads its packed sequences: its own `packed/` directory, another local
//! directory, or an S3 prefix, optionally through a block cache shared by all entries. Each
//! entry's file is found by the file name in the manifest.

use super::IndexError;
use super::manifest::{FileRecord, PACKED_DIRECTORY};
use crate::storage::MappedFile;
use crate::storage::ReadAt;
use crate::storage::cache::BlockCache;
use crate::storage::local::LocalFile;
use crate::storage::s3_store::{S3Config, S3Object};
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

/// Where a shard's packed sequences lie.
#[derive(Clone, Default)]
pub enum PackedLocation {
    /// The shard's own `packed/` directory.
    #[default]
    ShardDirectory,
    /// Another local directory with the packed sequence files, for example one that several
    /// shards share.
    Directory(PathBuf),
    /// An S3 prefix such as `s3://bucket/prefix/`; each file name is appended to it.
    S3 {
        /// The prefix, `s3://bucket/` followed by an optional key prefix.
        prefix: String,
        /// The store, the credentials and the limits.
        config: Box<S3Config>,
    },
}

/// How a shard opens its packed sequences.
#[derive(Clone, Default)]
pub struct PackedOptions {
    /// Where the packed sequences lie.
    pub location: PackedLocation,
    /// A block cache that every packed sequence is read through. With it, each entry's source
    /// stays open after its first read (one open file per entry for local files), so its blocks
    /// are found again.
    pub cache: Option<Arc<BlockCache>>,
}

/// The packed sequence sources of one open shard.
pub(crate) struct PackedSources {
    location: Location,
    cache: Option<Arc<BlockCache>>,
    /// Sources kept open per entry when they are cached or remote; empty otherwise.
    opened: Vec<OnceLock<Arc<dyn ReadAt>>>,
}

/// A resolved location.
enum Location {
    Local(PathBuf),
    S3 {
        prefix: String,
        config: Box<S3Config>,
    },
}

impl PackedSources {
    /// Resolves `options` for the shard in `directory` with `entries` entries.
    pub(crate) fn new(
        options: &PackedOptions,
        directory: &Path,
        entries: usize,
    ) -> Result<Self, IndexError> {
        let location = match &options.location {
            PackedLocation::ShardDirectory => Location::Local(directory.join(PACKED_DIRECTORY)),
            PackedLocation::Directory(path) => Location::Local(path.clone()),
            PackedLocation::S3 { prefix, config } => {
                let bucket = prefix
                    .strip_prefix("s3://")
                    .and_then(|rest| rest.split('/').next());
                if bucket.is_none_or(str::is_empty) {
                    return Err(IndexError::Input(format!(
                        "{prefix} is not an S3 prefix like s3://bucket/prefix/"
                    )));
                }
                let mut prefix = prefix.clone();
                if !prefix.ends_with('/') {
                    prefix.push('/');
                }
                let config = config.clone();
                Location::S3 { prefix, config }
            }
        };
        let keep = options.cache.is_some() || matches!(location, Location::S3 { .. });
        let opened = if keep {
            (0..entries).map(|_| OnceLock::new()).collect()
        } else {
            Vec::new()
        };
        Ok(Self {
            location,
            cache: options.cache.clone(),
            opened,
        })
    }

    /// The local path of a packed sequence file, or `None` when it lies on S3.
    pub(crate) fn local_path(&self, file: &FileRecord) -> Option<PathBuf> {
        match &self.location {
            Location::Local(directory) => Some(directory.join(&file.path)),
            Location::S3 { .. } => None,
        }
    }

    /// The source of entry `entry`'s packed sequence file, through the cache if there is one.
    pub(crate) fn open(
        &self,
        entry: usize,
        file: &FileRecord,
    ) -> Result<Arc<dyn ReadAt>, IndexError> {
        let Some(slot) = self.opened.get(entry) else {
            return self.open_uncached(file);
        };
        if let Some(source) = slot.get() {
            return Ok(Arc::clone(source));
        }
        let source: Arc<dyn ReadAt> = match (&self.cache, &self.location) {
            (Some(cache), Location::Local(directory)) => {
                Arc::new(cache.wrap(LocalFile::open(directory.join(&file.path))?))
            }
            (Some(cache), Location::S3 { prefix, config }) => {
                Arc::new(cache.wrap(open_s3(prefix, config, file)?))
            }
            (None, _) => self.open_uncached(file)?,
        };
        // Two threads may open the same entry at once; both get the source that was kept.
        Ok(Arc::clone(slot.get_or_init(|| source)))
    }

    /// The source of a packed sequence file without the cache, for example to verify it.
    pub(crate) fn open_uncached(&self, file: &FileRecord) -> Result<Arc<dyn ReadAt>, IndexError> {
        Ok(match &self.location {
            Location::Local(directory) => {
                Arc::new(MappedFile::open(&directory.join(&file.path), false)?)
            }
            Location::S3 { prefix, config } => Arc::new(open_s3(prefix, config, file)?),
        })
    }
}

fn open_s3(prefix: &str, config: &S3Config, file: &FileRecord) -> Result<S3Object, IndexError> {
    Ok(S3Object::open(
        config,
        &format!("{prefix}{}", file.path),
        Some(file.bytes),
    )?)
}
