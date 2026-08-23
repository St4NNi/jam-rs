//! Local-path and `file://` range readers.

use super::cache::BlockCache;
use super::error::{io_error, validate_range};
use super::metrics::MetricsCounter;
use super::{
    ByteRange, CacheIdentity, RangeReader, ResourceError, ResourceLocator, ResourceMetadata,
    ResourceOpenOptions, ResourceResult, ResourceScheme,
};
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

/// A local file exposed through the shared checked range interface.
pub struct LocalResource {
    locator: ResourceLocator,
    path: PathBuf,
    options: ResourceOpenOptions,
    cache: Arc<BlockCache>,
    metrics: Arc<MetricsCounter>,
    metadata: Mutex<Option<ResourceMetadata>>,
}

impl std::fmt::Debug for LocalResource {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("LocalResource")
            .field("locator", &self.locator.redacted())
            .field("path", &self.path)
            .field("options", &self.options)
            .finish_non_exhaustive()
    }
}

impl LocalResource {
    /// Open a local path or `file://` locator.
    pub fn open(input: impl AsRef<str>, options: ResourceOpenOptions) -> ResourceResult<Self> {
        let locator = ResourceLocator::parse(input.as_ref())?;
        Self::from_locator(locator, options)
    }

    pub fn from_locator(
        locator: ResourceLocator,
        options: ResourceOpenOptions,
    ) -> ResourceResult<Self> {
        if !matches!(
            locator.scheme(),
            ResourceScheme::Local | ResourceScheme::File
        ) {
            return Err(ResourceError::UnsupportedScheme(format!(
                "{:?}",
                locator.scheme()
            )));
        }
        let path = path_from_locator(&locator)?;
        // Opening is intentionally lazy with respect to file contents, but a
        // metadata call here gives callers an immediate and useful path error.
        std::fs::metadata(&path).map_err(|error| io_error(&locator, error))?;
        let cache = BlockCache::new(options.cache_block_bytes, options.max_cache_bytes)?;
        Ok(Self {
            locator,
            path,
            options,
            cache: Arc::new(cache),
            metrics: Arc::new(MetricsCounter::default()),
            metadata: Mutex::new(None),
        })
    }

    pub fn from_path(path: impl AsRef<Path>, options: ResourceOpenOptions) -> ResourceResult<Self> {
        let path = path.as_ref();
        Self::open(path.to_string_lossy(), options)
    }

    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    #[must_use]
    pub fn cache(&self) -> &BlockCache {
        &self.cache
    }

    fn current_metadata(&self) -> ResourceResult<ResourceMetadata> {
        let file_metadata =
            std::fs::metadata(&self.path).map_err(|error| io_error(&self.locator, error))?;
        let modified = file_metadata
            .modified()
            .ok()
            .and_then(|time| time.duration_since(std::time::UNIX_EPOCH).ok())
            .map(|duration| duration.as_nanos().to_string());
        Ok(ResourceMetadata {
            size: file_metadata.len(),
            etag: modified.clone(),
            last_modified: modified,
            accepts_ranges: true,
        })
    }

    fn identity(&self, metadata: &ResourceMetadata) -> CacheIdentity {
        CacheIdentity {
            redacted_locator: self.locator.redacted(),
            version: metadata
                .etag
                .as_deref()
                .or(metadata.last_modified.as_deref())
                .unwrap_or("unknown")
                .to_string(),
            size: metadata.size,
        }
    }

    fn read_uncached(&self, range: ByteRange) -> ResourceResult<Vec<u8>> {
        let length = usize::try_from(range.length).map_err(|_| ResourceError::Io {
            locator: self.locator.redacted(),
            message: "requested range is too large for this platform".to_string(),
        })?;
        let mut file = File::open(&self.path).map_err(|error| io_error(&self.locator, error))?;
        file.seek(SeekFrom::Start(range.offset))
            .map_err(|error| io_error(&self.locator, error))?;
        let mut bytes = vec![0; length];
        file.read_exact(&mut bytes)
            .map_err(|error| io_error(&self.locator, error))?;
        Ok(bytes)
    }

    fn read_cached(
        &self,
        range: ByteRange,
        metadata: &ResourceMetadata,
        identity: &CacheIdentity,
    ) -> ResourceResult<Vec<u8>> {
        if range.is_empty() {
            return Ok(Vec::new());
        }
        self.cache.prepare(identity)?;
        let mut output =
            Vec::with_capacity(
                usize::try_from(range.length).map_err(|_| ResourceError::Io {
                    locator: self.locator.redacted(),
                    message: "requested range is too large for this platform".to_string(),
                })?,
            );
        let requested_end = range.end()?;
        for block_offset in self.cache.block_offsets(range)? {
            let block_end = block_offset
                .saturating_add(self.cache.block_bytes())
                .min(metadata.size);
            let block_range = ByteRange::new(block_offset, block_end.saturating_sub(block_offset))?;
            let block = if let Some(bytes) = self.cache.get(identity, block_offset)? {
                let overlap_start = range.offset.max(block_offset);
                let overlap_end = requested_end.min(block_end);
                self.metrics.cache_hit();
                self.metrics
                    .cache_bytes(overlap_end.saturating_sub(overlap_start));
                bytes
            } else {
                let bytes = self.read_uncached(block_range)?;
                self.cache.insert(identity, block_offset, bytes.clone())?;
                bytes
            };
            let copy_start = range.offset.max(block_offset).saturating_sub(block_offset);
            let copy_end = requested_end.min(block_end).saturating_sub(block_offset);
            let copy_start = usize::try_from(copy_start).map_err(|_| ResourceError::Io {
                locator: self.locator.redacted(),
                message: "cache offset does not fit in platform usize".to_string(),
            })?;
            let copy_end = usize::try_from(copy_end).map_err(|_| ResourceError::Io {
                locator: self.locator.redacted(),
                message: "cache offset does not fit in platform usize".to_string(),
            })?;
            if copy_start > copy_end || copy_end > block.len() {
                return Err(ResourceError::Io {
                    locator: self.locator.redacted(),
                    message: "cached block has an invalid length".to_string(),
                });
            }
            output.extend_from_slice(&block[copy_start..copy_end]);
        }
        Ok(output)
    }
}

impl RangeReader for LocalResource {
    fn locator(&self) -> &ResourceLocator {
        &self.locator
    }

    fn metadata(&self) -> ResourceResult<ResourceMetadata> {
        self.metrics.metadata_request();
        let metadata = self.current_metadata()?;
        let mut cached = self.metadata.lock().map_err(|_| ResourceError::Io {
            locator: self.locator.redacted(),
            message: "metadata lock poisoned".to_string(),
        })?;
        *cached = Some(metadata.clone());
        Ok(metadata)
    }

    fn read_range(&self, range: ByteRange) -> ResourceResult<Vec<u8>> {
        self.metrics.range_request();
        let metadata = self.metadata()?;
        validate_range(&self.locator, range, metadata.size)?;
        if range.is_empty() {
            return Ok(Vec::new());
        }
        let identity = self.identity(&metadata);
        self.read_cached(range, &metadata, &identity)
    }

    fn stream(&self) -> ResourceResult<Box<dyn Read + Send>> {
        self.metrics.stream_request();
        let file = File::open(&self.path).map_err(|error| io_error(&self.locator, error))?;
        Ok(Box::new(BufReader::with_capacity(1024 * 1024, file)))
    }

    fn metrics(&self) -> super::ResourceMetrics {
        self.metrics.snapshot()
    }
}

fn path_from_locator(locator: &ResourceLocator) -> ResourceResult<PathBuf> {
    match locator.scheme() {
        ResourceScheme::Local => Ok(PathBuf::from(locator.raw())),
        ResourceScheme::File => {
            let raw_path = locator.raw().strip_prefix("file://").unwrap_or_default();
            let raw_path = raw_path.split_once('#').map_or(raw_path, |(path, _)| path);
            let raw_path = raw_path.split_once('?').map_or(raw_path, |(path, _)| path);
            let path = percent_decode(raw_path)?;
            let path = if let Some(path) = path.strip_prefix("localhost/") {
                format!("/{path}")
            } else {
                path
            };
            Ok(PathBuf::from(path))
        }
        _ => Err(ResourceError::UnsupportedScheme(format!(
            "{:?}",
            locator.scheme()
        ))),
    }
}

fn percent_decode(value: &str) -> ResourceResult<String> {
    let bytes = value.as_bytes();
    let mut decoded = Vec::with_capacity(bytes.len());
    let mut index = 0;
    while index < bytes.len() {
        if bytes[index] != b'%' {
            decoded.push(bytes[index]);
            index += 1;
            continue;
        }
        if index + 2 >= bytes.len() {
            return Err(ResourceError::InvalidLocator(
                "file URL contains an incomplete percent escape".to_string(),
            ));
        }
        let high = hex_value(bytes[index + 1]).ok_or_else(|| {
            ResourceError::InvalidLocator("file URL contains an invalid percent escape".to_string())
        })?;
        let low = hex_value(bytes[index + 2]).ok_or_else(|| {
            ResourceError::InvalidLocator("file URL contains an invalid percent escape".to_string())
        })?;
        decoded.push((high << 4) | low);
        index += 3;
    }
    String::from_utf8(decoded)
        .map_err(|_| ResourceError::InvalidLocator("file URL is not valid UTF-8".to_string()))
}

fn hex_value(value: u8) -> Option<u8> {
    match value {
        b'0'..=b'9' => Some(value - b'0'),
        b'a'..=b'f' => Some(value - b'a' + 10),
        b'A'..=b'F' => Some(value - b'A' + 10),
        _ => None,
    }
}
