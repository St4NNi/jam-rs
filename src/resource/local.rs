//! Local-path and `file://` range readers.

use super::cache::BlockCache;
use super::error::{io_error, validate_range};
use super::metrics::MetricsCounter;
use super::{
    ByteRange, CacheIdentity, RangeBytes, RangeReader, ResourceError, ResourceLocator,
    ResourceMetadata, ResourceOpenOptions, ResourceResult, ResourceScheme,
};
use memmap2::{Mmap, MmapOptions};
use std::collections::HashSet;
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
    mmap: Option<Mmap>,
    resident_pages: Mutex<HashSet<u64>>,
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
        let file = File::open(&path).map_err(|error| io_error(&locator, error))?;
        let file_size = file
            .metadata()
            .map_err(|error| io_error(&locator, error))?
            .len();
        let mmap = if file_size == 0 {
            None
        } else {
            // SAFETY: the read-only mapping owns its virtual-memory region for
            // the lifetime of this resource and the file is never mutated by
            // the resource implementation.
            Some(
                unsafe { MmapOptions::new().map(&file) }
                    .map_err(|error| io_error(&locator, error))?,
            )
        };
        let cache = BlockCache::new(options.cache_block_bytes, options.max_cache_bytes)?;
        let metrics = Arc::new(MetricsCounter::default());
        metrics.mapped_bytes(file_size);
        Ok(Self {
            locator,
            path,
            options,
            cache: Arc::new(cache),
            metrics,
            metadata: Mutex::new(None),
            mmap,
            resident_pages: Mutex::new(HashSet::new()),
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
            version: self.options.expected_sha256.as_ref().map_or_else(
                || {
                    metadata
                        .etag
                        .as_deref()
                        .or(metadata.last_modified.as_deref())
                        .unwrap_or("unknown")
                        .to_string()
                },
                |checksum| format!("sha256:{checksum}"),
            ),
            size: metadata.size,
        }
    }

    fn read_uncached(&self, range: ByteRange) -> ResourceResult<Vec<u8>> {
        self.metrics.requested_bytes(range.length);
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
        self.metrics
            .returned_bytes(u64::try_from(bytes.len()).unwrap_or(u64::MAX));
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
        self.cache.prepare_with_metrics(identity, &self.metrics)?;
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
            let block = if let Some(bytes) =
                self.cache
                    .get_with_metrics(identity, block_offset, &self.metrics)?
            {
                let overlap_start = range.offset.max(block_offset);
                let overlap_end = requested_end.min(block_end);
                self.metrics.cache_hit();
                self.metrics
                    .cache_bytes(overlap_end.saturating_sub(overlap_start));
                bytes
            } else {
                self.metrics.cache_miss();
                let bytes = self.read_uncached(block_range)?;
                self.cache.insert_with_metrics(
                    identity,
                    block_offset,
                    bytes.clone(),
                    &self.metrics,
                )?;
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
        self.metrics
            .decoded_bytes(u64::try_from(output.len()).unwrap_or(u64::MAX));
        Ok(output)
    }

    fn record_resident_pages(&self, range: ByteRange, file_size: u64) -> ResourceResult<()> {
        if range.is_empty() {
            return Ok(());
        }
        const PAGE_BYTES: u64 = 4096;
        let end = range.end()?;
        let first = range.offset / PAGE_BYTES;
        let last = (end - 1) / PAGE_BYTES;
        let mut pages = self.resident_pages.lock().map_err(|_| ResourceError::Io {
            locator: self.locator.redacted(),
            message: "mapped-page counter lock poisoned".to_string(),
        })?;
        let mut added = 0u64;
        for page in first..=last {
            if pages.insert(page) {
                let start = page.saturating_mul(PAGE_BYTES);
                added = added.saturating_add(file_size.saturating_sub(start).min(PAGE_BYTES));
            }
        }
        self.metrics.resident_bytes(added);
        Ok(())
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

    fn read_range_bytes(&self, range: ByteRange) -> ResourceResult<RangeBytes<'_>> {
        self.metrics.range_request();
        let metadata = self.metadata()?;
        validate_range(&self.locator, range, metadata.size)?;
        if range.is_empty() {
            return Ok(RangeBytes::Borrowed(&[]));
        }
        let mmap = self.mmap.as_ref().ok_or_else(|| ResourceError::Io {
            locator: self.locator.redacted(),
            message: "non-empty local resource has no memory map".to_string(),
        })?;
        let start = usize::try_from(range.offset).map_err(|_| ResourceError::Io {
            locator: self.locator.redacted(),
            message: "mapped range offset does not fit this platform".to_string(),
        })?;
        let end = usize::try_from(range.end()?).map_err(|_| ResourceError::Io {
            locator: self.locator.redacted(),
            message: "mapped range end does not fit this platform".to_string(),
        })?;
        let bytes = mmap.get(start..end).ok_or_else(|| ResourceError::Io {
            locator: self.locator.redacted(),
            message: "mapped range exceeds the local file".to_string(),
        })?;
        self.metrics.requested_bytes(range.length);
        self.metrics.returned_bytes(range.length);
        self.metrics.decoded_bytes(range.length);
        self.record_resident_pages(range, metadata.size)?;
        Ok(RangeBytes::Borrowed(bytes))
    }

    fn stream(&self) -> ResourceResult<Box<dyn Read + Send>> {
        self.metrics.stream_request();
        let file = File::open(&self.path).map_err(|error| io_error(&self.locator, error))?;
        Ok(Box::new(LocalReader {
            reader: BufReader::with_capacity(1024 * 1024, file),
            metrics: Arc::clone(&self.metrics),
        }))
    }

    fn metrics(&self) -> super::ResourceMetrics {
        self.metrics.snapshot()
    }
}

struct LocalReader {
    reader: BufReader<File>,
    metrics: Arc<MetricsCounter>,
}

impl Read for LocalReader {
    fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
        let count = self.reader.read(bytes)?;
        let count = u64::try_from(count).unwrap_or(u64::MAX);
        self.metrics.returned_bytes(count);
        self.metrics.decoded_bytes(count);
        Ok(usize::try_from(count).unwrap_or(usize::MAX))
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
