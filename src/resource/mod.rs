//! Uniform access to local and remote trace resources.

use serde::{Deserialize, Serialize};
use std::io::Read;
use std::path::PathBuf;
use thiserror::Error;

pub mod cache;
pub mod error;
pub mod local;
pub mod locator;
pub mod metrics;
pub mod object;
pub mod persistent_cache;
pub mod read_plan;
pub mod upload;

/// URI scheme understood by the trace resource layer.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ResourceScheme {
    Local,
    File,
    Http,
    Https,
    S3,
}

/// Parsed resource location. The unredacted value is never serialized.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct ResourceLocator {
    scheme: ResourceScheme,
    raw: String,
}

impl ResourceLocator {
    pub fn parse(input: &str) -> ResourceResult<Self> {
        locator::parse(input)
    }

    #[must_use]
    pub fn scheme(&self) -> ResourceScheme {
        self.scheme
    }

    #[must_use]
    pub fn redacted(&self) -> String {
        locator::redact(self)
    }

    pub(crate) fn from_parts(scheme: ResourceScheme, raw: String) -> Self {
        Self { scheme, raw }
    }

    pub(crate) fn raw(&self) -> &str {
        &self.raw
    }
}

impl std::str::FromStr for ResourceLocator {
    type Err = ResourceError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        Self::parse(input)
    }
}

impl Serialize for ResourceLocator {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&self.redacted())
    }
}

/// Half-open byte range `[offset, offset + length)`.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct ByteRange {
    pub offset: u64,
    pub length: u64,
}

impl ByteRange {
    pub fn new(offset: u64, length: u64) -> ResourceResult<Self> {
        offset
            .checked_add(length)
            .ok_or(ResourceError::RangeOverflow { offset, length })?;
        Ok(Self { offset, length })
    }

    pub fn end(self) -> ResourceResult<u64> {
        self.offset
            .checked_add(self.length)
            .ok_or(ResourceError::RangeOverflow {
                offset: self.offset,
                length: self.length,
            })
    }

    #[must_use]
    pub fn is_empty(self) -> bool {
        self.length == 0
    }
}

/// Stable metadata used for range validation and cache invalidation.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct ResourceMetadata {
    pub size: u64,
    pub etag: Option<String>,
    pub last_modified: Option<String>,
    pub accepts_ranges: bool,
}

/// Identity of cached bytes. `version` is an ETag, checksum, or last-modified token.
#[derive(Clone, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct CacheIdentity {
    pub redacted_locator: String,
    pub version: String,
    pub size: u64,
}

/// Per-resource counters. Implementations return snapshots, never live atomics.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct ResourceMetrics {
    pub metadata_requests: u64,
    pub head_requests: u64,
    pub get_requests: u64,
    pub range_requests: u64,
    pub stream_requests: u64,
    pub requested_bytes: u64,
    pub returned_bytes: u64,
    pub decoded_bytes: u64,
    pub mapped_bytes: u64,
    pub resident_bytes: u64,
    pub remote_bytes: u64,
    pub cache_bytes: u64,
    pub cache_hits: u64,
    pub cache_misses: u64,
    pub cache_evictions: u64,
    pub stale_cache_rejections: u64,
    pub retries: u64,
    pub full_object_fallbacks: u64,
    pub seed_buckets_read: u64,
    pub sequence_blocks_read: u64,
}

impl ResourceMetrics {
    #[must_use]
    pub fn saturating_add(self, other: Self) -> Self {
        Self {
            metadata_requests: self
                .metadata_requests
                .saturating_add(other.metadata_requests),
            head_requests: self.head_requests.saturating_add(other.head_requests),
            get_requests: self.get_requests.saturating_add(other.get_requests),
            range_requests: self.range_requests.saturating_add(other.range_requests),
            stream_requests: self.stream_requests.saturating_add(other.stream_requests),
            requested_bytes: self.requested_bytes.saturating_add(other.requested_bytes),
            returned_bytes: self.returned_bytes.saturating_add(other.returned_bytes),
            decoded_bytes: self.decoded_bytes.saturating_add(other.decoded_bytes),
            mapped_bytes: self.mapped_bytes.max(other.mapped_bytes),
            resident_bytes: self.resident_bytes.saturating_add(other.resident_bytes),
            remote_bytes: self.remote_bytes.saturating_add(other.remote_bytes),
            cache_bytes: self.cache_bytes.saturating_add(other.cache_bytes),
            cache_hits: self.cache_hits.saturating_add(other.cache_hits),
            cache_misses: self.cache_misses.saturating_add(other.cache_misses),
            cache_evictions: self.cache_evictions.saturating_add(other.cache_evictions),
            stale_cache_rejections: self
                .stale_cache_rejections
                .saturating_add(other.stale_cache_rejections),
            retries: self.retries.saturating_add(other.retries),
            full_object_fallbacks: self
                .full_object_fallbacks
                .saturating_add(other.full_object_fallbacks),
            seed_buckets_read: self
                .seed_buckets_read
                .saturating_add(other.seed_buckets_read),
            sequence_blocks_read: self
                .sequence_blocks_read
                .saturating_add(other.sequence_blocks_read),
        }
    }
}

/// Shared opening and cache limits for all resource schemes.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct ResourceOpenOptions {
    pub cache_dir: Option<PathBuf>,
    pub expected_sha256: Option<String>,
    pub cache_block_bytes: u64,
    pub max_cache_bytes: u64,
    pub request_timeout_seconds: u64,
    pub max_retries: u32,
    pub allow_full_download_fallback: bool,
}

impl Default for ResourceOpenOptions {
    fn default() -> Self {
        Self {
            cache_dir: None,
            expected_sha256: None,
            cache_block_bytes: 1024 * 1024,
            max_cache_bytes: 4 * 1024 * 1024 * 1024,
            request_timeout_seconds: 60,
            max_retries: 3,
            allow_full_download_fallback: true,
        }
    }
}

/// Synchronous range/stream contract used by archive and raw-sequence readers.
pub trait RangeReader: Send + Sync {
    fn locator(&self) -> &ResourceLocator;
    fn metadata(&self) -> ResourceResult<ResourceMetadata>;
    fn read_range(&self, range: ByteRange) -> ResourceResult<Vec<u8>>;
    fn read_range_bytes(&self, range: ByteRange) -> ResourceResult<RangeBytes<'_>> {
        self.read_range(range).map(RangeBytes::Owned)
    }
    fn stream(&self) -> ResourceResult<Box<dyn Read + Send>>;
    fn metrics(&self) -> ResourceMetrics;
}

pub enum RangeBytes<'a> {
    Borrowed(&'a [u8]),
    Owned(Vec<u8>),
}

impl AsRef<[u8]> for RangeBytes<'_> {
    fn as_ref(&self) -> &[u8] {
        match self {
            Self::Borrowed(bytes) => bytes,
            Self::Owned(bytes) => bytes,
        }
    }
}

impl<T: RangeReader + ?Sized> RangeReader for std::sync::Arc<T> {
    fn locator(&self) -> &ResourceLocator {
        (**self).locator()
    }

    fn metadata(&self) -> ResourceResult<ResourceMetadata> {
        (**self).metadata()
    }

    fn read_range(&self, range: ByteRange) -> ResourceResult<Vec<u8>> {
        (**self).read_range(range)
    }

    fn read_range_bytes(&self, range: ByteRange) -> ResourceResult<RangeBytes<'_>> {
        (**self).read_range_bytes(range)
    }

    fn stream(&self) -> ResourceResult<Box<dyn Read + Send>> {
        (**self).stream()
    }

    fn metrics(&self) -> ResourceMetrics {
        (**self).metrics()
    }
}

pub type ResourceResult<T> = Result<T, ResourceError>;

/// Errors are required to contain only redacted resource identifiers.
#[derive(Debug, Error)]
pub enum ResourceError {
    #[error("invalid resource locator: {0}")]
    InvalidLocator(String),
    #[error("unsupported resource scheme: {0}")]
    UnsupportedScheme(String),
    #[error("resource range overflows: offset={offset}, length={length}")]
    RangeOverflow { offset: u64, length: u64 },
    #[error("resource range exceeds size {size}: offset={offset}, length={length}")]
    RangeOutOfBounds { offset: u64, length: u64, size: u64 },
    #[error("resource changed while cached: {0}")]
    CacheIdentityChanged(String),
    #[error("resource does not support range reads: {0}")]
    RangeUnsupported(String),
    #[error("resource transport failed for {locator}: {message}")]
    Transport { locator: String, message: String },
    #[error("resource request returned HTTP status {status} for {locator}")]
    HttpStatus { locator: String, status: u16 },
    #[error("resource request timed out for {locator}")]
    Timeout { locator: String },
    #[error("resource I/O failed for {locator}: {message}")]
    Io { locator: String, message: String },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn byte_range_is_checked() {
        assert_eq!(ByteRange::new(4, 7).unwrap().end().unwrap(), 11);
        assert!(ByteRange::new(u64::MAX, 1).is_err());
    }

    #[test]
    fn locator_serialization_is_redacted() {
        let locator =
            ResourceLocator::parse("https://user:secret@example.org/a.jma?token=x").unwrap();
        let encoded = serde_json::to_string(&locator).unwrap();
        assert_eq!(encoded, "\"https://example.org/a.jma\"");
    }
}
