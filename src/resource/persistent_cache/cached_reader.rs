//! Read-through [`RangeReader`] wrapper backed by the persistent page cache.
//!
//! The wrapper deliberately requires a caller-supplied page descriptor.  The
//! persistent cache key includes the expected page checksum, so a lookup can
//! only be safe when the immutable section index (or equivalent) can provide
//! that checksum before bytes are fetched.

use super::{CacheLookup, CacheMissReason, PageKey, PersistentCacheError, PersistentPageCache};
use crate::resource::metrics::MetricsCounter;
use crate::resource::{
    ByteRange, RangeBytes, RangeReader, ResourceError, ResourceMetadata, ResourceMetrics,
    ResourceResult,
};
use sha2::{Digest, Sha256};
use std::io::Read;
use std::sync::Arc;

/// Exact section/page information required for a safe persistent lookup.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PageDescriptor {
    pub section_id: u32,
    pub page_id: u64,
    pub page_checksum: [u8; 32],
}

/// Resolves a requested exact range to the immutable page descriptor that
/// contains it.  Returning `None` bypasses the persistent cache for that
/// request and delegates directly to the wrapped reader.
pub trait PageKeyResolver: Send + Sync {
    fn resolve(&self, range: ByteRange) -> Option<PageDescriptor>;
}

impl<F> PageKeyResolver for F
where
    F: Fn(ByteRange) -> Option<PageDescriptor> + Send + Sync,
{
    fn resolve(&self, range: ByteRange) -> Option<PageDescriptor> {
        self(range)
    }
}

/// Generic read-through wrapper for immutable exact pages.
pub struct CachedRangeReader<R, P> {
    inner: R,
    cache: Arc<PersistentPageCache>,
    object_checksum: [u8; 32],
    resolver: P,
    metrics: Arc<MetricsCounter>,
}

impl<R, P> std::fmt::Debug for CachedRangeReader<R, P>
where
    R: std::fmt::Debug,
    P: std::fmt::Debug,
{
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("CachedRangeReader")
            .field("inner", &self.inner)
            .field("cache", &self.cache)
            .field("object_checksum", &"[redacted]")
            .field("resolver", &self.resolver)
            .finish_non_exhaustive()
    }
}

impl<R: RangeReader, P> CachedRangeReader<R, P> {
    /// Wrap `inner` with an immutable object identity and page resolver.
    pub fn new(
        inner: R,
        cache: Arc<PersistentPageCache>,
        object_checksum: [u8; 32],
        resolver: P,
    ) -> Self {
        Self {
            inner,
            cache,
            object_checksum,
            resolver,
            metrics: Arc::new(MetricsCounter::default()),
        }
    }

    #[must_use]
    pub fn cache(&self) -> &Arc<PersistentPageCache> {
        &self.cache
    }

    #[must_use]
    pub fn inner(&self) -> &R {
        &self.inner
    }

    #[must_use]
    pub fn object_checksum(&self) -> [u8; 32] {
        self.object_checksum
    }

    fn resource_error(&self, message: impl Into<String>) -> ResourceError {
        ResourceError::Io {
            locator: self.inner.locator().redacted(),
            message: message.into(),
        }
    }

    fn cache_error(&self, error: PersistentCacheError) -> ResourceError {
        self.resource_error(format!("persistent cache: {error}"))
    }

    fn identity_error(&self) -> ResourceError {
        ResourceError::CacheIdentityChanged(self.inner.locator().redacted())
    }

    fn checked_fetched_bytes(
        &self,
        range: ByteRange,
        expected_checksum: [u8; 32],
        bytes: Vec<u8>,
    ) -> ResourceResult<Vec<u8>> {
        let expected_len = usize::try_from(range.length)
            .map_err(|_| self.resource_error("fetched range length does not fit this platform"))?;
        if bytes.len() != expected_len || sha256(&bytes) != expected_checksum {
            self.metrics.stale_cache_rejection();
            return Err(self.identity_error());
        }
        Ok(bytes)
    }
}

impl<R, P> RangeReader for CachedRangeReader<R, P>
where
    R: RangeReader,
    P: PageKeyResolver,
{
    fn locator(&self) -> &crate::resource::ResourceLocator {
        self.inner.locator()
    }

    fn metadata(&self) -> ResourceResult<ResourceMetadata> {
        self.inner.metadata()
    }

    fn read_range(&self, range: ByteRange) -> ResourceResult<Vec<u8>> {
        if range.is_empty() {
            return self.inner.read_range(range);
        }
        let Some(descriptor) = self.resolver.resolve(range) else {
            return self.inner.read_range(range);
        };
        let key = PageKey::new(
            self.object_checksum,
            descriptor.section_id,
            descriptor.page_id,
            range,
            descriptor.page_checksum,
        )
        .map_err(|error| self.cache_error(error))?;

        match self
            .cache
            .get(&key)
            .map_err(|error| self.cache_error(error))?
        {
            CacheLookup::Hit(bytes) => {
                let bytes = self.checked_fetched_bytes(range, descriptor.page_checksum, bytes)?;
                self.metrics.range_request();
                self.metrics.requested_bytes(range.length);
                self.metrics.returned_bytes(range.length);
                self.metrics.decoded_bytes(range.length);
                self.metrics.cache_hit();
                self.metrics.cache_bytes(range.length);
                Ok(bytes)
            }
            CacheLookup::Miss(reason) => {
                self.metrics.cache_miss();
                if matches!(
                    reason,
                    CacheMissReason::StaleIdentity | CacheMissReason::Corrupt
                ) {
                    self.metrics.stale_cache_rejection();
                }
                let bytes = self.inner.read_range(range)?;
                let bytes = self.checked_fetched_bytes(range, descriptor.page_checksum, bytes)?;
                self.cache
                    .put(&key, &bytes)
                    .map_err(|error| self.cache_error(error))?;
                Ok(bytes)
            }
        }
    }

    fn read_range_bytes(&self, range: ByteRange) -> ResourceResult<RangeBytes<'_>> {
        self.read_range(range).map(RangeBytes::Owned)
    }

    fn stream(&self) -> ResourceResult<Box<dyn Read + Send>> {
        self.inner.stream()
    }

    fn metrics(&self) -> ResourceMetrics {
        self.inner.metrics().saturating_add(self.metrics.snapshot())
    }
}

fn sha256(bytes: &[u8]) -> [u8; 32] {
    Sha256::digest(bytes).into()
}

/// A resolver for one exact range, useful for callers whose index already
/// identifies one page per read.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ExactPageResolver {
    pub descriptor: PageDescriptor,
}

impl PageKeyResolver for ExactPageResolver {
    fn resolve(&self, _range: ByteRange) -> Option<PageDescriptor> {
        Some(self.descriptor)
    }
}

/// Convenience helper for a resolver that always returns one descriptor.
#[must_use]
pub const fn exact_page(descriptor: PageDescriptor) -> ExactPageResolver {
    ExactPageResolver { descriptor }
}
