//! HTTP, HTTPS, and S3 range readers.
//!
//! The implementation deliberately keeps the resource layer synchronous.  It
//! uses the ubiquitous `curl` command for transport so the core crate does not
//! need to own a TLS stack or cloud credential implementation.  Child stdout
//! is exposed directly for streaming; range and metadata requests validate all
//! returned lengths before handing bytes to callers.

use super::cache::BlockCache;
use super::error::{redact_message, retry_count, transport_error, validate_range};
use super::metrics::MetricsCounter;
use super::{
    ByteRange, CacheIdentity, RangeReader, ResourceError, ResourceLocator, ResourceMetadata,
    ResourceOpenOptions, ResourceResult, ResourceScheme,
};
use std::io::{self, Read};
use std::process::{Child, ChildStdout, Command, Stdio};
use std::sync::{Arc, Mutex};

const CURL: &str = "curl";

/// A remote object addressed by HTTP(S) or S3.
pub struct ObjectResource {
    locator: ResourceLocator,
    request_url: String,
    options: ResourceOpenOptions,
    cache: Arc<BlockCache>,
    metrics: Arc<MetricsCounter>,
    metadata: Mutex<Option<ResourceMetadata>>,
}

impl std::fmt::Debug for ObjectResource {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ObjectResource")
            .field("locator", &self.locator.redacted())
            .field("options", &self.options)
            .finish_non_exhaustive()
    }
}

impl ObjectResource {
    /// Open an HTTP(S) or S3 resource lazily.
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
            ResourceScheme::Http | ResourceScheme::Https | ResourceScheme::S3
        ) {
            return Err(ResourceError::UnsupportedScheme(format!(
                "{:?}",
                locator.scheme()
            )));
        }
        let cache = BlockCache::new(options.cache_block_bytes, options.max_cache_bytes)?;
        let request_url = request_url(&locator)?;
        Ok(Self {
            locator,
            request_url,
            options,
            cache: Arc::new(cache),
            metrics: Arc::new(MetricsCounter::default()),
            metadata: Mutex::new(None),
        })
    }

    #[must_use]
    pub fn cache(&self) -> &BlockCache {
        &self.cache
    }

    fn curl_args(&self) -> Vec<String> {
        vec![
            "--fail".to_string(),
            "--silent".to_string(),
            "--show-error".to_string(),
            "--location".to_string(),
            "--connect-timeout".to_string(),
            self.options.request_timeout_seconds.max(1).to_string(),
            "--max-time".to_string(),
            self.options.request_timeout_seconds.max(1).to_string(),
        ]
    }

    fn execute(&self, extra: &[String]) -> ResourceResult<Vec<u8>> {
        let attempts = retry_count(self.options.max_retries);
        let mut last_error = None;
        for attempt in 0..attempts {
            let mut command = Command::new(CURL);
            command.args(self.curl_args());
            command.args(extra);
            command.arg(&self.request_url);
            let output = command.output().map_err(|error| {
                transport_error(&self.locator, format!("unable to execute curl: {error}"))
            })?;
            if output.status.success() {
                return Ok(output.stdout);
            }
            if attempt + 1 < attempts {
                self.metrics.retry();
            }
            last_error = Some(format!("curl exited with status {}", output.status));
        }
        Err(transport_error(
            &self.locator,
            redact_message(&last_error.unwrap_or_else(|| "curl request failed".to_string())),
        ))
    }

    fn execute_headers(&self, extra: &[String]) -> ResourceResult<Vec<u8>> {
        let mut args = vec!["--dump-header".to_string(), "-".to_string()];
        args.extend_from_slice(extra);
        // Discard the response body.  `--head` is not used here because a few
        // object servers reject HEAD while accepting a range GET.
        args.extend(["--output".to_string(), "/dev/null".to_string()]);
        self.execute(&args)
    }

    fn fetch_metadata(&self) -> ResourceResult<ResourceMetadata> {
        let mut headers = self.execute_headers(&["--head".to_string()])?;
        let mut metadata = parse_metadata_headers(&headers);
        if metadata.size.is_none() {
            headers = self.execute_headers(&["--range".to_string(), "0-0".to_string()])?;
            metadata = parse_metadata_headers(&headers);
        }
        let size = metadata.size.ok_or_else(|| {
            transport_error(
                &self.locator,
                "remote response did not include Content-Length or Content-Range",
            )
        })?;
        Ok(ResourceMetadata {
            size,
            etag: metadata.etag,
            last_modified: metadata.last_modified,
            accepts_ranges: metadata.accepts_ranges,
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

    fn read_uncached(
        &self,
        range: ByteRange,
        metadata: &ResourceMetadata,
    ) -> ResourceResult<Vec<u8>> {
        let expected = usize::try_from(range.length).map_err(|_| ResourceError::Io {
            locator: self.locator.redacted(),
            message: "requested range is too large for this platform".to_string(),
        })?;
        if !metadata.accepts_ranges && !self.options.allow_full_download_fallback {
            return Err(ResourceError::RangeUnsupported(self.locator.redacted()));
        }
        let mut extra = Vec::new();
        if metadata.accepts_ranges {
            let end = range
                .end()?
                .checked_sub(1)
                .ok_or_else(|| ResourceError::Io {
                    locator: self.locator.redacted(),
                    message: "non-empty range has no end".to_string(),
                })?;
            extra.extend(["--range".to_string(), format!("{}-{end}", range.offset)]);
        }
        let bytes = self.execute(&extra)?;
        self.metrics
            .remote_bytes(u64::try_from(bytes.len()).unwrap_or(u64::MAX));
        if bytes.len() == expected {
            return Ok(bytes);
        }
        // Some servers advertise range support but return a complete object.
        // Slice it only when explicitly permitted by the opening options.
        if self.options.allow_full_download_fallback
            && bytes.len() == usize::try_from(metadata.size).unwrap_or(usize::MAX)
        {
            let start = usize::try_from(range.offset).map_err(|_| ResourceError::Io {
                locator: self.locator.redacted(),
                message: "range offset is too large for this platform".to_string(),
            })?;
            let end = start
                .checked_add(expected)
                .ok_or_else(|| ResourceError::Io {
                    locator: self.locator.redacted(),
                    message: "range end is too large for this platform".to_string(),
                })?;
            if end <= bytes.len() {
                return Ok(bytes[start..end].to_vec());
            }
        }
        Err(transport_error(
            &self.locator,
            format!(
                "range response length {} does not match requested {}",
                bytes.len(),
                expected
            ),
        ))
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
        let requested_end = range.end()?;
        let mut output =
            Vec::with_capacity(
                usize::try_from(range.length).map_err(|_| ResourceError::Io {
                    locator: self.locator.redacted(),
                    message: "requested range is too large for this platform".to_string(),
                })?,
            );
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
                let bytes = self.read_uncached(block_range, metadata)?;
                self.cache.insert(identity, block_offset, bytes.clone())?;
                bytes
            };
            let copy_start = usize::try_from(range.offset.max(block_offset) - block_offset)
                .map_err(|_| ResourceError::Io {
                    locator: self.locator.redacted(),
                    message: "cache offset does not fit in platform usize".to_string(),
                })?;
            let copy_end =
                usize::try_from(requested_end.min(block_end) - block_offset).map_err(|_| {
                    ResourceError::Io {
                        locator: self.locator.redacted(),
                        message: "cache offset does not fit in platform usize".to_string(),
                    }
                })?;
            if copy_start > copy_end || copy_end > block.len() {
                return Err(transport_error(
                    &self.locator,
                    "cached block has an invalid length",
                ));
            }
            output.extend_from_slice(&block[copy_start..copy_end]);
        }
        Ok(output)
    }
}

impl RangeReader for ObjectResource {
    fn locator(&self) -> &ResourceLocator {
        &self.locator
    }

    fn metadata(&self) -> ResourceResult<ResourceMetadata> {
        self.metrics.metadata_request();
        let mut cached = self.metadata.lock().map_err(|_| ResourceError::Io {
            locator: self.locator.redacted(),
            message: "metadata lock poisoned".to_string(),
        })?;
        if let Some(metadata) = cached.as_ref() {
            return Ok(metadata.clone());
        }
        let metadata = self.fetch_metadata()?;
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
        let mut command = Command::new(CURL);
        command
            .args(self.curl_args())
            .arg(&self.request_url)
            .stdout(Stdio::piped())
            .stderr(Stdio::null());
        let mut child = command.spawn().map_err(|error| {
            transport_error(&self.locator, format!("unable to execute curl: {error}"))
        })?;
        let stdout = child.stdout.take().ok_or_else(|| {
            transport_error(&self.locator, "curl did not provide a streaming stdout")
        })?;
        Ok(Box::new(ChildReader {
            child,
            stdout,
            metrics: Arc::clone(&self.metrics),
        }))
    }

    fn metrics(&self) -> super::ResourceMetrics {
        self.metrics.snapshot()
    }
}

struct ChildReader {
    child: Child,
    stdout: ChildStdout,
    metrics: Arc<MetricsCounter>,
}

impl Read for ChildReader {
    fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
        let count = self.stdout.read(bytes)?;
        self.metrics
            .remote_bytes(u64::try_from(count).unwrap_or(u64::MAX));
        if count == 0 {
            let status = self.child.wait()?;
            if !status.success() {
                return Err(io::Error::other("remote stream command failed"));
            }
        }
        Ok(count)
    }
}

impl Drop for ChildReader {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

#[derive(Default)]
struct ParsedHeaders {
    size: Option<u64>,
    etag: Option<String>,
    last_modified: Option<String>,
    accepts_ranges: bool,
}

fn parse_metadata_headers(bytes: &[u8]) -> ParsedHeaders {
    let text = String::from_utf8_lossy(bytes);
    let block = text
        .split("\r\n\r\n")
        .filter(|part| part.lines().any(|line| line.starts_with("HTTP/")))
        .last()
        .unwrap_or(text.as_ref());
    let mut parsed = ParsedHeaders::default();
    let mut status = 0u16;
    let mut content_range_size = None;
    for line in block.lines() {
        if let Some(value) = line.strip_prefix("HTTP/") {
            status = value
                .split_whitespace()
                .nth(1)
                .and_then(|value| value.parse().ok())
                .unwrap_or(0);
            continue;
        }
        let Some((name, value)) = line.split_once(':') else {
            continue;
        };
        let value = value.trim();
        match name.to_ascii_lowercase().as_str() {
            "content-length" => parsed.size = value.parse().ok(),
            "content-range" => {
                if let Some((_, total)) = value.rsplit_once('/') {
                    content_range_size = total.parse().ok();
                }
            }
            "accept-ranges" => parsed.accepts_ranges = value.eq_ignore_ascii_case("bytes"),
            "etag" => parsed.etag = Some(value.trim_matches('"').to_string()),
            "last-modified" => parsed.last_modified = Some(value.to_string()),
            _ => {}
        }
    }
    if let Some(size) = content_range_size {
        parsed.size = Some(size);
    }
    parsed.accepts_ranges |= status == 206;
    parsed
}

fn request_url(locator: &ResourceLocator) -> ResourceResult<String> {
    match locator.scheme() {
        ResourceScheme::Http | ResourceScheme::Https => Ok(locator.raw().to_string()),
        ResourceScheme::S3 => {
            let remainder = locator.raw().strip_prefix("s3://").unwrap_or_default();
            let (bucket, key) = remainder.split_once('/').ok_or_else(|| {
                ResourceError::InvalidLocator("S3 URL must contain a bucket and key".to_string())
            })?;
            let endpoint = std::env::var("JAM_S3_ENDPOINT")
                .or_else(|_| std::env::var("AWS_ENDPOINT_URL"))
                .or_else(|_| std::env::var("S3_ENDPOINT_URL"))
                .unwrap_or_else(|_| "https://s3.amazonaws.com".to_string())
                .trim_end_matches('/')
                .to_string();
            if endpoint.contains("{bucket}") {
                Ok(endpoint.replace("{bucket}", bucket) + "/" + key)
            } else {
                Ok(format!("{endpoint}/{bucket}/{key}"))
            }
        }
        _ => Err(ResourceError::UnsupportedScheme(format!(
            "{:?}",
            locator.scheme()
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::{parse_metadata_headers, request_url};
    use crate::resource::ResourceLocator;

    #[test]
    fn parses_final_redirect_headers() {
        let headers = b"HTTP/1.1 302 Found\r\nLocation: x\r\n\r\nHTTP/1.1 206 Partial Content\r\nContent-Range: bytes 0-0/42\r\nAccept-Ranges: bytes\r\n\r\n";
        let metadata = parse_metadata_headers(headers);
        assert_eq!(metadata.size, Some(42));
        assert!(metadata.accepts_ranges);
    }

    #[test]
    fn maps_s3_locator_to_virtual_path() {
        let locator = ResourceLocator::parse("s3://bucket/path/object.jma").unwrap();
        assert_eq!(
            request_url(&locator).unwrap(),
            "https://s3.amazonaws.com/bucket/path/object.jma"
        );
    }
}
