//! HTTP, HTTPS, and S3 range readers.
//!
//! The implementation deliberately keeps the resource layer synchronous.  It
//! uses the ubiquitous `curl` command for transport so the core crate does not
//! need to own a TLS stack or cloud credential implementation.  Child stdout
//! is exposed directly for streaming; range and metadata requests validate all
//! returned lengths before handing bytes to callers.

use super::cache::BlockCache;
use super::error::{
    http_status_error, is_retryable, retry_count, timeout_error, transport_error, validate_range,
};
use super::metrics::MetricsCounter;
use super::{
    ByteRange, CacheIdentity, RangeReader, ResourceError, ResourceLocator, ResourceMetadata,
    ResourceOpenOptions, ResourceResult, ResourceScheme,
};
use std::io::{self, Read};
use std::process::{Child, ChildStdout, Command, Stdio};
use std::sync::{Arc, Mutex};

const CURL: &str = "curl";
const MAX_RESPONSE_HEADER_BYTES: u64 = 1024 * 1024;

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

    fn execute(
        &self,
        extra: &[String],
        max_output_bytes: u64,
        requested_bytes: Option<u64>,
    ) -> ResourceResult<Vec<u8>> {
        let attempts = retry_count(self.options.max_retries);
        let mut last_error = None;
        for attempt in 0..attempts {
            let mut command = Command::new(CURL);
            command.args(self.curl_args());
            command.args(extra);
            command.arg(&self.request_url);
            let is_head = extra.iter().any(|argument| argument == "--head");
            if is_head {
                self.metrics.head_request();
            } else {
                self.metrics.get_request();
            }
            if let Some(bytes) = requested_bytes {
                self.metrics.requested_bytes(bytes);
            }
            let mut child = match command.stdout(Stdio::piped()).stderr(Stdio::null()).spawn() {
                Ok(child) => child,
                Err(error) => {
                    let resource_error =
                        transport_error(&self.locator, format!("unable to execute curl: {error}"));
                    if is_retryable(&resource_error) && attempt + 1 < attempts {
                        self.metrics.retry();
                        last_error = Some(resource_error);
                        continue;
                    }
                    return Err(resource_error);
                }
            };
            let mut stdout = match child.stdout.take() {
                Some(stdout) => stdout,
                None => {
                    let resource_error =
                        transport_error(&self.locator, "curl did not provide a readable stdout");
                    let _ = child.kill();
                    let _ = child.wait();
                    if is_retryable(&resource_error) && attempt + 1 < attempts {
                        self.metrics.retry();
                        last_error = Some(resource_error);
                        continue;
                    }
                    return Err(resource_error);
                }
            };
            let max_output_bytes = usize::try_from(max_output_bytes).unwrap_or(usize::MAX);
            let mut output = Vec::new();
            let mut buffer = [0_u8; 8192];
            let mut read_error = None;
            loop {
                let count = match stdout.read(&mut buffer) {
                    Ok(count) => count,
                    Err(error) => {
                        let _ = child.kill();
                        let _ = child.wait();
                        read_error = Some(if error.kind() == io::ErrorKind::TimedOut {
                            timeout_error(&self.locator)
                        } else {
                            transport_error(
                                &self.locator,
                                format!("unable to read curl response: {error}"),
                            )
                        });
                        break;
                    }
                };
                if count == 0 {
                    break;
                }
                if output.len().saturating_add(count) > max_output_bytes {
                    let _ = child.kill();
                    let _ = child.wait();
                    read_error = Some(transport_error(
                        &self.locator,
                        format!("response exceeds configured limit of {max_output_bytes} bytes"),
                    ));
                    break;
                }
                output.extend_from_slice(&buffer[..count]);
            }
            if let Some(resource_error) = read_error {
                if is_retryable(&resource_error) && attempt + 1 < attempts {
                    self.metrics.retry();
                    last_error = Some(resource_error);
                    continue;
                }
                return Err(resource_error);
            }
            let status = match child.wait() {
                Ok(status) => status,
                Err(error) => {
                    let resource_error =
                        transport_error(&self.locator, format!("unable to wait for curl: {error}"));
                    if is_retryable(&resource_error) && attempt + 1 < attempts {
                        self.metrics.retry();
                        last_error = Some(resource_error);
                        continue;
                    }
                    return Err(resource_error);
                }
            };
            if let Some(http_status) = response_status(&output).filter(|status| *status >= 400) {
                let resource_error = http_status_error(&self.locator, http_status);
                if is_retryable(&resource_error) && attempt + 1 < attempts {
                    self.metrics.retry();
                    last_error = Some(resource_error);
                    continue;
                }
                return Err(resource_error);
            }
            if status.success() {
                return Ok(output);
            }
            let resource_error = if status.code() == Some(28) {
                timeout_error(&self.locator)
            } else {
                transport_error(&self.locator, format!("curl exited with status {status}"))
            };
            if is_retryable(&resource_error) && attempt + 1 < attempts {
                self.metrics.retry();
                last_error = Some(resource_error);
                continue;
            }
            return Err(resource_error);
        }
        Err(last_error.unwrap_or_else(|| transport_error(&self.locator, "curl request failed")))
    }

    fn execute_headers(
        &self,
        extra: &[String],
        requested_bytes: Option<u64>,
    ) -> ResourceResult<Vec<u8>> {
        let mut args = vec!["--dump-header".to_string(), "-".to_string()];
        args.extend_from_slice(extra);
        // Discard the response body.  `--head` is not used here because a few
        // object servers reject HEAD while accepting a range GET.
        args.extend(["--output".to_string(), "/dev/null".to_string()]);
        self.execute(&args, MAX_RESPONSE_HEADER_BYTES, requested_bytes)
    }

    fn fetch_metadata(&self) -> ResourceResult<ResourceMetadata> {
        let (mut metadata, mut head_error) =
            match self.execute_headers(&["--head".to_string()], None) {
                Ok(headers) => (parse_metadata_headers(&headers), None),
                Err(error @ ResourceError::HttpStatus { status, .. })
                    if status == 403 || status == 404 || (500..600).contains(&status) =>
                {
                    return Err(error);
                }
                Err(error @ ResourceError::Timeout { .. }) => return Err(error),
                Err(error) => (ParsedHeaders::default(), Some(error)),
            };
        if metadata.size.is_none() {
            let headers =
                self.execute_headers(&["--range".to_string(), "0-0".to_string()], Some(1));
            match headers {
                Ok(headers) => metadata = parse_metadata_headers(&headers),
                Err(error) => {
                    if head_error.is_none() {
                        head_error = Some(error);
                    }
                }
            }
        }
        let size = metadata.size.ok_or_else(|| {
            head_error.unwrap_or_else(|| {
                transport_error(
                    &self.locator,
                    "remote response did not include Content-Length or Content-Range",
                )
            })
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
        let mut extra = vec!["--dump-header".to_string(), "-".to_string()];
        let requested_bytes = if metadata.accepts_ranges {
            u64::try_from(expected).unwrap_or(u64::MAX)
        } else {
            metadata.size
        };
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
        let fallback_limit = metadata.size.min(self.options.max_cache_bytes);
        let max_body = if metadata.accepts_ranges {
            u64::try_from(expected)
                .unwrap_or(u64::MAX)
                .max(fallback_limit)
        } else {
            if metadata.size > self.options.max_cache_bytes {
                return Err(transport_error(
                    &self.locator,
                    format!(
                        "full-object fallback of {} bytes exceeds max_cache_bytes {}",
                        metadata.size, self.options.max_cache_bytes
                    ),
                ));
            }
            metadata.size
        };
        let max_output = MAX_RESPONSE_HEADER_BYTES.saturating_add(max_body);
        let raw = self.execute(&extra, max_output, Some(requested_bytes))?;
        let response = parse_http_response(&raw, &self.locator)?;
        let body_length = u64::try_from(response.body.len()).unwrap_or(u64::MAX);
        self.metrics.returned_bytes(body_length);
        self.metrics.remote_bytes(body_length);
        match response.status {
            206 => {
                let expected_end =
                    range
                        .end()?
                        .checked_sub(1)
                        .ok_or_else(|| ResourceError::Io {
                            locator: self.locator.redacted(),
                            message: "non-empty range has no end".to_string(),
                        })?;
                let Some(content_range) = response.headers.content_range else {
                    return Err(transport_error(
                        &self.locator,
                        "206 response did not include Content-Range",
                    ));
                };
                if content_range.start != range.offset
                    || content_range.end != expected_end
                    || content_range.total != metadata.size
                {
                    return Err(transport_error(
                        &self.locator,
                        format!(
                            "206 Content-Range bytes {}-{}/{} does not match requested {}-{}/{}",
                            content_range.start,
                            content_range.end,
                            content_range.total,
                            range.offset,
                            expected_end,
                            metadata.size
                        ),
                    ));
                }
                if response.body.len() != expected {
                    return Err(transport_error(
                        &self.locator,
                        format!(
                            "206 response length {} does not match requested {}",
                            response.body.len(),
                            expected
                        ),
                    ));
                }
                Ok(response.body)
            }
            200 => {
                if !self.options.allow_full_download_fallback {
                    return Err(ResourceError::RangeUnsupported(self.locator.redacted()));
                }
                if metadata.size > self.options.max_cache_bytes {
                    return Err(transport_error(
                        &self.locator,
                        format!(
                            "full-object fallback of {} bytes exceeds max_cache_bytes {}",
                            metadata.size, self.options.max_cache_bytes
                        ),
                    ));
                }
                let object_size =
                    usize::try_from(metadata.size).map_err(|_| ResourceError::Io {
                        locator: self.locator.redacted(),
                        message: "full object is too large for this platform".to_string(),
                    })?;
                if response.body.len() != object_size {
                    return Err(transport_error(
                        &self.locator,
                        format!(
                            "full-object fallback length {} does not match object size {}",
                            response.body.len(),
                            object_size
                        ),
                    ));
                }
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
                if end > response.body.len() {
                    return Err(transport_error(
                        &self.locator,
                        "full-object fallback does not contain requested range",
                    ));
                }
                self.metrics.full_object_fallback();
                Ok(response.body[start..end].to_vec())
            }
            status => Err(transport_error(
                &self.locator,
                format!("unexpected HTTP status {status} for range request"),
            )),
        }
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
                let bytes = self.read_uncached(block_range, metadata)?;
                self.cache.insert_with_metrics(
                    identity,
                    block_offset,
                    bytes.clone(),
                    &self.metrics,
                )?;
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
        self.metrics
            .decoded_bytes(u64::try_from(output.len()).unwrap_or(u64::MAX));
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
        // Probe metadata before handing out a streaming handle.  This turns
        // HTTP authorization/not-found failures into structured resource
        // errors instead of hiding them behind a later parser I/O error.
        self.metadata()?;
        self.metrics.get_request();
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
        let count = u64::try_from(count).unwrap_or(u64::MAX);
        self.metrics.remote_bytes(count);
        self.metrics.returned_bytes(count);
        self.metrics.decoded_bytes(count);
        if count == 0 {
            let status = self.child.wait()?;
            if !status.success() {
                return Err(io::Error::other("remote stream command failed"));
            }
        }
        Ok(usize::try_from(count).unwrap_or(usize::MAX))
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
    status: u16,
    size: Option<u64>,
    etag: Option<String>,
    last_modified: Option<String>,
    accepts_ranges: bool,
    content_range: Option<ContentRange>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ContentRange {
    start: u64,
    end: u64,
    total: u64,
}

struct HttpResponse {
    status: u16,
    headers: ParsedHeaders,
    body: Vec<u8>,
}

fn response_status(bytes: &[u8]) -> Option<u16> {
    let mut cursor = 0;
    loop {
        let relative_start = find_http_status(bytes, cursor)?;
        let start = cursor + relative_start;
        let relative_end = find_bytes(&bytes[start..], b"\r\n\r\n")?;
        let header_end = start + relative_end + 4;
        let headers = parse_header_block(&bytes[start..header_end]);
        if headers.status == 0 {
            return None;
        }
        if headers.status == 100 || (300..400).contains(&headers.status) {
            cursor = header_end;
            continue;
        }
        return Some(headers.status);
    }
}

fn parse_http_response(bytes: &[u8], locator: &ResourceLocator) -> ResourceResult<HttpResponse> {
    let mut cursor = 0;
    loop {
        let Some(relative_start) = find_http_status(bytes, cursor) else {
            return Err(transport_error(
                locator,
                "curl response did not include an HTTP status header",
            ));
        };
        let start = cursor + relative_start;
        let Some(relative_end) = find_bytes(&bytes[start..], b"\r\n\r\n") else {
            return Err(transport_error(
                locator,
                "curl response headers were incomplete",
            ));
        };
        let header_end = start + relative_end + 4;
        let headers = parse_header_block(&bytes[start..header_end]);
        if headers.status == 0 {
            return Err(transport_error(
                locator,
                "curl response had an invalid HTTP status header",
            ));
        }
        // Curl emits intermediate redirect and informational response headers
        // before the final response when --location is enabled.
        if headers.status == 100 || (300..400).contains(&headers.status) {
            cursor = header_end;
            continue;
        }
        return Ok(HttpResponse {
            status: headers.status,
            headers,
            body: bytes[header_end..].to_vec(),
        });
    }
}

fn parse_metadata_headers(bytes: &[u8]) -> ParsedHeaders {
    let mut cursor = 0;
    let mut parsed = ParsedHeaders::default();
    while let Some(relative_start) = find_http_status(bytes, cursor) {
        let start = cursor + relative_start;
        let Some(relative_end) = find_bytes(&bytes[start..], b"\r\n\r\n") else {
            break;
        };
        let header_end = start + relative_end + 4;
        parsed = parse_header_block(&bytes[start..header_end]);
        cursor = header_end;
    }
    parsed
}

fn parse_header_block(block: &[u8]) -> ParsedHeaders {
    let text = String::from_utf8_lossy(block);
    let mut parsed = ParsedHeaders::default();
    for line in text.lines() {
        if let Some(value) = line.strip_prefix("HTTP/") {
            parsed.status = value
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
                if let Some(content_range) = parse_content_range(value) {
                    parsed.size = Some(content_range.total);
                    parsed.content_range = Some(content_range);
                }
            }
            "accept-ranges" => parsed.accepts_ranges = value.eq_ignore_ascii_case("bytes"),
            "etag" => parsed.etag = Some(value.trim_matches('"').to_string()),
            "last-modified" => parsed.last_modified = Some(value.to_string()),
            _ => {}
        }
    }
    parsed.accepts_ranges |= parsed.status == 206;
    parsed
}

fn parse_content_range(value: &str) -> Option<ContentRange> {
    let (unit, value) = value.split_once(|character: char| character.is_ascii_whitespace())?;
    let value = value.trim_start();
    if !unit.eq_ignore_ascii_case("bytes") {
        return None;
    }
    let value = value.trim_start();
    let (range, total) = value.split_once('/')?;
    let total = total.parse().ok()?;
    let (start, end) = range.split_once('-')?;
    Some(ContentRange {
        start: start.parse().ok()?,
        end: end.parse().ok()?,
        total,
    })
}

fn find_http_status(bytes: &[u8], from: usize) -> Option<usize> {
    bytes
        .get(from..)?
        .windows(5)
        .position(|window| window == b"HTTP/")
}

fn find_bytes(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    haystack
        .windows(needle.len())
        .position(|window| window == needle)
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
    use super::{parse_http_response, parse_metadata_headers, request_url};
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

    #[test]
    fn parses_response_body_separately_from_content_range() {
        let locator = ResourceLocator::parse("https://example.org/object").unwrap();
        let response = parse_http_response(
            b"HTTP/1.1 206 Partial Content\r\nContent-Length: 4\r\nContent-Range: bytes 3-6/16\r\n\r\n3456",
            &locator,
        )
        .unwrap();
        assert_eq!(response.status, 206);
        assert_eq!(response.body, b"3456");
        assert_eq!(
            response.headers.content_range,
            Some(super::ContentRange {
                start: 3,
                end: 6,
                total: 16,
            })
        );
    }

    #[test]
    fn parses_redirect_then_final_response() {
        let locator = ResourceLocator::parse("https://example.org/object").unwrap();
        let response = parse_http_response(
            b"HTTP/1.1 302 Found\r\nLocation: /object\r\n\r\nHTTP/1.1 200 OK\r\nContent-Length: 4\r\n\r\nbody",
            &locator,
        )
        .unwrap();
        assert_eq!(response.status, 200);
        assert_eq!(response.body, b"body");
    }
}
