//! Upload a completed result file to an HTTP(S) or S3-compatible resource.
//!
//! Uploads are deliberately one-shot from the caller's finalized local file.
//! A retry repeats that complete request; it never exposes a partially written
//! remote result as a successful upload.

use super::error::{http_status_error, is_retryable, retry_count, timeout_error, transport_error};
use super::{ResourceError, ResourceLocator, ResourceResult, ResourceScheme};
use std::fs;
use std::io::Read;
use std::path::Path;
use std::process::{Command, Stdio};

const CURL: &str = "curl";
const MAX_RESPONSE_HEADER_BYTES: usize = 1024 * 1024;

/// Options controlling a single-file result upload.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct UploadOptions {
    pub request_timeout_seconds: u64,
    pub max_retries: u32,
}

impl Default for UploadOptions {
    fn default() -> Self {
        Self {
            request_timeout_seconds: 60,
            max_retries: 3,
        }
    }
}

/// Summary returned after a complete upload request succeeds.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct UploadResult {
    pub bytes_uploaded: u64,
    pub attempts: u32,
    pub status: u16,
}

/// Upload a finalized local file to an HTTP(S) or S3-compatible object URL.
///
/// The local file is never modified.  A successful return means curl
/// received a successful HTTP response for the complete file.  HTTP 403 and
/// 404 responses are returned without retrying; retryable transport failures
/// and 5xx responses follow `options.max_retries`.
pub fn upload_file(
    locator: impl AsRef<str>,
    local_file: impl AsRef<Path>,
    options: UploadOptions,
) -> ResourceResult<UploadResult> {
    let locator = ResourceLocator::parse(locator.as_ref())?;
    if !matches!(
        locator.scheme(),
        ResourceScheme::Http | ResourceScheme::Https | ResourceScheme::S3
    ) {
        return Err(ResourceError::UnsupportedScheme(format!(
            "{:?}",
            locator.scheme()
        )));
    }
    let local_file = local_file.as_ref();
    let metadata = fs::metadata(local_file).map_err(|error| ResourceError::Io {
        locator: locator.redacted(),
        message: format!("unable to inspect upload file: {error}"),
    })?;
    if !metadata.is_file() {
        return Err(ResourceError::Io {
            locator: locator.redacted(),
            message: "upload source is not a regular file".to_string(),
        });
    }
    let bytes_uploaded = metadata.len();
    let request_url = request_url(&locator)?;
    let attempts = retry_count(options.max_retries);
    let mut last_error = None;

    for attempt in 0..attempts {
        let mut command = Command::new(CURL);
        let timeout = options.request_timeout_seconds.max(1).to_string();
        command.args([
            "--fail".to_string(),
            "--silent".to_string(),
            "--show-error".to_string(),
            "--location".to_string(),
            "--connect-timeout".to_string(),
            timeout.clone(),
            "--max-time".to_string(),
            timeout,
            "--dump-header".to_string(),
            "-".to_string(),
            "--output".to_string(),
            "/dev/null".to_string(),
            "--request".to_string(),
            "PUT".to_string(),
            "--upload-file".to_string(),
        ]);
        command.arg(local_file).arg(&request_url);
        let mut child = match command.stdout(Stdio::piped()).stderr(Stdio::null()).spawn() {
            Ok(child) => child,
            Err(error) => {
                let resource_error =
                    transport_error(&locator, format!("unable to execute curl: {error}"));
                if is_retryable(&resource_error) && attempt + 1 < attempts {
                    last_error = Some(resource_error);
                    continue;
                }
                return Err(resource_error);
            }
        };
        let mut headers = Vec::new();
        let max_output = MAX_RESPONSE_HEADER_BYTES;
        let read_error = match child.stdout.take() {
            Some(mut stdout) => {
                let mut buffer = [0_u8; 8192];
                loop {
                    let count = match stdout.read(&mut buffer) {
                        Ok(count) => count,
                        Err(error) => {
                            let _ = child.kill();
                            let _ = child.wait();
                            break Some(if error.kind() == std::io::ErrorKind::TimedOut {
                                timeout_error(&locator)
                            } else {
                                transport_error(
                                    &locator,
                                    format!("unable to read curl response: {error}"),
                                )
                            });
                        }
                    };
                    if count == 0 {
                        break None;
                    }
                    if headers.len().saturating_add(count) > max_output {
                        let _ = child.kill();
                        let _ = child.wait();
                        break Some(transport_error(
                            &locator,
                            format!(
                                "response headers exceed configured limit of {max_output} bytes"
                            ),
                        ));
                    }
                    headers.extend_from_slice(&buffer[..count]);
                }
            }
            None => {
                let _ = child.kill();
                let _ = child.wait();
                Some(transport_error(
                    &locator,
                    "curl did not provide readable response headers",
                ))
            }
        };
        if let Some(resource_error) = read_error {
            if is_retryable(&resource_error) && attempt + 1 < attempts {
                last_error = Some(resource_error);
                continue;
            }
            return Err(resource_error);
        }
        let status = match child.wait() {
            Ok(status) => status,
            Err(error) => {
                let resource_error =
                    transport_error(&locator, format!("unable to wait for curl: {error}"));
                if is_retryable(&resource_error) && attempt + 1 < attempts {
                    last_error = Some(resource_error);
                    continue;
                }
                return Err(resource_error);
            }
        };
        if let Some(http_status) = response_status(&headers).filter(|status| *status >= 400) {
            let resource_error = http_status_error(&locator, http_status);
            if is_retryable(&resource_error) && attempt + 1 < attempts {
                last_error = Some(resource_error);
                continue;
            }
            return Err(resource_error);
        }
        if status.success() {
            let response_status = response_status(&headers).unwrap_or(200);
            return Ok(UploadResult {
                bytes_uploaded,
                attempts: attempt.saturating_add(1),
                status: response_status,
            });
        }
        let resource_error = if status.code() == Some(28) {
            timeout_error(&locator)
        } else {
            transport_error(&locator, format!("curl upload exited with status {status}"))
        };
        if is_retryable(&resource_error) && attempt + 1 < attempts {
            last_error = Some(resource_error);
            continue;
        }
        return Err(resource_error);
    }

    Err(last_error.unwrap_or_else(|| transport_error(&locator, "curl upload failed")))
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

fn response_status(bytes: &[u8]) -> Option<u16> {
    let mut cursor = 0;
    loop {
        let relative_start = find_http_status(bytes, cursor)?;
        let start = cursor + relative_start;
        let relative_end = find_bytes(&bytes[start..], b"\r\n\r\n")?;
        let header_end = start + relative_end + 4;
        let status = parse_status(&bytes[start..header_end])?;
        if status == 100 || (300..400).contains(&status) {
            cursor = header_end;
            continue;
        }
        return Some(status);
    }
}

fn parse_status(block: &[u8]) -> Option<u16> {
    let text = std::str::from_utf8(block).ok()?;
    text.lines().find_map(|line| {
        let value = line.strip_prefix("HTTP/")?;
        value.split_whitespace().nth(1)?.parse().ok()
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
