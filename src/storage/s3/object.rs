//! One object on S3, read with exact ranges. Opening sends one ranged GET for the first byte,
//! which gives the size and a strong ETag. A HEAD request is not used: the S3 service of the
//! project answers it with 403 where the same signature works for a GET. Every read is one
//! ranged GET with `If-Match`, whose status, content range, ETag and length are checked, so a
//! changed or truncated object fails instead of returning other bytes.
//!
//! Each try has the request timeout of the config. A try that fails for a temporary reason (no
//! connection, a timeout, a broken body, status 408, 429, 500, 502, 503 or 504) is repeated
//! after a growing wait, up to the tries of the config. Anything else fails at once. At most
//! `max_in_flight` requests of one config run at the same time.

use super::config::{S3Config, Target};
use crate::storage::{ReadAt, StorageError, checked_end};
use std::io::{self, Read};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Duration;

/// Statuses after which the same request may succeed.
fn transient(status: u16) -> bool {
    matches!(status, 408 | 429 | 500 | 502 | 503 | 504)
}

/// How one try of a request failed.
enum Failure {
    /// Another try may fix it (connection, temporary status, broken body).
    Retry(io::Error),
    /// It will not change (another object version or range, a body of the wrong size).
    Final(io::Error),
}

/// An object `s3://bucket/key` with the size and ETag it had when it was opened.
pub struct S3Object {
    config: S3Config,
    target: Target,
    size: u64,
    etag: String,
    requests: AtomicU64,
}

impl S3Object {
    /// Opens `s3://bucket/key`. If `expected_size` is given, the object must have that size.
    pub fn open(
        config: &S3Config,
        location: &str,
        expected_size: Option<u64>,
    ) -> Result<Self, StorageError> {
        let (bucket, key) = location
            .strip_prefix("s3://")
            .and_then(|rest| rest.split_once('/'))
            .filter(|(bucket, key)| !bucket.is_empty() && !key.is_empty())
            .ok_or_else(|| StorageError::Location(location.to_string()))?;
        let target = config.target(bucket, key);
        let unusable = |reason: String| StorageError::Object {
            location: location.to_string(),
            reason,
        };
        let (size, etag) = with_tries(config, || first_byte(config, &target))
            .map_err(|error| unusable(error.to_string()))?;
        if expected_size.is_some_and(|expected| expected != size) {
            return Err(unusable(format!(
                "it has {size} bytes, not {}",
                expected_size.unwrap_or(0)
            )));
        }
        Ok(Self {
            config: config.clone(),
            target,
            size,
            etag,
            requests: AtomicU64::new(1),
        })
    }

    /// Requests sent so far, the one that opened the object included.
    pub fn requests(&self) -> u64 {
        self.requests.load(Ordering::Relaxed)
    }

    /// One try to read `buffer.len()` bytes at `offset`.
    fn try_range(&self, offset: u64, buffer: &mut [u8]) -> Result<(), Failure> {
        let last = offset + buffer.len() as u64 - 1;
        let _permit = self.config.in_flight().acquire();
        let response = self
            .config
            .get_range(&self.target, (offset, last), Some(&self.etag))
            .map_err(Failure::Retry)?;
        self.requests.fetch_add(1, Ordering::Relaxed);
        let header = |name: &str| response.header(name);
        let status = response.status();
        let whole = offset == 0 && buffer.len() as u64 == self.size;
        let range = format!("bytes {offset}-{last}/{}", self.size);
        let valid = (status == 206 || (status == 200 && whole))
            && (status != 206 || header("content-range") == Some(range.as_str()))
            && header("etag") == Some(self.etag.as_str())
            && header("content-length")
                .is_none_or(|value| value.parse::<usize>().ok() == Some(buffer.len()));
        if !valid {
            let error = io::Error::new(
                io::ErrorKind::InvalidData,
                format!("S3 returned another range or object version (status {status})"),
            );
            return Err(if transient(status) {
                Failure::Retry(error)
            } else {
                Failure::Final(error)
            });
        }
        let mut body = response.into_reader();
        body.read_exact(buffer).map_err(Failure::Retry)?;
        match body.read(&mut [0; 1]) {
            Ok(0) => Ok(()),
            Ok(_) => Err(Failure::Final(io::Error::new(
                io::ErrorKind::InvalidData,
                "S3 range body exceeds its declared span",
            ))),
            Err(error) => Err(Failure::Retry(error)),
        }
    }
}

impl ReadAt for S3Object {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        checked_end(offset, buffer.len(), self.size)?;
        if buffer.is_empty() {
            return Ok(());
        }
        with_tries(&self.config, || self.try_range(offset, buffer))
    }

    fn size(&self) -> u64 {
        self.size
    }
}

/// One try to learn the size and the ETag from the first byte of the object.
fn first_byte(config: &S3Config, target: &Target) -> Result<(u64, String), Failure> {
    let _permit = config.in_flight().acquire();
    let response = config
        .get_range(target, (0, 0), None)
        .map_err(Failure::Retry)?;
    let status = response.status();
    let invalid =
        |reason: &str| Failure::Final(io::Error::other(format!("{reason} (status {status})")));
    if transient(status) {
        return Err(Failure::Retry(io::Error::other(format!("status {status}"))));
    }
    let size = response
        .header("content-range")
        .and_then(|range| range.strip_prefix("bytes 0-0/"))
        .and_then(|total| total.parse::<u64>().ok());
    let etag = response
        .header("etag")
        .filter(|etag| !etag.trim().is_empty() && !etag.trim_start().starts_with("W/"));
    match (status, size, etag) {
        (206, Some(size), Some(etag)) => Ok((size, etag.to_string())),
        (206, _, _) => Err(invalid("the answer has no total size or no strong ETag")),
        _ => Err(invalid("the store did not answer a range request")),
    }
}

/// Runs `attempt` until it succeeds, fails for good or has used the tries of the config. The
/// wait between tries starts at `backoff` and triples each time.
fn with_tries<T>(
    config: &S3Config,
    mut attempt: impl FnMut() -> Result<T, Failure>,
) -> io::Result<T> {
    let limits = config.limits();
    let mut wait = limits.backoff;
    let mut tries = 0;
    loop {
        tries += 1;
        match attempt() {
            Ok(value) => return Ok(value),
            Err(Failure::Retry(_)) if tries < limits.tries => {
                std::thread::sleep(wait);
                wait = wait.saturating_mul(3).min(Duration::from_secs(60));
            }
            Err(Failure::Retry(error) | Failure::Final(error)) => return Err(error),
        }
    }
}
