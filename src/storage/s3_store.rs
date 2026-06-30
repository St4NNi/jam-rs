//! Objects on S3 and compatible stores, read with exact byte ranges through Apache OpenDAL.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{ReadAt, StorageError, checked_end};
use opendal::layers::{ConcurrentLimitLayer, RetryLayer, TimeoutLayer};
use opendal::services::S3;
use opendal::{ErrorKind, Operator};
use std::collections::HashMap;
use std::io;
use std::sync::{Arc, Mutex, PoisonError};
use std::time::Duration;
use tokio::runtime::Handle;

/// Limits on requests to the store.
#[derive(Clone, Copy, Debug)]
pub struct S3Limits {
    /// The longest one request may take, from its start to its last byte.
    pub request_timeout: Duration,
    /// Tries of one request, the first included. A read fails after the last one.
    pub tries: u32,
    /// The wait after the first failed try; it triples after each further one.
    pub backoff: Duration,
    /// Requests in flight at once, over all objects of one config. Some stores fail range reads
    /// when 128 run in parallel, so callers can set this lower.
    pub max_in_flight: usize,
}

impl Default for S3Limits {
    fn default() -> Self {
        Self {
            request_timeout: Duration::from_secs(60),
            tries: 4,
            backoff: Duration::from_millis(200),
            max_in_flight: 64,
        }
    }
}

/// The place of the store and the limits of its requests. Credentials come from the
/// environment or the AWS profile files. All objects of one config share its request limit.
#[derive(Clone)]
pub struct S3Config {
    region: String,
    endpoint: Option<String>,
    path_style: bool,
    limits: S3Limits,
    in_flight: ConcurrentLimitLayer,
    runtime: Handle,
    /// One operator per bucket, built at the first open.
    operators: Arc<Mutex<HashMap<String, Operator>>>,
    /// Sends requests without signing them (tests against a local store).
    #[cfg(test)]
    unsigned: bool,
}

impl S3Config {
    /// A config for `region`. With an `endpoint` (for example `https://s3.example.org`) the store
    /// is not AWS; `path_style` puts the bucket in the path. Reads wait on `runtime`.
    pub fn new(
        region: &str,
        endpoint: Option<&str>,
        path_style: bool,
        limits: S3Limits,
        runtime: Handle,
    ) -> Result<Self, StorageError> {
        if region.is_empty() {
            return Err(StorageError::Config("the region is empty".into()));
        }
        if limits.tries == 0 {
            return Err(StorageError::Config("at least one try is needed".into()));
        }
        if let Some(endpoint) = endpoint
            && !["http://", "https://"]
                .iter()
                .any(|scheme| endpoint.len() > scheme.len() && endpoint.starts_with(scheme))
        {
            return Err(StorageError::Config(format!(
                "{endpoint:?} is not an http or https endpoint"
            )));
        }
        opendal::install_default();
        Ok(Self {
            region: region.to_string(),
            endpoint: endpoint.map(str::to_string),
            path_style,
            limits,
            in_flight: ConcurrentLimitLayer::new(limits.max_in_flight.max(1)),
            runtime,
            operators: Arc::default(),
            #[cfg(test)]
            unsigned: false,
        })
    }

    /// The same config, sending requests without signatures.
    #[cfg(test)]
    pub(crate) fn unsigned(mut self) -> Self {
        self.unsigned = true;
        self
    }

    /// The operator of `bucket`, built once.
    fn operator(&self, bucket: &str) -> Result<Operator, StorageError> {
        let mut operators = self
            .operators
            .lock()
            .unwrap_or_else(PoisonError::into_inner);
        if let Some(operator) = operators.get(bucket) {
            return Ok(operator.clone());
        }
        let mut service = S3::default()
            .bucket(bucket)
            .region(&self.region)
            .disable_ec2_metadata();
        if let Some(endpoint) = &self.endpoint {
            service = service.endpoint(endpoint);
        }
        if !self.path_style {
            service = service.enable_virtual_host_style();
        }
        #[cfg(test)]
        if self.unsigned {
            service = service.skip_signature();
        }
        let limits = self.limits;
        let retry = RetryLayer::new()
            .with_max_times(limits.tries as usize - 1)
            .with_min_delay(limits.backoff)
            .with_factor(3.0)
            .with_max_delay(Duration::from_secs(60))
            .with_jitter();
        let timeout = TimeoutLayer::new()
            .with_timeout(limits.request_timeout)
            .with_io_timeout(limits.request_timeout);
        // Each try waits for a place in flight and is timed alone; the backoff runs outside.
        let operator = Operator::new(service)
            .map_err(|error| StorageError::Config(error.to_string()))?
            .layer(timeout)
            .layer(self.in_flight.clone())
            .layer(retry);
        operators.insert(bucket.to_string(), operator.clone());
        Ok(operator)
    }
}

/// An object `s3://bucket/key` with the size and ETag it had when it was opened. Every read is
/// one ranged GET with `If-Match`, so a changed object fails instead of returning other bytes.
pub struct S3Object {
    operator: Operator,
    runtime: Handle,
    key: String,
    size: u64,
    etag: String,
}

impl S3Object {
    /// Opens `s3://bucket/key`; if `expected_size` is given, the object must have that size.
    /// The size and ETag come from a GET of the first byte: some stores refuse HEAD requests.
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
        let operator = config.operator(bucket)?;
        let unusable = |source: io::Error| StorageError::Object {
            location: location.to_string(),
            source,
        };
        let metadata = config
            .runtime
            .block_on(async {
                let reader = operator.reader_with(key).await?;
                reader.read(0..1).await?;
                Ok::<_, opendal::Error>(reader.metadata().cloned())
            })
            .map_err(|error| unusable(error.into()))?
            .ok_or_else(|| unusable(invalid("the answer holds no size")))?;
        let size = metadata.content_length();
        let etag = metadata
            .etag()
            .filter(|etag| !etag.trim().is_empty() && !etag.trim_start().starts_with("W/"))
            .ok_or_else(|| unusable(invalid("the answer holds no strong ETag")))?;
        if let Some(expected) = expected_size
            && expected != size
        {
            return Err(unusable(invalid(&format!(
                "it has {size} bytes, not {expected}"
            ))));
        }
        Ok(Self {
            etag: etag.to_string(),
            operator,
            runtime: config.runtime.clone(),
            key: key.to_string(),
            size,
        })
    }
}

impl ReadAt for S3Object {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        let end = checked_end(offset, buffer.len(), self.size)?;
        if buffer.is_empty() {
            return Ok(());
        }
        // A range with an end needs no HEAD request for the object size.
        let read = self
            .operator
            .read_with(&self.key)
            .range(offset..end)
            .if_match(&self.etag);
        let bytes = self.runtime.block_on(read.into_future()).map_err(|error| {
            if error.kind() == ErrorKind::ConditionNotMatch {
                invalid("the S3 object changed since it was opened")
            } else {
                error.into()
            }
        })?;
        if bytes.len() != buffer.len() {
            return Err(invalid(&format!(
                "S3 returned {} bytes for a range of {}",
                bytes.len(),
                buffer.len()
            )));
        }
        buffer.copy_from_slice(&bytes.to_bytes());
        Ok(())
    }

    fn size(&self) -> u64 {
        self.size
    }
}

fn invalid(reason: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, reason.to_string())
}

#[cfg(test)]
#[path = "s3_store_tests.rs"]
mod tests;
