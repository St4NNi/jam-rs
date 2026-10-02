//! Objects on S3 and compatible stores, read with exact byte ranges through Apache OpenDAL.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::counts::{count_failed, count_s3};
use super::etags::{EtagBook, Latch, Start, changed_error, version_of};
use super::{ReadAt, StorageError, checked_end, lock};
use opendal::Operator;
use opendal::layers::{ConcurrentLimitLayer, RetryEvent, RetryLayer, TimeoutLayer};
use opendal::services::S3;
use std::collections::HashMap;
use std::io;
use std::ops::Range;
use std::path::Path;
use std::sync::{Arc, Mutex};
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
            max_in_flight: 256,
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
    /// The ETags of the objects read so far.
    etags: Arc<EtagBook>,
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
            etags: Arc::default(),
            #[cfg(test)]
            unsigned: false,
        })
    }

    /// The config of the standard AWS variables: `AWS_REGION` (default `us-east-1`) and
    /// `AWS_ENDPOINT_URL`, a store that is not AWS and takes the bucket in the path.
    pub fn from_environment(limits: S3Limits, runtime: Handle) -> Result<Self, StorageError> {
        let variable = |name| std::env::var(name).ok().filter(|value| !value.is_empty());
        let region = variable("AWS_REGION").unwrap_or_else(|| "us-east-1".into());
        let endpoint = variable("AWS_ENDPOINT_URL");
        Self::new(
            &region,
            endpoint.as_deref(),
            endpoint.is_some(),
            limits,
            runtime,
        )
    }

    /// Keeps the ETags of the objects in `directory` as well, so a later run knows which version
    /// its cached blocks belong to. Call it before the first read.
    pub fn keep_etags(&self, directory: &Path) -> Result<(), StorageError> {
        self.etags
            .keep(directory)
            .map_err(|source| StorageError::File {
                path: directory.to_path_buf(),
                source,
            })
    }

    /// Forgets the ETag of the object `s3://bucket/key`. It is learned again at the next read.
    pub fn forget_etag(&self, location: &str) {
        self.etags.forget(location);
    }

    /// How often an object was found replaced while it was read.
    pub fn changes(&self) -> u64 {
        self.etags.changes()
    }

    /// The same config, sending requests without signatures.
    #[cfg(test)]
    pub(crate) fn unsigned(mut self) -> Self {
        self.unsigned = true;
        self
    }

    /// The operator of `bucket`, built once.
    fn operator(&self, bucket: &str) -> Result<Operator, StorageError> {
        let mut operators = lock(&self.operators);
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
            .with_jitter()
            // Each retry follows one failed attempt; the last attempt is counted by the reader.
            .with_notify(|_: RetryEvent<'_>| count_failed());
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

/// An object `s3://bucket/key` of a known size. Every read is one ranged GET; the first learns
/// the ETag and later ones only succeed while the object has it. Opening sends no request, as
/// the size comes from the manifest and some stores refuse HEAD.
pub struct S3Object {
    operator: Operator,
    runtime: Handle,
    etags: Arc<EtagBook>,
    location: String,
    key: String,
    size: u64,
    /// Called with the version of the old ETag when the object was found replaced.
    on_change: Option<ChangeHook>,
}

/// What to do when an object is found replaced.
pub type ChangeHook = Arc<dyn Fn(u64) + Send + Sync>;

impl S3Object {
    /// Opens `s3://bucket/key`, an object of `size` bytes.
    pub fn open(config: &S3Config, location: &str, size: u64) -> Result<Self, StorageError> {
        let (bucket, key) = split_location(location)?;
        Ok(Self {
            operator: config.operator(bucket)?,
            runtime: config.runtime.clone(),
            etags: Arc::clone(&config.etags),
            location: location.to_string(),
            key: key.to_string(),
            size,
            on_change: None,
        })
    }

    /// The same object, calling `hook` when it is found replaced.
    pub fn on_change(mut self, hook: ChangeHook) -> Self {
        self.on_change = Some(hook);
        self
    }

    /// One GET of `range`. It sends the known ETag as `If-Match`; without one it learns the ETag
    /// from the answer.
    fn get(&self, range: Range<u64>) -> io::Result<Vec<u8>> {
        let known = self.etags.get(&self.location);
        let key = &self.key;
        let read = async {
            if let Some(etag) = &known {
                let read = (self.operator.read_with(key)).range(range).if_match(etag);
                return Ok((read.await?.to_vec(), None));
            }
            let reader = self.operator.reader(key).await?;
            let bytes = reader.read(range).await?.to_vec();
            let etag = (reader.metadata()).and_then(|meta| meta.etag().map(str::to_string));
            Ok::<_, opendal::Error>((bytes, etag))
        };
        let (bytes, seen) = self.runtime.block_on(read).map_err(|error| {
            count_failed();
            match error.kind() {
                opendal::ErrorKind::ConditionNotMatch => self.changed(known.as_deref()),
                _ => error.into(),
            }
        })?;
        count_s3(bytes.len() as u64);
        if known.is_none() {
            let seen = seen.ok_or_else(|| invalid(&format!("{} has no ETag", self.location)))?;
            if let Err(old) = self.etags.settle(&self.location, &seen) {
                return Err(self.changed(Some(&old)));
            }
        }
        Ok(bytes)
    }

    /// Takes `etag`, the ETag that a local copy of the object was read under, as the known one
    /// if none is known. Whether it is the known one.
    pub fn adopt(&self, etag: &str) -> bool {
        self.etags.settle(&self.location, etag).is_ok()
    }

    /// Forgets the known ETag and learns the one that the object has now, with a GET of its
    /// first byte.
    pub fn refresh(&self) -> io::Result<Option<String>> {
        self.etags.forget(&self.location);
        if self.size > 0 {
            self.get(0..1)?;
        }
        Ok(self.etags.get(&self.location))
    }

    /// Fails if the object no longer has `etag`, with a GET of its first byte that must match.
    pub fn verify(&self, etag: &str) -> io::Result<()> {
        if !self.adopt(etag) {
            return Err(self.changed(self.etags.get(&self.location).as_deref()));
        }
        self.get(0..1).map(drop)
    }

    /// Marks the check of this object as passed, after reads that its known ETag bound.
    pub(crate) fn pass_check(&self) {
        self.etags.pass_check(&self.location);
    }

    /// Asks for the check of this object.
    pub(crate) fn begin_check(&self) -> Start {
        self.etags.begin_check(&self.location)
    }

    /// Whether the ETag of the object is known.
    pub(crate) fn is_known(&self) -> bool {
        self.etags.version(&self.location).is_some()
    }

    /// How often an object of the config was found replaced.
    pub(crate) fn changes(&self) -> u64 {
        self.etags.changes()
    }

    /// Runs the check that `latch` waits for: a GET of the first byte. With a known ETag it runs
    /// beside the caller; else the caller waits, as the GET teaches the ETag.
    pub(crate) fn run_check(self: &Arc<Self>, latch: &Arc<Latch>) {
        let (this, latch) = (Arc::clone(self), Arc::clone(latch));
        let check = move || {
            let result = match this.size {
                0 => Ok(()),
                _ => this.get(0..1).map(drop),
            };
            this.etags.end_check(&this.location, &latch, &result);
        };
        match self.is_known() {
            true => drop(self.runtime.spawn_blocking(check)),
            false => check(),
        }
    }

    /// Forgets the ETag `old` of the object and gives the error of a replaced object.
    fn changed(&self, old: Option<&str>) -> io::Error {
        self.etags.changed(&self.location, old);
        if let (Some(hook), Some(old)) = (&self.on_change, old) {
            hook(version_of(old));
        }
        changed_error()
    }

    /// All bytes of the small object `s3://bucket/key`, read with one GET.
    pub fn read_whole(config: &S3Config, location: &str) -> Result<Vec<u8>, StorageError> {
        let (bucket, key) = split_location(location)?;
        let operator = config.operator(bucket)?;
        let read = config.runtime.block_on(operator.read(key));
        let bytes = read.map_err(|error| {
            count_failed();
            StorageError::Object {
                location: location.to_string(),
                source: error.into(),
            }
        })?;
        count_s3(bytes.len() as u64);
        Ok(bytes.to_vec())
    }
}

/// The bucket and key of `s3://bucket/key`.
fn split_location(location: &str) -> Result<(&str, &str), StorageError> {
    location
        .strip_prefix("s3://")
        .and_then(|rest| rest.split_once('/'))
        .filter(|(bucket, key)| !bucket.is_empty() && !key.is_empty())
        .ok_or_else(|| StorageError::Location(location.to_string()))
}

impl ReadAt for S3Object {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        let end = checked_end(offset, buffer.len(), self.size)?;
        if buffer.is_empty() {
            return Ok(());
        }
        // A range with an end needs no HEAD request for the object size.
        let bytes = self.get(offset..end)?;
        if bytes.len() != buffer.len() {
            return Err(invalid(&format!(
                "S3 returned {} bytes for a range of {}",
                bytes.len(),
                buffer.len()
            )));
        }
        buffer.copy_from_slice(&bytes);
        Ok(())
    }

    fn size(&self) -> u64 {
        self.size
    }

    fn version(&self) -> io::Result<u64> {
        if self.size == 0 {
            return Ok(0);
        }
        // The first byte teaches the ETag if no read has yet.
        if let Some(version) = self.etags.version(&self.location) {
            return Ok(version);
        }
        self.get(0..1)?;
        (self.etags.version(&self.location)).ok_or_else(changed_error)
    }
}

fn invalid(reason: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, reason.to_string())
}

#[cfg(test)]
#[path = "s3_tests.rs"]
mod tests;
