//! Raw FASTA/FASTQ resource access for trace candidates.
//!
//! The indexed path is preferred by the runner.  This module is the explicit
//! fallback for a catalog row without JMA, and uses the same resource
//! abstraction so local paths, file URLs, HTTP(S), S3, and compressed parser
//! input have identical redaction and metric behavior.

use crate::resource::local::LocalResource;
use crate::resource::object::ObjectResource;
use crate::resource::{
    ByteRange, RangeReader, ResourceError, ResourceLocator, ResourceMetadata, ResourceMetrics,
    ResourceOpenOptions, ResourceResult, ResourceScheme,
};
use needletail::Sequence;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::io::{self, Read};
use thiserror::Error;

/// A concrete resource handle usable by both the raw parser and JMA reader.
/// The enum avoids making callers choose a different generic reader type for
/// local versus remote resources.
pub enum AssemblyResource {
    Local(LocalResource),
    Object(ObjectResource),
}

impl std::fmt::Debug for AssemblyResource {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("AssemblyResource")
            .field("locator", &self.locator().redacted())
            .finish()
    }
}

impl RangeReader for AssemblyResource {
    fn locator(&self) -> &ResourceLocator {
        match self {
            Self::Local(resource) => resource.locator(),
            Self::Object(resource) => resource.locator(),
        }
    }

    fn metadata(&self) -> ResourceResult<ResourceMetadata> {
        match self {
            Self::Local(resource) => resource.metadata(),
            Self::Object(resource) => resource.metadata(),
        }
    }

    fn read_range(&self, range: ByteRange) -> ResourceResult<Vec<u8>> {
        match self {
            Self::Local(resource) => resource.read_range(range),
            Self::Object(resource) => resource.read_range(range),
        }
    }

    fn stream(&self) -> ResourceResult<Box<dyn Read + Send>> {
        match self {
            Self::Local(resource) => resource.stream(),
            Self::Object(resource) => resource.stream(),
        }
    }

    fn metrics(&self) -> ResourceMetrics {
        match self {
            Self::Local(resource) => resource.metrics(),
            Self::Object(resource) => resource.metrics(),
        }
    }
}

/// Open any catalog resource using the shared checked resource contract.
pub fn open_resource(
    locator: impl AsRef<str>,
    options: ResourceOpenOptions,
) -> Result<AssemblyResource, RawError> {
    let locator = ResourceLocator::parse(locator.as_ref())?;
    match locator.scheme() {
        ResourceScheme::Local | ResourceScheme::File => Ok(AssemblyResource::Local(
            LocalResource::from_locator(locator, options)?,
        )),
        ResourceScheme::Http | ResourceScheme::Https | ResourceScheme::S3 => Ok(
            AssemblyResource::Object(ObjectResource::from_locator(locator, options)?),
        ),
    }
}

/// One assembly contig.  The sequence is normalized to uppercase and keeps
/// ambiguity symbols as `N`; identifiers are stable FASTA/FASTQ record IDs.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct RawContig {
    pub id: String,
    pub sequence: Vec<u8>,
}

impl RawContig {
    #[must_use]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// Materialized raw assembly and the metrics observed while parsing it.
/// Materialization is bounded by the caller's catalog/resource policy; the
/// runner never retains formatted output or candidate matrices for all rows.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct RawAssembly {
    pub redacted_locator: String,
    pub contigs: Vec<RawContig>,
    pub metrics: ResourceMetrics,
}

impl RawAssembly {
    /// Parse one FASTA/FASTQ resource, including parser-supported compression.
    pub fn open(locator: impl AsRef<str>, options: ResourceOpenOptions) -> Result<Self, RawError> {
        let resource = open_resource(locator, options)?;
        Self::from_resource(resource)
    }

    /// Parse an already opened resource.  This is also used by tests with
    /// mock resource implementations routed through the shared enum.
    pub fn from_resource(resource: AssemblyResource) -> Result<Self, RawError> {
        let redacted_locator = resource.locator().redacted();
        let stream = resource.stream()?;
        let contigs = parse_stream(stream, &redacted_locator)?;
        let metrics = resource.metrics();
        Ok(Self {
            redacted_locator,
            contigs,
            metrics,
        })
    }

    #[must_use]
    pub fn total_bases(&self) -> u64 {
        self.contigs
            .iter()
            .map(|contig| u64::try_from(contig.sequence.len()).unwrap_or(u64::MAX))
            .fold(0u64, u64::saturating_add)
    }

    #[must_use]
    pub fn contig(&self, id: &str) -> Option<&RawContig> {
        self.contigs.iter().find(|contig| contig.id == id)
    }
}

/// Parse FASTA/FASTQ from a generic stream.  A truly empty input is a valid
/// empty assembly, while malformed non-empty input remains an explicit error.
pub fn parse_stream(
    stream: impl Read + Send + 'static,
    redacted_locator: &str,
) -> Result<Vec<RawContig>, RawError> {
    let mut reader = match needletail::parse_fastx_reader(stream) {
        Ok(reader) => reader,
        Err(error) if error.kind == needletail::errors::ParseErrorKind::EmptyFile => {
            return Ok(Vec::new());
        }
        Err(error) => {
            return Err(RawError::Parse {
                locator: redacted_locator.to_string(),
                message: error.to_string(),
            });
        }
    };
    let mut contigs = Vec::new();
    let mut identifiers = HashSet::new();
    while let Some(record) = reader.next() {
        let record = record.map_err(|error| RawError::Parse {
            locator: redacted_locator.to_string(),
            message: error.to_string(),
        })?;
        let id = std::str::from_utf8(record.id())
            .unwrap_or("unknown")
            .split_whitespace()
            .next()
            .unwrap_or("unknown")
            .to_string();
        if id == "unknown" || !identifiers.insert(id.clone()) {
            return Err(RawError::Parse {
                locator: redacted_locator.to_string(),
                message: format!("duplicate or invalid contig identifier {id:?}"),
            });
        }
        let sequence = record.normalize(false).to_vec();
        if sequence.is_empty() {
            return Err(RawError::Parse {
                locator: redacted_locator.to_string(),
                message: format!("contig {id:?} has an empty sequence"),
            });
        }
        contigs.push(RawContig { id, sequence });
    }
    Ok(contigs)
}

#[derive(Debug, Error)]
pub enum RawError {
    #[error("raw resource error: {0}")]
    Resource(#[from] ResourceError),
    #[error("raw assembly parse failed for {locator}: {message}")]
    Parse { locator: String, message: String },
    #[error("raw assembly I/O failed for {locator}: {message}")]
    Io { locator: String, message: String },
}

impl From<io::Error> for RawError {
    fn from(error: io::Error) -> Self {
        Self::Io {
            locator: "<stream>".to_string(),
            message: error.to_string(),
        }
    }
}
