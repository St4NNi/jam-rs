use std::io;
use thiserror::Error;

/// Why a source could not be opened or configured. Reads themselves fail with `io::Error`.
#[derive(Debug, Error)]
pub enum StorageError {
    /// A local file or the local disk level failed.
    #[error("storage I/O failed: {0}")]
    Io(#[from] io::Error),
    /// The location is not `s3://bucket/key`, or the bucket or key is empty.
    #[error("not a usable S3 location: {0}")]
    Location(String),
    /// The S3 settings are not valid, for example an empty endpoint.
    #[error("invalid S3 setting: {0}")]
    Config(String),
    /// The object cannot be read, or does not answer as expected: another size, no strong ETag.
    #[error("S3 object {location} cannot be used: {source}")]
    Object {
        /// The `s3://` location that was opened.
        location: String,
        /// What went wrong.
        #[source]
        source: io::Error,
    },
}
