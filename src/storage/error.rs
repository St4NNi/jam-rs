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
    /// No usable static credentials were found; roles, SSO and processes are not supported.
    #[error("no usable S3 credentials: {0}")]
    Credentials(String),
    /// The S3 settings are not valid, for example an empty endpoint.
    #[error("invalid S3 setting: {0}")]
    Config(String),
    /// The object does not answer as expected: another size, no strong ETag or no ranges.
    #[error("S3 object {location} cannot be used: {reason}")]
    Object {
        /// The `s3://` location that was opened.
        location: String,
        /// What is wrong, in plain words.
        reason: String,
    },
}
