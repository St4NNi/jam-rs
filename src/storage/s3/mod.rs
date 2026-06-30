//! Objects on S3 and compatible stores, read with signed byte ranges.

mod config;
mod credentials;
mod object;
mod sign;

pub use config::{S3Config, S3Limits};
pub use credentials::Credentials;
pub use object::S3Object;

