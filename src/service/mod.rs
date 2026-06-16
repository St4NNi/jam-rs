//! Serving a shard over TCP and searching remote shards: protocol, server, client and coordinator.
//!
//! A server opens one shard and answers searches on it. A coordinator is a server that also knows
//! its peers: it sends each search to every peer and to its own shard, and merges their frames as
//! they arrive. Every answer is a stream of frames that ends with one summary or failure per shard
//! and `Complete`.

pub mod frame;
pub mod records;
pub mod settings;

use crate::index::IndexError;
use crate::search::SearchError;
use thiserror::Error;

pub use frame::{Frame, MAX_FRAME_BYTES, PROTOCOL_VERSION, SearchRequest, ShardIdentity};
pub use records::{EventRecord, SummaryRecord};
pub use settings::SearchSettings;

/// Errors of serving and of remote searches.
#[derive(Debug, Error)]
pub enum ServiceError {
    /// Reading or writing a connection failed.
    #[error("{0}")]
    Io(#[from] std::io::Error),
    /// A frame could not be read or written.
    #[error("bad frame: {0}")]
    Frame(String),
    /// The other side speaks another protocol version.
    #[error("protocol version {theirs} differs from version {ours} of this jam")]
    Version { ours: u32, theirs: u32 },
    /// A frame came that does not fit the conversation.
    #[error("unexpected frame: {0}")]
    Unexpected(String),
    /// A search setting is invalid.
    #[error("invalid search setting: {0}")]
    Settings(String),
    /// The server refused the search.
    #[error("the server refused the search: {0}")]
    Refused(String),
    /// A peer does not fit the shards of this coordinator.
    #[error("peer {address}: {message}")]
    Peer { address: String, message: String },
    /// Opening the shard failed.
    #[error(transparent)]
    Index(#[from] IndexError),
    /// The search could not start.
    #[error(transparent)]
    Search(#[from] SearchError),
}
