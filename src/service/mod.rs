//! Serving a shard over iroh (QUIC) and searching remote shards: messages, server, coordinator and
//! client.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

mod client;
mod coordinator;
mod endpoint;
mod protocol;
mod server;
mod settings;
mod shard;

use crate::index::IndexError;
use crate::loss::LossError;
use crate::report::ReportError;
use crate::search::SearchError;
use iroh::EndpointId;
use iroh::endpoint::{BindError, ClosedStream, ConnectError, ConnectionError, InvalidSocketAddr};
use thiserror::Error;

pub use client::{Answer, search};
pub use coordinator::PEER_GRACE;
pub use endpoint::{bind, parse_address};
pub use protocol::{ALPN, MAX_MESSAGE_BYTES, Reply, SearchRequest};
pub use server::{Server, ServerLimits};
pub use settings::{CoverageSetting, SearchSettings};
pub use shard::ServedShard;

/// Errors of serving and of remote searches.
#[derive(Debug, Error)]
pub enum ServiceError {
    /// Reading or writing a stream failed.
    #[error("reading or writing a stream failed: {0}")]
    Io(#[from] std::io::Error),
    /// A message cannot be encoded or decoded.
    #[error("bad message: {0}")]
    Message(#[from] postcard::Error),
    /// The endpoint cannot be bound.
    #[error("cannot bind the endpoint: {0}")]
    Bind(#[from] BindError),
    /// A bind address is not usable.
    #[error("invalid bind address: {0}")]
    BindAddress(#[from] InvalidSocketAddr),
    /// The server cannot be reached.
    #[error("cannot connect: {0}")]
    Connect(#[from] ConnectError),
    /// The connection was lost.
    #[error("the connection failed: {0}")]
    Connection(#[from] ConnectionError),
    /// A stream was used after it was closed.
    #[error("the stream is already closed")]
    Closed(#[from] ClosedStream),
    /// The other side finished its stream before the last message.
    #[error("the answer ended without its end message")]
    EndedEarly,
    /// The server did not run the search.
    #[error("the server refused the search: {0}")]
    Refused(String),
    /// A peer or server address is not `NODE_ID@ADDRESS`.
    #[error("{0} is not NODE_ID@ADDRESS")]
    Address(String),
    /// The server stops and takes no more searches.
    #[error("the server is stopping")]
    Stopping(#[from] tokio::sync::AcquireError),
    /// A peer passed the deadline and its grace.
    #[error("no answer before the deadline")]
    Timeout(#[from] tokio::time::error::Elapsed),
    /// A peer failed; its shard gets a failed summary.
    #[error("peer {peer}: {source}")]
    Peer {
        peer: EndpointId,
        source: Box<ServiceError>,
    },
    /// The search thread panicked.
    #[error("the search thread panicked")]
    Panic,
    /// Reading the shard failed.
    #[error(transparent)]
    Index(#[from] IndexError),
    /// The loss budget of a request is invalid.
    #[error(transparent)]
    Loss(#[from] LossError),
    /// The search could not run, or its settings are invalid.
    #[error(transparent)]
    Search(#[from] SearchError),
    /// A record cannot be built.
    #[error(transparent)]
    Report(#[from] ReportError),
}

