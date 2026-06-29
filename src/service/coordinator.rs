//! The coordinator: one search on the own shard and on every peer at once, messages forwarded
//! as they arrive; a failed or slow peer becomes a failed summary of its shard alone.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::ServiceError;
use super::client;
use super::protocol::{Reply, SearchRequest};
use super::shard::ServedShard;
use crate::search::{SearchConfig, SearchControl};
use iroh::{Endpoint, EndpointAddr};
use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};
use tokio::sync::OwnedSemaphorePermit;
use tokio::sync::mpsc;

/// How long a peer may answer past the deadline before its shard counts as failed.
pub const PEER_GRACE: Duration = Duration::from_secs(5);

/// How long connecting to a peer may take; inside one cluster it takes milliseconds.
pub const PEER_CONNECT: Duration = Duration::from_secs(5);

/// One search as a server runs it.
pub struct ShardSearch {
    pub request: SearchRequest,
    /// The server's settings with the client's choices in place.
    pub config: SearchConfig,
    pub deadline: Option<Instant>,
    /// Set when the client left; the own shard's search then ends.
    pub stop: Arc<AtomicBool>,
}

/// The shards one search runs on: the own shard and, unless the request is local only, every
/// peer, dialled through `endpoint`.
pub struct Shards<'a> {
    pub endpoint: &'a Endpoint,
    pub own: &'a Arc<ServedShard>,
    pub peers: &'a [EndpointAddr],
}

/// Starts the search on every shard. Each task holds a clone of `replies`, so the channel closes
/// once every shard has sent its summary. `running` is held until the own shard's search ends.
pub fn start(
    shards: &Shards<'_>,
    search: &Arc<ShardSearch>,
    replies: &mpsc::Sender<Reply>,
    running: OwnedSemaphorePermit,
) {
    let peers = if search.request.local_only {
        &[]
    } else {
        shards.peers
    };
    for peer in peers {
        let (endpoint, peer) = (shards.endpoint.clone(), peer.clone());
        let (search, replies) = (search.clone(), replies.clone());
        tokio::spawn(async move { forward_peer(&endpoint, peer, &search, &replies).await });
    }
    let (shard, search, replies) = (shards.own.clone(), search.clone(), replies.clone());
    tokio::task::spawn_blocking(move || {
        let control = SearchControl {
            deadline: search.deadline,
            stop: &search.stop,
        };
        shard.search(&search.request, &search.config, control, &replies);
        drop(running);
    });
}

/// Forwards the messages of one peer. A peer that fails before its summary gets a failed summary,
/// named by its node id.
async fn forward_peer(
    endpoint: &Endpoint,
    peer: EndpointAddr,
    search: &ShardSearch,
    replies: &mpsc::Sender<Reply>,
) {
    let id = peer.id;
    if let Err(source) = relay(endpoint, peer, search, replies).await {
        let error = ServiceError::Peer {
            peer: id,
            source: Box::new(source),
        };
        let failed = Reply::failed(id.to_string(), &search.request, &error);
        // A closed channel means the client left.
        let _ = replies.send(failed).await;
    }
}

/// Sends the search to `peer` and passes its messages on. Each wait is bounded by the deadline
/// plus `PEER_GRACE`, and connecting also by `PEER_CONNECT`. An error after the peer's summary changes nothing, so it is dropped.
async fn relay(
    endpoint: &Endpoint,
    peer: EndpointAddr,
    search: &ShardSearch,
    replies: &mpsc::Sender<Reply>,
) -> Result<(), ServiceError> {
    let limit = search.deadline.map(|deadline| deadline + PEER_GRACE);
    let mut request = search.request.clone();
    request.local_only = true;
    request.deadline_ms = search.deadline.map(|deadline| {
        let left = deadline.saturating_duration_since(Instant::now());
        u64::try_from(left.as_millis()).unwrap_or(u64::MAX)
    });
    let connected = Instant::now() + PEER_CONNECT;
    let connect_limit = limit.map_or(connected, |limit| limit.min(connected));
    let connect = client::search(endpoint, peer, &request);
    let mut answer = within(Some(connect_limit), connect).await?;
    let mut summarized = false;
    loop {
        let reply = match within(limit, answer.next()).await {
            Ok(Some(reply)) => reply,
            Ok(None) => return Ok(()),
            Err(_) if summarized => return Ok(()),
            Err(error) => return Err(error),
        };
        summarized |= matches!(reply, Reply::Summary { .. });
        if replies.send(reply).await.is_err() {
            return Ok(());
        }
    }
}

/// Runs `work`, failing with `Timeout` when it is not done by `limit`.
async fn within<T>(
    limit: Option<Instant>,
    work: impl Future<Output = Result<T, ServiceError>>,
) -> Result<T, ServiceError> {
    match limit {
        Some(limit) => tokio::time::timeout_at(limit.into(), work).await?,
        None => work.await,
    }
}
