//! The server: it takes one search per QUIC stream, limits how many run, and streams the
//! messages of its own shard and its peers back.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::ServiceError;
use super::coordinator::{ShardSearch, Shards, start};
use super::protocol::{
    ALPN, MessageReader, MessageWriter, Reply, SearchRequest, framed, receive, send,
};
use super::shard::ServedShard;
use crate::search::SearchConfig;
use iroh::endpoint::Connection;
use iroh::protocol::{AcceptError, ProtocolHandler, Router};
use iroh::{Endpoint, EndpointAddr};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant};
use tokio::sync::{Semaphore, mpsc};

/// Messages held between the shards and the stream. A full queue makes the shards wait, so a
/// slow client slows the search instead of filling RAM.
const QUEUED_REPLIES: usize = 256;

/// The refusal of a search over the limits.
pub const BUSY: &str = "the server is busy";

/// The refusal of a search whose deadline passed while it waited for a running place.
pub const EXPIRED: &str = "the deadline passed while the search waited for a place";

/// How much work a server takes at once.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ServerLimits {
    /// Searches that run at once.
    pub searches: usize,
    /// Searches that wait for a place; a search over both limits is refused as busy.
    pub queue: usize,
}

impl Default for ServerLimits {
    fn default() -> Self {
        Self {
            searches: 4,
            queue: 32,
        }
    }
}

/// A server of one shard, and the coordinator of its peers' shards when it has peers.
pub struct Server {
    endpoint: Endpoint,
    shard: Arc<ServedShard>,
    peers: Vec<EndpointAddr>,
    /// The settings a client cannot choose, such as threads, and the defaults of the others.
    config: SearchConfig,
    running: Arc<Semaphore>,
    /// Running and waiting searches together.
    places: Semaphore,
}

impl std::fmt::Debug for Server {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let shard = self.shard.name();
        formatter
            .debug_struct("Server")
            .field("shard", &shard)
            .finish()
    }
}

impl Server {
    /// A server of `shard` on `endpoint` that also sends each search to `peers`.
    pub fn new(
        endpoint: Endpoint,
        shard: ServedShard,
        peers: Vec<EndpointAddr>,
        config: SearchConfig,
        limits: ServerLimits,
    ) -> Self {
        Self {
            endpoint,
            shard: Arc::new(shard),
            peers,
            config,
            running: Arc::new(Semaphore::new(limits.searches)),
            places: Semaphore::new(limits.searches + limits.queue),
        }
    }

    /// Starts answering on the endpoint. `Router::shutdown` stops it.
    pub fn spawn(self) -> Router {
        let endpoint = self.endpoint.clone();
        Router::builder(endpoint)
            .accept(ALPN, Arc::new(self))
            .spawn()
    }

    /// Answers the search on one stream. When the client stops the stream or leaves, the
    /// search stops.
    async fn answer(
        &self,
        mut writer: MessageWriter,
        mut reader: MessageReader,
    ) -> Result<(), ServiceError> {
        let gone = writer.get_ref().stopped();
        let request: SearchRequest = receive(&mut reader)
            .await?
            .ok_or(ServiceError::EndedEarly)?;
        let deadline = (request.deadline_ms).map(|ms| Instant::now() + Duration::from_millis(ms));
        let config = match request.settings.apply(&self.config) {
            Ok(config) => config,
            Err(error) => return refuse(writer, error.to_string()).await,
        };
        if let Some(note) = config.rule.identity_note()
            && !LOW_IDENTITY_LOGGED.swap(true, Ordering::Relaxed)
        {
            eprintln!("jam serve: {note}");
        }
        let Ok(_place) = self.places.try_acquire() else {
            return refuse(writer, BUSY.to_string()).await;
        };
        let search = Arc::new(ShardSearch {
            request,
            config,
            deadline,
            stop: Arc::new(AtomicBool::new(false)),
        });
        let outcome = tokio::select! {
            outcome = self.run(&search, &mut writer) => outcome,
            _ = gone => Ok(()),
        };
        // Ends the own shard's search if the client left; peers end as their messages find no
        // receiver.
        search.stop.store(true, Ordering::Relaxed);
        outcome
    }

    /// Waits for a running place until the deadline, starts the search on every shard and
    /// writes its messages. A search still waiting at its deadline is refused as expired.
    async fn run(
        &self,
        search: &Arc<ShardSearch>,
        writer: &mut MessageWriter,
    ) -> Result<(), ServiceError> {
        let place = self.running.clone().acquire_owned();
        let running = match search.deadline {
            Some(deadline) => match tokio::time::timeout_at(deadline.into(), place).await {
                Ok(running) => running?,
                Err(_) => {
                    let reason = EXPIRED.to_string();
                    send(writer, &Reply::Refused { reason }).await?;
                    writer.get_mut().finish()?;
                    return Ok(());
                }
            },
            None => place.await?,
        };
        let (replies, mut received) = mpsc::channel(QUEUED_REPLIES);
        let shards = Shards {
            endpoint: &self.endpoint,
            own: &self.shard,
            peers: &self.peers,
        };
        start(&shards, search, &replies, running);
        drop(replies);
        while let Some(reply) = received.recv().await {
            send(writer, &reply).await?;
        }
        send(writer, &Reply::End).await?;
        writer.get_mut().finish()?;
        Ok(())
    }
}

impl ProtocolHandler for Server {
    async fn accept(&self, connection: Connection) -> Result<(), AcceptError> {
        // The streams of one connection are answered in turn; a client opens one connection
        // per search.
        while let Ok((send_half, receive_half)) = connection.accept_bi().await {
            let (writer, reader) = framed(send_half, receive_half);
            if let Err(error) = self.answer(writer, reader).await {
                eprintln!("jam serve: {error}");
            }
        }
        Ok(())
    }
}

/// Whether the server has logged the note on identity limits below 0.8; it logs it once.
static LOW_IDENTITY_LOGGED: AtomicBool = AtomicBool::new(false);

async fn refuse(mut writer: MessageWriter, reason: String) -> Result<(), ServiceError> {
    send(&mut writer, &Reply::Refused { reason }).await?;
    writer.get_mut().finish()?;
    Ok(())
}
