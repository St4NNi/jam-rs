//! The server: it accepts connections, says hello, and answers each search with frames streamed
//! as the entries change their state. With peers it is the coordinator of their shards.

use super::ServiceError;
use super::admission::{Admission, NoPlace, ServerLimits};
use super::client::FINAL_FRAMES;
use super::connection::{Connection, POLL, watch_disconnect};
use super::coordinator::{Peer, check_peers, search_all};
use super::frame::{Frame, PROTOCOL_VERSION, SearchRequest, ShardIdentity, check_version};
use super::local::ServedShard;
use crate::search::{SearchConfig, SearchControl};
use std::io::ErrorKind;
use std::net::{TcpListener, TcpStream};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant};

/// How long a connection may wait for its next search before the server closes it.
pub const IDLE: Duration = Duration::from_secs(60);

/// How a server runs searches.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ServerConfig {
    pub limits: ServerLimits,
    /// The search settings a client cannot choose (threads and memory limits), and defaults.
    pub search: SearchConfig,
}

/// A server of one shard, and the coordinator of its peers' shards when it has peers.
pub struct Server {
    local: ServedShard,
    peers: Vec<Peer>,
    shards: Vec<ShardIdentity>,
    search: SearchConfig,
    admission: Admission,
}

impl Server {
    /// A server of `local` alone.
    pub fn new(local: ServedShard, config: ServerConfig) -> Self {
        Self {
            shards: vec![local.identity().clone()],
            local,
            peers: Vec::new(),
            search: config.search,
            admission: Admission::new(config.limits),
        }
    }

    /// A coordinator of `local` and the shards of `peers`, which must fit the local shard.
    pub fn with_peers(
        local: ServedShard,
        config: ServerConfig,
        peers: Vec<Peer>,
    ) -> Result<Self, ServiceError> {
        check_peers(local.identity(), &peers)?;
        let mut server = Self::new(local, config);
        server
            .shards
            .extend(peers.iter().map(|peer| peer.identity().clone()));
        server.peers = peers;
        Ok(server)
    }

    /// The shards the server answers for: its own first, then its peers'.
    pub fn shards(&self) -> &[ShardIdentity] {
        &self.shards
    }

    /// Accepts connections until `stop` is set, each served in its own thread. A connection over
    /// the limit is refused and closed.
    pub fn serve(&self, listener: TcpListener, stop: &AtomicBool) -> Result<(), ServiceError> {
        listener.set_nonblocking(true)?;
        std::thread::scope(|scope| {
            while !stop.load(Ordering::Relaxed) {
                let stream = match listener.accept() {
                    Ok((stream, _)) => stream,
                    Err(error) if error.kind() == ErrorKind::WouldBlock => {
                        std::thread::sleep(POLL);
                        continue;
                    }
                    Err(error) if error.kind() == ErrorKind::ConnectionAborted => continue,
                    Err(error) => return Err(error.into()),
                };
                stream.set_nonblocking(false)?;
                let Some(place) = self.admission.connection() else {
                    refuse(stream, self.admission.limits().connections);
                    continue;
                };
                scope.spawn(move || {
                    let _place = place;
                    if let Err(error) = self.serve_connection(stream) {
                        eprintln!("jam serve: a connection ended: {error}");
                    }
                });
            }
            Ok(())
        })
    }

    /// Serves one connection: the hello, then its searches one after another until the client
    /// closes it or stays idle for `IDLE`.
    pub fn serve_connection(&self, stream: TcpStream) -> Result<(), ServiceError> {
        let mut connection = Connection::new(stream)?;
        connection.set_write_deadline(Some(Instant::now() + IDLE));
        let hello = Frame::Hello {
            protocol_version: PROTOCOL_VERSION,
            shards: self.shards.clone(),
        };
        // A connection closed before the hello, such as a readiness probe, is no error.
        if connection.write(&hello).is_err() {
            return Ok(());
        }
        loop {
            connection.set_read_deadline(Some(Instant::now() + IDLE));
            if !connection.wait_for_frame()? {
                return Ok(());
            }
            let Some(frame) = connection.read()? else {
                return Ok(());
            };
            let Frame::Search {
                protocol_version,
                request,
            } = frame
            else {
                return Err(ServiceError::Unexpected("a client sent no search".into()));
            };
            if let Err(error) = check_version(protocol_version) {
                let message = error.to_string();
                connection.write(&Frame::Refused { message })?;
                return Err(error);
            }
            self.answer(&mut connection, &request)?;
        }
    }

    /// Answers one search. The search stops at its deadline, and when the client closes the
    /// connection or stops reading.
    fn answer(
        &self,
        connection: &mut Connection,
        request: &SearchRequest,
    ) -> Result<(), ServiceError> {
        let deadline = request
            .deadline_ms
            .map(|ms| Instant::now() + Duration::from_millis(ms));
        connection.set_write_deadline(deadline.map(|deadline| deadline + FINAL_FRAMES));
        let stop = Arc::new(AtomicBool::new(false));
        let done = AtomicBool::new(false);
        let watched = connection.stream().try_clone()?;
        std::thread::scope(|scope| {
            scope.spawn(|| watch_disconnect(&watched, &done, &stop));
            let outcome = self.run(connection, request, deadline, &stop);
            done.store(true, Ordering::Relaxed);
            outcome
        })
    }

    /// Takes a place for the search, runs it on every shard and ends with `Complete`, or sends
    /// `Refused`.
    fn run(
        &self,
        connection: &mut Connection,
        request: &SearchRequest,
        deadline: Option<Instant>,
        stop: &Arc<AtomicBool>,
    ) -> Result<(), ServiceError> {
        let config = match request.settings.apply(&self.search) {
            Ok(config) => config,
            Err(error) => return refuse_search(connection, error.to_string()),
        };
        let control = SearchControl {
            deadline,
            stop: stop.as_ref(),
        };
        let _place = match self.admission.search(control) {
            Ok(place) => place,
            Err(NoPlace::Busy) => {
                let ServerLimits {
                    searches, queue, ..
                } = self.admission.limits();
                let message = format!("busy: {searches} searches run and {queue} wait");
                return refuse_search(connection, message);
            }
            Err(NoPlace::Stopped) => {
                let message = "the search stopped while it waited for a place".into();
                return refuse_search(connection, message);
            }
        };
        let mut emit = |frame: Frame| connection.write(&frame);
        let complete = search_all(
            &self.local,
            &self.peers,
            request,
            &config,
            deadline,
            stop,
            &mut emit,
        )?;
        connection.write(&Frame::Complete { complete })
    }
}

fn refuse_search(connection: &mut Connection, message: String) -> Result<(), ServiceError> {
    connection.write(&Frame::Refused { message })
}

/// Tells a connection over the limit why it is closed, without waiting for it.
fn refuse(stream: TcpStream, connections: usize) {
    let Ok(mut connection) = Connection::new(stream) else {
        return;
    };
    connection.set_write_deadline(Some(Instant::now() + Duration::from_secs(1)));
    let message = format!("the server has its {connections} connections open");
    let _ = connection.write(&Frame::Refused { message });
}

#[cfg(test)]
#[path = "server_tests.rs"]
mod tests;
