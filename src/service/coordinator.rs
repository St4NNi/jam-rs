//! The coordinator: one search on the local shard and on every peer at once, with their frames
//! merged as they arrive. Every shard runs in its own thread, so a slow peer holds back only its
//! own frames; a lost or failing peer becomes `ShardFailed`, never a failed search.

use super::ServiceError;
use super::client::Client;
use super::frame::{Frame, SearchRequest, ShardIdentity};
use super::local::ServedShard;
use crate::search::{SearchConfig, SearchControl};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{SyncSender, sync_channel};
use std::time::{Duration, Instant};

/// The most shards one coordinator answers for, its own included.
pub const MAX_SHARDS: usize = 64;

/// How long the last frames of a peer may take after the deadline. It is shorter than the grace
/// of the coordinator's own client, so a stalled peer is cut before the client gives up.
pub const PEER_GRACE: Duration = Duration::from_secs(5);

/// Frames held between the shards and the writer to the client. A full queue makes the shards
/// wait, so a slow client slows the search instead of filling RAM.
const QUEUED_FRAMES: usize = 256;

/// A peer server and the shard it answered for at start.
#[derive(Clone, Debug)]
pub struct Peer {
    address: String,
    identity: ShardIdentity,
    timeout: Duration,
}

impl Peer {
    /// Connects to the server at `address` once and pins its shard. `timeout` bounds this and
    /// every later connection to the peer.
    pub fn connect(address: &str, timeout: Duration) -> Result<Self, ServiceError> {
        let error = |message: String| ServiceError::Peer {
            address: address.to_string(),
            message,
        };
        let client = Client::connect(address, timeout).map_err(|e| error(e.to_string()))?;
        let [identity] = client.shards() else {
            return Err(error(format!(
                "it answers for {} shards; a peer answers for one",
                client.shards().len()
            )));
        };
        Ok(Self {
            address: address.to_string(),
            identity: identity.clone(),
            timeout,
        })
    }

    /// The address of the peer.
    pub fn address(&self) -> &str {
        &self.address
    }

    /// The shard the peer answered for at start.
    pub fn identity(&self) -> &ShardIdentity {
        &self.identity
    }
}

/// Checks that the peers fit the local shard: the same seeds, no shard twice, and at most
/// `MAX_SHARDS` shards.
pub fn check_peers(local: &ShardIdentity, peers: &[Peer]) -> Result<(), ServiceError> {
    if peers.len() >= MAX_SHARDS {
        return Err(ServiceError::Settings(format!(
            "{} peers; a coordinator answers for at most {MAX_SHARDS} shards",
            peers.len()
        )));
    }
    let mut digests = vec![&local.manifest_sha256];
    for peer in peers {
        let error = |message: &str| ServiceError::Peer {
            address: peer.address.clone(),
            message: message.to_string(),
        };
        if !local.same_seeds(&peer.identity) {
            return Err(error("its shard has other seed parameters"));
        }
        if digests.contains(&&peer.identity.manifest_sha256) {
            return Err(error("its shard is already served"));
        }
        digests.push(&peer.identity.manifest_sha256);
    }
    Ok(())
}

/// Runs one search on the local shard (tag 0) and on every peer (tags 1 and up), and passes each
/// frame to `emit` as it arrives. Each shard ends with one `Summary` or `ShardFailed`; the result
/// tells whether every shard completed. When `emit` fails, `stop` is set and its error returned.
pub fn search_all(
    local: &ServedShard,
    peers: &[Peer],
    request: &SearchRequest,
    config: &SearchConfig,
    deadline: Option<Instant>,
    stop: &Arc<AtomicBool>,
    emit: &mut dyn FnMut(Frame) -> Result<(), ServiceError>,
) -> Result<bool, ServiceError> {
    let (sender, receiver) = sync_channel(QUEUED_FRAMES);
    // Per shard: `None` while it runs, then whether it completed.
    let mut ended: Vec<Option<bool>> = vec![None; 1 + peers.len()];
    let mut failure = None;
    std::thread::scope(|scope| {
        let local_sender = sender.clone();
        scope.spawn(move || {
            let control = SearchControl {
                deadline,
                stop: stop.as_ref(),
            };
            let send = move |frame| {
                // A closed queue means the client left and `stop` is set.
                let _ = local_sender.send(frame);
            };
            let query = request.query.as_bytes();
            local.search(query, request.circular, config, control, 0, &send);
        });
        for (place, peer) in peers.iter().enumerate() {
            let sender = sender.clone();
            let tag = place as u32 + 1;
            scope.spawn(move || search_peer(peer, tag, request, deadline, stop, &sender));
        }
        drop(sender);
        for frame in receiver {
            match &frame {
                Frame::Summary { shard, summary } => {
                    ended[*shard as usize] = Some(summary.is_complete());
                }
                Frame::ShardFailed { shard, .. } => ended[*shard as usize] = Some(false),
                _ => {}
            }
            if let Err(error) = emit(frame) {
                stop.store(true, Ordering::Relaxed);
                failure = Some(error);
                break;
            }
        }
    });
    if let Some(error) = failure {
        return Err(error);
    }
    // A shard thread that ended without its last frame still gets one.
    for (shard, state) in ended.iter_mut().enumerate() {
        if state.is_none() {
            *state = Some(false);
            emit(Frame::ShardFailed {
                shard: shard as u32,
                message: "the shard ended without a summary".into(),
            })?;
        }
    }
    Ok(ended.iter().all(|state| *state == Some(true)))
}

/// Sends the search to one peer and forwards its frames with the shard renamed to `tag`. Any
/// failure before the peer's summary becomes `ShardFailed`.
fn search_peer(
    peer: &Peer,
    tag: u32,
    request: &SearchRequest,
    deadline: Option<Instant>,
    stop: &Arc<AtomicBool>,
    sender: &SyncSender<Frame>,
) {
    let mut ended = false;
    let outcome = forward_peer(peer, tag, request, deadline, stop, sender, &mut ended);
    if let (Err(error), false) = (outcome, ended) {
        let _ = sender.send(Frame::ShardFailed {
            shard: tag,
            message: format!("{}: {error}", peer.address),
        });
    }
}

/// The work of `search_peer`; `ended` is set once the peer's last frame was forwarded.
fn forward_peer(
    peer: &Peer,
    tag: u32,
    request: &SearchRequest,
    deadline: Option<Instant>,
    stop: &Arc<AtomicBool>,
    sender: &SyncSender<Frame>,
    ended: &mut bool,
) -> Result<(), ServiceError> {
    let mut client = Client::connect(&peer.address, peer.timeout)?;
    if client.shards() != std::slice::from_ref(&peer.identity) {
        return Err(ServiceError::Unexpected(
            "the peer now answers for another shard".into(),
        ));
    }
    client.set_grace(PEER_GRACE);
    let mut request = request.clone();
    request.deadline_ms = deadline.map(|deadline| {
        let left = deadline.saturating_duration_since(Instant::now());
        u64::try_from(left.as_millis()).unwrap_or(u64::MAX)
    });
    for frame in client.search(&request, Some(stop.clone()))? {
        let frame = match frame? {
            Frame::Event { record, .. } => Frame::Event { shard: tag, record },
            Frame::Summary { summary, .. } => Frame::Summary {
                shard: tag,
                summary,
            },
            Frame::ShardFailed { message, .. } => Frame::ShardFailed {
                shard: tag,
                message: format!("{}: {message}", peer.address),
            },
            Frame::Complete { .. } => return Ok(()),
            _ => {
                return Err(ServiceError::Unexpected(
                    "a frame a peer does not send".into(),
                ));
            }
        };
        *ended |= matches!(frame, Frame::Summary { .. } | Frame::ShardFailed { .. });
        if sender.send(frame).is_err() {
            return Ok(());
        }
    }
    Ok(())
}
