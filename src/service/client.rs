//! The client side: one connection to a server, which sends searches and yields their frames.

use super::ServiceError;
use super::connection::Connection;
use super::frame::{Frame, PROTOCOL_VERSION, SearchRequest, ShardIdentity, check_version};
use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};

/// How long the last frames of a search may take after its deadline: the server stops at the
/// deadline and then still sends the summaries and `Complete`.
pub const FINAL_FRAMES: Duration = Duration::from_secs(10);

/// A connection to a server, after its hello.
pub struct Client {
    connection: Connection,
    shards: Vec<ShardIdentity>,
    grace: Duration,
}

impl Client {
    /// Connects to `address` (host and port) and reads the hello, both within `timeout`.
    pub fn connect(address: &str, timeout: Duration) -> Result<Self, ServiceError> {
        let mut connection = Connection::open(address, timeout)?;
        connection.set_read_deadline(Some(Instant::now() + timeout));
        let shards = match connection.read()? {
            Some(Frame::Hello {
                protocol_version,
                shards,
            }) => {
                check_version(protocol_version)?;
                shards
            }
            Some(Frame::Refused { message }) => return Err(ServiceError::Refused(message)),
            Some(_) => return Err(ServiceError::Unexpected("the server sent no hello".into())),
            None => return Err(ServiceError::Unexpected("the server closed at once".into())),
        };
        if shards.is_empty() || u32::try_from(shards.len()).is_err() {
            return Err(ServiceError::Unexpected(
                "the server names no shards".into(),
            ));
        }
        Ok(Self {
            connection,
            shards,
            grace: FINAL_FRAMES,
        })
    }

    /// The shards the server answers for; result frames name them by their place here.
    pub fn shards(&self) -> &[ShardIdentity] {
        &self.shards
    }

    /// Sets how long the last frames may take after the deadline (`FINAL_FRAMES` by default).
    pub fn set_grace(&mut self, grace: Duration) {
        self.grace = grace;
    }

    /// Sends one search and returns its frames. Reads end at the deadline plus the grace, and
    /// soon after `stop` is set; dropping the frames early leaves the connection unusable.
    pub fn search(
        &mut self,
        request: &SearchRequest,
        stop: Option<Arc<AtomicBool>>,
    ) -> Result<SearchFrames<'_>, ServiceError> {
        let deadline = request
            .deadline_ms
            .map(|ms| Instant::now() + Duration::from_millis(ms) + self.grace);
        self.connection.set_read_deadline(deadline);
        self.connection.set_write_deadline(deadline);
        self.connection.set_stop(stop);
        self.connection.write(&Frame::Search {
            protocol_version: PROTOCOL_VERSION,
            request: request.clone(),
        })?;
        Ok(SearchFrames {
            ended: vec![false; self.shards.len()],
            complete: vec![false; self.shards.len()],
            entries: self.shards.iter().map(|shard| shard.entries).collect(),
            connection: &mut self.connection,
            done: false,
        })
    }
}

/// The frames of one search, checked as they arrive: every shard ends once, with a summary or a
/// failure, and `Complete` comes last. `Refused` and broken frames end the search with an error.
pub struct SearchFrames<'c> {
    connection: &'c mut Connection,
    /// Per shard: whether it ended, and whether it ended complete.
    ended: Vec<bool>,
    complete: Vec<bool>,
    entries: Vec<u64>,
    done: bool,
}

impl SearchFrames<'_> {
    /// Checks one frame against what came before.
    fn check(&mut self, frame: &Frame) -> Result<(), ServiceError> {
        let open = |shard: u32| -> Result<usize, ServiceError> {
            let place = shard as usize;
            match self.ended.get(place) {
                Some(false) => Ok(place),
                _ => Err(ServiceError::Unexpected(format!(
                    "a frame of shard {shard}, which is unknown or has ended"
                ))),
            }
        };
        match frame {
            Frame::Event { shard, record } => {
                let place = open(*shard)?;
                if u64::from(record.entry) >= self.entries[place] {
                    return Err(ServiceError::Unexpected(format!(
                        "entry {} is not in shard {shard}",
                        record.entry
                    )));
                }
            }
            Frame::Summary { shard, summary } => {
                let place = open(*shard)?;
                self.ended[place] = true;
                self.complete[place] = summary.is_complete();
            }
            Frame::ShardFailed { shard, .. } => {
                let place = open(*shard)?;
                self.ended[place] = true;
            }
            Frame::Complete { complete } => {
                let all = self.complete.iter().all(|done| *done);
                if !self.ended.iter().all(|ended| *ended) || (*complete && !all) {
                    return Err(ServiceError::Unexpected(
                        "complete before every shard ended, or complete with a stopped shard"
                            .into(),
                    ));
                }
                self.done = true;
            }
            Frame::Refused { message } => return Err(ServiceError::Refused(message.clone())),
            Frame::Hello { .. } | Frame::Search { .. } => {
                return Err(ServiceError::Unexpected(
                    "a frame only a client sends".into(),
                ));
            }
        }
        Ok(())
    }
}

impl Iterator for SearchFrames<'_> {
    type Item = Result<Frame, ServiceError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let checked = match self.connection.read() {
            Ok(Some(frame)) => self.check(&frame).map(|()| frame),
            Ok(None) => Err(ServiceError::Unexpected(
                "the server closed the connection".into(),
            )),
            Err(error) => Err(error),
        };
        self.done |= checked.is_err();
        Some(checked)
    }
}
