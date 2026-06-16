//! The frames of the protocol and how they travel: a 4-byte little-endian length, then one JSON
//! object.
//!
//! A server first sends `Hello`. A client then sends `Search` frames, one at a time. The answer
//! to a search is either `Refused`, or result frames (`Event`, `Summary`, `ShardFailed`) followed
//! by `Complete`. Each shard of the answer ends with exactly one `Summary` or `ShardFailed`.

use super::ServiceError;
use super::records::{EventRecord, SummaryRecord};
use super::settings::SearchSettings;
use crate::index::{Parameters, Shard};
use serde::{Deserialize, Serialize};
use std::io::{Read, Write};

/// The version of the protocol. Both sides must speak the same one.
pub const PROTOCOL_VERSION: u32 = 1;

/// The largest frame, so a damaged length cannot claim unbounded memory. A search frame holds the
/// query, so this also bounds the query length.
pub const MAX_FRAME_BYTES: u32 = 64 << 20;

/// What a server tells about one shard it answers for.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ShardIdentity {
    /// The name of the shard directory.
    pub name: String,
    /// The number of entries.
    pub entries: u64,
    pub parameters: Parameters,
    /// The SHA-256 digest of `manifest.json`, which tells shards apart.
    pub manifest_sha256: String,
}

impl ShardIdentity {
    /// The identity of `shard` with the given name and manifest digest.
    pub fn of(shard: &Shard, name: String, manifest_sha256: String) -> Self {
        Self {
            name,
            entries: shard.manifest().entries.len() as u64,
            parameters: shard.manifest().parameters,
            manifest_sha256,
        }
    }

    /// Whether a search of both shards finds the same seeds: same k, minimizer window, entropy
    /// limit and bin size. Bucket bits may differ, as they follow the shard size.
    pub fn same_seeds(&self, other: &Self) -> bool {
        let (ours, theirs) = (&self.parameters, &other.parameters);
        (
            ours.k,
            ours.window,
            ours.min_entropy.to_bits(),
            ours.bin_bits,
        ) == (
            theirs.k,
            theirs.window,
            theirs.min_entropy.to_bits(),
            theirs.bin_bits,
        )
    }
}

/// One search as a client asks for it.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SearchRequest {
    pub query_name: String,
    /// The query bases.
    pub query: String,
    pub circular: bool,
    pub settings: SearchSettings,
    /// Milliseconds from receipt until the search stops; the server's clock starts it, so the
    /// clocks of client and server need not agree.
    pub deadline_ms: Option<u64>,
}

/// A frame of the protocol. Result frames name their shard by its place in `Hello::shards`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "frame", rename_all = "snake_case")]
pub enum Frame {
    /// The first frame of every connection, from the server: the shards it answers for, its own
    /// first and its peers after.
    Hello {
        protocol_version: u32,
        shards: Vec<ShardIdentity>,
    },
    /// A search, from the client.
    Search {
        protocol_version: u32,
        request: SearchRequest,
    },
    /// An entry of `shard` changed its state.
    Event { shard: u32, record: EventRecord },
    /// The search of `shard` ended; its summary.
    Summary { shard: u32, summary: SummaryRecord },
    /// The search of `shard` failed or its peer was lost; its events so far stand.
    ShardFailed { shard: u32, message: String },
    /// The server did not run the search, such as when it is busy or the settings are invalid.
    Refused { message: String },
    /// The last frame of a search: complete when every shard searched every possible entry.
    Complete { complete: bool },
}

/// Fails unless `theirs` is the protocol version of this jam.
pub fn check_version(theirs: u32) -> Result<(), ServiceError> {
    match theirs {
        PROTOCOL_VERSION => Ok(()),
        _ => Err(ServiceError::Version {
            ours: PROTOCOL_VERSION,
            theirs,
        }),
    }
}

/// Writes one frame and flushes it. A frame over `MAX_FRAME_BYTES` is not sent.
pub fn write_frame(writer: &mut impl Write, frame: &Frame) -> Result<(), ServiceError> {
    let body = serde_json::to_vec(frame).map_err(|error| ServiceError::Frame(error.to_string()))?;
    let length = u32::try_from(body.len())
        .ok()
        .filter(|length| *length <= MAX_FRAME_BYTES)
        .ok_or_else(|| {
            ServiceError::Frame(format!("a frame of {} bytes is too large", body.len()))
        })?;
    writer.write_all(&length.to_le_bytes())?;
    writer.write_all(&body)?;
    writer.flush()?;
    Ok(())
}

/// Reads the next frame, or `None` when the other side closed the connection between frames.
pub fn read_frame(reader: &mut impl Read) -> Result<Option<Frame>, ServiceError> {
    let mut length = [0u8; 4];
    loop {
        match reader.read(&mut length[..1]) {
            Ok(0) => return Ok(None),
            Ok(_) => break,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error.into()),
        }
    }
    reader.read_exact(&mut length[1..])?;
    let length = u32::from_le_bytes(length);
    if length > MAX_FRAME_BYTES {
        return Err(ServiceError::Frame(format!(
            "a frame of {length} bytes is too large"
        )));
    }
    let mut body = vec![0; length as usize];
    reader.read_exact(&mut body)?;
    serde_json::from_slice(&body)
        .map(Some)
        .map_err(|error| ServiceError::Frame(error.to_string()))
}

#[cfg(test)]
#[path = "frame_tests.rs"]
mod tests;
