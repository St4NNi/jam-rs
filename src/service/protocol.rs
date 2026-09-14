//! The messages of a remote search and how they travel over one QUIC stream: a length, then
//! postcard.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::ServiceError;
use super::settings::SearchSettings;
use crate::report::{
    Confirmed, Dropped, Located, Possible, Record, ShardRecord, Summary, Unfinished,
};
use futures_util::{SinkExt, StreamExt};
use iroh::endpoint::{RecvStream, SendStream};
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use tokio_util::bytes::Bytes;
use tokio_util::codec::{FramedRead, FramedWrite, LengthDelimitedCodec};

/// The ALPN of the protocol. Its last part is the version, so a client and a server of other
/// versions do not connect.
pub const ALPN: &[u8] = b"jam/1";

/// The largest message, so a damaged length cannot claim unbounded memory. It also bounds the
/// query length.
pub const MAX_MESSAGE_BYTES: usize = 64 << 20;

/// A search: the only message a client sends on its stream.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SearchRequest {
    /// The name the records of the search carry.
    pub query_name: String,
    /// The query bases.
    pub bases: Vec<u8>,
    /// The query is circular.
    pub circular: bool,
    /// The settings the client chose.
    pub settings: SearchSettings,
    /// Milliseconds from receipt until the search stops; the server's clock starts it.
    pub deadline_ms: Option<u64>,
    /// Search the server's own shard only. A coordinator sets it for its peers, so peers that
    /// name each other do not loop.
    pub local_only: bool,
}

/// A message from the server: the records of every shard, each shard ending with one summary,
/// then `End`; or one `Refused`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Reply {
    /// A possible entry of a shard.
    Possible { shard: String, record: Possible },
    /// A located entry of a shard.
    Located { shard: String, record: Located },
    /// A confirmed entry of a shard.
    Confirmed { shard: String, record: Confirmed },
    /// A dropped entry of a shard.
    Dropped { shard: String, record: Dropped },
    /// An entry of a shard that the search did not finish.
    Unfinished { shard: String, record: Unfinished },
    /// The summary of a shard, its last record.
    Summary { shard: String, record: Box<Summary> },
    /// The server did not run the search: it is busy, or the settings are invalid.
    Refused { reason: String },
    /// The last message of a search.
    End,
}

impl Reply {
    /// The message that carries `record` of `shard`.
    pub fn record(shard: String, record: Record) -> Self {
        match record {
            Record::Possible(record) => Self::Possible { shard, record },
            Record::Located(record) => Self::Located { shard, record },
            Record::Confirmed(record) => Self::Confirmed { shard, record },
            Record::Dropped(record) => Self::Dropped { shard, record },
            Record::Unfinished(record) => Self::Unfinished { shard, record },
            Record::Summary(record) => Self::Summary { shard, record },
        }
    }

    /// The failed summary of `shard`, for a search that could not run or a peer that failed.
    pub fn failed(shard: String, request: &SearchRequest, error: &ServiceError) -> Self {
        let length = request.bases.len() as u64;
        Self::record(shard, Record::failed(&request.query_name, length, error))
    }

    /// The record of this message with its shard, as `jam search` prints it; `None` for a
    /// refusal or the end.
    pub fn into_record(self) -> Option<ShardRecord> {
        let (shard, record) = match self {
            Self::Possible { shard, record } => (shard, Record::Possible(record)),
            Self::Located { shard, record } => (shard, Record::Located(record)),
            Self::Confirmed { shard, record } => (shard, Record::Confirmed(record)),
            Self::Dropped { shard, record } => (shard, Record::Dropped(record)),
            Self::Unfinished { shard, record } => (shard, Record::Unfinished(record)),
            Self::Summary { shard, record } => (shard, Record::Summary(record)),
            Self::Refused { .. } | Self::End => return None,
        };
        Some(ShardRecord { shard, record })
    }
}

/// The sending half of a stream, one message per frame.
pub type MessageWriter = FramedWrite<SendStream, LengthDelimitedCodec>;

/// The receiving half of a stream, one message per frame.
pub type MessageReader = FramedRead<RecvStream, LengthDelimitedCodec>;

/// Frames both halves of a stream.
pub fn framed(send: SendStream, receive: RecvStream) -> (MessageWriter, MessageReader) {
    let codec = LengthDelimitedCodec::builder()
        .max_frame_length(MAX_MESSAGE_BYTES)
        .new_codec();
    (
        FramedWrite::new(send, codec.clone()),
        FramedRead::new(receive, codec),
    )
}

/// Writes one message.
pub async fn send(
    writer: &mut MessageWriter,
    message: &impl Serialize,
) -> Result<(), ServiceError> {
    let bytes = postcard::to_allocvec(message)?;
    writer.send(Bytes::from(bytes)).await?;
    Ok(())
}

/// Reads the next message, or `None` when the other side finished its stream.
pub async fn receive<T: DeserializeOwned>(
    reader: &mut MessageReader,
) -> Result<Option<T>, ServiceError> {
    match reader.next().await {
        Some(frame) => Ok(Some(postcard::from_bytes(&frame?)?)),
        None => Ok(None),
    }
}
