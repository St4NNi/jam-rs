//! Remote searches: one search sent to a server, its messages read as they arrive.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::ServiceError;
use super::protocol::{ALPN, MessageReader, Reply, SearchRequest, framed, receive, send};
use iroh::endpoint::Connection;
use iroh::{Endpoint, EndpointAddr};

/// The messages of one search. Dropping it stops the stream, which stops the search.
pub struct Answer {
    reader: MessageReader,
    _connection: Connection,
}

/// Sends `request` to the server at `server` on its own connection and stream.
pub async fn search(
    endpoint: &Endpoint,
    server: EndpointAddr,
    request: &SearchRequest,
) -> Result<Answer, ServiceError> {
    let connection = endpoint.connect(server, ALPN).await?;
    let (send_half, receive_half) = connection.open_bi().await?;
    let (mut writer, reader) = framed(send_half, receive_half);
    send(&mut writer, request).await?;
    writer.get_mut().finish()?;
    Ok(Answer {
        reader,
        _connection: connection,
    })
}

impl Answer {
    /// The next record message, or `None` after the end. A refusal, or a stream that ends
    /// without its end message, is an error.
    pub async fn next(&mut self) -> Result<Option<Reply>, ServiceError> {
        match receive(&mut self.reader).await? {
            Some(Reply::End) => Ok(None),
            Some(Reply::Refused { reason }) => Err(ServiceError::Refused(reason)),
            Some(reply) => Ok(Some(reply)),
            None => Err(ServiceError::EndedEarly),
        }
    }
}
