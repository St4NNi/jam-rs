//! Small helpers that keep resource errors checked and redacted.

use super::{ByteRange, ResourceError, ResourceLocator, ResourceResult};
use std::io;

pub(crate) fn validate_range(
    locator: &ResourceLocator,
    range: ByteRange,
    size: u64,
) -> ResourceResult<()> {
    let end = range.end()?;
    if end > size {
        return Err(ResourceError::RangeOutOfBounds {
            offset: range.offset,
            length: range.length,
            size,
        });
    }
    usize::try_from(range.length).map_err(|_| ResourceError::Io {
        locator: locator.redacted(),
        message: "requested range is too large for this platform".to_string(),
    })?;
    Ok(())
}

pub(crate) fn io_error(locator: &ResourceLocator, error: io::Error) -> ResourceError {
    ResourceError::Io {
        locator: locator.redacted(),
        message: error.to_string(),
    }
}

pub(crate) fn transport_error(
    locator: &ResourceLocator,
    message: impl Into<String>,
) -> ResourceError {
    ResourceError::Transport {
        locator: locator.redacted(),
        message: redact_message(&message.into()),
    }
}

pub(crate) fn http_status_error(locator: &ResourceLocator, status: u16) -> ResourceError {
    ResourceError::HttpStatus {
        locator: locator.redacted(),
        status,
    }
}

pub(crate) fn timeout_error(locator: &ResourceLocator) -> ResourceError {
    ResourceError::Timeout {
        locator: locator.redacted(),
    }
}

pub(crate) fn is_retryable(error: &ResourceError) -> bool {
    matches!(
        error,
        ResourceError::Timeout { .. } | ResourceError::Transport { .. }
    ) || matches!(error, ResourceError::HttpStatus { status, .. } if (500..600).contains(status))
}

pub(crate) fn redact_message(message: &str) -> String {
    // Curl and HTTP stacks may include a complete request URL in diagnostics.
    // Keep the stable error useful while removing query strings and userinfo.
    message
        .split_whitespace()
        .map(redact_word)
        .collect::<Vec<_>>()
        .join(" ")
}

fn redact_word(word: &str) -> String {
    let without_suffix = word
        .split_once('#')
        .map_or(word, |(prefix, _)| prefix)
        .split_once('?')
        .map_or_else(
            || word.split_once('#').map_or(word, |(prefix, _)| prefix),
            |(prefix, _)| prefix,
        );
    let Some((scheme, remainder)) = without_suffix.split_once("://") else {
        return without_suffix.to_string();
    };
    let (authority, path) = remainder
        .split_once('/')
        .map_or((remainder, ""), |(authority, path)| (authority, path));
    let authority = authority
        .rsplit_once('@')
        .map_or(authority, |(_, host)| host);
    if path.is_empty() {
        format!("{scheme}://{authority}")
    } else {
        format!("{scheme}://{authority}/{path}")
    }
}

pub(crate) fn retry_count(max_retries: u32) -> u32 {
    max_retries.saturating_add(1)
}
