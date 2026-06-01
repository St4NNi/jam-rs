//! The errors of a search.

use super::control::Stopped;
use crate::index::IndexError;
use crate::loss::LossError;
use thiserror::Error;

/// Why a search or one of its stages failed.
#[derive(Debug, Error)]
pub enum SearchError {
    /// The query cannot be searched.
    #[error("invalid query: {0}")]
    InvalidQuery(&'static str),
    /// A setting of the search is out of range.
    #[error("invalid search setting: {0}")]
    InvalidConfig(&'static str),
    /// Reading the index failed.
    #[error(transparent)]
    Index(#[from] IndexError),
    /// The miss model does not fit the shard.
    #[error(transparent)]
    Loss(#[from] LossError),
    /// The search passed its deadline or was asked to stop.
    #[error(transparent)]
    Stopped(#[from] Stopped),
}
