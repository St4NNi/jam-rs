//! The errors of a search.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::align_stage::AlignStageError;
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
    /// Stage 3 could not decide an entry, or its settings are invalid.
    #[error(transparent)]
    Align(#[from] AlignStageError),
    /// The search passed its deadline or was asked to stop.
    #[error(transparent)]
    Stopped(#[from] Stopped),
}
