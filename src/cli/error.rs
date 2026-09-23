//! The errors of the commands, printed by `main` in one plain line.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::index::IndexError;
use crate::loss::LossError;
use crate::report::ReportError;
use crate::screen::ScreenError;
use crate::search::SearchError;
use crate::service::ServiceError;
use crate::storage::StorageError;
use std::path::PathBuf;
use thiserror::Error;

/// What can go wrong in a command.
#[derive(Debug, Error)]
pub enum CliError {
    /// An option or argument is wrong.
    #[error("{0}")]
    Usage(String),
    /// A file named on the command line cannot be read or created.
    #[error("cannot use {path}: {source}")]
    File {
        path: PathBuf,
        source: std::io::Error,
    },
    /// The query file is not valid FASTA or FASTQ.
    #[error("cannot read the query file {path}: {source}")]
    Query {
        path: PathBuf,
        source: needletail::errors::ParseError,
    },
    /// The query file holds no records.
    #[error("the query file {0} holds no records")]
    NoQueries(PathBuf),
    /// A shard cannot be opened.
    #[error("cannot open the shard {path}: {source}")]
    Shard { path: PathBuf, source: IndexError },
    /// Building or joining shards failed.
    #[error(transparent)]
    Index(#[from] IndexError),
    /// A search could not run.
    #[error("searching {query} failed: {source}")]
    Search { query: String, source: SearchError },
    /// The search settings are invalid.
    #[error(transparent)]
    Settings(#[from] SearchError),
    /// The loss budget or the miss model is invalid.
    #[error(transparent)]
    Loss(#[from] LossError),
    /// A record cannot be built.
    #[error(transparent)]
    Report(#[from] ReportError),
    /// A record cannot be written as JSON.
    #[error("cannot write a record: {0}")]
    Json(#[from] serde_json::Error),
    /// The output cannot be written.
    #[error("cannot write the output: {0}")]
    Output(#[from] std::io::Error),
    /// Searches ended before every possible entry was decided: a deadline passed, a stage
    /// failed, or a resource limit dropped an entry. PAF output then lacks those entries.
    #[error(
        "{0} searches did not complete (deadline, failure or resource limit); \
         the JSON summary lines tell which"
    )]
    Incomplete(usize),
    /// The screen settings are invalid or the shard cannot be screened.
    #[error(transparent)]
    Screen(#[from] ScreenError),
    /// Assemblies could not be screened; each has a failed record.
    #[error("{0} assemblies could not be screened")]
    FailedAssemblies(usize),
    /// The S3 store or the block cache cannot be set up.
    #[error(transparent)]
    Storage(#[from] StorageError),
    /// The runtime that waits for S3 reads cannot start.
    #[error("cannot start the runtime for S3 reads: {0}")]
    Runtime(std::io::Error),
    /// Serving or a remote search failed.
    #[error(transparent)]
    Service(#[from] ServiceError),
    /// The limit of open files cannot be raised to what a build needs.
    #[error(
        "the limit of open files is {limit}, but building a shard needs {needed}: raise the hard \
         limit (ulimit -Hn) and try again"
    )]
    OpenFiles { limit: u64, needed: u64 },
    /// The limit of open files cannot be read or raised.
    #[error("cannot raise the limit of open files: {0}")]
    FileLimit(std::io::Error),
    /// The search thread panicked.
    #[error("the search thread panicked")]
    Panic,
}
