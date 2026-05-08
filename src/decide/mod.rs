//! Fragments, the mosaic of an entry, query coverage, identity and the decision that reports an
//! entry.

mod coverage;

pub use coverage::QueryCoverage;

use crate::align::AlignmentError;
use thiserror::Error;

/// What can go wrong while building fragments, a mosaic or a decision.
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum DecideError {
    #[error("the query is empty")]
    EmptyQuery,
    #[error("the query window does not fit the query")]
    InvalidQueryWindow,
    #[error("the query segments of a fragment are invalid")]
    InvalidQuerySegments,
    #[error("the query segments of a fragment do not match its CIGAR")]
    QuerySpanMismatch,
    #[error("a coordinate is too large")]
    LengthOverflow,
    #[error("the query coverage is invalid: {0}")]
    InvalidCoverage(&'static str),
    #[error("the decision rule is invalid: {0}")]
    InvalidRule(&'static str),
    #[error(transparent)]
    Alignment(#[from] AlignmentError),
}
