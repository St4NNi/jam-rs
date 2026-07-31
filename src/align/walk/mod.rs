//! Walks: ungapped copies of dense and low-complexity k-mers, their alignments and their joins.

mod dense;
mod gapped;
mod join;
#[cfg(test)]
mod tests;
mod ungapped;

pub use dense::{DenseWalks, Walk, dense_walks, long_enough};
pub use gapped::{DenseOptimum, align_dense_walk};
pub use join::{JoinContext, JoinLimits, JoinScratch};
pub use ungapped::{gap_continuation, holds_walk, walk_alignment};
