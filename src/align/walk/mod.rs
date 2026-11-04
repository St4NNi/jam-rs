//! Walks: ungapped copies of dense and low-complexity k-mers, their alignments and joins.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

mod dense;
mod gapped;
mod join;
#[cfg(test)]
mod tests;
mod ungapped;

pub use dense::{DenseWalks, Walk, dense_walks, long_enough};
pub use gapped::{DenseResults, DenseWalk, KernelInput, align_dense_walk};
#[cfg(test)]
pub(crate) use gapped::{SAME_FOUND, SAME_OFF};
pub use join::{JoinContext, JoinLimits, JoinScratch};
pub use ungapped::{gap_continuation, holds_walk, walk_alignment};
