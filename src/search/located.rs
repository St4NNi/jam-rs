//! What stage 2 hands to stage 3 for one candidate entry.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use std::ops::Range;

/// A candidate entry after stage 2: its bounds from the seed hits and where to align.
#[derive(Clone, Debug, PartialEq)]
pub struct LocatedEntry {
    /// The entry's index in the shard.
    pub entry: u32,
    /// The highest query coverage, in bases, that the entry can still reach.
    pub coverage_bound: u64,
    /// The containment estimate again, from the seed hits with seed locations only (frequent
    /// seeds included, false seed hits removed).
    pub containment_estimate: f64,
    /// The alignment windows, by contig in ascending order.
    pub windows: Vec<ContigWindows>,
}

/// The alignment windows of one contig: sorted, non-overlapping ranges in contig coordinates.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ContigWindows {
    pub contig: u32,
    pub ranges: Vec<Range<u64>>,
}
