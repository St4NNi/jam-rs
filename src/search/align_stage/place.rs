//! Alignments placed as fragments, and the fragments one entry collects within its limit.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{AlignStageError, Limit, PreparedQuery};
use crate::align::{EditRun, Interval};
use crate::decide::Fragment;
use std::mem::size_of;

/// Whether the query has bases before and after an alignment with this query interval. A
/// circular query has room on both sides until an alignment spans all of it.
pub(super) fn query_room(interval: Interval, query: &PreparedQuery) -> (bool, bool) {
    if query.circular() {
        let room = interval.len() < query.length();
        (room, room)
    } else {
        (interval.start > 0, interval.end < query.length())
    }
}

/// The fragments of one entry and the bytes they hold, within a limit.
pub(super) struct Collected {
    list: Vec<Fragment>,
    bytes: usize,
    limit: usize,
}

/// A count of fragments and their bytes, to go back to when a round is aligned again.
#[derive(Clone, Copy)]
pub(super) struct Mark {
    length: usize,
    bytes: usize,
}

fn fragment_bytes(fragment: &Fragment) -> usize {
    size_of::<Fragment>()
        + fragment.query_segments.capacity() * size_of::<Interval>()
        + fragment.alignment.runs.capacity() * size_of::<EditRun>()
}

impl Collected {
    pub(super) fn new(limit: usize) -> Self {
        Self {
            list: Vec::new(),
            bytes: 0,
            limit,
        }
    }

    /// Adds a fragment, or fails when the entry's fragments would pass the limit.
    pub(super) fn push(&mut self, fragment: Fragment) -> Result<(), AlignStageError> {
        let bytes = self.bytes + fragment_bytes(&fragment);
        if bytes > self.limit {
            return Err(AlignStageError::TooLarge {
                what: "the fragments of one entry",
                needed: bytes as u64,
                limit: self.limit as u64,
                kind: Limit::FragmentBytes,
            });
        }
        self.bytes = bytes;
        self.list.push(fragment);
        Ok(())
    }

    pub(super) fn mark(&self) -> Mark {
        Mark {
            length: self.list.len(),
            bytes: self.bytes,
        }
    }

    /// Drops the fragments added after `mark`.
    pub(super) fn reset(&mut self, mark: Mark) {
        self.list.truncate(mark.length);
        self.bytes = mark.bytes;
    }

    pub(super) fn into_list(self) -> Vec<Fragment> {
        self.list
    }
}

#[cfg(test)]
mod tests {
    use super::super::AlignStageConfig;
    use super::*;
    use crate::align::{Alignment, Strand, parse_cigar};

    fn alignment(query: (u64, u64)) -> Alignment {
        let length = query.1 - query.0;
        Alignment {
            score: 1,
            strand: Strand::Forward,
            query_interval: Interval::new(query.0, query.1).unwrap(),
            contig_interval: Interval::new(0, length).unwrap(),
            runs: parse_cigar(&format!("{length}=")).unwrap(),
        }
    }

    fn interval(start: u64, end: u64) -> Interval {
        Interval::new(start, end).unwrap()
    }

    #[test]
    fn query_room_sides() {
        let bases = b"ACGTTGCAAC".repeat(10);
        let config = AlignStageConfig::default();
        let circular = PreparedQuery::new(&bases, true, &config).unwrap();
        let linear = PreparedQuery::new(&bases, false, &config).unwrap();
        assert_eq!(query_room(interval(0, 100), &linear), (false, false));
        assert_eq!(query_room(interval(5, 50), &linear), (true, true));
        assert_eq!(query_room(interval(0, 50), &circular), (true, true));
        assert_eq!(query_room(interval(50, 150), &circular), (false, false));
    }

    #[test]
    fn fragment_limit_rounds() {
        let bases = b"ACGTTGCAAC".repeat(10);
        let query = PreparedQuery::new(&bases, false, &AlignStageConfig::default()).unwrap();
        let fragment = Fragment::place(0, alignment((0, 50)), query.length(), false).unwrap();
        let one = fragment_bytes(&fragment);
        let mut collected = Collected::new(2 * one);
        collected.push(fragment.clone()).unwrap();
        let mark = collected.mark();
        collected.push(fragment.clone()).unwrap();
        assert!(collected.push(fragment.clone()).is_err());
        collected.reset(mark);
        collected.push(fragment).unwrap();
        assert_eq!(collected.into_list().len(), 2);
    }
}
