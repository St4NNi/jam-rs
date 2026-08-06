//! Fragments: alignments placed on the query and on a contig of an entry.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::DecideError;
use crate::align::{Alignment, Interval};
use serde::Serialize;

/// An alignment placed on the query and on one contig of an entry. Over the end of a circular
/// query it has two query segments, up to the query end and from the query start.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Fragment {
    /// The index of the contig within its entry.
    pub contig_index: u32,
    /// The query bases the alignment spans, in alignment order.
    pub query_segments: Vec<Interval>,
    /// The alignment, with its query interval moved to start inside the query.
    pub alignment: Alignment,
}

impl Fragment {
    /// The fragment of `alignment` on contig `contig_index`. A circular alignment may start one
    /// or more query lengths on; it is moved back to start inside the query, and an alignment
    /// over the query end gets two query segments.
    pub fn place(
        contig_index: u32,
        mut alignment: Alignment,
        query_length: u64,
        circular: bool,
    ) -> Result<Self, DecideError> {
        if query_length == 0 {
            return Err(DecideError::EmptyQuery);
        }
        let interval = alignment.query_interval;
        let shift = interval.start / query_length * query_length;
        let moved = Interval::new(interval.start - shift, interval.end - shift)?;
        let fits = match circular {
            false => moved.end <= query_length,
            true => moved.len() <= query_length,
        };
        if !fits {
            return Err(DecideError::QueryOutside);
        }
        alignment.query_interval = moved;
        Ok(Self {
            contig_index,
            query_segments: moved.wrap(query_length).collect(),
            alignment,
        })
    }
}

#[cfg(test)]
pub(super) mod tests {
    use super::*;
    use crate::align::{Strand, parse_cigar};

    /// A forward fragment with `cigar`, scored 2 per match and -3 or -1 per other column, from
    /// position 0 of the query and of the contig.
    pub(in crate::decide) fn fragment(
        contig_index: u32,
        segments: &[(u64, u64)],
        cigar: &str,
    ) -> Fragment {
        let mut alignment = Alignment {
            score: 0,
            strand: Strand::Forward,
            query_interval: Interval::new(0, 0).unwrap(),
            contig_interval: Interval::new(0, 0).unwrap(),
            runs: parse_cigar(cigar).unwrap(),
        };
        let counts = alignment.counts();
        alignment.score = 2 * counts.matches as i32
            - (3 * counts.substitutions + counts.insertions + counts.deletions) as i32;
        alignment.query_interval = Interval::new(0, counts.query_bases()).unwrap();
        alignment.contig_interval = Interval::new(0, counts.contig_bases()).unwrap();
        let segments = segments.iter();
        let segments = segments.map(|&(start, end)| Interval::new(start, end).unwrap());
        Fragment {
            contig_index,
            query_segments: segments.collect(),
            alignment,
        }
    }

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
    fn origin_splits_circular() {
        let over = Fragment::place(3, alignment((90, 120)), 100, true).unwrap();
        assert_eq!(over.query_segments, [interval(90, 100), interval(0, 20)]);
        assert_eq!(over.contig_index, 3);
        let moved = Fragment::place(0, alignment((110, 150)), 100, true).unwrap();
        assert_eq!(moved.alignment.query_interval, interval(10, 50));
        assert_eq!(moved.query_segments, [interval(10, 50)]);
        let whole = Fragment::place(0, alignment((30, 130)), 100, true).unwrap();
        assert_eq!(whole.query_segments, [interval(30, 100), interval(0, 30)]);
    }

    #[test]
    fn outside_query_rejected() {
        let linear = Fragment::place(0, alignment((90, 120)), 100, false);
        assert_eq!(linear, Err(DecideError::QueryOutside));
        let longer = Fragment::place(0, alignment((10, 111)), 100, true);
        assert_eq!(longer, Err(DecideError::QueryOutside));
        let empty = Fragment::place(0, alignment((0, 5)), 0, true);
        assert_eq!(empty, Err(DecideError::EmptyQuery));
        let inside = Fragment::place(0, alignment((5, 100)), 100, false).unwrap();
        assert_eq!(inside.query_segments, [interval(5, 100)]);
    }
}
