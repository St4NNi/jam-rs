//! Counted segments: the stretches of a fragment whose identity reaches the fragment limit.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::mosaic::project_segments;
use super::{DecideError, Fragment};
use crate::align::{Alignment, EditCounts, EditOperation, EditRun, Interval, Strand};

/// The limit is turned into whole numbers with this many steps, so that a stretch is tested
/// exactly: matches * (SCALE - step) >= (columns - matches) * step.
const SCALE: i128 = 1_000_000_000;

/// The limit in steps of 1 / `SCALE`.
fn limit_steps(limit: f64) -> i128 {
    (limit.clamp(0.0, 1.0) * SCALE as f64).round() as i128
}

/// Whether `counts` reach the identity `limit`, taken to nine decimals and compared exactly, as
/// the segments are.
pub(super) fn reaches(counts: EditCounts, limit: f64) -> bool {
    i128::from(counts.matches) * SCALE >= limit_steps(limit) * i128::from(counts.columns())
}

/// Totals of a fragment before one of its runs: score at the limit, columns, and the query and
/// contig bases spanned.
#[derive(Clone, Copy, Default)]
struct Totals {
    score: i128,
    columns: u64,
    query: u64,
    contig: u64,
}

/// The run totals of a fragment, one entry before each run and one after the last.
pub(super) struct RunOffsets(Vec<Totals>);

impl RunOffsets {
    /// Sums the runs once, with a match scoring 1 - `limit` per column and any other column
    /// -`limit`, as whole numbers.
    pub(super) fn new(runs: &[EditRun], limit: f64) -> Self {
        let step = limit_steps(limit);
        let mut totals = Vec::with_capacity(runs.len() + 1);
        let mut total = Totals::default();
        totals.push(total);
        for run in runs {
            let length = u64::from(run.length);
            let (query, contig) = match run.operation {
                EditOperation::Match | EditOperation::Substitution => (true, true),
                EditOperation::Insertion => (true, false),
                EditOperation::Deletion => (false, true),
            };
            total.score += match run.operation {
                EditOperation::Match => (SCALE - step) * i128::from(run.length),
                _ => -step * i128::from(run.length),
            };
            total.columns += length;
            total.query += if query { length } else { 0 };
            total.contig += if contig { length } else { 0 };
            totals.push(total);
        }
        Self(totals)
    }
}

/// A candidate of the scan: runs `first..last`, with the score total before and after them, and
/// the earlier candidate that the search for a join jumps to.
struct Candidate {
    first: usize,
    last: usize,
    before: i128,
    after: i128,
    /// The rightmost earlier candidate whose "before" is at most this one's.
    left: Option<usize>,
}

/// The maximal segments of `runs` with identity at least the limit of `offsets`, as ranges of
/// runs (Ruzzo and Tompa 1999; ties in score are joined). The jumps along `left` make the scan
/// linear in the number of runs.
pub(super) fn maximal_segments(runs: &[EditRun], offsets: &RunOffsets) -> Vec<(usize, usize)> {
    let totals = &offsets.0;
    let mut list: Vec<Candidate> = Vec::new();
    for (index, run) in runs.iter().enumerate() {
        if run.operation != EditOperation::Match {
            continue;
        }
        let mut item = Candidate {
            first: index,
            last: index + 1,
            before: totals[index].score,
            after: totals[index + 1].score,
            left: None,
        };
        let mut at = list.len().checked_sub(1);
        loop {
            // Candidates skipped by a jump all have a higher "before" than the one jumped from.
            while let Some(other) = at.filter(|&other| list[other].before > item.before) {
                at = list[other].left;
            }
            match at {
                Some(other) if list[other].after <= item.after => {
                    item.first = list[other].first;
                    item.before = list[other].before;
                    at = list[other].left;
                    list.truncate(other);
                }
                _ => break,
            }
        }
        item.left = at;
        list.push(item);
    }
    list.iter().map(|item| (item.first, item.last)).collect()
}

/// The part of `fragment` made of runs `first..last`, with intervals and score cut to it.
pub(super) fn slice_fragment(
    fragment: &Fragment,
    offsets: &RunOffsets,
    first: usize,
    last: usize,
) -> Result<Fragment, DecideError> {
    let alignment = &fragment.alignment;
    let (head, tail, whole) = (
        offsets.0[first],
        offsets.0[last],
        offsets.0[alignment.runs.len()],
    );
    let (query_start, query_length) = (head.query, tail.query - head.query);
    let (contig_start, contig_length) = (head.contig, tail.contig - head.contig);
    let query = alignment.query_interval.start + query_start;
    let contig = alignment.contig_interval;
    let contig_interval = match alignment.strand {
        Strand::Forward => Interval::new(
            contig.start + contig_start,
            contig.start + contig_start + contig_length,
        )?,
        Strand::Reverse => Interval::new(
            contig.end - contig_start - contig_length,
            contig.end - contig_start,
        )?,
    };
    let mut query_segments = Vec::new();
    project_segments(
        &fragment.query_segments,
        query_start,
        query_length,
        &mut query_segments,
    )?;
    let part = tail.columns - head.columns;
    let score = i128::from(alignment.score) * i128::from(part) / i128::from(whole.columns.max(1));
    Ok(Fragment {
        contig_index: fragment.contig_index,
        query_segments,
        alignment: Alignment {
            score: score as i32,
            strand: alignment.strand,
            query_interval: Interval::new(query, query + query_length)?,
            contig_interval,
            runs: alignment.runs[first..last].to_vec(),
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::parse_cigar;
    use crate::decide::fragment::tests::fragment;
    use rand::RngExt;

    fn segments(cigar: &str, limit: f64) -> Vec<(usize, usize)> {
        let runs = parse_cigar(cigar).unwrap();
        maximal_segments(&runs, &RunOffsets::new(&runs, limit))
    }

    #[test]
    fn plain_runs() {
        assert_eq!(segments("50=", 0.9), [(0, 1)]);
        assert_eq!(segments("50X", 0.9), []);
        // The mismatches cost more than the first stretch brings, so it stays alone.
        assert_eq!(segments("10=5X100=", 0.9), [(0, 1), (2, 3)]);
    }

    #[test]
    fn dense_stretch_joins() {
        // 90 matches and 5 mismatches over 100 columns: 100 / 106 is above 0.9 and joins.
        assert_eq!(segments("50=3X50=", 0.9), [(0, 3)]);
        // Too many mismatches in between keep the two stretches apart.
        assert_eq!(segments("20=30X20=", 0.9), [(0, 1), (2, 3)]);
        // The last stretch joins the first by a jump over the middle one.
        assert_eq!(segments("50=2X10=2X100=", 0.9), [(0, 5)]);
    }

    #[test]
    fn break_even_counts() {
        // Each 9 matches and 1 mismatch are exactly 0.9, so the whole stretch counts.
        assert_eq!(segments("9=1X9=1X9=", 0.9), [(0, 5)]);
        // A limit of 1 takes every exact run and nothing else.
        assert_eq!(segments("9=1X9=", 1.0), [(0, 1), (2, 3)]);
    }

    #[test]
    fn many_runs_linear() {
        // 116,279 blocks of 40= 3X: each block scores below 0, so no two join, and each new
        // candidate starts below all earlier ones.
        let runs = parse_cigar(&"40=3X".repeat(116_279)).unwrap();
        let found = maximal_segments(&runs, &RunOffsets::new(&runs, 0.95));
        assert_eq!(found.len(), 116_279);
        assert!(found.iter().all(|&(first, last)| last == first + 1));
    }

    /// The scan without jumps: each new candidate searches all earlier ones.
    fn plain_scan(runs: &[EditRun], offsets: &RunOffsets) -> Vec<(usize, usize)> {
        let mut list: Vec<(usize, usize, i128, i128)> = Vec::new();
        for (index, run) in runs.iter().enumerate() {
            if run.operation != EditOperation::Match {
                continue;
            }
            let (before, after) = (offsets.0[index].score, offsets.0[index + 1].score);
            let mut item = (index, index + 1, before, after);
            while let Some(at) = list.iter().rposition(|other| other.2 <= item.2) {
                if list[at].3 > item.3 {
                    break;
                }
                (item.0, item.2) = (list[at].0, list[at].2);
                list.truncate(at);
            }
            list.push(item);
        }
        list.iter().map(|item| (item.0, item.1)).collect()
    }

    #[test]
    fn jumps_match_plain() {
        let mut rng = crate::test_support::rng(11);
        for _ in 0..2000 {
            let mut cigar = String::new();
            for _ in 0..rng.random_range(1..30) {
                let other = ["X", "I", "D"][rng.random_range(0..3)];
                let (matches, edits) = (rng.random_range(1..40), rng.random_range(1..6));
                cigar += &format!("{matches}={edits}{other}");
            }
            let runs = parse_cigar(&cigar).unwrap();
            let limit = [0.0, 0.8, 0.9, 0.95, 1.0][rng.random_range(0..5)];
            let offsets = RunOffsets::new(&runs, limit);
            assert_eq!(
                maximal_segments(&runs, &offsets),
                plain_scan(&runs, &offsets)
            );
        }
    }

    #[test]
    fn slice_cuts_intervals() {
        let whole = fragment(2, &[(0, 12)], "4=2X6=");
        let offsets = RunOffsets::new(&whole.alignment.runs, 0.9);
        let part = slice_fragment(&whole, &offsets, 2, 3).unwrap();
        assert_eq!(part.query_segments, [Interval::new(6, 12).unwrap()]);
        assert_eq!(part.alignment.query_interval, Interval::new(6, 12).unwrap());
        assert_eq!(
            part.alignment.contig_interval,
            Interval::new(6, 12).unwrap()
        );
        assert_eq!(part.alignment.runs, parse_cigar("6=").unwrap());
        assert!(part.alignment.validate().is_ok());
        let mut reverse = whole.clone();
        reverse.alignment.strand = Strand::Reverse;
        let part = slice_fragment(&reverse, &offsets, 0, 1).unwrap();
        assert_eq!(
            part.alignment.contig_interval,
            Interval::new(8, 12).unwrap()
        );
        assert!(part.alignment.validate().is_ok());
    }
}
