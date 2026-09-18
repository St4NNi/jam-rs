//! Counted segments: the stretches of a fragment whose identity reaches the fragment limit.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::mosaic::project_segments;
use super::{DecideError, Fragment};
use crate::align::{Alignment, EditOperation, EditRun, Interval, Strand};

/// The limit is turned into whole numbers with this many steps, so that a segment is tested
/// exactly: matches * (SCALE - step) >= (columns - matches) * step.
const SCALE: i128 = 1_000_000_000;

/// A candidate of the scan: runs `first..last`, with the score total before and after them.
struct Candidate {
    first: usize,
    last: usize,
    before: i128,
    after: i128,
}

/// The maximal segments of `runs` with identity at least `limit`, as ranges of runs (Ruzzo and
/// Tompa 1999; ties in score are joined). A match scores 1 - limit per column and any other
/// column -limit, so a segment scores 0 or more when matches / columns >= limit.
pub(super) fn maximal_segments(runs: &[EditRun], limit: f64) -> Vec<(usize, usize)> {
    let step = (limit.clamp(0.0, 1.0) * SCALE as f64).round() as i128;
    if step >= SCALE {
        // Only exact stretches reach a limit of 1.
        let matches = runs.iter().enumerate();
        let matches = matches.filter(|(_, run)| run.operation == EditOperation::Match);
        return matches.map(|(index, _)| (index, index + 1)).collect();
    }
    let mut list: Vec<Candidate> = Vec::new();
    let mut total = 0i128;
    for (index, run) in runs.iter().enumerate() {
        let length = i128::from(run.length);
        if run.operation != EditOperation::Match {
            total -= step * length;
            continue;
        }
        let mut item = Candidate {
            first: index,
            last: index + 1,
            before: total,
            after: total + (SCALE - step) * length,
        };
        total = item.after;
        loop {
            let left = list.iter().rposition(|other| other.before <= item.before);
            match left {
                Some(at) if list[at].after <= item.after => {
                    item.first = list[at].first;
                    item.before = list[at].before;
                    list.truncate(at);
                }
                _ => break,
            }
        }
        list.push(item);
    }
    list.iter().map(|item| (item.first, item.last)).collect()
}

/// The part of `fragment` made of runs `first..last`, with intervals and score cut to it.
pub(super) fn slice_fragment(
    fragment: &Fragment,
    first: usize,
    last: usize,
) -> Result<Fragment, DecideError> {
    let alignment = &fragment.alignment;
    let bases = |runs: &[EditRun], query: bool| -> u64 {
        let mut total = 0;
        for run in runs {
            let length = u64::from(run.length);
            total += match (run.operation, query) {
                (EditOperation::Match | EditOperation::Substitution, _) => length,
                (EditOperation::Insertion, true) | (EditOperation::Deletion, false) => length,
                _ => 0,
            };
        }
        total
    };
    let (head, body) = (&alignment.runs[..first], &alignment.runs[first..last]);
    let (query_start, query_length) = (bases(head, true), bases(body, true));
    let (contig_start, contig_length) = (bases(head, false), bases(body, false));
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
    let columns: u64 = alignment.runs.iter().map(|run| u64::from(run.length)).sum();
    let part: u64 = body.iter().map(|run| u64::from(run.length)).sum();
    let score =
        (i128::from(alignment.score) * i128::from(part) / i128::from(columns.max(1))) as i32;
    Ok(Fragment {
        contig_index: fragment.contig_index,
        query_segments,
        alignment: Alignment {
            score,
            strand: alignment.strand,
            query_interval: Interval::new(query, query + query_length)?,
            contig_interval,
            runs: body.to_vec(),
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::parse_cigar;
    use crate::decide::fragment::tests::fragment;

    fn segments(cigar: &str, limit: f64) -> Vec<(usize, usize)> {
        maximal_segments(&parse_cigar(cigar).unwrap(), limit)
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
    }

    #[test]
    fn break_even_counts() {
        // Each 9 matches and 1 mismatch are exactly 0.9, so the whole stretch counts.
        let runs = parse_cigar("9=1X9=1X9=").unwrap();
        let found = maximal_segments(&runs, 0.9);
        assert_eq!(found, [(0, 5)]);
        // A limit of 1 takes every exact run and nothing else.
        assert_eq!(segments("9=1X9=", 1.0), [(0, 1), (2, 3)]);
    }

    #[test]
    fn slice_cuts_intervals() {
        let whole = fragment(2, &[(0, 12)], "4=2X6=");
        let part = slice_fragment(&whole, 2, 3).unwrap();
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
        let part = slice_fragment(&reverse, 0, 1).unwrap();
        assert_eq!(
            part.alignment.contig_interval,
            Interval::new(8, 12).unwrap()
        );
        assert!(part.alignment.validate().is_ok());
    }
}
