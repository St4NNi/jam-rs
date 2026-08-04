//! The mosaic: the fragments of one entry that cover the query, each query base once.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{DecideError, Fragment};
use crate::align::{EditCounts, EditOperation, Interval};
use serde::Serialize;
use std::cmp::Ordering;

/// A fragment of the mosaic with the query bases it alone adds.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct SelectedFragment {
    pub new_bases: u64,
    pub fragment: Fragment,
}

/// The fragments chosen to cover the query, best first, and the query bases they cover.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Mosaic {
    pub query_length: u64,
    /// Query bases that are aligned to a contig base (matches and substitutions) in some fragment.
    pub covered_bases: u64,
    pub covered_intervals: Vec<Interval>,
    pub gaps: Vec<Interval>,
    /// Fragments that cover query bases no better fragment covers.
    pub primary: Vec<SelectedFragment>,
    /// Fragments that add nothing new, such as other copies of the same region.
    pub alternatives: Vec<Fragment>,
}

/// A fragment with the query bases it pairs and its base counts, kept for sorting.
struct Prepared {
    fragment: Fragment,
    support: Vec<Interval>,
    support_bases: u64,
    counts: EditCounts,
}

/// Builds the mosaic of `fragments` for a query of `query_length` bases.
///
/// Fragments are taken by score, best first. A fragment is primary when it covers a query base
/// that no earlier fragment covers. The result does not depend on the input order, and identical
/// fragments count once.
pub(super) fn build_mosaic(
    query_length: u64,
    fragments: &[Fragment],
) -> Result<Mosaic, DecideError> {
    if query_length == 0 {
        return Err(DecideError::EmptyQuery);
    }
    let mut prepared = fragments
        .iter()
        .map(|fragment| {
            let support = covered_intervals(query_length, fragment)?;
            Ok(Prepared {
                fragment: fragment.clone(),
                support_bases: Interval::total_length(&support),
                support,
                counts: fragment.alignment.counts(),
            })
        })
        .collect::<Result<Vec<_>, DecideError>>()?;
    prepared.sort_by(compare_fragments);
    prepared.dedup_by(|left, right| left.fragment == right.fragment);

    let mut covered = Vec::new();
    let mut primary = Vec::new();
    let mut alternatives = Vec::new();
    for item in prepared {
        let new_bases = item.support_bases - overlap(&covered, &item.support);
        if new_bases == 0 {
            alternatives.push(item.fragment);
        } else {
            covered.extend_from_slice(&item.support);
            covered = Interval::union(covered);
            primary.push(SelectedFragment {
                new_bases,
                fragment: item.fragment,
            });
        }
    }
    Ok(Mosaic {
        query_length,
        covered_bases: Interval::total_length(&covered),
        gaps: gaps_between(&covered, query_length),
        covered_intervals: covered,
        primary,
        alternatives,
    })
}

fn compare_fragments(left: &Prepared, right: &Prepared) -> Ordering {
    let (a, b) = (&left.fragment, &right.fragment);
    let edits = |counts: EditCounts| counts.substitutions + counts.insertions + counts.deletions;
    b.alignment
        .score
        .cmp(&a.alignment.score)
        .then_with(|| right.support_bases.cmp(&left.support_bases))
        .then_with(|| right.counts.matches.cmp(&left.counts.matches))
        .then_with(|| edits(left.counts).cmp(&edits(right.counts)))
        .then_with(|| a.contig_index.cmp(&b.contig_index))
        .then_with(|| a.alignment.strand.cmp(&b.alignment.strand))
        .then_with(|| {
            a.alignment
                .contig_interval
                .cmp(&b.alignment.contig_interval)
        })
        .then_with(|| a.alignment.query_interval.cmp(&b.alignment.query_interval))
        .then_with(|| a.query_segments.cmp(&b.query_segments))
        .then_with(|| a.alignment.runs.cmp(&b.alignment.runs))
}

/// The query bases that a fragment pairs with a contig base: matches and substitutions. Bases
/// only in the query (`I`) are skipped over and bases only in the contig (`D`) use no query base.
fn covered_intervals(query_length: u64, fragment: &Fragment) -> Result<Vec<Interval>, DecideError> {
    fragment.alignment.validate()?;
    let segments = &fragment.query_segments;
    if segments.is_empty()
        || segments
            .iter()
            .any(|segment| segment.is_empty() || segment.end > query_length)
    {
        return Err(DecideError::InvalidQuerySegments);
    }
    let mut ordered = segments.clone();
    ordered.sort_unstable();
    if ordered.windows(2).any(|pair| pair[0].end > pair[1].start) {
        return Err(DecideError::InvalidQuerySegments);
    }
    let segment_bases = Interval::total_length(segments);
    if segment_bases != fragment.alignment.query_interval.len() {
        return Err(DecideError::QuerySpanMismatch);
    }

    let mut offset = 0u64;
    let mut support = Vec::new();
    for run in &fragment.alignment.runs {
        let length = u64::from(run.length);
        match run.operation {
            EditOperation::Match | EditOperation::Substitution => {
                project_segments(segments, offset, length, &mut support)?;
                offset += length;
            }
            EditOperation::Insertion => offset += length,
            EditOperation::Deletion => {}
        }
    }
    if offset != segment_bases {
        return Err(DecideError::QuerySpanMismatch);
    }
    Ok(Interval::union(support))
}

/// Maps `length` fragment query positions from `offset` on to intervals of the query segments.
fn project_segments(
    segments: &[Interval],
    mut offset: u64,
    mut length: u64,
    output: &mut Vec<Interval>,
) -> Result<(), DecideError> {
    for segment in segments {
        if offset >= segment.len() {
            offset -= segment.len();
            continue;
        }
        let count = length.min(segment.len() - offset);
        let start = segment.start + offset;
        output.push(Interval::new(start, start + count)?);
        length -= count;
        offset = 0;
        if length == 0 {
            return Ok(());
        }
    }
    Err(DecideError::QuerySpanMismatch)
}

/// The bases of the sorted, disjoint `support` that `covered`, also sorted and disjoint, holds.
fn overlap(covered: &[Interval], support: &[Interval]) -> u64 {
    let mut bases = 0;
    for piece in support {
        let first = covered.partition_point(|interval| interval.end <= piece.start);
        for interval in covered[first..]
            .iter()
            .take_while(|interval| interval.start < piece.end)
        {
            bases += interval.end.min(piece.end) - interval.start.max(piece.start);
        }
    }
    bases
}

fn gaps_between(covered: &[Interval], query_length: u64) -> Vec<Interval> {
    let mut gaps = Vec::new();
    let mut start = 0;
    for interval in covered {
        if start < interval.start {
            gaps.push(Interval {
                start,
                end: interval.start,
            });
        }
        start = interval.end;
    }
    if start < query_length {
        gaps.push(Interval {
            start,
            end: query_length,
        });
    }
    gaps
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{Alignment, Strand, parse_cigar};

    /// A fragment on a forward strand; the contig interval just has the right length.
    fn fragment(contig_index: u32, score: i32, segments: &[(u64, u64)], cigar: &str) -> Fragment {
        let runs = parse_cigar(cigar).unwrap();
        let mut alignment = Alignment {
            score,
            strand: Strand::Forward,
            query_interval: Interval::new(0, 0).unwrap(),
            contig_interval: Interval::new(0, 0).unwrap(),
            runs,
        };
        let counts = alignment.counts();
        let query_bases = counts.matches + counts.substitutions + counts.insertions;
        let contig_bases = counts.matches + counts.substitutions + counts.deletions;
        alignment.query_interval = Interval::new(0, query_bases).unwrap();
        alignment.contig_interval = Interval::new(0, contig_bases).unwrap();
        Fragment {
            contig_index,
            query_segments: segments
                .iter()
                .map(|&(start, end)| Interval::new(start, end).unwrap())
                .collect(),
            alignment,
        }
    }

    fn equal_fragment(contig_index: u32, score: i32, start: u64, end: u64) -> Fragment {
        fragment(
            contig_index,
            score,
            &[(start, end)],
            &format!("{}=", end - start),
        )
    }

    fn interval(start: u64, end: u64) -> Interval {
        Interval::new(start, end).unwrap()
    }

    #[test]
    fn bases_only_in_the_query_are_gaps() {
        // 2= 1I 2=: the middle query base has no contig base.
        let mosaic = build_mosaic(5, &[fragment(0, 7, &[(0, 5)], "2=1I2=")]).unwrap();
        assert_eq!(mosaic.covered_bases, 4);
        assert_eq!(mosaic.gaps, [interval(2, 3)]);
    }

    #[test]
    fn bases_only_in_the_contig_use_no_query_base() {
        let mosaic = build_mosaic(4, &[fragment(0, 7, &[(0, 4)], "2=3D2=")]).unwrap();
        assert_eq!(mosaic.covered_bases, 4);
        assert!(mosaic.gaps.is_empty());
    }

    #[test]
    fn later_fragments_count_only_new_bases() {
        let input = [
            equal_fragment(2, 7, 1, 5),
            equal_fragment(1, 8, 4, 10),
            equal_fragment(0, 10, 0, 6),
        ];
        let mosaic = build_mosaic(10, &input).unwrap();
        assert_eq!(mosaic.covered_bases, 10);
        assert_eq!(mosaic.primary.len(), 2);
        assert_eq!(mosaic.primary[0].fragment.contig_index, 0);
        assert_eq!(mosaic.primary[1].new_bases, 4);
        assert_eq!(mosaic.alternatives.len(), 1);
        let mut reversed = input.to_vec();
        reversed.reverse();
        assert_eq!(build_mosaic(10, &reversed).unwrap(), mosaic);
    }

    #[test]
    fn duplicates_count_once() {
        let mut a = equal_fragment(0, 100, 0, 30);
        let mut b = equal_fragment(1, 100, 30, 60);
        b.alignment.strand = Strand::Reverse;
        let combined = build_mosaic(101, &[a.clone(), b.clone(), a.clone(), b]).unwrap();
        assert_eq!(combined.covered_bases, 60);
        assert_eq!(combined.primary.len(), 2);
        assert!(combined.alternatives.is_empty());
        assert_eq!(combined.query_length, 101);
        a.contig_index = 5;
        let one = build_mosaic(101, &[a]).unwrap();
        assert_eq!(one.covered_bases, 30);
    }

    #[test]
    fn a_fragment_over_the_query_end_covers_both_ends() {
        let wrapped = fragment(0, 8, &[(8, 10), (0, 2)], "4=");
        let mosaic = build_mosaic(10, &[wrapped]).unwrap();
        assert_eq!(mosaic.covered_intervals, [interval(0, 2), interval(8, 10)]);
        assert_eq!(mosaic.gaps, [interval(2, 8)]);
    }

    #[test]
    fn query_only_bases_next_to_the_query_end_stay_uncovered() {
        let wrapped = fragment(0, 8, &[(8, 10), (0, 3)], "1=1I1=1I1=");
        let mosaic = build_mosaic(10, &[wrapped]).unwrap();
        assert_eq!(
            mosaic.covered_intervals,
            [interval(0, 1), interval(2, 3), interval(8, 9)]
        );
    }

    #[test]
    fn bad_fragments_are_errors() {
        assert_eq!(build_mosaic(0, &[]), Err(DecideError::EmptyQuery));
        let outside = fragment(0, 1, &[(5, 12)], "7=");
        assert_eq!(
            build_mosaic(10, &[outside]),
            Err(DecideError::InvalidQuerySegments)
        );
        let overlap = fragment(0, 1, &[(0, 4), (3, 6)], "7=");
        assert_eq!(
            build_mosaic(10, &[overlap]),
            Err(DecideError::InvalidQuerySegments)
        );
        let short = fragment(0, 1, &[(0, 4)], "5=");
        assert_eq!(
            build_mosaic(10, &[short]),
            Err(DecideError::QuerySpanMismatch)
        );
    }
}
