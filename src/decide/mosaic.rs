//! The mosaic: the fragments of one entry that cover the query, each query base once.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{DecideError, Fragment};
use crate::align::{EditCounts, EditOperation, Interval};
use serde::Serialize;
use std::cmp::Ordering;
use std::collections::BTreeMap;

/// A fragment of the mosaic with the query bases it alone adds.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct SelectedFragment {
    /// Covered query bases that no better fragment covers.
    pub new_bases: u64,
    /// The fragment.
    pub fragment: Fragment,
}

/// The fragments chosen to cover the query, best first, and the query bases they cover.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Mosaic {
    /// The number of query bases.
    pub query_length: u64,
    /// Query bases that are aligned to a contig base (matches and substitutions) in some fragment.
    pub covered_bases: u64,
    /// The covered query bases as sorted, disjoint intervals.
    pub covered_intervals: Vec<Interval>,
    /// The query bases that no fragment covers, sorted.
    pub gaps: Vec<Interval>,
    /// Fragments that cover query bases no better fragment covers.
    pub primary: Vec<SelectedFragment>,
    /// Fragments that add nothing new, such as other copies of the same region.
    pub alternatives: Vec<Fragment>,
    /// Base counts of the primary fragments, each query base counted once: in the first primary
    /// fragment that spans it.
    pub counts: EditCounts,
}

/// A fragment with the query bases it pairs and its base counts, kept for sorting.
struct Prepared {
    fragment: Fragment,
    /// The query bases the fragment pairs with contig bases, sorted and disjoint.
    support: Vec<Interval>,
    support_bases: u64,
    counts: EditCounts,
}

/// A primary fragment with the covered query bases it adds and its counts over them.
struct Added {
    fragment: Fragment,
    bases: Vec<Interval>,
    counts: EditCounts,
}

/// Builds the mosaic of `fragments`, highest identity first: a fragment counts over the query
/// bases it adds outside `low_complexity`. When all are below `min_identity`, those whose added
/// bases are below it are left out, so the rest reach it.
pub(super) fn build_mosaic(
    query_length: u64,
    fragments: &[Fragment],
    low_complexity: &[Interval],
    min_identity: f64,
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

    // Low-complexity bases start as taken, so no fragment adds them or counts over them.
    let mut covered = IntervalSet::default();
    low_complexity
        .iter()
        .for_each(|&piece| covered.insert(piece));
    // The query bases that earlier primary fragments span, insertions too, for the counts.
    let mut spanned = IntervalSet(covered.0.clone());
    let mut added = Vec::new();
    let mut alternatives = Vec::new();
    for item in prepared {
        let mut bases = Vec::new();
        for &piece in &item.support {
            covered.outside(piece, &mut bases);
        }
        if bases.is_empty() {
            alternatives.push(item.fragment);
            continue;
        }
        let counts = new_counts(&item.fragment, &spanned)?;
        item.fragment
            .query_segments
            .iter()
            .for_each(|&piece| spanned.insert(piece));
        bases.iter().for_each(|&piece| covered.insert(piece));
        added.push(Added {
            fragment: item.fragment,
            bases,
            counts,
        });
    }
    let total = added.iter().fold(EditCounts::default(), |mut total, part| {
        total += part.counts;
        total
    });
    if !total.reaches_identity(min_identity) {
        added.retain(|part| part.counts.reaches_identity(min_identity));
    }
    Ok(finish_mosaic(query_length, added, alternatives))
}

/// The mosaic of the primary fragments in `added`, whose added bases are disjoint.
fn finish_mosaic(query_length: u64, added: Vec<Added>, alternatives: Vec<Fragment>) -> Mosaic {
    let mut counts = EditCounts::default();
    let mut covered = Vec::new();
    let mut primary = Vec::with_capacity(added.len());
    for part in added {
        counts += part.counts;
        covered.extend_from_slice(&part.bases);
        primary.push(SelectedFragment {
            new_bases: Interval::total_length(&part.bases),
            fragment: part.fragment,
        });
    }
    let covered = Interval::union(covered);
    Mosaic {
        query_length,
        covered_bases: Interval::total_length(&covered),
        gaps: gaps_between(&covered, query_length),
        covered_intervals: covered,
        primary,
        alternatives,
        counts,
    }
}

/// Sorted, disjoint query intervals in a map from start to end, so that adding an interval and
/// finding the parts of one inside take logarithmic time.
#[derive(Default)]
struct IntervalSet(BTreeMap<u64, u64>);

impl IntervalSet {
    /// Adds `piece`, joined with the intervals it overlaps or touches.
    fn insert(&mut self, piece: Interval) {
        if piece.is_empty() {
            return;
        }
        let (mut start, mut end) = (piece.start, piece.end);
        if let Some((&left, &left_end)) = self.0.range(..start).next_back()
            && left_end >= start
        {
            (start, end) = (left, end.max(left_end));
        }
        let joined: Vec<u64> = self.0.range(start..=end).map(|(&key, _)| key).collect();
        for key in joined {
            end = end.max(self.0.remove(&key).unwrap_or(end));
        }
        self.0.insert(start, end);
    }

    /// The parts of `piece` inside the set, in order.
    fn inside(&self, piece: Interval) -> impl Iterator<Item = Interval> + '_ {
        let left = self.0.range(..piece.start).next_back();
        let left = left.filter(|&(_, &end)| end > piece.start);
        let rest = self.0.range(piece.start..piece.end.max(piece.start));
        left.into_iter()
            .chain(rest)
            .map(move |(&start, &end)| Interval {
                start: start.max(piece.start),
                end: end.min(piece.end),
            })
    }

    /// The bases of `piece` inside the set.
    fn overlap(&self, piece: Interval) -> u64 {
        self.inside(piece).map(Interval::len).sum()
    }

    /// Adds the parts of `piece` outside the set to `output`, in order.
    fn outside(&self, piece: Interval, output: &mut Vec<Interval>) {
        let mut start = piece.start;
        for part in self.inside(piece) {
            if start < part.start {
                output.push(Interval {
                    start,
                    end: part.start,
                });
            }
            start = part.end;
        }
        if start < piece.end {
            output.push(Interval {
                start,
                end: piece.end,
            });
        }
    }
}

/// The base counts of `fragment` over the query bases outside `spanned`, from its CIGAR. A
/// match, substitution or insertion counts at its own query base; a deletion counts at the query
/// base before it, or after it at the fragment start.
fn new_counts(fragment: &Fragment, spanned: &IntervalSet) -> Result<EditCounts, DecideError> {
    let segments = &fragment.query_segments;
    let mut bases = Vec::new();
    let mut new_in = |offset: u64, length: u64| -> Result<u64, DecideError> {
        bases.clear();
        project_segments(segments, offset, length, &mut bases)?;
        let old: u64 = bases.iter().map(|&piece| spanned.overlap(piece)).sum();
        Ok(length - old)
    };
    let mut counts = EditCounts::default();
    let mut offset = 0u64;
    for run in &fragment.alignment.runs {
        let length = u64::from(run.length);
        match run.operation {
            EditOperation::Match => counts.matches += new_in(offset, length)?,
            EditOperation::Substitution => counts.substitutions += new_in(offset, length)?,
            EditOperation::Insertion => counts.insertions += new_in(offset, length)?,
            EditOperation::Deletion => {
                if new_in(offset.saturating_sub(1), 1)? == 1 {
                    counts.deletions += length;
                }
                continue;
            }
        }
        offset += length;
    }
    Ok(counts)
}

/// Best first: higher identity, higher score, more support, more matches, fewer edits. The
/// remaining fields only break ties, so that the mosaic does not depend on the input order.
fn compare_fragments(left: &Prepared, right: &Prepared) -> Ordering {
    let (a, b) = (&left.fragment, &right.fragment);
    let edits = |counts: EditCounts| counts.substitutions + counts.insertions + counts.deletions;
    // Identities compared exactly, as matches times the other's columns.
    let identity =
        |one: EditCounts, other: EditCounts| u128::from(one.matches) * u128::from(other.columns());
    identity(right.counts, left.counts)
        .cmp(&identity(left.counts, right.counts))
        .then_with(|| b.alignment.score.cmp(&a.alignment.score))
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
pub(super) fn project_segments(
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

/// The query bases outside the sorted, disjoint `covered` intervals.
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
    use crate::align::Strand;
    use crate::decide::fragment::tests::fragment;

    /// The shared test fragment with another score.
    fn scored(contig_index: u32, score: i32, segments: &[(u64, u64)], cigar: &str) -> Fragment {
        let mut scored = fragment(contig_index, segments, cigar);
        scored.alignment.score = score;
        scored
    }

    fn equal_fragment(contig_index: u32, score: i32, start: u64, end: u64) -> Fragment {
        scored(
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
    fn insertions_leave_gaps() {
        // 2= 1I 2=: the middle query base has no contig base.
        let mosaic = build_mosaic(5, &[scored(0, 7, &[(0, 5)], "2=1I2=")], &[], 0.0).unwrap();
        assert_eq!(mosaic.covered_bases, 4);
        assert_eq!(mosaic.gaps, [interval(2, 3)]);
    }

    #[test]
    fn deletions_take_nothing() {
        let mosaic = build_mosaic(4, &[scored(0, 7, &[(0, 4)], "2=3D2=")], &[], 0.0).unwrap();
        assert_eq!(mosaic.covered_bases, 4);
        assert!(mosaic.gaps.is_empty());
    }

    #[test]
    fn later_count_new() {
        let input = [
            equal_fragment(2, 7, 1, 5),
            equal_fragment(1, 8, 4, 10),
            equal_fragment(0, 10, 0, 6),
        ];
        let mosaic = build_mosaic(10, &input, &[], 0.0).unwrap();
        assert_eq!(mosaic.covered_bases, 10);
        assert_eq!(mosaic.primary.len(), 2);
        assert_eq!(mosaic.primary[0].fragment.contig_index, 0);
        assert_eq!(mosaic.primary[1].new_bases, 4);
        assert_eq!(mosaic.alternatives.len(), 1);
        let mut reversed = input.to_vec();
        reversed.reverse();
        assert_eq!(build_mosaic(10, &reversed, &[], 0.0).unwrap(), mosaic);
    }

    #[test]
    fn conserved_copy_primary() {
        let conserved = scored(0, 100, &[(0, 1000)], "1000=");
        let diverged = scored(1, 200, &[(0, 1000)], "900=100X");
        let input = [diverged, conserved.clone()];
        let mosaic = build_mosaic(1000, &input, &[], 0.95).unwrap();
        assert_eq!(mosaic.primary.len(), 1);
        assert_eq!(mosaic.primary[0].fragment.contig_index, 0);
        assert_eq!(mosaic.alternatives[0].contig_index, 1);
        // Without a limit the higher identity still comes first.
        let mosaic = build_mosaic(1000, &input, &[], 0.0).unwrap();
        assert_eq!(mosaic.primary[0].fragment.contig_index, 0);
        // At equal identity the higher score wins.
        let longer = scored(2, 300, &[(0, 1000)], "1000=");
        let mosaic = build_mosaic(1000, &[conserved, longer], &[], 0.0).unwrap();
        assert_eq!(mosaic.primary[0].fragment.contig_index, 2);
    }

    #[test]
    fn counts_new_bases() {
        let first = scored(0, 2000, &[(0, 500)], "500=");
        let count = |cigar: &str| {
            let second = scored(1, 1000, &[(0, 1000)], cigar);
            build_mosaic(1000, &[first.clone(), second], &[], 0.0)
                .unwrap()
                .counts
        };
        // Bases 0 to 500 count in the first fragment, with the gaps of the second there.
        let counts = count("498=2I3D500=");
        assert_eq!(
            (counts.matches, counts.insertions, counts.deletions),
            (1000, 0, 0)
        );
        let counts = count("600=4D400=");
        assert_eq!((counts.matches, counts.deletions), (1000, 4));
        let counts = count("600=5I395=");
        assert_eq!((counts.matches, counts.insertions), (995, 5));
        let counts = count("450=50X450=50X");
        assert_eq!((counts.matches, counts.substitutions), (950, 50));
    }

    #[test]
    fn duplicates_count_once() {
        let mut a = equal_fragment(0, 100, 0, 30);
        let mut b = equal_fragment(1, 100, 30, 60);
        b.alignment.strand = Strand::Reverse;
        let combined = build_mosaic(101, &[a.clone(), b.clone(), a.clone(), b], &[], 0.0).unwrap();
        assert_eq!(combined.covered_bases, 60);
        assert_eq!(combined.primary.len(), 2);
        assert!(combined.alternatives.is_empty());
        assert_eq!(combined.query_length, 101);
        a.contig_index = 5;
        let one = build_mosaic(101, &[a], &[], 0.0).unwrap();
        assert_eq!(one.covered_bases, 30);
    }

    #[test]
    fn low_complexity_uncounted() {
        let masked = [interval(20, 50)];
        let whole = scored(0, 100, &[(0, 100)], "60=1X39=");
        let inside = scored(1, 300, &[(25, 45)], "20=");
        let mosaic = build_mosaic(100, &[whole, inside], &masked, 0.0).unwrap();
        assert_eq!(mosaic.covered_bases, 70);
        assert_eq!(
            mosaic.covered_intervals,
            [interval(0, 20), interval(50, 100)]
        );
        assert_eq!(
            (mosaic.counts.matches, mosaic.counts.substitutions),
            (69, 1)
        );
        assert_eq!(mosaic.alternatives.len(), 1);
    }

    #[test]
    fn wrap_covers_ends() {
        let wrapped = scored(0, 8, &[(8, 10), (0, 2)], "4=");
        let mosaic = build_mosaic(10, &[wrapped], &[], 0.0).unwrap();
        assert_eq!(mosaic.covered_intervals, [interval(0, 2), interval(8, 10)]);
        assert_eq!(mosaic.gaps, [interval(2, 8)]);
    }

    #[test]
    fn wrap_insertions_uncovered() {
        let wrapped = scored(0, 8, &[(8, 10), (0, 3)], "1=1I1=1I1=");
        let mosaic = build_mosaic(10, &[wrapped], &[], 0.0).unwrap();
        assert_eq!(
            mosaic.covered_intervals,
            [interval(0, 1), interval(2, 3), interval(8, 9)]
        );
        assert_eq!((mosaic.counts.matches, mosaic.counts.insertions), (3, 2));
    }

    #[test]
    fn bad_fragments_fail() {
        assert_eq!(build_mosaic(0, &[], &[], 0.0), Err(DecideError::EmptyQuery));
        let outside = scored(0, 1, &[(5, 12)], "7=");
        assert_eq!(
            build_mosaic(10, &[outside], &[], 0.0),
            Err(DecideError::InvalidQuerySegments)
        );
        let overlap = scored(0, 1, &[(0, 4), (3, 6)], "7=");
        assert_eq!(
            build_mosaic(10, &[overlap], &[], 0.0),
            Err(DecideError::InvalidQuerySegments)
        );
        let short = scored(0, 1, &[(0, 4)], "5=");
        assert_eq!(
            build_mosaic(10, &[short], &[], 0.0),
            Err(DecideError::QuerySpanMismatch)
        );
    }
}
