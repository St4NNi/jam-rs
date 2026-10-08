//! The mosaic: the fragments of one entry that cover the query, each query base once.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::segments::reaches;
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

/// A primary fragment, by its place in the sorted list, with the covered query bases it adds and
/// its counts over them.
struct Added {
    index: usize,
    bases: Vec<Interval>,
    counts: EditCounts,
}

/// Builds the mosaic of `fragments` from the first of four selections (see the steps below) that
/// covers `required_bases` at `min_identity`, `normal_bases` of them outside `low_complexity`.
/// Masked bases count toward coverage and identity. Fails when `stopped` returns true.
pub(super) fn build_mosaic(
    query_length: u64,
    fragments: Vec<Fragment>,
    low_complexity: &[Interval],
    min_identity: f64,
    required_bases: u64,
    normal_bases: u64,
    stopped: &dyn Fn() -> bool,
) -> Result<Mosaic, DecideError> {
    if query_length == 0 {
        return Err(DecideError::EmptyQuery);
    }
    let mut prepared = fragments
        .into_iter()
        .map(|fragment| {
            if stopped() {
                return Err(DecideError::Stopped);
            }
            let support = covered_intervals(query_length, &fragment)?;
            Ok(Prepared {
                support_bases: Interval::total_length(&support),
                support,
                counts: fragment.alignment.counts(),
                fragment,
            })
        })
        .collect::<Result<Vec<_>, DecideError>>()?;
    prepared.sort_by(compare_fragments);
    prepared.dedup_by(|left, right| left.fragment == right.fragment);
    let mut masked = IntervalSet::default();
    low_complexity
        .iter()
        .for_each(|&piece| masked.insert(piece));

    // 1: all fragments, highest identity first, each over the bases it adds.
    let all = 0..prepared.len();
    let (mut added, alternatives) = add_in_order(&prepared, all, stopped)?;
    if !reaches(total_counts(&added), min_identity) {
        // 2: without the fragments whose added bases are below the limit.
        added.retain(|part| reaches(part.counts, min_identity));
        // 3: the one fragment that covers the most bases at the limit alone.
        let enough = |added: &[Added]| {
            let bases = added.iter().map(|part| part.bases.as_slice());
            let (all, normal) = base_totals(bases, &masked);
            all >= required_bases && normal >= normal_bases
        };
        if !enough(&added) {
            let found = best_single(&prepared, &masked, min_identity, normal_bases, stopped)?;
            if let Some(single) = found
                && enough(std::slice::from_ref(&single))
            {
                return Ok(finish_mosaic(
                    query_length,
                    prepared,
                    vec![single],
                    Vec::new(),
                ));
            }
            // 4: the fragments of the one contig that covers the most bases at the limit.
            let found = best_contig(&prepared, min_identity, required_bases, &enough, stopped)?;
            if let Some((contig, rest)) = found {
                return Ok(finish_mosaic(query_length, prepared, contig, rest));
            }
        }
    }
    Ok(finish_mosaic(query_length, prepared, added, alternatives))
}

/// Adds the fragments in order: one is primary when it covers a query base that no earlier one
/// covers, and counts over the bases no earlier one spans; else it is an alternative.
fn add_in_order(
    prepared: &[Prepared],
    order: impl Iterator<Item = usize>,
    stopped: &dyn Fn() -> bool,
) -> Result<(Vec<Added>, Vec<usize>), DecideError> {
    let mut covered = IntervalSet::default();
    // The query bases that earlier primary fragments span, insertions too, for the counts.
    let mut spanned = IntervalSet::default();
    let mut added = Vec::new();
    let mut alternatives = Vec::new();
    for index in order {
        let item = &prepared[index];
        if stopped() {
            return Err(DecideError::Stopped);
        }
        let mut bases = Vec::new();
        for &piece in &item.support {
            covered.outside(piece, &mut bases);
        }
        if bases.is_empty() {
            alternatives.push(index);
            continue;
        }
        let counts = new_counts(&item.fragment, &spanned)?;
        item.fragment
            .query_segments
            .iter()
            .for_each(|&piece| spanned.insert(piece));
        bases.iter().for_each(|&piece| covered.insert(piece));
        added.push(Added {
            index,
            bases,
            counts,
        });
    }
    Ok((added, alternatives))
}

/// The fragment that alone covers the most query bases with identity at least `min_identity`,
/// among those with `normal_bases` outside `masked`; the first in order on ties.
fn best_single(
    prepared: &[Prepared],
    masked: &IntervalSet,
    min_identity: f64,
    normal_bases: u64,
    stopped: &dyn Fn() -> bool,
) -> Result<Option<Added>, DecideError> {
    let mut best: Option<Added> = None;
    for (index, item) in prepared.iter().enumerate() {
        if stopped() {
            return Err(DecideError::Stopped);
        }
        let counts = item.counts;
        if !reaches(counts, min_identity) {
            continue;
        }
        let bases = item.support.clone();
        let (length, normal) = base_totals([bases.as_slice()], masked);
        if normal >= normal_bases
            && best
                .as_ref()
                .is_none_or(|old| length > Interval::total_length(&old.bases))
        {
            best = Some(Added {
                index,
                bases,
                counts,
            });
        }
    }
    Ok(best)
}

/// A contig selection: the added pieces and the indices of the fragments they came from.
type ContigPick = (Vec<Added>, Vec<usize>);

/// The fragments of the one contig that alone cover the most bases, enough by `enough`, at
/// `min_identity`; the higher identity wins a tie. Each contig is tried as
/// selections 1 and 2 try all fragments. Contigs with fewer paired bases than needed are skipped.
fn best_contig(
    prepared: &[Prepared],
    min_identity: f64,
    required_bases: u64,
    enough: &dyn Fn(&[Added]) -> bool,
    stopped: &dyn Fn() -> bool,
) -> Result<Option<ContigPick>, DecideError> {
    let mut contigs: BTreeMap<u32, (u64, Vec<usize>)> = BTreeMap::new();
    for (index, item) in prepared.iter().enumerate() {
        let entry = contigs.entry(item.fragment.contig_index).or_default();
        entry.0 += item.support_bases;
        entry.1.push(index);
    }
    let mut best: Option<ContigPick> = None;
    let mut best_key = (0, EditCounts::default());
    for (span, order) in contigs.into_values() {
        if span < required_bases {
            continue;
        }
        let (mut added, alternatives) = add_in_order(prepared, order.into_iter(), stopped)?;
        if !reaches(total_counts(&added), min_identity) {
            added.retain(|part| reaches(part.counts, min_identity));
        }
        let counts = total_counts(&added);
        let bases: u64 = added.iter().map(|p| Interval::total_length(&p.bases)).sum();
        if !enough(&added) || !reaches(counts, min_identity) {
            continue;
        }
        let better = best.is_none()
            || bases > best_key.0
            || (bases == best_key.0
                && u128::from(counts.matches) * u128::from(best_key.1.columns())
                    > u128::from(best_key.1.matches) * u128::from(counts.columns()));
        if better {
            best_key = (bases, counts);
            best = Some((added, alternatives));
        }
    }
    Ok(best)
}

/// The number of bases in `pieces`, and of those outside `masked`.
fn base_totals<'a>(
    pieces: impl IntoIterator<Item = &'a [Interval]>,
    masked: &IntervalSet,
) -> (u64, u64) {
    let (mut all, mut normal) = (0, 0);
    for piece in pieces.into_iter().flatten() {
        let mut outside = Vec::new();
        masked.outside(*piece, &mut outside);
        all += piece.end - piece.start;
        normal += Interval::total_length(&outside);
    }
    (all, normal)
}

/// The query bases of `covered` outside `low_complexity`.
pub(super) fn normal_covered(covered: &[Interval], low_complexity: &[Interval]) -> u64 {
    let mut masked = IntervalSet::default();
    low_complexity
        .iter()
        .for_each(|&piece| masked.insert(piece));
    base_totals([covered], &masked).1
}

/// The summed counts of `added`.
fn total_counts(added: &[Added]) -> EditCounts {
    added.iter().fold(EditCounts::default(), |mut total, part| {
        total += part.counts;
        total
    })
}

/// The mosaic of the primary fragments in `added`, whose added bases are disjoint, and the
/// alternatives, both by place in `prepared`.
fn finish_mosaic(
    query_length: u64,
    prepared: Vec<Prepared>,
    added: Vec<Added>,
    alternatives: Vec<usize>,
) -> Mosaic {
    let mut fragments: Vec<Option<Fragment>> = prepared
        .into_iter()
        .map(|item| Some(item.fragment))
        .collect();
    let counts = total_counts(&added);
    let mut covered = Vec::new();
    let mut primary = Vec::with_capacity(added.len());
    for part in added {
        covered.extend_from_slice(&part.bases);
        if let Some(fragment) = fragments[part.index].take() {
            primary.push(SelectedFragment {
                new_bases: Interval::total_length(&part.bases),
                fragment,
            });
        }
    }
    let alternatives = alternatives.into_iter();
    let alternatives = alternatives.filter_map(|index| fragments[index].take());
    let covered = Interval::union(covered);
    Mosaic {
        query_length,
        covered_bases: Interval::total_length(&covered),
        gaps: gaps_between(&covered, query_length),
        covered_intervals: covered,
        primary,
        alternatives: alternatives.collect(),
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
        let mosaic = build_mosaic(
            5,
            vec![scored(0, 7, &[(0, 5)], "2=1I2=")],
            &[],
            0.0,
            1,
            1,
            &|| false,
        )
        .unwrap();
        assert_eq!(mosaic.covered_bases, 4);
        assert_eq!(mosaic.gaps, [interval(2, 3)]);
    }

    #[test]
    fn deletions_take_nothing() {
        let mosaic = build_mosaic(
            4,
            vec![scored(0, 7, &[(0, 4)], "2=3D2=")],
            &[],
            0.0,
            1,
            1,
            &|| false,
        )
        .unwrap();
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
        let mosaic = build_mosaic(10, input.to_vec(), &[], 0.0, 1, 1, &|| false).unwrap();
        assert_eq!(mosaic.covered_bases, 10);
        assert_eq!(mosaic.primary.len(), 2);
        assert_eq!(mosaic.primary[0].fragment.contig_index, 0);
        assert_eq!(mosaic.primary[1].new_bases, 4);
        assert_eq!(mosaic.alternatives.len(), 1);
        let mut reversed = input.to_vec();
        reversed.reverse();
        assert_eq!(
            build_mosaic(10, reversed.to_vec(), &[], 0.0, 1, 1, &|| false).unwrap(),
            mosaic
        );
    }

    #[test]
    fn conserved_copy_primary() {
        let conserved = scored(0, 100, &[(0, 1000)], "1000=");
        let diverged = scored(1, 200, &[(0, 1000)], "900=100X");
        let input = [diverged, conserved.clone()];
        let mosaic = build_mosaic(1000, input.to_vec(), &[], 0.95, 1, 1, &|| false).unwrap();
        assert_eq!(mosaic.primary.len(), 1);
        assert_eq!(mosaic.primary[0].fragment.contig_index, 0);
        assert_eq!(mosaic.alternatives[0].contig_index, 1);
        // Without a limit the higher identity still comes first.
        let mosaic = build_mosaic(1000, input.to_vec(), &[], 0.0, 1, 1, &|| false).unwrap();
        assert_eq!(mosaic.primary[0].fragment.contig_index, 0);
        // At equal identity the higher score wins.
        let longer = scored(2, 300, &[(0, 1000)], "1000=");
        let mosaic =
            build_mosaic(1000, vec![conserved, longer], &[], 0.0, 1, 1, &|| false).unwrap();
        assert_eq!(mosaic.primary[0].fragment.contig_index, 2);
    }

    #[test]
    fn counts_new_bases() {
        let first = scored(0, 2000, &[(0, 500)], "500=");
        let count = |cigar: &str| {
            let second = scored(1, 1000, &[(0, 1000)], cigar);
            build_mosaic(1000, vec![first.clone(), second], &[], 0.0, 1, 1, &|| false)
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
        let combined = build_mosaic(
            101,
            vec![a.clone(), b.clone(), a.clone(), b],
            &[],
            0.0,
            1,
            1,
            &|| false,
        )
        .unwrap();
        assert_eq!(combined.covered_bases, 60);
        assert_eq!(combined.primary.len(), 2);
        assert!(combined.alternatives.is_empty());
        assert_eq!(combined.query_length, 101);
        a.contig_index = 5;
        let one = build_mosaic(101, vec![a], &[], 0.0, 1, 1, &|| false).unwrap();
        assert_eq!(one.covered_bases, 30);
    }

    #[test]
    fn low_complexity_covered() {
        let masked = [interval(20, 50)];
        let whole = scored(0, 100, &[(0, 100)], "60=1X39=");
        let inside = scored(1, 300, &[(25, 45)], "20=");
        let mosaic = build_mosaic(100, vec![whole, inside], &masked, 0.0, 1, 1, &|| false).unwrap();
        assert_eq!(mosaic.covered_bases, 100);
        assert_eq!(mosaic.covered_intervals, [interval(0, 100)]);
        assert_eq!(
            (mosaic.counts.matches, mosaic.counts.substitutions),
            (99, 1)
        );
        // The exact fragment inside the masked bases ranks first and now adds its bases.
        assert_eq!(mosaic.primary.len(), 2);
        assert!(mosaic.alternatives.is_empty());
    }

    #[test]
    fn wrap_covers_ends() {
        let wrapped = scored(0, 8, &[(8, 10), (0, 2)], "4=");
        let mosaic = build_mosaic(10, vec![wrapped], &[], 0.0, 1, 1, &|| false).unwrap();
        assert_eq!(mosaic.covered_intervals, [interval(0, 2), interval(8, 10)]);
        assert_eq!(mosaic.gaps, [interval(2, 8)]);
    }

    #[test]
    fn wrap_insertions_uncovered() {
        let wrapped = scored(0, 8, &[(8, 10), (0, 3)], "1=1I1=1I1=");
        let mosaic = build_mosaic(10, vec![wrapped], &[], 0.0, 1, 1, &|| false).unwrap();
        assert_eq!(
            mosaic.covered_intervals,
            [interval(0, 1), interval(2, 3), interval(8, 9)]
        );
        assert_eq!((mosaic.counts.matches, mosaic.counts.insertions), (3, 2));
    }

    #[test]
    fn bad_fragments_fail() {
        assert_eq!(
            build_mosaic(0, vec![], &[], 0.0, 1, 1, &|| false),
            Err(DecideError::EmptyQuery)
        );
        let outside = scored(0, 1, &[(5, 12)], "7=");
        assert_eq!(
            build_mosaic(10, vec![outside], &[], 0.0, 1, 1, &|| false),
            Err(DecideError::InvalidQuerySegments)
        );
        let overlap = scored(0, 1, &[(0, 4), (3, 6)], "7=");
        assert_eq!(
            build_mosaic(10, vec![overlap], &[], 0.0, 1, 1, &|| false),
            Err(DecideError::InvalidQuerySegments)
        );
        let short = scored(0, 1, &[(0, 4)], "5=");
        assert_eq!(
            build_mosaic(10, vec![short], &[], 0.0, 1, 1, &|| false),
            Err(DecideError::QuerySpanMismatch)
        );
    }
}
