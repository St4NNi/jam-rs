//! The paired bases of the alignments already made from one alignment window, so that a chain or
//! a dense walk that only finds one of them again is not aligned a second time.
//!
//! Positions are oriented: window positions are read on the alignment's strand (the reverse
//! complement of the window for the reverse strand), and query positions run on past the origin
//! of a circular query without wrapping. The diagonal of a pair is its window position minus its
//! query position.

#[cfg(test)]
use crate::align::Walk;
use crate::align::{Alignment, Chain, EditOperation, Strand};
use std::collections::HashMap;

/// The paired stretches of one alignment: (query position, oriented window position, length).
pub(super) struct AlignedPath {
    reverse: bool,
    paired: Vec<(u64, u64, u64)>,
}

impl AlignedPath {
    /// The path of `alignment`, which was made from a window of `window_length` bases that starts
    /// at contig position `offset`.
    pub(super) fn of(alignment: &Alignment, window_length: u64, offset: u64) -> Self {
        let reverse = alignment.strand == Strand::Reverse;
        let contig = alignment.contig_interval;
        let mut window_at = if reverse {
            window_length - (contig.end - offset)
        } else {
            contig.start - offset
        };
        let mut query_at = alignment.query_interval.start;
        let mut paired = Vec::new();
        for run in &alignment.runs {
            let length = u64::from(run.length);
            match run.operation {
                EditOperation::Deletion => window_at += length,
                EditOperation::Insertion => query_at += length,
                EditOperation::Match | EditOperation::Substitution => {
                    paired.push((query_at, window_at, length));
                    query_at += length;
                    window_at += length;
                }
            }
        }
        Self { reverse, paired }
    }

    #[cfg(test)]
    /// The paired stretches as walks, with query positions wrapped into a query of
    /// `query_length` bases.
    pub(super) fn walks(&self, query_length: u64) -> impl Iterator<Item = Walk> + '_ {
        self.paired
            .iter()
            .map(move |&(query, window, length)| Walk {
                reverse: self.reverse,
                contig: window as u32,
                query: (query % query_length) as u32,
                length: length as u32,
            })
    }
}

/// Paired stretches by strand and diagonal, as (query start, length).
type Segments = HashMap<(bool, i64), Vec<(u64, u64)>>;

/// The paired stretches of the alignments made from one window, by strand and diagonal, so a
/// chain is compared only with the stretches on its own diagonals. Touching stretches of one
/// diagonal are kept as one.
#[derive(Default)]
pub(super) struct PathIndex {
    /// Stretches of every alignment, dense walks included.
    diagonals: Segments,
    /// The same stretches of kernel alignments only, without the dense walks.
    aligned: Segments,
}

fn add_segments(segments: &mut Segments, path: &AlignedPath) {
    for &(query, window, length) in &path.paired {
        let list = segments
            .entry((path.reverse, window as i64 - query as i64))
            .or_default();
        match list.last_mut() {
            Some(last) if last.0 + last.1 == query => last.1 += length,
            _ => list.push((query, length)),
        }
    }
}

#[cfg(test)]
/// The bases covered by the union of the half-open spans.
fn covered_bases(mut spans: Vec<(i64, i64)>) -> u64 {
    spans.sort_unstable();
    let (mut covered, mut reached) = (0u64, i64::MIN);
    for (from, to) in spans {
        let from = from.max(reached);
        if to > from {
            covered += (to - from) as u64;
            reached = to;
        }
    }
    covered
}

#[cfg(test)]
/// The query lengths a circular position may be shifted by to meet a stored stretch.
fn shifts(circular: bool) -> &'static [i64] {
    if circular { &[0, 1, -1] } else { &[0] }
}

impl PathIndex {
    /// Adds the path of a kernel alignment.
    pub(super) fn add(&mut self, path: &AlignedPath) {
        add_segments(&mut self.diagonals, path);
        add_segments(&mut self.aligned, path);
    }

    #[cfg(test)]
    /// Adds an ungapped dense walk.
    pub(super) fn add_walk(&mut self, walk: &Walk) {
        let path = AlignedPath {
            reverse: walk.reverse,
            paired: vec![(walk.query.into(), walk.contig.into(), walk.length.into())],
        };
        add_segments(&mut self.diagonals, &path);
    }

    /// Whether an anchor lies on a paired stretch; a circular query position may also count one
    /// query length on.
    fn covers(&self, reverse: bool, window: u64, query: u64, query_length: u64) -> bool {
        [query, query + query_length].iter().any(|&query| {
            self.diagonals
                .get(&(reverse, window as i64 - query as i64))
                .is_some_and(|list| {
                    list.iter()
                        .any(|&(start, length)| (start..start + length).contains(&query))
                })
        })
    }

    /// A chain is redundant when nine tenths of its anchors lie on stretches already aligned. An
    /// anchor of another copy (another diagonal) does not count, so alternative alignments of
    /// repeats inside an aligned region are kept.
    pub(super) fn redundant(&self, chain: &Chain, query_length: u64) -> bool {
        let covered = chain
            .anchors
            .iter()
            .filter(|&&(window, query)| {
                self.covers(chain.reverse, window.into(), query.into(), query_length)
            })
            .count();
        covered * 10 >= chain.anchors.len() * 9
    }

    /// The anchors of `chain` that no stretch covers, as a chain of their own. `None` when none
    /// of them is covered (nothing was gained) or fewer than `min_anchors` are left.
    pub(super) fn uncovered(
        &self,
        chain: &Chain,
        query_length: u64,
        min_anchors: u32,
    ) -> Option<Chain> {
        let kept: Vec<usize> = (0..chain.anchors.len())
            .filter(|&index| {
                let (window, query) = chain.anchors[index];
                !self.covers(chain.reverse, window.into(), query.into(), query_length)
            })
            .collect();
        if kept.len() == chain.anchors.len() || kept.len() < min_anchors as usize {
            return None;
        }
        let base = chain.offsets[kept[0]];
        Some(Chain {
            reverse: chain.reverse,
            anchors: kept.iter().map(|&index| chain.anchors[index]).collect(),
            offsets: kept
                .iter()
                .map(|&index| chain.offsets[index] - base)
                .collect(),
            score: chain.score,
        })
    }

    #[cfg(test)]
    /// Whether a dense walk lies on stretches of its own diagonal: nine tenths of it for a linear
    /// query, all of it for a circular one, so a shifted copy that adds contig bases stays.
    pub(super) fn covers_walk(&self, walk: &Walk, query_length: u64, circular: bool) -> bool {
        let (start, length) = (i64::from(walk.query), i64::from(walk.length));
        let diagonal = i64::from(walk.contig) - start;
        // Overlaps with every shifted copy, moved back to the walk's query positions and merged
        // once, so a walk base covered by two copies counts once.
        let mut overlaps: Vec<(i64, i64)> = Vec::new();
        for &shift in shifts(circular) {
            let moved = shift * query_length as i64;
            let Some(list) = self.diagonals.get(&(walk.reverse, diagonal - moved)) else {
                continue;
            };
            let (low, high) = (start + moved, start + length + moved);
            overlaps.extend(
                list.iter()
                    .map(|&(at, run)| {
                        let (at, end) = (at as i64, (at + run) as i64);
                        (at.max(low) - moved, end.min(high) - moved)
                    })
                    .filter(|(from, to)| from < to),
            );
        }
        let required = if circular { 10 } else { 9 };
        covered_bases(overlaps).min(length as u64) * 10 >= length as u64 * required
    }

    #[cfg(test)]
    /// Whether kernel alignments within `band` diagonals of a dense walk already pair nine tenths
    /// of its query bases and nine tenths of its window bases. A gapped alignment of such a walk
    /// would only find that region again. Dense walks do not count: the gapped join between them
    /// is what an alignment of a walk looks for.
    pub(super) fn aligned_around(
        &self,
        walk: &Walk,
        band: i64,
        query_length: u64,
        circular: bool,
    ) -> bool {
        let (start, length) = (i64::from(walk.query), i64::from(walk.length));
        let window = i64::from(walk.contig);
        let diagonal = window - start;
        let (mut queries, mut windows) = (Vec::new(), Vec::new());
        for &shift in shifts(circular) {
            let moved = shift * query_length as i64;
            let (low, high) = (start + moved, start + length + moved);
            for line in diagonal - moved - band..=diagonal - moved + band {
                let Some(list) = self.aligned.get(&(walk.reverse, line)) else {
                    continue;
                };
                for &(at, run) in list {
                    let (at, end) = (at as i64, (at + run) as i64);
                    if at.max(low) < end.min(high) {
                        queries.push((at.max(low) - moved, end.min(high) - moved));
                    }
                    let (from, to) = ((at + line).max(window), (end + line).min(window + length));
                    if from < to {
                        windows.push((from, to));
                    }
                }
            }
        }
        let nine_tenths = |spans| covered_bases(spans).min(length as u64) * 10 >= length as u64 * 9;
        nine_tenths(queries) && nine_tenths(windows)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{Interval, parse_cigar};

    fn alignment(strand: Strand, query: (u64, u64), contig: (u64, u64), cigar: &str) -> Alignment {
        let alignment = Alignment {
            score: 1,
            strand,
            query_interval: Interval::new(query.0, query.1).unwrap(),
            contig_interval: Interval::new(contig.0, contig.1).unwrap(),
            runs: parse_cigar(cigar).unwrap(),
        };
        alignment.validate().unwrap();
        alignment
    }

    fn chain(reverse: bool, anchors: &[(u32, u32)]) -> Chain {
        Chain {
            reverse,
            anchors: anchors.to_vec(),
            offsets: anchors.iter().map(|&(_, q)| q - anchors[0].1).collect(),
            score: 100,
        }
    }

    #[test]
    fn insertions_move_the_query_and_deletions_move_the_window() {
        // Window bases 100..130 of a window starting at contig position 1,000.
        let found = alignment(Strand::Forward, (10, 35), (1100, 1130), "10=5I5=10D5=");
        let path = AlignedPath::of(&found, 500, 1000);
        assert_eq!(path.paired, [(10, 100, 10), (25, 110, 5), (30, 125, 5)]);
        // On the reverse strand the window is read from its end.
        let reverse = alignment(Strand::Reverse, (0, 20), (1400, 1420), "20=");
        let path = AlignedPath::of(&reverse, 500, 1000);
        assert_eq!(path.paired, [(0, 80, 20)]);
        let walks: Vec<Walk> = AlignedPath::of(&found, 500, 1000).walks(28).collect();
        assert_eq!(walks[2].query, 2);
    }

    #[test]
    fn a_chain_on_an_aligned_diagonal_is_redundant_and_its_rest_is_kept() {
        let mut paths = PathIndex::default();
        let found = alignment(Strand::Forward, (0, 100), (0, 100), "100=");
        paths.add(&AlignedPath::of(&found, 1000, 0));
        let inside = chain(false, &[(10, 10), (40, 40), (80, 80)]);
        assert!(paths.redundant(&inside, 1000));
        let other_copy = chain(false, &[(510, 10), (540, 40), (580, 80)]);
        assert!(!paths.redundant(&other_copy, 1000));
        let longer = chain(
            false,
            &[(50, 50), (90, 90), (300, 300), (400, 400), (500, 500)],
        );
        let rest = paths.uncovered(&longer, 1000, 3).unwrap();
        assert_eq!(rest.anchors, [(300, 300), (400, 400), (500, 500)]);
        assert_eq!(rest.offsets, [0, 100, 200]);
        assert!(paths.uncovered(&other_copy, 1000, 3).is_none());
    }

    #[test]
    fn circular_positions_meet_stretches_one_query_length_on() {
        let mut paths = PathIndex::default();
        // An alignment over the origin of a query of 100 bases: query 90..130.
        let found = alignment(Strand::Forward, (90, 130), (190, 230), "40=");
        paths.add(&AlignedPath::of(&found, 1000, 0));
        // Query position 10 is 110 one query length on, on window position 210.
        assert!(paths.redundant(&chain(false, &[(205, 5), (210, 10), (215, 15)]), 100));
        let walk = Walk {
            reverse: false,
            contig: 200,
            query: 0,
            length: 30,
        };
        assert!(paths.covers_walk(&walk, 100, true));
        assert!(!paths.covers_walk(&walk, 100, false));
    }

    #[test]
    fn walks_count_for_coverage_but_not_for_alignments_around() {
        let walk = Walk {
            reverse: true,
            contig: 100,
            query: 0,
            length: 200,
        };
        let mut paths = PathIndex::default();
        paths.add_walk(&walk);
        assert!(paths.covers_walk(&walk, 1000, false));
        assert!(!paths.aligned_around(&walk, 64, 1000, false));
        let mut paths = PathIndex::default();
        // A kernel alignment two diagonals off pairs the same bases.
        let found = alignment(Strand::Reverse, (0, 200), (698, 898), "200=");
        paths.add(&AlignedPath::of(&found, 1000, 0));
        assert!(paths.aligned_around(&walk, 64, 1000, false));
        assert!(!paths.aligned_around(&walk, 1, 1000, false));
    }
}
