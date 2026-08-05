//! The contig rule: which contigs of an entry count for the decision.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::Fragment;
use crate::align::{EditCounts, Interval};
use std::collections::{HashMap, HashSet};

/// The fragments of the contigs that count: `min_coverage` of their bases lie in fragments that
/// reach `min_identity` together. A contig that only shares a gene or a repeat with the query
/// is left out, so it cannot lower the identity. `contig_length` gives a contig length.
pub(super) fn counted_contig_fragments(
    fragments: Vec<Fragment>,
    min_coverage: f64,
    min_identity: f64,
    contig_length: impl Fn(u32) -> u64,
) -> Vec<Fragment> {
    let mut spans: HashMap<u32, Vec<Interval>> = HashMap::new();
    let mut pairing: HashMap<u32, EditCounts> = HashMap::new();
    for fragment in &fragments {
        let contig = fragment.contig_index;
        spans
            .entry(contig)
            .or_default()
            .push(fragment.alignment.contig_interval);
        *pairing.entry(contig).or_default() += fragment.alignment.counts();
    }
    let counted: HashSet<u32> = pairing
        .into_iter()
        .filter(|&(contig, counts)| {
            let union = Interval::union(spans.remove(&contig).unwrap_or_default());
            let aligned = Interval::total_length(&union);
            aligned as f64 >= min_coverage * contig_length(contig) as f64
                && counts.reaches_identity(min_identity)
        })
        .map(|(contig, _)| contig)
        .collect();
    fragments
        .into_iter()
        .filter(|fragment| counted.contains(&fragment.contig_index))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{Alignment, Interval, Strand, parse_cigar};

    fn fragment(contig: u32, contig_start: u64, cigar: &str) -> Fragment {
        let runs = parse_cigar(cigar).unwrap();
        let counts = Alignment {
            score: 0,
            strand: Strand::Forward,
            query_interval: Interval::new(0, 0).unwrap(),
            contig_interval: Interval::new(0, 0).unwrap(),
            runs: runs.clone(),
        }
        .counts();
        let query = counts.matches + counts.substitutions + counts.insertions;
        let on_contig = counts.matches + counts.substitutions + counts.deletions;
        Fragment {
            contig_index: contig,
            query_segments: vec![Interval::new(0, query).unwrap()],
            alignment: Alignment {
                score: 0,
                strand: Strand::Forward,
                query_interval: Interval::new(0, query).unwrap(),
                contig_interval: Interval::new(contig_start, contig_start + on_contig).unwrap(),
                runs,
            },
        }
    }

    fn kept(fragments: &[Fragment], coverage: f64, identity: f64) -> Vec<u32> {
        let lengths = [3700u64, 10_000];
        let kept = counted_contig_fragments(fragments.to_vec(), coverage, identity, |contig| {
            lengths[contig as usize]
        });
        kept.iter().map(|fragment| fragment.contig_index).collect()
    }

    #[test]
    fn a_contig_that_is_mostly_not_the_query_is_left_out() {
        // Contig 0 is 81 percent paired. Contig 1 shares 600 of 10,000 bases.
        let fragments = [fragment(0, 400, "3000="), fragment(1, 0, "600=")];
        assert_eq!(kept(&fragments, 0.0, 0.0), [0, 1]);
        assert_eq!(kept(&fragments, 0.5, 0.0), [0]);
        assert!(kept(&fragments, 0.9, 0.0).is_empty());
    }

    #[test]
    fn overlapping_fragments_of_a_contig_count_their_union() {
        let fragments = [fragment(0, 0, "2000="), fragment(0, 1000, "2000=")];
        // The union is 3,000 of 3,700 bases (81 percent), not 4,000.
        assert_eq!(kept(&fragments, 0.8, 0.0), [0, 0]);
        assert!(kept(&fragments, 0.82, 0.0).is_empty());
    }

    #[test]
    fn the_identity_is_taken_over_all_fragments_of_a_contig() {
        let fragments = [fragment(0, 0, "1800=200X"), fragment(0, 2000, "1000=")];
        // 2,800 matches in 3,000 columns.
        assert_eq!(kept(&fragments, 0.0, 0.93), [0, 0]);
        assert!(kept(&fragments, 0.0, 0.94).is_empty());
    }
}
