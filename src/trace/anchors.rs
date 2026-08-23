//! Orientation-aware positional anchors from exact JMA seed matches.
//!
//! JMA returns occurrences for an exact `(hash, canonical_kmer)` query.  This
//! module adds query positions and resolves the relative strand from the two
//! canonicalization orientations.  Highly repetitive seeds are omitted as a
//! whole, and all final truncation is deterministic so thread count cannot
//! change the retained anchor set.

use crate::jma::{SeedOccurrence, SeedQuery};
use crate::trace::model::Strand;
use crate::trace::seeds::QuerySeed;
use serde::{Deserialize, Serialize};

/// A query seed and all occurrences returned for its exact JMA lookup.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SeedOccurrenceGroup {
    pub seed: QuerySeed,
    pub k: u8,
    pub occurrences: Vec<SeedOccurrence>,
}

impl SeedOccurrenceGroup {
    #[must_use]
    pub fn query(&self) -> SeedQuery {
        self.seed.query(self.k)
    }
}

/// One orientation-aware positional correspondence between query and target.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct Anchor {
    pub query_position: u64,
    /// Forward-stored target coordinate of the k-mer start.
    pub target_position: u64,
    pub contig_id: u32,
    pub strand: Strand,
    pub k: u8,
    pub hash: u64,
    pub canonical_kmer: u64,
    pub query_reverse: bool,
    pub target_reverse: bool,
}

impl Anchor {
    #[must_use]
    pub const fn query_interval(self) -> (u64, u64) {
        (self.query_position, self.query_position + self.k as u64)
    }

    #[must_use]
    pub const fn target_interval(self) -> (u64, u64) {
        (self.target_position, self.target_position + self.k as u64)
    }

    #[must_use]
    pub const fn orientation_matches(self) -> bool {
        matches!(self.strand, Strand::Forward)
    }
}

/// Bounded deterministic output from anchor generation.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct AnchorSet {
    pub anchors: Vec<Anchor>,
    pub repetitive_seeds: u64,
    pub skipped_occurrences: u64,
    pub truncated_anchors: u64,
}

impl AnchorSet {
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.anchors.is_empty()
    }
}

/// Generate bounded anchors from exact seed occurrence groups.
///
/// A group with more than `max_occurrences_per_seed` records is considered
/// repetitive and contributes no anchors.  A zero limit therefore produces an
/// empty set.  Input occurrence order is irrelevant.
#[must_use]
pub fn generate_anchors(
    groups: &[SeedOccurrenceGroup],
    max_occurrences_per_seed: u32,
    max_anchors: u32,
) -> AnchorSet {
    let mut anchors = Vec::new();
    let mut repetitive_seeds = 0u64;
    let mut skipped_occurrences = 0u64;

    for group in groups {
        let occurrence_count = group.occurrences.len();
        if occurrence_count > usize::try_from(max_occurrences_per_seed).unwrap_or(usize::MAX) {
            repetitive_seeds += 1;
            skipped_occurrences += occurrence_count as u64;
            continue;
        }
        for occurrence in &group.occurrences {
            anchors.push(Anchor {
                query_position: group.seed.position,
                target_position: occurrence.position,
                contig_id: occurrence.contig_id,
                strand: anchor_strand(group.seed.reverse, occurrence.reverse),
                k: group.k,
                hash: group.seed.hash,
                canonical_kmer: group.seed.canonical_kmer,
                query_reverse: group.seed.reverse,
                target_reverse: occurrence.reverse,
            });
        }
    }

    anchors.sort_by_key(anchor_sort_key);
    anchors.dedup();
    let max_anchors = usize::try_from(max_anchors).unwrap_or(usize::MAX);
    let truncated_anchors = anchors.len().saturating_sub(max_anchors) as u64;
    anchors.truncate(max_anchors);

    AnchorSet {
        anchors,
        repetitive_seeds,
        skipped_occurrences,
        truncated_anchors,
    }
}

/// Alias for callers that describe the operation as building anchors.
#[must_use]
pub fn build_anchors(
    groups: &[SeedOccurrenceGroup],
    max_occurrences_per_seed: u32,
    max_anchors: u32,
) -> AnchorSet {
    generate_anchors(groups, max_occurrences_per_seed, max_anchors)
}

/// Resolve relative orientation from canonicalization flags.
#[must_use]
pub const fn anchor_strand(query_reverse: bool, target_reverse: bool) -> Strand {
    if query_reverse == target_reverse {
        Strand::Forward
    } else {
        Strand::Reverse
    }
}

fn anchor_sort_key(anchor: &Anchor) -> (u32, u8, u64, u64, u64, u64, u8, bool, bool) {
    (
        anchor.contig_id,
        strand_rank(anchor.strand),
        anchor.query_position,
        anchor.target_position,
        anchor.hash,
        anchor.canonical_kmer,
        anchor.k,
        anchor.query_reverse,
        anchor.target_reverse,
    )
}

const fn strand_rank(strand: Strand) -> u8 {
    match strand {
        Strand::Forward => 0,
        Strand::Reverse => 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seed(position: u64, reverse: bool) -> QuerySeed {
        QuerySeed {
            position,
            hash: position + 1,
            canonical_kmer: position + 9,
            reverse,
        }
    }

    fn occurrence(contig_id: u32, position: u64, reverse: bool) -> SeedOccurrence {
        SeedOccurrence {
            contig_id,
            position,
            reverse,
        }
    }

    #[test]
    fn relative_orientation_is_xor_of_canonical_flags() {
        assert_eq!(anchor_strand(false, false), Strand::Forward);
        assert_eq!(anchor_strand(true, true), Strand::Forward);
        assert_eq!(anchor_strand(false, true), Strand::Reverse);
        assert_eq!(anchor_strand(true, false), Strand::Reverse);
    }

    #[test]
    fn repetitive_groups_are_suppressed() {
        let groups = [SeedOccurrenceGroup {
            seed: seed(4, false),
            k: 31,
            occurrences: vec![occurrence(0, 10, false), occurrence(0, 20, false)],
        }];
        let anchors = generate_anchors(&groups, 1, 20);
        assert!(anchors.anchors.is_empty());
        assert_eq!(anchors.repetitive_seeds, 1);
        assert_eq!(anchors.skipped_occurrences, 2);
    }

    #[test]
    fn output_order_and_anchor_limit_are_deterministic() {
        let groups = [
            SeedOccurrenceGroup {
                seed: seed(8, false),
                k: 31,
                occurrences: vec![occurrence(1, 40, false)],
            },
            SeedOccurrenceGroup {
                seed: seed(2, false),
                k: 31,
                occurrences: vec![occurrence(0, 20, false), occurrence(0, 10, false)],
            },
        ];
        let anchors = generate_anchors(&groups, 4, 2);
        assert_eq!(anchors.anchors.len(), 2);
        assert_eq!(anchors.anchors[0].query_position, 2);
        assert_eq!(anchors.anchors[0].target_position, 10);
        assert_eq!(anchors.truncated_anchors, 1);
    }

    #[test]
    fn duplicate_occurrences_do_not_duplicate_anchors() {
        let group = SeedOccurrenceGroup {
            seed: seed(4, false),
            k: 31,
            occurrences: vec![occurrence(0, 10, false), occurrence(0, 10, false)],
        };
        let anchors = generate_anchors(&[group], 4, 20);
        assert_eq!(anchors.anchors.len(), 1);
    }
}
