//! Stage 3 on planted copies with substitutions only: from one located window inside a copy,
//! the counted fragments cover the whole copy at the identity limit or above.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::tests::build;
use super::*;
use crate::align::Interval;
use crate::decide::Fragment;
use crate::search::ContigWindows;
use crate::seeds::{complement, reverse_complement};
use crate::test_support::{no_stop, random_bases, rng};
use rand::RngExt;
use rand::rngs::StdRng;

/// Query bases before and after the copy.
const FLANK: usize = 300;

/// Where the substitutions of a copy lie.
#[derive(Clone, Copy, Debug)]
enum Placement {
    /// 4% of the bases, uniformly at random.
    Uniform,
    /// Clusters of 6 substitutions within 12 bases, 5% of the bases in all.
    Clustered,
    /// One every 20 bases.
    Periodic,
    /// Exactly 5% of the bases (the identity limit), uniformly at random.
    Limit,
    /// Exactly this many bases, uniformly at random.
    Count(usize),
}

/// Substitution positions in a copy of `length` bases. The first and last 10 bases stay equal,
/// so that every prefix and every suffix of the copy scores above zero.
fn substitutions(placement: Placement, length: usize, random: &mut StdRng) -> Vec<usize> {
    let inner = 10..length - 10;
    let mut picked: Vec<usize> = match placement {
        Placement::Uniform | Placement::Limit | Placement::Count(_) => {
            let count = match placement {
                Placement::Uniform => length * 4 / 100,
                Placement::Count(count) => count,
                _ => length * 5 / 100,
            };
            let mut all: Vec<usize> = inner.collect();
            for at in 0..count.min(all.len()) {
                let other = random.random_range(at..all.len());
                all.swap(at, other);
            }
            all.truncate(count);
            all
        }
        Placement::Clustered => {
            let clusters = length / 20 / 6;
            let mut all = Vec::new();
            for cluster in 0..clusters {
                // One cluster per stretch of 120 bases: 6 of its first 12 bases.
                let start = 10 + cluster * (length - 20) / clusters.max(1);
                let mut bases: Vec<usize> = (start..start + 12).collect();
                for at in 0..6 {
                    let other = random.random_range(at..12);
                    bases.swap(at, other);
                }
                all.extend(bases[..6].iter().filter(|&&at| at < length - 10));
            }
            all
        }
        Placement::Periodic => (10..length).step_by(20).collect(),
    };
    picked.sort_unstable();
    picked.dedup();
    picked
}

/// A copy of `bases` with a different base at each of `positions`.
fn substituted(bases: &[u8], positions: &[usize]) -> Vec<u8> {
    let mut copy = bases.to_vec();
    for &at in positions {
        copy[at] = match copy[at] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
    }
    copy
}

/// One planted case: the query, its copy in a contig and the contig range of one located bin.
struct Case {
    query: Vec<u8>,
    contig: Vec<u8>,
    window: std::ops::Range<u64>,
    substitutions: usize,
}

/// The query is flank, interval, flank; the contig holds the copy between the complements of
/// the flanks, so no column next to the copy matches. The located bin lies around the middle.
fn case(length: usize, placement: Placement, reverse: bool, seed: u64) -> Case {
    let mut random = rng(seed);
    let (before, interval, after) = (
        random_bases(&mut random, FLANK),
        random_bases(&mut random, length),
        random_bases(&mut random, FLANK),
    );
    let positions = substitutions(placement, length, &mut random);
    let copy = substituted(&interval, &positions);
    let opposite = |bases: &[u8]| {
        bases
            .iter()
            .map(|&base| complement(base))
            .collect::<Vec<_>>()
    };
    let forward = [opposite(&before), copy, opposite(&after)].concat();
    let contig = if reverse {
        reverse_complement(&forward)
    } else {
        forward
    };
    let middle = (FLANK + length / 2) as u64;
    let bin = middle >> 10 << 10;
    Case {
        query: [before, interval, after].concat(),
        window: bin..(bin + 1024).min(contig.len() as u64),
        contig,
        substitutions: positions.len(),
    }
}

/// The query bases that the counted fragments cover.
fn counted_cover(fragments: &[Fragment], rule: &DecisionRule) -> Vec<Interval> {
    let counted = fragments
        .iter()
        .filter(|fragment| rule.counts_fragment(&fragment.alignment));
    Interval::union(counted.flat_map(|f| f.query_segments.clone()).collect())
}

#[test]
fn substituted_copies_covered() {
    let rule = DecisionRule::default();
    let config = AlignStageConfig::default();
    let placements = [
        Placement::Uniform,
        Placement::Clustered,
        Placement::Periodic,
        Placement::Limit,
    ];
    let mut cases = Vec::new();
    for (index, &length) in [400usize, 1000, 2500, 20_000].iter().enumerate() {
        for (kind, &placement) in placements.iter().enumerate() {
            for reverse in [false, true] {
                let seed = 1000 + (index * 8 + kind * 2) as u64 + u64::from(reverse);
                cases.push((
                    length,
                    placement,
                    reverse,
                    case(length, placement, reverse, seed),
                ));
            }
        }
    }
    let directory = tempfile::tempdir().unwrap();
    let entries: Vec<_> = cases
        .iter()
        .enumerate()
        .map(|(at, (.., case))| {
            (
                format!("e{at}"),
                vec![("c".to_string(), case.contig.clone())],
            )
        })
        .collect();
    let shard = build(directory.path(), &entries);
    let workspace = AlignStageWorkspace::default();
    for (entry, (length, placement, reverse, case)) in cases.iter().enumerate() {
        let identity = 1.0 - case.substitutions as f64 / *length as f64;
        assert!(
            identity >= rule.min_fragment_identity,
            "{placement:?} {length}"
        );
        let located = LocatedEntry {
            entry: entry as u32,
            coverage_bound: case.query.len() as u64,
            containment_estimate: 1.0,
            alignment_windows: vec![ContigWindows {
                contig: 0,
                ranges: vec![case.window.clone()],
            }],
        };
        let prepared = PreparedQuery::new(&case.query, false, &config).unwrap();
        let packed = shard.packed(entry as u32).unwrap();
        let (fragments, _) = entry_fragments(
            &packed,
            &prepared,
            &located,
            &config,
            &rule,
            &workspace,
            no_stop(),
        )
        .unwrap();
        let copy = Interval::new(FLANK as u64, (FLANK + length) as u64).unwrap();
        assert_eq!(
            counted_cover(&fragments, &rule),
            [copy],
            "{placement:?} length {length} reverse {reverse}"
        );
    }
}

/// Short copies at 90% and 95% with the rule at that identity: one fragment counts whole over
/// the whole copy, and the decision confirms the entry.
#[test]
fn short_copies_confirmed() {
    let config = AlignStageConfig::default();
    let mut cases = Vec::new();
    for length in [150usize, 250, 390] {
        for identity in [0.9, 0.95] {
            for reverse in [false, true] {
                let count = (length as f64 * (1.0 - identity) + 1e-6) as usize;
                let seed = 2000 + (length + count) as u64 * 2 + u64::from(reverse);
                let planted = case(length, Placement::Count(count), reverse, seed);
                cases.push((length, identity, reverse, planted));
            }
        }
    }
    let directory = tempfile::tempdir().unwrap();
    let entries: Vec<_> = cases
        .iter()
        .enumerate()
        .map(|(at, (.., case))| {
            (
                format!("e{at}"),
                vec![("c".to_string(), case.contig.clone())],
            )
        })
        .collect();
    let shard = build(directory.path(), &entries);
    let workspace = AlignStageWorkspace::default();
    for (entry, (length, identity, reverse, case)) in cases.iter().enumerate() {
        let rule = DecisionRule {
            min_query_coverage: format!("{length}bp").parse().unwrap(),
            min_identity: *identity,
            min_fragment_identity: *identity,
            ..DecisionRule::default()
        };
        let located = LocatedEntry {
            entry: entry as u32,
            coverage_bound: case.query.len() as u64,
            containment_estimate: 1.0,
            alignment_windows: vec![ContigWindows {
                contig: 0,
                ranges: vec![case.window.clone()],
            }],
        };
        let prepared = PreparedQuery::new(&case.query, false, &config).unwrap();
        let packed = shard.packed(entry as u32).unwrap();
        let control = no_stop();
        let found = entry_fragments(
            &packed, &prepared, &located, &config, &rule, &workspace, control,
        );
        let (fragments, _) = found.unwrap();
        let copy = Interval::new(FLANK as u64, (FLANK + length) as u64).unwrap();
        let whole = fragments.iter().filter(|fragment| {
            rule.counts_fragment(&fragment.alignment) && fragment.query_segments == [copy]
        });
        let label = format!("length {length} identity {identity} reverse {reverse}");
        assert_eq!(whole.count(), 1, "{label}");
        let query_length = case.query.len() as u64;
        let decision = rule.decide(query_length, fragments, &[], || false, |_| 0);
        assert!(decision.unwrap().qualifies, "{label}");
    }
}
