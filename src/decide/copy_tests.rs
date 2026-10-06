//! The decision confirms an entry that holds one copy of enough query bases at the identity
//! limit, whatever other fragments the entry holds.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::fragment::tests::fragment;
use super::{DecisionRule, Fragment};
use crate::align::Interval;
use crate::test_support::rng;
use rand::RngExt;
use rand::rngs::StdRng;

/// The CIGAR of `columns` columns with a substitution at each place in `changes`.
fn substitutions(columns: u64, changes: &[u64]) -> String {
    let mut runs: Vec<(char, u64)> = Vec::new();
    for column in 0..columns {
        let operation = if changes.binary_search(&column).is_ok() {
            'X'
        } else {
            '='
        };
        match runs.last_mut() {
            Some(last) if last.0 == operation => last.1 += 1,
            _ => runs.push((operation, 1)),
        }
    }
    runs.iter()
        .map(|(operation, count)| format!("{count}{operation}"))
        .collect()
}

/// A substitution-only fragment on `contig` over query bases `start..start + columns`.
fn copy(contig: u32, start: u64, columns: u64, changes: &[u64]) -> Fragment {
    let cigar = substitutions(columns, changes);
    fragment(contig, &[(start, start + columns)], &cigar)
}

/// 50 substitutions spread over 1,000 columns: identity 0.95 exactly.
fn spread() -> Vec<u64> {
    (0..50).map(|block| block * 20 + 19).collect()
}

fn half_rule() -> DecisionRule {
    DecisionRule {
        min_query_coverage: "0.5".parse().unwrap(),
        ..DecisionRule::default()
    }
}

fn confirmed(
    rule: &DecisionRule,
    length: u64,
    fragments: Vec<Fragment>,
    masked: &[Interval],
) -> bool {
    let decision = rule
        .decide(length, fragments, masked, || false, |_| 0)
        .unwrap();
    decision.qualifies
}

#[test]
fn copy_alone() {
    let rule = half_rule();
    assert!(confirmed(
        &rule,
        2000,
        vec![copy(0, 500, 1000, &spread())],
        &[]
    ));
    // One substitution more is below the limit.
    let mut more = spread();
    more.push(0);
    more.sort_unstable();
    assert!(!confirmed(
        &rule,
        2000,
        vec![copy(0, 500, 1000, &more)],
        &[]
    ));
}

#[test]
fn worse_overlap() {
    let rule = half_rule();
    let changes: Vec<u64> = (0..150).map(|at| at * 10).collect();
    for start in [0, 300, 500, 900] {
        // 0.9 over 1,500 bases, overlapping the copy in different ways.
        let worse = copy(1, start, 1500.min(2000 - start), &changes);
        let fragments = vec![worse, copy(0, 500, 1000, &spread())];
        assert!(confirmed(&rule, 2000, fragments, &[]), "{start}");
    }
}

#[test]
fn lower_elsewhere() {
    // 0.85 fragments elsewhere pull the identity of all fragments together below 0.95.
    let rule = DecisionRule {
        min_fragment_identity: 0.8,
        ..half_rule()
    };
    let changes: Vec<u64> = (0..45).map(|at| at * 20 + 3).collect();
    let mut fragments = vec![copy(0, 0, 1000, &spread())];
    for contig in 1..4 {
        fragments.push(copy(contig, 700 + 300 * contig as u64, 300, &changes));
    }
    assert!(confirmed(&rule, 2000, fragments, &[]));
}

#[test]
fn clustered_substitutions() {
    // Ruzzo-Tompa segments would split this copy into 450= and 500= and lose 50 bases.
    let clustered: Vec<u64> = (450..500).collect();
    let rule = half_rule();
    let alone = vec![copy(0, 0, 1000, &clustered)];
    assert!(confirmed(&rule, 2000, alone, &[]));
    let changes: Vec<u64> = (0..100).map(|at| at * 10).collect();
    let fragments = vec![copy(1, 0, 1000, &changes), copy(0, 0, 1000, &clustered)];
    assert!(confirmed(&rule, 2000, fragments, &[]));
}

#[test]
fn exact_limits() {
    for (limit, changed) in [(0.8, 200), (0.9, 100), (0.95, 50), (0.97, 30), (0.99, 10)] {
        let rule = DecisionRule {
            min_identity: limit,
            min_fragment_identity: limit,
            ..half_rule()
        };
        let clustered: Vec<u64> = (300..300 + changed).collect();
        assert!(confirmed(
            &rule,
            2000,
            vec![copy(0, 0, 1000, &clustered)],
            &[]
        ));
        let more: Vec<u64> = (300..301 + changed).collect();
        assert!(!confirmed(&rule, 2000, vec![copy(0, 0, 1000, &more)], &[]));
    }
}

#[test]
fn masked_better_fragment() {
    // The other fragment ranks first by its exact masked part (masked columns count toward
    // identity, as in the fragment limits), and covers the copy's bases.
    let masked = [Interval::new(0, 3000).unwrap()];
    let other = copy(1, 0, 4000, &(3000..3100).collect::<Vec<_>>());
    let rule = DecisionRule {
        min_query_coverage: "1000bp".parse().unwrap(),
        ..DecisionRule::default()
    };
    let fragments = vec![other, copy(0, 3000, 1000, &spread())];
    let decision = rule
        .decide(4000, fragments, &masked, || false, |_| 0)
        .unwrap();
    assert!(decision.qualifies);
    assert_eq!(decision.mosaic.covered_bases, 1000);
}

#[test]
fn masked_matches_count() {
    // A copy exactly at the limit stays confirmed when the mask holds only its matches.
    let rule = half_rule();
    let changes: Vec<u64> = (0..50).map(|at| 100 + at * 18).collect();
    let masked = [Interval::new(0, 100).unwrap()];
    assert!(confirmed(
        &rule,
        1800,
        vec![copy(0, 0, 1000, &changes)],
        &masked
    ));
}

/// Random substitution places: `count` of `columns`, sorted.
fn random_changes(random: &mut StdRng, columns: u64, count: u64) -> Vec<u64> {
    let mut changes: Vec<u64> = (0..columns).collect();
    for at in 0..count as usize {
        let other = random.random_range(at..columns as usize);
        changes.swap(at, other);
    }
    changes.truncate(count as usize);
    changes.sort_unstable();
    changes
}

#[test]
fn random_entries() {
    let mut random = rng(31);
    let length = 5000;
    let mut checked = 0;
    for round in 0..400 {
        let rule = DecisionRule {
            min_query_coverage: "0.3".parse().unwrap(),
            min_fragment_identity: [0.95, 0.8][round % 2],
            ..DecisionRule::default()
        };
        let columns = random.random_range(1500..3000);
        let start = random.random_range(0..length - columns);
        let count = random.random_range(0..=columns / 20);
        let changes = random_changes(&mut random, columns, count);
        let mut masked = Vec::new();
        let mut at = 0;
        for _ in 0..random.random_range(0..3) {
            let from = random.random_range(at..length);
            let to = random.random_range(from..length.min(from + 400) + 1);
            masked.push(Interval::new(from, to).unwrap());
            at = to + 1;
            if at >= length {
                break;
            }
        }
        // The copy's bases and matches outside the masked bases.
        let inside = |base: u64| {
            masked
                .iter()
                .any(|piece| piece.start <= base && base < piece.end)
        };
        let free = (start..start + columns).filter(|&base| !inside(base));
        let bases = free.clone().count() as u64;
        let changed = free
            .filter(|base| changes.contains(&(base - start)))
            .count() as u64;
        let mut fragments = vec![copy(0, start, columns, &changes)];
        for contig in 1..random.random_range(1..10) {
            let other = random.random_range(100..3000);
            let from = random.random_range(0..length - other);
            let share = random.random_range(0..=other * 3 / 10);
            fragments.push(copy(
                contig,
                from,
                other,
                &random_changes(&mut random, other, share),
            ));
        }
        if bases >= 1500 && (bases - changed) * 20 >= bases * 19 {
            checked += 1;
            let decision = rule
                .decide(length, fragments, &masked, || false, |_| 0)
                .unwrap();
            assert!(decision.qualifies, "round {round}");
        }
    }
    assert!(checked > 100, "{checked}");
}

#[test]
fn split_copy_outranked() {
    // A and B are adjacent on contig 0 at 0.95. C on contig 1 is better and takes part of A's
    // bases first, so A adds 450 matches in 500 bases and all fragments together reach 0.94.
    let rule = DecisionRule {
        min_query_coverage: "2000bp".parse().unwrap(),
        ..DecisionRule::default()
    };
    let a = fragment(0, &[(0, 1000)], "500=50X450=");
    let b = fragment(0, &[(1000, 2000)], "500=50X450=");
    let c = fragment(1, &[(0, 500)], "240=20X240=");
    assert!(confirmed(&rule, 4000, vec![a.clone(), b.clone()], &[]));
    let decision = rule
        .decide(4000, vec![a, b, c], &[], || false, |_| 0)
        .unwrap();
    assert!(decision.qualifies);
    assert_eq!(decision.mosaic.covered_bases, 2000);
    assert!(
        decision
            .mosaic
            .primary
            .iter()
            .all(|part| part.fragment.contig_index == 0)
    );
}

#[test]
fn contig_with_extras() {
    let mut random = rng(47);
    let length = 4000;
    let rule = DecisionRule {
        min_query_coverage: "1500bp".parse().unwrap(),
        ..DecisionRule::default()
    };
    for round in 0..300 {
        // One contig: 2 or 3 adjacent pieces, each at 0.95 or better, covering 1,500 or more.
        let mut start = random.random_range(0..length - 3000);
        let mut own = Vec::new();
        for _ in 0..random.random_range(2..4) {
            let columns = random.random_range(760..1000);
            let count = columns / 20 - random.random_range(0..2);
            let changes = match random.random_range(0..2) {
                0 => random_changes(&mut random, columns, count),
                _ => (columns - count..columns).collect(),
            };
            own.push(copy(0, start, columns, &changes));
            start += columns;
        }
        assert!(confirmed(&rule, length, own.clone(), &[]), "own {round}");
        let first = own[0].alignment.query_interval.start;
        let mut all = own;
        // A better fragment over the start of the copy, as in the case above.
        let columns = random.random_range(300..700);
        let spread = random_changes(&mut random, columns, columns / 25);
        all.push(copy(9, first, columns, &spread));
        for contig in 1..random.random_range(2..8) {
            let columns = random.random_range(100..1500);
            let from = random.random_range(0..length - columns);
            let count = random.random_range(0..=columns / 25);
            let changes = random_changes(&mut random, columns, count);
            all.push(copy(contig, from, columns, &changes));
        }
        assert!(confirmed(&rule, length, all, &[]), "extras {round}");
    }
}
