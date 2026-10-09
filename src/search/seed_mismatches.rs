//! Which lookups a search makes for each query k-mer, so that seeds with mismatches or an indel
//! are hit too.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::decide::DecisionRule;
use crate::loss::{MaxMiss, MissModel};

/// The extra CPU seconds per query k-mer of levels 1, 2 and 3, summed over threads. Level 1 is
/// measured on 1,000 assemblies at 80% identity, with the extra alignment windows. Level 2 is
/// measured as 5.5 seconds per 1,000 k-mers. Level 3 stays within 10 seconds up to 1,100 k-mers.
const SECONDS_PER_KMER: [f64; 3] = [0.001, 0.0055, 0.009];

/// The extra CPU seconds per query above which a level is stepped down. Only low identities reach it.
const SAFETY_SECONDS: f64 = 60.0;

/// The seed hits that a whole copy is expected to give before a level is enough. A copy then
/// shows a seed hit with a chance of about 95%.
const EXPECTED_HITS: f64 = 3.0;

/// The most query k-mers for which level 2 also looks up the k-mers one indel away.
const MAX_INDEL_KMERS: u64 = 1_100;

/// The chance that a k-mer holds at most `mismatches` substitutions at identity `identity`.
fn at_most(k: usize, mismatches: usize, identity: f64) -> f64 {
    let (mut choose, mut sum) = (1.0, 0.0);
    for j in 0..=mismatches {
        sum += choose * (1.0 - identity).powi(j as i32) * identity.powi((k - j) as i32);
        choose *= (k - j) as f64 / (j + 1) as f64;
    }
    sum
}

/// The lookup level of a query: 0 exact, 1 or 2 mismatches, 3 level 2 plus indels (up to
/// `MAX_INDEL_KMERS`). It is the cheapest level where a whole copy at the rule's identity is
/// expected to show `EXPECTED_HITS` seed hits, else 2; over `SAFETY_SECONDS` it steps down.
pub fn lookup_level(
    models: &[MissModel; 3],
    rule: &DecisionRule,
    length: u64,
    circular: bool,
) -> usize {
    let (k, w) = (models[0].k(), models[0].w());
    let kmers = if circular {
        length
    } else {
        (length + 1).saturating_sub(k as u64)
    };
    let per_window = 2.0 / (w as f64 + 1.0);
    let hits = |level: usize| kmers as f64 * per_window * at_most(k, level, rule.min_identity);
    let mut level = (0..2).find(|&r| hits(r) >= EXPECTED_HITS).unwrap_or(2);
    if level == 2 && kmers <= MAX_INDEL_KMERS {
        level = 3;
    }
    while level > 0 && kmers as f64 * SECONDS_PER_KMER[level - 1] > SAFETY_SECONDS {
        level -= 1;
    }
    level
}

/// The proven miss bound of lookup level `level` for a query of `length` bases. Level 3 has the
/// bound of level 2.
fn level_bound(models: &[MissModel; 3], rule: &DecisionRule, length: u64, level: usize) -> f64 {
    let required = rule.min_query_coverage.required_bases(length).min(length);
    let bounded = models[0].qualifying_identity(rule.min_identity, required);
    models[level.min(2)].miss_upper_bound(required, bounded)
}

/// The longest query length checked by `target_length`.
const LONGEST_QUERY: f64 = 1e9;

/// The smallest length from which every longer query reaches `--max-miss`, or `None` up to 1e9.
/// The bound is not monotone in the length, as the level changes with it. So lengths are
/// sampled 1% apart, and the first reaching length after the last failing sample is bisected.
pub fn target_length(models: &[MissModel; 3], rule: &DecisionRule, budget: MaxMiss) -> Option<u64> {
    let reaches = |length: u64| {
        let level = lookup_level(models, rule, length, false);
        level_bound(models, rule, length, level) <= budget.chance()
    };
    let first = models[0].k() as u64;
    let lengths = std::iter::successors(Some(first as f64), |&n| Some((n * 1.01).max(n + 1.0)))
        .take_while(|&n| n <= LONGEST_QUERY)
        .map(|n| n as u64);
    let mut failing = None;
    for length in lengths {
        if !reaches(length) {
            failing = Some(length);
        }
    }
    let failing = failing?;
    if failing as f64 * 1.01 > LONGEST_QUERY {
        return None;
    }
    // The samples above `failing` reach; bisect (failing, failing * 1.01 + 1].
    let (mut low, mut high) = (failing, (failing as f64 * 1.01) as u64 + 1);
    while high - low > 1 {
        let middle = low + (high - low) / 2;
        if reaches(middle) {
            high = middle;
        } else {
            low = middle;
        }
    }
    Some(high)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seeds::SeedScheme;

    fn models() -> [MissModel; 3] {
        [0, 1, 2].map(|mismatches| MissModel::new(SeedScheme::default(), mismatches).unwrap())
    }

    fn choice(identity: f64, length: u64) -> usize {
        let rule = DecisionRule {
            min_identity: identity,
            ..DecisionRule::default()
        };
        lookup_level(&models(), &rule, length, false)
    }

    #[test]
    fn level_follows_length() {
        // Level 3 is level 2 with the indel lookups, up to 1,100 k-mers.
        for (identity, length, level) in [
            (0.95, 250, 3),
            (0.95, 500, 1),
            (0.95, 600, 0),
            (0.95, 10_000, 0),
            (0.9, 250, 3),
            (0.9, 500, 3),
            (0.9, 1_000, 1),
            (0.9, 2_200, 0),
            (0.85, 1_000, 3),
            (0.85, 2_000, 1),
            (0.85, 8_900, 0),
            (0.8, 1_500, 2),
            (0.8, 5_000, 2),
            (0.8, 6_000, 1),
            (0.8, 40_000, 1),
            (0.8, 41_000, 0),
            (0.8, 1_000_000, 0),
        ] {
            assert_eq!(choice(identity, length), level, "{identity} {length}");
        }
    }

    /// The safety limit steps a level down at low identity: 0.6 wants level 2 for 12 kb.
    #[test]
    fn safety_steps_down() {
        assert_eq!(choice(0.6, 12_000), 1);
        assert_eq!(choice(0.6, 70_000), 0);
    }

    #[test]
    fn target_length_found() {
        let models = models();
        let from = |identity: f64| {
            let rule = DecisionRule {
                min_identity: identity,
                ..DecisionRule::default()
            };
            target_length(&models, &rule, MaxMiss::default()).unwrap() as f64
        };
        let near = |at: f64, expected: f64| (at / expected - 1.0).abs() < 0.05;
        assert!(near(from(0.95), 7_200.0));
        assert!(near(from(0.9), 31_000.0));
        assert!(near(from(0.85), 132_000.0));
        assert!(near(from(0.8), 604_000.0));
    }
}
