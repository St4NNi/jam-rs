//! Which lookups a search makes for each query k-mer, so that seeds with mismatches or an indel
//! are hit too.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::decide::DecisionRule;
use crate::loss::{MaxMiss, MissModel};

/// The seed hits that a qualifying copy is expected to give on the required bases before a plan
/// is enough. A copy then shows a seed hit with a chance of about 99.97%.
const EXPECTED_HITS: f64 = 8.0;

/// The seed hits the plan by hits asks for from `HIGH_IDENTITY` up, instead of `EXPECTED_HITS`.
/// Chosen by a pre-registered experiment on separate tuning queries: at 95% to 100% identity 4
/// hits lost no truth pair on tuning and cost a quarter of the CPU of 8 on genes; the few test
/// losses were mosaics of pieces under 124 bases at the coverage limit. Below 90% 8 hits keep
/// partial divergent copies that 4 loses.
const HIGH_IDENTITY_HITS: f64 = 4.0;

/// The identity from which the plan by hits asks for `HIGH_IDENTITY_HITS`.
const HIGH_IDENTITY: f64 = 0.9;

/// The most query k-mers for which the k-mers with level-2 lookups also look up the k-mers one
/// indel away.
const MAX_INDEL_KMERS: u64 = 1_100;

/// The identity and the share of required bases at which `full_two_limit` is set.
const LIMIT_IDENTITY: f64 = 0.8;
const LIMIT_COVERAGE: f64 = 0.5;

/// Query positions per stratum of the sampled lookups.
const STRATUM: u64 = 32;

/// Salts of the hashes that pick the sampled positions.
const STRATUM_SALT: u64 = 0x5bd1_e995_9e37_79b9;
const POSITION_SALT: u64 = 0x2545_f491_4f6c_dd1d;

/// The lookups of one query. Every query k-mer gets the lookups of level `full`; about `sampled`
/// of the `starts` query starts, spread by `sampled_in`, also get those of level `full + 1`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LookupPlan {
    /// The level of every query k-mer: 0 exact, 1 or 2 mismatches.
    pub full: usize,
    /// The query starts expected to get the lookups of level `full + 1`.
    pub sampled: u64,
    /// The query starts: the query length on a circular query, else the length less k - 1.
    pub starts: u64,
    /// Whether the k-mers with level-2 lookups also look up the k-mers one indel away.
    pub indels: bool,
    /// Whether the proven bound of level `full` meets the miss target.
    pub proven: bool,
    /// The seed hits that a copy at the identity limit is expected to give on the required bases.
    pub expected_hits: f64,
}

impl LookupPlan {
    /// The plan of level `level` on every k-mer, 3 being level 2 with the indel lookups.
    pub fn full(level: usize) -> Self {
        Self {
            full: level.min(2),
            sampled: 0,
            starts: 0,
            indels: level >= 3,
            proven: false,
            expected_hits: 0.0,
        }
    }

    /// The highest level of the plan: 0 exact, 1 or 2 mismatches, 3 also one indel.
    pub fn level(&self) -> usize {
        match self.indels {
            true => 3,
            false => self.full + usize::from(self.sampled > 0),
        }
    }

    /// The query positions of stratum `stratum` (positions 32 s to 32 s + 31, as bits) that get
    /// the lookups of level `full + 1`. With x = 32 times the sampled share, the stratum takes
    /// floor(x) positions, and one more with chance x - floor(x), by a hash of the stratum. They
    /// are those with the lowest hashes of their position. So the choice is the same on every
    /// run, follows no period of the query, and a larger share picks a superset.
    pub fn sampled_in(&self, stratum: u64) -> u32 {
        if self.sampled == 0 {
            return 0;
        }
        if self.sampled >= self.starts {
            return u32::MAX;
        }
        let target = self.sampled as f64 / self.starts as f64 * STRATUM as f64;
        let extra = unit(mix(stratum ^ STRATUM_SALT)) < target.fract();
        let count = target as usize + usize::from(extra);
        let mut order: [(u64, u32); STRATUM as usize] = std::array::from_fn(|at| {
            (
                mix((stratum * STRATUM + at as u64) ^ POSITION_SALT),
                at as u32,
            )
        });
        order.sort_unstable();
        order[..count]
            .iter()
            .fold(0, |bits, &(_, at)| bits | 1 << at)
    }

    /// Whether query position `position` gets the lookups of level `full + 1`.
    pub fn is_sampled(&self, position: u64) -> bool {
        self.sampled_in(position / STRATUM) >> (position % STRATUM) & 1 == 1
    }
}

/// The splitmix64 finalizer.
fn mix(value: u64) -> u64 {
    let mut z = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// A hash as a number from 0 up to 1.
fn unit(hash: u64) -> f64 {
    (hash >> 11) as f64 / (1u64 << 53) as f64
}

/// The chance that a k-mer holds at most `mismatches` substitutions at identity `identity`.
fn at_most(k: usize, mismatches: usize, identity: f64) -> f64 {
    let (mut choose, mut sum) = (1.0, 0.0);
    for j in 0..=mismatches {
        sum += choose * (1.0 - identity).powi(j as i32) * identity.powi((k - j) as i32);
        choose *= (k - j) as f64 / (j + 1) as f64;
    }
    sum
}

/// The most query k-mers on which a full level 2 counts as cheap: those at which the plan by
/// hits needs level 2 on every k-mer at 80% identity and coverage 0.5,
/// 8 / (0.5 * 2 / (w + 1) * P(at most 2 substitutions at 80%)), about 8,226 with k 25 and w 100.
/// From 80% identity and coverage 0.5 up, the plan by hits never holds more level-2 k-mers.
pub fn full_two_limit(k: usize, w: usize) -> u64 {
    let per_kmer = 2.0 / (w as f64 + 1.0) * at_most(k, 2, LIMIT_IDENTITY);
    (EXPECTED_HITS / (LIMIT_COVERAGE * per_kmer)).ceil() as u64
}

/// The lookups of a query. First the proof step: the cheapest full level whose proven bound
/// meets `budget`, among levels 0 and 1 and level 2 up to `full_two_limit` k-mers. Otherwise the
/// plan by hits: the cheapest plan at which a copy at the identity limit is expected to give
/// `EXPECTED_HITS` seed hits on the required bases, adding level 1, then level 2, on just enough
/// sampled k-mers (`LookupPlan::sampled_in`); full level 2 if even that is not enough. The k-mers with level-2
/// lookups also get indel lookups on queries of at most `MAX_INDEL_KMERS` k-mers.
pub fn lookup_plan(
    models: &[MissModel; 3],
    rule: &DecisionRule,
    budget: MaxMiss,
    length: u64,
    circular: bool,
) -> LookupPlan {
    let (k, w) = (models[0].k(), models[0].w());
    let starts = if circular {
        length
    } else {
        (length + 1).saturating_sub(k as u64)
    };
    let required = rule.min_query_coverage.required_bases(length).min(length);
    let per_start = required as f64 * 2.0 / (w as f64 + 1.0);
    let chance = |r: usize| at_most(k, r, rule.min_identity);
    let plan = |full: usize, sampled: u64, proven: bool| {
        let added = match sampled {
            0 => 0.0,
            _ => sampled as f64 / starts as f64 * (chance(full + 1) - chance(full)),
        };
        LookupPlan {
            full,
            sampled,
            starts,
            indels: (full == 2 || sampled > 0 && full == 1) && starts <= MAX_INDEL_KMERS,
            proven,
            expected_hits: per_start * (chance(full) + added),
        }
    };
    let cheap = |r: usize| r < 2 || starts <= full_two_limit(k, w);
    let proven =
        (0..3).find(|&r| cheap(r) && level_bound(models, rule, length, r) <= budget.chance());
    if let Some(r) = proven {
        return plan(r, 0, true);
    }
    let target = match rule.min_identity >= HIGH_IDENTITY {
        true => HIGH_IDENTITY_HITS,
        false => EXPECTED_HITS,
    };
    let needed = target / per_start;
    match (0..3).find(|&r| chance(r) >= needed) {
        Some(0) => plan(0, 0, false),
        Some(r) => {
            let share = (needed - chance(r - 1)) / (chance(r) - chance(r - 1));
            let sampled = ((share * starts as f64).ceil() as u64).clamp(1, starts);
            // Level 1 is cheap (76 lookups per k-mer), so it goes on every k-mer; only level 2
            // is sampled.
            match sampled == starts || r == 1 {
                true => plan(r, 0, false),
                false => plan(r - 1, sampled, false),
            }
        }
        None => plan(2, 0, false),
    }
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

/// The smallest length from which the proof step meets `--max-miss` on every longer query, or
/// `None` up to 1e9. It is not monotone in the length, as level 2 stops being cheap. So lengths
/// are sampled 1% apart, and the first reaching length after the last failing sample is bisected.
pub fn target_length(models: &[MissModel; 3], rule: &DecisionRule, budget: MaxMiss) -> Option<u64> {
    let reaches = |length: u64| lookup_plan(models, rule, budget, length, false).proven;
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
    use crate::decide::QueryCoverage;
    use crate::seeds::SeedScheme;

    fn models() -> [MissModel; 3] {
        [0, 1, 2].map(|mismatches| MissModel::new(SeedScheme::default(), mismatches).unwrap())
    }

    fn rule(identity: f64) -> DecisionRule {
        DecisionRule {
            min_identity: identity,
            ..DecisionRule::default()
        }
    }

    fn plan(identity: f64, length: u64) -> LookupPlan {
        lookup_plan(
            &models(),
            &rule(identity),
            MaxMiss::default(),
            length,
            false,
        )
    }

    /// Full level, sampled share in percent (rounded), indels and proof, by identity and length.
    #[test]
    fn plan_follows_length() {
        for (identity, length, full, share, indels, proven) in [
            (0.95, 250, 2, 0, true, false),
            (0.95, 1_000, 1, 0, false, false),
            (0.95, 1_300, 2, 0, false, true),
            (0.95, 3_000, 1, 0, false, true),
            (0.95, 10_000, 0, 0, false, true),
            (0.9, 1_000, 1, 50, true, false),
            (0.9, 2_000, 1, 0, false, false),
            (0.9, 3_000, 1, 0, false, false),
            (0.9, 5_000, 2, 0, false, true),
            (0.9, 50_000, 0, 0, false, true),
            (0.85, 10_000, 1, 0, false, false),
            (0.85, 30_000, 1, 0, false, true),
            (0.8, 5_000, 2, 0, false, false),
            (0.8, 10_000, 1, 75, false, false),
            (0.8, 50_000, 1, 0, false, false),
            (0.8, 100_000, 1, 0, false, true),
            (0.8, 1_000_000, 0, 0, false, true),
        ] {
            let plan = plan(identity, length);
            let at = format!("{identity} {length}");
            assert_eq!(plan.full, full, "{at}");
            let percent = (100.0 * plan.sampled as f64 / plan.starts as f64).round();
            assert_eq!(percent as u64, share, "{at}");
            assert_eq!((plan.indels, plan.proven), (indels, proven), "{at}");
            let target = match identity >= HIGH_IDENTITY {
                true => HIGH_IDENTITY_HITS,
                false => EXPECTED_HITS,
            };
            if !proven && plan.full < 2 {
                assert!(plan.expected_hits >= target - 0.01, "{at}");
            }
        }
    }

    /// Level 2 proves 0.85 up to `full_two_limit` k-mers, one more is not cheap; indels stop
    /// after 1,100 k-mers.
    #[test]
    fn plan_boundaries() {
        let limit = full_two_limit(25, 100);
        assert_eq!(limit, 8_226);
        let at = |identity, kmers: u64| plan(identity, kmers + 24);
        assert!(at(0.85, limit).proven && at(0.85, limit).full == 2);
        assert!(!at(0.85, limit + 1).proven);
        assert!(at(0.8, 1_100).indels && !at(0.8, 1_101).indels);
        // At 100% identity a copy of 124 bases misses no seed: level 0 from 247 bases.
        assert!(plan(1.0, 247).proven && plan(1.0, 247).full == 0);
        assert!(!plan(1.0, 246).proven);
        assert!(plan(0.99, 1_000).proven);
    }

    /// The plan follows the coverage rule: a small share or a fixed number of bases needs more
    /// lookups, and when even level 2 gives no 8 hits it is full level 2.
    #[test]
    fn plan_follows_coverage() {
        let with = |coverage: QueryCoverage, length: u64| {
            let rule = DecisionRule {
                min_query_coverage: coverage,
                ..rule(0.8)
            };
            lookup_plan(&models(), &rule, MaxMiss::default(), length, false)
        };
        let bases = with(QueryCoverage::Bases(1_000), 300_000);
        assert_eq!((bases.full, bases.sampled, bases.proven), (2, 0, false));
        assert!(bases.expected_hits < 2.0);
        let share = QueryCoverage::Share {
            numerator: 5,
            decimal_places: 2,
        };
        let small = with(share, 100_000);
        assert!(small.level() > plan(0.8, 100_000).level());
        // A very small target is not met at 3 kb; the plan by hits applies.
        let tight = lookup_plan(
            &models(),
            &rule(0.95),
            MaxMiss::new(1e-12).unwrap(),
            3_000,
            false,
        );
        assert!(!tight.proven && tight.full == 0 && tight.expected_hits >= EXPECTED_HITS);
    }

    #[test]
    fn target_length_found() {
        let models = models();
        let from = |identity: f64| target_length(&models, &rule(identity), MaxMiss::default());
        assert_eq!(from(0.95), Some(1_239));
        assert_eq!(from(0.9), Some(3_223));
        assert_eq!(from(0.85), Some(23_559));
        assert_eq!(from(0.8), Some(82_583));
    }

    fn sampled(plan: &LookupPlan, positions: u64) -> Vec<u64> {
        (0..positions).filter(|&at| plan.is_sampled(at)).collect()
    }

    /// The same plan picks the same positions, about the share of them, and a larger share
    /// picks a superset.
    #[test]
    fn sampling_nested() {
        let mut plan = plan(0.8, 10_000);
        let first = sampled(&plan, plan.starts);
        assert_eq!(first, sampled(&plan, plan.starts));
        let share = first.len() as f64 / plan.starts as f64;
        assert!((share / (plan.sampled as f64 / plan.starts as f64) - 1.0).abs() < 0.02);
        plan.sampled += plan.starts / 10;
        let more = sampled(&plan, plan.starts);
        assert!(first.iter().all(|at| more.binary_search(at).is_ok()));
        plan.sampled = plan.starts;
        assert_eq!(sampled(&plan, 100), (0..100).collect::<Vec<_>>());
    }

    /// Useful seeds every d positions from offset o, for every period and offset, and their
    /// mirror on the other strand: each layout keeps close to the share of sampled positions.
    /// A lattice of every second position kept none of them at period 2 and offset 0.
    #[test]
    fn sampling_periodic() {
        let mut plan = plan(0.8, 50_000);
        for share in [0.25, 0.5, 0.75] {
            plan.sampled = (share * plan.starts as f64) as u64;
            for period in [2, 3, 4, 7, 16, 32, 64, 101] {
                for offset in 0..period {
                    let useful: Vec<u64> = (offset..plan.starts).step_by(period as usize).collect();
                    let hit = |at: &u64| plan.is_sampled(*at);
                    let forward = useful.iter().filter(|at| hit(at)).count();
                    let mirror = useful
                        .iter()
                        .filter(|&&at| hit(&(plan.starts - 1 - at)))
                        .count();
                    // Within four standard deviations of independent picks.
                    let spread = 4.0 * (share * (1.0 - share) / useful.len() as f64).sqrt();
                    for count in [forward, mirror] {
                        let got = count as f64 / useful.len() as f64;
                        assert!(
                            (got - share).abs() < spread,
                            "{share} {period} {offset} {got}"
                        );
                    }
                }
            }
        }
    }
}
