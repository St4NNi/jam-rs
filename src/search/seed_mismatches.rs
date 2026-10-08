//! Whether a search also looks up every query k-mer one substitution away, so that seeds with
//! one mismatch are hit too.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::decide::DecisionRule;
use crate::loss::{MaxMiss, MissModel};
use serde::{Deserialize, Serialize};

/// The mismatches a query k-mer may have against a seed and still hit it.
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize, clap::ValueEnum)]
#[serde(rename_all = "snake_case")]
pub enum SeedMismatches {
    /// Chosen per query from the identity, --max-miss and the query length
    #[default]
    Auto,
    /// A seed is hit only by an equal query k-mer
    #[value(name = "0")]
    Zero,
    /// A seed is also hit by a query k-mer with one mismatch
    #[value(name = "1")]
    One,
    /// A seed is also hit by a query k-mer with two mismatches; on queries of up to about
    /// 1.1 kb also by one with one indel and at most one mismatch
    #[value(name = "2")]
    Two,
}

/// The extra CPU seconds per query k-mer that the lookups cost, summed over threads: the
/// extra lookups and the extra alignment windows, measured on 1,000 assemblies at 80% identity.
const SECONDS_PER_KMER: f64 = 0.001;

/// The most extra CPU seconds per query that `Auto` spends on the lookups.
const MAX_EXTRA_SECONDS: f64 = 10.0;

/// The most query k-mers for which `Auto` takes 2: two mismatches cost about 5.5 CPU seconds per
/// 1,000 query k-mers (measured), which keeps 1,800 k-mers within `MAX_EXTRA_SECONDS`.
const MAX_TWO_KMERS: u64 = 1_800;

/// The identity below which `Auto` may take 2.
const TWO_BELOW_IDENTITY: f64 = 0.9;

impl SeedMismatches {
    /// The mismatches (0, 1 or 2) for a query of `length` bases. `Auto` takes 1 when exact seeds
    /// miss above the budget, one mismatch lowers it and the CPU limit holds; 2 instead when also
    /// the identity, `MAX_TWO_KMERS` and the bounds of one and two mismatches allow it.
    pub fn resolve(
        self,
        models: &[MissModel; 3],
        rule: &DecisionRule,
        budget: MaxMiss,
        length: u64,
        circular: bool,
    ) -> usize {
        match self {
            Self::Zero => return 0,
            Self::One => return 1,
            Self::Two => return 2,
            Self::Auto => {}
        }
        let [exact, one_mismatch, two_mismatch] = models;
        let required = rule.min_query_coverage.required_bases(length).min(length);
        let identity = rule.min_identity;
        let bounded = exact.qualifying_identity(identity, required);
        let missed = exact.miss_upper_bound(required, bounded);
        let lowered = one_mismatch.miss_upper_bound(required, bounded) < missed;
        let kmers = if circular {
            length
        } else {
            (length + 1).saturating_sub(exact.k() as u64)
        };
        let small = kmers as f64 * SECONDS_PER_KMER <= MAX_EXTRA_SECONDS;
        if !(missed > budget.chance() && lowered && small) {
            return 0;
        }
        let one = one_mismatch.miss_upper_bound(required, bounded);
        let two = two_mismatch.miss_upper_bound(required, bounded);
        let lower = identity < TWO_BELOW_IDENTITY && kmers <= MAX_TWO_KMERS;
        if lower && one > budget.chance() && two < one {
            2
        } else {
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seeds::SeedScheme;

    fn choice(identity: f64, length: u64) -> usize {
        choice_with(identity, length, MaxMiss::default().chance())
    }

    fn models() -> [MissModel; 3] {
        [0, 1, 2].map(|mismatches| MissModel::new(SeedScheme::default(), mismatches).unwrap())
    }

    fn choice_with(identity: f64, length: u64, budget: f64) -> usize {
        let rule = DecisionRule {
            min_identity: identity,
            ..DecisionRule::default()
        };
        let budget = MaxMiss::new(budget).unwrap();
        let mode = SeedMismatches::Auto;
        mode.resolve(&models(), &rule, budget, length, false)
    }

    #[test]
    fn auto_takes_two() {
        // 1,000 bases hold 976 k-mers; 2,000 hold 1,976, past the limit of 1,800.
        assert_eq!(choice(0.8, 1_000), 2);
        assert_eq!(choice(0.8, 1_824), 2);
        assert_eq!(choice(0.8, 1_826), 1);
        assert_eq!(choice(0.8, 5_000), 1);
        assert_eq!(choice(0.89, 1_000), 2);
        // Not at the default identity, nor from 0.9.
        assert_eq!(choice(0.95, 1_000), 1);
        assert_eq!(choice(0.9, 1_000), 1);
        assert_eq!(choice(0.99, 1_000), 1);
    }

    /// When one mismatch already keeps the bound within the budget, 2 is not taken.
    #[test]
    fn auto_keeps_one() {
        // At 0.89 and 1,000 bases the bound with one mismatch is 0.37, with two 0.078.
        assert_eq!(choice_with(0.89, 1_000, 0.0001), 2);
        assert_eq!(choice_with(0.89, 1_000, 0.1), 2);
        assert_eq!(choice_with(0.89, 1_000, 0.5), 1);
    }

    #[test]
    fn auto_follows_identity() {
        // At 80% the exact seeds miss short pieces often; past the CPU limit the lookups stay off.
        assert_eq!(choice(0.8, 1_000), 2);
        assert_eq!(choice(0.8, 5_000), 1);
        assert_eq!(choice(0.8, 20_000), 0);
        assert_eq!(choice(0.8, 100_000), 0);
        // At 95% the proven bound of exact seeds reaches the budget from 29 blocks of 124 bases,
        // a required piece of 3,596 bases, so a query of 7,192 bases.
        assert_eq!(choice(0.95, 1_000), 1);
        assert_eq!(choice(0.95, 7_190), 1);
        assert_eq!(choice(0.95, 7_192), 0);
        assert_eq!(choice(0.95, 100_000), 0);
        // Below the identities where seeds can be hit at all, no seeding helps.
        assert_eq!(choice(0.0, 1_000), 0);
    }

    /// Where the CPU limit allows one mismatch, auto takes the mode with the lower proven bound,
    /// so that bound fits the budget whenever either mode's bound does.
    #[test]
    fn auto_takes_lower() {
        let models = models();
        let budget = MaxMiss::default().chance();
        for identity in [0.8, 0.85, 0.9, 0.95, 0.98, 0.99] {
            for length in (200u64..=10_000).step_by(97) {
                let required = length.div_ceil(2);
                let bounds = models.map(|model| model.miss_upper_bound(required, identity));
                let chosen = bounds[choice(identity, length)];
                assert!(
                    chosen <= bounds[0].min(bounds[1]).max(budget),
                    "{identity} {length}"
                );
            }
        }
    }

    #[test]
    fn forced_choice_holds() {
        let models = models();
        let rule = DecisionRule::default();
        let budget = MaxMiss::default();
        let resolve =
            |mode: SeedMismatches, length| mode.resolve(&models, &rule, budget, length, false);
        assert_eq!(resolve(SeedMismatches::Zero, 1_000), 0);
        assert_eq!(resolve(SeedMismatches::One, 1_000_000), 1);
        assert_eq!(resolve(SeedMismatches::Two, 1_000_000), 2);
    }
}
