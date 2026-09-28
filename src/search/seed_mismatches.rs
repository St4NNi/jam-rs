//! Whether a search also looks up every query k-mer one substitution away, so that seeds with
//! one mismatch are hit too.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::decide::DecisionRule;
use crate::loss::{LossBudget, MissModel};
use serde::{Deserialize, Serialize};

/// The mismatches a query k-mer may have against a seed and still hit it.
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize, clap::ValueEnum)]
#[serde(rename_all = "snake_case")]
pub enum SeedMismatches {
    /// Chosen per query from the identity, the loss budget and the query length
    #[default]
    Auto,
    /// A seed is hit only by an equal query k-mer
    #[value(name = "0")]
    Zero,
    /// A seed is also hit by a query k-mer with one mismatch
    #[value(name = "1")]
    One,
}

/// The extra CPU seconds per query k-mer that the lookups cost, summed over threads: the
/// extra lookups and the extra alignment windows, measured on 1,000 assemblies at 80% identity.
const SECONDS_PER_KMER: f64 = 0.001;

/// The most extra CPU seconds per query that `Auto` spends on the lookups.
const MAX_EXTRA_SECONDS: f64 = 10.0;

impl SeedMismatches {
    /// The mismatches (0 or 1) for a query of `length` bases. `Auto` takes 1 when the proven miss
    /// bound of exact seeds for the required coverage exceeds the budget, one mismatch lowers it,
    /// and query k-mers times `SECONDS_PER_KMER` fit `MAX_EXTRA_SECONDS`.
    pub fn resolve(
        self,
        exact: &MissModel,
        one_mismatch: &MissModel,
        rule: &DecisionRule,
        budget: LossBudget,
        length: u64,
        circular: bool,
    ) -> usize {
        match self {
            Self::Zero => return 0,
            Self::One => return 1,
            Self::Auto => {}
        }
        let required = rule.min_query_coverage.required_bases(length).min(length);
        let identity = rule.min_identity;
        let missed = exact.miss_upper_bound(required, identity);
        let lowered = one_mismatch.miss_upper_bound(required, identity) < missed;
        let kmers = if circular {
            length
        } else {
            (length + 1).saturating_sub(exact.k() as u64)
        };
        let small = kmers as f64 * SECONDS_PER_KMER <= MAX_EXTRA_SECONDS;
        usize::from(missed > budget.chance() && lowered && small)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seeds::SeedScheme;

    fn choice(identity: f64, length: u64) -> usize {
        let scheme = SeedScheme::default();
        let exact = MissModel::new(scheme).unwrap();
        let one_mismatch = MissModel::one_mismatch(scheme).unwrap();
        let rule = DecisionRule {
            min_identity: identity,
            ..DecisionRule::default()
        };
        let budget = LossBudget::default();
        SeedMismatches::Auto.resolve(&exact, &one_mismatch, &rule, budget, length, false)
    }

    #[test]
    fn auto_follows_identity() {
        // At 80% the exact seeds miss short pieces often; past the CPU limit the lookups stay off.
        assert_eq!(choice(0.8, 1_000), 1);
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
        let scheme = SeedScheme::default();
        let (exact, one) = (
            MissModel::new(scheme).unwrap(),
            MissModel::one_mismatch(scheme).unwrap(),
        );
        let budget = LossBudget::default().chance();
        for identity in [0.8, 0.85, 0.9, 0.95, 0.98, 0.99] {
            for length in (200u64..=10_000).step_by(97) {
                let required = length.div_ceil(2);
                let bounds = [&exact, &one].map(|model| model.miss_upper_bound(required, identity));
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
        let model = MissModel::new(SeedScheme::default()).unwrap();
        let rule = DecisionRule::default();
        let budget = LossBudget::default();
        let resolve = |mode: SeedMismatches, length| {
            mode.resolve(&model, &model, &rule, budget, length, false)
        };
        assert_eq!(resolve(SeedMismatches::Zero, 1_000), 0);
        assert_eq!(resolve(SeedMismatches::One, 1_000_000), 1);
    }
}
