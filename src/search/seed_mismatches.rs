//! Which lookups a search makes for each query k-mer, so that seeds with mismatches or an indel
//! are hit too.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::decide::DecisionRule;
use crate::loss::{MaxMiss, MissModel};

/// The CPU seconds per query that a search may spend on the lookups, summed over threads. It is
/// fixed: 30 seconds cover a read of 2 kb at 80% identity.
pub const CPU_BUDGET_SECONDS: f64 = 30.0;

/// The extra CPU seconds per query k-mer of levels 1, 2 and 3, summed over threads. Level 1 is
/// measured on 1,000 assemblies at 80% identity, with the extra alignment windows. Level 2 is
/// measured as 5.5 seconds per 1,000 k-mers. Level 3 stays within 10 seconds up to 1,100 k-mers.
const SECONDS_PER_KMER: [f64; 3] = [0.001, 0.0055, 0.009];

/// The lookup level for a query of `length` bases: 0 exact, 1 one mismatch, 2 two mismatches,
/// 3 level 2 plus the k-mers one indel away. Only levels that fit `CPU_BUDGET_SECONDS` count.
/// It takes the cheapest with a bound within `budget`, else the lowest bound; level 3 replaces 2.
pub fn lookup_level(
    models: &[MissModel; 3],
    rule: &DecisionRule,
    budget: MaxMiss,
    length: u64,
    circular: bool,
) -> usize {
    let kmers = if circular {
        length
    } else {
        (length + 1).saturating_sub(models[0].k() as u64)
    };
    let bound = |level: usize| level_bound(models, rule, length, level);
    let fits = |level: usize| {
        level == 0 || kmers as f64 * SECONDS_PER_KMER[level - 1] <= CPU_BUDGET_SECONDS
    };
    let target = (0..4).find(|&level| fits(level) && bound(level) <= budget.chance());
    // Level 3 has the bound of level 2 and also finds indels, so it replaces level 2 where it fits.
    let upgrade = |level: usize| if level == 2 && fits(3) { 3 } else { level };
    if let Some(level) = target {
        return upgrade(level);
    }
    let mut best = 0;
    for level in 1..4 {
        if fits(level) && bound(level) < bound(best) {
            best = level;
        }
    }
    upgrade(best)
}

/// The proven miss bound of lookup level `level` for a query of `length` bases. Level 3 has the
/// bound of level 2.
fn level_bound(models: &[MissModel; 3], rule: &DecisionRule, length: u64, level: usize) -> f64 {
    let required = rule.min_query_coverage.required_bases(length).min(length);
    let bounded = models[0].qualifying_identity(rule.min_identity, required);
    models[level.min(2)].miss_upper_bound(required, bounded)
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
        lookup_level(&models(), &rule, MaxMiss::default(), length, false)
    }

    #[test]
    fn level_follows_length() {
        // At 95% exact seeds meet the target from about 7.2 kb.
        assert_eq!(choice(0.95, 10_000), 0);
        assert_eq!(choice(0.95, 3_000), 1);
        assert_eq!(choice(0.95, 1_000), 3);
        // At 80%, 1,976 k-mers cost about 18 s at level 3.
        assert_eq!(choice(0.8, 2_000), 3);
        assert_eq!(choice(0.8, 8_000), 1);
        assert_eq!(choice(0.8, 10_000), 1);
        // Level 1 would meet the target but costs about 100 s.
        assert_eq!(choice(0.8, 100_000), 0);
        assert_eq!(choice(0.8, 1_000_000), 0);
        // Below the identities where seeds can be hit at all, no lookup helps.
        assert_eq!(choice(0.0, 1_000), 0);
    }

    /// The chosen bound is the lowest among the allowed levels, or within the budget.
    #[test]
    fn level_takes_lowest() {
        let models = models();
        let budget = MaxMiss::default().chance();
        for identity in [0.8, 0.85, 0.9, 0.95, 0.98, 0.99] {
            for length in (200u64..=10_000).step_by(97) {
                let required = length.div_ceil(2);
                let bounds = models.map(|model| model.miss_upper_bound(required, identity));
                let chosen = bounds[choice(identity, length).min(2)];
                assert!(chosen <= bounds[0].max(budget), "{identity} {length}");
            }
        }
    }
}
