//! The search settings a client sends with a search; threads and memory limits belong to the
//! server.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::ServiceError;
use crate::decide::DecisionRule;
use crate::loss::LossBudget;
use crate::search::{SearchConfig, SearchStage, SeedMismatches};
use serde::{Deserialize, Serialize};

/// The settings of one search that a client may choose.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SearchSettings {
    /// The last stage that runs.
    pub last_stage: SearchStage,
    /// The decision rule of the search.
    pub rule: DecisionRule,
    /// The identity below which a chain's anchors are too sparse for it to be aligned.
    pub min_chain_identity: f64,
    /// The largest accepted chance of missing a qualifying entry.
    pub loss_budget: f64,
    /// Whether query k-mers one substitution away are looked up too.
    pub seed_mismatches: SeedMismatches,
    /// The alignment score of a match.
    pub match_score: i32,
    /// The alignment score of a mismatch.
    pub mismatch_score: i32,
    /// The score added once for each gap.
    pub gap_open_score: i32,
    /// The score added for each gap base.
    pub gap_extend_score: i32,
}

impl SearchSettings {
    /// The settings a client chose in `config`.
    pub fn of(config: &SearchConfig) -> Self {
        let (rule, scores) = (&config.rule, &config.align.band.alignment);
        Self {
            last_stage: config.last_stage,
            rule: *rule,
            min_chain_identity: config.align.min_chain_identity,
            loss_budget: config.loss_budget.chance(),
            seed_mismatches: config.seed_mismatches,
            match_score: scores.match_score,
            mismatch_score: scores.mismatch_score,
            gap_open_score: scores.gap_open_score,
            gap_extend_score: scores.gap_extend_score,
        }
    }

    /// The server's `base` configuration with these settings in place, checked.
    pub fn apply(&self, base: &SearchConfig) -> Result<SearchConfig, ServiceError> {
        let mut config = *base;
        config.last_stage = self.last_stage;
        config.rule = self.rule;
        config.align.min_chain_identity = self.min_chain_identity;
        config.loss_budget = LossBudget::new(self.loss_budget)?;
        config.seed_mismatches = self.seed_mismatches;
        let scores = &mut config.align.band.alignment;
        scores.match_score = self.match_score;
        scores.mismatch_score = self.mismatch_score;
        scores.gap_open_score = self.gap_open_score;
        scores.gap_extend_score = self.gap_extend_score;
        config.validate()?;
        Ok(config)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn settings_keep_limits() {
        let mut client = SearchConfig {
            last_stage: SearchStage::Locate,
            threads: 3,
            ..SearchConfig::default()
        };
        client.rule.min_query_coverage = "0.7".parse().unwrap();
        client.rule.min_identity = 0.9;
        client.loss_budget = LossBudget::new(0.001).unwrap();
        client.seed_mismatches = SeedMismatches::One;
        client.align.band.alignment.mismatch_score = -4;
        let settings = SearchSettings::of(&client);
        let text = serde_json::to_string(&settings).unwrap();
        let read: SearchSettings = serde_json::from_str(&text).unwrap();
        assert_eq!(read, settings);

        let base = SearchConfig {
            threads: 8,
            ..SearchConfig::default()
        };
        let applied = read.apply(&base).unwrap();
        assert_eq!(applied.threads, 8);
        assert_eq!(applied.last_stage, SearchStage::Locate);
        assert_eq!(applied.rule, client.rule);
        assert_eq!(applied.loss_budget, client.loss_budget);
        assert_eq!(applied.seed_mismatches, SeedMismatches::One);
        assert_eq!(applied.align.band.alignment, client.align.band.alignment);
    }

    #[test]
    fn invalid_settings_fail() {
        let mut settings = SearchSettings::of(&SearchConfig::default());
        settings.loss_budget = 0.0;
        assert!(settings.apply(&SearchConfig::default()).is_err());
        let mut settings = SearchSettings::of(&SearchConfig::default());
        settings.rule.min_identity = 2.0;
        assert!(settings.apply(&SearchConfig::default()).is_err());
    }
}
