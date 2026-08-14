//! The search settings a client sends with a search; threads and memory limits belong to the
//! server.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::ServiceError;
use crate::decide::QueryCoverage;
use crate::loss::LossBudget;
use crate::report::Stage;
use crate::search::{SearchConfig, SearchStage};
use serde::{Deserialize, Serialize};

/// The settings of one search that a client may choose.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SearchSettings {
    /// The last stage that runs.
    pub level: Stage,
    /// The least query coverage of the mosaic.
    pub min_query_coverage: CoverageSetting,
    /// The least identity over all counted fragments together.
    pub min_identity: f64,
    /// The least identity of one fragment to count.
    pub min_fragment_identity: f64,
    /// The least number of alignment columns of one fragment to count.
    pub min_fragment_columns: u64,
    /// The least share of a contig's bases that its fragments span for the contig to count.
    pub min_contig_coverage: f64,
    /// The least identity of the fragments of one contig together for the contig to count.
    pub min_contig_identity: f64,
    /// The identity below which a chain's anchors are too sparse for it to be aligned.
    pub min_chain_identity: f64,
    /// The largest accepted chance of missing a qualifying entry.
    pub loss_budget: f64,
    /// The alignment score of a match.
    pub match_score: i32,
    /// The alignment score of a mismatch.
    pub mismatch_score: i32,
    /// The score added once for each gap.
    pub gap_open_score: i32,
    /// The score added for each gap base.
    pub gap_extend_score: i32,
}

/// The least query coverage: an exact decimal share of the query, or a number of bases.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CoverageSetting {
    /// A share of the query, `numerator` divided by 10 to the power `decimal_places`.
    Share { numerator: u64, decimal_places: u32 },
    /// A number of query bases.
    Bases(u64),
}

impl SearchSettings {
    /// The settings a client chose in `config`.
    pub fn of(config: &SearchConfig) -> Self {
        let (rule, scores) = (&config.rule, &config.align.band.alignment);
        Self {
            level: config.level.into(),
            min_query_coverage: rule.min_query_coverage.into(),
            min_identity: rule.min_identity,
            min_fragment_identity: rule.min_fragment_identity,
            min_fragment_columns: rule.min_fragment_columns,
            min_contig_coverage: rule.min_contig_coverage,
            min_contig_identity: rule.min_contig_identity,
            min_chain_identity: config.align.min_chain_identity,
            loss_budget: config.loss_budget.chance(),
            match_score: scores.match_score,
            mismatch_score: scores.mismatch_score,
            gap_open_score: scores.gap_open_score,
            gap_extend_score: scores.gap_extend_score,
        }
    }

    /// The server's `base` configuration with these settings in place, checked.
    pub fn apply(&self, base: &SearchConfig) -> Result<SearchConfig, ServiceError> {
        let mut config = *base;
        config.level = match self.level {
            Stage::Sketch => SearchStage::Sketch,
            Stage::Locate => SearchStage::Locate,
            Stage::Align => SearchStage::Align,
        };
        let rule = &mut config.rule;
        rule.min_query_coverage = self.min_query_coverage.into();
        rule.min_identity = self.min_identity;
        rule.min_fragment_identity = self.min_fragment_identity;
        rule.min_fragment_columns = self.min_fragment_columns;
        rule.min_contig_coverage = self.min_contig_coverage;
        rule.min_contig_identity = self.min_contig_identity;
        config.align.min_chain_identity = self.min_chain_identity;
        config.loss_budget = LossBudget::new(self.loss_budget)?;
        let scores = &mut config.align.band.alignment;
        scores.match_score = self.match_score;
        scores.mismatch_score = self.mismatch_score;
        scores.gap_open_score = self.gap_open_score;
        scores.gap_extend_score = self.gap_extend_score;
        config.validate()?;
        Ok(config)
    }
}

impl From<QueryCoverage> for CoverageSetting {
    fn from(coverage: QueryCoverage) -> Self {
        match coverage {
            QueryCoverage::Share {
                numerator,
                decimal_places,
            } => Self::Share {
                numerator,
                decimal_places,
            },
            QueryCoverage::Bases(bases) => Self::Bases(bases),
        }
    }
}

impl From<CoverageSetting> for QueryCoverage {
    fn from(coverage: CoverageSetting) -> Self {
        match coverage {
            CoverageSetting::Share {
                numerator,
                decimal_places,
            } => Self::Share {
                numerator,
                decimal_places,
            },
            CoverageSetting::Bases(bases) => Self::Bases(bases),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn settings_keep_limits() {
        let mut client = SearchConfig {
            level: SearchStage::Locate,
            threads: 3,
            ..SearchConfig::default()
        };
        client.rule.min_query_coverage = "0.7".parse().unwrap();
        client.rule.min_identity = 0.9;
        client.loss_budget = LossBudget::new(0.001).unwrap();
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
        assert_eq!(applied.level, SearchStage::Locate);
        assert_eq!(applied.rule, client.rule);
        assert_eq!(applied.loss_budget, client.loss_budget);
        assert_eq!(applied.align.band.alignment, client.align.band.alignment);
    }

    #[test]
    fn invalid_settings_fail() {
        let mut settings = SearchSettings::of(&SearchConfig::default());
        settings.loss_budget = 0.0;
        assert!(settings.apply(&SearchConfig::default()).is_err());
        let mut settings = SearchSettings::of(&SearchConfig::default());
        settings.min_identity = 2.0;
        assert!(settings.apply(&SearchConfig::default()).is_err());
    }
}
