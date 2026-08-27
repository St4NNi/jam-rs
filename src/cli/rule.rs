//! The decision options that `jam search` and `jam screen` share.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::decide::{DecisionRule, QueryCoverage};
use clap::Args;

/// The options of the decision rule; their defaults are those of `DecisionRule`.
#[derive(Debug, Args)]
pub(super) struct RuleArgs {
    /// Least query coverage: a share such as 0.5, or a base count such as 800bp
    #[arg(short, long, default_value_t = DecisionRule::default().min_query_coverage)]
    coverage: QueryCoverage,
    /// Least identity over all counted fragments together
    #[arg(short, long, default_value_t = DecisionRule::default().min_identity)]
    identity: f64,
    /// Least identity of one fragment to count
    #[arg(long, default_value_t = DecisionRule::default().min_fragment_identity)]
    fragment_identity: f64,
    /// Least alignment columns of one fragment to count
    #[arg(long, default_value_t = DecisionRule::default().min_fragment_columns)]
    fragment_columns: u64,
    /// Contig rule: least share of a contig's bases that must align (0 turns the rule off)
    #[arg(long, default_value_t = DecisionRule::default().min_contig_coverage)]
    contig_coverage: f64,
    /// Contig rule: least identity of the fragments of a contig
    #[arg(long, default_value_t = DecisionRule::default().min_contig_identity)]
    contig_identity: f64,
}

impl RuleArgs {
    /// The decision rule of the options, not yet checked.
    pub(super) fn rule(&self) -> DecisionRule {
        DecisionRule {
            min_query_coverage: self.coverage,
            min_identity: self.identity,
            min_fragment_identity: self.fragment_identity,
            min_fragment_columns: self.fragment_columns,
            min_contig_coverage: self.contig_coverage,
            min_contig_identity: self.contig_identity,
        }
    }
}
