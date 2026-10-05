//! The decision options that `jam search` and `jam screen` share.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::decide::{DecisionRule, QueryCoverage};
use clap::Args;

/// The options of the decision rule; their defaults are those of `DecisionRule`.
#[derive(Debug, Args)]
pub(super) struct RuleArgs {
    /// Least share of the query that the alignments must cover: a share such as 0.5, or a base
    /// count such as 800bp
    #[arg(short, long, default_value_t = DecisionRule::default().min_query_coverage)]
    coverage: QueryCoverage,
    /// Least identity over the covered query bases
    #[arg(short, long, default_value_t = DecisionRule::default().min_identity)]
    identity: f64,
    /// Least identity of a stretch of an alignment to count toward coverage [default: the
    /// --identity value]
    #[arg(long)]
    stretch_identity: Option<f64>,
    /// Least alignment columns of a stretch to count toward coverage
    #[arg(long, default_value_t = DecisionRule::default().min_fragment_columns)]
    stretch_length: u64,
    /// A contig counts only if the alignments cover at least this share of the contig itself
    /// (0 turns this off)
    #[arg(long, default_value_t = DecisionRule::default().min_contig_coverage)]
    contig_coverage: f64,
    /// A contig counts only if its alignments reach at least this identity (0 turns this off)
    #[arg(long, default_value_t = DecisionRule::default().min_contig_identity)]
    contig_identity: f64,
}

impl RuleArgs {
    /// The decision rule of the options, not yet checked.
    pub(super) fn rule(&self) -> DecisionRule {
        DecisionRule {
            min_query_coverage: self.coverage,
            min_identity: self.identity,
            min_fragment_identity: self.stretch_identity.unwrap_or(self.identity),
            min_fragment_columns: self.stretch_length,
            min_contig_coverage: self.contig_coverage,
            min_contig_identity: self.contig_identity,
        }
    }
}
