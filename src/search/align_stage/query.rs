//! The query as stage 3 uses it: its bases and its anchor tables, built once and shared by all
//! workers.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{AlignStageConfig, AlignStageError};
use crate::align::{Interval, QueryTable};

/// The query bases in upper case and the tables of its k-mers. Built once per query; all
/// workers share it.
pub struct PreparedQuery {
    /// The k-mers that pass the entropy limit, for anchors and chains.
    pub(super) table: QueryTable,
    /// The low-complexity k-mers, left to dense walks; `None` when the query has none.
    pub(super) low: Option<QueryTable>,
    /// Low-complexity query stretches, sorted and disjoint; they count toward no coverage in the
    /// decision. Empty from `new`; the search sets them.
    pub low_complexity: Vec<Interval>,
}

impl PreparedQuery {
    /// Prepares `query`, linear or circular, with the anchor settings of `config`. Bases are
    /// taken in upper case; other bytes such as IUPAC codes stay, as packed sequences keep them.
    pub fn new(
        query: &[u8],
        circular: bool,
        config: &AlignStageConfig,
    ) -> Result<Self, AlignStageError> {
        config.validate()?;
        if query.is_empty() {
            return Err(AlignStageError::EmptyQuery);
        }
        let bases = query.to_ascii_uppercase();
        let (table, low) = QueryTable::split(&bases, circular, config.anchors, config.chains)?;
        Ok(Self {
            table,
            low: (low.max_occurrences() > 0).then_some(low),
            low_complexity: Vec::new(),
        })
    }

    /// Drops the low-complexity k-mers, so that their copies give no dense walks.
    pub fn skip_low_complexity_walks(&mut self) {
        self.low = None;
    }

    /// The query bases, in upper case.
    pub fn bases(&self) -> &[u8] {
        self.table.query()
    }

    /// The number of query bases.
    pub fn length(&self) -> u64 {
        u64::from(self.table.query_length())
    }

    /// Whether the query is circular.
    pub fn circular(&self) -> bool {
        self.table.circular()
    }

    /// Checks that the query was prepared with the anchor settings of `config`.
    pub(super) fn check(&self, config: &AlignStageConfig) -> Result<(), AlignStageError> {
        if self.table.anchors() != config.anchors || self.table.chains() != config.chains {
            return Err(AlignStageError::InvalidConfig(
                "the query was prepared with other anchor settings",
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn queries_upper_case() {
        let config = AlignStageConfig::default();
        let query = PreparedQuery::new(b"acgtRYacgt", true, &config).unwrap();
        assert_eq!(query.bases(), b"ACGTRYACGT");
        assert_eq!(query.length(), 10);
        assert!(query.circular());
        assert!(query.low.is_none());
        assert!(PreparedQuery::new(b"", false, &config).is_err());
        let homopolymer = PreparedQuery::new(&[b'A'; 100], false, &config).unwrap();
        assert!(homopolymer.low.is_some());
        homopolymer.check(&config).unwrap();
        let mut other = config;
        other.anchors.a = 17;
        assert!(homopolymer.check(&other).is_err());
    }
}
