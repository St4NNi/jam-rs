//! The decision rule: which fragments count, and when an entry is reported.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::contigs::counted_contig_fragments;
use super::mosaic::build_mosaic;
use super::{DecideError, Fragment, Mosaic, QueryCoverage};
use crate::align::{Alignment, EditCounts};
use serde::{Deserialize, Serialize};

/// The limits that decide whether an entry is reported.
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
pub struct DecisionRule {
    /// The least query coverage of the mosaic. A share must be at least 0.05, as a large plasmid
    /// can show as a few small pieces.
    pub min_query_coverage: QueryCoverage,
    /// The least identity over all counted fragments together.
    pub min_identity: f64,
    /// The least identity of one fragment to count. 0 lets every fragment count.
    pub min_fragment_identity: f64,
    /// The least number of alignment columns of one fragment to count.
    pub min_fragment_columns: u64,
    /// The least share of a contig's bases that its fragments span for the contig to count.
    /// 0 counts every contig.
    pub min_contig_coverage: f64,
    /// The least identity of the fragments of one contig together for the contig to count.
    /// 0 counts every contig.
    pub min_contig_identity: f64,
}

impl Default for DecisionRule {
    fn default() -> Self {
        Self {
            min_query_coverage: QueryCoverage::Share {
                numerator: 5,
                decimal_places: 1,
            },
            min_identity: 0.95,
            min_fragment_identity: 0.0,
            min_fragment_columns: 100,
            min_contig_coverage: 0.0,
            min_contig_identity: 0.0,
        }
    }
}

/// The result of applying a rule to the fragments of one entry.
#[derive(Clone, Debug, PartialEq)]
pub struct Decision {
    /// The mosaic of the counted fragments.
    pub mosaic: Mosaic,
    /// Base counts summed over the counted fragments.
    pub counts: EditCounts,
    /// Fragments left out by the contig rule or the fragment limits.
    pub dropped_fragments: usize,
    /// Whether the entry is reported.
    pub qualifies: bool,
}

impl Decision {
    /// Matches over columns of the counted fragments; 0 when none counts.
    pub fn identity(&self) -> f64 {
        self.counts.identity()
    }
}

impl DecisionRule {
    /// Checks that the limits make sense.
    pub fn validate(&self) -> Result<(), DecideError> {
        self.min_query_coverage.validate()?;
        let fraction = |value: f64| value.is_finite() && (0.0..=1.0).contains(&value);
        let shares = [
            self.min_identity,
            self.min_fragment_identity,
            self.min_contig_coverage,
            self.min_contig_identity,
        ];
        if !shares.into_iter().all(fraction) {
            return Err(DecideError::InvalidRule(
                "a share or identity is not from 0 to 1",
            ));
        }
        if self.min_fragment_columns == 0 {
            return Err(DecideError::InvalidRule("fragment columns is zero"));
        }
        Ok(())
    }

    /// Whether a fragment with this alignment counts: it needs `min_fragment_columns` columns
    /// and `min_fragment_identity`, compared exactly as `matches >= limit * columns`.
    pub fn counts_fragment(&self, alignment: &Alignment) -> bool {
        let counts = alignment.counts();
        counts.columns() >= self.min_fragment_columns
            && counts.reaches_identity(self.min_fragment_identity)
    }

    /// Builds the mosaic of the counted fragments of one entry and decides. The identity sums
    /// matches and columns over all counted fragments, each copy too, so a conserved contig can
    /// offset a diverged one. `contig_length` gives a contig length for the contig rule.
    pub fn decide(
        &self,
        query_length: u64,
        fragments: Vec<Fragment>,
        contig_length: impl Fn(u32) -> u64,
    ) -> Result<Decision, DecideError> {
        let before = fragments.len();
        let fragments = if self.min_contig_coverage > 0.0 || self.min_contig_identity > 0.0 {
            counted_contig_fragments(
                fragments,
                self.min_contig_coverage,
                self.min_contig_identity,
                contig_length,
            )
        } else {
            fragments
        };
        let counted: Vec<Fragment> = fragments
            .into_iter()
            .filter(|fragment| self.counts_fragment(&fragment.alignment))
            .collect();
        let dropped_fragments = before - counted.len();
        let mosaic = build_mosaic(query_length, &counted)?;
        let mut counts = EditCounts::default();
        let used = mosaic
            .primary
            .iter()
            .map(|selected| &selected.fragment)
            .chain(&mosaic.alternatives);
        for fragment in used {
            counts += fragment.alignment.counts();
        }
        let qualifies = counts.columns() > 0
            && counts.reaches_identity(self.min_identity)
            && self
                .min_query_coverage
                .accepts(mosaic.covered_bases, query_length);
        Ok(Decision {
            mosaic,
            counts,
            dropped_fragments,
            qualifies,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::decide::fragment::tests::fragment;

    fn rule(coverage: &str) -> DecisionRule {
        DecisionRule {
            min_query_coverage: coverage.parse().unwrap(),
            min_identity: 0.8,
            ..DecisionRule::default()
        }
    }

    fn matches_and_columns(decision: &Decision) -> (u64, u64) {
        (decision.counts.matches, decision.counts.columns())
    }

    #[test]
    fn mixed_identity_combined() {
        let conserved = fragment(0, &[(0, 3000)], "1425=150X1425=");
        let diverged = fragment(1, &[(3000, 5000)], "720=280X720=280X");
        let both = rule("0.5")
            .decide(10_000, vec![conserved.clone(), diverged.clone()], |_| 0)
            .unwrap();
        assert!(both.qualifies);
        assert_eq!(matches_and_columns(&both), (2850 + 1440, 5000));
        assert_eq!(both.mosaic.covered_bases, 5000);
        for alone in [conserved.clone(), diverged.clone()] {
            assert!(
                !rule("0.5")
                    .decide(10_000, vec![alone], |_| 0)
                    .unwrap()
                    .qualifies
            );
        }
        // A limit on each fragment drops the 72 percent contig and the coverage falls short.
        let strict = DecisionRule {
            min_fragment_identity: 0.8,
            ..rule("0.5")
        };
        let decision = strict
            .decide(10_000, vec![conserved, diverged], |_| 0)
            .unwrap();
        assert!(!decision.qualifies);
        assert_eq!(decision.dropped_fragments, 1);
    }

    #[test]
    fn identity_below_fails() {
        let low = fragment(1, &[(1000, 2000)], "700=300X");
        let below = fragment(0, &[(0, 1000)], "899=101X");
        let decision = rule("1")
            .decide(2000, vec![below, low.clone()], |_| 0)
            .unwrap();
        assert_eq!(matches_and_columns(&decision), (1599, 2000));
        assert!(!decision.qualifies);
        let exact = fragment(0, &[(0, 1000)], "900=100X");
        assert!(
            rule("1")
                .decide(2000, vec![exact, low], |_| 0)
                .unwrap()
                .qualifies
        );
    }

    #[test]
    fn identity_column_weighted() {
        let short = fragment(1, &[(9000, 10_000)], "990=10X");
        // 75 and 99 percent average 87, but 7,740 of 10,000 columns match.
        let long = fragment(0, &[(0, 9000)], "6750=2250X");
        let decision = rule("1")
            .decide(10_000, vec![long, short.clone()], |_| 0)
            .unwrap();
        assert!(!decision.qualifies);
        assert!((decision.identity() - 0.774).abs() < 1e-12);
        let long = fragment(0, &[(0, 9000)], "7020=1980X");
        assert!(
            rule("1")
                .decide(10_000, vec![long, short], |_| 0)
                .unwrap()
                .qualifies
        );
    }

    #[test]
    fn copies_count_once() {
        let first = fragment(0, &[(0, 5000)], "4500=500X");
        let copy = fragment(1, &[(0, 5000)], "3600=1400X");
        let decision = rule("0.5")
            .decide(10_000, vec![first.clone(), copy, first.clone()], |_| 0)
            .unwrap();
        assert_eq!(matches_and_columns(&decision), (8100, 10_000));
        assert_eq!(decision.mosaic.covered_bases, 5000);
        assert_eq!(decision.mosaic.primary.len(), 1);
        assert_eq!(decision.mosaic.alternatives.len(), 1);
        assert!(decision.qualifies);
        // A 65 percent copy stays in, so the identity falls to 77.5 percent.
        let weak = fragment(1, &[(0, 5000)], "3250=1750X");
        let decision = rule("0.5")
            .decide(10_000, vec![first.clone(), weak.clone()], |_| 0)
            .unwrap();
        assert!(!decision.qualifies);
        assert!((decision.identity() - 0.775).abs() < 1e-12);
        let strict = DecisionRule {
            min_fragment_identity: 0.7,
            ..rule("0.5")
        };
        assert!(
            strict
                .decide(10_000, vec![first, weak], |_| 0)
                .unwrap()
                .qualifies
        );
    }

    #[test]
    fn weak_piece_needed() {
        let conserved = fragment(0, &[(0, 4000)], "3960=40X");
        let divergent = fragment(1, &[(4000, 5000)], "650=350X");
        for coverage in ["0.5", "1"] {
            let decision = rule(coverage)
                .decide(5000, vec![conserved.clone(), divergent.clone()], |_| 0)
                .unwrap();
            assert!(decision.qualifies);
            assert_eq!(decision.dropped_fragments, 0);
            assert_eq!(decision.mosaic.covered_bases, 5000);
            assert!((decision.identity() - 0.922).abs() < 1e-12);
        }
        let strict = DecisionRule {
            min_fragment_identity: 0.7,
            ..rule("0.5")
        };
        assert!(
            !strict
                .decide(10_000, vec![conserved, divergent], |_| 0)
                .unwrap()
                .qualifies
        );
    }

    #[test]
    fn indels_count_columns() {
        let gapped = fragment(0, &[(0, 1000)], "800=50I100=50D50=");
        let decision = rule("1").decide(1000, vec![gapped.clone()], |_| 0).unwrap();
        assert_eq!(matches_and_columns(&decision), (950, 1050));
        // The 50 query-only bases are not covered, and the 50 contig-only bases use no query base.
        assert_eq!(decision.mosaic.covered_bases, 950);
        assert!(!decision.qualifies);
        assert!(
            rule("0.95")
                .decide(1000, vec![gapped], |_| 0)
                .unwrap()
                .qualifies
        );
    }

    #[test]
    fn coverage_limits_exact() {
        for (limit, required) in [("0.50", 501), ("0.75", 751), ("0.90", 901), ("1.00", 1001)] {
            for (bases, expected) in [(required, true), (required - 1, false)] {
                let exact = fragment(0, &[(0, bases)], &format!("{bases}="));
                let decision = rule(limit).decide(1001, vec![exact], |_| 0).unwrap();
                assert_eq!(decision.qualifies, expected, "{limit} {bases}");
            }
        }
        let wrapped = fragment(0, &[(900, 1001), (0, 400)], "501=");
        let decision = rule("0.5").decide(1001, vec![wrapped], |_| 0).unwrap();
        assert!(decision.qualifies);
        assert_eq!(decision.mosaic.covered_bases, 501);
    }

    #[test]
    fn base_count_decides() {
        let piece = fragment(0, &[(0, 800)], "800=");
        assert!(
            rule("800bp")
                .decide(100_000, vec![piece.clone()], |_| 0)
                .unwrap()
                .qualifies
        );
        assert!(
            !rule("801bp")
                .decide(100_000, vec![piece], |_| 0)
                .unwrap()
                .qualifies
        );
    }

    #[test]
    fn fragment_limits() {
        let rule = rule("0.5");
        let counts = |cigar: &str, length| {
            rule.counts_fragment(&fragment(0, &[(0, length)], cigar).alignment)
        };
        assert!(!counts("99=", 99));
        assert!(counts("100=", 100));
        assert!(counts("70=30X", 100));
        assert!(counts("65=35X", 100));
        let strict = DecisionRule {
            min_fragment_identity: 0.7,
            ..rule
        };
        assert!(!strict.counts_fragment(&fragment(0, &[(0, 100)], "69=31X").alignment));
        // Four decimals would print 0.79996 as 0.8000, but the limit is compared exactly.
        let rounded = fragment(0, &[(0, 100_000)], "79996=20004X");
        let decision = rule.decide(100_000, vec![rounded], |_| 0).unwrap();
        assert_eq!(format!("{:.4}", decision.identity()), "0.8000");
        assert!(!decision.qualifies);
    }

    #[test]
    fn contig_rule_identity() {
        // Contig 0 is the query. Contig 1 shares 400 of 10,000 bases at 75 percent.
        let query_copy = fragment(0, &[(0, 3000)], "2970=30X");
        let shared_part = fragment(1, &[(3000, 3400)], "300=100X");
        let lengths = |contig: u32| [3000, 10_000][contig as usize];
        let fragments = || vec![query_copy.clone(), shared_part.clone()];
        let open = rule("0.5");
        let decision = open.decide(3400, fragments(), lengths).unwrap();
        assert_eq!(matches_and_columns(&decision), (3270, 3400));
        let by_contig = DecisionRule {
            min_contig_coverage: 0.5,
            min_identity: 0.98,
            ..open
        };
        let decision = by_contig.decide(3400, fragments(), lengths).unwrap();
        assert!(decision.qualifies);
        assert_eq!(matches_and_columns(&decision), (2970, 3000));
        assert_eq!(decision.dropped_fragments, 1);
        let without = DecisionRule {
            min_contig_coverage: 0.0,
            ..by_contig
        };
        assert!(
            !without
                .decide(3400, fragments(), lengths)
                .unwrap()
                .qualifies
        );
        let by_identity = DecisionRule {
            min_contig_identity: 0.9,
            ..without
        };
        let decision = by_identity.decide(3400, fragments(), lengths).unwrap();
        assert_eq!(decision.dropped_fragments, 1);
    }

    #[test]
    fn rule_limits_checked() {
        assert!(DecisionRule::default().validate().is_ok());
        for coverage in ["0.05", "0.1", "0.4999", "0.5", "0.75", "1", "1bp", "5000bp"] {
            assert!(rule(coverage).validate().is_ok(), "{coverage}");
        }
        for coverage in ["0.0499", "0", "0bp"] {
            assert!(rule(coverage).validate().is_err(), "{coverage}");
        }
        let invalid = [
            DecisionRule {
                min_identity: f64::NAN,
                ..DecisionRule::default()
            },
            DecisionRule {
                min_fragment_identity: 1.5,
                ..DecisionRule::default()
            },
            DecisionRule {
                min_contig_coverage: 1.5,
                ..DecisionRule::default()
            },
            DecisionRule {
                min_contig_identity: -0.1,
                ..DecisionRule::default()
            },
            DecisionRule {
                min_fragment_columns: 0,
                ..DecisionRule::default()
            },
        ];
        assert!(invalid.iter().all(|rule| rule.validate().is_err()));
    }
}
