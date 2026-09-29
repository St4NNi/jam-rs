//! The decision rule: which fragments count, and when an entry is reported.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::contigs::counted_contig_fragments;
use super::mosaic::build_mosaic;
use super::segments::{RunOffsets, maximal_segments, slice_fragment};
use super::{DecideError, Fragment, Mosaic, QueryCoverage};
use crate::align::{Alignment, EditCounts};
use serde::{Deserialize, Serialize};

/// The limits that decide whether an entry is reported.
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
pub struct DecisionRule {
    /// The least query coverage of the mosaic. A share must be at least 0.05, as a large plasmid
    /// can show as a few small pieces.
    pub min_query_coverage: QueryCoverage,
    /// The least identity over the query bases the mosaic covers, each fragment counted only
    /// over the query bases it adds.
    pub min_identity: f64,
    /// The least identity of a stretch of a fragment to count. It is `min_identity` by default,
    /// as in LexicMap and megablast. A fragment counts over its maximal segments at this
    /// identity. 0 lets every whole fragment count.
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
            min_fragment_identity: 0.95,
            min_fragment_columns: 40,
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
    /// Base counts over the query bases the mosaic covers, as in `Mosaic::counts`.
    pub counts: EditCounts,
    /// Fragments left out by the contig rule or the fragment limits.
    pub dropped_fragments: usize,
    /// Whether the entry is reported.
    pub qualifies: bool,
}

impl Decision {
    /// Matches over columns of the covered query bases; 0 when none counts.
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

    /// Adds the counted stretches of a fragment to `output`: the whole fragment when it meets the
    /// fragment limits, else its maximal segments at the fragment identity with enough columns.
    fn counted_parts(
        &self,
        fragment: &Fragment,
        output: &mut Vec<Fragment>,
    ) -> Result<(), DecideError> {
        if self.counts_fragment(&fragment.alignment) {
            output.push(fragment.clone());
            return Ok(());
        }
        let runs = &fragment.alignment.runs;
        let offsets = RunOffsets::new(runs, self.min_fragment_identity);
        for (first, last) in maximal_segments(runs, &offsets) {
            let part = slice_fragment(fragment, &offsets, first, last)?;
            if part.alignment.counts().columns() >= self.min_fragment_columns {
                output.push(part);
            }
        }
        Ok(())
    }

    /// Builds the mosaic of the counted fragments of one entry and decides. Each primary fragment
    /// counts only over the query bases it adds, so a copy that adds few bases does not lower the
    /// identity. `contig_length` gives a contig length for the contig rule.
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
        let mut counted = Vec::new();
        let mut kept = 0;
        for fragment in &fragments {
            let found = counted.len();
            self.counted_parts(fragment, &mut counted)?;
            kept += usize::from(counted.len() > found);
        }
        let dropped_fragments = before - kept;
        let mosaic = build_mosaic(query_length, &counted, self.min_identity)?;
        let counts = mosaic.counts;
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
            min_fragment_identity: 0.0,
            ..DecisionRule::default()
        }
    }

    fn matches_and_columns(decision: &Decision) -> (u64, u64) {
        (decision.counts.matches, decision.counts.columns())
    }

    #[test]
    fn mixed_identity_combined() {
        let conserved = fragment(0, &[(0, 3000)], "1425=150X1425=");
        let diverged = fragment(1, &[(3000, 5000)], &"36=14X".repeat(40));
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

    /// `count` blocks of `matches` matches, one mismatch between neighbours.
    fn spread(matches: u32, count: usize) -> String {
        vec![format!("{matches}="); count].join("1X")
    }

    #[test]
    fn default_fragment_limit() {
        // 960 matches in 999 columns, and 850 matches in 1000 columns over the other half.
        let good = fragment(0, &[(0, 999)], &spread(24, 40));
        let weak = fragment(1, &[(1000, 1902)], &format!("{}3X", spread(17, 50)));
        let default = DecisionRule {
            min_query_coverage: "0.4".parse().unwrap(),
            ..DecisionRule::default()
        };
        assert_eq!(default.min_fragment_identity, default.min_identity);
        let decision = default
            .decide(2000, vec![good.clone(), weak.clone()], |_| 0)
            .unwrap();
        assert!(decision.qualifies);
        assert_eq!(decision.dropped_fragments, 1);
        assert_eq!(decision.mosaic.covered_bases, 999);
        assert!(decision.identity() >= 0.95);
        let more = DecisionRule {
            min_query_coverage: "0.6".parse().unwrap(),
            ..default
        };
        assert!(
            !more
                .decide(2000, vec![good, weak], |_| 0)
                .unwrap()
                .qualifies
        );
    }

    #[test]
    fn lower_fragment_limit() {
        let good = fragment(0, &[(0, 999)], &spread(24, 40));
        let weak = fragment(1, &[(1000, 1899)], &spread(17, 50));
        let lower = DecisionRule {
            min_query_coverage: "0.9".parse().unwrap(),
            min_identity: 0.9,
            min_fragment_identity: 0.8,
            ..DecisionRule::default()
        };
        let decision = lower.decide(2000, vec![good, weak], |_| 0).unwrap();
        assert_eq!(decision.dropped_fragments, 0);
        assert_eq!(decision.mosaic.covered_bases, 999 + 899);
        assert!(decision.qualifies);
    }

    #[test]
    fn dense_stretch_counts() {
        // 1,280 matches in 1,319 columns (0.97), then 480 matches in 539 columns (0.89): the
        // fragment as a whole is below 0.95, its first stretch is not.
        let cigar = format!("{}1X{}", spread(32, 40), spread(8, 60));
        let mixed = fragment(0, &[(0, 1859)], &cigar);
        assert!(!DecisionRule::default().counts_fragment(&mixed.alignment));
        let decision = DecisionRule {
            min_query_coverage: "0.5".parse().unwrap(),
            ..DecisionRule::default()
        }
        .decide(1859, vec![mixed], |_| 0)
        .unwrap();
        assert!(decision.qualifies);
        assert_eq!(decision.dropped_fragments, 0);
        assert_eq!(decision.mosaic.covered_bases, 1319);
        assert_eq!(matches_and_columns(&decision), (1280, 1319));
    }

    #[test]
    fn whole_fragment_kept() {
        // 539 of 542 columns match (99.4 percent) over exactly half of the query.
        let fragment = fragment(0, &[(0, 542)], "39=3X500=");
        let decision = DecisionRule::default()
            .decide(1084, vec![fragment], |_| 0)
            .unwrap();
        assert!(decision.qualifies);
        assert_eq!(decision.mosaic.covered_bases, 542);
        assert_eq!(matches_and_columns(&decision), (539, 542));
    }

    /// 19 matches and 1 mismatch, `blocks` times, then 20 matches.
    fn spread_mismatches(blocks: usize) -> String {
        format!("{}20=", "19=1X".repeat(blocks))
    }

    #[test]
    fn mosaic_identity_checked() {
        // Each fragment is above 0.95, but B adds bases 5,000 to 6,000 at 0.9 when it comes
        // after A1 and A2. Counted first, B gives 10,546 matches over 11,000 bases.
        let a1 = fragment(0, &[(0, 5000)], &spread_mismatches(249));
        let a2 = fragment(1, &[(6000, 11_000)], &spread_mismatches(249));
        let b = fragment(2, &[(4000, 7900)], "1900=100X1900=");
        assert!(a1.alignment.score > b.alignment.score);
        let rule = DecisionRule::default();
        let decision = rule
            .decide(22_000, vec![a1.clone(), a2.clone(), b.clone()], |_| 0)
            .unwrap();
        assert!(decision.qualifies);
        assert_eq!(decision.mosaic.covered_bases, 11_000);
        assert_eq!(decision.mosaic.primary[0].fragment.contig_index, 2);
        assert!(decision.counts.reaches_identity(0.95));
        // At 0.959 the counts fall short; leaving out the pieces below the limit then loses
        // coverage, so the entry fails.
        let strict = DecisionRule {
            min_identity: 0.959,
            min_fragment_identity: 0.95,
            ..rule
        };
        let decision = strict.decide(22_000, vec![a1, a2, b], |_| 0).unwrap();
        assert!(!decision.qualifies);
        assert!(decision.counts.reaches_identity(0.959));
    }

    #[test]
    fn weak_addition_dropped() {
        // The 0.86 fragment adds 600 bases at 0.75 to an exact one, 0.906 in all; without it
        // the identity holds and 1,000 of 1,600 bases are still covered.
        let exact = fragment(0, &[(0, 1000)], "1000=");
        let weak = fragment(1, &[(500, 1600)], "500=150X450=");
        let rule = DecisionRule {
            min_query_coverage: "0.6".parse().unwrap(),
            min_identity: 0.95,
            min_fragment_identity: 0.8,
            ..DecisionRule::default()
        };
        let decision = rule
            .decide(1600, vec![exact.clone(), weak.clone()], |_| 0)
            .unwrap();
        assert!(decision.qualifies);
        assert_eq!(decision.mosaic.covered_bases, 1000);
        assert_eq!(decision.mosaic.primary.len(), 1);
        assert_eq!(matches_and_columns(&decision), (1000, 1000));
        let more = DecisionRule {
            min_query_coverage: "0.7".parse().unwrap(),
            ..rule
        };
        assert!(
            !more
                .decide(1600, vec![exact, weak], |_| 0)
                .unwrap()
                .qualifies
        );
    }

    #[test]
    fn many_segments_fast() {
        // 116,279 blocks of 40= 3X: each 40= counts alone, as 116,279 primary fragments.
        let blocks = 116_279;
        let long = fragment(0, &[(0, 43 * blocks)], &"40=3X".repeat(blocks as usize));
        let start = std::time::Instant::now();
        let decision = DecisionRule::default()
            .decide(43 * blocks, vec![long], |_| 0)
            .unwrap();
        assert!(decision.qualifies);
        assert_eq!(decision.mosaic.primary.len(), blocks as usize);
        assert_eq!(decision.mosaic.covered_bases, 40 * blocks);
        assert!(start.elapsed().as_secs() < 20, "{:?}", start.elapsed());
    }

    #[test]
    fn fragment_below_limit() {
        let weak = fragment(0, &[(0, 902)], &format!("{}3X", spread(17, 50)));
        let decision = DecisionRule {
            min_query_coverage: "0.5".parse().unwrap(),
            ..DecisionRule::default()
        }
        .decide(902, vec![weak], |_| 0)
        .unwrap();
        assert!(!decision.qualifies);
        assert_eq!(decision.dropped_fragments, 1);
        assert_eq!(decision.mosaic.covered_bases, 0);
    }

    #[test]
    fn identity_below_fails() {
        let low = fragment(1, &[(1000, 2000)], "700=300X");
        let below = fragment(0, &[(0, 1000)], "899=101X");
        let decision = rule("1")
            .decide(2000, vec![below, low.clone()], |_| 0)
            .unwrap();
        // 1,599 of 2,000 is below 0.8, and without the 0.7 fragment half the query is left.
        assert_eq!(matches_and_columns(&decision), (899, 1000));
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
        let decision = rule("0.5")
            .decide(10_000, vec![long.clone(), short.clone()], |_| 0)
            .unwrap();
        assert!(!decision.qualifies);
        // The long fragment is left out as it is below the limit; the short one alone remains.
        assert_eq!(matches_and_columns(&decision), (990, 1000));
        let decision = rule("0.1").decide(10_000, vec![long, short.clone()], |_| 0);
        assert!(decision.unwrap().qualifies);
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
        assert_eq!(matches_and_columns(&decision), (4500, 5000));
        assert_eq!(decision.mosaic.covered_bases, 5000);
        assert_eq!(decision.mosaic.primary.len(), 1);
        assert_eq!(decision.mosaic.alternatives.len(), 1);
        assert!(decision.qualifies);
        // A 65 percent copy is an alternative, so it does not lower the identity.
        let weak = fragment(1, &[(0, 5000)], "3250=1750X");
        let decision = rule("0.5")
            .decide(10_000, vec![first, weak], |_| 0)
            .unwrap();
        assert!(decision.qualifies);
        assert!((decision.identity() - 0.9).abs() < 1e-12);
    }

    #[test]
    fn alternative_left_out() {
        // A perfect fragment and a 370 column copy inside it at 81 percent.
        let perfect = fragment(0, &[(0, 1000)], "1000=");
        let weak = fragment(1, &[(200, 570)], "300=70X");
        let decision = DecisionRule {
            min_fragment_identity: 0.0,
            ..DecisionRule::default()
        }
        .decide(1000, vec![perfect, weak], |_| 0)
        .unwrap();
        assert!(decision.qualifies);
        assert_eq!(matches_and_columns(&decision), (1000, 1000));
        assert_eq!(decision.mosaic.alternatives.len(), 1);
    }

    #[test]
    fn conserved_copy_primary() {
        // The diverged copy is longer and scores higher, but the conserved one comes first.
        let conserved = fragment(0, &[(0, 1000)], "1000=");
        let diverged = fragment(1, &[(0, 1500)], "700=50X700=50X");
        assert!(diverged.alignment.score > conserved.alignment.score);
        let decision = DecisionRule {
            min_fragment_identity: 0.0,
            ..DecisionRule::default()
        }
        .decide(2000, vec![diverged, conserved], |_| 0)
        .unwrap();
        assert!(decision.qualifies);
        assert_eq!(decision.mosaic.primary[0].fragment.contig_index, 0);
        assert_eq!(decision.mosaic.primary[1].new_bases, 500);
        assert_eq!(decision.mosaic.covered_bases, 1500);
        // The diverged copy counts only over query bases 1000 to 1500.
        assert_eq!(matches_and_columns(&decision), (1450, 1500));
    }

    #[test]
    fn paralog_counted_little() {
        // A 76 percent paralog that adds 10 query bases counts only over those 10 bases.
        let source = fragment(0, &[(0, 1000)], "1000=");
        let paralog = fragment(1, &[(0, 1010)], "250X760=");
        let decision = DecisionRule {
            min_fragment_identity: 0.0,
            ..DecisionRule::default()
        }
        .decide(1010, vec![paralog, source], |_| 0)
        .unwrap();
        assert!(decision.qualifies);
        assert_eq!(decision.mosaic.primary.len(), 2);
        assert_eq!(matches_and_columns(&decision), (1010, 1010));
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
        assert!(!counts("39=", 39));
        assert!(counts("40=", 40));
        assert!(counts("70=30X", 100));
        assert!(counts("65=35X", 100));
        let strict = DecisionRule {
            min_fragment_identity: 0.7,
            ..rule
        };
        assert!(!strict.counts_fragment(&fragment(0, &[(0, 100)], "69=31X").alignment));
        // Four decimals would print 0.79996 as 0.8000, but the limit is compared exactly.
        let rounded = fragment(0, &[(0, 100_000)], "79996=20004X");
        assert_eq!(format!("{:.4}", rounded.alignment.identity()), "0.8000");
        let decision = rule.decide(100_000, vec![rounded], |_| 0).unwrap();
        assert!(!decision.qualifies);
        assert_eq!(decision.mosaic.covered_bases, 0);
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
        // Without it, the shared part is left out as its bases are below the identity limit.
        let decision = without.decide(3400, fragments(), lengths).unwrap();
        assert!(decision.qualifies);
        assert_eq!(matches_and_columns(&decision), (2970, 3000));
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
