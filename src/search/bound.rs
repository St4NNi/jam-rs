//! The coverage bound of an entry: the most query bases it can still cover, given where its seed
//! hits lie on the query.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::loss::{LossBudget, MissModel};

/// Computes coverage bounds for one query against one shard.
#[derive(Clone, Copy, Debug)]
pub struct CoverageBound<'a> {
    model: &'a MissModel,
    budget: LossBudget,
    /// The lowest identity at which an entry still qualifies.
    identity: f64,
    query_length: u64,
    circular: bool,
    /// The query coverage, in bases, that an entry must reach.
    required_bases: u64,
    /// Query starts that count as possibly hit for every entry, ascending.
    masked_starts: &'a [u32],
}

impl<'a> CoverageBound<'a> {
    /// A bound for a query of `query_length` bases, spending `budget` per entry. The ascending
    /// `masked_starts` count as possibly hit in every entry, as low-complexity stretches have no
    /// query k-mers.
    pub fn new(
        model: &'a MissModel,
        budget: LossBudget,
        identity: f64,
        query_length: u64,
        circular: bool,
        required_bases: u64,
        masked_starts: &'a [u32],
    ) -> Self {
        Self {
            model,
            budget,
            identity,
            query_length,
            circular,
            required_bases,
            masked_starts,
        }
    }

    /// The highest query coverage, in bases, of an entry with seed hits (or possible seed hits)
    /// at the ascending query starts `positions`: the bases their k-mers cover, plus at most the
    /// allowance `a` of each stretch without a seed hit.
    pub fn bound(&self, positions: &[u32]) -> u64 {
        // Assumes a stretch holds at most one present piece: the decision counts only long
        // fragments, and the part of an entry that matches a stretch aligns as one piece.
        let stretches = if self.masked_starts.is_empty() {
            self.stretches(positions)
        } else {
            let mut all = [positions, self.masked_starts].concat();
            all.sort_unstable();
            self.stretches(&all)
        };
        let allowance = self.allowance(stretches.len() as u64);
        let left_out: u64 = stretches
            .iter()
            .map(|&bases| bases.saturating_sub(allowance))
            .sum();
        self.query_length - left_out
    }

    /// The shortest length `a` with `stretches * P_miss(a) <= B` (a union bound over the
    /// stretches, B the stage's share of the loss budget), or the query length if none fits.
    fn allowance(&self, stretches: u64) -> u64 {
        let fits = |length: u64| {
            stretches as f64 * self.model.miss_chance(length, self.identity) <= self.budget.chance()
        };
        if stretches == 0 || !fits(self.query_length) {
            return self.query_length;
        }
        // The miss chance never rises with the length, so the lengths that fit form a tail.
        let (mut low, mut high) = (0, self.query_length);
        while high - low > 1 {
            let middle = low + (high - low) / 2;
            if fits(middle) {
                high = middle;
            } else {
                low = middle;
            }
        }
        high
    }

    /// Whether an entry with this bound can still reach the required query coverage.
    pub fn keeps(&self, bound: u64) -> bool {
        bound >= self.required_bases
    }

    /// The lengths of the stretches of query bases that no k-mer with a seed hit covers. On a
    /// circular query the stretch over the end joins both ends.
    fn stretches(&self, positions: &[u32]) -> Vec<u64> {
        let k = self.model.k() as u64;
        let length = self.query_length;
        let (Some(&first), Some(&last)) = (positions.first(), positions.last()) else {
            return vec![length];
        };
        let mut stretches: Vec<u64> = positions
            .windows(2)
            .map(|pair| u64::from(pair[1]).saturating_sub(u64::from(pair[0]) + k))
            .filter(|&bases| bases > 0)
            .collect();
        // A k-mer of a circular query may run past the end and cover the first bases again.
        let (first, last_end) = (u64::from(first), u64::from(last) + k);
        let ends = if self.circular {
            vec![(first + length).saturating_sub(last_end)]
        } else {
            vec![first, length.saturating_sub(last_end)]
        };
        stretches.extend(ends.into_iter().filter(|&bases| bases > 0));
        stretches
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::miss_model as model;

    fn bound(length: u64, circular: bool) -> CoverageBound<'static> {
        CoverageBound::new(
            model(),
            LossBudget::default(),
            0.95,
            length,
            circular,
            0,
            &[],
        )
    }

    #[test]
    fn stretch_lengths() {
        let linear = bound(1_000, false);
        assert_eq!(linear.stretches(&[]), vec![1_000]);
        assert_eq!(
            linear.stretches(&[100, 110, 500, 900]),
            vec![369, 379, 100, 79]
        );
        let circular = bound(1_000, true);
        assert_eq!(
            circular.stretches(&[100, 110, 500, 900]),
            vec![369, 379, 179]
        );
        // The k-mer at 990 runs over the end and covers bases 0 to 10.
        assert_eq!(circular.stretches(&[50, 990]), vec![919, 39]);
    }

    #[test]
    fn allowance_fits_budget() {
        let linear = bound(10_000, false);
        let one = linear.allowance(1);
        assert!(model().miss_chance(one, 0.95) <= 1e-4);
        assert!(model().miss_chance(one - 1, 0.95) > 1e-4);
        assert!(linear.allowance(40) > one);
        assert_eq!(linear.allowance(0), 10_000);
        assert_eq!(bound(10_000, false).bound(&[]), one);
    }

    #[test]
    fn stretches_keep_allowance() {
        let linear = bound(10_000, false);
        // 40 k-mers with seed hits, 40 bases apart, then 7,600 bases without one.
        let positions: Vec<u32> = (0..40).map(|i| 61 * i).collect();
        let allowance = linear.allowance(40);
        let covered = 40 * 21;
        assert_eq!(
            linear.bound(&positions),
            covered + 39 * 40 + allowance.min(7_600)
        );
        assert!(allowance < 7_600);
    }

    #[test]
    fn circular_ends_join() {
        // Stretches of 250 and 300 bases at the ends: each is shorter than the allowance, but
        // together they are not.
        let positions: Vec<u32> = (250..=5_679).collect();
        let linear = bound(6_000, false);
        assert!(linear.allowance(2) > 300);
        assert_eq!(linear.bound(&positions), 6_000);
        let circular = bound(6_000, true);
        assert_eq!(circular.stretches(&positions), vec![550]);
        let allowance = circular.allowance(1);
        assert!(allowance < 550);
        assert_eq!(circular.bound(&positions), 6_000 - 550 + allowance);
    }

    #[test]
    fn masked_starts_count() {
        let positions: Vec<u32> = (0..40).map(|i| 61 * i).collect();
        let plain = bound(10_000, false);
        let masked_starts: Vec<u32> = (2_400..10_000).step_by(21).chain([9_979]).collect();
        let masked = CoverageBound::new(
            model(),
            LossBudget::default(),
            0.95,
            10_000,
            false,
            0,
            &masked_starts,
        );
        assert!(plain.bound(&positions) < 10_000);
        assert_eq!(masked.bound(&positions), 10_000);
        assert_eq!(
            masked.bound(&[]),
            10_000 - 2_400u64.saturating_sub(masked.allowance(1))
        );
    }
}
