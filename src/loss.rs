//! The proven bound on the chance that a query piece present in an entry shows no seed hit, and
//! the largest total chance of missing a qualifying entry that a search accepts.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

//! The bound holds for one continuous piece without indels, under independent substitutions at
//! rate `1 - identity`. Repeats, indels or clustered changes can miss more often.

use crate::seeds::{SeedError, SeedScheme};
use thiserror::Error;

/// Bad input to the loss model.
#[derive(Clone, Debug, Error, PartialEq)]
pub enum LossError {
    /// The seed scheme has a value outside its range.
    #[error(transparent)]
    Seeds(#[from] SeedError),
    /// The --max-miss value is not a chance above 0 and below 1.
    #[error("--max-miss must be a chance above 0 and below 1")]
    InvalidMaxMiss,
}

/// The proven bound on the chance that a piece of a given length and identity shows no seed hit.
#[derive(Clone, Copy, Debug)]
pub struct MissModel {
    k: usize,
    w: usize,
    /// Mismatches a query k-mer may have against a seed and still hit it (0, 1 or 2).
    mismatches: usize,
}

impl MissModel {
    /// The model for the k and w of a seed scheme, for seeds hit with at most `mismatches`
    /// mismatches (0, 1 or 2; more counts as 2).
    pub fn new(scheme: SeedScheme, mismatches: usize) -> Result<Self, LossError> {
        scheme.validate()?;
        Ok(Self {
            k: usize::from(scheme.k),
            w: usize::from(scheme.w),
            mismatches: mismatches.min(2),
        })
    }

    /// Mismatches a query k-mer may have against a seed and still hit it (0, 1 or 2).
    pub fn mismatches(&self) -> usize {
        self.mismatches
    }

    /// A proven upper bound on the chance that a piece present at `identity` shows no seed hit,
    /// under independent substitutions: each disjoint block of w + k - 1 bases holds a whole window,
    /// so a seed (if it passes the entropy limit), and these seeds are missed independently.
    pub fn miss_upper_bound(&self, stretch_bases: u64, identity: f64) -> f64 {
        if identity.is_nan() || identity <= 0.0 {
            return 1.0;
        }
        let (k, p) = (self.k as i32, identity.min(1.0));
        let mut missed = 1.0 - p.powi(k);
        if self.mismatches > 0 {
            missed -= f64::from(k) * (1.0 - p) * p.powi(k - 1);
        }
        if self.mismatches > 1 {
            missed -= f64::from(k * (k - 1) / 2) * (1.0 - p).powi(2) * p.powi(k - 2);
        }
        let blocks = stretch_bases / (self.w + self.k - 1) as u64;
        missed.clamp(0.0, 1.0).powf(blocks as f64)
    }

    /// The identity p* = 1 - (1 - p) R / (R - k + 1) for a copy of R bases that qualifies at p:
    /// its at most (1 - p) R random substitutions change each seed base with at most that chance.
    /// 0, so no bound, when R < k.
    pub fn qualifying_identity(&self, identity: f64, required: u64) -> f64 {
        let span = required.saturating_sub(self.k as u64 - 1);
        if span == 0 || identity.is_nan() {
            return 0.0;
        }
        let changed = (1.0 - identity) * required as f64 / span as f64;
        (1.0 - changed).clamp(0.0, 1.0)
    }

    /// A proven upper bound on the chance that a qualifying copy shows fewer than `hits` distinct
    /// hit starts, for copies of `required` to `longest` bases at `identity`. A copy of L bases
    /// holds n = floor(L / (w + k - 1)) disjoint blocks with one seed each, and its
    /// s = floor((1 - identity) L) substitutions fall on uniform positions. Fewer than `hits` hit
    /// starts need fewer than `hits` of these n seeds with at most `mismatches` substitutions. That
    /// chance is counted exactly up to `EXACT_BLOCKS` blocks, and bounded by a union over the sets
    /// of failing seeds above. The bound is the largest over L. A longer copy splits into pieces of
    /// `required` to 2 `required` - 1 bases, one of which qualifies on its own, so L stops there.
    pub fn hit_miss_bound(&self, identity: f64, required: u64, longest: u64, hits: u64) -> f64 {
        if hits == 0 {
            return 0.0;
        }
        if identity.is_nan() || identity <= 0.0 || required == 0 {
            return 1.0;
        }
        let span = (self.w + self.k - 1) as u64;
        let rate = 1.0 - identity.min(1.0);
        let last = longest.min(2 * required - 1).max(required);
        let factorials = ln_factorials(last.min((EXACT_BLOCKS + 1) * span));
        let (mut worst, mut previous) = (0.0f64, None);
        for length in required..=last {
            // The small term keeps a product such as 0.05 * 740 from rounding down below 37.
            let changes = ((rate * length as f64 + 1e-9).floor() as u64).min(length);
            let blocks = length / span;
            // At the same blocks and substitutions a longer copy spreads them wider: only the
            // shortest length of each pair can be the largest.
            if previous == Some((blocks, changes)) {
                continue;
            }
            previous = Some((blocks, changes));
            let union = self.union_hit_miss(length, changes, blocks, hits);
            if union <= worst {
                continue;
            }
            worst = match blocks <= EXACT_BLOCKS {
                true => self.exact_hit_miss(length, changes, blocks, hits, &factorials),
                false => union,
            }
            .max(worst);
        }
        worst.min(1.0)
    }

    /// The chance that one seed of a copy of `length` bases with `changes` uniform substitutions
    /// keeps at most `mismatches` of them.
    fn kept_chance(&self, length: u64, changes: u64) -> f64 {
        let (k, l, s) = (self.k as u64, length as f64, changes as f64);
        (0..=self.mismatches as u64)
            .filter(|&u| u <= changes && k - u <= length - changes)
            .map(|u| {
                let mut chance = choose(k, u);
                for i in 0..u {
                    chance *= (s - i as f64) / (l - i as f64);
                }
                for i in 0..k - u {
                    chance *= (l - s - i as f64) / (l - (u + i) as f64);
                }
                chance
            })
            .sum()
    }

    /// The union bound C(n, hits - 1) f^(n - hits + 1), f the chance that one seed fails: some
    /// n - hits + 1 seeds all fail, and seeds on disjoint bases fail at most independently.
    fn union_hit_miss(&self, length: u64, changes: u64, blocks: u64, hits: u64) -> f64 {
        if blocks < hits {
            return 1.0;
        }
        let failed = (1.0 - self.kept_chance(length, changes)).clamp(0.0, 1.0);
        (choose(blocks, hits - 1) * failed.powf((blocks - hits + 1) as f64)).min(1.0)
    }

    /// The exact chance that fewer than `hits` of `blocks` seeds keep at most `mismatches` of the
    /// `changes` substitutions placed uniformly on `length` bases: the placements are counted
    /// seed by seed, with the rest falling on the other bases.
    fn exact_hit_miss(
        &self,
        length: u64,
        changes: u64,
        blocks: u64,
        hits: u64,
        factorials: &[f64],
    ) -> f64 {
        let (k, s, hits) = (self.k, changes as usize, hits as usize);
        let ways: Vec<f64> = (0..=k as u64).map(|u| choose(k as u64, u)).collect();
        // counts[j][t]: placements of t substitutions on the seeds so far, j of which keep at
        // most `mismatches`, for j below `hits`.
        let mut counts = vec![vec![0.0f64; s + 1]; hits];
        counts[0][0] = 1.0;
        for _ in 0..blocks {
            let mut next = vec![vec![0.0f64; s + 1]; hits];
            for (j, row) in counts.iter().enumerate() {
                for (t, &count) in row.iter().enumerate().filter(|(_, c)| **c > 0.0) {
                    for (u, &way) in ways.iter().enumerate().take(s - t + 1) {
                        let kept = j + usize::from(u <= self.mismatches);
                        if kept < hits {
                            next[kept][t + u] += count * way;
                        }
                    }
                }
            }
            counts = next;
        }
        let rest = length - k as u64 * blocks;
        let ln_choose = |n: u64, r: u64| {
            factorials[n as usize] - factorials[r as usize] - factorials[(n - r) as usize]
        };
        let all = ln_choose(length, changes);
        let chance: f64 = (counts.iter())
            .flat_map(|row| row.iter().enumerate())
            .filter(|&(t, _)| changes - t as u64 <= rest)
            .map(|(t, &count)| count * (ln_choose(rest, changes - t as u64) - all).exp())
            .sum();
        chance.min(1.0)
    }

    /// The minimizer window of the model, in k-mers.
    pub fn w(&self) -> usize {
        self.w
    }

    /// The k-mer size of the model.
    pub fn k(&self) -> usize {
        self.k
    }
}

/// The blocks up to which `MissModel::hit_miss_bound` counts placements exactly.
const EXACT_BLOCKS: u64 = 16;

/// C(n, r) as a float.
fn choose(n: u64, r: u64) -> f64 {
    (0..r).fold(1.0, |ways, i| ways * (n - i) as f64 / (i + 1) as f64)
}

/// ln(n!) for n from 0 to `up_to`.
fn ln_factorials(up_to: u64) -> Vec<f64> {
    let mut sums = vec![0.0f64; up_to as usize + 1];
    for n in 1..sums.len() {
        sums[n] = sums[n - 1] + (n as f64).ln();
    }
    sums
}

/// The largest total chance of missing a qualifying entry (`--max-miss`): the seed lookups and the
/// coverage pruning together. It holds under the assumptions of `MissModel` only.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaxMiss {
    chance: f64,
}

impl Default for MaxMiss {
    /// 0.001, which is 0.1%.
    fn default() -> Self {
        Self { chance: 0.001 }
    }
}

impl MaxMiss {
    /// A limit of `chance`, which must be above 0 and below 1.
    pub fn new(chance: f64) -> Result<Self, LossError> {
        if chance.is_nan() || chance <= 0.0 || chance >= 1.0 {
            return Err(LossError::InvalidMaxMiss);
        }
        Ok(Self { chance })
    }

    /// The chance of the limit.
    pub fn chance(self) -> f64 {
        self.chance
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The proven bound: q^floor(L / 124), q the chance that one seed has more mismatches than
    /// the lookups allow.
    #[test]
    fn block_upper_bound() {
        let scheme = SeedScheme::default();
        let (exact, one) = (MissModel::new(scheme, 0), MissModel::new(scheme, 1));
        let (exact, one) = (exact.unwrap(), one.unwrap());
        let q0 = 1.0 - 0.95f64.powi(25);
        let q1 = q0 - 25.0 * 0.05 * 0.95f64.powi(24);
        assert_eq!(exact.miss_upper_bound(123, 0.95), 1.0);
        assert!((exact.miss_upper_bound(750, 0.95) - q0.powi(6)).abs() < 1e-12);
        assert!((one.miss_upper_bound(750, 0.95) - q1.powi(6)).abs() < 1e-12);
        // The values of the review at L = 750.
        assert!((exact.miss_upper_bound(750, 0.95) - 0.142_372_25).abs() < 1e-8);
        assert!((one.miss_upper_bound(750, 0.95) - 0.002_092_00).abs() < 1e-8);
        assert_eq!(exact.miss_upper_bound(124, 1.0), 0.0);
        assert_eq!(exact.miss_upper_bound(1000, f64::NAN), 1.0);
    }

    /// The bound for two mismatches: q_2 = 1 - sum over j <= 2 of C(25, j) (1 - p)^j p^(25 - j).
    #[test]
    fn two_mismatch_bound() {
        let two = MissModel::new(SeedScheme::default(), 2).unwrap();
        assert_eq!(two.mismatches(), 2);
        for p in [0.8f64, 0.85, 0.95] {
            let hit: f64 = [(1.0, 0), (25.0, 1), (300.0, 2)]
                .iter()
                .map(|&(ways, j)| ways * (1.0 - p).powi(j) * p.powi(25 - j))
                .sum();
            let q2 = 1.0 - hit;
            assert!(
                (two.miss_upper_bound(750, p) - q2.powi(6)).abs() < 1e-12,
                "{p}"
            );
        }
        assert_eq!(two.miss_upper_bound(123, 0.8), 1.0);
        assert_eq!(two.miss_upper_bound(124, 1.0), 0.0);
        let one = MissModel::new(SeedScheme::default(), 1).unwrap();
        assert!(two.miss_upper_bound(750, 0.8) < one.miss_upper_bound(750, 0.8));
    }

    /// At R = 1,116 and p = 0.95, p* = 1 - 0.05 * 1,116 / 1,092 = 0.948 90.
    #[test]
    fn qualifying_identity_lower() {
        let model = MissModel::new(SeedScheme::default(), 0).unwrap();
        let lower = model.qualifying_identity(0.95, 1_116);
        assert!((lower - 0.948_901).abs() < 1e-6, "{lower}");
        assert_eq!(model.qualifying_identity(1.0, 500), 1.0);
        // Up to k - 1 required bases there is no bound: the chance stays 1.
        assert_eq!(model.qualifying_identity(0.95, 24), 0.0);
        assert_eq!(
            model.miss_upper_bound(500, model.qualifying_identity(0.95, 24)),
            1.0
        );
        // With 25 bases, 1.25 changes may fall on the one seed span: no bound either.
        assert_eq!(model.qualifying_identity(0.95, 25), 0.0);
    }

    /// The values of the design study for m = 1,299, R = 650 and identity 0.95: the existing bound
    /// q_r(p*)^5, and the exact bound for fewer than h hit starts over copies of 650 to 1,299
    /// bases, largest at 740 bases with 37 substitutions.
    #[test]
    fn hit_miss_design_values() {
        let model = |r| MissModel::new(SeedScheme::default(), r).unwrap();
        let lower = model(0).qualifying_identity(0.95, 650);
        let close = |got: f64, want: f64| (got / want - 1.0).abs() < 1e-5;
        assert!(close(model(1).miss_upper_bound(650, lower), 0.007_437_94));
        assert!(close(
            model(2).miss_upper_bound(650, lower),
            0.000_050_001_8
        ));
        for (r, hits, want) in [
            (1, 1, 0.003_621_28),
            (1, 2, 0.047_690_7),
            (1, 3, 0.237_266),
            (2, 1, 0.000_005_966_23),
            (2, 2, 0.000_454_314),
            (2, 3, 0.010_998_3),
        ] {
            let got = model(r).hit_miss_bound(0.95, 650, 1_299, hits);
            assert!(close(got, want), "{r} {hits} {got}");
            let at_740 = model(r).hit_miss_bound(0.95, 740, 740, hits);
            assert!(close(at_740, want), "{r} {hits} {at_740}");
        }
        // At L = 650 alone, s = 32: the values of the first table column.
        let at_650 = model(2).hit_miss_bound(0.95, 650, 650, 2);
        assert!(close(at_650, 0.000_324_190), "{at_650}");
    }

    /// Skipping lengths with the same blocks and substitutions finds the same largest value as
    /// every length; the union bound is never below the exact count; one more hit start or a
    /// lower identity never lowers the bound.
    #[test]
    fn hit_miss_checks() {
        let two = MissModel::new(SeedScheme::default(), 2).unwrap();
        let factorials = ln_factorials(2_200);
        let mut every = 0.0f64;
        for length in 300..=599u64 {
            let changes = (0.1 * length as f64 + 1e-9).floor() as u64;
            let blocks = length / 124;
            let exact = two.exact_hit_miss(length, changes, blocks, 2, &factorials);
            assert!(exact <= two.union_hit_miss(length, changes, blocks, 2) + 1e-12);
            every = every.max(exact);
        }
        let skipped = two.hit_miss_bound(0.9, 300, 10_000, 2);
        assert!((skipped / every - 1.0).abs() < 1e-12, "{skipped} {every}");
        assert!(
            two.hit_miss_bound(0.9, 2_000, 4_000, 2) < two.hit_miss_bound(0.9, 1_000, 2_000, 2)
        );
        assert!(two.hit_miss_bound(0.95, 650, 1_299, 3) > two.hit_miss_bound(0.95, 650, 1_299, 2));
        assert!(two.hit_miss_bound(0.9, 650, 1_299, 2) > two.hit_miss_bound(0.95, 650, 1_299, 2));
        assert_eq!(two.hit_miss_bound(0.95, 650, 1_299, 0), 0.0);
        assert_eq!(two.hit_miss_bound(0.95, 200, 300, 2), 1.0);
        assert_eq!(two.hit_miss_bound(1.0, 248, 300, 2), 0.0);
        // A long query takes the union bound above 16 blocks, far below any target.
        let long = two.hit_miss_bound(0.95, 5_000, 10_000, 2);
        assert!(long > 0.0 && long < 1e-20, "{long}");
    }

    #[test]
    fn bad_seeding_rejected() {
        for (k, w) in [(0, 100), (33, 100), (25, 0)] {
            let scheme = SeedScheme {
                k,
                w,
                ..SeedScheme::default()
            };
            let error = MissModel::new(scheme, 0).unwrap_err();
            assert!(matches!(error, LossError::Seeds(_)), "{k} {w}");
        }
    }

    #[test]
    fn default_max_miss() {
        assert_eq!(MaxMiss::default().chance(), 0.001);
    }

    #[test]
    fn max_miss_checked() {
        for bad in [0.0, 1.0, -0.1, f64::NAN, 2.0] {
            assert_eq!(MaxMiss::new(bad).unwrap_err(), LossError::InvalidMaxMiss);
        }
    }
}
