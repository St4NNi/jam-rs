//! The proven bound on the chance that a query piece present in an entry shows no seed hit, and
//! the loss budget that the shortcuts of a search stay within.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

//! The bound holds for one continuous piece without indels, under independent substitutions at
//! rate `1 - identity`. Repeats, indels or clustered changes can miss more often.

use crate::seeds::{SeedError, SeedScheme};
use thiserror::Error;

/// The table estimate of the miss chance, no longer used; it goes next.
#[allow(dead_code)]
mod estimate;

/// Bad input to the loss model.
#[derive(Clone, Debug, Error, PartialEq)]
pub enum LossError {
    /// The seed scheme has a value outside its range.
    #[error(transparent)]
    Seeds(#[from] SeedError),
    /// The loss budget is not a chance above 0 and below 1.
    #[error("the loss budget must be a chance above 0 and below 1")]
    InvalidBudget,
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

    /// The k-mer size of the model.
    pub fn k(&self) -> usize {
        self.k
    }

    /// The minimizer window of the model, in k-mers.
    pub fn w(&self) -> usize {
        self.w
    }
}

/// The chance of losing a qualifying entry that the shortcuts of one search may spend. It holds
/// under the assumptions of `MissModel` only.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LossBudget {
    chance: f64,
}

impl Default for LossBudget {
    /// 0.0001, which is 0.01%.
    fn default() -> Self {
        Self { chance: 0.0001 }
    }
}

impl LossBudget {
    /// A budget of `chance`, which must be above 0 and below 1.
    pub fn new(chance: f64) -> Result<Self, LossError> {
        if chance.is_nan() || chance <= 0.0 || chance >= 1.0 {
            return Err(LossError::InvalidBudget);
        }
        Ok(Self { chance })
    }

    /// The chance of the budget.
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
    fn default_budget() {
        assert_eq!(LossBudget::default().chance(), 0.0001);
    }

    #[test]
    fn budget_is_chance() {
        for bad in [0.0, 1.0, -0.1, f64::NAN, 2.0] {
            assert_eq!(LossBudget::new(bad).unwrap_err(), LossError::InvalidBudget);
        }
    }
}
