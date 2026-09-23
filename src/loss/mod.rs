//! The chance that a query stretch present in an entry shows no seed hit, and the loss budget that
//! bounds every shortcut of a search.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

//! The chance is exact for one seed layout under independent substitutions at rate `1 - identity`;
//! the model averages it over the layouts of random sequences in a table over length and identity.

use crate::seeds::{SeedError, SeedScheme};
use rand_core::{Rng, SeedableRng};
use rand_pcg::Pcg64Mcg;
use rayon::prelude::*;
use thiserror::Error;

/// Bad input to the loss model.
#[derive(Clone, Debug, Error, PartialEq)]
pub enum LossError {
    /// The seed scheme has a value outside its range.
    #[error(transparent)]
    Seeds(#[from] SeedError),
    /// The model was asked to average over no layouts.
    #[error("the miss model needs at least one layout per stretch length")]
    NoLayouts,
    /// The loss budget is not a chance above 0 and below 1.
    #[error("the loss budget must be a chance above 0 and below 1")]
    InvalidBudget,
}

/// The chance that a stretch of a given length and identity shows no seed hit.
#[derive(Clone, Debug)]
pub struct MissModel {
    k: usize,
    w: usize,
    lengths: Vec<f64>,
    identities: Vec<f64>,
    /// Mismatches a query k-mer may have against a seed and still hit it (0 or 1).
    mismatches: usize,
    /// Natural logarithm of the chance, row by length, column by identity.
    ln_miss: Vec<f64>,
}

impl MissModel {
    /// The model for the k and w of a seed scheme, with the table averaged over 1,000 layouts
    /// per length. The table of the default k and w is stored, so it is read, not averaged.
    pub fn new(scheme: SeedScheme) -> Result<Self, LossError> {
        Self::stored_or_averaged(scheme, 0)
    }

    /// As `new`, for queries that also look up every k-mer with one substitution: a seed is hit
    /// when the query holds it with at most one mismatch.
    pub fn one_mismatch(scheme: SeedScheme) -> Result<Self, LossError> {
        Self::stored_or_averaged(scheme, 1)
    }

    /// As `new`, averaging over `layouts` random layouts per stretch length, for seeds hit with
    /// at most `mismatches` mismatches (0 or 1).
    pub fn with_layouts(
        scheme: SeedScheme,
        layouts: usize,
        mismatches: usize,
    ) -> Result<Self, LossError> {
        scheme.validate()?;
        if layouts == 0 {
            return Err(LossError::NoLayouts);
        }
        let (k, w) = (usize::from(scheme.k), usize::from(scheme.w));
        let lengths = table_lengths(k);
        let identities = identity_grid();
        let rows: Vec<Vec<f64>> = lengths
            .par_iter()
            .map(|&length| average_row(length as usize, k, w, layouts, mismatches, &identities))
            .collect();
        let ln_miss = rows
            .into_iter()
            .flatten()
            .map(|chance| chance.max(1e-300).ln())
            .collect();
        Ok(Self {
            k,
            w,
            lengths,
            identities,
            mismatches: mismatches.min(1),
            ln_miss,
        })
    }

    /// The stored table for the default k and w, else a new average over `LAYOUTS` layouts.
    fn stored_or_averaged(scheme: SeedScheme, mismatches: usize) -> Result<Self, LossError> {
        let default = SeedScheme::default();
        if (scheme.k, scheme.w) != (default.k, default.w) {
            return Self::with_layouts(scheme, LAYOUTS, mismatches);
        }
        scheme.validate()?;
        let stored = match mismatches {
            0 => DEFAULT_LN_MISS,
            _ => ONE_MISMATCH_TABLE,
        };
        let (values, _) = stored.as_chunks::<{ size_of::<f64>() }>();
        let ln_miss = values
            .iter()
            .map(|&bytes| f64::from_le_bytes(bytes))
            .collect();
        Ok(Self {
            k: usize::from(scheme.k),
            w: usize::from(scheme.w),
            lengths: table_lengths(usize::from(scheme.k)),
            identities: identity_grid(),
            mismatches: mismatches.min(1),
            ln_miss,
        })
    }

    /// Mismatches a query k-mer may have against a seed and still hit it.
    pub fn mismatches(&self) -> usize {
        self.mismatches
    }

    /// The chance that a stretch of `stretch_bases` bases present at `identity` shows no seed hit,
    /// interpolated in the logarithm between table points. Beyond the table the answer errs on the
    /// safe side: the longest length and highest identity hold, and a lower identity gives 1.
    pub fn miss_chance(&self, stretch_bases: u64, identity: f64) -> f64 {
        if stretch_bases < self.k as u64 || identity.is_nan() || identity < self.identities[0] {
            return 1.0;
        }
        if identity >= 1.0 {
            return 0.0;
        }
        let (row, row_weight) = locate(&self.lengths, stretch_bases as f64);
        let (column, column_weight) = locate(&self.identities, identity);
        let width = self.identities.len();
        let at = |r: usize, c: usize| self.ln_miss[r * width + c];
        let along = |r: usize| {
            at(r, column) * (1.0 - column_weight)
                + at(r, (column + 1).min(width - 1)) * column_weight
        };
        let next_row = (row + 1).min(self.lengths.len() - 1);
        let ln = along(row) * (1.0 - row_weight) + along(next_row) * row_weight;
        ln.exp().min(1.0)
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

/// The largest chance of missing a qualifying entry that all shortcuts of one search may add up to.
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

/// The stretch lengths of the table, in bases. Lengths below the k-mer size are left out.
const LENGTHS: [u64; 21] = [
    25, 35, 50, 70, 100, 140, 200, 280, 400, 560, 800, 1000, 1250, 1600, 2000, 2500, 3200, 4000,
    5000, 6500, 8000,
];

/// The lengths of the table for k-mers of `k` bases.
fn table_lengths(k: usize) -> Vec<f64> {
    LENGTHS
        .into_iter()
        .filter(|&length| length >= k as u64)
        .map(|length| length as f64)
        .collect()
}

/// The identities of the table: 0.700 to 0.990 in steps of 0.005, then 0.995 and 0.999.
fn identity_grid() -> Vec<f64> {
    let mut grid: Vec<f64> = (0..=58)
        .map(|step| 0.70 + 0.005 * f64::from(step))
        .collect();
    grid.extend([0.995, 0.999]);
    grid
}

/// The layouts a table averages over per stretch length.
const LAYOUTS: usize = 1000;

/// The table of the default k and w, averaged over `LAYOUTS` layouts per length, as
/// little-endian f64 values row by row; a test checks it against a new average.
const DEFAULT_LN_MISS: &[u8] = include_bytes!("default_table.bin");

/// As `DEFAULT_LN_MISS`, for seeds hit with at most one mismatch.
const ONE_MISMATCH_TABLE: &[u8] = include_bytes!("one_mismatch_table.bin");

/// The largest k the model supports, as the seed scheme does.
const MAX_K: usize = 32;

/// The average miss chance over random layouts of one stretch length, for every identity.
fn average_row(
    length: usize,
    k: usize,
    w: usize,
    layouts: usize,
    mismatches: usize,
    identities: &[f64],
) -> Vec<f64> {
    // PCG gives the same numbers in every version, so the table never changes between builds.
    let mut random = Pcg64Mcg::seed_from_u64(0x6A09_E667_F3BC_C908 ^ length as u64);
    let mut sums = vec![0.0f64; identities.len()];
    for _ in 0..layouts {
        let seeds = random_layout(length, k, w, &mut random);
        for (sum, &identity) in sums.iter_mut().zip(identities) {
            *sum += if mismatches == 0 {
                layout_miss_chance(&seeds, k, identity)
            } else {
                one_mismatch_chance(&seeds, k, identity)
            };
        }
    }
    sums.iter().map(|sum| sum / layouts as f64).collect()
}

/// The base positions where the seeds of a random stretch end (the last base of the seed k-mer),
/// in increasing order. The hashes of the k-mers of a random sequence are independent and uniform.
fn random_layout(length: usize, k: usize, w: usize, random: &mut Pcg64Mcg) -> Vec<usize> {
    let kmers = length + 1 - k;
    let hashes: Vec<u64> = (0..kmers).map(|_| random.next_u64()).collect();
    let mut seeds = Vec::new();
    if kmers < w {
        let smallest = (0..kmers).min_by_key(|&i| hashes[i]).unwrap_or(0);
        seeds.push(smallest + k - 1);
        return seeds;
    }
    // A monotone queue of k-mer positions with rising hashes gives the minimum of each minimizer window.
    let mut queue = std::collections::VecDeque::new();
    for position in 0..kmers {
        while queue
            .back()
            .is_some_and(|&back| hashes[back] > hashes[position])
        {
            queue.pop_back();
        }
        queue.push_back(position);
        if queue[0] + w <= position {
            queue.pop_front();
        }
        if position + 1 >= w && seeds.last() != Some(&(queue[0] + k - 1)) {
            seeds.push(queue[0] + k - 1);
        }
    }
    seeds
}

/// The exact chance that every seed has a mutated base, for seeds ending at `seed_ends`. The state
/// is the number of unmutated bases before the current seed end; a seed is hit when it reaches
/// `k`, so only states below `k` are kept.
fn layout_miss_chance(seed_ends: &[usize], k: usize, identity: f64) -> f64 {
    let mut powers = [1.0f64; MAX_K + 1];
    for t in 1..=k {
        powers[t] = powers[t - 1] * identity;
    }
    let mut state = [0.0f64; MAX_K];
    state[0] = 1.0;
    // One position before the first base.
    let mut previous_end = -1i64;
    for &end in seed_ends {
        let gap = (end as i64 - previous_end) as usize;
        let total: f64 = state[..k].iter().sum();
        let mut next = [0.0f64; MAX_K];
        // A mutation inside the gap leaves `run` unmutated bases at the seed end.
        for (run, slot) in next.iter_mut().enumerate().take(gap.min(k)) {
            *slot = total * powers[run] * (1.0 - identity);
        }
        // No mutation in the gap extends the earlier run.
        for run in 0..k.saturating_sub(gap) {
            next[run + gap] += state[run] * powers[gap.min(k)];
        }
        state = next;
        previous_end = end as i64;
    }
    state[..k].iter().sum()
}

/// The exact chance that every seed, ending at `seed_ends`, has two or more mutated bases. State
/// (a, b): unmutated bases after the last mutation and bases after the second-last one; only
/// b < k misses a seed, so only those are kept, and `free` holds the mass with no mutation yet.
fn one_mismatch_chance(seed_ends: &[usize], k: usize, identity: f64) -> f64 {
    let (p, q) = (identity, 1.0 - identity);
    let mut powers = [1.0f64; MAX_K + 1];
    for t in 1..=k {
        powers[t] = powers[t - 1] * p;
    }
    let at = |a: usize, b: usize| a * MAX_K + b;
    let mut state = vec![0.0f64; MAX_K * MAX_K];
    let mut free = 1.0f64;
    let mut previous_end = -1i64;
    for &end in seed_ends {
        let gap = (end as i64 - previous_end) as usize;
        let total = free + state.iter().sum::<f64>();
        let mut next = vec![0.0f64; MAX_K * MAX_K];
        // Two or more mutations in the gap: the last two alone set the state.
        for t in 0..k {
            for b in t + 1..k {
                let u = b - t - 1;
                if t + u + 2 <= gap {
                    next[at(t, b)] += total * q * q * powers[t + u];
                }
            }
        }
        if gap < k {
            for a in 0..k {
                for b in a + 1..k {
                    let mass = state[at(a, b)];
                    if mass == 0.0 {
                        continue;
                    }
                    // Exactly one mutation: the old last mutation becomes the second-last.
                    if a + gap < k {
                        for t in 0..gap {
                            next[at(t, a + gap)] += mass * q * powers[gap - 1];
                        }
                    }
                    // No mutation: both distances grow.
                    if b + gap < k {
                        next[at(a + gap, b + gap)] += mass * powers[gap];
                    }
                }
            }
        }
        free = 0.0;
        state = next;
        previous_end = end as i64;
    }
    free + state.iter().sum::<f64>()
}

/// The index of the grid point at or below `value` and the share of the way to the next point.
/// Values beyond the grid end stay on its last point.
fn locate(grid: &[f64], value: f64) -> (usize, f64) {
    let last = grid.len() - 1;
    if value >= grid[last] {
        return (last, 0.0);
    }
    let below = grid
        .partition_point(|&point| point <= value)
        .saturating_sub(1);
    let weight = (value - grid[below]) / (grid[below + 1] - grid[below]);
    (below, weight)
}

#[cfg(test)]
mod tests;
