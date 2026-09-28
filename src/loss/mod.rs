//! The modeled chance that a query stretch present in an entry shows no seed hit, and the loss
//! budget that the shortcuts of a search aim to stay within.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

//! The chance is exact for one seed layout under independent substitutions at rate `1 - identity`
//! and no indels. The table holds an upper value of its mean over the layouts of random sequences,
//! so it is a model estimate: repeats, indels or clustered changes can miss more often.

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

/// The modeled chance that a stretch of a given length and identity shows no seed hit.
#[derive(Clone, Debug)]
pub struct MissModel {
    k: usize,
    w: usize,
    lengths: Vec<f64>,
    identities: Vec<f64>,
    /// Mismatches a query k-mer may have against a seed and still hit it (0 or 1).
    mismatches: usize,
    /// Natural logarithm of the chance, row by length, column by identity, never rising along
    /// either.
    ln_miss: Vec<f64>,
}

impl MissModel {
    /// The model for the k and w of a seed scheme, with the table taken over 1,000 layouts per
    /// length. The table of the default k and w is stored, so it is read, not computed.
    pub fn new(scheme: SeedScheme) -> Result<Self, LossError> {
        Self::stored_or_averaged(scheme, 0)
    }

    /// As `new`, for queries that also look up every k-mer with one substitution: a seed is hit
    /// when the query holds it with at most one mismatch.
    pub fn one_mismatch(scheme: SeedScheme) -> Result<Self, LossError> {
        Self::stored_or_averaged(scheme, 1)
    }

    /// As `new`, taking `layouts` random layouts per stretch length, for seeds hit with at most
    /// `mismatches` mismatches (0 or 1).
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
            .map(|&length| {
                let row = average_row(length as usize, k, w, layouts, mismatches, &identities);
                row.iter()
                    .map(|&(mean, error)| upper_value(mean, error))
                    .collect()
            })
            .collect();
        let ln_miss = never_rising(rows, identities.len());
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

    /// The modeled chance that a stretch of `stretch_bases` bases present at `identity` shows no
    /// seed hit: the table point at or below both values, as the chance falls with each of them.
    /// At identity 1 a stretch holding a whole minimizer window always has a seed.
    pub fn miss_chance(&self, stretch_bases: u64, identity: f64) -> f64 {
        let length = stretch_bases as f64;
        if identity.is_nan() || identity < self.identities[0] || length < self.lengths[0] {
            return 1.0;
        }
        if identity >= 1.0 {
            let window_bases = (self.w + self.k - 1) as u64;
            return if stretch_bases >= window_bases {
                0.0
            } else {
                1.0
            };
        }
        let at_or_below = |grid: &[f64], value: f64| grid.partition_point(|&p| p <= value) - 1;
        let row = at_or_below(&self.lengths, length);
        let column = at_or_below(&self.identities, identity);
        self.ln_miss[row * self.identities.len() + column]
            .exp()
            .min(1.0)
    }

    /// A proven upper bound on the chance that a stretch present at `identity` shows no seed hit,
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

/// The modeled chance of losing a qualifying entry that the coverage bounds of one search may
/// spend. It holds under the assumptions of `MissModel` only, so it is an estimate, not a limit.
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

/// The lengths of the table for k-mers of `k` bases: 25 bases, then 10% longer each step up to
/// 8,000 (61 lengths), rounded. Lengths below the k-mer size are left out.
fn table_lengths(k: usize) -> Vec<f64> {
    std::iter::successors(Some(25.0f64), |length| Some(length * 1.1))
        .map(f64::round)
        .take_while(|&length| length <= 8000.0)
        .filter(|&length| length >= k as f64)
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

/// The layouts a table takes per stretch length.
const LAYOUTS: usize = 1000;

/// The table of the default k and w, taken over `LAYOUTS` layouts per length, as little-endian
/// f64 values row by row; a test checks it against a new computation.
const DEFAULT_LN_MISS: &[u8] = include_bytes!("default_table.bin");

/// As `DEFAULT_LN_MISS`, for seeds hit with at most one mismatch.
const ONE_MISMATCH_TABLE: &[u8] = include_bytes!("one_mismatch_table.bin");

/// The largest k the model supports, as the seed scheme does.
const MAX_K: usize = 32;

/// The upper value of a sample mean: three standard errors above it, at most 1.
fn upper_value(mean: f64, error: f64) -> f64 {
    (mean + 3.0 * error).min(1.0)
}

/// The logarithms of the rows, each value lowered to the least value at a shorter length or a
/// lower identity. The true chance never rises along either, so these stay upper values.
fn never_rising(rows: Vec<Vec<f64>>, width: usize) -> Vec<f64> {
    let mut ln_miss: Vec<f64> = (rows.into_iter().flatten())
        .map(|chance| chance.max(1e-300).ln())
        .collect();
    for at in width..ln_miss.len() {
        ln_miss[at] = ln_miss[at].min(ln_miss[at - width]);
    }
    for at in 1..ln_miss.len() {
        if at % width != 0 {
            ln_miss[at] = ln_miss[at].min(ln_miss[at - 1]);
        }
    }
    ln_miss
}

/// The mean miss chance over random layouts of one stretch length and its standard error, for
/// every identity. A single layout gives no error.
fn average_row(
    length: usize,
    k: usize,
    w: usize,
    layouts: usize,
    mismatches: usize,
    identities: &[f64],
) -> Vec<(f64, f64)> {
    // PCG gives the same numbers in every version, so the table never changes between builds.
    let mut random = Pcg64Mcg::seed_from_u64(0x6A09_E667_F3BC_C908 ^ length as u64);
    let mut sums = vec![(0.0f64, 0.0f64); identities.len()];
    for _ in 0..layouts {
        let seeds = random_layout(length, k, w, &mut random);
        for ((sum, squares), &identity) in sums.iter_mut().zip(identities) {
            let chance = if mismatches == 0 {
                layout_miss_chance(&seeds, k, identity)
            } else {
                one_mismatch_chance(&seeds, k, identity)
            };
            *sum += chance;
            *squares += chance * chance;
        }
    }
    let n = layouts as f64;
    (sums.iter())
        .map(|&(sum, squares)| {
            let mean = sum / n;
            let variance = (squares - n * mean * mean).max(0.0) / (n - 1.0).max(1.0);
            (mean, (variance / n).sqrt())
        })
        .collect()
}

/// The ascending last bases of the seeds of a random stretch inside a longer random contig: its
/// k-mers and `w - 1` on each side get independent uniform hashes, and seeds count only if their
/// k-mer lies inside. So a stretch shorter than one window may have no seed.
fn random_layout(length: usize, k: usize, w: usize, random: &mut Pcg64Mcg) -> Vec<usize> {
    let inside = length + 1 - k;
    let kmers = inside + 2 * (w - 1);
    let hashes: Vec<u64> = (0..kmers).map(|_| random.next_u64()).collect();
    let mut seeds: Vec<usize> = Vec::new();
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
        // The seed's k-mer start inside the stretch, if it lies there.
        let start = queue[0].checked_sub(w - 1).filter(|&start| start < inside);
        if let Some(end) = start.map(|start| start + k - 1)
            && position + 1 >= w
            && seeds.last() != Some(&end)
        {
            seeds.push(end);
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

#[cfg(test)]
mod tests;
