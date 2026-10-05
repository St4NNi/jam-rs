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
    /// The present bases that a stretch without a seed hit may hold, for this query and these
    /// settings; it never depends on the hits, so a further hit can only raise a bound.
    allowance: u64,
}

impl<'a> CoverageBound<'a> {
    /// A bound for a query of `query_length` bases, spending `budget` per entry.
    pub fn new(
        model: &'a MissModel,
        budget: LossBudget,
        identity: f64,
        query_length: u64,
        circular: bool,
        required_bases: u64,
    ) -> Self {
        let mut bound = Self {
            model,
            budget,
            identity,
            query_length,
            circular,
            required_bases,
            allowance: query_length,
        };
        bound.allowance = bound.proven_allowance();
        bound
    }

    /// The highest query coverage, in bases, of an entry with seed hits (or possible seed hits)
    /// at the ascending query starts `positions`: the bases their k-mers cover, plus at most the
    /// allowance `a` of each stretch without a seed hit.
    pub fn bound(&self, positions: &[u32]) -> u64 {
        // Assumes a stretch holds at most one present piece, as the model of the miss chance
        // does; several pieces separated by unrelated bases are not covered.
        let left_out: u64 = (self.stretches(positions).iter())
            .map(|&bases| bases.saturating_sub(self.allowance))
            .sum();
        self.query_length - left_out
    }

    /// The proven bound on the chance that a present piece of the required length shows no seed
    /// hit: the chance to lose a qualifying entry before any bound sees it.
    pub fn seed_miss_bound(&self) -> f64 {
        let required = self.required_bases.min(self.query_length);
        self.model.miss_upper_bound(required, self.identity)
    }

    /// The shortest length `a` with `n * U(a + 1) <= B`, U the proven miss bound and n the fixed
    /// starts of a seedless stretch of a + 1 bases (L on a circular query, else L - a), or L.
    /// It holds when each stretch meets at most one coherent present piece.
    fn proven_allowance(&self) -> u64 {
        let length = self.query_length;
        let fits = |allowance: u64| {
            let starts = if self.circular {
                length
            } else {
                length.saturating_sub(allowance).max(1)
            };
            let miss = self.model.miss_upper_bound(allowance + 1, self.identity);
            starts as f64 * miss <= self.budget.chance()
        };
        smallest_fit(length, fits)
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

/// The least length up to `length` that `fits`, or `length`. No term of a fit rises with the
/// length, so the lengths that fit form a tail.
fn smallest_fit(length: u64, fits: impl Fn(u64) -> bool) -> u64 {
    if !fits(length) {
        return length;
    }
    let (mut low, mut high) = (0, length);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::search::QueryKmers;
    use crate::seeds::{SeedScheme, Seeder};
    use crate::test_support::miss_model as model;
    use rand::RngExt;

    fn bound(length: u64, circular: bool) -> CoverageBound<'static> {
        CoverageBound::new(model(), LossBudget::default(), 0.95, length, circular, 0)
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
        // The test model has k 21 and w 20.
        let fits = |allowance: u64| {
            let miss = model().miss_upper_bound(allowance + 1, 0.95);
            (10_000 - allowance) as f64 * miss <= 1e-4
        };
        assert!(fits(linear.allowance) && !fits(linear.allowance - 1));
        assert_eq!(linear.bound(&[]), linear.allowance);
        // A longer query can hold more long stretches, so it needs a longer allowance.
        assert!(bound(1_000_000, false).allowance > linear.allowance);
    }

    #[test]
    fn stretches_keep_allowance() {
        let linear = bound(10_000, false);
        // 40 k-mers with seed hits, 40 bases apart, then 7,600 bases without one.
        let positions: Vec<u32> = (0..40).map(|i| 61 * i).collect();
        let covered = 40 * 21;
        assert_eq!(
            linear.bound(&positions),
            covered + 39 * 40 + linear.allowance.min(7_600)
        );
        assert!(linear.allowance < 7_600);
    }

    #[test]
    fn circular_ends_join() {
        // Stretches of three quarters of the allowance at both ends: each is shorter than the
        // allowance, but together they are not.
        let length = 20_000;
        let circular = bound(length, true);
        let end = circular.allowance * 3 / 4;
        let positions: Vec<u32> = (end as u32..=(length - end - 21) as u32).collect();
        let linear = bound(length, false);
        assert!(linear.allowance > end);
        assert_eq!(linear.bound(&positions), length);
        assert_eq!(circular.stretches(&positions), vec![2 * end]);
        assert_eq!(
            circular.bound(&positions),
            length - 2 * end + circular.allowance
        );
    }

    /// The model of the default k 25 and w 100.
    fn default_model() -> &'static MissModel {
        static MODEL: std::sync::OnceLock<MissModel> = std::sync::OnceLock::new();
        MODEL.get_or_init(|| MissModel::new(SeedScheme::default()).unwrap())
    }

    /// Hit starts 1 to 3,075 leave stretches of 1 and 6,900 bases; a hit at 0 removes the first.
    /// With an allowance from the number of stretches, the extra hit lowered the bound from
    /// 5,034 to 4,911, below the 5,000 required bases.
    #[test]
    fn extra_hit_kept() {
        let budget = LossBudget::new(0.00005).unwrap();
        let bound = CoverageBound::new(default_model(), budget, 0.95, 10_000, false, 5_000);
        let without: Vec<u32> = (1..=3_075).collect();
        let with: Vec<u32> = (0..=3_075).collect();
        assert!(bound.bound(&with) >= bound.bound(&without));
        assert_eq!(
            bound.keeps(bound.bound(&with)),
            bound.keeps(bound.bound(&without))
        );
        // Random hit sets: one more hit never lowers the bound.
        let mut random = crate::test_support::rng(4);
        for _ in 0..200 {
            let mut hits: Vec<u32> = (0..random.random_range(0..40))
                .map(|_| random.random_range(0..9_976))
                .collect();
            hits.sort_unstable();
            hits.dedup();
            let before = bound.bound(&hits);
            hits.push(random.random_range(0..9_976));
            hits.sort_unstable();
            hits.dedup();
            assert!(bound.bound(&hits) >= before);
        }
    }

    /// The query and the contig are nine copies of a 124-base sequence with one minimizer. At
    /// identity 1 the allowance was one k-mer, so this exact match got a bound of 803 bases.
    #[test]
    fn exact_match_kept() {
        let unit = b"AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTATTTGAAGGAGTCTAGCAGCCGCAGT";
        let sequence = unit.repeat(9);
        let seeder = Seeder::new(SeedScheme::default()).unwrap();
        let mut seeds = Vec::new();
        seeder.seeds(&sequence, |kmer| seeds.push(kmer.hash));
        let query = QueryKmers::new(&seeder, &sequence, false);
        let hit: Vec<u32> = (0u32..)
            .zip(query.distinct())
            .filter(|(_, kmer)| seeds.contains(&kmer.hash))
            .map(|(seed, _)| seed)
            .collect();
        let positions = query.sorted_positions(&hit);
        let length = sequence.len() as u64;
        let budget = LossBudget::default();
        let bound = CoverageBound::new(default_model(), budget, 1.0, length, false, length);
        assert_eq!(bound.bound(&positions), length);
    }

    /// The bound of a query, at identity 0.95 and budget 0.00005, with the default k and w.
    fn proven(model: &'static MissModel, length: u64, circular: bool) -> CoverageBound<'static> {
        let budget = LossBudget::new(0.00005).unwrap();
        CoverageBound::new(model, budget, 0.95, length, circular, 0)
    }

    #[test]
    fn proven_allowances() {
        let exact = default_model();
        let one = MissModel::one_mismatch(SeedScheme::default()).unwrap();
        let one: &'static MissModel = Box::leak(Box::new(one));
        // A linear query of L bases holds L - a starts of a + 1 bases, a circular one L.
        assert_eq!(proven(exact, 5_000, false).allowance, 4_978);
        assert_eq!(proven(exact, 5_000, true).allowance, 5_000);
        assert_eq!(proven(exact, 10_000, false).allowance, 6_943);
        assert_eq!(proven(exact, 10_000, true).allowance, 7_315);
        assert_eq!(proven(exact, 50_000, false).allowance, 7_935);
        assert_eq!(proven(one, 5_000, false).allowance, 2_231);
        assert_eq!(proven(one, 10_000, false).allowance, 2_355);
        assert_eq!(proven(one, 50_000, false).allowance, 2_603);
        // One more hit never lowers a proven bound either, so a key collision cannot drop an entry.
        let bound = proven(one, 10_000, false);
        let mut random = crate::test_support::rng(5);
        let mut hits: Vec<u32> = Vec::new();
        for _ in 0..100 {
            let before = bound.bound(&hits);
            hits.push(random.random_range(0..9_976));
            hits.sort_unstable();
            hits.dedup();
            assert!(bound.bound(&hits) >= before);
        }
    }

    /// A circular exact copy of one 124-base window: its only seed leaves 99 uncovered bases,
    /// so an allowance of 75 bases gave a bound of 100.
    #[test]
    fn circular_copy_kept() {
        let unit = b"AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTATTTGAAGGAGTCTAGCAGCCGCAGT";
        let seeder = Seeder::new(SeedScheme::default()).unwrap();
        let mut seeds = Vec::new();
        seeder.seeds(unit, |kmer| seeds.push(kmer.hash));
        let query = QueryKmers::new(&seeder, unit, true);
        let hit: Vec<u32> = (0u32..)
            .zip(query.distinct())
            .filter(|(_, kmer)| seeds.contains(&kmer.hash))
            .map(|(seed, _)| seed)
            .collect();
        let positions = query.sorted_positions(&hit);
        let budget = LossBudget::default();
        let bound = CoverageBound::new(default_model(), budget, 1.0, 124, true, 124);
        assert_eq!(bound.bound(&positions), 124);
    }
}
