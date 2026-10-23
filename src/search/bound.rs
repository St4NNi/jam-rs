//! The coverage bound of an entry: the most query bases it can still cover, given where its seed
//! hits lie on the query.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::loss::{MaxMiss, MissModel};

/// Computes coverage bounds for one query against one shard.
#[derive(Clone, Copy, Debug)]
pub struct CoverageBound<'a> {
    model: &'a MissModel,
    /// The identity at which the proven bound holds for a qualifying copy: below the identity
    /// limit, see `MissModel::qualifying_identity`.
    identity: f64,
    query_length: u64,
    circular: bool,
    /// The query coverage, in bases, that an entry must reach.
    required_bases: u64,
    /// The present bases that a stretch without a seed hit may hold, for this query and these
    /// settings; it never depends on the hits, so a further hit can only raise a bound. The query
    /// length turns the pruning off.
    allowance: u64,
    /// The chance that the pruning may spend: --max-miss less the seed miss bound.
    pruning_share: f64,
    /// The distinct hit starts that an entry needs.
    min_hits: u64,
    /// The proven chance that a present piece of the required length shows fewer than
    /// `min_hits` hit starts.
    seed_miss: f64,
}

impl<'a> CoverageBound<'a> {
    /// A bound for a query of `query_length` bases. The seed miss bound and the pruning together
    /// stay within `max_miss`; when the seed miss bound alone reaches it, the pruning is off.
    pub fn new(
        model: &'a MissModel,
        max_miss: MaxMiss,
        identity: f64,
        query_length: u64,
        circular: bool,
        required_bases: u64,
    ) -> Self {
        let identity = model.qualifying_identity(identity, required_bases.min(query_length));
        Self {
            model,
            identity,
            query_length,
            circular,
            required_bases,
            allowance: query_length,
            pruning_share: 0.0,
            min_hits: 1,
            seed_miss: model.miss_upper_bound(required_bases.min(query_length), identity),
        }
        .spend(max_miss)
    }

    /// The bound when an entry needs `min_hits` distinct hit starts, `seed_miss` being the proven
    /// chance that a present piece of the required length shows fewer. The pruning gets the rest
    /// of `max_miss`.
    pub fn with_min_hits(mut self, min_hits: u64, seed_miss: f64, max_miss: MaxMiss) -> Self {
        (self.min_hits, self.seed_miss, self.allowance) = (min_hits, seed_miss, self.query_length);
        self.spend(max_miss)
    }

    /// Gives the pruning what `max_miss` leaves after the seed miss bound; none turns it off.
    fn spend(mut self, max_miss: MaxMiss) -> Self {
        self.pruning_share = max_miss.chance() - self.seed_miss;
        if self.pruning_share > 0.0 {
            self.allowance = self.proven_allowance();
        }
        self
    }

    /// The distinct hit starts that an entry needs.
    pub fn min_hits(&self) -> u64 {
        self.min_hits
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

    /// The proven bound on the chance that a present piece of the required length shows fewer hit
    /// starts than an entry needs (no seed hit when one is enough): the chance to lose a
    /// qualifying entry before any bound sees it.
    pub fn seed_miss_bound(&self) -> f64 {
        self.seed_miss
    }

    /// Whether the coverage bound can drop an entry: false when the allowance is the whole query.
    pub fn prunes(&self) -> bool {
        self.allowance < self.query_length
    }

    /// The modeled total chance of missing a qualifying entry: the seed miss bound, plus the
    /// chance the pruning spends, `n * U(a + 1)`, when it is on.
    pub fn total_miss(&self) -> f64 {
        let pruning = match self.prunes() {
            true => self.starts(self.allowance) as f64 * self.stretch_miss(self.allowance),
            false => 0.0,
        };
        (self.seed_miss_bound() + pruning).min(1.0)
    }

    /// The starts of a seedless stretch of a + 1 bases: L on a circular query, else L - a.
    fn starts(&self, allowance: u64) -> u64 {
        match self.circular {
            true => self.query_length,
            false => self.query_length.saturating_sub(allowance).max(1),
        }
    }

    /// The proven chance that a present piece of a + 1 bases shows no seed hit.
    fn stretch_miss(&self, allowance: u64) -> f64 {
        self.model.miss_upper_bound(allowance + 1, self.identity)
    }

    /// The shortest length `a` with `n * U(a + 1) <= B`, or L: B the pruning share, U the proven
    /// bound and n the starts of a seedless stretch of a + 1 bases. It holds when each stretch
    /// meets at most one coherent present piece.
    fn proven_allowance(&self) -> u64 {
        let fits = |allowance: u64| {
            let miss = self.stretch_miss(allowance);
            self.starts(allowance) as f64 * miss <= self.pruning_share
        };
        smallest_fit(self.query_length, fits)
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
        CoverageBound::new(model(), MaxMiss::default(), 0.95, length, circular, length)
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
            (10_000 - allowance) as f64 * miss <= 1e-4 - linear.seed_miss_bound()
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
        MODEL.get_or_init(|| MissModel::new(SeedScheme::default(), 0).unwrap())
    }

    /// Hit starts 1 to 3,075 leave stretches of 1 and 6,900 bases; a hit at 0 removes the first.
    /// With an allowance from the number of stretches, the extra hit lowered the bound from
    /// 5,034 to 4,911, below the 5,000 required bases.
    #[test]
    fn extra_hit_kept() {
        let budget = MaxMiss::new(0.00005).unwrap();
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
        let budget = MaxMiss::default();
        let bound = CoverageBound::new(default_model(), budget, 1.0, length, false, length);
        assert_eq!(bound.bound(&positions), length);
    }

    /// The bound of a query that must be covered whole, at identity 0.95 and --max-miss
    /// 0.00005, with the default k and w.
    fn proven(model: &'static MissModel, length: u64, circular: bool) -> CoverageBound<'static> {
        let max_miss = MaxMiss::new(0.00005).unwrap();
        CoverageBound::new(model, max_miss, 0.95, length, circular, length)
    }

    /// When the seed miss bound alone reaches --max-miss, the pruning is off: every bound is the
    /// whole query, and the total is the seed miss bound.
    #[test]
    fn seed_miss_unpruned() {
        let max_miss = MaxMiss::default();
        // Half of 1,000 bases holds 4 blocks of 124 bases: U = (1 - 0.95^25)^4, about 0.27.
        let off = CoverageBound::new(default_model(), max_miss, 0.95, 1_000, false, 500);
        assert!(!off.prunes());
        assert_eq!(off.bound(&[]), 1_000);
        assert!(off.seed_miss_bound() > 0.2);
        assert_eq!(off.total_miss(), off.seed_miss_bound());
        // A whole 10 kb query leaves room: the pruning is on and the total stays within.
        let on = CoverageBound::new(default_model(), max_miss, 0.95, 10_000, false, 10_000);
        assert!(on.prunes());
        assert!(on.total_miss() <= max_miss.chance());
        assert!(on.total_miss() > on.seed_miss_bound());
    }

    /// SecY (1,299 bases, 650 required) at identity 0.95 with level 2 and two hit starts: the
    /// seed part 0.000454314 plus the pruning of the 0.0001 allowance of 991 bases is 0.000494765
    /// (the design study); at --max-miss 0.001 the pruning gets the rest, an allowance of 867
    /// bases, for a total of 0.000865524.
    #[test]
    fn two_hit_totals() {
        let two = MissModel::new(SeedScheme::default(), 2).unwrap();
        let two: &'static MissModel = Box::leak(Box::new(two));
        let seed = two.hit_miss_bound(0.95, 650, 1_299, 2);
        let at = |max_miss| CoverageBound::new(two, max_miss, 0.95, 1_299, false, 650);
        let close = |got: f64, want: f64| (got / want - 1.0).abs() < 1e-5;
        let old = at(MaxMiss::new(0.0001).unwrap());
        assert_eq!((old.allowance, old.min_hits()), (991, 1));
        let combined = seed + old.total_miss() - old.seed_miss_bound();
        assert!(close(combined, 0.000_494_765), "{combined}");
        let budget = MaxMiss::new(0.001).unwrap();
        let new = at(budget).with_min_hits(2, seed, budget);
        assert_eq!((new.allowance, new.min_hits()), (867, 2));
        assert_eq!(new.seed_miss_bound(), seed);
        assert!(
            close(new.total_miss(), 0.000_865_524),
            "{}",
            new.total_miss()
        );
        // When the seed part alone reaches the target the pruning is off.
        let tight = at(budget).with_min_hits(2, 0.001, budget);
        assert!(!tight.prunes() && tight.total_miss() == 0.001);
    }

    #[test]
    fn proven_allowances() {
        let exact = default_model();
        let one = MissModel::new(SeedScheme::default(), 1).unwrap();
        let one: &'static MissModel = Box::leak(Box::new(one));
        // A linear query of L bases holds L - a starts of a + 1 bases, a circular one L.
        assert_eq!(proven(exact, 5_000, false).allowance, 4_982);
        assert_eq!(proven(exact, 5_000, true).allowance, 5_000);
        assert_eq!(proven(exact, 10_000, false).allowance, 6_943);
        assert_eq!(proven(exact, 10_000, true).allowance, 7_439);
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
        let budget = MaxMiss::default();
        let bound = CoverageBound::new(default_model(), budget, 1.0, 124, true, 124);
        assert_eq!(bound.bound(&positions), 124);
    }
}
