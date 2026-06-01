//! The coverage bound of an entry: the query length minus the stretches without a seed hit that
//! the loss budget lets count as absent.

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
        Self {
            model,
            budget,
            identity,
            query_length,
            circular,
            required_bases,
        }
    }

    /// The query coverage, in bases, that an entry must reach.
    pub fn required_bases(&self) -> u64 {
        self.required_bases
    }

    /// The highest query coverage, in bases, of an entry with seed hits (or possible seed hits)
    /// at the query starts `positions`, which must be ascending.
    pub fn bound(&self, positions: &[u32]) -> u64 {
        let stretches = self.stretches(positions);
        let absent = self
            .budget
            .absent_stretches(self.model, self.identity, &stretches);
        let absent_bases: u64 = stretches
            .iter()
            .zip(absent)
            .filter(|(_, absent)| *absent)
            .map(|(bases, _)| bases)
            .sum();
        self.query_length - absent_bases
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
    use std::sync::OnceLock;

    fn model() -> &'static MissModel {
        static MODEL: OnceLock<MissModel> = OnceLock::new();
        MODEL.get_or_init(|| MissModel::with_layouts(21, 20, 50).unwrap())
    }

    fn bound(length: u64, circular: bool) -> CoverageBound<'static> {
        CoverageBound::new(model(), LossBudget::default(), 0.95, length, circular, 0)
    }

    #[test]
    fn stretches_lie_between_and_around_the_kmers_with_seed_hits() {
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
    fn long_stretches_count_as_absent_and_short_ones_do_not() {
        let linear = bound(10_000, false);
        assert_eq!(linear.bound(&[]), 0);
        // A stretch of 40 bases may well hold the query without a seed hit; 8,000 bases cannot.
        let positions: Vec<u32> = (0..40).map(|i| 61 * i).collect();
        let end = u64::from(*positions.last().unwrap()) + 21;
        assert_eq!(linear.bound(&positions), end);
    }

    #[test]
    fn the_stretch_over_the_end_joins_both_ends() {
        // Stretches of 250 and 300 bases at the ends: each alone may hold the query without a
        // seed hit, but not both together.
        let positions: Vec<u32> = (0..90).map(|i| 250 + 61 * i).collect();
        assert_eq!(bound(6_000, false).bound(&positions), 6_000);
        let circular = bound(6_000, true);
        assert_eq!(circular.stretches(&positions).last(), Some(&550));
        assert_eq!(circular.bound(&positions), 6_000 - 550);
    }
}
