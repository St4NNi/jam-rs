//! The query k-mers of one search, grouped by hash, so each distinct k-mer is looked up once.

use crate::seeds::Seeder;
use rayon::prelude::*;
use std::ops::Range;

/// One distinct query k-mer and where it occurs on the query.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct QuerySeed {
    /// jamhash of the canonical k-mer, salted below the entropy limit, as the index stores it.
    pub hash: u64,
    /// Its starts on the query, as a range into `QueryKmers::positions`.
    positions: Range<u32>,
}

/// The query k-mers of one query. On a circular query the k-mers over the end are included, so
/// every start from 0 to the query length minus 1 can hold a k-mer.
#[derive(Clone, Debug)]
pub struct QueryKmers {
    length: u64,
    circular: bool,
    /// The minimizer window of the seed scheme, in k-mers.
    window: u64,
    kmers: u64,
    seeds: Vec<QuerySeed>,
    positions: Vec<u32>,
}

impl QueryKmers {
    /// Hashes every k-mer of `query` with the seed scheme of the shard. The query must have at
    /// least k bases and fewer than 2^32.
    pub fn new(seeder: &Seeder, query: &[u8], circular: bool) -> Self {
        let k = usize::from(seeder.scheme().k);
        let mut scanned = query.to_vec();
        if circular {
            scanned.extend_from_slice(&query[..(k - 1).min(query.len())]);
        }
        let mut kmers: Vec<(u64, u32)> = Vec::with_capacity(query.len());
        seeder.query_kmers(&scanned, |kmer| {
            if kmer.position < query.len() as u64 {
                kmers.push((kmer.hash, kmer.position as u32));
            }
        });
        kmers.par_sort_unstable();
        let mut seeds: Vec<QuerySeed> = Vec::new();
        for (index, &(hash, _)) in kmers.iter().enumerate() {
            match seeds.last_mut() {
                Some(last) if last.hash == hash => last.positions.end += 1,
                _ => seeds.push(QuerySeed {
                    hash,
                    positions: index as u32..index as u32 + 1,
                }),
            }
        }
        Self {
            length: query.len() as u64,
            circular,
            window: u64::from(seeder.scheme().window),
            kmers: kmers.len() as u64,
            seeds,
            positions: kmers.into_iter().map(|(_, position)| position).collect(),
        }
    }

    /// Bases of the query.
    pub fn length(&self) -> u64 {
        self.length
    }

    /// Whether the query is circular.
    pub fn circular(&self) -> bool {
        self.circular
    }

    /// The number of query k-mers, counting repeats.
    pub fn kmers(&self) -> u64 {
        self.kmers
    }

    /// The distinct query k-mers, by ascending hash.
    pub fn seeds(&self) -> &[QuerySeed] {
        &self.seeds
    }

    /// The query starts of the distinct k-mer `seed`, ascending.
    pub fn positions(&self, seed: u32) -> &[u32] {
        let range = &self.seeds[seed as usize].positions;
        &self.positions[range.start as usize..range.end as usize]
    }

    /// The query starts with a seed hit, divided by the seeds a full exact match would show (query
    /// k-mers times the seed density 2 / (w + 1)), at most 1.
    pub fn containment_estimate(&self, hit_positions: u64) -> f64 {
        let expected = self.kmers as f64 * 2.0 / (self.window as f64 + 1.0);
        if expected <= 0.0 {
            return 0.0;
        }
        (hit_positions as f64 / expected).min(1.0)
    }

    /// The query starts of all the given distinct k-mers, ascending.
    pub fn sorted_positions(&self, seeds: &[u32]) -> Vec<u32> {
        let mut positions: Vec<u32> = seeds
            .iter()
            .flat_map(|&seed| self.positions(seed))
            .copied()
            .collect();
        positions.sort_unstable();
        positions
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seeds::SeedScheme;

    fn seeder() -> Seeder {
        Seeder::new(SeedScheme {
            k: 5,
            window: 3,
            min_entropy: 0.0,
        })
        .unwrap()
    }

    #[test]
    fn repeated_kmers_form_one_seed_with_every_start() {
        let query = b"ACGTTGCAACGTTGCA";
        let kmers = QueryKmers::new(&seeder(), query, false);
        assert_eq!(kmers.kmers(), 12);
        let all: Vec<u32> = (0..kmers.seeds().len() as u32).collect();
        assert_eq!(kmers.sorted_positions(&all), (0..12).collect::<Vec<_>>());
        let first = (0..kmers.seeds().len() as u32)
            .find(|&seed| kmers.positions(seed).contains(&0))
            .unwrap();
        assert_eq!(kmers.positions(first), &[0, 7, 8]);
    }

    #[test]
    fn a_circular_query_has_a_kmer_at_every_start() {
        let query = b"ACGTTGCAAGGT";
        let linear = QueryKmers::new(&seeder(), query, false);
        let circular = QueryKmers::new(&seeder(), query, true);
        assert_eq!(linear.kmers(), 8);
        assert_eq!(circular.kmers(), 12);
        let all: Vec<u32> = (0..circular.seeds().len() as u32).collect();
        assert_eq!(circular.sorted_positions(&all), (0..12).collect::<Vec<_>>());
    }

    #[test]
    fn bases_other_than_acgt_leave_no_kmer() {
        let query = b"ACGTTNGCAAGGT";
        let kmers = QueryKmers::new(&seeder(), query, false);
        let all: Vec<u32> = (0..kmers.seeds().len() as u32).collect();
        assert_eq!(kmers.sorted_positions(&all), vec![0, 6, 7, 8]);
    }
}
