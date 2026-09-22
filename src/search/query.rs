//! The query k-mers of one search, grouped by hash, so each distinct k-mer is looked up once.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::seeds::Seeder;
use rayon::prelude::*;
use std::ops::Range;

/// One distinct query k-mer and where it occurs on the query.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DistinctKmer {
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
    /// The minimizer window of the seed scheme, in k-mers.
    window: u64,
    kmers: u64,
    distinct: Vec<DistinctKmer>,
    positions: Vec<u32>,
}

impl QueryKmers {
    /// Hashes every k-mer of `query` with the seed scheme of the shard. The query must have at
    /// least k bases and fewer than 2^32.
    pub fn new(seeder: &Seeder, query: &[u8], circular: bool) -> Self {
        Self::hashed(seeder, query, circular, false)
    }

    /// As `new`, adding every k-mer one substitution away from a query k-mer at its start, so
    /// seeds with one mismatch are found too.
    pub fn with_neighbours(seeder: &Seeder, query: &[u8], circular: bool) -> Self {
        Self::hashed(seeder, query, circular, true)
    }

    fn hashed(seeder: &Seeder, query: &[u8], circular: bool, neighbours: bool) -> Self {
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
        let exact = kmers.len() as u64;
        if neighbours {
            seeder.query_neighbours(&scanned, |kmer| {
                if kmer.position < query.len() as u64 {
                    kmers.push((kmer.hash, kmer.position as u32));
                }
            });
        }
        kmers.par_sort_unstable();
        let mut distinct: Vec<DistinctKmer> = Vec::new();
        for (index, &(hash, _)) in kmers.iter().enumerate() {
            match distinct.last_mut() {
                Some(last) if last.hash == hash => last.positions.end += 1,
                _ => distinct.push(DistinctKmer {
                    hash,
                    positions: index as u32..index as u32 + 1,
                }),
            }
        }
        Self {
            length: query.len() as u64,
            window: u64::from(seeder.scheme().w),
            kmers: exact,
            distinct,
            positions: kmers.into_iter().map(|(_, position)| position).collect(),
        }
    }

    /// Bases of the query.
    pub fn length(&self) -> u64 {
        self.length
    }

    /// The number of query k-mers, counting repeats.
    pub fn kmers(&self) -> u64 {
        self.kmers
    }

    /// The distinct query k-mers, by ascending hash.
    pub fn distinct(&self) -> &[DistinctKmer] {
        &self.distinct
    }

    /// The query starts of the distinct k-mer `seed`, ascending.
    pub fn positions(&self, seed: u32) -> &[u32] {
        let range = &self.distinct[seed as usize].positions;
        &self.positions[range.start as usize..range.end as usize]
    }

    /// The containment estimate: the share of the query held by the entry, as the query starts
    /// with a seed hit over the seed hits of a full exact match (query k-mers times 2 / (w + 1)),
    /// at most 1. Mutations remove seed hits, so they lower it.
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
            w: 3,
            min_entropy: 0.0,
        })
        .unwrap()
    }

    #[test]
    fn repeats_share_kmer() {
        let query = b"ACGTTGCAACGTTGCA";
        let kmers = QueryKmers::new(&seeder(), query, false);
        assert_eq!(kmers.kmers(), 12);
        let all: Vec<u32> = (0..kmers.distinct().len() as u32).collect();
        assert_eq!(kmers.sorted_positions(&all), (0..12).collect::<Vec<_>>());
        let first = (0..kmers.distinct().len() as u32)
            .find(|&seed| kmers.positions(seed).contains(&0))
            .unwrap();
        assert_eq!(kmers.positions(first), &[0, 7, 8]);
    }

    #[test]
    fn circular_kmers_wrap() {
        let query = b"ACGTTGCAAGGT";
        let linear = QueryKmers::new(&seeder(), query, false);
        let circular = QueryKmers::new(&seeder(), query, true);
        assert_eq!(linear.kmers(), 8);
        assert_eq!(circular.kmers(), 12);
        let all: Vec<u32> = (0..circular.distinct().len() as u32).collect();
        assert_eq!(circular.sorted_positions(&all), (0..12).collect::<Vec<_>>());
    }

    #[test]
    fn unknown_bases_skipped() {
        let query = b"ACGTTNGCAAGGT";
        let kmers = QueryKmers::new(&seeder(), query, false);
        let all: Vec<u32> = (0..kmers.distinct().len() as u32).collect();
        assert_eq!(kmers.sorted_positions(&all), vec![0, 6, 7, 8]);
    }
}
