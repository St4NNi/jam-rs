//! The query k-mers of one search, grouped by hash, so each distinct k-mer is looked up once.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::seed_mismatches::LookupPlan;
use crate::seeds::{Kmer, Seeder};
use rayon::prelude::*;
use std::convert::Infallible;
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
    lookups: u64,
    distinct: Vec<DistinctKmer>,
    positions: Vec<u32>,
}

impl QueryKmers {
    /// Hashes every k-mer of `query` with the seed scheme of the shard. The query must have at
    /// least k bases and fewer than 2^32.
    pub fn new(seeder: &Seeder, query: &[u8], circular: bool) -> Self {
        Self::planned(seeder, query, circular, &LookupPlan::full(0))
    }

    /// As `new`, adding every k-mer one substitution away from a query k-mer at its start, so
    /// seeds with one mismatch are found too.
    pub fn with_neighbours(seeder: &Seeder, query: &[u8], circular: bool) -> Self {
        Self::planned(seeder, query, circular, &LookupPlan::full(1))
    }

    /// As `with_neighbours`, adding every k-mer two substitutions away too.
    pub fn with_two_neighbours(seeder: &Seeder, query: &[u8], circular: bool) -> Self {
        Self::planned(seeder, query, circular, &LookupPlan::full(2))
    }

    /// As `with_two_neighbours`, adding every k-mer one indel and at most one substitution away.
    pub fn with_indel_neighbours(seeder: &Seeder, query: &[u8], circular: bool) -> Self {
        Self::planned(seeder, query, circular, &LookupPlan::full(3))
    }

    /// Hashes every k-mer of `query` with the lookups of `plan`: those of its full level at every
    /// query position, and those of the next level at its sampled positions.
    pub fn planned(seeder: &Seeder, query: &[u8], circular: bool, plan: &LookupPlan) -> Self {
        let kept = Self::filtered(seeder, query, circular, plan, |_| Ok::<(), Infallible>(()));
        match kept {
            Ok(kmers) => kmers,
            Err(never) => match never {},
        }
    }

    /// As `planned`, making the lookups in batches of query positions, each passed to `keep`,
    /// which may drop lookups (those the index does not hold), so only the kept ones are held.
    pub fn filtered<E>(
        seeder: &Seeder,
        query: &[u8],
        circular: bool,
        plan: &LookupPlan,
        mut keep: impl FnMut(&mut Vec<(u64, u32)>) -> Result<(), E>,
    ) -> Result<Self, E> {
        let k = usize::from(seeder.scheme().k);
        let mut scanned = query.to_vec();
        if circular {
            scanned.extend_from_slice(&query[..(k - 1).min(query.len())]);
        }
        let mut exact = 0;
        seeder.query_kmers(&scanned, |kmer| {
            exact += u64::from(kmer.position < query.len() as u64);
        });
        let per_position = LOOKUPS_PER_POSITION[plan.level()];
        let batch = (BATCH_LOOKUPS / per_position).max(1) as usize;
        let (mut kmers, mut lookups) = (Vec::new(), 0);
        for first in (0..query.len()).step_by(batch) {
            let positions = first..(first + batch).min(query.len());
            let mut made: Vec<(u64, u32)> = (positions.into_par_iter())
                .flat_map_iter(|position| lookups_at(seeder, &scanned, plan, position))
                .collect();
            lookups += made.len() as u64;
            keep(&mut made)?;
            kmers.append(&mut made);
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
        Ok(Self {
            length: query.len() as u64,
            window: u64::from(seeder.scheme().w),
            kmers: exact,
            lookups,
            distinct,
            positions: kmers.into_iter().map(|(_, position)| position).collect(),
        })
    }

    /// The lookups made, counting repeats: query k-mers, their neighbours and indel variants.
    pub fn lookups(&self) -> u64 {
        self.lookups
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

    /// The query starts of all the given distinct k-mers, ascending. A start that several
    /// k-mers share (a k-mer and its one-substitution neighbours) is listed once.
    pub fn sorted_positions(&self, seeds: &[u32]) -> Vec<u32> {
        let mut positions: Vec<u32> = seeds
            .iter()
            .flat_map(|&seed| self.positions(seed))
            .copied()
            .collect();
        positions.sort_unstable();
        positions.dedup();
        positions
    }
}

/// The most lookups that one batch of `QueryKmers::filtered` makes, about 64 MiB.
const BATCH_LOOKUPS: u64 = 1 << 22;

/// The lookups per query position, by the highest level of a plan, at k 25: 1, 1 + 75,
/// 1 + 75 + 2,700, and that plus 116 indel variants with their 75 substitutions each.
const LOOKUPS_PER_POSITION: [u64; 4] = [1, 76, 2_776, 2_776 + 116 * 76];

/// The lookups that `plan` makes at query position `position`, each on the bases from there alone
/// (k, or k + 1 for the indel lookups, fewer at the end), so a position gets exactly the lookups
/// that enumerating the whole query gives it. An indel lookup also starts at the last k - 1 bases.
fn lookups_at(
    seeder: &Seeder,
    scanned: &[u8],
    plan: &LookupPlan,
    position: usize,
) -> Vec<(u64, u32)> {
    let k = usize::from(seeder.scheme().k);
    let bases = &scanned[position..(position + k + 1).min(scanned.len())];
    let level = plan.full + usize::from(plan.is_sampled(position as u64));
    let mut made = Vec::new();
    let mut push = |kmer: Kmer| {
        if kmer.position == 0 {
            made.push((kmer.hash, position as u32));
        }
    };
    if bases.len() >= k {
        seeder.query_kmers(&bases[..k], &mut push);
        if level >= 1 {
            seeder.query_neighbours(&bases[..k], &mut push);
        }
        if level >= 2 {
            seeder.query_neighbours_two(&bases[..k], &mut push);
        }
    }
    if level >= 2 && plan.indels {
        seeder.query_neighbours_indel_at(bases, true, &mut push);
    }
    made
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seeds::SeedScheme;
    use crate::test_support::{random_bases, rng};

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

    #[test]
    fn shared_starts_once() {
        let query = b"ACGTTGCAAGGTCATG";
        let kmers = QueryKmers::with_neighbours(&seeder(), query, false);
        let all: Vec<u32> = (0..kmers.distinct().len() as u32).collect();
        let total: usize = all.iter().map(|&seed| kmers.positions(seed).len()).sum();
        assert!(total > kmers.kmers() as usize);
        assert_eq!(kmers.sorted_positions(&all), (0..12).collect::<Vec<_>>());
    }

    #[test]
    fn indel_lookups_added() {
        let seeder = seeder();
        let pairs = |kmers: &QueryKmers| -> usize {
            let all = 0..kmers.distinct().len() as u32;
            all.map(|seed| kmers.positions(seed).len()).sum()
        };
        // A k-mer of 15 query bases with the eighth of 16 removed: one indel, many substitutions.
        let wide = Seeder::new(SeedScheme {
            k: 15,
            w: 10,
            min_entropy: 0.0,
        })
        .unwrap();
        let query = random_bases(&mut rng(4), 40);
        let mut removed = Vec::new();
        wide.query_kmers(&[&query[..7], &query[8..16]].concat(), |kmer| {
            removed.push(kmer.hash)
        });
        let kmers = QueryKmers::with_indel_neighbours(&wide, &query, false);
        assert!(kmers.distinct().iter().any(|kmer| kmer.hash == removed[0]));
        // Two neighbours are substitutions only: 1 + 15 + 90 lookups per k-mer of 5 bases.
        let long = random_bases(&mut rng(3), 1200);
        let two = QueryKmers::with_two_neighbours(&seeder, &long, false);
        assert_eq!(pairs(&two), 1196 * 106);
        let indel = QueryKmers::with_indel_neighbours(&seeder, &long, false);
        assert!(pairs(&indel) > 1196 * 106);
    }

    /// The pairs of query start and hash, sorted.
    fn pairs(kmers: &QueryKmers) -> Vec<(u32, u64)> {
        let mut pairs: Vec<(u32, u64)> = (0..kmers.distinct().len() as u32)
            .flat_map(|seed| kmers.positions(seed).iter().map(move |&at| (at, seed)))
            .map(|(at, seed)| (at, kmers.distinct()[seed as usize].hash))
            .collect();
        pairs.sort_unstable();
        pairs
    }

    /// A sampled position gets exactly the lookups that the whole query gives it at the higher
    /// level, the terminal indel lookups of the last k - 1 bases too; the others those of the
    /// full level.
    #[test]
    fn sampled_like_full() {
        let seeder = seeder();
        let query = random_bases(&mut rng(7), 300);
        let mut plan = LookupPlan {
            full: 1,
            sampled: 150,
            starts: 296,
            indels: true,
            proven: false,
            expected_hits: 0.0,
        };
        // A share near half that samples position 296, which holds only k - 1 bases.
        while !plan.is_sampled(296) {
            plan.sampled += 1;
        }
        let planned = pairs(&QueryKmers::planned(&seeder, &query, false, &plan));
        let high = pairs(&QueryKmers::with_indel_neighbours(&seeder, &query, false));
        let low = pairs(&QueryKmers::with_neighbours(&seeder, &query, false));
        let mut expected: Vec<(u32, u64)> = (high.iter())
            .filter(|(at, _)| plan.is_sampled(u64::from(*at)))
            .chain(
                low.iter()
                    .filter(|(at, _)| !plan.is_sampled(u64::from(*at))),
            )
            .copied()
            .collect();
        expected.sort_unstable();
        assert_eq!(planned, expected);
        assert!(planned.iter().any(|&(at, _)| at == 296));
    }

    /// Lookups dropped batch by batch leave the same k-mers as dropping them at the end.
    #[test]
    fn filtered_batches_keep() {
        let seeder = seeder();
        let query = random_bases(&mut rng(8), 3_000);
        let plan = LookupPlan::full(2);
        let even = |pairs: &mut Vec<(u64, u32)>| {
            pairs.retain(|&(hash, _)| hash % 2 == 0);
            Ok::<(), ()>(())
        };
        let kept = QueryKmers::filtered(&seeder, &query, false, &plan, even).unwrap();
        let all = QueryKmers::planned(&seeder, &query, false, &plan);
        let mut expected = pairs(&all);
        expected.retain(|&(_, hash)| hash % 2 == 0);
        assert_eq!(pairs(&kept), expected);
        assert_eq!((kept.kmers(), kept.lookups()), (all.kmers(), all.lookups()));
    }
}
