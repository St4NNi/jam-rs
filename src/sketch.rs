use crate::hasher::NoHashHasher;
use serde::{Deserialize, Serialize};
use sourmash::sketch::{minhash::KmerMinHash, Sketch as SourmashSketch};
use std::{
    collections::{BTreeSet, HashMap},
    hash::BuildHasherDefault,
};

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
#[serde(transparent)]
pub struct Stats((u8, u8));

impl Stats {
    pub fn new(size_class: u8, gc_class: u8) -> Self {
        Stats((size_class, gc_class))
    }

    pub fn from_seq(seq: &[u8]) -> Self {
        let mut gc_count = 0;
        for base in seq {
            if *base == b'G' || *base == b'C' {
                gc_count += 1;
            }
        }
        Stats((
            Self::get_size_class(seq.len()),
            Self::get_gc_class(gc_count, seq.len()),
        ))
    }

    // 255 size classes a 2000 bp
    // 0 - 512.000 bp classes
    pub fn get_size_class(size: usize) -> u8 {
        if size > 2000 * u8::MAX as usize {
            return u8::MAX;
        }
        (size / 2000) as u8
    }

    // GC class are 100 classes from 0 - 100%
    pub fn get_gc_class(gc_size: usize, size: usize) -> u8 {
        let result = (gc_size * 100 / size * 100) / 100;
        result as u8
    }

    #[inline(always)]
    pub fn compare(&self, other: &Self, gc_bounds: Option<(u8, u8)>) -> bool {
        // Expect that gc_upper_range and gc_lower_range are smaller than 100
        if let Some((gc_lower_range, gc_upper_range)) = gc_bounds {
            // If size class of self is larger than other
            if self.0 .0 >= other.0 .0 {
                // Example: Self = 55, other = 50, gc_upper_range = 5, gc_lower_range = 5
                // Self + 5 >= other
                // Self - 5 >= other --> possible hit
                if self.0 .1 + gc_upper_range >= other.0 .1
                    && self.0 .1 - gc_lower_range <= other.0 .1
                {
                    // Other is within range of self -> it is a subset
                    return true;
                }
            }
        } else {
            return self.0 .0 >= other.0 .0;
        }
        false
    }
}

impl PartialOrd for Stats {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.0.cmp(&other.0))
    }
}

#[derive(Debug, Serialize, Deserialize, Clone, Default, PartialEq)]
pub struct Sketch {
    pub name: String, // Name of file or sequence
    pub hashes: HashMap<u64, Option<Stats>, BuildHasherDefault<NoHashHasher>>, // Hashes with stats
    pub num_kmers: usize, // Number of kmers (collected)
    pub max_kmers: usize, // Max number of kmers (budget)
    pub kmer_size: u8, // Kmer size
}

impl Sketch {
    pub fn new(name: String, num_kmers: usize, max_kmers: usize, kmer_size: u8) -> Self {
        Sketch {
            name,
            max_kmers,
            num_kmers,
            kmer_size,
            hashes: HashMap::with_capacity_and_hasher(
                1_000_000,
                BuildHasherDefault::<NoHashHasher>::default(),
            ),
        }
    }
}

impl Sketch {
    pub fn into_sourmash(self, max_hash: u64) -> SourmashSketch {
        let sketch = KmerMinHash::builder()
            .ksize(self.kmer_size as u32)
            .num(self.hashes.len() as u32)
            .max_hash(max_hash)
            .mins(
                self.hashes
                    .into_keys()
                    .collect::<BTreeSet<u64>>()
                    .into_iter()
                    .collect::<Vec<u64>>(),
            )
            .build();

        SourmashSketch::MinHash(sketch)
    }
}

#[cfg(test)]
mod tests {
    use super::Stats;

    #[test]
    fn test_stats() {
        let stat = Stats::from_seq(b"ATGC");
        assert_eq!(stat, Stats::new(0, 50));

        let stata = Stats::new(2, 50);
        let statb = Stats::new(2, 50);
        assert!(stata.compare(&statb, Some((5, 5))));

        let stata = Stats::new(2, 20);
        let statb = Stats::new(4, 50);
        // Should be false because size class of statb is larger than stata
        // and gc class of statb is not within range of stata
        assert!(!stata.compare(&statb, Some((5, 5))));
        // Should fail only because of size class
        assert!(!stata.compare(&statb, Some((30, 5))));
        // Should succeed because of size class and gc class
        assert!(statb.compare(&stata, Some((30, 30))));
    }
}
