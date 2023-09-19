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

    pub fn get_size_class(size: usize) -> u8 {
        if size > 2000 * u8::MAX as usize {
            return u8::MAX;
        }
        (size / 2000) as u8
    }

    pub fn get_gc_class(gc_size: usize, size: usize) -> u8 {
        (gc_size * 256 / size * 256) as u8
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
