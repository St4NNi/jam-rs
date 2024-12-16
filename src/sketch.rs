use itertools::Itertools;
use serde::{Deserialize, Serialize};
use sourmash::sketch::{minhash::KmerMinHash, Sketch as SourmashSketch};
use std::collections::BTreeSet;

#[derive(Debug, Serialize, Deserialize, Clone, Default)]
pub struct Sketch {
    pub name: String,          // Name of file or sequence
    pub hashes: BTreeSet<u64>, // Hashes with stats
    pub num_kmers: usize,      // Number of kmers (collected)
    pub kmer_size: u8,         // Kmer size
}

impl Sketch {
    pub fn new(name: String, num_kmers: usize, kmer_size: u8) -> Self {
        Sketch {
            name,
            num_kmers,
            kmer_size,
            hashes: BTreeSet::new(),
        }
    }
}

impl Sketch {
    pub fn into_sourmash(self, max_hash: u64) -> SourmashSketch {
        let sketch = KmerMinHash::builder()
            .ksize(self.kmer_size as u32)
            .num(self.hashes.len() as u32)
            .max_hash(max_hash)
            .mins(self.hashes.into_iter().sorted().collect::<Vec<u64>>())
            .build();
        SourmashSketch::MinHash(sketch)
    }
}
