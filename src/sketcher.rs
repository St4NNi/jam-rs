use crate::hasher::NoHashHasher;
use needletail::Sequence;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::hash::BuildHasherDefault;

#[derive(Debug, Serialize, Deserialize)]
pub struct Sketch {
    kmer_budget: usize,
    pub hashes: HashSet<u64, BuildHasherDefault<NoHashHasher>>,
    pub lowest_hash: u64,
    pub max_kmers: usize,
}

impl Sketch {
    fn push(&mut self, kmer: &[u8]) {
        let hash = xxhash_rust::xxh3::xxh3_64(kmer);
        if self.kmer_budget > 0 {
            self.kmer_budget -= 1;
            self.hashes.insert(hash);
            if hash < self.lowest_hash {
                self.lowest_hash = hash;
            }
        } else if hash < self.lowest_hash {
            self.hashes.insert(hash);
            self.lowest_hash = hash;
        }
    }
}

pub struct Sketcher {
    kmer_length: u8,
    current_sketch: Sketch,
    num_kmers: usize,
}

impl Sketcher {
    pub fn new(kmer_length: u8, budget: u64) -> Self {
        Sketcher {
            kmer_length,
            current_sketch: Sketch {
                kmer_budget: budget as usize,
                hashes: HashSet::default(),
                lowest_hash: u64::MAX,
                max_kmers: 0,
            },
            num_kmers: 0,
        }
    }
}

impl Sketcher {
    pub fn process<'seq, 'a, 'inner>(&'a mut self, seq: &'seq dyn Sequence<'inner>)
    where
        'a: 'seq,
        'seq: 'inner,
    {
        let rc = seq.reverse_complement();
        for (_, kmer, _) in seq.normalize(false).canonical_kmers(self.kmer_length, &rc) {
            self.current_sketch.push(kmer);
            self.num_kmers += 1;
        }
    }

    pub fn finalize(self) -> Sketch {
        let mut sketch = self.current_sketch;
        sketch.max_kmers = self.num_kmers;
        sketch
    }
}
