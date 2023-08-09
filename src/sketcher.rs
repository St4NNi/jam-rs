use needletail::Sequence;
use serde::{Deserialize, Serialize};
use std::{
    collections::{BinaryHeap, HashSet},
    hash::BuildHasherDefault,
};

use crate::hasher::NoHashHasher;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Sketch {
    pub name: String,
    pub hashes: HashSet<u64, BuildHasherDefault<NoHashHasher>>,
    pub num_kmers: usize,
    pub max_kmers: usize,
    pub kmer_size: u8,
}

impl From<WipSketch> for Sketch {
    fn from(wip: WipSketch) -> Self {
        let num_kmers = wip.hashes.len();
        Sketch {
            name: wip.name,
            hashes: wip.hashes.into_iter().collect(),
            num_kmers,
            max_kmers: 0,
            kmer_size: wip.kmer_size,
        }
    }
}

pub struct WipSketch {
    kmer_budget: usize,
    hashes: BinaryHeap<u64>,
    name: String,
    kmer_size: u8,
}

impl WipSketch {
    fn push(&mut self, kmer: &[u8]) {
        let hash = xxhash_rust::xxh3::xxh3_64(kmer);
        match self.hashes.peek() {
            Some(largest) => {
                if hash < *largest {
                    self.hashes.push(hash);
                    self.kmer_budget -= 1;
                    if self.kmer_budget == 0 {
                        self.hashes.pop();
                    }
                }
            }
            None => {
                self.hashes.push(hash);
                self.kmer_budget -= 1;
            }
        }
    }
}

pub struct Sketcher {
    kmer_length: u8,
    current_sketch: WipSketch,
    num_kmers: usize,
}

impl Sketcher {
    pub fn new(kmer_length: u8, budget: u64, name: String) -> Self {
        Sketcher {
            kmer_length,
            current_sketch: WipSketch {
                kmer_budget: budget as usize,
                hashes: BinaryHeap::with_capacity(budget as usize),
                name,
                kmer_size: kmer_length,
            },
            num_kmers: 0,
        }
    }
}

impl Sketcher {
    // This is more or less derived from the `process` method in `finch-rs`:
    // https://github.com/onecodex/finch-rs/blob/master/lib/src/sketch_schemes/mash.rs
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
        let mut sketch: Sketch = self.current_sketch.into();
        sketch.max_kmers = self.num_kmers;
        sketch
    }
}
