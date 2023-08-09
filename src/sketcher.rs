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
    pub heap: BinaryHeap<u64>,
    pub hashes: HashSet<u64, BuildHasherDefault<NoHashHasher>>,
    pub num_kmers: usize,
    pub max_kmers: usize,
    pub kmer_size: u8,
    pub kmer_budget: u64,
}

impl Sketch {
    fn push(&mut self, kmer: &[u8]) {
        let hash = xxhash_rust::xxh3::xxh3_64(kmer);
        let add = match self.heap.peek() {
            Some(largest) => hash < *largest,
            None => true,
        };
        if add {
            if self.hashes.insert(hash) {
                self.heap.push(hash);
                self.kmer_budget -= 1;
                if self.kmer_budget == 0 {
                    self.heap.pop();
                    self.kmer_budget += 1;
                }
            }
        }
    }
}

pub struct Sketcher {
    kmer_length: u8,
    current_sketch: Sketch,
}

impl Sketcher {
    pub fn new(kmer_length: u8, budget: u64, name: String) -> Self {
        Sketcher {
            kmer_length,
            current_sketch: Sketch {
                kmer_budget: budget,
                heap: BinaryHeap::with_capacity(budget as usize),
                hashes: HashSet::with_capacity_and_hasher(
                    budget as usize,
                    BuildHasherDefault::<NoHashHasher>::default(),
                ),
                name,
                kmer_size: kmer_length,
                num_kmers: 0,
                max_kmers: 0,
            },
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
            self.current_sketch.max_kmers += 1;
        }
    }

    pub fn finalize(self) -> Sketch {
        let mut sketch: Sketch = self.current_sketch.into();
        sketch.heap.clear();
        sketch.heap.shrink_to_fit();
        sketch.num_kmers = sketch.hashes.len();
        sketch
    }
}
