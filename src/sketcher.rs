use needletail::Sequence;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Sketch {
    pub name: String,
    kmer_budget: usize,
    pub hashes: BTreeSet<u64>,
    pub lowest_hash: u64,
    pub max_kmers: usize,
    pub kmer_size: u8,
}

impl Sketch {
    fn push(&mut self, kmer: &[u8]) {
        let hash = xxhash_rust::xxh3::xxh3_64(kmer);
        if self.kmer_budget > 0 {
            self.kmer_budget -= 1;
            self.hashes.insert(hash);
        } else if hash < self.lowest_hash {
            self.hashes.insert(hash);
            self.lowest_hash = hash;
            self.hashes.pop_last();
        }
    }
}

pub struct Sketcher {
    kmer_length: u8,
    current_sketch: Sketch,
    num_kmers: usize,
}

impl Sketcher {
    pub fn new(kmer_length: u8, budget: u64, name: String) -> Self {
        Sketcher {
            kmer_length,
            current_sketch: Sketch {
                kmer_budget: budget as usize,
                hashes: BTreeSet::new(),
                lowest_hash: u64::MAX,
                max_kmers: 0,
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
        let mut sketch = self.current_sketch;
        sketch.max_kmers = self.num_kmers;
        sketch
    }
}
