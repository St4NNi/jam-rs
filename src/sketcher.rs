use crate::{
    cli::HashAlgorithms,
    hash_functions::Function,
    signature::Signature,
    sketch::Sketch,
};
use needletail::{parser::SequenceRecord, Sequence};
use std::{
    collections::BinaryHeap, fs::File
};

pub enum Storage {
    Sourmash(File),
    Lmdb(heed::Env),
}

#[derive(Debug, Default)]
struct SketchHelper {
    pub max_hash: u64,
    hit_counter: u64,
    kmer_seq_counter: u64,
    pub nmax: u64,
    pub heap: BinaryHeap<u64>,
}

impl SketchHelper {
    pub fn new(max_hash: u64, nmax: Option<u64>) -> Self {
        SketchHelper {
            nmax: nmax.unwrap_or(u64::MAX),
            hit_counter: 0,
            kmer_seq_counter: 0,
            max_hash,
            heap: BinaryHeap::new(),
        }
    }

    pub fn push(&mut self, hash: u64) {
        // Increase the local sequence counter in any case
        self.kmer_seq_counter += 1;
        if hash < self.max_hash {
            self.hit_counter += 1;
            self.heap.push(hash);
            if self.heap.len() > self.nmax as usize {
                self.heap.pop();
            }
        }
    }

    pub fn reset(&mut self) {
        let nmax = self.nmax;
        *self = Self::default();
        self.nmax = nmax;
    }

    pub fn into_sketch(&mut self, name: String, kmer_size: u8) -> Sketch {
        let mut sketch = Sketch::new(
            name,
            self.heap.len(),
            kmer_size,
        );
        sketch.hashes = self.heap.drain().collect();
        self.reset();
        sketch
    }
}

pub struct Sketcher<'a> {
    name: String,
    kmer_length: u8,
    helper: SketchHelper,
    completed_sketches: Vec<Sketch>,
    singleton: bool,
    function: Function<'a>,
    algorithm: HashAlgorithms,    
}

impl<'a> Sketcher<'a> {
    pub fn new(
        kmer_length: u8,
        name: String,
        singleton: bool,
        max_hash: u64,
        nmax: Option<u64>,
        function: Function<'a>,
        algorithm: HashAlgorithms,
    ) -> Self {
        Sketcher {
            name,
            kmer_length,
            helper: SketchHelper::new(max_hash, nmax),
            singleton,
            completed_sketches: Vec::new(),
            function,
            algorithm,
        }
    }
}

impl Sketcher<'_> {
    // This is more or less derived from the `process` method in `finch-rs`:
    // https://github.com/onecodex/finch-rs/blob/master/lib/src/sketch_schemes/mash.rs
    pub fn process<'seq, 'a, 'inner>(&'a mut self, seq: &'seq SequenceRecord<'inner>)
    where
        'a: 'seq,
        'seq: 'inner,
    {
        let name = seq.id();
        let seq = seq.normalize(false);
        if self.kmer_length <= 31 {
            let func_small = self.function.get_small().unwrap();
            for (_, kmer, _) in seq.bit_kmers(self.kmer_length, true) {
                self.helper.push(func_small(kmer.0));
            }
        } else {
            let func_large = self.function.get_large().unwrap();
            let rc = seq.reverse_complement();
            for (_, kmer, _) in seq.canonical_kmers(self.kmer_length, &rc) {
                self.helper.push(func_large(kmer));
            }
        }
        if self.singleton {
            self.completed_sketches.push(
                self.helper
                    .into_sketch(String::from_utf8_lossy(name).to_string(), self.kmer_length),
            );
        }
    }

    pub fn finish(self) -> Signature {
        let max_hash = self.helper.max_hash;
        let file_name = self.name.to_string();
        let algorithm = self.algorithm.clone();
        let kmer_size = self.kmer_length;
        let mut sketches = self.completed_sketches;
        let mut helper = self.helper;
        sketches.push(helper.into_sketch(self.name, self.kmer_length));
        Signature {
            file_name,
            sketches,
            max_hash,
            algorithm,
            kmer_size,
        }
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_sketch_helper() {
//         let mut helper = SketchHelper::new(1, 100, None, None);
//         helper.initialize_record(Some(Stats::new(0, 0)));
//         helper.push(1);
//         helper.push(2);
//         helper.push(3);
//         assert_eq!(
//             helper.into_sketch("sketch".to_string(), 1),
//             Sketch {
//                 name: "sketch".to_string(),
//                 hashes: HashMap::from_iter(vec![(1, Some(Stats::new(0, 0)))]),
//                 num_kmers: 1,
//                 max_kmers: 3,
//                 kmer_size: 1
//             }
//         );
//     }
// }
