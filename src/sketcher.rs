use crate::{
    cli::HashAlgorithms,
    hash_functions::Function,
    hasher::NoHashHasher,
    signature::Signature,
    sketch::{Sketch, Stats},
};
use needletail::{parser::SequenceRecord, Sequence};
use std::{
    collections::{BinaryHeap, HashMap},
    hash::BuildHasherDefault,
};

#[derive(Debug, Default)]
struct SketchHelper {
    pub kmer_budget: u64,
    pub max_hash: u64,
    counter: u64,
    seq_counter: u64,
    pub nmax: u64,
    pub global_heap: BinaryHeap<u64>,
    pub local_heap: Option<(u64, BinaryHeap<u64>)>,
    pub hashes: HashMap<u64, Option<Stats>, BuildHasherDefault<NoHashHasher>>,
    current_stat: Option<Stats>,
}

impl SketchHelper {
    pub fn new(kmer_budget: u64, max_hash: u64, nmin: Option<u64>, nmax: Option<u64>) -> Self {
        let local_heap = if let Some(nmin) = nmin {
            Some((nmin, BinaryHeap::with_capacity((nmin * 2) as usize)))
        } else {
            None
        };

        SketchHelper {
            kmer_budget,
            nmax: nmax.unwrap_or_else(|| u64::MAX),
            counter: 0,
            seq_counter: 0,
            max_hash: max_hash,
            global_heap: BinaryHeap::with_capacity(1_000_000 as usize),
            local_heap,
            hashes: HashMap::with_capacity_and_hasher(
                kmer_budget as usize,
                BuildHasherDefault::<NoHashHasher>::default(),
            ),
            current_stat: None,
        }
    }

    pub fn initialize_record(&mut self, stats: Option<Stats>) {
        self.current_stat = stats;
    }

    pub fn push(&mut self, hash: u64) {
        // Increase the local sequence counter
        self.seq_counter += 1;

        // Check if the hash is smaller than the global max hash cutoff
        if hash < self.max_hash {
            if self.global_heap.len() < self.nmax as usize
                || self.global_heap.len() < self.kmer_budget as usize
            {
                self.global_heap.push(hash);
                self.hashes.insert(hash, self.current_stat.clone());
            } else {
                let mut max = self.global_heap.peek_mut().unwrap();
                if hash < *max {
                    *max = hash;
                    self.hashes.insert(hash, self.current_stat.clone());
                    self.hashes.remove(&max);
                }
            }
        }

        // If there is a local_heap
        if let Some((nmin, local_heap)) = &mut self.local_heap {
            if local_heap.len() < *nmin as usize {
                local_heap.push(hash);
            } else {
                let mut max = local_heap.peek_mut().unwrap();
                if hash < *max {
                    *max = hash;
                }
            }
        }
    }

    pub fn next_record(&mut self) {
        self.counter += self.seq_counter;
        if let Some((_, local_heap)) = &mut self.local_heap {
            self.hashes
                .extend(local_heap.drain().map(|x| (x, self.current_stat.clone())));
        }
        self.seq_counter = 0;
    }

    pub fn reset(&mut self) {
        self.global_heap.clear();
        self.hashes.clear();
        self.counter = 0;
    }

    pub fn into_sketch(&mut self, name: String, kmer_size: u8) -> Sketch {
        let mut sketch = Sketch::new(name, self.hashes.len(), self.counter as usize, kmer_size);
        sketch.hashes = self.hashes.drain().collect();
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
    stats: bool,
    function: Function<'a>,
    algorithm: HashAlgorithms,
}

impl<'a> Sketcher<'a> {
    pub fn new(
        kmer_length: u8,
        name: String,
        singleton: bool,
        stats: bool,
        budget: u64,
        max_hash: u64,
        nmin: Option<u64>,
        nmax: Option<u64>,
        function: Function<'a>,
        algorithm: HashAlgorithms,
    ) -> Self {
        Sketcher {
            name,
            kmer_length,
            helper: SketchHelper::new(budget, max_hash, nmin, nmax),
            singleton,
            completed_sketches: Vec::new(),
            function,
            stats,
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
        let stats = if self.stats {
            Some(Stats::from_seq(seq.sequence()))
        } else {
            None
        };
        if self.kmer_length <= 31 {
            let func_small = self.function.get_small().unwrap();
            let seq_normalized = seq.normalize(true);
            self.helper.initialize_record(stats);

            for (_, kmer, _) in seq_normalized.bit_kmers(self.kmer_length, true) {
                self.helper.push(func_small(kmer.0));
            }

            self.helper.next_record();

            if self.singleton {
                self.completed_sketches.push(self.helper.into_sketch(
                    String::from_utf8_lossy(seq.id()).to_string(),
                    self.kmer_length,
                ));
            }
        } else {
            let func_large = self.function.get_large().unwrap();
            self.helper.initialize_record(stats);
            let rc = seq.reverse_complement();
            for (_, kmer, _) in seq.normalize(false).canonical_kmers(self.kmer_length, &rc) {
                self.helper.push(func_large(kmer));
            }
            self.helper.next_record();
            if self.singleton {
                self.completed_sketches.push(self.helper.into_sketch(
                    String::from_utf8_lossy(seq.id()).to_string(),
                    self.kmer_length,
                ));
            }
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
