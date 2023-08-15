use crate::hasher::NoHashHasher;
use needletail::Sequence;
use serde::{Deserialize, Serialize};
use std::{
    collections::{BinaryHeap, HashMap},
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

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Sketch {
    pub name: String, // Name of file or sequence
    pub hashes: HashMap<u64, Stats, BuildHasherDefault<NoHashHasher>>, // Hashes with stats
    pub num_kmers: usize, // Number of kmers (collected)
    pub max_kmers: usize, // Max number of kmers (budget)
    pub kmer_size: u8, // Kmer size
}

#[derive(Debug, Default)]
struct SketchHelper {
    pub kmer_budget: u64,
    pub max_hash: u64,
    pub binary_heap: BinaryHeap<u64>,
    pub hashes: HashMap<u64, Stats, BuildHasherDefault<NoHashHasher>>,
}

impl SketchHelper {
    pub fn new(kmer_budget: u64, max_hash: Option<u64>) -> Self {
        SketchHelper {
            kmer_budget,
            max_hash: max_hash.unwrap_or(u64::MAX),
            binary_heap: BinaryHeap::with_capacity(kmer_budget as usize),
            hashes: HashMap::with_capacity_and_hasher(
                kmer_budget as usize,
                BuildHasherDefault::<NoHashHasher>::default(),
            ),
        }
    }

    pub fn push(&mut self, hash: u64, stats: Option<Stats>) {
        if hash < self.max_hash {
            if self.binary_heap.len() < self.kmer_budget as usize {
                self.binary_heap.push(hash);
                if let Some(stats) = stats {
                    self.hashes.insert(hash, stats);
                }
            } else {
                let mut max = self.binary_heap.peek_mut().unwrap();
                if hash < *max {
                    *max = hash;
                    if let Some(stats) = stats {
                        self.hashes.insert(hash, stats);
                    }
                    self.hashes.remove(&max);
                }
            }
        }
    }
}

pub struct Sketcher<'a> {
    kmer_length: u8,
    helper: SketchHelper,
    current_sketch: Sketch,
    completed_sketches: Vec<Sketch>,
    singleton: bool,
    stats: bool,
    function: Function<'a>,
}

pub enum Function<'a> {
    Large(&'a dyn Fn(&[u8]) -> u64),
    Small(&'a dyn Fn(u64) -> u64),
}

impl Function<'_> {
    pub fn get_large(&self) -> Option<&dyn Fn(&[u8]) -> u64> {
        match self {
            Function::Large(f) => Some(f),
            _ => None,
        }
    }
    pub fn get_small(&self) -> Option<&dyn Fn(u64) -> u64> {
        match self {
            Function::Small(f) => Some(f),
            _ => None,
        }
    }
}

impl<'a> Sketcher<'_> {
    pub fn new(
        kmer_length: u8,
        name: String,
        singleton: bool,
        stats: bool,
        budget: Option<u64>,
        max_hash: Option<u64>,
        nmin: Option<u64>,
        nmax: Option<u64>,
        function: Function,
    ) -> Self {
        let budget = budget.unwrap_or(u64::MAX);

        Sketcher {
            kmer_length,
            helper: SketchHelper::new(budget, max_hash),
            current_sketch: Sketch {
                hashes: HashMap::with_capacity_and_hasher(
                    budget as usize,
                    BuildHasherDefault::<NoHashHasher>::default(),
                ),
                name,
                kmer_size: kmer_length,
                num_kmers: 0,
                max_kmers: 0,
            },
            singleton,
            completed_sketches: Vec::with_capacity(budget as usize),
            function,
            stats,
        }
    }
}

impl Sketcher<'_> {
    // This is more or less derived from the `process` method in `finch-rs`:
    // https://github.com/onecodex/finch-rs/blob/master/lib/src/sketch_schemes/mash.rs
    pub fn process_small<'seq, 'a, 'inner>(&'a mut self, seq: &'seq dyn Sequence<'inner>)
    where
        'a: 'seq,
        'seq: 'inner,
    {
        let stats = if self.stats {
            Some(Stats::from_seq(seq.sequence()))
        } else {
            None
        };
        let func_small = self.function.get_small().unwrap();
        let seq = seq.normalize(true);

        for (_, kmer, _) in seq.bit_kmers(self.kmer_length, true) {
            self.helper.push(func_small(kmer.0), stats.clone());
        }
    }

    // This is more or less derived from the `process` method in `finch-rs`:
    // https://github.com/onecodex/finch-rs/blob/master/lib/src/sketch_schemes/mash.rs
    /// To process larger kmers > 31 bases
    pub fn process_large<'seq, 'a, 'inner>(&'a mut self, seq: &'seq dyn Sequence<'inner>)
    where
        'a: 'seq,
        'seq: 'inner,
    {
        let func_large = self.function.get_large().unwrap();
        let stats = if self.stats {
            Some(Stats::from_seq(seq.sequence()))
        } else {
            None
        };
        let rc = seq.reverse_complement();
        for (_, kmer, is_rev_complement) in
            seq.normalize(false).canonical_kmers(self.kmer_length, &rc)
        {
            self.helper.push(func_large(kmer), stats.clone());
        }
    }

    pub fn finalize(self) -> Sketch {
        todo!()
    }
}
