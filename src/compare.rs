use crate::sketcher;
use serde::{Deserialize, Serialize};
use std::fmt::{self, Display, Formatter};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct CompareResult {
    pub from_name: String,
    pub to_name: String,
    pub num_common: usize,
    pub num_kmers: usize,
}

impl Display for CompareResult {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.from_name,
            self.to_name,
            self.num_common,
            self.num_kmers,
            self.num_common as f64 / self.num_kmers as f64 * 100.0
        )
    }
}

pub struct MultiComp {
    from: Vec<sketcher::Sketch>,
    to: Vec<sketcher::Sketch>,
    results: Vec<CompareResult>,
}

impl MultiComp {
    pub fn new(from: Vec<sketcher::Sketch>, to: Vec<sketcher::Sketch>) -> Self {
        MultiComp {
            from,
            to,
            results: Vec::new(),
        }
    }

    pub fn compare(&mut self) {
        for origin in &self.from {
            for target in &self.to {
                let mut comparator = Comparator::new(origin, target);
                comparator.compare();
                self.results.push(comparator.finalize());
            }
        }
    }

    pub fn finalize(self) -> Vec<CompareResult> {
        self.results
    }
}

pub struct Comparator<'a> {
    larger: &'a sketcher::Sketch,
    smaller: &'a sketcher::Sketch,
    num_kmers: usize,
    num_common: usize,
}

impl<'a> Comparator<'a> {
    pub fn new(sketch_a: &'a sketcher::Sketch, sketch_b: &'a sketcher::Sketch) -> Self {
        let (larger, smaller) = if sketch_a.hashes.len() > sketch_b.hashes.len() {
            (sketch_a, sketch_b)
        } else {
            (sketch_b, sketch_a)
        };
        Comparator {
            larger,
            smaller,
            num_kmers: 0,
            num_common: 0,
        }
    }

    #[inline]
    pub fn compare(&mut self) {
        for hash in &self.smaller.hashes {
            self.num_kmers += 1;
            if *hash < self.larger.lowest_hash {
                continue;
            }
            if self.larger.hashes.contains(hash) {
                self.num_common += 1;
            }
        }
    }

    pub fn finalize(self) -> CompareResult {
        CompareResult {
            from_name: self.larger.name.clone(),
            to_name: self.smaller.name.clone(),
            num_kmers: self.num_kmers,
            num_common: self.num_common,
        }
    }

    #[allow(dead_code)]
    pub fn reset(&mut self) {
        self.num_kmers = 0;
        self.num_common = 0;
    }
}
