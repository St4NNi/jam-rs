use crate::sketcher;

pub struct Comparator {
    larger: sketcher::Sketch,
    smaller: sketcher::Sketch,
    num_kmers: usize,
    num_common: usize,
}

impl Comparator {
    pub fn new(sketch_a: sketcher::Sketch, sketch_b: sketcher::Sketch) -> Self {
        let (larger, smaller) = if sketch_a.max_kmers > sketch_b.max_kmers {
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

    pub fn finalize(self) -> (usize, usize) {
        (self.num_kmers, self.num_common)
    }

    pub fn reset(&mut self) {
        self.num_kmers = 0;
        self.num_common = 0;
    }
}
