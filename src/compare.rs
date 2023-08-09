use crate::sketcher;
use anyhow::Result;
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};
use std::{
    fmt::{self, Display, Formatter},
    ops::DerefMut,
    sync::Mutex,
};

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
    threads: usize,
}

impl MultiComp {
    pub fn new(from: Vec<sketcher::Sketch>, to: Vec<sketcher::Sketch>, threads: usize) -> Self {
        MultiComp {
            from,
            to,
            results: Vec::new(),
            threads,
        }
    }

    pub fn compare(&mut self) -> Result<()> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;

        let results = Mutex::new(Vec::new());

        pool.install(|| {
            self.from.par_iter().try_for_each(|origin| {
                self.to.par_iter().try_for_each(|target| {
                    let mut comparator = Comparator::new(origin, target);
                    comparator.compare();
                    results
                        .lock()
                        .unwrap()
                        .deref_mut()
                        .push(comparator.finalize());
                    Ok::<(), anyhow::Error>(())
                })
            })
        })?;

        self.results = results.into_inner().unwrap();
        Ok(())
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
