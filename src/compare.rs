use crate::sketcher;
use anyhow::anyhow;
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
    pub reverse: bool,
}

impl Display for CompareResult {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        if self.reverse {
            write!(
                f,
                "{}\t{}\t{}\t{}\t{}",
                self.to_name,
                self.from_name,
                self.num_common,
                self.num_kmers,
                self.num_common as f64 / self.num_kmers as f64 * 100.0
            )?;
            return Ok(());
        } else {
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
}

pub struct MultiComp {
    from: Vec<sketcher::Sketch>,
    to: Vec<sketcher::Sketch>,
    results: Vec<CompareResult>,
    threads: usize,
    kmer_size: u8,
    cutoff: f64,
}

impl MultiComp {
    pub fn new(
        from: Vec<sketcher::Sketch>,
        to: Vec<sketcher::Sketch>,
        threads: usize,
        cutoff: f64,
    ) -> Result<Self> {
        let kmer_size = from
            .first()
            .ok_or_else(|| anyhow!("Empty from list"))?
            .kmer_size
            .clone();

        Ok(MultiComp {
            from,
            to,
            results: Vec::new(),
            threads,
            kmer_size,
            cutoff,
        })
    }

    pub fn compare(&mut self) -> Result<()> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;

        let results = Mutex::new(Vec::new());

        pool.install(|| {
            self.from.par_iter().try_for_each(|origin| {
                self.to.par_iter().try_for_each(|target| {
                    if target.kmer_size != self.kmer_size || origin.kmer_size != self.kmer_size {
                        return Err(anyhow!(
                            "Kmer sizes do not match expected: {} got: {}",
                            self.kmer_size,
                            origin.kmer_size
                        ));
                    }
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
            .into_iter()
            .filter(|e| e.num_common as f64 / e.num_kmers as f64 * 100.0 > self.cutoff)
            .collect()
    }
}

pub struct Comparator<'a> {
    larger: &'a sketcher::Sketch,
    smaller: &'a sketcher::Sketch,
    num_kmers: usize,
    num_common: usize,
    reverse: bool,
}

impl<'a> Comparator<'a> {
    pub fn new(sketch_a: &'a sketcher::Sketch, sketch_b: &'a sketcher::Sketch) -> Self {
        let (larger, smaller, reverse) = if sketch_a.hashes.len() > sketch_b.hashes.len() {
            (sketch_a, sketch_b, false)
        } else {
            (sketch_b, sketch_a, true)
        };
        Comparator {
            larger,
            smaller,
            num_kmers: 0,
            num_common: 0,
            reverse,
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
            reverse: self.reverse,
        }
    }

    #[allow(dead_code)]
    pub fn reset(&mut self) {
        self.num_kmers = 0;
        self.num_common = 0;
    }
}
