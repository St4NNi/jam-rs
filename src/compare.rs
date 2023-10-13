use crate::signature::Signature;
use crate::sketch::Sketch;
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
    pub estimated_containment: f64,
}

impl Display for CompareResult {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        if self.reverse {
            write!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}",
                self.to_name,
                self.from_name,
                self.num_common,
                self.num_kmers,
                self.num_common as f64 / self.num_kmers as f64 * 100.0, // Percent
                self.estimated_containment
            )?;
            Ok(())
        } else {
            write!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}",
                self.from_name,
                self.to_name,
                self.num_common,
                self.num_kmers,
                self.num_common as f64 / self.num_kmers as f64 * 100.0,
                self.estimated_containment
            )
        }
    }
}

pub struct MultiComp {
    from: Vec<Sketch>,
    to: Vec<Sketch>,
    results: Vec<CompareResult>,
    threads: usize,
    kmer_size: u8,
    cutoff: f64,
}

impl MultiComp {
    pub fn new(
        mut from: Vec<Signature>,
        mut to: Vec<Signature>,
        threads: usize,
        cutoff: f64,
    ) -> Result<Self> {
        let kmer_size = from
            .first()
            .ok_or_else(|| anyhow!("Empty from list"))?
            .kmer_size;

        Ok(MultiComp {
            from: from.iter_mut().map(|e| e.collapse()).collect(),
            to: to.iter_mut().map(|e| e.collapse()).collect(),
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
                            "Kmer sizes do not match, expected: {}, got: {}",
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
    larger: &'a Sketch,
    smaller: &'a Sketch,
    num_kmers: usize,
    num_common: usize,
    reverse: bool,
}

impl<'a> Comparator<'a> {
    pub fn new(sketch_a: &'a Sketch, sketch_b: &'a Sketch) -> Self {
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
        for (hash, _) in &self.smaller.hashes {
            self.num_kmers += 1;
            if self.larger.hashes.contains_key(hash) {
                self.num_common += 1;
            }
        }
    }

    pub fn finalize(self) -> CompareResult {
        // Eg 0.1
        let larger_fraction = self.larger.num_kmers as f64 / self.larger.max_kmers as f64;
        // Eg 1.0
        let smaller_fraction = self.smaller.num_kmers as f64 / self.smaller.max_kmers as f64;
        // How much smaller is the smaller sketch
        let fraction = if larger_fraction < smaller_fraction {
            smaller_fraction / larger_fraction
        } else {
            larger_fraction / smaller_fraction
        };
        let estimated_containment =
            self.num_common as f64 / self.num_kmers as f64 * fraction * 100.0;

        CompareResult {
            from_name: self.larger.name.clone(),
            to_name: self.smaller.name.clone(),
            num_kmers: self.num_kmers,
            num_common: self.num_common,
            reverse: self.reverse,
            estimated_containment,
        }
    }

    #[allow(dead_code)]
    pub fn reset(&mut self) {
        self.num_kmers = 0;
        self.num_common = 0;
    }
}
