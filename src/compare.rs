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
    pub option_num_skipped: Option<usize>,
    pub reverse: bool,
    pub estimated_containment: f64,
}

impl Display for CompareResult {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        if self.reverse {
            write!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.to_name,
                self.from_name,
                self.num_common,
                self.num_kmers,
                self.num_common as f64 / self.num_kmers as f64 * 100.0, // Percent
                self.estimated_containment,
                self.option_num_skipped.unwrap_or(0)
            )?;
            Ok(())
        } else {
            write!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.from_name,
                self.to_name,
                self.num_common,
                self.num_kmers,
                self.num_common as f64 / self.num_kmers as f64 * 100.0,
                self.estimated_containment,
                self.option_num_skipped.unwrap_or(0)
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
    use_stats: bool,
    gc_bounds: Option<(u8, u8)>,
}

impl MultiComp {
    pub fn new(
        mut from: Vec<Signature>,
        mut to: Vec<Signature>,
        threads: usize,
        cutoff: f64,
        use_stats: bool,
        gc_bounds: Option<(u8, u8)> 
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
            use_stats,
            gc_bounds,
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
                    let mut comparator = Comparator::new(origin, target, self.use_stats, self.gc_bounds);
                    comparator.compare()?;
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
    num_skipped: usize,
    reverse: bool,
    use_stats: bool,
    gc_bounds: Option<(u8, u8)>,
}

impl<'a> Comparator<'a> {
    pub fn new(sketch_a: &'a Sketch, sketch_b: &'a Sketch, use_stats: bool, gc_bounds: Option<(u8, u8)>) -> Self {
        let (larger, smaller, reverse) = if sketch_a.hashes.len() > sketch_b.hashes.len() {
            // DATABASE, INPUT -> Reverse = false
            (sketch_a, sketch_b, false)
        } else {
            // INPUT, DATABASE -> Reverse = true
            (sketch_b, sketch_a, true)
        };
        Comparator {
            larger,
            smaller,
            num_kmers: 0,
            num_common: 0,
            num_skipped: 0,
            reverse,
            use_stats,
            gc_bounds,
        }
    }

    // Stats handling:
    // GC & Size for the original contig are stored in the Stats struct
    // This comparison is always in relation to the query sketch
    // If reverse is true, the query sketch is the larger sketch
    #[inline]
    pub fn compare(&mut self) -> Result<()> {
        if self.use_stats {
            for (hash, stats) in &self.smaller.hashes {
                let smaller_stats = stats.as_ref().ok_or_else(|| anyhow!("Missing stats"))?;
                self.num_kmers += 1;
                if let Some(stats) = self.larger.hashes.get(hash) {
                    let larger_stats = stats.as_ref().ok_or_else(|| anyhow!("Missing stats"))?;
                    if self.reverse {
                        if !larger_stats.compare(smaller_stats, self.gc_bounds) {
                            self.num_skipped += 1;
                        }
                    } else {
                        if !smaller_stats.compare(larger_stats, self.gc_bounds) {
                            self.num_skipped += 1;
                        }
                    }
                    self.num_common += 1;
                };
            }
        } else {
            for (hash, _) in &self.smaller.hashes {
                self.num_kmers += 1;
                if self.larger.hashes.contains_key(hash) {
                    self.num_common += 1;
                };
            }
        }
        Ok(())
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
            option_num_skipped: if self.use_stats {
                Some(self.num_skipped)
            } else {
                None
            },
            reverse: self.reverse,
            estimated_containment,
        }
    }

    #[allow(dead_code)]
    pub fn reset(&mut self) {
        self.num_kmers = 0;
        self.num_common = 0;
        self.num_skipped = 0;
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    #[test]
    fn test_comp_without_stats() {
        let mut hashmap = HashMap::default();
        hashmap.extend([(1, None), (2, None), (3, None)]);
        let sketch_a = crate::sketch::Sketch {
            name: "a".to_string(),
            hashes: hashmap,
            num_kmers: 3,
            max_kmers: 10,
            kmer_size: 21,
        };
        let mut hashmap2 = HashMap::default();
        hashmap2.extend([(1, None), (2, None), (4, None)]);
        let sketch_b = crate::sketch::Sketch {
            name: "b".to_string(),
            hashes: hashmap2,
            num_kmers: 3,
            max_kmers: 10,
            kmer_size: 21,
        };

        let mut comp = super::Comparator::new(&sketch_a, &sketch_b, false, None);
        comp.compare().unwrap();
        let result = comp.finalize();
        assert_eq!(result.num_kmers, 3);
        assert_eq!(result.num_common, 2);
        assert_eq!(result.estimated_containment, 66.66666666666666);
        assert_eq!(result.option_num_skipped, None);
    }
}
