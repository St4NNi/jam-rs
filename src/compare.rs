use crate::file_io::ShortSketchInfo;
use crate::signature::Signature;
use crate::sketch::Sketch;
use anyhow::anyhow;
use anyhow::Result;
use byteorder::BigEndian;
use heed::types::SerdeBincode;
use heed::types::U32;
use heed::types::U64;
use heed::DatabaseFlags;
use heed::EnvFlags;
use indicatif::ParallelProgressIterator;
use indicatif::ProgressBar;
use indicatif::ProgressDrawTarget;
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};
use std::cmp::max;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::RwLock;
use std::{
    fmt::{self, Display, Formatter},
    ops::DerefMut,
    sync::Mutex,
};

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
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
                "{}\t{}\t{}\t{}\t{:.2}",
                self.to_name,
                self.from_name,
                self.num_common,
                self.num_kmers,
                self.estimated_containment,
            )?;
            Ok(())
        } else {
            write!(
                f,
                "{}\t{}\t{}\t{}\t{:.2}",
                self.from_name,
                self.to_name,
                self.num_common,
                self.num_kmers,
                self.estimated_containment,
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
}

impl<'a> Comparator<'a> {
    pub fn new(sketch_a: &'a Sketch, sketch_b: &'a Sketch) -> Self {
        let (larger, smaller, reverse) = if sketch_a.hashes.len() >= sketch_b.hashes.len() {
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
        }
    }

    // Stats handling:
    // GC & Size for the original contig are stored in the Stats struct
    // This comparison is always in relation to the query sketch
    // If reverse is true, the query sketch is the larger sketch
    #[inline]
    pub fn compare(&mut self) -> Result<()> {
        self.num_kmers = max(self.larger.num_kmers, self.smaller.num_kmers);

        let mut larger = self.larger.hashes.iter();
        let mut smaller = self.smaller.hashes.iter();

        let mut larger_item = larger.next();
        let mut smaller_item = smaller.next();

        loop {
            match (larger_item, smaller_item) {
                (Some(l), Some(s)) => {
                    if l == s {
                        self.num_common += 1;
                        smaller_item = smaller.next();
                        larger_item = larger.next();
                    } else if l < s {
                        smaller_item = smaller.next();
                    } else {
                        larger_item = larger.next();
                    }
                }
                (Some(_), None) => {
                    larger_item = larger.next();
                }
                (None, Some(_)) => {
                    smaller_item = smaller.next();
                }
                (None, None) => break,
            }
        }

        Ok(())
    }

    pub fn finalize(self) -> CompareResult {
        // Eg 0.1
        let larger_fraction = self.larger.num_kmers as f64 / self.larger.hashes.len() as f64;
        // Eg 1.0
        let smaller_fraction = self.smaller.num_kmers as f64 / self.smaller.hashes.len() as f64;
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
        self.num_skipped = 0;
    }
}

pub struct LmdbComparator {
    pub signatures: Vec<Signature>,
    pub lmdb_env: heed::Env,
    pub threads: usize,
    pub cutoff: f64,
    pub infos: Arc<RwLock<HashMap<u32, ShortSketchInfo>>>,
    pub kmer_size: u8,
    pub fscale: Option<u64>,
    pub silent: bool,
}

impl LmdbComparator {
    pub fn new(lmdb_env: PathBuf, threads: usize, cutoff: f64, silent: bool) -> Result<Self> {
        let lmdb_env = unsafe {
            heed::EnvOpenOptions::new()
                .flags(EnvFlags::READ_ONLY | EnvFlags::NO_LOCK | EnvFlags::NO_SUB_DIR)
                .map_size(10 * 1024 * 1024 * 1024)
                .max_dbs(2)
                .open(lmdb_env)
                .unwrap()
        };

        let txn = lmdb_env.read_txn()?;

        let sigs_db = lmdb_env
            .open_database::<U32<BigEndian>, SerdeBincode<ShortSketchInfo>>(&txn, Some("sigs"))?
            .ok_or_else(|| anyhow!("Database sigs not found"))?;

        let infos = RwLock::new(HashMap::new());

        let mut kmer_size = None;
        let mut fscale = None;
        for sig in sigs_db.iter(&txn)? {
            let (key, value) = sig?;
            if let Some(kmer_size) = kmer_size {
                if kmer_size != value.kmer_size {
                    return Err(anyhow!("Kmer sizes do not match"));
                }
            } else {
                kmer_size = Some(value.kmer_size);
            }

            if fscale.is_some() {
                if fscale != value.fscale {
                    return Err(anyhow!("Fscale sizes do not match"));
                }
            } else {
                fscale = value.fscale;
            }

            infos.write().expect("poisoned lock").insert(key, value);
        }

        txn.commit()?;

        Ok(LmdbComparator {
            signatures: vec![],
            lmdb_env,
            threads,
            cutoff,
            infos: Arc::new(infos),
            kmer_size: kmer_size.unwrap(),
            fscale,
            silent,
        })
    }

    pub fn set_signatures(&mut self, signatures: Vec<Signature>) {
        self.signatures = signatures;
    }

    pub fn compare(&self) -> Result<Vec<CompareResult>> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;

        let results = Mutex::new(Vec::new());

            let pb = ProgressBar::new(self.signatures.len() as u64);
            pb.set_style(
                indicatif::ProgressStyle::default_bar()
                    .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")?
                    .progress_chars("##-"),
            );
        if self.silent {
            pb.set_draw_target(ProgressDrawTarget::hidden())
        }
        let infos = self.infos.clone();

        pool.install(|| {
            self.signatures
                .par_iter()
                .progress_with(pb)
                .try_for_each(|origin| {
                    origin.sketches.par_iter().try_for_each(|target| {
                        let txn = self.lmdb_env.read_txn()?;

                        let hashes = self
                            .lmdb_env
                            .database_options()
                            .types::<U64<BigEndian>, U32<BigEndian>>()
                            .name("hashes")
                            .flags(DatabaseFlags::DUP_SORT)
                            .open(&txn)?
                            .ok_or_else(|| anyhow!("Database hashes not found"))?;
                        let mut result_map = HashMap::new();

                        for hash in target.hashes.iter() {
                            if let Some(key) = hashes.get_duplicates(&txn, hash)? {
                                for item in key {
                                    let (_, sketch) = item?;
                                    let entry = result_map.entry(sketch).or_insert(0);
                                    *entry += 1u64;
                                }
                            };
                        }

                        let mut final_results = vec![];
                        for (idx, num_common) in result_map {
                            let read_infos = infos.read().expect("poisoned lock");
                            let infos = read_infos.get(&idx).expect("Key not found");
                            let num_kmers = if target.hashes.len() < infos.num_hashes {
                                target.hashes.len()
                            } else {
                                infos.num_hashes
                            };
                            let estimated_containment =
                                num_common as f64 / num_kmers as f64 * 100.0;
                            final_results.push(CompareResult {
                                from_name: target.name.clone(),
                                to_name: infos.file_name.clone(),
                                num_kmers,
                                num_common: num_common as usize,
                                reverse: false,
                                estimated_containment,
                            })
                        }

                        results
                            .lock()
                            .unwrap()
                            .extend(final_results.into_iter().filter(|e| {
                                e.num_common as f64 / e.num_kmers as f64 * 100.0 > self.cutoff
                            }));

                        Ok::<(), anyhow::Error>(())
                    })
                })
        })?;
        Ok(results.into_inner().expect("poisoned lock"))
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use crate::compare::CompareResult;

    #[test]
    fn test_comp_without_stats() {
        let mut bheap1 = BTreeSet::default();
        bheap1.extend([1, 2, 3]);
        let sketch_a = crate::sketch::Sketch {
            name: "a".to_string(),
            hashes: bheap1,
            num_kmers: 3,
            kmer_size: 21,
        };
        let mut bheap2 = BTreeSet::default();
        bheap2.extend([1, 2, 4]);
        let sketch_b = crate::sketch::Sketch {
            name: "b".to_string(),
            hashes: bheap2,
            num_kmers: 3,
            kmer_size: 21,
        };

        let mut comp = super::Comparator::new(&sketch_a, &sketch_b);
        comp.compare().unwrap();
        let result = comp.finalize();
        assert_eq!(result.num_kmers, 3);
        assert_eq!(result.num_common, 2);
        assert_eq!(result.estimated_containment, 66.66666666666666);

        let constructed_result = CompareResult {
            from_name: "a".to_string(),
            to_name: "b".to_string(),
            num_kmers: 3,
            num_common: 2,
            reverse: false,
            estimated_containment: 66.66666666666666,
        };
        assert_eq!(result, constructed_result);
    }
}
