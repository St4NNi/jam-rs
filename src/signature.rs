use crate::{cli::HashAlgorithms, sketch::Sketch};
use serde::{Deserialize, Serialize};
use sourmash::signature::{Signature as SourmashSignature, SigsTrait};
use std::collections::BTreeSet;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Signature {
    pub file_name: String,
    pub sketches: Vec<Sketch>,
    pub algorithm: HashAlgorithms,
    pub kmer_size: u8,
    pub max_hash: u64,
}

impl From<Signature> for SourmashSignature {
    fn from(val: Signature) -> Self {
        SourmashSignature::builder()
            .hash_function(format!("{:?}", val.algorithm))
            .filename(Some(val.file_name))
            .email("".to_string())
            .license("CC0".to_string())
            .name(None)
            .signatures(
                val.sketches
                    .into_iter()
                    .map(|sketch| sketch.into_sourmash(val.max_hash))
                    .collect(),
            )
            .build()
    }
}

impl From<SourmashSignature> for Signature {
    fn from(sourmash_signature: SourmashSignature) -> Self {
        let mut sketches = Vec::new();
        let mut max_hash = None;
        let mut kmer_size = None;
        for sketch in sourmash_signature.sketches() {
            match sketch {
                sourmash::sketch::Sketch::MinHash(mash) => {
                    if let Some(max_hash) = max_hash {
                        if max_hash != mash.max_hash() {
                            panic!("Max hash of sketches is not equal");
                        }
                    } else {
                        max_hash = Some(mash.max_hash());
                    }

                    if let Some(kmer_size) = kmer_size {
                        if kmer_size != mash.ksize() as u8 {
                            panic!("Kmer size of sketches is not equal");
                        }
                    } else {
                        kmer_size = Some(mash.ksize() as u8);
                    }

                    let mut sketch = Sketch::new(
                        sourmash_signature.filename(),
                        mash.mins().len(),
                        mash.ksize() as u8,
                    );
                    sketch.hashes = mash.mins().into_iter().collect::<BTreeSet<u64>>();
                    sketches.push(sketch);
                }
                sourmash::sketch::Sketch::LargeMinHash(mash) => {
                    if let Some(max_hash) = max_hash {
                        if max_hash != mash.max_hash() {
                            panic!("Max hash of sketches is not equal");
                        }
                    } else {
                        max_hash = Some(mash.max_hash());
                    }

                    if let Some(kmer_size) = kmer_size {
                        if kmer_size != mash.ksize() as u8 {
                            panic!("Kmer size of sketches is not equal");
                        }
                    } else {
                        kmer_size = Some(mash.ksize() as u8);
                    }

                    let mut sketch = Sketch::new(
                        sourmash_signature.filename(),
                        mash.mins().len(),
                        mash.ksize() as u8,
                    );
                    sketch.hashes = mash.mins().into_iter().collect::<BTreeSet<u64>>();
                    sketches.push(sketch);
                }
                sourmash::sketch::Sketch::HyperLogLog(_) => {
                    unimplemented!("HyperLogLog sketches are not supported")
                }
            }
        }
        Signature {
            file_name: sourmash_signature.filename(),
            sketches,
            algorithm: HashAlgorithms::Murmur3,
            kmer_size: kmer_size.expect("No sketch with kmer_size found"),
            max_hash: max_hash.expect("No sketch with max hash found"),
        }
    }
}

impl Signature {
    pub fn collapse(&mut self) -> Sketch {
        let mut sketch = Sketch::new(self.file_name.to_string(), 0, self.kmer_size);
        for old_sketch in self.sketches.drain(..) {
            sketch.hashes.extend(old_sketch.hashes);
            sketch.num_kmers += old_sketch.num_kmers;
        }
        sketch
    }
}
