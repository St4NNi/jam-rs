use crate::{cli::HashAlgorithms, sketch::Sketch};
use serde::{Deserialize, Serialize};
use sourmash::signature::Signature as SourmashSignature;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Signature {
    pub file_name: String,
    pub sketches: Vec<Sketch>,
    pub algorithm: HashAlgorithms,
    pub kmer_size: u8,
    pub max_hash: u64,
}

impl Into<SourmashSignature> for Signature {
    fn into(self) -> SourmashSignature {
        SourmashSignature::builder()
            .hash_function(format!("{:?}", self.algorithm))
            .filename(Some(self.file_name))
            .email("".to_string())
            .license("CC0".to_string())
            .name(None)
            .signatures(
                self.sketches
                    .into_iter()
                    .map(|sketch| sketch.into_sourmash(self.max_hash))
                    .collect(),
            )
            .build()
    }
}

impl Signature {
    pub fn collapse(&mut self) -> Sketch {
        let mut sketch = Sketch::new(self.file_name.to_string(), 0, 0, self.kmer_size);
        for old_sketch in self.sketches.drain(..) {
            sketch.hashes.extend(old_sketch.hashes);
            sketch.num_kmers += old_sketch.num_kmers;
            sketch.max_kmers += old_sketch.max_kmers;
        }
        sketch
    }
}
