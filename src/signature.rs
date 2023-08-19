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
