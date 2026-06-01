//! Small generated shards for the search tests.

use crate::index::{BuildOptions, EntryInput, Parameters, Shard, build_shard};
use crate::loss::MissModel;
use std::path::Path;
use std::sync::OnceLock;

/// The shard parameters of the tests: small minimizer windows and bins, so short queries have
/// enough seeds.
pub fn parameters() -> Parameters {
    Parameters {
        k: 21,
        window: 20,
        bin_bits: 8,
        ..Parameters::default()
    }
}

/// The miss model of `parameters`, built once.
pub fn model() -> &'static MissModel {
    static MODEL: OnceLock<MissModel> = OnceLock::new();
    MODEL.get_or_init(|| MissModel::with_layouts(21, 20, 200).unwrap())
}

/// A small xorshift generator.
pub struct Random(pub u64);

impl Random {
    pub fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    pub fn bases(&mut self, length: usize) -> Vec<u8> {
        (0..length)
            .map(|_| b"ACGT"[(self.next() % 4) as usize])
            .collect()
    }

    /// A copy of `sequence` with substitutions, each base changed with chance `1 - identity`.
    pub fn mutate(&mut self, sequence: &[u8], identity: f64) -> Vec<u8> {
        sequence
            .iter()
            .map(|&base| {
                let roll = (self.next() % 1_000_000) as f64 / 1_000_000.0;
                if roll < 1.0 - identity {
                    let other = b"ACGT".iter().filter(|&&b| b != base);
                    *other.clone().nth((self.next() % 3) as usize).unwrap()
                } else {
                    base
                }
            })
            .collect()
    }
}

/// Builds a shard in `directory/shard` with one entry per item, each a list of contigs.
pub fn build(directory: &Path, entries: &[Vec<Vec<u8>>], bucket_bits: Option<u8>) -> Shard {
    let inputs: Vec<EntryInput> = entries
        .iter()
        .enumerate()
        .map(|(index, contigs)| {
            let path = directory.join(format!("entry{index}.fa"));
            let mut text = Vec::new();
            for (contig, sequence) in contigs.iter().enumerate() {
                text.extend_from_slice(format!(">entry{index}_contig{contig}\n").as_bytes());
                text.extend_from_slice(sequence);
                text.push(b'\n');
            }
            std::fs::write(&path, text).unwrap();
            EntryInput {
                name: format!("entry{index}"),
                path,
            }
        })
        .collect();
    let output = directory.join("shard");
    let options = BuildOptions {
        parameters: parameters(),
        bucket_bits,
        threads: 2,
        temp_dir: directory.to_path_buf(),
        entry_per_record: false,
    };
    build_shard(&inputs, &output, &options).unwrap();
    Shard::open(&output).unwrap()
}
