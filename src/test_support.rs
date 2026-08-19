//! Test data shared by all tests: seeded random bases and bytes, mutated copies, test shards.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::index::{BuildOptions, EntryInput, Parameters, Shard, build_shard};
use crate::loss::MissModel;
use crate::search::SearchControl;
use crate::seeds::SeedScheme;
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};
use std::path::Path;
use std::sync::OnceLock;
use std::sync::atomic::AtomicBool;
use tokio::runtime::{Handle, Runtime};

/// A random generator that gives the same values for the same seed on every run.
pub fn rng(seed: u64) -> StdRng {
    StdRng::seed_from_u64(seed)
}

/// `length` random bases, each of `A`, `C`, `G` and `T` with the same chance.
pub fn random_bases(rng: &mut StdRng, length: usize) -> Vec<u8> {
    (0..length)
        .map(|_| b"ACGT"[rng.random_range(0..4)])
        .collect()
}

/// `length` random bytes.
pub fn random_bytes(rng: &mut StdRng, length: usize) -> Vec<u8> {
    let mut bytes = vec![0; length];
    rng.fill(&mut bytes[..]);
    bytes
}

/// The handle of one tokio runtime shared by all tests.
pub fn runtime() -> Handle {
    static RUNTIME: OnceLock<Runtime> = OnceLock::new();
    RUNTIME
        .get_or_init(|| Runtime::new().unwrap())
        .handle()
        .clone()
}

/// A copy of `bases` in which each base becomes another base with chance `1 - identity`.
pub fn with_substitutions(rng: &mut StdRng, bases: &[u8], identity: f64) -> Vec<u8> {
    bases
        .iter()
        .map(|&base| {
            if rng.random::<f64>() < identity {
                return base;
            }
            let others: Vec<u8> = b"ACGT"
                .iter()
                .copied()
                .filter(|&other| other != base)
                .collect();
            others[rng.random_range(0..others.len())]
        })
        .collect()
}

/// The FASTA text of the contigs of entry `entry`, named `entry{entry}_contig{index}`.
pub fn fasta_text(entry: usize, contigs: &[Vec<u8>]) -> Vec<u8> {
    let mut text = Vec::new();
    for (index, bases) in contigs.iter().enumerate() {
        text.extend_from_slice(format!(">entry{entry}_contig{index} description\n").as_bytes());
        for line in bases.chunks(70) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    text
}

/// Writes each entry (a list of contigs) to `directory/entry{index}.fa` and returns the inputs,
/// named `entry{index}`.
pub fn write_entries(directory: &Path, entries: &[Vec<Vec<u8>>]) -> Vec<EntryInput> {
    let inputs = entries.iter().enumerate().map(|(index, contigs)| {
        let path = directory.join(format!("entry{index}.fa"));
        std::fs::write(&path, fasta_text(index, contigs)).unwrap();
        EntryInput::from_path(path)
    });
    inputs.collect()
}

/// The parameters of the test shards: k 21, w 20 and bins of 256 bases, so short queries have
/// enough seeds.
pub fn parameters() -> Parameters {
    Parameters {
        scheme: SeedScheme {
            k: 21,
            w: 20,
            ..SeedScheme::default()
        },
        bin_bits: 8,
        ..Parameters::default()
    }
}

/// The build options of the test shards, with `parameters`.
pub fn build_options(directory: &Path, threads: usize) -> BuildOptions {
    BuildOptions {
        parameters: parameters(),
        bucket_bits: None,
        threads,
        temporary_directory: directory.to_path_buf(),
        entry_per_record: false,
    }
}

/// Builds and opens the shard `directory/shard` of `entries`, each a list of contigs.
pub fn test_shard(directory: &Path, entries: &[Vec<Vec<u8>>]) -> Shard {
    let inputs = write_entries(directory, entries);
    let output = directory.join("shard");
    build_shard(&inputs, &output, &build_options(directory, 2)).unwrap();
    Shard::open(&output).unwrap()
}

/// The miss model of the test shards, built once per test process. 50 layouts per length keep
/// the build short and are enough for tests.
pub fn miss_model() -> &'static MissModel {
    static MODEL: OnceLock<MissModel> = OnceLock::new();
    MODEL.get_or_init(|| MissModel::with_layouts(parameters().scheme, 50).unwrap())
}

/// A search control without a deadline whose stop flag is never set.
pub fn no_stop() -> SearchControl<'static> {
    static NEVER: AtomicBool = AtomicBool::new(false);
    SearchControl {
        deadline: None,
        stop: &NEVER,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seeds_repeat_bases() {
        assert_eq!(random_bases(&mut rng(7), 50), random_bases(&mut rng(7), 50));
        assert_ne!(random_bases(&mut rng(7), 50), random_bases(&mut rng(8), 50));
    }

    #[test]
    fn substitutions_follow_identity() {
        let mut rng = rng(1);
        let bases = random_bases(&mut rng, 100_000);
        let copy = with_substitutions(&mut rng, &bases, 0.9);
        let same = bases.iter().zip(&copy).filter(|(a, b)| a == b).count();
        assert!((89_000..91_000).contains(&same), "{same}");
        assert_eq!(with_substitutions(&mut rng, &bases, 1.0), bases);
    }
}
