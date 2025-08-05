use std::{
    cmp::max,
    sync::{
        Arc,
        atomic::{AtomicU64, Ordering},
    },
    u64,
};

//use arc_swap::{ArcSwap, ArcSwapAny};
use jam_rs::hash_functions::{double_fold, jamhash};
use rand::{Rng, rng};

// a429664bff768316, c6d35be8d0acb457

#[inline]
pub fn murmur3_u64(kmer: u64) -> u64 {
    fastmurmur3::murmur3_x64_128(&kmer.to_be_bytes(), 42) as u64
}

#[inline]
pub fn xxhash3(kmer: u64) -> u64 {
    xxhash_rust::xxh3::xxh3_64(&kmer.to_be_bytes())
}

fn compute_u64_avalanche() {
    let mut rng = rng();

    let mut bit_flips_jam = vec![0; 64 * 64];
    let mut bit_flips_murmur = vec![0; 64 * 64];
    let mut bit_flips_xxhash = vec![0; 64 * 64];

    let mut jam_values = vec![];
    let mut murmur_values = vec![];
    let mut xxhash_values = vec![];

    for iter in 0..1000 {
        for _ in 0..10000 {
            let base_val: u64 = rng.random();
            let jam_hash = jamhash(base_val);
            for flip_pos in 0..64 {
                let delta_val = base_val ^ (1 << flip_pos);
                let delta_hash = jamhash(delta_val);
                for test_pos in 0..64 {
                    let flipped = ((jam_hash ^ delta_hash) >> test_pos) & 1;
                    bit_flips_jam[test_pos * 64 + flip_pos] += flipped as usize;
                }
            }

            let murmur_hash = murmur3_u64(base_val);

            for flip_pos in 0..64 {
                let delta_val = base_val ^ (1 << flip_pos);
                let delta_hash = murmur3_u64(delta_val);
                for test_pos in 0..64 {
                    let flipped = ((murmur_hash ^ delta_hash) >> test_pos) & 1;
                    bit_flips_murmur[test_pos * 64 + flip_pos] += flipped as usize;
                }
            }

            let xxhash = xxhash3(base_val);
            for flip_pos in 0..64 {
                let delta_val = base_val ^ (1 << flip_pos);
                let delta_hash = xxhash3(delta_val);
                for test_pos in 0..64 {
                    let flipped = ((xxhash ^ delta_hash) >> test_pos) & 1;
                    bit_flips_xxhash[test_pos * 64 + flip_pos] += flipped as usize;
                }
            }
        }

        jam_values.push(analyze_bitflips("jamhash", bit_flips_jam.clone()));
        murmur_values.push(analyze_bitflips("murmur3", bit_flips_murmur.clone()));
        xxhash_values.push(analyze_bitflips("xxhash3", bit_flips_xxhash.clone()));

        bit_flips_jam.fill(0);
        bit_flips_murmur.fill(0);
        bit_flips_xxhash.fill(0);

        if iter % 10 == 0 {
            println!(
                "Iteration {iter}: Jam: {:.4}, Murmur: {:.4}, XXHash: {:.4}",
                jam_values.last().unwrap().0,
                murmur_values.last().unwrap().0,
                xxhash_values.last().unwrap().0
            );
        }
    }

    // Write results to CSV files
    std::fs::create_dir_all("out").unwrap();
    std::fs::write(
        "out/avalanche-jam.csv",
        jam_values
            .iter()
            .map(|v| format!("{:.4},{:.4}", v.0, v.1))
            .collect::<Vec<_>>()
            .join("\n"),
    )
    .unwrap();
    std::fs::write(
        "out/avalanche-murmur.csv",
        murmur_values
            .iter()
            .map(|v| format!("{:.4},{:.4}", v.0, v.1))
            .collect::<Vec<_>>()
            .join("\n"),
    )
    .unwrap();
    std::fs::write(
        "out/avalanche-xxhash.csv",
        xxhash_values
            .iter()
            .map(|v| format!("{:.4},{:.4}", v.0, v.1))
            .collect::<Vec<_>>()
            .join("\n"),
    )
    .unwrap();
}

pub fn analyze_bitflips(algorithm_name: &str, bit_flips: Vec<usize>) -> (f64, f64) {
    let avg_bit_flips: f64 = bit_flips.iter().sum::<usize>() as f64 / bit_flips.len() as f64;

    let variance = bit_flips
        .iter()
        .map(|x| {
            let diff = *x as f64 - avg_bit_flips as f64;
            diff * diff
        })
        .sum::<f64>()
        / bit_flips.len() as f64;

    let stddev = variance.sqrt();

    (avg_bit_flips, stddev)
}

fn main() {
    compute_u64_avalanche();
}
