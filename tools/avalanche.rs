use rand::{Rng, rng};

// Adapted from https://github.com/orlp/foldhash/blob/master/benches/avalanche.rs
fn compute_u64_avalanche<F: Fn(u64) -> u64>(
    num_hashers: usize,
    iters_per_hasher: usize,
    hash_fn: F,
) -> Vec<f64> {
    let mut rng = rng();
    let mut worst_bias = vec![0.5f64; 64 * 64];
    for _ in 0..num_hashers {
        let mut bit_flips = vec![0; 64 * 64];
        for _ in 0..iters_per_hasher {
            let base_val: u64 = rng.random();
            let base_hash = hash_fn(base_val);
            for flip_pos in 0..64 {
                let delta_val = base_val ^ (1 << flip_pos);
                let delta_hash = hash_fn(delta_val);

                for test_pos in 0..64 {
                    let flipped = ((base_hash ^ delta_hash) >> test_pos) & 1;
                    bit_flips[test_pos * 64 + flip_pos] += flipped as usize;
                }
            }
        }

        for i in 0..64 * 64 {
            let flip_frac = bit_flips[i] as f64 / iters_per_hasher as f64;
            if (flip_frac - 0.5).abs() > (worst_bias[i] - 0.5).abs() {
                worst_bias[i] = flip_frac;
            }
        }
    }

    worst_bias
}

fn write_avalanche_csv<F: Fn(u64) -> u64>(name: &str, hash_fn: F) {
    println!("calculating avalanche properties of {name}");
    let strings: Vec<String> = compute_u64_avalanche(10000, 1000, hash_fn)
        .into_iter()
        .map(|b| format!("{b}"))
        .collect();
    std::fs::create_dir_all("out").unwrap();
    std::fs::write(format!("out/avalanche-{name}.csv"), strings.join(",")).unwrap();
}

#[inline]
pub fn murmur3_u64(kmer: u64) -> u64 {
    fastmurmur3::murmur3_x64_128(&kmer.to_be_bytes(), 42) as u64
}

#[inline]
pub fn xxhash3(kmer: u64) -> u64 {
    xxhash_rust::xxh3::xxh3_64(&kmer.to_be_bytes())
}

fn main() {
    write_avalanche_csv("jamhash", jam_rs::hash_functions::jamhash);
    write_avalanche_csv("murmur3", murmur3_u64);
    write_avalanche_csv("xxhash3", xxhash3);
}
