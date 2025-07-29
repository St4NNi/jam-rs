use std::{sync::Arc, u64};

use arc_swap::{ArcSwap, ArcSwapAny};
use rand::{Rng, rng};

/// Cumulative Distribution Function for the Uniform Distribution.
fn cdf_uniform(x: u64) -> f64 {
    // Wish we had f128s. Gonna be issues here.
    (x as f64) / (std::u64::MAX as f64)
}

/// Compute the Kolmogorov-Smirnov test.
///
/// ECDF: Experimental Cumulative Distribution Function. The distribution represented by the
/// samples.
///
/// TCDF: Theoretical Cumulative Distribution Function. The theoretical distribution to be tested
/// against; in this case the uniform distribution.
fn ks(samples: &[u64]) -> f64 {
    let n = samples.len() as f64;
    let mut last_ecdf = 0.0f64;
    let mut ks = std::f64::MIN;
    for (i, x) in samples.iter().enumerate() {
        let tcdf = (i as f64) / n;
        let next_ecdf = cdf_uniform(*x);
        let d1 = (last_ecdf - tcdf).abs();
        let d2 = (tcdf - next_ecdf).abs();
        ks = ks.max(d1.max(d2));
        last_ecdf = next_ecdf;
    }
    ks
}

fn main() {
    let samples = (0..1_000_000).collect::<Vec<_>>();

    let mut results = vec![0u64; samples.len()];

    let const1 = 0x01c3defadcba8000u64;
    let const2 = 0xb85b352f48103b89u64;

    let current_best_ks = {
        for x in 0..samples.len() {
            results[x] = jam_rs::hash_functions::double_fold(samples[x], const1, const2);
        }

        results.sort();

        ks(&results).to_bits()
    };

    let current_best_deformity = test_max_deformity(const1, const2);

    let values = Arc::new(ArcSwap::from_pointee((
        current_best_ks,
        current_best_deformity,
    )));

    let num_threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4);

    // 04feab043089c9a1, 4c6102befcb61250, 430a4fbd481e7756, f640d38a7b661335

    let mut threads = vec![];
    for _ in 0..num_threads {
        let samples = samples.clone();
        let results = results.clone();
        let values = values.clone();

        threads.push(std::thread::spawn(move || {
            find_best_constants(samples, results, values);
        }));
    }

    for thread in threads {
        thread.join().expect("Thread panicked");
    }
}

fn find_best_constants(
    samples: Vec<u64>,
    mut results: Vec<u64>,
    values: Arc<ArcSwapAny<Arc<(u64, u64)>>>,
) {
    let mut rng = rng();
    loop {
        let const1 = rng.random::<u64>();
        let const2 = rng.random::<u64>();

        for x in 0..samples.len() {
            results[x] = jam_rs::hash_functions::double_fold(samples[x], const1, const2);
        }

        results.sort();

        let result = ks(&results);
        let best = values.load();

        let best_ks = best.0.clone();
        let best_deformity = best.1.clone();
        drop(best);
        if best_ks < result.to_bits() {
            continue; // Skip if the current result is worse than the best found so far.
        } else {
            let calc_deformity = test_max_deformity(const1, const2);

            if calc_deformity >= best_deformity {
                continue; // Skip if the deformity is worse than the best found so far.
            }

            let old_arc = values.rcu(|best| {
                if best.0 > result.to_bits() && best.1 > calc_deformity {
                    Arc::new((result.to_bits(), calc_deformity))
                } else {
                    best.clone()
                }
            });

            if old_arc.0 < best_ks && old_arc.1 < best_deformity {
                continue;
            }

            println!(
                "Found new KS / deformity: {result} / {calc_deformity} for constants {const1:016x}, {const2:016x}"
            );
        }
    }
}

fn unrolled_64bits(num: u64, nums: &mut [u64; 64]) {
    for i in 0..64 {
        if num & (1u64 << i) != 0 {
            nums[i] += 1;
        }
    }
}

const RANGE: std::ops::Range<u64> = 0..100_000_000u64;
const RANGE_LEN: u64 = RANGE.end - RANGE.start;

fn test_max_deformity(const1: u64, const2: u64) -> u64 {
    let mut jamhash_bits = [0u64; 64];

    for x in RANGE {
        let ah = jam_rs::hash_functions::double_fold(x, const1, const2);
        unrolled_64bits(ah, &mut jamhash_bits);
    }

    let mut max_deviation_jamhash = 0u64;
    for x in 0..64 {
        max_deviation_jamhash = ((jamhash_bits[x] as i128 - (RANGE_LEN / 2) as i128).abs() as u64)
            .max(max_deviation_jamhash);
    }

    max_deviation_jamhash
}
