use num_bigint::{BigInt, ToBigInt};
use num_rational::BigRational;
use num_traits::{One, Pow, ToPrimitive, Zero};
use std::time::Instant;

fn hamming_distance(x: u64, y: u64) -> u32 {
    (x ^ y).count_ones()
}

fn triplet_distance(a: u64, b: u64, c: u64) -> u32 {
    hamming_distance(a, b) + hamming_distance(b, c) + hamming_distance(a, c)
}

// High-precision pi - 3 using Machin's formula for better accuracy
fn calculate_pi_minus_3(precision_digits: usize) -> BigRational {
    // Using Machin's formula: π/4 = 4*arctan(1/5) - arctan(1/239)
    // We'll compute this to high precision using power series

    fn arctan_series(x_inv: i64, terms: usize) -> BigRational {
        let x_inv_big = BigRational::from(BigInt::from(x_inv));
        let x_inv_squared = &x_inv_big * &x_inv_big;

        let mut result = BigRational::from(BigInt::one()) / &x_inv_big;
        let mut term = result.clone();

        for n in 1..terms {
            // Next term: (-1)^n * x^(2n+1) / (2n+1)
            term = &term / &x_inv_squared;
            let sign = if n % 2 == 0 { 1 } else { -1 };
            let denominator = 2 * n + 1;
            let next_term = &term * BigRational::from(BigInt::from(sign))
                / BigRational::from(BigInt::from(denominator));
            result = result + next_term;
        }

        result
    }

    let terms = precision_digits / 2 + 100; // More terms for higher precision

    // π/4 = 4*arctan(1/5) - arctan(1/239)
    let pi_over_4 =
        BigRational::from(BigInt::from(4)) * arctan_series(5, terms) - arctan_series(239, terms);
    let pi = BigRational::from(BigInt::from(4)) * pi_over_4;

    // Return π - 3
    pi - BigRational::from(BigInt::from(3))
}

fn generate_constants(range_size: usize) -> Vec<u64> {
    println!(
        "Generating {} constants with pure Rust precision...",
        range_size - 1
    );

    // Calculate π - 3 with high precision
    let start_calc = Instant::now();
    let a = calculate_pi_minus_3(1000); // 1000 decimal digits should be plenty
    println!(
        "Calculated π - 3 in {:.2}s",
        start_calc.elapsed().as_secs_f64()
    );

    let mut constants = Vec::with_capacity(range_size - 1);

    for i in 1..range_size {
        // Calculate 2^(64*i) as BigInt
        let power_exponent = 64 * i;
        let power_of_two = BigInt::from(2).pow(power_exponent);

        // a * 2^(64*i)
        let product = &a * BigRational::from(power_of_two);

        // Get the integer part (floor)
        let floored = product.to_integer();

        // Convert to u64 (this handles modulo 2^64 automatically via wrapping)
        let b = floored.to_u64().unwrap_or_else(|| {
            // For very large numbers, take modulo 2^64 manually
            let mask = (BigInt::one() << 64) - BigInt::one();
            let reduced = floored & mask;
            reduced.to_u64().unwrap_or(0)
        });

        constants.push(b);

        if i % 100 == 0 {
            println!("Generated {} constants...", i);
        }
    }

    constants
}

fn find_best_triplet_exhaustive(
    constants: &[u64],
) -> (u32, (usize, usize, usize), (u64, u64, u64)) {
    let n = constants.len();
    let total_triplets = (n * (n - 1) * (n - 2)) / 6;
    println!("Searching {} triplets exhaustively...", total_triplets);

    let mut max_distance = 0u32;
    let mut best_indices = (0, 0, 0);
    let mut best_values = (0u64, 0u64, 0u64);
    let mut processed = 0usize;

    let start_time = Instant::now();

    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                let distance = triplet_distance(constants[i], constants[j], constants[k]);

                if distance > max_distance {
                    max_distance = distance;
                    best_indices = (i, j, k);
                    best_values = (constants[i], constants[j], constants[k]);
                }

                processed += 1;
                if processed % 100_000 == 0 {
                    let elapsed = start_time.elapsed().as_secs_f64();
                    let rate = processed as f64 / elapsed;
                    let eta = (total_triplets - processed) as f64 / rate;
                    println!(
                        "Processed {:>10} triplets ({:>6.2}%) | Best: {} | Rate: {:>8.0}/sec | ETA: {:>6.1}s",
                        processed,
                        100.0 * processed as f64 / total_triplets as f64,
                        max_distance,
                        rate,
                        eta
                    );
                }
            }
        }
    }

    (max_distance, best_indices, best_values)
}

fn find_best_triplet_sampled(
    constants: &[u64],
    sample_size: usize,
) -> (u32, (usize, usize, usize), (u64, u64, u64)) {
    use rand::prelude::*;

    println!("Sampling {} random triplets...", sample_size);

    let mut rng = thread_rng();
    let mut max_distance = 0u32;
    let mut best_indices = (0, 0, 0);
    let mut best_values = (0u64, 0u64, 0u64);

    let start_time = Instant::now();

    for i in 0..sample_size {
        // Generate three unique random indices efficiently
        let n = constants.len();
        let idx1 = rng.gen_range(0..n);
        let mut idx2 = rng.gen_range(0..n - 1);
        if idx2 >= idx1 {
            idx2 += 1;
        }
        let mut idx3 = rng.gen_range(0..n - 2);
        if idx3 >= idx1.min(idx2) {
            idx3 += 1;
        }
        if idx3 >= idx1.max(idx2) {
            idx3 += 1;
        }

        let distance = triplet_distance(constants[idx1], constants[idx2], constants[idx3]);

        if distance > max_distance {
            max_distance = distance;
            best_indices = (idx1, idx2, idx3);
            best_values = (constants[idx1], constants[idx2], constants[idx3]);
        }

        if i % 10_000 == 0 && i > 0 {
            let elapsed = start_time.elapsed().as_secs_f64();
            let rate = i as f64 / elapsed;
            let eta = (sample_size - i) as f64 / rate;
            println!(
                "Processed {:>8} samples | Best: {} | Rate: {:>8.0}/sec | ETA: {:>6.1}s",
                i, max_distance, rate, eta
            );
        }
    }

    (max_distance, best_indices, best_values)
}

fn main() {
    println!("Pure Rust High-Performance Bit Difference Triplet Finder\n");

    // Configuration
    let range_size = 500; // Change this: 100, 200, 500, 1000
    let sample_size: Option<usize> = None; // Some(100_000) for sampling, None for exhaustive

    // Generate constants
    let start_time = Instant::now();
    let constants = generate_constants(range_size);
    println!(
        "Generated {} constants in {:.2}s\n",
        constants.len(),
        start_time.elapsed().as_secs_f64()
    );

    // Calculate total possible triplets
    let n = constants.len();
    let total_triplets = (n * (n - 1) * (n - 2)) / 6;
    println!("Total possible triplets: {}", total_triplets);

    // Find best triplet
    let search_start = Instant::now();
    let (max_distance, best_indices, best_values) = match sample_size {
        Some(samples) if samples < total_triplets => find_best_triplet_sampled(&constants, samples),
        _ => find_best_triplet_exhaustive(&constants),
    };

    let search_time = search_start.elapsed().as_secs_f64();
    println!("\nSearch completed in {:.2}s\n", search_time);

    // Display results
    println!("BEST TRIPLET FOUND:");
    println!("Indices: {:?}", best_indices);
    println!("Values:");
    for (&idx, &val) in best_indices.iter().zip(best_values.iter()) {
        println!("  constants[{}] = 0x{:016x}", idx, val);
    }

    println!("\nPairwise Hamming distances:");
    let d01 = hamming_distance(best_values.0, best_values.1);
    let d12 = hamming_distance(best_values.1, best_values.2);
    let d02 = hamming_distance(best_values.0, best_values.2);

    println!("  {} ↔ {}: {} bits", best_indices.0, best_indices.1, d01);
    println!("  {} ↔ {}: {} bits", best_indices.1, best_indices.2, d12);
    println!("  {} ↔ {}: {} bits", best_indices.0, best_indices.2, d02);
    println!("Total distance: {} bits", max_distance);

    println!("\nBit patterns:");
    for (&idx, &val) in best_indices.iter().zip(best_values.iter()) {
        println!("constants[{}]: {:064b}", idx, val);
    }

    println!("\nStatistics:");
    println!("Maximum possible distance per pair: 64 bits");
    println!("Maximum possible triplet distance: 192 bits");
    println!(
        "Achieved distance: {} bits ({:.1}% of maximum)",
        max_distance,
        max_distance as f64 / 192.0 * 100.0
    );

    // Verification: Print first few constants for comparison with Python
    println!("\nFirst 5 constants (hex) for verification:");
    for (i, &val) in constants.iter().take(5).enumerate() {
        println!("  constants[{}] = 0x{:016x}", i, val);
    }
}
