use jamhash::jamhash_u64;
use rand::{Rng, rng};

#[allow(dead_code)]
fn cdf_uniform(x: u64) -> f64 {
    (x as f64) / (std::u64::MAX as f64)
}

#[allow(dead_code)]
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

fn compute_u64_avalanche() -> (f64, f64) {
    let mut rng = rng();

    let mut bit_flips = vec![0; 64 * 64];
    for _ in 0..10000 {
        let base_val: u64 = rng.random();
        let base_hash = jamhash_u64(base_val);
        for flip_pos in 0..64 {
            let delta_val = base_val ^ (1 << flip_pos);
            let delta_hash = jamhash_u64(delta_val);
            for test_pos in 0..64 {
                let flipped = ((base_hash ^ delta_hash) >> test_pos) & 1;
                bit_flips[test_pos * 64 + flip_pos] += flipped as usize;
            }
        }
    }

    let avg_bit_flips: f64 = bit_flips.iter().sum::<usize>() as f64 / bit_flips.len() as f64;

    let variance = bit_flips
        .iter()
        .map(|x| {
            let diff = *x as f64 - avg_bit_flips;
            diff * diff
        })
        .sum::<f64>()
        / bit_flips.len() as f64;

    let stddev = variance.sqrt();

    (avg_bit_flips, stddev)
}

fn main() {
    println!("Analyzing jamhash avalanche properties...");
    println!("Running 10 iterations of 10,000 samples each...\n");

    let mut results = Vec::new();
    for i in 0..10 {
        let (avg, stddev) = compute_u64_avalanche();
        println!(
            "Iteration {}: avg={:.4} (target: 5000), stddev={:.4}",
            i + 1,
            avg,
            stddev
        );
        results.push((avg, stddev));
    }

    let overall_avg: f64 = results.iter().map(|r| r.0).sum::<f64>() / results.len() as f64;
    let overall_stddev: f64 = results.iter().map(|r| r.1).sum::<f64>() / results.len() as f64;

    println!("\n=== Summary ===");
    println!("Overall average: {:.4} (target: 5000)", overall_avg);
    println!("Overall stddev: {:.4}", overall_stddev);
    println!(
        "Deviation from ideal: {:.4}%",
        ((overall_avg - 5000.0).abs() / 5000.0) * 100.0
    );

    if (overall_avg - 5000.0).abs() < 50.0 && overall_stddev < 100.0 {
        println!("\n✓ jamhash has excellent avalanche properties!");
    } else {
        println!("\n⚠ jamhash avalanche properties may need review");
    }
}
