// kolmogorov-smirnov from: https://github.com/tmmcguire/hashers/blob/master/examples/kolmogorov-smirnov.rs
// See
// - https://www.itl.nist.gov/div898/handbook/eda/section3/eda35g.htm
// - https://onlinecourses.science.psu.edu/stat414/node/322/

/// Hash a sequence of values, returning the hashes sorted.
#[inline]
fn do_hashes_bytes(fcn: fn(&[u8]) -> u64, data: &[Vec<u8>]) -> Vec<u64> {
    let mut res: Vec<u64> = data.iter().map(|elt| fcn(elt)).collect();
    res.sort();
    res
}

#[inline]
fn do_hashes_u64(fcn: fn(u64) -> u64, data: &[u64]) -> Vec<u64> {
    let mut res: Vec<u64> = data.iter().map(|elt| fcn(*elt)).collect();
    res.sort();
    res
}

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

fn print_ks(hash: &str, d: f64) {
    assert!(d < 0.005); // 0.5% confidence interval that the distribution is not uniform.
    println!("{:10} {: <10.10}", hash, d);
}

#[inline]
pub fn murmur3_old(kmer: &[u8]) -> u64 {
    murmurhash3::murmurhash3_x64_128(kmer, 42).0
}

#[inline]
pub fn murmur3_new(kmer: &[u8]) -> u64 {
    fastmurmur3::murmur3_x64_128(kmer, 42) as u64
}

#[test]
fn run_ks() {
    let samples = (100_000_000_000..100_000_100_000u64).collect::<Vec<_>>();

    let samples_bytes = samples
        .iter()
        .map(|x| x.to_be_bytes().to_vec())
        .collect::<Vec<_>>();
    print_ks(
        "xxhash3",
        ks(&do_hashes_bytes(
            jam_rs::hash_functions::xxhash3,
            samples_bytes.as_slice(),
        )),
    );
    print_ks(
        "ahash",
        ks(&do_hashes_u64(jam_rs::hash_functions::ahash, &samples)),
    );
    print_ks(
        "murmur3_old",
        ks(&do_hashes_bytes(murmur3_old, &samples_bytes)),
    );
    print_ks(
        "murmur3_new",
        ks(&do_hashes_bytes(murmur3_new, &samples_bytes)),
    );
}

#[test]
fn test_bit_distribution() {
    let samples = (100_000_000_000..100_010_000_000u64).collect::<Vec<_>>();

    let samples_bytes = samples
        .iter()
        .map(|x| x.to_be_bytes().to_vec())
        .collect::<Vec<_>>();

    let mut xxhash3_bits = [0u64; 64];
    let mut ahash_bits = [0u64; 64];
    let mut murmur3_old_bits = [0u64; 64];
    let mut murmur3_new_bits = [0u64; 64];

    for x in 0..samples.len() {
        let xx = jam_rs::hash_functions::xxhash3(samples_bytes[x].as_slice());
        unrolled_64bits(xx, &mut xxhash3_bits);
        let ah = jam_rs::hash_functions::ahash(samples[x]);
        unrolled_64bits(ah, &mut ahash_bits);
        let mo = murmur3_old(samples_bytes[x].as_slice());
        unrolled_64bits(mo, &mut murmur3_old_bits);
        let mn = murmur3_new(samples_bytes[x].as_slice());
        unrolled_64bits(mn, &mut murmur3_new_bits);
    }

    println!("bit|xxhash3|ahash|murmur3_old|murmur3_new");
    for x in 0..64 {
        let xxhash_bit = xxhash3_bits[x] as f64 / 10_000_000f64;
        let ahash_bit = ahash_bits[x] as f64 / 10_000_000f64;
        let murmur3_old_bit = murmur3_old_bits[x] as f64 / 10_000_000f64;
        let murmur3_new_bit = murmur3_new_bits[x] as f64 / 10_000_000f64;
        assert!(xxhash_bit > 0.49);
        assert!(xxhash_bit < 0.51);
        assert!(ahash_bit > 0.49);
        assert!(ahash_bit < 0.51);
        assert!(murmur3_old_bit > 0.49);
        assert!(murmur3_old_bit < 0.51);
        assert!(murmur3_new_bit > 0.49);
        assert!(murmur3_new_bit < 0.51);
        println!(
            "{}|{}|{}|{}|{}",
            x, xxhash_bit, ahash_bit, murmur3_old_bit, murmur3_new_bit
        );
    }
}

fn unrolled_64bits(num: u64, nums: &mut [u64; 64]) {
    for i in 0..64 {
        if num & (1u64 << i) != 0 {
            nums[i] += 1;
        }
    }
}
