use rand::{Rng, rng};
use roaring::RoaringTreemap;

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
    println!("{hash:10} {d: <10.10}");
}

#[inline]
pub fn murmur3_old(kmer: &[u8]) -> u64 {
    murmurhash3::murmurhash3_x64_128(kmer, 42).0
}

#[inline]
pub fn murmur3_new(kmer: &[u8]) -> u64 {
    fastmurmur3::murmur3_x64_128(kmer, 42) as u64
}

#[inline]
pub fn xxhash3(kmer: &[u8]) -> u64 {
    xxhash_rust::xxh3::xxh3_64(kmer)
}

const RANGE: std::ops::Range<u64> = 0..100_000_000u64;
const RANGE_LEN: usize = RANGE.end as usize - RANGE.start as usize;

// Skip
#[test]
#[ignore]
fn run_collision_analysis() {
    let mut hll_custom = RoaringTreemap::new();
    let mut hll_murmur3 = RoaringTreemap::new();
    let mut hll_xxhash = RoaringTreemap::new();

    let full_range = 0..4_398_046_511_104u64; // 2^42, the range of u64s.

    for x in full_range.take(1_000_000) {
        let xx = xxhash3(&x.to_be_bytes());
        let mo = murmur3_new(&x.to_be_bytes());
        let ah = jam_rs::hash_functions::jamhash(x);
        if !hll_custom.insert(ah) {
            panic!("Custom hash collision detected for value: {x}");
        };
        if !hll_murmur3.insert(mo) {
            panic!("Murmur3 hash collision detected for value: {x}");
        };
        if !hll_xxhash.insert(xx) {
            panic!("XXHash3 hash collision detected for value: {x}");
        };
    }
}

#[test]
#[ignore]
fn find_constants() {
    let mut rng = rng();
    let samples = (0..1_000_000).collect::<Vec<_>>();

    let mut results = vec![0u64; samples.len()];

    // 04feab043089c9a1, 4c6102befcb61250, 430a4fbd481e7756, f640d38a7b661335
    let mut smallest_ks = std::f64::MAX;
    loop {
        let const1 = rng.random::<u64>();
        let const2 = rng.random::<u64>();

        for x in 0..samples.len() {
            results[x] = jam_rs::hash_functions::double_fold(samples[x], const1, const2);
        }

        results.sort();
        let result = ks(&results);
        if result < smallest_ks {
            smallest_ks = result;
            println!("New best constants: {const1:016x}, {const2:016x}, with KS={smallest_ks}");
        }
    }
}

#[test]
fn run_ks() {
    let samples = RANGE.collect::<Vec<_>>();

    let samples_bytes = samples
        .iter()
        .map(|x| x.to_be_bytes().to_vec())
        .collect::<Vec<_>>();
    print_ks(
        "xxhash3",
        ks(&do_hashes_bytes(xxhash3, samples_bytes.as_slice())),
    );
    print_ks(
        "jamhash",
        ks(&do_hashes_u64(jam_rs::hash_functions::jamhash, &samples)),
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
    let mut xxhash3_bits = [0u64; 64];
    let mut jamhash_bits = [0u64; 64];
    let mut murmur3_old_bits = [0u64; 64];
    let mut murmur3_new_bits = [0u64; 64];

    for x in RANGE {
        let xx = xxhash3(x.to_be_bytes().as_slice());
        unrolled_64bits(xx, &mut xxhash3_bits);
        let ah = jam_rs::hash_functions::jamhash(x);
        unrolled_64bits(ah, &mut jamhash_bits);
        let mo = murmur3_old(x.to_be_bytes().as_slice());
        unrolled_64bits(mo, &mut murmur3_old_bits);
        let mn = murmur3_new(x.to_be_bytes().as_slice());
        unrolled_64bits(mn, &mut murmur3_new_bits);
    }

    let mut max_deviation_xxhash = 0i128;
    let mut max_deviation_jamhash = 0i128;
    let mut max_deviation_murmur3_old = 0i128;
    let mut max_deviation_murmur3_new = 0i128;
    for x in 0..64 {
        max_deviation_xxhash = (xxhash3_bits[x] as i128 - RANGE_LEN as i128)
            .abs()
            .max(max_deviation_xxhash);
        max_deviation_jamhash = (jamhash_bits[x] as i128 - RANGE_LEN as i128)
            .abs()
            .max(max_deviation_jamhash);
        max_deviation_murmur3_old = (murmur3_old_bits[x] as i128 - RANGE_LEN as i128)
            .abs()
            .max(max_deviation_murmur3_old);
        max_deviation_murmur3_new = (murmur3_new_bits[x] as i128 - RANGE_LEN as i128)
            .abs()
            .max(max_deviation_murmur3_new);
    }

    let max_deviation_xxhash = max_deviation_xxhash - (RANGE_LEN as i128 / 2);
    let max_deviation_jamhash = max_deviation_jamhash - (RANGE_LEN as i128 / 2);
    let max_deviation_murmur3_old = max_deviation_murmur3_old - (RANGE_LEN as i128 / 2);
    let max_deviation_murmur3_new = max_deviation_murmur3_new - (RANGE_LEN as i128 / 2);

    println!(
        "Max deviation from 50%: xxhash3: {max_deviation_xxhash}, jamhash: {max_deviation_jamhash}, murmur3_old: {max_deviation_murmur3_old}, murmur3_new: {max_deviation_murmur3_new}"
    );

    println!("bit|xxhash3|jamhash|murmur3_old|murmur3_new");
    for x in 0..64 {
        let xxhash_bit = xxhash3_bits[x] as f64 / RANGE_LEN as f64;
        let jamhash_bit = jamhash_bits[x] as f64 / RANGE_LEN as f64;
        let murmur3_old_bit = murmur3_old_bits[x] as f64 / RANGE_LEN as f64;
        let murmur3_new_bit = murmur3_new_bits[x] as f64 / RANGE_LEN as f64;
        assert!(xxhash_bit > 0.49);
        assert!(xxhash_bit < 0.51);
        assert!(jamhash_bit > 0.49);
        assert!(jamhash_bit < 0.51);
        assert!(murmur3_old_bit > 0.49);
        assert!(murmur3_old_bit < 0.51);
        assert!(murmur3_new_bit > 0.49);
        assert!(murmur3_new_bit < 0.51);
        println!("{x}|{xxhash_bit}|{jamhash_bit}|{murmur3_old_bit}|{murmur3_new_bit}");
    }
}

fn unrolled_64bits(num: u64, nums: &mut [u64; 64]) {
    for i in 0..64 {
        if num & (1u64 << i) != 0 {
            nums[i] += 1;
        }
    }
}
