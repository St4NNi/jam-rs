//! A list of hash functions to compare
//!
//! Constants chosen by testing different digits of pi;
use crate::cli::HashAlgorithms;
const KEY1: u64 = 0xe121_19c4_114f_22a7; // = 0x4528_21e6_38d0_1377 ^ 0xa409_3822_299f_31d0;
const KEY2: u32 = 0x60e5; //(0xbe54_66cf_34e9_0c6c ^ 0x082e_fa98_ec4e_6c89) & 63;

// Standard xxhash function for all sizes
#[inline]
pub fn xxhash3(kmer: &[u8]) -> u64 {
    xxhash_rust::xxh3::xxh3_64(kmer)
}

// Standard xxhash function for all sizes
#[inline]
pub fn xxhash3_u64(kmer: u64) -> u64 {
    xxhash_rust::xxh3::xxh3_64(&kmer.to_be_bytes())
}

// Specialized hash function for kmers < 32
// Simplified version of ahash-fallback from the ahash crate
#[inline]
pub fn ahash(kmer: u64) -> u64 {
    let temp = (kmer ^ KEY1) as u128 * 6364136223846793005_u128;
    let temp2 = ((temp & 0xffff_ffff_ffff_ffff) as u64) ^ ((temp >> 64) as u64); // XOR the lower 64 bits with the upper 64 bits.
    return temp2.rotate_left(KEY2);
}

// Faster version of murmur3 with equivalent output
#[inline]
pub fn murmur3(kmer: &[u8]) -> u64 {
    fastmurmur3::murmur3_x64_128(kmer, 42) as u64
}

#[inline]
pub fn murmur3_u64(kmer: u64) -> u64 {
    fastmurmur3::murmur3_x64_128(&kmer.to_be_bytes(), 42) as u64
}

/// Stores a function pointer to a hash function
#[derive(Clone)]
pub enum Function<'a> {
    Large(&'a (dyn Fn(&[u8]) -> u64 + Send + Sync)),
    Small(&'a (dyn Fn(u64) -> u64 + Send + Sync)),
}

impl Function<'_> {
    pub fn get_large(&self) -> Option<&dyn Fn(&[u8]) -> u64> {
        match self {
            Function::Large(f) => Some(f),
            _ => None,
        }
    }
    pub fn get_small(&self) -> Option<&dyn Fn(u64) -> u64> {
        match self {
            Function::Small(f) => Some(f),
            _ => None,
        }
    }

    pub fn from_alg(algo: HashAlgorithms, kmer_size: u8) -> Self {
        if kmer_size < 32 {
            match algo {
                HashAlgorithms::Ahash => Function::Small(&ahash),
                HashAlgorithms::Murmur3 => Function::Small(&murmur3_u64),
                HashAlgorithms::Xxhash => Function::Small(&xxhash3_u64),
                HashAlgorithms::Default => Function::Small(&ahash),
            }
        } else {
            match algo {
                HashAlgorithms::Murmur3 => Function::Large(&murmur3),
                HashAlgorithms::Xxhash | HashAlgorithms::Default => Function::Large(&xxhash3),
                _ => panic!("Hash function not supported for kmer size > 32"),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_xxhash3() {
        assert_eq!(xxhash3(b"AAAAAAAAAAA"), 0x92994E9987384EE2);
    }

    #[test]
    fn test_ahash() {
        assert_eq!(ahash(0xAAAAAAAAAAAAAAA), 6369629604220809163);
    }

    #[test]
    fn test_murmur3() {
        assert_eq!(murmur3(b"AAAAAAAAAAA"), 7773142420371383521);
    }

    #[test]
    fn test_xxhash3_u64() {
        assert_eq!(xxhash3_u64(0xAAAAAAAAAAAAAAA), 5855080426738543665);
    }

    #[test]
    fn test_murmur3_u64() {
        assert_eq!(murmur3_u64(0xAAAAAAAAAAAAAAA), 442865051503200633);
    }

    #[test]
    fn function_test() {
        let f = Function::from_alg(HashAlgorithms::Ahash, 21);
        assert_eq!(f.get_small().unwrap()(0xAAAAAAAAAAAAAAA), 6369629604220809163);
        let f = Function::from_alg(HashAlgorithms::Murmur3, 21);
        assert_eq!(f.get_small().unwrap()(0xAAAAAAAAAAAAAAA), 442865051503200633);
        let f = Function::from_alg(HashAlgorithms::Xxhash, 21);
        assert_eq!(f.get_small().unwrap()(0xAAAAAAAAAAAAAAA), 5855080426738543665);
        let f = Function::from_alg(HashAlgorithms::Default, 21);
        assert_eq!(f.get_small().unwrap()(0xAAAAAAAAAAAAAAA), 6369629604220809163);
        let f = Function::from_alg(HashAlgorithms::Murmur3, 32);
        assert_eq!(f.get_large().unwrap()(b"AAAAAAAAAAA"), 7773142420371383521);
        let f = Function::from_alg(HashAlgorithms::Xxhash, 32);
        assert_eq!(f.get_large().unwrap()(b"AAAAAAAAAAA"), 10563560822279786210);
        let f = Function::from_alg(HashAlgorithms::Default, 32);
        assert_eq!(f.get_large().unwrap()(b"AAAAAAAAAAA"), 10563560822279786210);
    }
}
