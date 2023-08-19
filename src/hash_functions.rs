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
        let default = if kmer_size < 32 {
            Function::Small(&ahash)
        } else {
            Function::Large(&xxhash3)
        };

        match algo {
            HashAlgorithms::Ahash => Function::Small(&ahash),
            HashAlgorithms::Murmur3 => Function::Large(&murmur3),
            HashAlgorithms::Xxhash => Function::Large(&xxhash3),
            HashAlgorithms::Default => default,
        }
    }
}
