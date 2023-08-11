//! A list of hash functions to compare
const KEY1: u64 = 0x4528_21e6_38d0_1377 ^ 0x4;
const KEY2: u64 = 0xbe54_66cf_34e9_0c6c ^ 0x2;

// Standard xxhash function for all sizes
#[inline]
pub fn xxhash3(kmer: &[u8]) -> u64 {
    xxhash_rust::xxh3::xxh3_64(kmer)
}

// Specialized hash function for kmers < 32
#[inline]
pub fn ahash(kmer: u64) -> u64 {
    let temp = (kmer ^ KEY1) as u128 * 6364136223846793005_u128;
    let temp2 = ((temp & 0xffff_ffff_ffff_ffff) as u64) ^ ((temp >> 64) as u64); // XOR the lower 64 bits with the upper 64 bits.
    return temp2.rotate_left((KEY2 & 63) as u32);
}
