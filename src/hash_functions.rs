//! Constants chosen by testing different digits of pi;
const KEY1: u64 = 0xe121_19c4_114f_22a7; // = 0x4528_21e6_38d0_1377 ^ 0xa409_3822_299f_31d0;
const KEY2: u32 = 0x60e5; //(0xbe54_66cf_34e9_0c6c ^ 0x082e_fa98_ec4e_6c89) & 63;

// Specialized hash function for kmers < 32
// Simplified version of ahash-fallback from the ahash crate
#[inline]
pub fn ahash(kmer: u64) -> u64 {
    let temp = (kmer ^ KEY1) as u128 * 6364136223846793005_u128;
    let temp2 = ((temp & 0xffff_ffff_ffff_ffff) as u64) ^ ((temp >> 64) as u64); // XOR the lower 64 bits with the upper 64 bits.
    temp2.rotate_left(KEY2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ahash() {
        assert_eq!(ahash(0xAAAAAAAAAAAAAAA), 6369629604220809163);
    }
}
