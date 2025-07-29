//! Constants chosen by testing different digits of pi;
const KEY1: u64 = 0xe121_19c4_114f_22a7; // = 0x4528_21e6_38d0_1377 ^ 0xa409_3822_299f_31d0;
const PRNG_CONSTANT: u128 = 6364136223846793005_u128; // Well known constant from Donald Knuth's MMIX PRNG.
const MASK: u128 = 0xffff_ffff_ffff_ffff; // Mask for the lower 64 bits of a u128.

// Specialized hash function for bitkmers < 32 and a constant u64 input
// Inspired by the ahash fallback algorithm. https://github.com/tkaitchuck/aHash/wiki/AHash-fallback-algorithm
#[inline]
pub const fn jamhash(kmer: u64) -> u64 {
    let temp = ((kmer ^ KEY1) as u128).wrapping_mul(PRNG_CONSTANT); // XOR the input with a constant, then multiply by the PRNG constant.
    ((temp & MASK) as u64) ^ ((temp >> 64) as u64) // XOR the lower 64 bits with the upper 64 bits.
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jamhash() {
        assert_eq!(jamhash(0xAAAAAAAAAAAAAAA), 338655621601766339);
    }

    #[test]
    fn test_empty() {
        // This is expected, since it is the XOR constant XORed with itself.
        assert_eq!(jamhash(0xe121_19c4_114f_22a7), 0);
    }
}
