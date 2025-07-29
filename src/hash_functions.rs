//! Constants chosen by testing different digits of pi;
const KEY1: u64 = 0xe12119c4114f22a7; // = 0x4528_21e6_38d0_1377 ^ 0xa409_3822_299f_31d0;
const KEY2: u128 = 0xd1310ba698dfb5ac;
const PRNG_CONSTANT: u128 = 6364136223846793005_u128; // Well known constant from Donald Knuth's MMIX PRNG.
const MASK: u128 = 0xffff_ffff_ffff_ffff; // Mask for the lower 64 bits of a u128.

// Specialized hash function for bitkmers < 32 and a constant u64 input
// Inspired by the ahash fallback algorithm. https://github.com/tkaitchuck/aHash/wiki/AHash-fallback-algorithm
// and foldhash: https://github.com/orlp/foldhash
#[inline]
pub fn jamhash(kmer: u64) -> u64 {
    let temp = ((kmer ^ KEY1) as u128).wrapping_mul(PRNG_CONSTANT); // XOR the input with a constant, then multiply by the PRNG constant.
    let temp2 = ((temp & MASK) as u64) ^ ((temp >> 64) as u64); // XOR the lower 64 bits with the upper 64 bits.
    let temp3 = (temp2 as u128).wrapping_mul(KEY2);
    ((temp3 & MASK) as u64) ^ ((temp3 >> 64) as u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jamhash() {
        assert_eq!(jamhash(0xAAAAAAAAAAAAAAA), 5926757574219312021);
    }

    #[test]
    fn test_empty() {
        // This is expected, since it is the XOR constant XORed with itself.
        assert_eq!(jamhash(0xe121_19c4_114f_22a7), 0);
    }
}
