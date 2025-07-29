//! Constants chosen by testing different digits of pi;
const KEY1: u64 = 0x4327DE7F11A64EB7; // 0xd1310ba698dfb5ac ^ 0x9216d5d98979fb1b;
const KEY2: u64 = 0xe12119c4114f22a7;
//const KEY2: u128 = 0x4327DE7F11A64EB7; // 0xd1310ba698dfb5ac ^ 0x9216d5d98979fb1b;
const PRNG_CONSTANT: u128 = 6364136223846793005_u128; // Well known constant from Donald Knuth's MMIX PRNG.
const MASK: u128 = 0xffff_ffff_ffff_ffff; // Mask for the lower 64 bits of a u128.

// Specialized hash function for bitkmers < 32 and a constant u64 input
// Inspired by the ahash fallback algorithm. https://github.com/tkaitchuck/aHash/wiki/AHash-fallback-algorithm
// and foldhash: https://github.com/orlp/foldhash
#[inline]
pub fn jamhash(kmer: u64) -> u64 {
    let temp = ((kmer ^ KEY1) as u128).wrapping_mul(PRNG_CONSTANT); // XOR the input with a constant, then multiply by the PRNG constant.
    let temp2 = ((temp & MASK) as u64) ^ ((temp >> 64) as u64); // XOR the lower 64 bits with the upper 64 bits.
    let temp3 = ((temp2 ^ KEY2) as u128).wrapping_mul(PRNG_CONSTANT); // Second round to fully mix the bits
    ((temp3 & MASK) as u64) ^ ((temp3 >> 64) as u64) // XOR the lower 64 bits with the upper 64 bits again.
}

#[inline]
pub fn double_fold(input: u64, const_1: u64, const_2: u64) -> u64 {
    let first = fold_multiply(input, const_1);
    fold_multiply(first, const_2)
}

#[inline]
fn fold_multiply(input: u64, const_1: u64) -> u64 {
    let temp = ((input) as u128).wrapping_mul(const_1 as u128);
    ((temp & MASK) as u64) ^ ((temp >> 64) as u64)
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
