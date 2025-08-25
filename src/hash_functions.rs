//! Constants chosen by testing different digits of pi;
//! see tools/piconstants.py
const CONST1: u64 = 0xb8e1afed6a267e96;
const CONST2: u64 = 0x082efa98ec4e6c89;

// Specialized hash function for bitkmers < 32 and a constant u64 input
// Inspired by the ahash fallback algorithm. https://github.com/tkaitchuck/aHash/wiki/AHash-fallback-algorithm
// and foldhash: https://github.com/orlp/foldhash
#[inline]
pub fn jamhash(kmer: u64) -> u64 {
    let part1 = kmer.rotate_left(16);
    let part2 = kmer.rotate_left(48);
    let fold1 = fold_multiply(part1, CONST1);
    let fold2 = fold_multiply(part2, CONST2);
    fold_multiply(fold1, fold2)
}

#[inline(always)]
pub fn double_fold(input: u64, const_1: u64, const_2: u64) -> u64 {
    let first = fold_multiply(input, const_1);
    fold_multiply(first, const_2)
}

#[inline(always)]
fn fold_multiply(input: u64, const_1: u64) -> u64 {
    let temp = (input as u128).wrapping_mul(const_1 as u128);
    (temp as u64) ^ ((temp >> 64) as u64)
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
