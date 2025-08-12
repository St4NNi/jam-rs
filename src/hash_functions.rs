//! Constants chosen by testing different digits of pi;
//! see tools/piconstants.py
const CONST1: u64 = 0xe8ddc98368c78e76; // 621
const CONST2: u64 = 0x9de2168e9de8b1d1; // 632
const CONST3: u64 = 0x160ffe7183376388; // 770

// Specialized hash function for bitkmers < 32 and a constant u64 input
// Inspired by the ahash fallback algorithm. https://github.com/tkaitchuck/aHash/wiki/AHash-fallback-algorithm
// and foldhash: https://github.com/orlp/foldhash
#[inline]
pub fn jamhash(kmer: u64) -> u64 {
    let kmer_const1 = kmer ^ CONST1;
    let kmer_const2 = kmer_const1 ^ CONST2;
    let fold1 = fold_multiply(kmer_const1, kmer_const2);
    fold_multiply(fold1, CONST3)
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
