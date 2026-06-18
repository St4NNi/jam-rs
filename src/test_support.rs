//! Test data shared by all tests: seeded random bases, copies with substitutions, and reverse
//! complements. Every test that needs random data uses these, with a fixed seed.

use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

/// A random generator that gives the same values for the same seed on every run.
pub fn rng(seed: u64) -> StdRng {
    StdRng::seed_from_u64(seed)
}

/// `length` random bases, each of `A`, `C`, `G` and `T` with the same chance.
pub fn random_bases(rng: &mut StdRng, length: usize) -> Vec<u8> {
    (0..length)
        .map(|_| b"ACGT"[rng.random_range(0..4)])
        .collect()
}

/// A copy of `bases` in which each base becomes another base with chance `1 - identity`.
pub fn with_substitutions(rng: &mut StdRng, bases: &[u8], identity: f64) -> Vec<u8> {
    bases
        .iter()
        .map(|&base| {
            if rng.random::<f64>() < identity {
                return base;
            }
            let others: Vec<u8> = b"ACGT"
                .iter()
                .copied()
                .filter(|&other| other != base)
                .collect();
            others[rng.random_range(0..others.len())]
        })
        .collect()
}

/// The reverse complement of `bases`.
pub fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .map(|&base| crate::align::complement(base))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_same_seed_gives_the_same_bases() {
        assert_eq!(random_bases(&mut rng(7), 50), random_bases(&mut rng(7), 50));
        assert_ne!(random_bases(&mut rng(7), 50), random_bases(&mut rng(8), 50));
    }

    #[test]
    fn substitutions_follow_the_identity() {
        let mut rng = rng(1);
        let bases = random_bases(&mut rng, 100_000);
        let copy = with_substitutions(&mut rng, &bases, 0.9);
        let same = bases.iter().zip(&copy).filter(|(a, b)| a == b).count();
        assert!((89_000..91_000).contains(&same), "{same}");
        assert_eq!(with_substitutions(&mut rng, &bases, 1.0), bases);
    }

    #[test]
    fn reverse_complement_reverses_and_complements() {
        assert_eq!(reverse_complement(b"AACGTN"), b"NACGTT");
    }
}
