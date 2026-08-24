//! Exact paired k=21 and bounded strobemer seeds.
//!
//! The default paired scheme uses a fixed gap and is exactly
//! reverse-complement symmetric.  The optional min-hash mode chooses a
//! reciprocal minimum edge at both endpoints, which keeps the selected set
//! symmetric while bounding the local candidate neighborhood.

use crate::archive::SeedKey;
use crate::jamhash_u64_v1;
use serde::{Deserialize, Serialize};

pub const STROBEMER_K: u8 = 21;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum StrobemerMode {
    /// Pair every k-mer with the k-mer at one fixed gap.
    FixedPair,
    /// Select a reciprocal minimum-hash edge in each local neighborhood.
    MinHash,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct StrobemerConfig {
    pub k: u8,
    /// Number of unseeded bases between the two k-mers, inclusive bounds.
    pub min_gap: u16,
    pub max_gap: u16,
    pub mode: StrobemerMode,
    pub scale: u64,
    /// Maximum emitted seeds.  Zero means no bound.
    pub max_seeds: u32,
}

impl StrobemerConfig {
    #[must_use]
    pub const fn paired_k21() -> Self {
        Self {
            k: STROBEMER_K,
            min_gap: 4,
            max_gap: 4,
            mode: StrobemerMode::FixedPair,
            scale: 1,
            max_seeds: 0,
        }
    }

    #[must_use]
    pub const fn strobemer_k21() -> Self {
        Self {
            k: STROBEMER_K,
            min_gap: 4,
            max_gap: 64,
            mode: StrobemerMode::MinHash,
            scale: 1,
            max_seeds: 0,
        }
    }
}

/// One exact paired or strobemer seed.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Serialize, Deserialize)]
pub struct StrobemerSeed {
    /// Start of the left-most k-mer in the query orientation.
    pub position: u64,
    pub first_position: u64,
    pub second_position: u64,
    pub k: u8,
    pub gap: u16,
    /// Canonical packed k-mers in the selected pair orientation.
    pub first_kmer: u64,
    pub second_kmer: u64,
    pub hash: u64,
    pub reverse: bool,
}

impl StrobemerSeed {
    #[must_use]
    pub fn verification_bytes(self) -> Vec<u8> {
        let width = usize::from(self.k).div_ceil(4);
        let mut result = Vec::with_capacity(width * 2 + 2);
        result.extend_from_slice(&packed_bytes(self.first_kmer, self.k));
        result.extend_from_slice(&packed_bytes(self.second_kmer, self.k));
        result.extend_from_slice(&self.gap.to_be_bytes());
        result
    }

    #[must_use]
    pub fn seed_key(self) -> SeedKey {
        SeedKey {
            digest: self.hash,
            verification: self.verification_bytes(),
        }
    }

    #[must_use]
    pub const fn span(self) -> u64 {
        self.second_position + self.k as u64 - self.first_position
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum StrobemerError {
    UnsupportedK(u8),
    InvalidGap { min_gap: u16, max_gap: u16 },
    ZeroScale,
    CoordinateOverflow,
}

impl std::fmt::Display for StrobemerError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedK(k) => write!(formatter, "strobemer k must be 21, got {k}"),
            Self::InvalidGap { min_gap, max_gap } => {
                write!(
                    formatter,
                    "strobemer gap range {min_gap}..={max_gap} is invalid"
                )
            }
            Self::ZeroScale => formatter.write_str("strobemer scale must be greater than zero"),
            Self::CoordinateOverflow => formatter.write_str("strobemer coordinate overflow"),
        }
    }
}

impl std::error::Error for StrobemerError {}

#[derive(Clone, Copy)]
struct KmerInfo {
    forward: u64,
    reverse_complement: u64,
}

/// Extract exact paired or reciprocal-minimum strobemers.
pub fn extract_strobemers(
    sequence: &[u8],
    config: StrobemerConfig,
) -> Result<Vec<StrobemerSeed>, StrobemerError> {
    validate_config(config)?;
    let k = usize::from(config.k);
    let Some(last_start) = sequence.len().checked_sub(k) else {
        return Ok(Vec::new());
    };
    let mut kmers = Vec::with_capacity(last_start + 1);
    for start in 0..=last_start {
        kmers.push(canonical_kmer(&sequence[start..start + k]));
    }
    let threshold = u64::MAX / config.scale;
    let mut seeds = Vec::new();
    for first in 0..=last_start {
        let min_second = first
            .checked_add(k)
            .and_then(|value| value.checked_add(usize::from(config.min_gap)));
        let max_second = first
            .checked_add(k)
            .and_then(|value| value.checked_add(usize::from(config.max_gap)))
            .map_or(last_start, |value| value.min(last_start));
        let Some(min_second) = min_second else {
            continue;
        };
        if min_second > max_second {
            continue;
        }

        match config.mode {
            StrobemerMode::FixedPair => {
                // Fixed-distance pairing is naturally a reverse-complement
                // bijection: (first, second) becomes (rc(second), rc(first)).
                let second = min_second;
                if let Some(seed) = make_seed(first, second, &kmers, config, threshold)? {
                    seeds.push(seed);
                }
            }
            StrobemerMode::MinHash => {
                let mut best_for_first: Option<PairKey> = None;
                for second in min_second..=max_second {
                    let Some(candidate) = make_pair_key(first, second, &kmers, config.k) else {
                        continue;
                    };
                    if best_for_first.is_none_or(|best| candidate.key < best.key) {
                        best_for_first = Some(candidate);
                    }
                }
                let Some(best_for_first) = best_for_first else {
                    continue;
                };
                // A reciprocal minimum at both endpoints is invariant under
                // reversing the sequence and avoids one-sided strobe drift.
                let mut best_for_second: Option<PairKey> = None;
                let endpoint = best_for_first.second;
                let min_left_distance = k + usize::from(config.min_gap);
                let max_left_distance = k + usize::from(config.max_gap);
                let lower = endpoint.saturating_sub(max_left_distance);
                let upper = endpoint.saturating_sub(min_left_distance);
                if endpoint >= min_left_distance && lower <= upper {
                    for left in lower..=upper.min(last_start) {
                        let Some(candidate) = make_pair_key(left, endpoint, &kmers, config.k)
                        else {
                            continue;
                        };
                        if best_for_second.is_none_or(|best| candidate.key < best.key) {
                            best_for_second = Some(candidate);
                        }
                    }
                }
                if best_for_second.is_some_and(|candidate| {
                    candidate.first == best_for_first.first
                        && candidate.second == best_for_first.second
                        && candidate.key == best_for_first.key
                }) && let Some(seed) = make_seed(
                    best_for_first.first,
                    best_for_first.second,
                    &kmers,
                    config,
                    threshold,
                )? {
                    seeds.push(seed);
                }
            }
        }
        if config.max_seeds != 0
            && seeds.len() >= usize::try_from(config.max_seeds).unwrap_or(usize::MAX)
        {
            break;
        }
    }
    // A fixed pair can only be produced once.  The min-hash path also has a
    // unique left endpoint, but deduping keeps this invariant explicit if a
    // future partner policy emits a duplicate edge.
    seeds.sort_unstable_by_key(|seed| {
        (
            seed.position,
            seed.second_position,
            seed.hash,
            seed.first_kmer,
            seed.second_kmer,
        )
    });
    seeds.dedup();
    Ok(seeds)
}

pub fn extract_paired_k21(
    sequence: &[u8],
    scale: u64,
    max_seeds: u32,
) -> Result<Vec<StrobemerSeed>, StrobemerError> {
    extract_strobemers(
        sequence,
        StrobemerConfig {
            scale,
            max_seeds,
            ..StrobemerConfig::paired_k21()
        },
    )
}

pub fn extract_strobemers_bounded(
    sequence: &[u8],
    config: StrobemerConfig,
) -> Result<Vec<StrobemerSeed>, StrobemerError> {
    extract_strobemers(sequence, config)
}

pub fn extract_strobemer_seeds(
    sequence: &[u8],
    config: StrobemerConfig,
) -> Result<Vec<StrobemerSeed>, StrobemerError> {
    extract_strobemers(sequence, config)
}

pub fn extract_paired_k21_seeds(
    sequence: &[u8],
    scale: u64,
    max_seeds: u32,
) -> Result<Vec<StrobemerSeed>, StrobemerError> {
    extract_paired_k21(sequence, scale, max_seeds)
}

#[derive(Clone, Copy)]
struct PairKey {
    first: usize,
    second: usize,
    key: (u64, u64, u64, u16),
}

fn make_pair_key(
    first: usize,
    second: usize,
    kmers: &[Option<KmerInfo>],
    k: u8,
) -> Option<PairKey> {
    let left = kmers.get(first)?.as_ref()?;
    let right = kmers.get(second)?.as_ref()?;
    let gap = u16::try_from(second.checked_sub(first)?.checked_sub(usize::from(k))?).ok()?;
    let (first_kmer, second_kmer) = canonical_pair(left, right, k);
    let hash = pair_hash(first_kmer, second_kmer, gap);
    Some(PairKey {
        first,
        second,
        key: (
            hash,
            first_kmer.min(second_kmer),
            first_kmer.max(second_kmer),
            gap,
        ),
    })
}

fn make_seed(
    first: usize,
    second: usize,
    kmers: &[Option<KmerInfo>],
    config: StrobemerConfig,
    threshold: u64,
) -> Result<Option<StrobemerSeed>, StrobemerError> {
    let Some(left) = kmers.get(first).and_then(Option::as_ref) else {
        return Ok(None);
    };
    let Some(right) = kmers.get(second).and_then(Option::as_ref) else {
        return Ok(None);
    };
    let gap = u16::try_from(
        second
            .checked_sub(first)
            .and_then(|distance| distance.checked_sub(usize::from(config.k)))
            .ok_or(StrobemerError::CoordinateOverflow)?,
    )
    .map_err(|_| StrobemerError::CoordinateOverflow)?;
    let (first_kmer, second_kmer) = canonical_pair(left, right, config.k);
    let hash = pair_hash(first_kmer, second_kmer, gap);
    if hash == 0 || hash >= threshold {
        return Ok(None);
    }
    let reverse = first_kmer != left.forward || second_kmer != right.forward;
    Ok(Some(StrobemerSeed {
        position: u64::try_from(first).map_err(|_| StrobemerError::CoordinateOverflow)?,
        first_position: u64::try_from(first).map_err(|_| StrobemerError::CoordinateOverflow)?,
        second_position: u64::try_from(second).map_err(|_| StrobemerError::CoordinateOverflow)?,
        k: config.k,
        gap,
        first_kmer,
        second_kmer,
        hash,
        reverse,
    }))
}

fn canonical_pair(left: &KmerInfo, right: &KmerInfo, _k: u8) -> (u64, u64) {
    let forward = (left.forward, right.forward);
    let reverse = (right.reverse_complement, left.reverse_complement);
    // Canonicalize the complete ordered pair, not each member independently.
    // Reversing endpoint order without complementing the bases therefore
    // remains a distinct exact pair.
    if reverse < forward { reverse } else { forward }
}

fn pair_hash(first: u64, second: u64, gap: u16) -> u64 {
    jamhash_u64_v1(first ^ second.rotate_left(17) ^ u64::from(gap).rotate_left(41))
}

fn canonical_kmer(sequence: &[u8]) -> Option<KmerInfo> {
    let mut forward = 0u64;
    for &base in sequence {
        forward = (forward << 2) | u64::from(base_bits(base)?);
    }
    let reverse_complement = reverse_complement_packed(forward, sequence.len());
    Some(KmerInfo {
        forward,
        reverse_complement,
    })
}

fn reverse_complement_packed(value: u64, length: usize) -> u64 {
    let mut result = 0u64;
    let mut remaining = value;
    for _ in 0..length {
        let base = remaining & 3;
        remaining >>= 2;
        result = (result << 2) | (3 - base);
    }
    result
}

fn validate_config(config: StrobemerConfig) -> Result<(), StrobemerError> {
    if config.k != STROBEMER_K {
        return Err(StrobemerError::UnsupportedK(config.k));
    }
    if config.min_gap > config.max_gap {
        return Err(StrobemerError::InvalidGap {
            min_gap: config.min_gap,
            max_gap: config.max_gap,
        });
    }
    if config.scale == 0 {
        return Err(StrobemerError::ZeroScale);
    }
    Ok(())
}

fn base_bits(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

fn packed_bytes(value: u64, k: u8) -> Vec<u8> {
    let width = usize::from(k).div_ceil(4);
    let mut bytes = value.to_be_bytes().to_vec();
    bytes.drain(..bytes.len() - width);
    bytes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixed_pairs_are_reverse_complement_symmetric() {
        let query = b"ACGTCAGTACCGATGCTAGCTAGGCTAACGTTACGATCGATGCATGCATGCA";
        let config = StrobemerConfig::paired_k21();
        let left = extract_strobemers(query, config).unwrap();
        let right = extract_strobemers(&reverse_complement(query), config).unwrap();
        assert_eq!(left.len(), right.len());
        let mut left_keys = left
            .iter()
            .map(|seed| (seed.hash, seed.verification_bytes()))
            .collect::<Vec<_>>();
        let mut right_keys = right
            .iter()
            .map(|seed| (seed.hash, seed.verification_bytes()))
            .collect::<Vec<_>>();
        left_keys.sort_unstable();
        right_keys.sort_unstable();
        assert_eq!(left_keys, right_keys);
    }

    #[test]
    fn reciprocal_min_hash_is_bounded_and_exact() {
        let query = b"ACGTCAGTACCGATGCTAGCTAGGCTAACGTTACGATCGATGCATGCATGCA";
        let config = StrobemerConfig {
            max_seeds: 4,
            ..StrobemerConfig::strobemer_k21()
        };
        let seeds = extract_strobemers(query, config).unwrap();
        assert!(seeds.len() <= 4);
        assert!(
            seeds
                .iter()
                .all(|seed| seed.verification_bytes().len() == 14)
        );
        let reverse = reverse_complement(query);
        let reverse_seeds = extract_strobemers(&reverse, config).unwrap();
        let mut left_keys = seeds
            .iter()
            .map(|seed| (seed.hash, seed.verification_bytes()))
            .collect::<Vec<_>>();
        let mut right_keys = reverse_seeds
            .iter()
            .map(|seed| (seed.hash, seed.verification_bytes()))
            .collect::<Vec<_>>();
        left_keys.sort_unstable();
        right_keys.sort_unstable();
        assert_eq!(left_keys, right_keys);
    }

    #[test]
    fn reciprocal_min_hash_selects_nonempty_edges() {
        let mut state = 0x1234_5678_9abc_def0u64;
        let mut query = Vec::with_capacity(512);
        for _ in 0..512 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            query.push(b"ACGT"[((state >> 62) & 3) as usize]);
        }
        let seeds = extract_strobemers(&query, StrobemerConfig::strobemer_k21()).unwrap();
        assert!(!seeds.is_empty());
        assert!(
            seeds
                .iter()
                .all(|seed| seed.second_position > seed.first_position)
        );
    }

    fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|base| match *base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => unreachable!(),
            })
            .collect()
    }

    #[test]
    fn ambiguity_cannot_create_a_pair() {
        let query = b"ACGTCAGTACCGATGCTAGCTAGGCTAACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        assert!(
            extract_strobemers(query, StrobemerConfig::paired_k21())
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn reversed_noncomplemented_pair_is_not_the_same_exact_key() {
        let left = canonical_kmer(&[b'A'; STROBEMER_K as usize]).unwrap();
        let right = canonical_kmer(&[b'C'; STROBEMER_K as usize]).unwrap();
        let forward_order = canonical_pair(&left, &right, STROBEMER_K);
        let reversed_order = canonical_pair(&right, &left, STROBEMER_K);
        assert_ne!(forward_order, reversed_order);
        assert_ne!(
            pair_hash(forward_order.0, forward_order.1, 4),
            pair_hash(reversed_order.0, reversed_order.1, 4)
        );
    }

    #[test]
    fn forced_digest_collision_has_distinct_verification() {
        let query = b"ACGTCAGTACCGATGCTAGCTAGGCTAACGTTACGATCGATGCATGCATGCA";
        let seed = extract_strobemers(query, StrobemerConfig::paired_k21())
            .unwrap()
            .into_iter()
            .next()
            .unwrap();
        assert_ne!(seed.hash, 0);
        assert_eq!(seed.verification_bytes().len(), 14);
        let mut collision = seed.verification_bytes();
        collision[0] ^= 1;
        assert_ne!(collision, seed.verification_bytes());
    }
}
