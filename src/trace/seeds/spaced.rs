//! Reverse-complement-symmetric spaced exact seeds.
//!
//! The spaced scheme has a 31-base span and 21 informative bases.  The
//! selected bases are packed and canonicalized independently of the skipped
//! positions, so a digest collision can never be treated as an exact seed
//! without comparing this packed verification value.

use crate::archive::SeedKey;
use crate::jamhash_u64_v1;
use serde::{Deserialize, Serialize};

pub const SPACED_SPAN: u8 = 31;
pub const SPACED_WEIGHT: u8 = 21;

/// A symmetric 31-span/21-weight mask.  Bit zero is the left-most base.
///
/// The mask intentionally includes both span ends and is not merely a
/// contiguous k=21 window.  Its reverse is identical, which makes canonical
/// selection and query/target orientation independent.
pub const SPACED31_WEIGHT21_MASK: u32 = 0x6aaf_faab;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SpacedSeedConfig {
    pub span: u8,
    pub weight: u8,
    pub mask: u32,
    pub scale: u64,
    /// Maximum emitted seeds.  Zero means no bound.
    pub max_seeds: u32,
}

impl SpacedSeedConfig {
    #[must_use]
    pub const fn span31_weight21() -> Self {
        Self {
            span: SPACED_SPAN,
            weight: SPACED_WEIGHT,
            mask: SPACED31_WEIGHT21_MASK,
            scale: 1,
            max_seeds: 0,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Serialize, Deserialize)]
pub struct SpacedSeed {
    pub position: u64,
    pub span: u8,
    pub weight: u8,
    pub mask: u32,
    pub hash: u64,
    /// Canonical packed selected bases, in canonical sequence orientation.
    pub packed_selected: u64,
    pub reverse: bool,
}

impl SpacedSeed {
    #[must_use]
    pub fn verification_bytes(self) -> Vec<u8> {
        packed_bytes(self.packed_selected, self.weight)
    }

    #[must_use]
    pub fn seed_key(self) -> SeedKey {
        SeedKey {
            digest: self.hash,
            verification: self.verification_bytes(),
        }
    }

    pub fn selected_offsets(self) -> impl Iterator<Item = u8> {
        selected_offsets(self.mask, self.span)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SpacedSeedError {
    InvalidSpan(u8),
    InvalidWeight { span: u8, weight: u8 },
    InvalidMask { span: u8, mask: u32 },
    NonSymmetricMask(u32),
    ZeroScale,
    CoordinateOverflow,
}

impl std::fmt::Display for SpacedSeedError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidSpan(span) => write!(formatter, "spaced span must be 21..=31, got {span}"),
            Self::InvalidWeight { span, weight } => {
                write!(
                    formatter,
                    "spaced weight {weight} is invalid for span {span}"
                )
            }
            Self::InvalidMask { span, mask } => {
                write!(formatter, "spaced mask {mask:#x} exceeds span {span}")
            }
            Self::NonSymmetricMask(mask) => {
                write!(
                    formatter,
                    "spaced mask {mask:#x} is not reverse-complement symmetric"
                )
            }
            Self::ZeroScale => formatter.write_str("spaced seed scale must be greater than zero"),
            Self::CoordinateOverflow => formatter.write_str("spaced seed coordinate overflow"),
        }
    }
}

impl std::error::Error for SpacedSeedError {}

/// Extract exact spaced seeds from both represented strands.
pub fn extract_spaced_seeds(
    sequence: &[u8],
    config: SpacedSeedConfig,
) -> Result<Vec<SpacedSeed>, SpacedSeedError> {
    validate_config(config)?;
    let span = usize::from(config.span);
    if sequence.len() < span {
        return Ok(Vec::new());
    }
    let threshold = u64::MAX / config.scale;
    let mut seeds = Vec::new();
    for start in 0..=sequence.len() - span {
        let window = &sequence[start..start + span];
        let Some((forward, reverse_complement)) = packed_selected_pair(window, config) else {
            continue;
        };
        let (packed_selected, reverse) = if reverse_complement < forward {
            (reverse_complement, true)
        } else {
            (forward, false)
        };
        let hash = jamhash_u64_v1(packed_selected);
        if hash == 0 || hash >= threshold {
            continue;
        }
        seeds.push(SpacedSeed {
            position: u64::try_from(start).map_err(|_| SpacedSeedError::CoordinateOverflow)?,
            span: config.span,
            weight: config.weight,
            mask: config.mask,
            hash,
            packed_selected,
            reverse,
        });
        if config.max_seeds != 0
            && seeds.len() >= usize::try_from(config.max_seeds).unwrap_or(usize::MAX)
        {
            break;
        }
    }
    Ok(seeds)
}

pub fn extract_spaced_seeds_bounded(
    sequence: &[u8],
    config: SpacedSeedConfig,
) -> Result<Vec<SpacedSeed>, SpacedSeedError> {
    extract_spaced_seeds(sequence, config)
}

pub fn extract_spaced_seed_level(
    sequence: &[u8],
    config: SpacedSeedConfig,
) -> Result<Vec<SpacedSeed>, SpacedSeedError> {
    extract_spaced_seeds(sequence, config)
}

fn validate_config(config: SpacedSeedConfig) -> Result<(), SpacedSeedError> {
    if !(21..=31).contains(&config.span) {
        return Err(SpacedSeedError::InvalidSpan(config.span));
    }
    if config.weight == 0 || config.weight > config.span {
        return Err(SpacedSeedError::InvalidWeight {
            span: config.span,
            weight: config.weight,
        });
    }
    if config.span < 32 && config.mask >> u32::from(config.span) != 0 {
        return Err(SpacedSeedError::InvalidMask {
            span: config.span,
            mask: config.mask,
        });
    }
    if config.mask.count_ones() != u32::from(config.weight) {
        return Err(SpacedSeedError::InvalidWeight {
            span: config.span,
            weight: config.weight,
        });
    }
    if reverse_mask(config.mask, config.span) != config.mask {
        return Err(SpacedSeedError::NonSymmetricMask(config.mask));
    }
    if config.scale == 0 {
        return Err(SpacedSeedError::ZeroScale);
    }
    Ok(())
}

fn packed_selected_pair(window: &[u8], config: SpacedSeedConfig) -> Option<(u64, u64)> {
    let offsets = selected_offsets(config.mask, config.span).collect::<Vec<_>>();
    let mut forward = 0u64;
    for &offset in &offsets {
        forward = (forward << 2) | u64::from(base_bits(window[usize::from(offset)])?);
    }

    // A reverse-complement window visits selected positions in reverse order.
    // The symmetric mask means this is also the canonical packed order for
    // the reverse strand; no host-endian conversion is involved.
    let mut reverse = 0u64;
    for &offset in offsets.iter().rev() {
        let base = base_bits(window[usize::from(offset)])?;
        reverse = (reverse << 2) | u64::from(3 - base);
    }
    Some((forward, reverse))
}

pub fn selected_offsets(mask: u32, span: u8) -> impl Iterator<Item = u8> {
    (0..span).filter(move |offset| mask & (1u32 << u32::from(*offset)) != 0)
}

fn reverse_mask(mask: u32, span: u8) -> u32 {
    let mut reversed = 0u32;
    for offset in 0..span {
        if mask & (1u32 << u32::from(offset)) != 0 {
            reversed |= 1u32 << u32::from(span - 1 - offset);
        }
    }
    reversed
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

fn packed_bytes(value: u64, weight: u8) -> Vec<u8> {
    let width = usize::from(weight).div_ceil(4);
    let mut bytes = value.to_be_bytes().to_vec();
    bytes.drain(..bytes.len() - width);
    bytes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_mask_is_symmetric_and_has_requested_weight() {
        let config = SpacedSeedConfig::span31_weight21();
        assert_eq!(config.mask.count_ones(), 21);
        assert_eq!(reverse_mask(config.mask, 31), config.mask);
        assert_eq!(selected_offsets(config.mask, 31).count(), 21);
    }

    #[test]
    fn reverse_complement_has_same_seed_count() {
        let query = b"ACGTCAGTACCGATGCTAGCTAGGCTAACGTTACGATCGATGCATGCATGCA";
        let mut reverse = query.to_vec();
        reverse.reverse();
        for base in &mut reverse {
            *base = match *base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => unreachable!(),
            };
        }
        let config = SpacedSeedConfig::span31_weight21();
        let left = extract_spaced_seeds(query, config).unwrap();
        let right = extract_spaced_seeds(&reverse, config).unwrap();
        assert_eq!(left.len(), right.len());
    }

    #[test]
    fn ambiguity_is_not_exact_evidence() {
        let config = SpacedSeedConfig::span31_weight21();
        let query = b"ACGTCAGTACCGATGCTAGCTAGGCTAACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        let seeds = extract_spaced_seeds(query, config).unwrap();
        assert!(seeds.is_empty());
    }

    #[test]
    fn collision_verification_material_differs_from_digest() {
        let config = SpacedSeedConfig::span31_weight21();
        let seed = extract_spaced_seeds(b"ACGTCAGTACCGATGCTAGCTAGGCTAACGTTACGATCGATGCA", config)
            .unwrap()
            .into_iter()
            .next()
            .unwrap();
        assert_eq!(seed.verification_bytes().len(), 6);
        assert_ne!(seed.hash, 0);
        assert_ne!(
            seed.hash.to_be_bytes().as_slice(),
            seed.verification_bytes()
        );
    }
}
