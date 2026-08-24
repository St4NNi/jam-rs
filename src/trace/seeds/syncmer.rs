//! Orientation-symmetric syncmer seed extraction.
//!
//! Syncmers are kept in their own module because they are an experimental
//! trace seed scheme, rather than part of the FracMinHash levels used by the
//! current runner.  A syncmer is still an exact seed: the digest is only a
//! lookup hint and the canonical packed k-mer is retained for verification.

use crate::archive::SeedKey;
use crate::jamhash_u64_v1;
use serde::{Deserialize, Serialize};

/// The two bounded syncmer selection rules used by the trace experiments.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum SyncmerMode {
    /// Select when the minimum canonical s-mer is at either k-mer boundary.
    Closed,
    /// Select the boundary that is first in the canonical k-mer orientation.
    /// This keeps the lower-density open rule reverse-complement symmetric.
    Open,
}

/// Parameters for one syncmer level.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SyncmerConfig {
    pub k: u8,
    pub s: u8,
    pub mode: SyncmerMode,
    pub scale: u64,
    /// Maximum number of selected seeds.  Zero means no bound.
    pub max_seeds: u32,
}

impl SyncmerConfig {
    /// A conservative k=31 primary syncmer configuration.
    #[must_use]
    pub const fn k31(mode: SyncmerMode) -> Self {
        Self {
            k: 31,
            s: 11,
            mode,
            scale: 1,
            max_seeds: 0,
        }
    }

    /// A denser k=21 rescue syncmer configuration.
    #[must_use]
    pub const fn k21_rescue(mode: SyncmerMode) -> Self {
        Self {
            k: 21,
            s: 9,
            mode,
            scale: 1,
            max_seeds: 0,
        }
    }
}

/// One selected exact syncmer.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Serialize, Deserialize)]
pub struct SyncmerSeed {
    pub position: u64,
    pub k: u8,
    pub s: u8,
    pub mode: SyncmerMode,
    pub hash: u64,
    /// Canonical packed k-mer, not the selected s-mer.  This is the exact
    /// verification value required before an anchor can be accepted.
    pub canonical_kmer: u64,
    pub reverse: bool,
    /// Position of the selected s-mer in the forward query window.
    pub selected_offset: u8,
}

impl SyncmerSeed {
    /// Adapt a contiguous syncmer to the common exact k-mer query type.
    #[must_use]
    pub const fn query_seed(self) -> super::QuerySeed {
        super::QuerySeed {
            position: self.position,
            hash: self.hash,
            canonical_kmer: self.canonical_kmer,
            reverse: self.reverse,
        }
    }

    /// Return the fixed-width verification bytes in big-endian order.
    ///
    /// The width is derived from k and is stable across architectures.  The
    /// selected mode and s-mer position are scheme metadata, not lookup data.
    #[must_use]
    pub fn verification_bytes(self) -> Vec<u8> {
        packed_bytes(self.canonical_kmer, self.k)
    }

    /// Construct the backend-neutral exact lookup key.  The digest remains a
    /// candidate filter; the packed canonical k-mer is mandatory evidence.
    #[must_use]
    pub fn seed_key(self) -> SeedKey {
        SeedKey {
            digest: self.hash,
            verification: self.verification_bytes(),
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SyncmerError {
    UnsupportedK(u8),
    InvalidS { k: u8, s: u8 },
    ZeroScale,
    AmbiguousBase { position: u64 },
    CoordinateOverflow,
}

impl std::fmt::Display for SyncmerError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedK(k) => write!(formatter, "syncmer k must be 21 or 31, got {k}"),
            Self::InvalidS { k, s } => write!(formatter, "syncmer s={s} is invalid for k={k}"),
            Self::ZeroScale => formatter.write_str("syncmer scale must be greater than zero"),
            Self::AmbiguousBase { position } => {
                write!(formatter, "ambiguous base at syncmer position {position}")
            }
            Self::CoordinateOverflow => formatter.write_str("syncmer coordinate overflow"),
        }
    }
}

impl std::error::Error for SyncmerError {}

/// Extract selected syncmers in query order.
pub fn extract_syncmers(
    sequence: &[u8],
    config: SyncmerConfig,
) -> Result<Vec<SyncmerSeed>, SyncmerError> {
    validate_config(config)?;
    let k = usize::from(config.k);
    let s = usize::from(config.s);
    if sequence.len() < k {
        return Ok(Vec::new());
    }

    let mut seeds = Vec::new();
    let threshold = u64::MAX / config.scale;
    for start in 0..=sequence.len() - k {
        let window = &sequence[start..start + k];
        let Some((canonical_kmer, reverse)) = canonical_packed(window) else {
            continue;
        };
        let mut minimum = u64::MAX;
        let mut minimum_at_left = false;
        let mut minimum_at_right = false;
        for offset in 0..=k - s {
            let Some((canonical_smer, _)) = canonical_packed(&window[offset..offset + s]) else {
                // The complete k-mer was valid, so this can only be reached
                // if a future base codec changes its validity rules.
                continue;
            };
            if canonical_smer < minimum {
                minimum = canonical_smer;
                minimum_at_left = offset == 0;
                minimum_at_right = offset == k - s;
            } else if canonical_smer == minimum {
                minimum_at_left |= offset == 0;
                minimum_at_right |= offset == k - s;
            }
        }
        let selected = match config.mode {
            SyncmerMode::Closed => minimum_at_left || minimum_at_right,
            SyncmerMode::Open => {
                if reverse {
                    minimum_at_right
                } else {
                    minimum_at_left
                }
            }
        };
        if minimum == u64::MAX || !selected {
            continue;
        }
        let minimum_offset = if config.mode == SyncmerMode::Closed {
            if minimum_at_left { 0 } else { k - s }
        } else if reverse && minimum_at_right {
            k - s
        } else {
            0
        };
        let hash = jamhash_u64_v1(canonical_kmer);
        if hash == 0 || hash >= threshold {
            continue;
        }
        seeds.push(SyncmerSeed {
            position: u64::try_from(start).map_err(|_| SyncmerError::CoordinateOverflow)?,
            k: config.k,
            s: config.s,
            mode: config.mode,
            hash,
            canonical_kmer,
            reverse,
            selected_offset: u8::try_from(minimum_offset)
                .map_err(|_| SyncmerError::CoordinateOverflow)?,
        });
        if config.max_seeds != 0
            && seeds.len() >= usize::try_from(config.max_seeds).unwrap_or(usize::MAX)
        {
            break;
        }
    }
    Ok(seeds)
}

/// Alias used by the benchmark harness for a bounded extraction round.
pub fn extract_syncmers_bounded(
    sequence: &[u8],
    config: SyncmerConfig,
) -> Result<Vec<SyncmerSeed>, SyncmerError> {
    extract_syncmers(sequence, config)
}

/// Descriptive alias used by seed-comparison callers.
pub fn extract_syncmer_seeds(
    sequence: &[u8],
    config: SyncmerConfig,
) -> Result<Vec<SyncmerSeed>, SyncmerError> {
    extract_syncmers(sequence, config)
}

pub fn extract_syncmer_query_seeds(
    sequence: &[u8],
    config: SyncmerConfig,
) -> Result<Vec<super::QuerySeed>, SyncmerError> {
    Ok(extract_syncmers(sequence, config)?
        .into_iter()
        .map(SyncmerSeed::query_seed)
        .collect())
}

/// Extract both closed and canonical-orientation open levels with one pass
/// over the query.  The cap applies independently to each level.
pub fn extract_syncmer_levels(
    sequence: &[u8],
    k: u8,
    s: u8,
    scale: u64,
    max_seeds: u32,
) -> Result<(Vec<SyncmerSeed>, Vec<SyncmerSeed>), SyncmerError> {
    let closed = extract_syncmers(
        sequence,
        SyncmerConfig {
            k,
            s,
            mode: SyncmerMode::Closed,
            scale,
            max_seeds,
        },
    )?;
    let open = extract_syncmers(
        sequence,
        SyncmerConfig {
            k,
            s,
            mode: SyncmerMode::Open,
            scale,
            max_seeds,
        },
    )?;
    Ok((closed, open))
}

fn validate_config(config: SyncmerConfig) -> Result<(), SyncmerError> {
    if config.k != 21 && config.k != 31 {
        return Err(SyncmerError::UnsupportedK(config.k));
    }
    if config.s < 3 || config.s > config.k {
        return Err(SyncmerError::InvalidS {
            k: config.k,
            s: config.s,
        });
    }
    if config.scale == 0 {
        return Err(SyncmerError::ZeroScale);
    }
    Ok(())
}

fn canonical_packed(sequence: &[u8]) -> Option<(u64, bool)> {
    if sequence.len() > 31 {
        return None;
    }
    let mut forward = 0u64;
    for &base in sequence {
        forward = (forward << 2) | u64::from(base_bits(base)?);
    }
    let reverse = reverse_complement_packed(forward, sequence.len());
    if reverse < forward {
        Some((reverse, true))
    } else {
        Some((forward, false))
    }
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
    fn closed_and_open_rules_are_bounded_to_supported_k() {
        let query = b"ACGTTGCATGTCAGTAGGCATCAGTACCGATGCTAGCTAGGCTAACGTTACGATCGATGCA";
        let closed = extract_syncmers(query, SyncmerConfig::k31(SyncmerMode::Closed)).unwrap();
        let open = extract_syncmers(query, SyncmerConfig::k31(SyncmerMode::Open)).unwrap();
        assert!(closed.len() >= open.len());
        assert!(closed.iter().all(|seed| seed.k == 31));
        assert!(open.iter().all(|seed| seed.k == 31));
    }

    #[test]
    fn reverse_complement_preserves_open_selection_count() {
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
        let left = extract_syncmers(query, SyncmerConfig::k21_rescue(SyncmerMode::Open)).unwrap();
        let right =
            extract_syncmers(&reverse, SyncmerConfig::k21_rescue(SyncmerMode::Open)).unwrap();
        assert_eq!(left.len(), right.len());
    }

    #[test]
    fn ambiguity_is_skipped_without_making_a_false_seed() {
        let query = b"ACGTCAGTACCGATGCTAGCTAGGCTAACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        let seeds = extract_syncmers(query, SyncmerConfig::k31(SyncmerMode::Closed)).unwrap();
        assert!(seeds.iter().all(|seed| seed.position + 31 <= 29));
    }

    #[test]
    fn exact_verification_material_is_not_just_the_digest() {
        let seed = extract_syncmers(
            b"ACGTCAGTACCGATGCTAGCTAGGCTAACGTTACGATCGATGCA",
            SyncmerConfig::k31(SyncmerMode::Closed),
        )
        .unwrap()
        .into_iter()
        .next()
        .unwrap();
        assert_eq!(seed.verification_bytes().len(), 8);
        assert_ne!(seed.hash, 0);
    }
}
