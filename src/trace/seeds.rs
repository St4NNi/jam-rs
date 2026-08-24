//! Positional, orientation-aware seeds used by the trace search.
//!
//! Trace seeds intentionally use the same packed k-mer representation and
//! `jamhash_u64_v1` identity as the existing sketch path.  Hash zero is not a
//! valid retained seed: it is skipped both here and by the JMA builder.  The
//! exact packed canonical k-mer is retained alongside the hash so a JMA lookup
//! can verify a seed without treating a hash collision as evidence.

use crate::jamhash_u64_v1;
use crate::jma::SeedQuery;
use crate::trace::config::SeedSensitivity;
use crate::trace::model::BaseInterval;
use needletail::Sequence;
use serde::{Deserialize, Serialize};
use thiserror::Error;

pub mod gear;
pub mod spaced;
pub mod strobemer;
pub mod syncmer;

/// Stable hash identity for trace seeds.
pub const HASH_ID: &str = "jamhash_u64_v1";
/// Primary trace seeds are fixed at k=31 by the frozen workflow contract.
pub const PRIMARY_K: u8 = 31;
/// Rescue trace seeds are fixed at k=21 by the frozen workflow contract.
pub const RESCUE_K: u8 = 21;

/// A positional query seed.  `reverse` is true when the canonical packed
/// k-mer is the reverse complement of the sequence at `position`.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Serialize, Deserialize)]
pub struct QuerySeed {
    pub position: u64,
    pub hash: u64,
    pub canonical_kmer: u64,
    pub reverse: bool,
}

impl QuerySeed {
    #[must_use]
    pub const fn query(self, k: u8) -> SeedQuery {
        SeedQuery {
            k,
            hash: self.hash,
            canonical_kmer: self.canonical_kmer,
        }
    }
}

/// One retained FracMinHash density level for a single k-mer size.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SeedLevel {
    pub k: u8,
    pub scale: u64,
    pub seeds: Vec<QuerySeed>,
    /// Hash-zero windows are counted for diagnostics but never emitted.
    pub skipped_hash_zero: u64,
    /// Valid non-zero windows rejected by the scale threshold.
    pub skipped_by_density: u64,
}

impl SeedLevel {
    #[must_use]
    pub fn len(&self) -> usize {
        self.seeds.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.seeds.is_empty()
    }

    /// The exact JMA queries represented by this level, in positional order.
    pub fn queries(&self) -> impl Iterator<Item = SeedQuery> + '_ {
        self.seeds.iter().copied().map(|seed| seed.query(self.k))
    }
}

/// All requested seed levels for one query sequence.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SeedSketch {
    pub sequence_length: u64,
    pub levels: Vec<SeedLevel>,
}

/// Compatibility-friendly name for callers that call the result a seed set.
pub type SeedSet = SeedSketch;

impl SeedSketch {
    #[must_use]
    pub fn level(&self, k: u8, scale: u64) -> Option<&SeedLevel> {
        self.levels
            .iter()
            .find(|level| level.k == k && level.scale == scale)
    }

    #[must_use]
    pub fn primary(&self) -> Option<&SeedLevel> {
        self.levels.iter().find(|level| level.k == PRIMARY_K)
    }

    #[must_use]
    pub fn rescue(&self) -> Option<&SeedLevel> {
        self.levels.iter().find(|level| level.k == RESCUE_K)
    }
}

/// Errors from trace seed configuration or packed-k-mer extraction.
#[derive(Debug, Error, Eq, PartialEq)]
pub enum SeedError {
    #[error("trace seed k must be 21 or 31, got {0}")]
    UnsupportedK(u8),
    #[error("trace seed scale must be greater than zero")]
    ZeroScale,
    #[error("trace seed levels must use unique k/scale pairs")]
    DuplicateLevel,
    #[error("trace seed level k={k} has scale {scale} below a denser level scale {denser_scale}")]
    NonNestedDensity {
        k: u8,
        scale: u64,
        denser_scale: u64,
    },
    #[error("trace seed interval is outside the query: start={start}, end={end}, length={length}")]
    IntervalOutOfBounds { start: u64, end: u64, length: u64 },
    #[error("trace seed interval coordinates overflowed")]
    CoordinateOverflow,
}

/// Compute the retained-hash threshold used by the existing `.jam` sketch
/// implementation.  Smaller scales retain a denser, nested set of hashes.
#[must_use]
pub const fn retention_threshold(scale: u64) -> u64 {
    u64::MAX / scale
}

/// Return whether a hash belongs to a FracMinHash level.
///
/// Hash zero is deliberately excluded before applying the scale threshold.
#[must_use]
pub const fn retained_hash(hash: u64, scale: u64) -> bool {
    scale != 0 && hash != 0 && hash < retention_threshold(scale)
}

/// Extract one fixed-k trace seed level from a DNA sequence.
pub fn extract_seed_level(
    sequence: &[u8],
    config: SeedSensitivity,
) -> Result<SeedLevel, SeedError> {
    extract_seed_level_in_intervals(
        sequence,
        config,
        &[BaseInterval {
            start: 0,
            end: sequence.len() as u64,
        }],
        0,
    )
}

/// Extract one seed level only from the requested query intervals.
///
/// Intervals use zero-based, half-open query coordinates.  Each interval is
/// expanded by `flank_bases`, clipped to the query, and overlapping expanded
/// intervals are merged before hashing.  Consequently a rescue round can
/// pass unresolved gaps without re-emitting duplicate seeds or changing their
/// global query positions.  A k-mer is emitted only when its complete window
/// lies inside the merged interval.
pub fn extract_seed_level_in_intervals(
    sequence: &[u8],
    config: SeedSensitivity,
    intervals: &[BaseInterval],
    flank_bases: u64,
) -> Result<SeedLevel, SeedError> {
    validate_config(config)?;
    let ranges = expanded_seed_intervals(intervals, sequence.len() as u64, flank_bases)?;
    extract_seed_level_in_ranges(sequence, config, &ranges)
}

fn extract_seed_level_in_ranges(
    sequence: &[u8],
    config: SeedSensitivity,
    ranges: &[BaseInterval],
) -> Result<SeedLevel, SeedError> {
    let mut level = SeedLevel {
        k: config.k,
        scale: config.scale,
        seeds: Vec::new(),
        skipped_hash_zero: 0,
        skipped_by_density: 0,
    };

    // `normalize(false)` makes U and lower-case bases consistent with the
    // existing sketch path.  `bit_kmers` skips windows containing N or any
    // other non-ACGT symbol and reports the position within the sliced range.
    // Hashing each merged range separately keeps gap rescue selective while
    // adding the range start back to every position.
    for range in ranges {
        let start = usize::try_from(range.start).map_err(|_| SeedError::CoordinateOverflow)?;
        let end = usize::try_from(range.end).map_err(|_| SeedError::CoordinateOverflow)?;
        let normalized = sequence[start..end].normalize(false);
        for (position, kmer, reverse) in normalized.bit_kmers(config.k, true) {
            let position = range
                .start
                .checked_add(position as u64)
                .ok_or(SeedError::CoordinateOverflow)?;
            let hash = jamhash_u64_v1(kmer.0);
            if hash == 0 {
                level.skipped_hash_zero += 1;
                continue;
            }
            if hash >= retention_threshold(config.scale) {
                level.skipped_by_density += 1;
                continue;
            }
            level.seeds.push(QuerySeed {
                position,
                hash,
                canonical_kmer: kmer.0,
                reverse,
            });
        }
    }

    Ok(level)
}

/// Extract several nested seed levels from selected query intervals.
pub fn extract_seed_levels_in_intervals(
    sequence: &[u8],
    configs: &[SeedSensitivity],
    intervals: &[BaseInterval],
    flank_bases: u64,
) -> Result<SeedSketch, SeedError> {
    let mut requested = configs.to_vec();
    for config in &requested {
        validate_config(*config)?;
    }
    requested.sort_by_key(|config| (k_order(config.k), config.scale));
    for pair in requested.windows(2) {
        if pair[0].k == pair[1].k && pair[0].scale == pair[1].scale {
            return Err(SeedError::DuplicateLevel);
        }
    }

    let ranges = expanded_seed_intervals(intervals, sequence.len() as u64, flank_bases)?;
    let mut levels = Vec::with_capacity(requested.len());
    for config in requested {
        levels.push(extract_seed_level_in_ranges(sequence, config, &ranges)?);
    }

    for pair in levels.windows(2) {
        if pair[0].k == pair[1].k && pair[0].scale > pair[1].scale {
            return Err(SeedError::NonNestedDensity {
                k: pair[1].k,
                scale: pair[1].scale,
                denser_scale: pair[0].scale,
            });
        }
    }

    Ok(SeedSketch {
        sequence_length: sequence.len() as u64,
        levels,
    })
}

/// Merge and clip gap-rescue intervals after adding checked flanks.
///
/// Empty intervals are ignored.  Invalid or out-of-bounds intervals are
/// rejected instead of silently changing the requested query coordinates.
pub fn expanded_seed_intervals(
    intervals: &[BaseInterval],
    sequence_length: u64,
    flank_bases: u64,
) -> Result<Vec<BaseInterval>, SeedError> {
    let mut expanded = Vec::with_capacity(intervals.len());
    for interval in intervals {
        if interval.start > interval.end || interval.end > sequence_length {
            return Err(SeedError::IntervalOutOfBounds {
                start: interval.start,
                end: interval.end,
                length: sequence_length,
            });
        }
        if interval.start == interval.end {
            continue;
        }
        let start = interval.start.saturating_sub(flank_bases);
        let end = interval
            .end
            .checked_add(flank_bases)
            .ok_or(SeedError::CoordinateOverflow)?
            .min(sequence_length);
        expanded.push(BaseInterval { start, end });
    }
    expanded.sort_by_key(|interval| (interval.start, interval.end));

    let mut merged = Vec::with_capacity(expanded.len());
    for interval in expanded {
        let Some(last) = merged.last_mut() else {
            merged.push(interval);
            continue;
        };
        if interval.start <= last.end {
            last.end = last.end.max(interval.end);
        } else {
            merged.push(interval);
        }
    }
    Ok(merged)
}

/// Extract deterministic nested levels.  Levels are ordered by k=31 before
/// k=21 and, within one k, from denser (smaller scale) to sparser (larger
/// scale).  Consequently every later level is a subset of an earlier level
/// with the same k.
pub fn extract_seed_levels(
    sequence: &[u8],
    configs: &[SeedSensitivity],
) -> Result<SeedSketch, SeedError> {
    let mut requested = configs.to_vec();
    for config in &requested {
        validate_config(*config)?;
    }
    requested.sort_by_key(|config| (k_order(config.k), config.scale));
    for pair in requested.windows(2) {
        if pair[0].k == pair[1].k && pair[0].scale == pair[1].scale {
            return Err(SeedError::DuplicateLevel);
        }
    }

    let mut levels = Vec::with_capacity(requested.len());
    for config in requested {
        levels.push(extract_seed_level(sequence, config)?);
    }

    // A sorted same-k sequence of scales is nested by construction.  Keep the
    // check explicit because this is a format/provenance invariant and should
    // not be silently weakened if level ordering changes later.
    for pair in levels.windows(2) {
        if pair[0].k == pair[1].k && pair[0].scale > pair[1].scale {
            return Err(SeedError::NonNestedDensity {
                k: pair[1].k,
                scale: pair[1].scale,
                denser_scale: pair[0].scale,
            });
        }
    }

    Ok(SeedSketch {
        sequence_length: sequence.len() as u64,
        levels,
    })
}

/// Extract the levels used by a sensitivity configuration.
pub fn extract_for_sensitivity(
    sequence: &[u8],
    primary: SeedSensitivity,
    rescue: Option<SeedSensitivity>,
) -> Result<SeedSketch, SeedError> {
    let mut configs = Vec::with_capacity(2);
    configs.push(primary);
    if let Some(rescue) = rescue {
        configs.push(rescue);
    }
    extract_seed_levels(sequence, &configs)
}

/// Alias emphasizing that the query is a plasmid sequence.
pub fn extract_query_seeds(
    sequence: &[u8],
    config: SeedSensitivity,
) -> Result<SeedLevel, SeedError> {
    extract_seed_level(sequence, config)
}

fn validate_config(config: SeedSensitivity) -> Result<(), SeedError> {
    if config.k != PRIMARY_K && config.k != RESCUE_K {
        return Err(SeedError::UnsupportedK(config.k));
    }
    if config.scale == 0 {
        return Err(SeedError::ZeroScale);
    }
    Ok(())
}

const fn k_order(k: u8) -> u8 {
    match k {
        PRIMARY_K => 0,
        RESCUE_K => 1,
        _ => 2,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hash_zero_is_never_retained() {
        assert!(!retained_hash(0, 1));
        assert!(!retained_hash(0, 500));
        assert!(!retained_hash(1, 0));
    }

    #[test]
    fn extracted_hash_zero_windows_are_counted_but_not_emitted() {
        for (k, window_count) in [(PRIMARY_K, 34u64), (RESCUE_K, 44u64)] {
            let sequence = vec![b'A'; 64];
            let level = extract_seed_level(
                &sequence,
                SeedSensitivity {
                    k,
                    scale: 1,
                    max_occurrences: 1,
                },
            )
            .unwrap();
            assert!(level.seeds.is_empty());
            assert_eq!(level.skipped_hash_zero, window_count);
            assert_eq!(level.skipped_by_density, 0);
        }
    }

    #[test]
    fn levels_are_ordered_dense_to_sparse() {
        let sketch = extract_seed_levels(
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            &[
                SeedSensitivity {
                    k: 31,
                    scale: 500,
                    max_occurrences: 1,
                },
                SeedSensitivity {
                    k: 31,
                    scale: 100,
                    max_occurrences: 1,
                },
            ],
        )
        .unwrap();
        assert_eq!(sketch.levels[0].scale, 100);
        assert_eq!(sketch.levels[1].scale, 500);
        assert!(
            sketch.levels[1]
                .seeds
                .iter()
                .all(|seed| sketch.levels[0].seeds.contains(seed))
        );
    }

    #[test]
    fn extraction_skips_windows_with_ambiguous_bases() {
        let level = extract_seed_level(
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNACGTACGT",
            SeedSensitivity {
                k: 31,
                scale: 1,
                max_occurrences: 1,
            },
        )
        .unwrap();
        assert!(level.seeds.iter().all(|seed| seed.position < 44));
    }

    #[test]
    fn invalid_k_and_scale_are_rejected() {
        let invalid_k = extract_seed_level(
            b"ACGT",
            SeedSensitivity {
                k: 15,
                scale: 1,
                max_occurrences: 1,
            },
        );
        assert_eq!(invalid_k, Err(SeedError::UnsupportedK(15)));
        let invalid_scale = extract_seed_level(
            b"ACGT",
            SeedSensitivity {
                k: 31,
                scale: 0,
                max_occurrences: 1,
            },
        );
        assert_eq!(invalid_scale, Err(SeedError::ZeroScale));
    }

    #[test]
    fn duplicate_levels_are_rejected() {
        let config = SeedSensitivity {
            k: PRIMARY_K,
            scale: 1,
            max_occurrences: 1,
        };
        assert_eq!(
            extract_seed_levels(b"ACGTACGTACGTACGTACGTACGTACGTACGT", &[config, config]),
            Err(SeedError::DuplicateLevel)
        );
    }
}
