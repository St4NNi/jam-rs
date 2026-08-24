//! Positional seed construction for JMA v1 archives.
//!
//! JMA uses the released `jamhash_u64_v1` identity and the same FracMinHash
//! threshold convention as the existing `.jam` path. Hash zero is always
//! excluded: `jamhash_u64_v1(0) == 0`, and retaining it would make every
//! zero-valued k-mer an artificial common seed.

use crate::core_utils::passes_entropy_filter;
use crate::jamhash_u64_v1;
use crate::jma::writer::{SeedLevelData, SeedRecord, SeedSection};
use crate::jma::{ContigId, JmaError, JmaResult, SeedLevel, SeedOccurrence, SeedQuery};
use needletail::Sequence;
use std::collections::BTreeSet;

/// Default primary and rescue k-mer identities used by the local builder.
/// The archive itself accepts an open list of scheme descriptors.
pub const PRIMARY_K: u8 = 31;
pub const RESCUE_K: u8 = 21;

/// Input sequence for positional seed construction.
#[derive(Clone, Copy, Debug)]
pub struct SeedInput<'a> {
    pub contig_id: ContigId,
    pub sequence: &'a [u8],
}

/// Seed construction parameters.  `scale` has the same meaning as the
/// existing `.jam` `fscale`: retained hashes satisfy
/// `hash < u64::MAX / scale`.
#[derive(Clone, Copy, Debug)]
pub struct SeedBuildConfig {
    pub k31_scale: u64,
    pub k21_scale: Option<u64>,
    pub min_entropy: Option<f64>,
}

impl Default for SeedBuildConfig {
    fn default() -> Self {
        Self {
            k31_scale: 200,
            k21_scale: Some(500),
            min_entropy: None,
        }
    }
}

/// Unique retained hashes for one seed level.  This is the embedded sample
/// sketch evidence used by candidate retrieval; positional records remain in
/// the corresponding `SeedSection`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EmbeddedSampleSketch {
    pub k: u8,
    pub scale: u64,
    pub hashes: Vec<u64>,
}

/// Result of building the positional sections and their sample-sketch view.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SeedBuildResult {
    pub sections: Vec<SeedSection>,
    pub sample_sketches: Vec<EmbeddedSampleSketch>,
}

/// Builds k=31 primary and optional k=21 rescue sections from contig
/// sequences.  Input order is not part of the encoded result: records are
/// sorted by the JMA writer, while this function assigns no new contig IDs.
pub fn build_seed_sections(
    inputs: &[SeedInput<'_>],
    config: SeedBuildConfig,
) -> JmaResult<SeedBuildResult> {
    validate_config(config)?;

    let mut sections = Vec::with_capacity(2);
    let mut sample_sketches = Vec::with_capacity(2);
    let primary = build_seed_section(inputs, PRIMARY_K, config.k31_scale, config.min_entropy)?;
    sample_sketches.push(primary.sample_sketch.clone());
    sections.push(primary.section);

    if let Some(scale) = config.k21_scale {
        let rescue = build_seed_section(inputs, RESCUE_K, scale, config.min_entropy)?;
        sample_sketches.push(rescue.sample_sketch.clone());
        sections.push(rescue.section);
    }

    Ok(SeedBuildResult {
        sections,
        sample_sketches,
    })
}

/// Builds one positional seed section and its unique retained-hash sketch.
pub fn build_seed_section(
    inputs: &[SeedInput<'_>],
    k: u8,
    scale: u64,
    min_entropy: Option<f64>,
) -> JmaResult<SeedBuildSectionResult> {
    if !(21..=32).contains(&k) {
        return Err(JmaError::CorruptSection(format!(
            "seed k must be between 21 and 32, got k={k}"
        )));
    }
    if scale == 0 {
        return Err(JmaError::CorruptSection(
            "seed scale must be greater than zero".to_string(),
        ));
    }
    if let Some(entropy) = min_entropy
        && (!entropy.is_finite() || !(0.0..=2.0).contains(&entropy))
    {
        return Err(JmaError::CorruptSection(format!(
            "minimum entropy {entropy} is invalid for k={k}"
        )));
    }

    let threshold = u64::MAX / scale;
    let mut records = Vec::new();
    let mut hashes = BTreeSet::new();
    for input in inputs {
        for (position, kmer, reverse) in input.sequence.bit_kmers(k, true) {
            let hash = jamhash_u64_v1(kmer.0);
            if hash == 0 || hash >= threshold {
                continue;
            }
            if let Some(min_entropy) = min_entropy
                && !passes_entropy_filter(kmer.0, k, min_entropy)
            {
                continue;
            }
            hashes.insert(hash);
            records.push(SeedRecord {
                query: SeedQuery {
                    k,
                    hash,
                    canonical_kmer: kmer.0,
                },
                occurrence: SeedOccurrence {
                    contig_id: input.contig_id,
                    position: u64::try_from(position).map_err(|_| JmaError::OffsetOverflow)?,
                    reverse,
                },
            });
        }
    }

    Ok(SeedBuildSectionResult {
        section: SeedSection {
            k,
            levels: vec![SeedLevelData {
                level: SeedLevel { k, scale },
                records,
            }],
        },
        sample_sketch: EmbeddedSampleSketch {
            k,
            scale,
            hashes: hashes.into_iter().collect(),
        },
    })
}

/// One section plus the unique retained-hash sketch associated with it.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SeedBuildSectionResult {
    pub section: SeedSection,
    pub sample_sketch: EmbeddedSampleSketch,
}

fn validate_config(config: SeedBuildConfig) -> JmaResult<()> {
    if config.k31_scale == 0 {
        return Err(JmaError::CorruptSection(
            "k=31 seed scale must be greater than zero".to_string(),
        ));
    }
    if config.k21_scale == Some(0) {
        return Err(JmaError::CorruptSection(
            "k=21 seed scale must be greater than zero".to_string(),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sequence() -> Vec<u8> {
        // A deterministic non-periodic sequence gives enough distinct kmers
        // for both levels without relying on a biological fixture.
        b"ACGTTGCATGTCAGTAGGCATCAGTACCGATGCTAGCTAGGCTAACGTTACGATCGATGCA".to_vec()
    }

    #[test]
    fn builds_both_levels_with_positions_and_sketches() {
        let sequence = sequence();
        let result = build_seed_sections(
            &[SeedInput {
                contig_id: 7,
                sequence: &sequence,
            }],
            SeedBuildConfig {
                k31_scale: 1,
                k21_scale: Some(1),
                min_entropy: None,
            },
        )
        .unwrap();

        assert_eq!(result.sections.len(), 2);
        assert_eq!(result.sections[0].k, 31);
        assert_eq!(result.sections[1].k, 21);
        assert!(!result.sections[0].levels[0].records.is_empty());
        assert!(!result.sections[1].levels[0].records.is_empty());
        assert_eq!(result.sample_sketches[0].k, 31);
        assert_eq!(
            result.sample_sketches[0].hashes.len(),
            result.sample_sketches[0]
                .hashes
                .iter()
                .collect::<BTreeSet<_>>()
                .len()
        );
        assert!(
            result.sections[0].levels[0]
                .records
                .iter()
                .all(|record| record.occurrence.contig_id == 7)
        );
    }

    #[test]
    fn zero_hash_is_never_retained() {
        // k=31 all-A has packed value zero and therefore jamhash zero.
        let sequence = vec![b'A'; 31];
        let result = build_seed_section(
            &[SeedInput {
                contig_id: 0,
                sequence: &sequence,
            }],
            PRIMARY_K,
            1,
            None,
        )
        .unwrap();
        assert!(result.section.levels[0].records.is_empty());
        assert!(result.sample_sketch.hashes.is_empty());
    }

    #[test]
    fn invalid_scales_and_kmers_are_rejected() {
        let sequence = sequence();
        assert!(
            build_seed_sections(
                &[SeedInput {
                    contig_id: 0,
                    sequence: &sequence,
                }],
                SeedBuildConfig {
                    k31_scale: 0,
                    ..SeedBuildConfig::default()
                },
            )
            .is_err()
        );
        assert!(build_seed_section(&[], 19, 1, None).is_err());
    }
}
