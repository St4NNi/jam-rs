//! Frozen sensitivity and alignment configuration shared by trace stages.

use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SensitivityProfile {
    Fast,
    Balanced,
    Sensitive,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SeedSensitivity {
    pub k: u8,
    pub scale: u64,
    pub max_occurrences: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct AlignmentScoring {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open_score: i32,
    pub gap_extend_score: i32,
    pub band_width: u32,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SensitivityConfig {
    pub profile: SensitivityProfile,
    pub primary: SeedSensitivity,
    pub rescue: Option<SeedSensitivity>,
    pub max_candidates: u32,
    pub max_anchors_per_candidate: u32,
    pub max_chains_per_candidate: u32,
    pub min_chain_anchors: u32,
    pub max_alignment_window_bases: u64,
    pub max_concurrent_candidates: u32,
    pub alignment: AlignmentScoring,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct TraceAlgorithmParameters {
    pub hash_id: String,
    pub seed_selection: String,
    pub sensitivity: SensitivityConfig,
    pub max_chain_predecessors: u32,
    pub max_query_gap: u64,
    pub max_target_gap: u64,
    pub chain_gap_penalty: i64,
    pub alignment_mode: String,
    pub x_drop: Option<i32>,
    pub band_widening: bool,
    pub coverage_mode: String,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct TraceAlgorithmMetadata {
    pub id: String,
    pub version: u16,
    pub parameters: TraceAlgorithmParameters,
}

impl TraceAlgorithmMetadata {
    #[must_use]
    pub fn for_sensitivity(sensitivity: SensitivityConfig) -> Self {
        let max_chain_predecessors = sensitivity.max_anchors_per_candidate.min(256);
        let max_gap = sensitivity.max_alignment_window_bases;
        Self {
            id: super::TRACE_ALGORITHM_ID.to_string(),
            version: super::TRACE_ALGORITHM_VERSION,
            parameters: TraceAlgorithmParameters {
                hash_id: "jamhash_u64_v1".to_string(),
                seed_selection: "all_canonical_kmers_fracminhash".to_string(),
                sensitivity,
                max_chain_predecessors,
                max_query_gap: max_gap,
                max_target_gap: max_gap,
                chain_gap_penalty: 1,
                alignment_mode: "local_affine_gap_fixed_band".to_string(),
                x_drop: None,
                band_widening: false,
                coverage_mode: "nonredundant_supported_query_union".to_string(),
            },
        }
    }
}

impl SensitivityConfig {
    #[must_use]
    pub fn for_profile(profile: SensitivityProfile) -> Self {
        match profile {
            SensitivityProfile::Fast => Self {
                profile,
                primary: SeedSensitivity {
                    k: 31,
                    scale: 500,
                    max_occurrences: 64,
                },
                rescue: None,
                max_candidates: 25,
                max_anchors_per_candidate: 25_000,
                max_chains_per_candidate: 4,
                min_chain_anchors: 3,
                max_alignment_window_bases: 250_000,
                max_concurrent_candidates: 4,
                alignment: AlignmentScoring {
                    match_score: 2,
                    mismatch_score: -3,
                    gap_open_score: -5,
                    gap_extend_score: -1,
                    band_width: 64,
                },
            },
            SensitivityProfile::Balanced => Self {
                profile,
                primary: SeedSensitivity {
                    k: 31,
                    scale: 200,
                    max_occurrences: 128,
                },
                rescue: Some(SeedSensitivity {
                    k: 21,
                    scale: 500,
                    max_occurrences: 128,
                }),
                max_candidates: 100,
                max_anchors_per_candidate: 100_000,
                max_chains_per_candidate: 8,
                min_chain_anchors: 3,
                max_alignment_window_bases: 1_000_000,
                max_concurrent_candidates: 4,
                alignment: AlignmentScoring {
                    match_score: 2,
                    mismatch_score: -3,
                    gap_open_score: -5,
                    gap_extend_score: -1,
                    band_width: 128,
                },
            },
            SensitivityProfile::Sensitive => Self {
                profile,
                primary: SeedSensitivity {
                    k: 31,
                    scale: 100,
                    max_occurrences: 256,
                },
                rescue: Some(SeedSensitivity {
                    k: 21,
                    scale: 200,
                    max_occurrences: 256,
                }),
                max_candidates: 250,
                max_anchors_per_candidate: 250_000,
                max_chains_per_candidate: 16,
                min_chain_anchors: 2,
                max_alignment_window_bases: 4_000_000,
                max_concurrent_candidates: 2,
                alignment: AlignmentScoring {
                    match_score: 2,
                    mismatch_score: -3,
                    gap_open_score: -5,
                    gap_extend_score: -1,
                    band_width: 256,
                },
            },
        }
    }

    pub fn validate(&self) -> Result<(), SensitivityError> {
        if self.primary.k != 31 {
            return Err(SensitivityError::PrimaryK(self.primary.k));
        }
        if self.rescue.is_some_and(|rescue| rescue.k != 21) {
            return Err(SensitivityError::RescueK(
                self.rescue.map_or(0, |rescue| rescue.k),
            ));
        }
        if self.primary.scale == 0 || self.rescue.is_some_and(|rescue| rescue.scale == 0) {
            return Err(SensitivityError::ZeroScale);
        }
        if self.max_candidates == 0
            || self.max_anchors_per_candidate == 0
            || self.max_chains_per_candidate == 0
            || self.max_alignment_window_bases == 0
            || self.max_concurrent_candidates == 0
        {
            return Err(SensitivityError::ZeroLimit);
        }
        Ok(())
    }
}

impl Default for SensitivityConfig {
    fn default() -> Self {
        Self::for_profile(SensitivityProfile::Balanced)
    }
}

#[derive(Debug, Error)]
pub enum SensitivityError {
    #[error("primary trace seed k must be 31, got {0}")]
    PrimaryK(u8),
    #[error("rescue trace seed k must be 21, got {0}")]
    RescueK(u8),
    #[error("trace seed scale must be greater than zero")]
    ZeroScale,
    #[error("trace resource limits must be greater than zero")]
    ZeroLimit,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shipped_profiles_validate() {
        for profile in [
            SensitivityProfile::Fast,
            SensitivityProfile::Balanced,
            SensitivityProfile::Sensitive,
        ] {
            SensitivityConfig::for_profile(profile).validate().unwrap();
        }
    }

    #[test]
    fn balanced_uses_fixed_primary_and_rescue_k() {
        let config = SensitivityConfig::default();
        assert_eq!(config.primary.k, 31);
        assert_eq!(config.rescue.unwrap().k, 21);
    }

    #[test]
    fn algorithm_metadata_resolves_every_fixed_parameter() {
        let metadata = TraceAlgorithmMetadata::for_sensitivity(SensitivityConfig::default());
        assert_eq!(metadata.id, super::super::TRACE_ALGORITHM_ID);
        assert_eq!(metadata.version, super::super::TRACE_ALGORITHM_VERSION);
        assert_eq!(metadata.parameters.hash_id, "jamhash_u64_v1");
        assert_eq!(
            metadata.parameters.seed_selection,
            "all_canonical_kmers_fracminhash"
        );
        assert_eq!(metadata.parameters.max_chain_predecessors, 256);
        assert_eq!(metadata.parameters.chain_gap_penalty, 1);
        assert_eq!(
            metadata.parameters.alignment_mode,
            "local_affine_gap_fixed_band"
        );
        assert_eq!(metadata.parameters.x_drop, None);
        assert!(!metadata.parameters.band_widening);
        assert_eq!(
            metadata.parameters.coverage_mode,
            "nonredundant_supported_query_union"
        );
    }
}
