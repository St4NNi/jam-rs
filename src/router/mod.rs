//! Shared contracts for the immutable collection witness router.

use crate::jamhash_u64_v1;
use serde::{Deserialize, Serialize};
use thiserror::Error;

pub mod build;
pub mod format;
pub mod postings;
pub mod reader;
pub mod search;
pub mod source;
pub mod witness;
pub mod writer;

pub const ROUTER_FORMAT_NAME: &str = "JAM Witness Router";
pub const ROUTER_FORMAT_VERSION: u16 = 1;
pub const WITNESS_K: u8 = 21;
const K21_PACKED_LIMIT: u64 = 1_u64 << (2 * WITNESS_K);

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum HashAlgorithmId {
    JamhashU64V1,
}

impl HashAlgorithmId {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::JamhashU64V1 => "jamhash_u64_v1",
        }
    }
}

/// One dense k=21 key set whose hash derives every available coarser tier.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct WitnessScheme {
    pub scheme_id: u32,
    pub k: u8,
    pub base_scale: u32,
    pub available_scales: Vec<u32>,
    pub hash_id: HashAlgorithmId,
    pub zero_excluded: bool,
}

impl WitnessScheme {
    pub fn validate(&self) -> Result<(), RouterContractError> {
        if self.k != WITNESS_K {
            return Err(RouterContractError::UnsupportedWitnessK(self.k));
        }
        if self.base_scale == 0 || self.available_scales.is_empty() {
            return Err(RouterContractError::InvalidScales);
        }
        if self.available_scales[0] != self.base_scale
            || self.available_scales.contains(&0)
            || self
                .available_scales
                .windows(2)
                .any(|pair| pair[0] >= pair[1])
        {
            return Err(RouterContractError::InvalidScales);
        }
        if self.hash_id != HashAlgorithmId::JamhashU64V1 || !self.zero_excluded {
            return Err(RouterContractError::UnsupportedHashScheme);
        }
        Ok(())
    }

    pub fn includes_hash(&self, hash: u64, scale: u32) -> Result<bool, RouterContractError> {
        self.validate()?;
        if self.available_scales.binary_search(&scale).is_err() {
            return Err(RouterContractError::UnavailableScale(scale));
        }
        Ok(hash != 0 && hash < u64::MAX / u64::from(scale))
    }

    pub fn retained_scales(&self, hash: u64) -> Result<Vec<u32>, RouterContractError> {
        self.validate()?;
        Ok(self
            .available_scales
            .iter()
            .copied()
            .filter(|scale| hash != 0 && hash < u64::MAX / u64::from(*scale))
            .collect())
    }
}

/// Exact k=21 identity. A digest match without `packed` equality is never evidence.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct WitnessKey {
    pub packed: u64,
    pub jamhash: u64,
}

impl WitnessKey {
    pub fn from_packed(packed: u64) -> Result<Self, RouterContractError> {
        if packed >= K21_PACKED_LIMIT {
            return Err(RouterContractError::InvalidPackedKmer(packed));
        }
        let jamhash = jamhash_u64_v1(packed);
        if jamhash == 0 {
            return Err(RouterContractError::ExcludedZeroHash);
        }
        Ok(Self { packed, jamhash })
    }

    pub fn checked(packed: u64, jamhash: u64) -> Result<Self, RouterContractError> {
        let key = Self::from_packed(packed)?;
        if key.jamhash != jamhash {
            return Err(RouterContractError::HashMismatch);
        }
        Ok(key)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum WitnessClass {
    Rare,
    ModeratelyCommon,
    CollectionCommon,
    Suppressed,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct QueryWitness {
    pub key: WitnessKey,
    pub query_position: u64,
    pub query_reverse: bool,
    pub query_window_ids: Vec<u32>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SharedWitness {
    pub key: WitnessKey,
    pub query_position: u64,
    pub query_reverse: bool,
    pub query_window_id: u32,
    pub document_frequency: u32,
    pub witness_tier: u32,
    pub witness_class: WitnessClass,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct PositionalWitnessOccurrence {
    pub witness: SharedWitness,
    pub contig_id: u32,
    pub position: u64,
    pub reverse: bool,
    pub scheme_id: u32,
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct CandidateWindowEvidence {
    pub matched_query_windows: u32,
    pub total_eligible_query_windows: u32,
    pub longest_supported_window_run: u32,
    pub rare_witness_count: u32,
    pub common_witness_count: u32,
    pub total_shared_witness_count: u32,
}

impl CandidateWindowEvidence {
    #[must_use]
    pub fn support_fraction(&self) -> f64 {
        if self.total_eligible_query_windows == 0 {
            0.0
        } else {
            f64::from(self.matched_query_windows) / f64::from(self.total_eligible_query_windows)
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RoutedCandidate {
    pub metagenome_id: String,
    pub rare_shared_witnesses: u32,
    pub common_shared_witnesses: u32,
    pub supported_query_windows: u32,
    pub weighted_witness_sum: f64,
    pub estimated_query_containment: f64,
    pub candidate_tier: u32,
    pub window_evidence: CandidateWindowEvidence,
    pub shared_witnesses: Vec<SharedWitness>,
    pub positional_witnesses: Vec<PositionalWitnessOccurrence>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum WitnessHandoffMode {
    SampleOnly,
    PositionBearing,
    Hybrid,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct RouterCatalogBinding {
    pub router_checksum: [u8; 32],
    pub catalog_checksum: [u8; 32],
    pub witness_scheme_checksum: [u8; 32],
}

#[derive(Debug, Error, Eq, PartialEq)]
pub enum RouterContractError {
    #[error("witness k must be 21, got {0}")]
    UnsupportedWitnessK(u8),
    #[error("witness scales must start at the base scale and be strictly increasing")]
    InvalidScales,
    #[error("witness scale {0} is not available")]
    UnavailableScale(u32),
    #[error("the witness scheme must use jamhash_u64_v1 and exclude hash zero")]
    UnsupportedHashScheme,
    #[error("packed k=21 value is out of range: {0}")]
    InvalidPackedKmer(u64),
    #[error("jamhash zero is excluded from witness indexes")]
    ExcludedZeroHash,
    #[error("stored jamhash does not match the exact packed kmer")]
    HashMismatch,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scheme() -> WitnessScheme {
        WitnessScheme {
            scheme_id: 1,
            k: WITNESS_K,
            base_scale: 20,
            available_scales: vec![20, 50, 100, 200, 500],
            hash_id: HashAlgorithmId::JamhashU64V1,
            zero_excluded: true,
        }
    }

    #[test]
    fn one_hash_derives_nested_tiers_without_duplicate_keys() {
        let threshold_100 = u64::MAX / 100;
        assert_eq!(
            scheme().retained_scales(threshold_100 - 1).unwrap(),
            vec![20, 50, 100]
        );
    }

    #[test]
    fn exact_key_validation_fails_closed() {
        let key = WitnessKey::from_packed(7).unwrap();
        assert_eq!(WitnessKey::checked(key.packed, key.jamhash).unwrap(), key);
        assert_eq!(
            WitnessKey::checked(key.packed, key.jamhash ^ 1),
            Err(RouterContractError::HashMismatch)
        );
    }
}
