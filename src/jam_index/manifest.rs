//! Checked manifest contracts for one logical Jam Index dataset.

use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use std::path::{Component, Path};
use thiserror::Error;

pub const JAM_INDEX_MANIFEST_SCHEMA: &str = "jam-index-manifest-v2";
pub const JAM_INDEX_FORMAT_VERSION: u16 = 2;
pub const BASELINE_SIGNATURE_POLICY_ID: &str = "contig-minhash-v1";
pub const SPATIAL_256_SIGNATURE_POLICY_ID: &str = "contig-spatial-min-256-v1";
pub const SPATIAL_256_TWO_SIGNATURE_POLICY_ID: &str = "contig-spatial-min-256x2-v1";
pub const SPATIAL_256_ADAPTIVE_512_POLICY_ID: &str = "contig-spatial-min-256-adaptive-512-v1";
pub const SPATIAL_256_ADAPTIVE_768_POLICY_ID: &str = "contig-spatial-min-256-adaptive-768-v1";
pub const SPATIAL_256_ADAPTIVE_1024_POLICY_ID: &str = "contig-spatial-min-256-adaptive-1024-v1";

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct ContigSignatureBudget {
    pub minimum: u32,
    pub maximum: u32,
    pub bases_per_signature: u64,
}

impl ContigSignatureBudget {
    pub fn validate(self) -> Result<(), JamIndexManifestError> {
        if self.minimum == 0 || self.maximum < self.minimum || self.bases_per_signature == 0 {
            return Err(JamIndexManifestError::InvalidSelectionPolicy);
        }
        Ok(())
    }

    #[must_use]
    pub fn budget_for_bases(self, bases: u64) -> u32 {
        let proportional =
            bases.saturating_add(self.bases_per_signature - 1) / self.bases_per_signature;
        u32::try_from(proportional)
            .unwrap_or(u32::MAX)
            .clamp(self.minimum, self.maximum)
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct ScreenSelectionPolicy {
    pub policy_id: String,
    pub k: u8,
    pub hash_id: String,
    pub zero_excluded: bool,
    pub contig_budget: ContigSignatureBudget,
    pub whole_metagenome_budget: u32,
    pub query_window_bases: u32,
    pub adaptive_second_minimum_bases: Option<u64>,
}

impl ScreenSelectionPolicy {
    #[must_use]
    pub fn default_signatures() -> Self {
        Self {
            policy_id: BASELINE_SIGNATURE_POLICY_ID.to_string(),
            k: 21,
            hash_id: "jamhash_u64_v1".to_string(),
            zero_excluded: true,
            contig_budget: ContigSignatureBudget {
                minimum: 16,
                maximum: 256,
                bases_per_signature: 1_024,
            },
            whole_metagenome_budget: 512,
            query_window_bases: 256,
            adaptive_second_minimum_bases: None,
        }
    }

    #[must_use]
    pub fn spatial_256(whole_metagenome_budget: u32) -> Self {
        Self::spatial(
            SPATIAL_256_SIGNATURE_POLICY_ID,
            256,
            whole_metagenome_budget,
        )
    }

    #[must_use]
    pub fn spatial_256_two(whole_metagenome_budget: u32) -> Self {
        Self::spatial(
            SPATIAL_256_TWO_SIGNATURE_POLICY_ID,
            128,
            whole_metagenome_budget,
        )
    }

    pub fn spatial_256_adaptive(
        threshold: u64,
        whole_metagenome_budget: u32,
    ) -> Result<Self, JamIndexManifestError> {
        let policy_id = match threshold {
            512 => SPATIAL_256_ADAPTIVE_512_POLICY_ID,
            768 => SPATIAL_256_ADAPTIVE_768_POLICY_ID,
            1_024 => SPATIAL_256_ADAPTIVE_1024_POLICY_ID,
            _ => return Err(JamIndexManifestError::InvalidSelectionPolicy),
        };
        let mut policy = Self::spatial(policy_id, 256, whole_metagenome_budget);
        policy.adaptive_second_minimum_bases = Some(threshold);
        Ok(policy)
    }

    fn spatial(policy_id: &str, segment_bases: u64, whole_metagenome_budget: u32) -> Self {
        Self {
            policy_id: policy_id.to_string(),
            k: 21,
            hash_id: "jamhash_u64_v1".to_string(),
            zero_excluded: true,
            contig_budget: ContigSignatureBudget {
                minimum: 1,
                maximum: 262_144,
                bases_per_signature: segment_bases,
            },
            whole_metagenome_budget,
            query_window_bases: 256,
            adaptive_second_minimum_bases: None,
        }
    }

    #[must_use]
    pub fn spatial_segment_bases(&self) -> Option<u32> {
        match self.policy_id.as_str() {
            SPATIAL_256_SIGNATURE_POLICY_ID
            | SPATIAL_256_TWO_SIGNATURE_POLICY_ID
            | SPATIAL_256_ADAPTIVE_512_POLICY_ID
            | SPATIAL_256_ADAPTIVE_768_POLICY_ID
            | SPATIAL_256_ADAPTIVE_1024_POLICY_ID => Some(256),
            _ => None,
        }
    }

    #[must_use]
    pub fn spatial_signatures_per_segment(&self, contig_bases: u64) -> Option<u32> {
        if self
            .adaptive_second_minimum_bases
            .is_some_and(|threshold| contig_bases >= threshold)
        {
            return Some(2);
        }
        match self.policy_id.as_str() {
            SPATIAL_256_SIGNATURE_POLICY_ID => Some(1),
            SPATIAL_256_TWO_SIGNATURE_POLICY_ID => Some(2),
            SPATIAL_256_ADAPTIVE_512_POLICY_ID
            | SPATIAL_256_ADAPTIVE_768_POLICY_ID
            | SPATIAL_256_ADAPTIVE_1024_POLICY_ID => Some(1),
            _ => None,
        }
    }

    #[must_use]
    pub fn contig_signature_budget(&self, bases: u64) -> u32 {
        if let Some(per_segment) = self.spatial_signatures_per_segment(bases) {
            let segments = bases.saturating_add(255) / 256;
            u32::try_from(segments.saturating_mul(u64::from(per_segment)))
                .unwrap_or(u32::MAX)
                .clamp(self.contig_budget.minimum, self.contig_budget.maximum)
        } else {
            self.contig_budget.budget_for_bases(bases)
        }
    }

    #[must_use]
    pub fn smaller_signatures() -> Self {
        Self {
            policy_id: "contig-minhash-small-v1".to_string(),
            k: 21,
            hash_id: "jamhash_u64_v1".to_string(),
            zero_excluded: true,
            contig_budget: ContigSignatureBudget {
                minimum: 4,
                maximum: 64,
                bases_per_signature: 4_096,
            },
            whole_metagenome_budget: 256,
            query_window_bases: 256,
            adaptive_second_minimum_bases: None,
        }
    }

    pub fn validate(&self) -> Result<(), JamIndexManifestError> {
        self.contig_budget.validate()?;
        if self.policy_id.trim().is_empty()
            || self.k != 21
            || self.hash_id != "jamhash_u64_v1"
            || !self.zero_excluded
            || self.whole_metagenome_budget == 0
            || self.query_window_bases < u32::from(self.k)
        {
            return Err(JamIndexManifestError::InvalidSelectionPolicy);
        }
        if let Some(segment_bases) = self.spatial_segment_bases()
            && (self.contig_budget.minimum != 1
                || self.contig_budget.bases_per_signature
                    != u64::from(segment_bases)
                        / u64::from(self.spatial_signatures_per_segment(0).unwrap_or(1))
                || !matches!(self.whole_metagenome_budget, 512 | 1_024))
        {
            return Err(JamIndexManifestError::InvalidSelectionPolicy);
        }
        let adaptive_id = matches!(
            self.policy_id.as_str(),
            SPATIAL_256_ADAPTIVE_512_POLICY_ID
                | SPATIAL_256_ADAPTIVE_768_POLICY_ID
                | SPATIAL_256_ADAPTIVE_1024_POLICY_ID
        );
        if adaptive_id != self.adaptive_second_minimum_bases.is_some()
            || self.adaptive_second_minimum_bases.is_some_and(|threshold| {
                !matches!(threshold, 512 | 768 | 1_024)
                    || self.policy_id
                        != match threshold {
                            512 => SPATIAL_256_ADAPTIVE_512_POLICY_ID,
                            768 => SPATIAL_256_ADAPTIVE_768_POLICY_ID,
                            _ => SPATIAL_256_ADAPTIVE_1024_POLICY_ID,
                        }
            })
        {
            return Err(JamIndexManifestError::InvalidSelectionPolicy);
        }
        Ok(())
    }

    #[must_use]
    pub fn estimated_signature_count(&self, contig_lengths: &[u64]) -> u64 {
        contig_lengths
            .iter()
            .map(|length| u64::from(self.contig_signature_budget(*length)))
            .fold(u64::from(self.whole_metagenome_budget), u64::saturating_add)
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct JamIndexPart {
    pub part_id: u32,
    pub directory: String,
    pub screen_file: String,
    pub data_file: String,
    pub metagenome_count: u32,
    pub contig_count: u64,
    pub total_bases: u64,
    pub estimated_signature_count: u64,
    pub screen_jam_bytes: u64,
    pub contig_posting_bytes: u64,
    pub source_reference_bytes: u64,
    pub metagenome_directory_bytes: u64,
    pub contig_length_bytes: u64,
    pub exceptional_length_bytes: u64,
    pub string_table_bytes: u64,
    pub screen_sha256: String,
    pub data_sha256: String,
}

impl JamIndexPart {
    pub fn validate(&self) -> Result<(), JamIndexManifestError> {
        if !valid_relative_path(&self.directory)
            || !valid_file_name(&self.screen_file)
            || !self.screen_file.ends_with(".jam")
            || !valid_file_name(&self.data_file)
            || self.metagenome_count == 0
            || self.screen_sha256.len() != 64
            || self.data_sha256.len() != 64
            || !self
                .screen_sha256
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit())
            || !self
                .data_sha256
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit())
        {
            return Err(JamIndexManifestError::InvalidPart(self.part_id));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct JamIndexManifest {
    pub schema_version: String,
    pub format_version: u16,
    pub selection_policy: ScreenSelectionPolicy,
    pub source_manifest_sha256: String,
    pub total_metagenomes: u64,
    pub total_contigs: u64,
    pub total_bases: u64,
    pub estimated_signature_count: u64,
    pub parts: Vec<JamIndexPart>,
}

impl JamIndexManifest {
    #[must_use]
    pub fn empty(selection_policy: ScreenSelectionPolicy) -> Self {
        Self {
            schema_version: JAM_INDEX_MANIFEST_SCHEMA.to_string(),
            format_version: JAM_INDEX_FORMAT_VERSION,
            selection_policy,
            source_manifest_sha256: String::new(),
            total_metagenomes: 0,
            total_contigs: 0,
            total_bases: 0,
            estimated_signature_count: 0,
            parts: Vec::new(),
        }
    }

    #[must_use]
    pub fn next_part_id(&self) -> u32 {
        self.parts
            .last()
            .map_or(0, |part| part.part_id.saturating_add(1))
    }

    pub fn append_parts(
        &mut self,
        source_manifest_sha256: String,
        parts: impl IntoIterator<Item = JamIndexPart>,
    ) -> Result<(), JamIndexManifestError> {
        if source_manifest_sha256.len() != 64
            || !source_manifest_sha256
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit())
        {
            return Err(JamIndexManifestError::InvalidSourceChecksum);
        }
        let expected_first = self.next_part_id();
        let appended = parts.into_iter().collect::<Vec<_>>();
        if appended
            .iter()
            .enumerate()
            .any(|(index, part)| part.part_id != expected_first.saturating_add(index as u32))
        {
            return Err(JamIndexManifestError::NonAppendOnlyParts);
        }
        self.parts.extend(appended);
        self.source_manifest_sha256 = source_manifest_sha256;
        self.recompute_totals();
        self.validate()
    }

    pub fn validate(&self) -> Result<(), JamIndexManifestError> {
        if self.schema_version != JAM_INDEX_MANIFEST_SCHEMA
            || self.format_version != JAM_INDEX_FORMAT_VERSION
        {
            return Err(JamIndexManifestError::UnsupportedManifest);
        }
        self.selection_policy.validate()?;
        if !self.parts.is_empty()
            && (self.source_manifest_sha256.len() != 64
                || !self
                    .source_manifest_sha256
                    .bytes()
                    .all(|byte| byte.is_ascii_hexdigit()))
        {
            return Err(JamIndexManifestError::InvalidSourceChecksum);
        }
        let mut directories = BTreeSet::new();
        for (index, part) in self.parts.iter().enumerate() {
            part.validate()?;
            if part.part_id != index as u32 || !directories.insert(part.directory.clone()) {
                return Err(JamIndexManifestError::NonAppendOnlyParts);
            }
        }
        let expected = totals(&self.parts);
        if expected
            != (
                self.total_metagenomes,
                self.total_contigs,
                self.total_bases,
                self.estimated_signature_count,
            )
        {
            return Err(JamIndexManifestError::TotalMismatch);
        }
        Ok(())
    }

    fn recompute_totals(&mut self) {
        let totals = totals(&self.parts);
        self.total_metagenomes = totals.0;
        self.total_contigs = totals.1;
        self.total_bases = totals.2;
        self.estimated_signature_count = totals.3;
    }
}

fn totals(parts: &[JamIndexPart]) -> (u64, u64, u64, u64) {
    parts.iter().fold((0, 0, 0, 0), |totals, part| {
        (
            totals.0 + u64::from(part.metagenome_count),
            totals.1.saturating_add(part.contig_count),
            totals.2.saturating_add(part.total_bases),
            totals.3.saturating_add(part.estimated_signature_count),
        )
    })
}

fn valid_relative_path(value: &str) -> bool {
    let path = Path::new(value);
    !value.is_empty()
        && !path.is_absolute()
        && path
            .components()
            .all(|component| matches!(component, Component::Normal(_)))
}

fn valid_file_name(value: &str) -> bool {
    valid_relative_path(value) && Path::new(value).components().count() == 1
}

#[derive(Debug, Error, Eq, PartialEq)]
pub enum JamIndexManifestError {
    #[error("unsupported Jam Index manifest schema or format version")]
    UnsupportedManifest,
    #[error("invalid Jam Index screen selection policy")]
    InvalidSelectionPolicy,
    #[error("invalid Jam Index part {0}")]
    InvalidPart(u32),
    #[error("Jam Index parts are not a strictly append-only sequence")]
    NonAppendOnlyParts,
    #[error("invalid Jam Index source manifest checksum")]
    InvalidSourceChecksum,
    #[error("Jam Index manifest totals do not match its parts")]
    TotalMismatch,
}
