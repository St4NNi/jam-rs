//! Existing `.jam` candidate search for the trace workflow.
//!
//! This stage intentionally reports sketch evidence only.  A candidate is a
//! metagenome whose sampled hashes overlap the queried plasmid; it is not a
//! plasmid-origin or presence call.  Uniform and catalog-specific bias scores
//! remain distinguishable, and the uniform E-value is omitted in bias mode.

use crate::core_utils::passes_entropy_filter;
use crate::jamhash_u64_v1;
use crate::query::{QueryEngine, QueryError};
use crate::reader::ReaderError;
use crate::trace::model::CandidateResult;
use needletail::Sequence;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::path::{Path, PathBuf};
use thiserror::Error;

/// Independent candidate filters.  Containment values use explicit
/// denominators: query/plasmid hashes and indexed metagenome hashes.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CandidateSearchConfig {
    pub min_shared_hashes: u32,
    pub min_plasmid_containment: f64,
    pub min_metagenome_containment: f64,
    pub top_candidates: usize,
}

impl Default for CandidateSearchConfig {
    fn default() -> Self {
        Self {
            min_shared_hashes: 3,
            min_plasmid_containment: 0.0,
            min_metagenome_containment: 0.0,
            top_candidates: 100,
        }
    }
}

impl CandidateSearchConfig {
    pub fn validate(&self) -> Result<(), CandidateError> {
        if self.top_candidates == 0 {
            return Err(CandidateError::InvalidConfig(
                "top_candidates must be greater than zero".to_string(),
            ));
        }
        if self.min_shared_hashes == 0 {
            return Err(CandidateError::InvalidConfig(
                "min_shared_hashes must be greater than zero".to_string(),
            ));
        }
        for (name, value) in [
            ("min_plasmid_containment", self.min_plasmid_containment),
            (
                "min_metagenome_containment",
                self.min_metagenome_containment,
            ),
        ] {
            if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                return Err(CandidateError::InvalidConfig(format!(
                    "{name} must be finite and between 0 and 1"
                )));
            }
        }
        Ok(())
    }
}

/// The score identity used to rank a candidate.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CandidateScoreMode {
    Uniform,
    Bias,
    Witness,
}

impl CandidateScoreMode {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Uniform => "uniform",
            Self::Bias => "bias",
            Self::Witness => "witness",
        }
    }
}

/// A ranked candidate with both raw containment directions and optional bias
/// evidence.  `candidate` always uses hit-count denominators, even when the
/// search was bias assisted.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RankedCandidate {
    pub candidate: CandidateResult,
    pub score_mode: CandidateScoreMode,
    pub bias_weighted_plasmid_containment: Option<f64>,
    pub uniform_hash_e_value: Option<f64>,
}

impl RankedCandidate {
    #[must_use]
    pub fn metagenome_id(&self) -> &str {
        &self.candidate.metagenome_id
    }
}

/// Result of one plasmid query against the existing `.jam` index.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CandidateSearchResult {
    pub query_hashes: u64,
    pub hashes_found: u64,
    pub candidates: Vec<RankedCandidate>,
    pub score_mode: CandidateScoreMode,
}

/// Candidate search engine backed by an existing memory-mapped `.jam` file.
pub struct CandidateSearcher {
    engine: QueryEngine,
    database: PathBuf,
}

/// Short alias for callers that refer to this stage as candidate search.
pub type CandidateSearch = CandidateSearcher;

impl std::fmt::Debug for CandidateSearcher {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("CandidateSearcher")
            .field("database", &self.database)
            .field("kmer_size", &self.engine.kmer_size())
            .field("threshold", &self.engine.threshold())
            .field("bias", &self.engine.has_bias_table())
            .finish()
    }
}

impl CandidateSearcher {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, CandidateError> {
        let database = path.as_ref().to_path_buf();
        let engine = QueryEngine::open(&database)?;
        Ok(Self { engine, database })
    }

    #[must_use]
    pub fn database(&self) -> &Path {
        &self.database
    }

    #[must_use]
    pub fn kmer_size(&self) -> u8 {
        self.engine.kmer_size()
    }

    #[must_use]
    pub fn score_mode(&self) -> CandidateScoreMode {
        if self.engine.has_bias_table() {
            CandidateScoreMode::Bias
        } else {
            CandidateScoreMode::Uniform
        }
    }

    /// Search one FASTA/FASTQ record from a path.  Multiple records are
    /// rejected so a caller cannot accidentally turn a plasmid catalog into a
    /// combined query.
    pub fn search_path(
        &self,
        path: impl AsRef<Path>,
        config: &CandidateSearchConfig,
    ) -> Result<CandidateSearchResult, CandidateError> {
        let path = path.as_ref();
        let mut reader = needletail::parse_fastx_file(path).map_err(|error| {
            CandidateError::Query(QueryError::Parse {
                path: path.display().to_string(),
                message: error.to_string(),
            })
        })?;
        let Some(record) = reader.next() else {
            return self.search_sequence(&[], config);
        };
        let record = record.map_err(|error| {
            CandidateError::Query(QueryError::Parse {
                path: path.display().to_string(),
                message: error.to_string(),
            })
        })?;
        let sequence = record.normalize(false).to_vec();
        if reader.next().is_some() {
            return Err(CandidateError::MultipleQueryRecords {
                path: path.to_path_buf(),
            });
        }
        self.search_sequence(&sequence, config)
    }

    /// Search one already parsed plasmid sequence.
    pub fn search_sequence(
        &self,
        sequence: &[u8],
        config: &CandidateSearchConfig,
    ) -> Result<CandidateSearchResult, CandidateError> {
        config.validate()?;
        let hashes = self.query_hashes(sequence);
        let query_size = hashes.len();
        let result = self.engine.query(&hashes);
        let score_mode = self.score_mode();
        let mut candidates = result
            .matches
            .into_iter()
            .filter_map(|matched| {
                let reference_hashes = self.engine.reader().sample_size(matched.sample_id)?;
                let query_containment = if query_size == 0 {
                    0.0
                } else {
                    f64::from(matched.hit_count) / query_size as f64
                };
                let reference_containment = if reference_hashes == 0 {
                    0.0
                } else {
                    f64::from(matched.hit_count) / reference_hashes as f64
                };
                if matched.hit_count < config.min_shared_hashes
                    || query_containment < config.min_plasmid_containment
                    || reference_containment < config.min_metagenome_containment
                {
                    return None;
                }
                let metagenome_id = self
                    .engine
                    .reader()
                    .sample_name(matched.sample_id)
                    .unwrap_or("unknown")
                    .to_string();
                Some(RankedCandidate {
                    candidate: CandidateResult {
                        metagenome_id,
                        shared_hashes: u64::from(matched.hit_count),
                        plasmid_hashes: query_size as u64,
                        metagenome_hashes: reference_hashes,
                        plasmid_containment: query_containment,
                        metagenome_containment: reference_containment,
                        rank: 0,
                        score_mode: score_mode.as_str().to_string(),
                        bias_weighted_plasmid_containment: (score_mode == CandidateScoreMode::Bias)
                            .then_some(matched.containment),
                        uniform_hash_e_value: (score_mode == CandidateScoreMode::Uniform)
                            .then_some(matched.e_value),
                    },
                    score_mode,
                    bias_weighted_plasmid_containment: (score_mode == CandidateScoreMode::Bias)
                        .then_some(matched.containment),
                    uniform_hash_e_value: (score_mode == CandidateScoreMode::Uniform)
                        .then_some(matched.e_value),
                })
            })
            .collect::<Vec<_>>();

        candidates.sort_by(|left, right| {
            let left_score = left
                .bias_weighted_plasmid_containment
                .unwrap_or(left.candidate.plasmid_containment);
            let right_score = right
                .bias_weighted_plasmid_containment
                .unwrap_or(right.candidate.plasmid_containment);
            right_score
                .total_cmp(&left_score)
                .then_with(|| {
                    right
                        .candidate
                        .shared_hashes
                        .cmp(&left.candidate.shared_hashes)
                })
                .then_with(|| {
                    right
                        .candidate
                        .metagenome_containment
                        .total_cmp(&left.candidate.metagenome_containment)
                })
                .then_with(|| {
                    left.candidate
                        .metagenome_id
                        .cmp(&right.candidate.metagenome_id)
                })
        });
        candidates.truncate(config.top_candidates);
        for (rank, candidate) in candidates.iter_mut().enumerate() {
            candidate.candidate.rank = u32::try_from(rank + 1).unwrap_or(u32::MAX);
        }
        Ok(CandidateSearchResult {
            query_hashes: query_size as u64,
            hashes_found: result.hashes_found as u64,
            candidates,
            score_mode,
        })
    }

    fn query_hashes(&self, sequence: &[u8]) -> Vec<u64> {
        let k = self.engine.kmer_size();
        let threshold = self
            .engine
            .bias_table()
            .as_ref()
            .filter(|table| table.is_soft_filter())
            .map_or(self.engine.threshold(), |table| {
                u64::MAX / table.min_fscale()
            });
        let min_entropy = self.engine.reader().min_entropy();
        let bias = self.engine.bias_table();
        let mut hashes = HashSet::new();
        for (_, kmer, _) in sequence.bit_kmers(k, true) {
            let hash = jamhash_u64_v1(kmer.0);
            // Hash zero is excluded by both the existing database builder and
            // query path.  Keeping it out here prevents artificial matches.
            if hash == 0 || hash >= threshold {
                continue;
            }
            if min_entropy > 0.0 && !passes_entropy_filter(kmer.0, k, min_entropy) {
                continue;
            }
            if bias
                .as_ref()
                .is_some_and(|table| !table.passes_filter(hash))
            {
                continue;
            }
            hashes.insert(hash);
        }
        hashes.into_iter().collect()
    }
}

#[derive(Debug, Error)]
pub enum CandidateError {
    #[error("candidate query failed: {0}")]
    Query(#[from] QueryError),
    #[error("candidate database failed: {0}")]
    Database(#[from] ReaderError),
    #[error("invalid candidate search configuration: {0}")]
    InvalidConfig(String),
    #[error("plasmid query contains more than one FASTA/FASTQ record: {path}")]
    MultipleQueryRecords { path: PathBuf },
}
