//! Parallel, append-only construction of local Jam Index parts.

use super::manifest::{JamIndexManifest, JamIndexPart, ScreenSelectionPolicy};
use super::part::{JamIndexPartReader, MetagenomeSource, write_part};
use crate::provenance;
use crate::reader::JamReader;
use crate::writer::{HashSampleInput, build_hash_samples};
use needletail::parse_fastx_file;
use rayon::prelude::*;
use sha2::{Digest, Sha256};
use std::collections::BTreeSet;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use thiserror::Error;

pub const JAM_INDEX_MANIFEST_FILE: &str = "manifest.json";

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct JamIndexBuildConfig {
    pub selection_policy: ScreenSelectionPolicy,
    pub max_part_bases: u64,
    pub max_part_signatures: u64,
    pub parallel_parts: usize,
    pub source_manifest_sha256: String,
}

impl Default for JamIndexBuildConfig {
    fn default() -> Self {
        Self {
            selection_policy: ScreenSelectionPolicy::default_signatures(),
            max_part_bases: 1_000_000_000,
            max_part_signatures: 1_000_000,
            parallel_parts: 1,
            source_manifest_sha256: "0".repeat(64),
        }
    }
}

impl JamIndexBuildConfig {
    pub fn validate(&self) -> Result<(), JamIndexBuildError> {
        self.selection_policy
            .validate()
            .map_err(|error| JamIndexBuildError::InvalidConfig(error.to_string()))?;
        if self.max_part_bases == 0
            || self.max_part_signatures == 0
            || self.parallel_parts == 0
            || !valid_checksum(&self.source_manifest_sha256)
        {
            return Err(JamIndexBuildError::InvalidConfig(
                "part bounds, parallelism, and source checksum must be valid".to_string(),
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct JamIndexBuildStats {
    pub new_parts: u32,
    pub total_parts: u32,
    pub new_metagenomes: u64,
    pub total_metagenomes: u64,
    pub new_contigs: u64,
    pub total_contigs: u64,
    pub new_bases: u64,
    pub total_bases: u64,
    pub screen_jam_bytes: u64,
    pub contig_posting_bytes: u64,
    pub source_reference_bytes: u64,
}

#[derive(Clone, Debug)]
struct SourceEstimate {
    source: MetagenomeSource,
    contig_count: u64,
    total_bases: u64,
    estimated_signatures: u64,
}

#[derive(Clone, Debug)]
struct PartPlan {
    part_id: u32,
    sources: Vec<SourceEstimate>,
    contig_count: u64,
    total_bases: u64,
    estimated_signatures: u64,
}

pub fn build_jam_index(
    output: impl AsRef<Path>,
    sources: &[MetagenomeSource],
    config: &JamIndexBuildConfig,
) -> Result<JamIndexBuildStats, JamIndexBuildError> {
    config.validate()?;
    let output = output.as_ref();
    fs::create_dir(output)?;
    fs::create_dir(output.join("parts"))?;
    let estimates = inspect_sources(sources, &config.selection_policy)?;
    let plans = split_parts(&estimates, config, 0)?;
    let parts = build_parts(output, &plans, config)?;
    let mut manifest = JamIndexManifest::empty(config.selection_policy.clone());
    manifest.append_parts(config.source_manifest_sha256.clone(), parts.clone())?;
    write_manifest_atomic(output, &manifest)?;
    Ok(stats(&parts, &manifest))
}

pub fn append_jam_index(
    output: impl AsRef<Path>,
    sources: &[MetagenomeSource],
    config: &JamIndexBuildConfig,
) -> Result<JamIndexBuildStats, JamIndexBuildError> {
    config.validate()?;
    let output = output.as_ref();
    let mut manifest = load_manifest(output)?;
    if manifest.selection_policy != config.selection_policy {
        return Err(JamIndexBuildError::PolicyMismatch);
    }
    reject_existing_metagenomes(output, &manifest, sources)?;
    let estimates = inspect_sources(sources, &config.selection_policy)?;
    let plans = split_parts(&estimates, config, manifest.next_part_id())?;
    let parts = build_parts(output, &plans, config)?;
    let chained_source_checksum = chained_checksum(
        &manifest.source_manifest_sha256,
        &config.source_manifest_sha256,
    );
    let old_metagenomes = manifest.total_metagenomes;
    let old_contigs = manifest.total_contigs;
    let old_bases = manifest.total_bases;
    manifest.append_parts(chained_source_checksum, parts.clone())?;
    write_manifest_atomic(output, &manifest)?;
    let mut result = stats(&parts, &manifest);
    result.new_metagenomes = manifest.total_metagenomes.saturating_sub(old_metagenomes);
    result.new_contigs = manifest.total_contigs.saturating_sub(old_contigs);
    result.new_bases = manifest.total_bases.saturating_sub(old_bases);
    Ok(result)
}

pub fn load_manifest(root: impl AsRef<Path>) -> Result<JamIndexManifest, JamIndexBuildError> {
    let path = root.as_ref().join(JAM_INDEX_MANIFEST_FILE);
    let manifest: JamIndexManifest =
        serde_json::from_reader(std::io::BufReader::new(fs::File::open(&path)?))?;
    manifest.validate()?;
    Ok(manifest)
}

fn inspect_sources(
    sources: &[MetagenomeSource],
    policy: &ScreenSelectionPolicy,
) -> Result<Vec<SourceEstimate>, JamIndexBuildError> {
    if sources.is_empty() {
        return Err(JamIndexBuildError::InvalidConfig(
            "at least one metagenome is required".to_string(),
        ));
    }
    let mut ids = BTreeSet::new();
    sources
        .iter()
        .map(|source| {
            if source.metagenome_id.trim().is_empty()
                || !source.sequence_path.is_file()
                || !ids.insert(source.metagenome_id.clone())
            {
                return Err(JamIndexBuildError::InvalidSource(
                    source.metagenome_id.clone(),
                ));
            }
            let mut reader = parse_fastx_file(&source.sequence_path).map_err(|error| {
                JamIndexBuildError::Parse {
                    path: source.sequence_path.clone(),
                    message: error.to_string(),
                }
            })?;
            let mut lengths = Vec::new();
            while let Some(record) = reader.next() {
                let record = record.map_err(|error| JamIndexBuildError::Parse {
                    path: source.sequence_path.clone(),
                    message: error.to_string(),
                })?;
                lengths.push(
                    u64::try_from(record.seq().len()).map_err(|_| JamIndexBuildError::Overflow)?,
                );
            }
            if lengths.is_empty() {
                return Err(JamIndexBuildError::InvalidSource(
                    source.metagenome_id.clone(),
                ));
            }
            Ok(SourceEstimate {
                source: source.clone(),
                contig_count: u64::try_from(lengths.len())
                    .map_err(|_| JamIndexBuildError::Overflow)?,
                total_bases: lengths.iter().copied().fold(0u64, u64::saturating_add),
                estimated_signatures: policy.estimated_signature_count(&lengths),
            })
        })
        .collect()
}

fn split_parts(
    sources: &[SourceEstimate],
    config: &JamIndexBuildConfig,
    first_part_id: u32,
) -> Result<Vec<PartPlan>, JamIndexBuildError> {
    let mut plans = Vec::new();
    let mut current = PartPlan {
        part_id: first_part_id,
        sources: Vec::new(),
        contig_count: 0,
        total_bases: 0,
        estimated_signatures: 0,
    };
    for source in sources {
        let exceeds = !current.sources.is_empty()
            && (current.total_bases.saturating_add(source.total_bases) > config.max_part_bases
                || current
                    .estimated_signatures
                    .saturating_add(source.estimated_signatures)
                    > config.max_part_signatures);
        if exceeds {
            plans.push(current);
            current = PartPlan {
                part_id: first_part_id
                    .checked_add(
                        u32::try_from(plans.len()).map_err(|_| JamIndexBuildError::Overflow)?,
                    )
                    .ok_or(JamIndexBuildError::Overflow)?,
                sources: Vec::new(),
                contig_count: 0,
                total_bases: 0,
                estimated_signatures: 0,
            };
        }
        current.contig_count = current.contig_count.saturating_add(source.contig_count);
        current.total_bases = current.total_bases.saturating_add(source.total_bases);
        current.estimated_signatures = current
            .estimated_signatures
            .saturating_add(source.estimated_signatures);
        current.sources.push(source.clone());
    }
    if !current.sources.is_empty() {
        plans.push(current);
    }
    Ok(plans)
}

fn build_parts(
    output: &Path,
    plans: &[PartPlan],
    config: &JamIndexBuildConfig,
) -> Result<Vec<JamIndexPart>, JamIndexBuildError> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.parallel_parts)
        .build()
        .map_err(|error| JamIndexBuildError::InvalidConfig(error.to_string()))?;
    let mut parts = pool.install(|| {
        plans
            .par_iter()
            .map(|plan| build_part(output, plan, &config.selection_policy))
            .collect::<Result<Vec<_>, JamIndexBuildError>>()
    })?;
    parts.sort_by_key(|part| part.part_id);
    Ok(parts)
}

fn build_part(
    root: &Path,
    plan: &PartPlan,
    policy: &ScreenSelectionPolicy,
) -> Result<JamIndexPart, JamIndexBuildError> {
    let parts_root = root.join("parts");
    let staging = tempfile::Builder::new()
        .prefix(&format!(".part-{:06}-", plan.part_id))
        .tempdir_in(&parts_root)?;
    let data_path = staging.path().join("part.bin");
    let screen_path = staging.path().join("screen.jam");
    let sources = plan
        .sources
        .iter()
        .map(|source| source.source.clone())
        .collect::<Vec<_>>();
    let result = write_part(&data_path, &sources, policy)?;
    if result.total_bases != plan.total_bases || result.contig_count != plan.contig_count {
        return Err(JamIndexBuildError::PlanMismatch(plan.part_id));
    }
    let screen_samples = result
        .screen_samples
        .iter()
        .map(|sample| HashSampleInput {
            sample_name: sample.metagenome_id.clone(),
            hashes: sample.hashes.clone(),
        })
        .collect::<Vec<_>>();
    let screen_stats = build_hash_samples(&screen_path, &screen_samples, 21, 1)?;
    let screen = JamReader::open(&screen_path)?;
    let data = JamIndexPartReader::open(&data_path)?;
    if screen.sample_names()
        != data
            .metagenomes()
            .iter()
            .map(|metagenome| metagenome.metagenome_id.as_str())
            .collect::<Vec<_>>()
        || screen_stats.sample_count != result.metagenome_count
        || screen_stats.total_entries != result.posting_count
    {
        return Err(JamIndexBuildError::PlanMismatch(plan.part_id));
    }
    let final_name = format!("part-{:06}", plan.part_id);
    let final_path = parts_root.join(&final_name);
    if final_path.exists() {
        return Err(JamIndexBuildError::ExistingPart(plan.part_id));
    }
    let staging_path = staging.keep();
    fs::rename(&staging_path, &final_path)?;
    Ok(JamIndexPart {
        part_id: plan.part_id,
        directory: format!("parts/{final_name}"),
        screen_file: "screen.jam".to_string(),
        data_file: "part.bin".to_string(),
        metagenome_count: result.metagenome_count,
        contig_count: result.contig_count,
        total_bases: result.total_bases,
        estimated_signature_count: result.estimated_signature_count,
        screen_jam_bytes: screen_stats.file_size,
        contig_posting_bytes: result.contig_posting_bytes,
        source_reference_bytes: result.source_reference_bytes,
        screen_sha256: provenance::sha256_file(&final_path.join("screen.jam"))?,
        data_sha256: provenance::sha256_file(&final_path.join("part.bin"))?,
    })
}

fn reject_existing_metagenomes(
    root: &Path,
    manifest: &JamIndexManifest,
    sources: &[MetagenomeSource],
) -> Result<(), JamIndexBuildError> {
    let requested = sources
        .iter()
        .map(|source| source.metagenome_id.as_str())
        .collect::<BTreeSet<_>>();
    if requested.len() != sources.len() {
        return Err(JamIndexBuildError::DuplicateMetagenome);
    }
    for part in &manifest.parts {
        let reader = JamIndexPartReader::open(root.join(&part.directory).join(&part.data_file))?;
        if reader
            .metagenomes()
            .iter()
            .any(|metagenome| requested.contains(metagenome.metagenome_id.as_str()))
        {
            return Err(JamIndexBuildError::DuplicateMetagenome);
        }
    }
    Ok(())
}

fn write_manifest_atomic(
    root: &Path,
    manifest: &JamIndexManifest,
) -> Result<(), JamIndexBuildError> {
    manifest.validate()?;
    let mut temporary = tempfile::NamedTempFile::new_in(root)?;
    serde_json::to_writer_pretty(&mut temporary, manifest)?;
    temporary.write_all(b"\n")?;
    temporary.as_file_mut().sync_all()?;
    temporary
        .persist(root.join(JAM_INDEX_MANIFEST_FILE))
        .map_err(|error| JamIndexBuildError::Io(error.error))?;
    Ok(())
}

fn stats(parts: &[JamIndexPart], manifest: &JamIndexManifest) -> JamIndexBuildStats {
    JamIndexBuildStats {
        new_parts: u32::try_from(parts.len()).unwrap_or(u32::MAX),
        total_parts: u32::try_from(manifest.parts.len()).unwrap_or(u32::MAX),
        new_metagenomes: parts
            .iter()
            .map(|part| u64::from(part.metagenome_count))
            .sum(),
        total_metagenomes: manifest.total_metagenomes,
        new_contigs: parts.iter().map(|part| part.contig_count).sum(),
        total_contigs: manifest.total_contigs,
        new_bases: parts.iter().map(|part| part.total_bases).sum(),
        total_bases: manifest.total_bases,
        screen_jam_bytes: parts.iter().map(|part| part.screen_jam_bytes).sum(),
        contig_posting_bytes: parts.iter().map(|part| part.contig_posting_bytes).sum(),
        source_reference_bytes: parts.iter().map(|part| part.source_reference_bytes).sum(),
    }
}

fn chained_checksum(previous: &str, next: &str) -> String {
    let mut digest = Sha256::new();
    digest.update(b"jam-index-source-manifest-chain-v1");
    digest.update(previous.as_bytes());
    digest.update(next.as_bytes());
    format!("{:x}", digest.finalize())
}

fn valid_checksum(value: &str) -> bool {
    value.len() == 64 && value.bytes().all(|byte| byte.is_ascii_hexdigit())
}

#[derive(Debug, Error)]
pub enum JamIndexBuildError {
    #[error("Jam Index builder I/O failed: {0}")]
    Io(#[from] std::io::Error),
    #[error("Jam Index builder JSON failed: {0}")]
    Json(#[from] serde_json::Error),
    #[error("Jam Index manifest failed: {0}")]
    Manifest(#[from] super::manifest::JamIndexManifestError),
    #[error("Jam Index part failed: {0}")]
    Part(#[from] super::part::JamIndexPartError),
    #[error("Jam Index screen writer failed: {0}")]
    ScreenWrite(#[from] crate::writer::CompactError),
    #[error("Jam Index screen reader failed: {0}")]
    ScreenRead(#[from] crate::reader::ReaderError),
    #[error("Jam Index provenance failed: {0}")]
    Provenance(#[from] anyhow::Error),
    #[error("invalid Jam Index build config: {0}")]
    InvalidConfig(String),
    #[error("invalid Jam Index source {0}")]
    InvalidSource(String),
    #[error("Jam Index source parse failed for {path}: {message}")]
    Parse { path: PathBuf, message: String },
    #[error("Jam Index append selection policy does not match the existing dataset")]
    PolicyMismatch,
    #[error("Jam Index metagenome already exists or is duplicated")]
    DuplicateMetagenome,
    #[error("Jam Index part {0} already exists")]
    ExistingPart(u32),
    #[error("Jam Index part {0} does not match its split plan")]
    PlanMismatch(u32),
    #[error("Jam Index build coordinate overflow")]
    Overflow,
}
