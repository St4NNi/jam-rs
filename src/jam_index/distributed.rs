//! Deterministic, restartable Jam Index construction stages.

use super::builder::{JAM_INDEX_MANIFEST_FILE, JamIndexBuildStats, write_manifest_atomic};
use super::manifest::{JamIndexManifest, JamIndexPart, ScreenSelectionPolicy};
use super::part::{
    JamIndexPartError, JamIndexPartReader, MetagenomeSource, PublishedMetagenomeSource,
    StagedMetagenomeSource, merge_part_fragments, write_part_staged,
};
use crate::provenance;
use crate::reader::{JamReader, ReaderError};
use crate::writer::{CompactError, HashSampleInput, build_hash_samples};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use thiserror::Error;

pub const DISTRIBUTED_PLAN_SCHEMA: &str = "jam-index-build-plan-v1";
pub const FRAGMENT_MANIFEST_SCHEMA: &str = "jam-index-fragment-v1";
pub const MERGED_PART_MANIFEST_SCHEMA: &str = "jam-index-merged-part-v1";

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct IndexPlanSource {
    pub metagenome_id: String,
    pub source_path: PathBuf,
    pub source_size: u64,
    pub source_sha256: String,
    pub estimated_bases: u64,
    pub estimated_signatures: u64,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct IndexBuildFragmentPlan {
    pub fragment_id: u32,
    pub part_id: u32,
    pub estimated_bases: u64,
    pub estimated_signatures: u64,
    pub sources: Vec<IndexPlanSource>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct IndexBuildPartPlan {
    pub part_id: u32,
    pub estimated_bases: u64,
    pub estimated_signatures: u64,
    pub fragments: Vec<IndexBuildFragmentPlan>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct IndexBuildPlan {
    pub schema: String,
    pub source_manifest_sha256: String,
    pub selection_policy: ScreenSelectionPolicy,
    pub target_parts: u32,
    pub fragments_per_part: u32,
    pub estimated_expansion: u64,
    pub parts: Vec<IndexBuildPartPlan>,
}

impl IndexBuildPlan {
    pub fn validate(&self) -> Result<(), DistributedIndexError> {
        self.selection_policy
            .validate()
            .map_err(|error| DistributedIndexError::InvalidPlan(error.to_string()))?;
        if self.schema != DISTRIBUTED_PLAN_SCHEMA
            || !valid_checksum(&self.source_manifest_sha256)
            || self.target_parts == 0
            || self.fragments_per_part == 0
            || self.estimated_expansion == 0
            || self.parts.is_empty()
            || usize::try_from(self.target_parts).ok() != Some(self.parts.len())
        {
            return Err(DistributedIndexError::InvalidPlan(
                "invalid plan header".to_string(),
            ));
        }
        let mut part_ids = BTreeSet::new();
        let mut fragment_ids = BTreeSet::new();
        let mut metagenome_ids = BTreeSet::new();
        for part in &self.parts {
            if !part_ids.insert(part.part_id) || part.fragments.is_empty() {
                return Err(DistributedIndexError::InvalidPlan(
                    "duplicate or empty part".to_string(),
                ));
            }
            let mut part_bases = 0u64;
            let mut part_signatures = 0u64;
            for fragment in &part.fragments {
                if fragment.part_id != part.part_id
                    || !fragment_ids.insert(fragment.fragment_id)
                    || fragment.sources.is_empty()
                {
                    return Err(DistributedIndexError::InvalidPlan(
                        "invalid build fragment".to_string(),
                    ));
                }
                let bases = fragment
                    .sources
                    .iter()
                    .map(|source| source.estimated_bases)
                    .sum::<u64>();
                let signatures = fragment
                    .sources
                    .iter()
                    .map(|source| source.estimated_signatures)
                    .sum::<u64>();
                if bases != fragment.estimated_bases
                    || signatures != fragment.estimated_signatures
                    || fragment.sources.iter().any(|source| {
                        source.metagenome_id.trim().is_empty()
                            || source.source_size == 0
                            || !valid_checksum(&source.source_sha256)
                            || !metagenome_ids.insert(source.metagenome_id.clone())
                    })
                {
                    return Err(DistributedIndexError::InvalidPlan(
                        "invalid fragment sources".to_string(),
                    ));
                }
                part_bases = part_bases.saturating_add(bases);
                part_signatures = part_signatures.saturating_add(signatures);
            }
            if part_bases != part.estimated_bases || part_signatures != part.estimated_signatures {
                return Err(DistributedIndexError::InvalidPlan(
                    "part estimates do not match fragments".to_string(),
                ));
            }
        }
        Ok(())
    }

    pub fn fragment(&self, fragment_id: u32) -> Option<&IndexBuildFragmentPlan> {
        self.parts
            .iter()
            .flat_map(|part| &part.fragments)
            .find(|fragment| fragment.fragment_id == fragment_id)
    }

    pub fn part(&self, part_id: u32) -> Option<&IndexBuildPartPlan> {
        self.parts.iter().find(|part| part.part_id == part_id)
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct IndexFragmentManifest {
    pub schema: String,
    pub fragment_id: u32,
    pub part_id: u32,
    pub source_manifest_sha256: String,
    pub metagenome_ids: Vec<String>,
    pub source_results: Vec<IndexFragmentSourceResult>,
    pub metagenome_count: u32,
    pub contig_count: u64,
    pub total_bases: u64,
    pub estimated_signature_count: u64,
    pub screen_sha256: String,
    pub data_sha256: String,
    pub screen_bytes: u64,
    pub data_bytes: u64,
    pub contig_signature_histogram: BTreeMap<u32, u64>,
    pub single_contig_mappings: u64,
    pub overflow_mappings: u64,
    pub overflow_contigs: u64,
    pub maximum_overflow_count: u32,
    pub signature_run_count: u32,
    pub signature_run_record_limit: u64,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct IndexFragmentSourceResult {
    pub metagenome_id: String,
    pub staged_path: PathBuf,
    pub published_path: PathBuf,
    pub source_size: u64,
    pub source_sha256: String,
    pub contig_count: u32,
    pub total_bases: u64,
    pub validation_status: String,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct MergedPartManifest {
    pub schema: String,
    pub part: JamIndexPart,
    pub fragment_ids: Vec<u32>,
    pub fragment_manifest_sha256: Vec<String>,
    pub contig_signature_histogram: BTreeMap<u32, u64>,
    pub single_contig_mappings: u64,
    pub overflow_mappings: u64,
    pub overflow_contigs: u64,
    pub maximum_overflow_count: u32,
    pub signature_run_count: u32,
    pub signature_run_record_limit: u64,
}

pub fn plan_index(
    mut sources: Vec<IndexPlanSource>,
    source_manifest_sha256: String,
    selection_policy: ScreenSelectionPolicy,
    target_parts: usize,
    fragments_per_part: usize,
    estimated_expansion: u64,
) -> Result<IndexBuildPlan, DistributedIndexError> {
    selection_policy
        .validate()
        .map_err(|error| DistributedIndexError::InvalidPlan(error.to_string()))?;
    if sources.is_empty()
        || target_parts == 0
        || fragments_per_part == 0
        || estimated_expansion == 0
        || !valid_checksum(&source_manifest_sha256)
    {
        return Err(DistributedIndexError::InvalidPlan(
            "planner inputs are invalid".to_string(),
        ));
    }
    let mut ids = BTreeSet::new();
    for source in &mut sources {
        if source.metagenome_id.trim().is_empty()
            || source.source_size == 0
            || !valid_checksum(&source.source_sha256)
            || !ids.insert(source.metagenome_id.clone())
        {
            return Err(DistributedIndexError::InvalidPlan(
                "planner source is invalid".to_string(),
            ));
        }
        source.estimated_bases = source.source_size.saturating_mul(estimated_expansion);
        source.estimated_signatures =
            selection_policy.estimated_signature_count(&[source.estimated_bases]);
    }
    sources.sort_by(|left, right| {
        right
            .estimated_bases
            .cmp(&left.estimated_bases)
            .then_with(|| right.estimated_signatures.cmp(&left.estimated_signatures))
            .then_with(|| left.metagenome_id.cmp(&right.metagenome_id))
    });
    let part_count = target_parts.min(sources.len());
    let mut assigned = vec![Vec::<IndexPlanSource>::new(); part_count];
    let mut part_loads = vec![(0u64, 0u64); part_count];
    for source in sources {
        let part = (0..part_count)
            .min_by_key(|part| (part_loads[*part].0, part_loads[*part].1, *part))
            .ok_or(DistributedIndexError::Overflow)?;
        part_loads[part].0 = part_loads[part].0.saturating_add(source.estimated_bases);
        part_loads[part].1 = part_loads[part]
            .1
            .saturating_add(source.estimated_signatures);
        assigned[part].push(source);
    }
    let mut next_fragment_id = 0u32;
    let mut parts = Vec::with_capacity(part_count);
    for (part_index, mut part_sources) in assigned.into_iter().enumerate() {
        part_sources.sort_by(|left, right| {
            right
                .estimated_bases
                .cmp(&left.estimated_bases)
                .then_with(|| left.metagenome_id.cmp(&right.metagenome_id))
        });
        let fragment_count = fragments_per_part.min(part_sources.len());
        let mut fragment_sources = vec![Vec::<IndexPlanSource>::new(); fragment_count];
        let mut fragment_loads = vec![(0u64, 0u64); fragment_count];
        for source in part_sources {
            let fragment = (0..fragment_count)
                .min_by_key(|fragment| {
                    (
                        fragment_loads[*fragment].0,
                        fragment_loads[*fragment].1,
                        *fragment,
                    )
                })
                .ok_or(DistributedIndexError::Overflow)?;
            fragment_loads[fragment].0 = fragment_loads[fragment]
                .0
                .saturating_add(source.estimated_bases);
            fragment_loads[fragment].1 = fragment_loads[fragment]
                .1
                .saturating_add(source.estimated_signatures);
            fragment_sources[fragment].push(source);
        }
        let part_id = u32::try_from(part_index).map_err(|_| DistributedIndexError::Overflow)?;
        let mut fragments = Vec::with_capacity(fragment_count);
        for (fragment, mut sources) in fragment_sources.into_iter().enumerate() {
            sources.sort_by(|left, right| left.metagenome_id.cmp(&right.metagenome_id));
            let (estimated_bases, estimated_signatures) = fragment_loads[fragment];
            fragments.push(IndexBuildFragmentPlan {
                fragment_id: next_fragment_id,
                part_id,
                estimated_bases,
                estimated_signatures,
                sources,
            });
            next_fragment_id = next_fragment_id
                .checked_add(1)
                .ok_or(DistributedIndexError::Overflow)?;
        }
        parts.push(IndexBuildPartPlan {
            part_id,
            estimated_bases: part_loads[part_index].0,
            estimated_signatures: part_loads[part_index].1,
            fragments,
        });
    }
    let plan = IndexBuildPlan {
        schema: DISTRIBUTED_PLAN_SCHEMA.to_string(),
        source_manifest_sha256,
        selection_policy,
        target_parts: u32::try_from(parts.len()).map_err(|_| DistributedIndexError::Overflow)?,
        fragments_per_part: u32::try_from(fragments_per_part)
            .map_err(|_| DistributedIndexError::Overflow)?,
        estimated_expansion,
        parts,
    };
    plan.validate()?;
    Ok(plan)
}

pub fn load_plan(path: impl AsRef<Path>) -> Result<IndexBuildPlan, DistributedIndexError> {
    let plan: IndexBuildPlan =
        serde_json::from_reader(std::io::BufReader::new(fs::File::open(path)?))?;
    plan.validate()?;
    Ok(plan)
}

pub fn write_plan_atomic(
    path: impl AsRef<Path>,
    plan: &IndexBuildPlan,
) -> Result<(), DistributedIndexError> {
    plan.validate()?;
    write_json_atomic(path.as_ref(), plan)
}

pub fn build_fragment(
    plan: &IndexBuildPlan,
    fragment_id: u32,
    staged_sources: &BTreeMap<String, MetagenomeSource>,
    output: impl AsRef<Path>,
) -> Result<IndexFragmentManifest, DistributedIndexError> {
    plan.validate()?;
    let fragment = plan
        .fragment(fragment_id)
        .ok_or(DistributedIndexError::UnknownFragment(fragment_id))?;
    let sources = fragment
        .sources
        .iter()
        .map(|source| {
            let staged = staged_sources.get(&source.metagenome_id).ok_or_else(|| {
                DistributedIndexError::MissingSource(source.metagenome_id.clone())
            })?;
            if staged.metagenome_id != source.metagenome_id
                || fs::metadata(&staged.sequence_path)?.len() != source.source_size
            {
                return Err(DistributedIndexError::SourceIdentity(
                    source.metagenome_id.clone(),
                ));
            }
            Ok(StagedMetagenomeSource {
                metagenome_id: source.metagenome_id.clone(),
                staged_sequence_path: staged.sequence_path.clone(),
                published_sequence_path: source.source_path.clone(),
                expected_source_size: source.source_size,
                expected_source_sha256: Some(parse_checksum(&source.source_sha256)?),
            })
        })
        .collect::<Result<Vec<_>, DistributedIndexError>>()?;
    let output = output.as_ref();
    if output.exists() {
        return Err(DistributedIndexError::ExistingOutput(output.to_path_buf()));
    }
    let parent = output.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(parent)?;
    let staging = tempfile::Builder::new()
        .prefix(&format!(".fragment-{fragment_id:06}-"))
        .tempdir_in(parent)?;
    let data_path = staging.path().join("part.bin");
    let screen_path = staging.path().join("screen.jam");
    let result = write_part_staged(&data_path, &sources, &plan.selection_policy)?;
    let screen_samples = result
        .screen_samples
        .iter()
        .map(|sample| HashSampleInput {
            sample_name: sample.metagenome_id.clone(),
            hashes: sample.hashes.clone(),
        })
        .collect::<Vec<_>>();
    let screen_stats = build_hash_samples(&screen_path, &screen_samples, 21, 1)?;
    let reader = JamIndexPartReader::open(&data_path)?;
    let mut source_results = Vec::with_capacity(fragment.sources.len());
    for planned in &fragment.sources {
        let observed = reader
            .metagenomes()
            .iter()
            .find(|metagenome| metagenome.metagenome_id == planned.metagenome_id)
            .ok_or_else(|| DistributedIndexError::MissingSource(planned.metagenome_id.clone()))?;
        if hex_checksum(&observed.source_sha256) != planned.source_sha256 {
            return Err(DistributedIndexError::SourceIdentity(
                planned.metagenome_id.clone(),
            ));
        }
        let staged = staged_sources
            .get(&planned.metagenome_id)
            .ok_or_else(|| DistributedIndexError::MissingSource(planned.metagenome_id.clone()))?;
        source_results.push(IndexFragmentSourceResult {
            metagenome_id: planned.metagenome_id.clone(),
            staged_path: staged.sequence_path.clone(),
            published_path: planned.source_path.clone(),
            source_size: planned.source_size,
            source_sha256: planned.source_sha256.clone(),
            contig_count: observed.contig_count,
            total_bases: observed.total_bases,
            validation_status: "passed".to_string(),
        });
    }
    let manifest = IndexFragmentManifest {
        schema: FRAGMENT_MANIFEST_SCHEMA.to_string(),
        fragment_id,
        part_id: fragment.part_id,
        source_manifest_sha256: plan.source_manifest_sha256.clone(),
        metagenome_ids: fragment
            .sources
            .iter()
            .map(|source| source.metagenome_id.clone())
            .collect(),
        source_results,
        metagenome_count: result.metagenome_count,
        contig_count: result.contig_count,
        total_bases: result.total_bases,
        estimated_signature_count: result.estimated_signature_count,
        screen_sha256: provenance::sha256_file(&screen_path)
            .map_err(|error| DistributedIndexError::Provenance(error.to_string()))?,
        data_sha256: provenance::sha256_file(&data_path)
            .map_err(|error| DistributedIndexError::Provenance(error.to_string()))?,
        screen_bytes: screen_stats.file_size,
        data_bytes: result.data_file_bytes,
        contig_signature_histogram: result.contig_signature_histogram,
        single_contig_mappings: result.single_contig_mappings,
        overflow_mappings: result.overflow_mappings,
        overflow_contigs: result.overflow_contigs,
        maximum_overflow_count: result.maximum_overflow_count,
        signature_run_count: result.signature_run_count,
        signature_run_record_limit: result.signature_run_record_limit,
    };
    write_json_atomic(&staging.path().join("fragment.json"), &manifest)?;
    let staging_path = staging.keep();
    fs::rename(staging_path, output)?;
    Ok(manifest)
}

pub fn merge_part(
    plan: &IndexBuildPlan,
    part_id: u32,
    fragments_root: impl AsRef<Path>,
    output: impl AsRef<Path>,
) -> Result<MergedPartManifest, DistributedIndexError> {
    plan.validate()?;
    let part_plan = plan
        .part(part_id)
        .ok_or(DistributedIndexError::UnknownPart(part_id))?;
    let output = output.as_ref();
    if output.exists() {
        return Err(DistributedIndexError::ExistingOutput(output.to_path_buf()));
    }
    let parent = output.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(parent)?;
    let staging = tempfile::Builder::new()
        .prefix(&format!(".part-{part_id:06}-"))
        .tempdir_in(parent)?;
    let mut fragment_paths = Vec::new();
    let mut fragment_ids = Vec::new();
    let mut fragment_manifest_sha256 = Vec::new();
    let mut estimated_signature_count = 0u64;
    let mut published = BTreeMap::new();
    for fragment in &part_plan.fragments {
        let root = fragments_root
            .as_ref()
            .join(format!("fragment-{:06}", fragment.fragment_id));
        let manifest_path = root.join("fragment.json");
        let manifest: IndexFragmentManifest =
            serde_json::from_reader(std::io::BufReader::new(fs::File::open(&manifest_path)?))?;
        if manifest.schema != FRAGMENT_MANIFEST_SCHEMA
            || manifest.fragment_id != fragment.fragment_id
            || manifest.part_id != part_id
            || manifest.source_manifest_sha256 != plan.source_manifest_sha256
            || manifest.metagenome_ids
                != fragment
                    .sources
                    .iter()
                    .map(|source| source.metagenome_id.clone())
                    .collect::<Vec<_>>()
            || manifest.source_results.len() != fragment.sources.len()
            || manifest
                .source_results
                .iter()
                .any(|source| source.validation_status != "passed")
            || provenance::sha256_file(&root.join("screen.jam"))
                .map_err(|error| DistributedIndexError::Provenance(error.to_string()))?
                != manifest.screen_sha256
            || provenance::sha256_file(&root.join("part.bin"))
                .map_err(|error| DistributedIndexError::Provenance(error.to_string()))?
                != manifest.data_sha256
        {
            return Err(DistributedIndexError::FragmentIdentity(
                fragment.fragment_id,
            ));
        }
        fragment_paths.push((root.join("screen.jam"), root.join("part.bin")));
        estimated_signature_count =
            estimated_signature_count.saturating_add(manifest.estimated_signature_count);
        fragment_ids.push(fragment.fragment_id);
        fragment_manifest_sha256.push(
            provenance::sha256_file(&manifest_path)
                .map_err(|error| DistributedIndexError::Provenance(error.to_string()))?,
        );
        for source in &fragment.sources {
            published.insert(
                source.metagenome_id.clone(),
                PublishedMetagenomeSource {
                    metagenome_id: source.metagenome_id.clone(),
                    sequence_path: source.source_path.clone(),
                    source_size: source.source_size,
                    source_sha256: parse_checksum(&source.source_sha256)?,
                },
            );
        }
    }
    let data_path = staging.path().join("part.bin");
    let screen_path = staging.path().join("screen.jam");
    let result = merge_part_fragments(&data_path, &fragment_paths, &published)?;
    if result.estimated_signature_count != estimated_signature_count {
        return Err(DistributedIndexError::PartIdentity(part_id));
    }
    let samples = result
        .screen_samples
        .iter()
        .map(|sample| HashSampleInput {
            sample_name: sample.metagenome_id.clone(),
            hashes: sample.hashes.clone(),
        })
        .collect::<Vec<_>>();
    let screen_stats = build_hash_samples(&screen_path, &samples, 21, 1)?;
    let screen = JamReader::open(&screen_path)?;
    let data = JamIndexPartReader::open(&data_path)?;
    if screen.sample_names()
        != data
            .metagenomes()
            .iter()
            .map(|metagenome| metagenome.metagenome_id.as_str())
            .collect::<Vec<_>>()
    {
        return Err(DistributedIndexError::PartIdentity(part_id));
    }
    let descriptor = JamIndexPart {
        part_id,
        directory: format!("parts/part-{part_id:06}"),
        screen_file: "screen.jam".to_string(),
        data_file: "part.bin".to_string(),
        metagenome_count: result.metagenome_count,
        contig_count: result.contig_count,
        total_bases: result.total_bases,
        estimated_signature_count: result.estimated_signature_count,
        screen_jam_bytes: screen_stats.file_size,
        contig_posting_bytes: result.contig_posting_bytes,
        source_reference_bytes: result.source_reference_bytes,
        screen_sha256: provenance::sha256_file(&screen_path)
            .map_err(|error| DistributedIndexError::Provenance(error.to_string()))?,
        data_sha256: provenance::sha256_file(&data_path)
            .map_err(|error| DistributedIndexError::Provenance(error.to_string()))?,
    };
    let manifest = MergedPartManifest {
        schema: MERGED_PART_MANIFEST_SCHEMA.to_string(),
        part: descriptor,
        fragment_ids,
        fragment_manifest_sha256,
        contig_signature_histogram: result.contig_signature_histogram,
        single_contig_mappings: result.single_contig_mappings,
        overflow_mappings: result.overflow_mappings,
        overflow_contigs: result.overflow_contigs,
        maximum_overflow_count: result.maximum_overflow_count,
        signature_run_count: result.signature_run_count,
        signature_run_record_limit: result.signature_run_record_limit,
    };
    write_json_atomic(&staging.path().join("part.json"), &manifest)?;
    let staging_path = staging.keep();
    fs::rename(staging_path, output)?;
    Ok(manifest)
}

pub fn finalize_index(
    plan: &IndexBuildPlan,
    output: impl AsRef<Path>,
) -> Result<JamIndexBuildStats, DistributedIndexError> {
    plan.validate()?;
    let output = output.as_ref();
    if output.join(JAM_INDEX_MANIFEST_FILE).exists() {
        return Err(DistributedIndexError::ExistingOutput(
            output.join(JAM_INDEX_MANIFEST_FILE),
        ));
    }
    let mut parts = Vec::with_capacity(plan.parts.len());
    for planned in &plan.parts {
        let root = output.join(format!("parts/part-{:06}", planned.part_id));
        let manifest: MergedPartManifest = serde_json::from_reader(std::io::BufReader::new(
            fs::File::open(root.join("part.json"))?,
        ))?;
        if manifest.schema != MERGED_PART_MANIFEST_SCHEMA
            || manifest.part.part_id != planned.part_id
            || provenance::sha256_file(&root.join(&manifest.part.screen_file))
                .map_err(|error| DistributedIndexError::Provenance(error.to_string()))?
                != manifest.part.screen_sha256
            || provenance::sha256_file(&root.join(&manifest.part.data_file))
                .map_err(|error| DistributedIndexError::Provenance(error.to_string()))?
                != manifest.part.data_sha256
        {
            return Err(DistributedIndexError::PartIdentity(planned.part_id));
        }
        parts.push(manifest.part);
    }
    parts.sort_by_key(|part| part.part_id);
    let mut manifest = JamIndexManifest::empty(plan.selection_policy.clone());
    manifest.append_parts(plan.source_manifest_sha256.clone(), parts.clone())?;
    write_manifest_atomic(output, &manifest)?;
    Ok(JamIndexBuildStats {
        new_parts: u32::try_from(parts.len()).unwrap_or(u32::MAX),
        total_parts: u32::try_from(parts.len()).unwrap_or(u32::MAX),
        new_metagenomes: manifest.total_metagenomes,
        total_metagenomes: manifest.total_metagenomes,
        new_contigs: manifest.total_contigs,
        total_contigs: manifest.total_contigs,
        new_bases: manifest.total_bases,
        total_bases: manifest.total_bases,
        screen_jam_bytes: parts.iter().map(|part| part.screen_jam_bytes).sum(),
        contig_posting_bytes: parts.iter().map(|part| part.contig_posting_bytes).sum(),
        source_reference_bytes: parts.iter().map(|part| part.source_reference_bytes).sum(),
    })
}

fn write_json_atomic<T: Serialize>(path: &Path, value: &T) -> Result<(), DistributedIndexError> {
    let parent = path.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(parent)?;
    let mut temporary = tempfile::NamedTempFile::new_in(parent)?;
    serde_json::to_writer_pretty(&mut temporary, value)?;
    temporary.write_all(b"\n")?;
    temporary.as_file_mut().sync_all()?;
    temporary
        .persist(path)
        .map_err(|error| DistributedIndexError::Io(error.error))?;
    Ok(())
}

pub(crate) fn parse_checksum(value: &str) -> Result<[u8; 32], DistributedIndexError> {
    if !valid_checksum(value) {
        return Err(DistributedIndexError::InvalidPlan(
            "invalid source checksum".to_string(),
        ));
    }
    let mut checksum = [0u8; 32];
    for (index, byte) in checksum.iter_mut().enumerate() {
        *byte = u8::from_str_radix(&value[index * 2..index * 2 + 2], 16).map_err(|_| {
            DistributedIndexError::InvalidPlan("invalid source checksum".to_string())
        })?;
    }
    Ok(checksum)
}

fn hex_checksum(value: &[u8; 32]) -> String {
    value.iter().map(|byte| format!("{byte:02x}")).collect()
}

fn valid_checksum(value: &str) -> bool {
    value.len() == 64 && value.bytes().all(|byte| byte.is_ascii_hexdigit())
}

#[derive(Debug, Error)]
pub enum DistributedIndexError {
    #[error("distributed Jam Index I/O failed: {0}")]
    Io(#[from] std::io::Error),
    #[error("distributed Jam Index JSON failed: {0}")]
    Json(#[from] serde_json::Error),
    #[error("invalid distributed Jam Index plan: {0}")]
    InvalidPlan(String),
    #[error("unknown build fragment {0}")]
    UnknownFragment(u32),
    #[error("unknown planned part {0}")]
    UnknownPart(u32),
    #[error("staged source is missing for {0}")]
    MissingSource(String),
    #[error("source identity changed for {0}")]
    SourceIdentity(String),
    #[error("fragment {0} identity validation failed")]
    FragmentIdentity(u32),
    #[error("part {0} identity validation failed")]
    PartIdentity(u32),
    #[error("distributed output already exists: {0}")]
    ExistingOutput(PathBuf),
    #[error("distributed index coordinate overflow")]
    Overflow,
    #[error("distributed provenance failed: {0}")]
    Provenance(String),
    #[error(transparent)]
    Part(#[from] JamIndexPartError),
    #[error(transparent)]
    ScreenWrite(#[from] CompactError),
    #[error(transparent)]
    ScreenRead(#[from] ReaderError),
    #[error(transparent)]
    Manifest(#[from] super::manifest::JamIndexManifestError),
    #[error(transparent)]
    Build(#[from] super::builder::JamIndexBuildError),
}
