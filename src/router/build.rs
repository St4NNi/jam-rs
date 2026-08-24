//! Deterministic collection-router construction from local JMA objects.

use crate::archive::SeedSchemeId;
use crate::jma::reader::JmaReader;
use crate::resource::ResourceOpenOptions;
use crate::resource::local::LocalResource;
use crate::router::format::{MetagenomeEntry, PostingHeader};
use crate::router::postings::{PositionPosting, encode_document_ids_wire, encode_positions_wire};
use crate::router::writer::{RouterBuildInput, RouterKeyInput, RouterWriterOptions, write_router};
use crate::router::{HashAlgorithmId, WITNESS_K, WitnessKey, WitnessScheme};
use std::collections::BTreeMap;
use std::path::Path;
use thiserror::Error;

#[derive(Clone, Debug)]
pub struct RouterCollectionBuildConfig {
    pub base_scale: u32,
    pub available_scales: Vec<u32>,
    pub rare_document_frequency: u32,
    pub common_document_frequency: u32,
    pub positions_per_metagenome: usize,
    pub writer: RouterWriterOptions,
}

impl Default for RouterCollectionBuildConfig {
    fn default() -> Self {
        Self {
            base_scale: 20,
            available_scales: vec![20, 50, 100, 200, 500],
            rare_document_frequency: 4,
            common_document_frequency: 16,
            positions_per_metagenome: 4,
            writer: RouterWriterOptions::default(),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct RouterCollectionBuildStats {
    pub metagenomes: u32,
    pub keys: u64,
    pub posting_ids: u64,
    pub positional_occurrences: u64,
    pub bytes: u64,
}

#[derive(Default)]
struct Aggregate {
    documents: Vec<u32>,
    positions: Vec<PositionPosting>,
}

pub fn build_router_from_local_jma(
    metagenomes: &[MetagenomeEntry],
    output: impl AsRef<Path>,
    config: &RouterCollectionBuildConfig,
) -> Result<RouterCollectionBuildStats, RouterBuildError> {
    build_router_from_local_jma_sources(metagenomes, metagenomes, output, config)
}

pub fn build_router_from_local_jma_sources(
    witness_sources: &[MetagenomeEntry],
    metagenomes: &[MetagenomeEntry],
    output: impl AsRef<Path>,
    config: &RouterCollectionBuildConfig,
) -> Result<RouterCollectionBuildStats, RouterBuildError> {
    if metagenomes.is_empty() {
        return Err(RouterBuildError::InvalidConfig(
            "metagenome catalog is empty".to_string(),
        ));
    }
    if witness_sources.len() != metagenomes.len()
        || witness_sources
            .iter()
            .zip(metagenomes)
            .any(|(source, bound)| source.id != bound.id)
    {
        return Err(RouterBuildError::InvalidConfig(
            "witness and bound metagenome catalogs differ in ID order".to_string(),
        ));
    }
    if config.positions_per_metagenome == 0 {
        return Err(RouterBuildError::InvalidConfig(
            "positions_per_metagenome must be greater than zero".to_string(),
        ));
    }
    let mut scales = config.available_scales.clone();
    scales.sort_unstable();
    scales.dedup();
    if scales.first().copied() != Some(config.base_scale) {
        return Err(RouterBuildError::InvalidConfig(
            "available scales must start at base_scale".to_string(),
        ));
    }

    let mut aggregates = BTreeMap::<WitnessKey, Aggregate>::new();
    let mut shared_scheme_id = None;
    for (document_index, metagenome) in witness_sources.iter().enumerate() {
        let observed_checksum = crate::provenance::sha256_file(Path::new(&metagenome.object_uri))
            .map_err(|error| RouterBuildError::Resource(error.to_string()))?;
        if observed_checksum != hex(&metagenome.checksum) {
            return Err(RouterBuildError::Resource(format!(
                "JMA checksum mismatch for {}",
                metagenome.id
            )));
        }
        let resource = LocalResource::from_path(
            &metagenome.object_uri,
            ResourceOpenOptions {
                expected_sha256: Some(hex(&metagenome.checksum)),
                allow_full_download_fallback: false,
                ..ResourceOpenOptions::default()
            },
        )
        .map_err(|error| RouterBuildError::Resource(error.to_string()))?;
        let reader = JmaReader::from_resource(resource)
            .map_err(|error| RouterBuildError::Jma(error.to_string()))?;
        let descriptor = reader
            .seed_schemes()
            .find(|descriptor| {
                descriptor.span == u16::from(WITNESS_K)
                    && descriptor.density_parameter == config.base_scale
            })
            .copied()
            .ok_or_else(|| RouterBuildError::MissingWitnessScheme {
                metagenome_id: metagenome.id.clone(),
                scale: config.base_scale,
            })?;
        if shared_scheme_id.is_some_and(|scheme_id| scheme_id != descriptor.scheme_id) {
            return Err(RouterBuildError::SchemeMismatch);
        }
        shared_scheme_id = Some(descriptor.scheme_id);
        let records = reader
            .export_seed_records(SeedSchemeId(descriptor.scheme_id))
            .map_err(|error| RouterBuildError::Jma(error.to_string()))?;
        let document_id =
            u32::try_from(document_index).map_err(|_| RouterBuildError::CountOverflow)?;
        for record in records {
            let key = WitnessKey::checked(record.query.canonical_kmer, record.query.hash)
                .map_err(|error| RouterBuildError::Jma(error.to_string()))?;
            let aggregate = aggregates.entry(key).or_default();
            aggregate.documents.push(document_id);
            aggregate.positions.extend(
                record
                    .occurrences
                    .into_iter()
                    .take(config.positions_per_metagenome)
                    .map(|occurrence| PositionPosting {
                        metagenome_id: document_id,
                        contig_id: occurrence.contig_id,
                        position: occurrence.position,
                        reverse: occurrence.reverse,
                    }),
            );
        }
    }

    let scheme_id = shared_scheme_id.ok_or(RouterBuildError::SchemeMismatch)?;
    let scheme = WitnessScheme {
        scheme_id,
        k: WITNESS_K,
        base_scale: config.base_scale,
        available_scales: scales,
        hash_id: HashAlgorithmId::JamhashU64V1,
        zero_excluded: true,
    };
    scheme
        .validate()
        .map_err(|error| RouterBuildError::InvalidConfig(error.to_string()))?;
    let mut keys = Vec::with_capacity(aggregates.len());
    let mut stats = RouterCollectionBuildStats {
        metagenomes: u32::try_from(metagenomes.len())
            .map_err(|_| RouterBuildError::CountOverflow)?,
        ..RouterCollectionBuildStats::default()
    };
    for (key, mut aggregate) in aggregates {
        aggregate.documents.sort_unstable();
        aggregate.documents.dedup();
        aggregate.positions.sort_unstable();
        aggregate.positions.dedup();
        let document_frequency = u32::try_from(aggregate.documents.len())
            .map_err(|_| RouterBuildError::CountOverflow)?;
        let mut flags = 0u32;
        if document_frequency > config.common_document_frequency {
            flags |= PostingHeader::FLAG_COMMON;
        }
        let positions = if document_frequency <= config.rare_document_frequency {
            flags |= PostingHeader::FLAG_POSITION_BEARING;
            Some(
                encode_positions_wire(&aggregate.positions)
                    .map_err(|error| RouterBuildError::Posting(error.to_string()))?,
            )
        } else {
            None
        };
        stats.posting_ids = stats
            .posting_ids
            .saturating_add(u64::from(document_frequency));
        if positions.is_some() {
            stats.positional_occurrences = stats
                .positional_occurrences
                .saturating_add(u64::try_from(aggregate.positions.len()).unwrap_or(u64::MAX));
        }
        keys.push(RouterKeyInput {
            key,
            document_frequency,
            flags,
            posting: encode_document_ids_wire(&aggregate.documents)
                .map_err(|error| RouterBuildError::Posting(error.to_string()))?,
            positions,
        });
    }
    stats.keys = u64::try_from(keys.len()).map_err(|_| RouterBuildError::CountOverflow)?;
    let input = RouterBuildInput {
        metagenomes: metagenomes.to_vec(),
        schemes: vec![scheme],
        keys,
    };
    write_router(output.as_ref(), &input, config.writer)
        .map_err(|error| RouterBuildError::Format(error.to_string()))?;
    stats.bytes = output
        .as_ref()
        .metadata()
        .map_err(|error| RouterBuildError::Resource(error.to_string()))?
        .len();
    Ok(stats)
}

fn hex(bytes: &[u8]) -> String {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut output = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        output.push(DIGITS[(byte >> 4) as usize] as char);
        output.push(DIGITS[(byte & 0x0f) as usize] as char);
    }
    output
}

#[derive(Debug, Error)]
pub enum RouterBuildError {
    #[error("invalid router build configuration: {0}")]
    InvalidConfig(String),
    #[error("JMA resource failed: {0}")]
    Resource(String),
    #[error("JMA seed export failed: {0}")]
    Jma(String),
    #[error("router posting encode failed: {0}")]
    Posting(String),
    #[error("router format failed: {0}")]
    Format(String),
    #[error("metagenome {metagenome_id} has no k=21 scale-{scale} witness scheme")]
    MissingWitnessScheme { metagenome_id: String, scale: u32 },
    #[error("JMA witness scheme identifiers differ across the catalog")]
    SchemeMismatch,
    #[error("router collection count overflow")]
    CountOverflow,
}
