//! Adapters from the in-memory JMA seed records to the compact section.

use crate::jma::seeds::{
    SeedBinding, SeedBuildConfig, encode_seed_collection, encode_seed_section,
};
use crate::jma::writer::{SeedLevelData, SeedSection};
use crate::jma::{JmaError, JmaResult};

/// Encode one existing JMA seed level using the compact key/header/payload
/// layout. The caller owns archive/catalog binding values.
pub fn encode_seed_level(
    section: &SeedSection,
    level: &SeedLevelData,
    binding: SeedBinding,
    config: SeedBuildConfig,
) -> JmaResult<Vec<u8>> {
    if section.k != level.level.k
        || section.k != binding.k
        || level.level.scale != u64::from(binding.scale)
    {
        return Err(JmaError::CorruptSection(
            "compact seed level does not match its binding".to_string(),
        ));
    }
    encode_seed_section(&level.records, binding, config)
}

/// Encode all levels of one k-mer section. Each level remains independently
/// bound and independently addressable by its caller.
pub fn encode_seed_levels(
    section: &SeedSection,
    bindings: &[(SeedBinding, SeedBuildConfig)],
) -> JmaResult<Vec<Vec<u8>>> {
    if section.levels.len() != bindings.len() {
        return Err(JmaError::CorruptSection(
            "compact seed binding count does not match levels".to_string(),
        ));
    }
    section
        .levels
        .iter()
        .zip(bindings.iter().copied())
        .map(|(level, (binding, config))| encode_seed_level(section, level, binding, config))
        .collect()
}

/// Encode all levels of one JMA section into the single collection payload
/// required by `SECTION_SEEDS`.
pub fn encode_seed_collection_levels(
    section: &SeedSection,
    bindings: &[(SeedBinding, SeedBuildConfig)],
) -> JmaResult<Vec<u8>> {
    let blobs = encode_seed_levels(section, bindings)?;
    let bound = bindings
        .iter()
        .map(|(binding, _)| *binding)
        .zip(blobs)
        .collect::<Vec<_>>();
    encode_seed_collection(&bound)
}

/// Encode every level from all JMA seed sections into one `SECTION_SEEDS`
/// collection. Bindings are consumed in deterministic section/level order;
/// the collection directory then sorts them by scheme identity.
pub fn encode_seed_collection_sections(
    sections: &[SeedSection],
    bindings: &[(SeedBinding, SeedBuildConfig)],
) -> JmaResult<Vec<u8>> {
    let level_count = sections
        .iter()
        .map(|section| section.levels.len())
        .sum::<usize>();
    if level_count != bindings.len() {
        return Err(JmaError::CorruptSection(
            "compact seed binding count does not match all levels".to_string(),
        ));
    }
    let mut bound = Vec::with_capacity(level_count);
    let mut binding_index = 0usize;
    for section in sections {
        for level in &section.levels {
            let (binding, config) = bindings.get(binding_index).copied().ok_or_else(|| {
                JmaError::CorruptSection("compact seed binding is unavailable".to_string())
            })?;
            bound.push((binding, encode_seed_level(section, level, binding, config)?));
            binding_index += 1;
        }
    }
    encode_seed_collection(&bound)
}
