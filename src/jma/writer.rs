//! Deterministic writer for the self-contained JMA format 1 object.

use crate::jma::contigs::encode_contigs;
use crate::jma::format::{
    ArchiveMetadataFields, SECTION_CONTIGS, SECTION_FLAG_REQUIRED, SECTION_METADATA, SECTION_SEEDS,
    SECTION_SEQUENCE_DIRECTORY, SECTION_SEQUENCE_PAYLOAD, SECTION_SKETCH, SUPERBLOCK_SIZE,
    SectionDescriptor, checksum_object, encode_header, encode_metadata, encode_section_directory,
};
use crate::jma::index::{
    encode_seed_index, encode_shared_sequence_directory, rebase_shared_sequence_directory,
    seed_directory_prefix_length,
};
use crate::jma::{
    ArchiveHeader, ContigMetadata, JMA_FORMAT_VERSION, JmaError, JmaResult, SeedLevel,
    SeedOccurrence, SeedQuery,
};
use crate::sequence::EncodedSequenceBlock;
use std::collections::BTreeSet;
use std::io::Write;

/// One exact seed occurrence retained by an in-memory archive build.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedRecord {
    pub query: SeedQuery,
    pub occurrence: SeedOccurrence,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SeedLevelData {
    pub level: SeedLevel,
    pub records: Vec<SeedRecord>,
}

/// In-memory seed data. The writer turns every level into an open descriptor
/// and page directory; there is no fixed two-level header restriction.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SeedSection {
    pub k: u8,
    pub levels: Vec<SeedLevelData>,
}

/// Inputs to the deterministic one-object writer.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ArchiveParts {
    pub flags: u32,
    pub source_sha256: [u8; 32],
    pub contigs: Vec<ContigMetadata>,
    pub sequence_blocks: Vec<EncodedSequenceBlock>,
    pub seed_sections: Vec<SeedSection>,
}

pub type JmaArchive = ArchiveParts;

/// Encodes one deterministic, self-contained JMA object.
pub fn encode_archive(parts: &ArchiveParts) -> JmaResult<Vec<u8>> {
    encode_archive_with_min_entropy(parts, None)
}

pub fn encode_archive_with_min_entropy(
    parts: &ArchiveParts,
    min_entropy: Option<f64>,
) -> JmaResult<Vec<u8>> {
    validate_contig_table(&parts.contigs)?;
    validate_complete_coverage(&parts.sequence_blocks, &parts.contigs)?;
    let contig_count = u32::try_from(parts.contigs.len()).map_err(|_| JmaError::OffsetOverflow)?;
    let total_bases = parts.contigs.iter().try_fold(0u64, |total, contig| {
        total
            .checked_add(contig.length)
            .ok_or(JmaError::OffsetOverflow)
    })?;

    let mut metadata = ArchiveMetadataFields::for_source(parts.source_sha256);
    metadata.min_entropy = min_entropy;
    let metadata_payload = encode_metadata(&metadata)?;
    let contig_payload = encode_contigs(&parts.contigs)?;
    let seed_payload = encode_seed_index(&parts.seed_sections)?;
    let (sequence_directory_payload, sequence_payload) =
        encode_shared_sequence_directory(&parts.sequence_blocks, &parts.contigs)?;
    let sketch_payload = encode_sketch(&parts.seed_sections)?;

    let mut payloads = vec![
        (SECTION_METADATA, metadata_payload),
        (SECTION_CONTIGS, contig_payload),
        (SECTION_SKETCH, sketch_payload),
        (SECTION_SEEDS, seed_payload),
        (SECTION_SEQUENCE_DIRECTORY, sequence_directory_payload),
        (SECTION_SEQUENCE_PAYLOAD, sequence_payload),
    ];
    let section_count = u32::try_from(payloads.len()).map_err(|_| JmaError::OffsetOverflow)?;
    let directory_offset = u64::try_from(SUPERBLOCK_SIZE).map_err(|_| JmaError::OffsetOverflow)?;
    let directory_length = u64::try_from(
        payloads
            .len()
            .checked_mul(crate::jma::format::SECTION_ENTRY_SIZE)
            .ok_or(JmaError::OffsetOverflow)?,
    )
    .map_err(|_| JmaError::OffsetOverflow)?;
    let mut next_offset = directory_offset
        .checked_add(directory_length)
        .ok_or(JmaError::OffsetOverflow)?;
    let sequence_payload_offset =
        payloads
            .iter()
            .take(5)
            .try_fold(next_offset, |offset, (_, payload)| {
                offset
                    .checked_add(
                        u64::try_from(payload.len()).map_err(|_| JmaError::OffsetOverflow)?,
                    )
                    .ok_or(JmaError::OffsetOverflow)
            })?;
    payloads[4].1 = rebase_shared_sequence_directory(&payloads[4].1, sequence_payload_offset)?;
    let mut descriptors = Vec::with_capacity(payloads.len());
    for (kind, payload) in &payloads {
        let length = u64::try_from(payload.len()).map_err(|_| JmaError::OffsetOverflow)?;
        let section_checksum = if *kind == SECTION_SEEDS {
            let prefix = seed_directory_prefix_length(payload)?;
            crate::jma::format::checksum(payload.get(..prefix).ok_or_else(|| {
                JmaError::CorruptSection("seed directory prefix is truncated".to_string())
            })?)
        } else {
            crate::jma::format::checksum(payload)
        };
        descriptors.push(SectionDescriptor {
            kind: *kind,
            flags: SECTION_FLAG_REQUIRED,
            offset: next_offset,
            length,
            uncompressed_length: length,
            checksum: section_checksum,
        });
        next_offset = next_offset
            .checked_add(length)
            .ok_or(JmaError::OffsetOverflow)?;
    }

    let archive_header = ArchiveHeader {
        format_version: JMA_FORMAT_VERSION,
        flags: parts.flags,
        contig_count,
        total_bases,
        source_sha256: parts.source_sha256,
        seed_levels: parts
            .seed_sections
            .iter()
            .flat_map(|section| section.levels.iter().map(|level| level.level))
            .collect(),
        algorithm_id: Some("jam-seed-chain-align-v1".to_string()),
        algorithm_version: Some(1),
        min_entropy,
    };
    let directory = encode_section_directory(&descriptors)?;
    let header = encode_header(
        &archive_header,
        section_count,
        directory_offset,
        directory_length,
        next_offset,
        crate::jma::format::checksum(&directory),
    )?;
    let capacity = usize::try_from(next_offset).map_err(|_| JmaError::OffsetOverflow)?;
    let mut output = Vec::with_capacity(capacity);
    output.extend_from_slice(&header);
    output.extend_from_slice(&directory);
    for (_, payload) in payloads {
        output.extend_from_slice(&payload);
    }
    if u64::try_from(output.len()).map_err(|_| JmaError::OffsetOverflow)? != next_offset {
        return Err(JmaError::CorruptSection(
            "encoded object length does not match superblock".to_string(),
        ));
    }
    // Finalize the two non-recursive checksum fields in deterministic order.
    let object_digest = checksum_object(&output);
    output[96..128].copy_from_slice(&object_digest);
    output[128..160].fill(0);
    let superblock_digest = crate::jma::format::checksum(&output[..SUPERBLOCK_SIZE]);
    output[128..160].copy_from_slice(&superblock_digest);
    Ok(output)
}

pub fn write_archive<W: Write>(writer: &mut W, parts: &ArchiveParts) -> JmaResult<()> {
    let bytes = encode_archive(parts)?;
    writer
        .write_all(&bytes)
        .map_err(|error| JmaError::CorruptSection(format!("archive write failed: {error}")))
}

pub fn write_archive_with_min_entropy<W: Write>(
    writer: &mut W,
    parts: &ArchiveParts,
    min_entropy: Option<f64>,
) -> JmaResult<()> {
    let bytes = encode_archive_with_min_entropy(parts, min_entropy)?;
    writer
        .write_all(&bytes)
        .map_err(|error| JmaError::CorruptSection(format!("archive write failed: {error}")))
}

fn encode_sketch(sections: &[SeedSection]) -> JmaResult<Vec<u8>> {
    let mut hashes = BTreeSet::new();
    for section in sections {
        validate_seed_section(section)?;
        for level in &section.levels {
            for record in &level.records {
                if record.query.hash != 0 {
                    hashes.insert(record.query.hash);
                }
            }
        }
    }
    let mut bytes = Vec::with_capacity(16 + hashes.len().saturating_mul(8));
    bytes.extend_from_slice(&1u32.to_le_bytes());
    bytes.extend_from_slice(&0u32.to_le_bytes());
    bytes.extend_from_slice(
        &u64::try_from(hashes.len())
            .map_err(|_| JmaError::OffsetOverflow)?
            .to_le_bytes(),
    );
    for hash in hashes {
        bytes.extend_from_slice(&hash.to_le_bytes());
    }
    Ok(bytes)
}

fn validate_contig_table(contigs: &[ContigMetadata]) -> JmaResult<()> {
    let mut previous = None;
    for contig in contigs {
        if previous.is_some_and(|previous_id| contig.id <= previous_id) {
            return Err(JmaError::CorruptSection(
                "contig identifiers must be strictly increasing".to_string(),
            ));
        }
        previous = Some(contig.id);
    }
    Ok(())
}

fn validate_complete_coverage(
    blocks: &[EncodedSequenceBlock],
    contigs: &[ContigMetadata],
) -> JmaResult<()> {
    for contig in contigs {
        let mut next = 0u64;
        let mut contig_blocks = blocks
            .iter()
            .filter(|block| block.record.contig_id == contig.id)
            .collect::<Vec<_>>();
        contig_blocks.sort_by_key(|block| block.record.base_start);
        for block in contig_blocks {
            if block.record.base_start != next {
                return Err(JmaError::CorruptSection(format!(
                    "sequence blocks leave a gap or overlap in contig {} at {}",
                    contig.id, next
                )));
            }
            next = block.record.base_end().ok_or(JmaError::OffsetOverflow)?;
        }
        if next != contig.length {
            return Err(JmaError::CorruptSection(format!(
                "sequence blocks cover {next} of contig {}, expected {}",
                contig.id, contig.length
            )));
        }
    }
    Ok(())
}

fn validate_seed_section(section: &SeedSection) -> JmaResult<()> {
    if !(21..=32).contains(&section.k) {
        return Err(JmaError::CorruptSection(format!(
            "seed span must be between 21 and 32, got k={}",
            section.k
        )));
    }
    let mut identities = BTreeSet::new();
    for level in &section.levels {
        if level.level.k != section.k
            || level.level.scale == 0
            || !identities.insert(level.level.scale)
        {
            return Err(JmaError::CorruptSection(
                "seed level identity or density parameter is invalid".to_string(),
            ));
        }
        for record in &level.records {
            if record.query.k != section.k
                || !crate::jma::format::valid_canonical_kmer(section.k, record.query.canonical_kmer)
            {
                return Err(JmaError::CorruptSection(
                    "seed record query does not match its section".to_string(),
                ));
            }
        }
    }
    Ok(())
}
