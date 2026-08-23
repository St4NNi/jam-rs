//! Checksum-bound sidecar indexes for JMA v1 archives.
//!
//! A sidecar contains byte ranges for the encoded sequence blocks and for
//! fixed-size seed-record buckets.  It is deliberately separate from the JMA
//! binary format: existing archives remain readable through the eager reader,
//! while indexed readers can fetch only the blocks needed by a query.

use crate::jma::contigs::decode_contigs;
use crate::jma::format::{
    HEADER_SIZE, ParsedHeader, SECTION_CONTIGS, SECTION_SEEDS_K21, SECTION_SEEDS_K31,
    SECTION_SEQUENCES, SectionDescriptor, checked_end, checked_usize, checksum, get_u32, get_u64,
};
use crate::jma::header::{parse_header, parse_section_directory};
use crate::jma::sequence::decode_sequence_block;
use crate::jma::writer::decode_seed_section;
use crate::jma::{ArchiveHeader, ContigMetadata, JmaError, JmaResult, SeedLevel};
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::{Path, PathBuf};

/// Version of the JSON sidecar schema.
pub const INDEX_FORMAT_VERSION: u16 = 1;
/// Number of high-order hash bits used to group fixed-size seed records.
pub const HASH_PREFIX_BITS: u8 = 16;
/// Fixed encoded size of one seed occurrence record.
pub const SEED_RECORD_SIZE: u64 = 32;
const SEQUENCE_SECTION_HEADER_SIZE: u64 = 8;
const SEED_SECTION_HEADER_SIZE: u64 = 12;
const SEED_LEVEL_HEADER_SIZE: u64 = 16;
const SEQUENCE_BLOCK_HEADER_SIZE: u64 = 36;
pub(crate) const MAX_INDEX_BYTES: u64 = 256 * 1024 * 1024;
pub(crate) const MAX_INDEXED_RANGE_BYTES: u64 = 256 * 1024 * 1024;

type ArchiveLayout<'a> = (
    &'a [u8],
    ParsedHeader,
    Vec<SectionDescriptor>,
    Vec<ContigMetadata>,
);

/// A checksum-bound range for one encoded sequence block.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SequenceBlockIndex {
    pub contig_id: u32,
    pub start: u64,
    pub base_length: u64,
    /// Absolute byte offset in the JMA archive, including the block header.
    pub offset: u64,
    pub length: u64,
    pub checksum: [u8; 32],
}

/// Structural metadata for one encoded seed density level.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SeedLevelIndex {
    pub k: u8,
    pub scale: u64,
    pub section_offset: u64,
    pub section_length: u64,
    pub level_header_offset: u64,
    pub record_offset: u64,
    pub record_count: u64,
    pub record_length: u64,
    pub checksum: [u8; 32],
}

/// A contiguous fixed-record range selected by the high-order hash prefix.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SeedBucketIndex {
    pub k: u8,
    pub scale: u64,
    pub hash_prefix: u16,
    /// Absolute byte offset in the JMA archive.
    pub offset: u64,
    pub length: u64,
    pub record_count: u64,
    pub checksum: [u8; 32],
}

/// JSON sidecar metadata for one JMA v1 archive.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct JmaSidecarIndex {
    pub index_version: u16,
    pub archive_size: u64,
    /// SHA-256 of the complete JMA archive, recorded for manifest and cache
    /// identity. Indexed opening does not reread the complete archive merely
    /// to recompute this value.
    pub archive_sha256: [u8; 32],
    pub archive_header_sha256: [u8; 32],
    pub archive_format_version: u16,
    pub archive_source_sha256: [u8; 32],
    pub algorithm_id: Option<String>,
    pub algorithm_version: Option<u16>,
    pub hash_prefix_bits: u8,
    pub sequence_blocks: Vec<SequenceBlockIndex>,
    pub seed_levels: Vec<SeedLevelIndex>,
    pub seed_buckets: Vec<SeedBucketIndex>,
}

/// Returns the deterministic sidecar path for a JMA archive.
#[must_use]
pub fn sidecar_path<P: AsRef<Path>>(archive: P) -> PathBuf {
    let archive = archive.as_ref();
    let mut name = archive.as_os_str().to_os_string();
    name.push(".idx.json");
    PathBuf::from(name)
}

/// Builds a deterministic sidecar index from complete encoded archive bytes.
///
/// This is used by the local archive builder, which already materializes the
/// encoded output before writing it.  It is also useful for tests and for
/// callers that generate an archive through [`crate::jma::writer`].
pub fn build_index(bytes: &[u8]) -> JmaResult<JmaSidecarIndex> {
    let (header_bytes, parsed, sections, contigs) = parse_archive_layout(bytes)?;
    let sequence_descriptor = unique_section(&sections, SECTION_SEQUENCES)?;
    let sequence_payload = section_bytes(bytes, sequence_descriptor)?;
    let sequence_blocks = scan_sequence_blocks(sequence_payload, sequence_descriptor)?;

    let mut seed_levels = Vec::new();
    let mut seed_buckets = Vec::new();
    for kind in [SECTION_SEEDS_K21, SECTION_SEEDS_K31] {
        if let Some(descriptor) = sections.iter().find(|entry| entry.kind == kind) {
            let payload = section_bytes(bytes, descriptor)?;
            let (levels, buckets) = scan_seed_section(payload, descriptor)?;
            seed_levels.extend(levels);
            seed_buckets.extend(buckets);
        }
    }

    // Validate all generated block metadata against the contig table before
    // exposing it to an indexed reader.
    validate_sequence_index(&sequence_blocks, sequence_descriptor, &contigs)?;
    validate_seed_index(&seed_levels, &seed_buckets, &sections, &parsed.archive)?;

    Ok(JmaSidecarIndex {
        index_version: INDEX_FORMAT_VERSION,
        archive_size: u64::try_from(bytes.len()).map_err(|_| JmaError::OffsetOverflow)?,
        archive_sha256: checksum(bytes),
        archive_header_sha256: checksum(header_bytes),
        archive_format_version: parsed.archive.format_version,
        archive_source_sha256: parsed.archive.source_sha256,
        algorithm_id: parsed.archive.algorithm_id,
        algorithm_version: parsed.archive.algorithm_version,
        hash_prefix_bits: HASH_PREFIX_BITS,
        sequence_blocks,
        seed_levels,
        seed_buckets,
    })
}

/// Serializes a sidecar with stable field order and a trailing newline.
pub fn encode_index(index: &JmaSidecarIndex) -> JmaResult<Vec<u8>> {
    if index.index_version != INDEX_FORMAT_VERSION {
        return Err(JmaError::UnsupportedVersion(index.index_version));
    }
    let mut bytes = serde_json::to_vec_pretty(index)
        .map_err(|error| JmaError::CorruptSection(format!("cannot encode JMA sidecar: {error}")))?;
    bytes.push(b'\n');
    Ok(bytes)
}

/// Parses a sidecar without performing archive-specific binding checks.
pub fn decode_index(bytes: &[u8]) -> JmaResult<JmaSidecarIndex> {
    if bytes.is_empty() {
        return Err(JmaError::CorruptSection("JMA sidecar is empty".to_string()));
    }
    let index: JmaSidecarIndex = serde_json::from_slice(bytes).map_err(|error| {
        JmaError::CorruptSection(format!("cannot decode JMA sidecar JSON: {error}"))
    })?;
    if index.index_version != INDEX_FORMAT_VERSION {
        return Err(JmaError::UnsupportedVersion(index.index_version));
    }
    if index.hash_prefix_bits != HASH_PREFIX_BITS {
        return Err(JmaError::CorruptSection(format!(
            "unsupported JMA sidecar hash prefix width {}",
            index.hash_prefix_bits
        )));
    }
    Ok(index)
}

/// Writes the sidecar associated with an already encoded archive path.
///
/// Only the explicitly derived `archive.jma.idx.json` path is touched.  The
/// archive itself and other sibling files are never replaced.
pub fn write_sidecar<P: AsRef<Path>>(archive: P, bytes: &[u8]) -> JmaResult<PathBuf> {
    let index = build_index(bytes)?;
    let encoded = encode_index(&index)?;
    let path = sidecar_path(archive);
    fs::write(&path, encoded).map_err(|error| {
        JmaError::CorruptSection(format!(
            "cannot write JMA sidecar {}: {error}",
            path.display()
        ))
    })?;
    Ok(path)
}

/// Reads and indexes an existing archive file, writing its deterministic
/// sidecar next to it.
pub fn write_sidecar_for_archive<P: AsRef<Path>>(archive: P) -> JmaResult<PathBuf> {
    let archive = archive.as_ref();
    let bytes = fs::read(archive).map_err(|error| {
        JmaError::CorruptSection(format!(
            "cannot read JMA archive {}: {error}",
            archive.display()
        ))
    })?;
    write_sidecar(archive, &bytes)
}

/// Validates sidecar identity and all index structure against an archive's
/// fixed header, directory, and contig table.  Payload checksums are checked
/// when an indexed bucket or sequence block is fetched.
pub fn validate_against_archive(
    index: &JmaSidecarIndex,
    archive_size: u64,
    header_bytes: &[u8],
    archive: &ArchiveHeader,
    sections: &[SectionDescriptor],
    contigs: &[ContigMetadata],
) -> JmaResult<()> {
    if index.index_version != INDEX_FORMAT_VERSION {
        return Err(JmaError::UnsupportedVersion(index.index_version));
    }
    if index.hash_prefix_bits != HASH_PREFIX_BITS {
        return Err(JmaError::CorruptSection(format!(
            "unsupported JMA sidecar hash prefix width {}",
            index.hash_prefix_bits
        )));
    }
    if index.archive_size != archive_size {
        return Err(JmaError::CorruptSection(format!(
            "JMA sidecar archive size {} does not match {}",
            index.archive_size, archive_size
        )));
    }
    if index.archive_format_version != archive.format_version {
        return Err(JmaError::CorruptSection(
            "JMA sidecar archive format version does not match archive".to_string(),
        ));
    }
    if index.archive_header_sha256 != checksum(header_bytes) {
        return Err(JmaError::ChecksumMismatch(
            "JMA sidecar archive header".to_string(),
        ));
    }
    if index.archive_source_sha256 != archive.source_sha256 {
        return Err(JmaError::CorruptSection(
            "JMA sidecar source checksum does not match archive".to_string(),
        ));
    }
    if index.algorithm_id != archive.algorithm_id
        || index.algorithm_version != archive.algorithm_version
    {
        return Err(JmaError::CorruptSection(
            "JMA sidecar algorithm identity does not match archive".to_string(),
        ));
    }

    let sequence_descriptor = unique_section(sections, SECTION_SEQUENCES)?;
    validate_sequence_index(&index.sequence_blocks, sequence_descriptor, contigs)?;
    validate_seed_index(&index.seed_levels, &index.seed_buckets, sections, archive)
}

/// Finds a bucket for a query hash using the sidecar's fixed prefix rule.
#[must_use]
pub fn bucket_for(
    index: &JmaSidecarIndex,
    k: u8,
    scale: u64,
    hash: u64,
) -> Option<&SeedBucketIndex> {
    let prefix = (hash >> (64 - u32::from(HASH_PREFIX_BITS))) as u16;
    index
        .seed_buckets
        .iter()
        .find(|bucket| bucket.k == k && bucket.scale == scale && bucket.hash_prefix == prefix)
}

#[must_use]
pub fn level_for(index: &JmaSidecarIndex, k: u8, scale: u64) -> Option<&SeedLevelIndex> {
    index
        .seed_levels
        .iter()
        .find(|level| level.k == k && level.scale == scale)
}

fn parse_archive_layout(bytes: &[u8]) -> JmaResult<ArchiveLayout<'_>> {
    if bytes.len() < HEADER_SIZE {
        return Err(JmaError::CorruptSection(format!(
            "JMA archive is truncated: {} bytes, expected at least {HEADER_SIZE}",
            bytes.len()
        )));
    }
    let header_bytes = &bytes[..HEADER_SIZE];
    let parsed = parse_header(header_bytes)?;
    if parsed.section_count > 1 << 20 {
        return Err(JmaError::CorruptSection(
            "JMA section count exceeds the v1 limit".to_string(),
        ));
    }
    let archive_size = u64::try_from(bytes.len()).map_err(|_| JmaError::OffsetOverflow)?;
    let directory_end = checked_end(
        parsed.section_directory_offset,
        parsed.section_directory_length,
    )?;
    if directory_end > archive_size {
        return Err(JmaError::CorruptSection(
            "JMA section directory exceeds archive size".to_string(),
        ));
    }
    let directory_start = checked_usize(parsed.section_directory_offset, "directory offset")?;
    let directory_length = checked_usize(parsed.section_directory_length, "directory length")?;
    let directory_end_usize = directory_start
        .checked_add(directory_length)
        .ok_or(JmaError::OffsetOverflow)?;
    let directory = bytes
        .get(directory_start..directory_end_usize)
        .ok_or_else(|| JmaError::CorruptSection("JMA directory is truncated".to_string()))?;
    let sections = parse_section_directory(
        directory,
        parsed.section_count,
        parsed.section_directory_offset,
        archive_size,
    )?;
    validate_section_kinds(&sections)?;
    let contig_descriptor = unique_section(&sections, SECTION_CONTIGS)?;
    let contig_payload = section_bytes(bytes, contig_descriptor)?;
    let contigs = decode_contigs(contig_payload, parsed.archive.contig_count)?;
    let total_bases = contigs.iter().try_fold(0u64, |total, contig| {
        total
            .checked_add(contig.length)
            .ok_or(JmaError::OffsetOverflow)
    })?;
    if total_bases != parsed.archive.total_bases {
        return Err(JmaError::CorruptSection(
            "contig table total does not match archive header".to_string(),
        ));
    }
    Ok((header_bytes, parsed, sections, contigs))
}

fn section_bytes<'a>(bytes: &'a [u8], descriptor: &SectionDescriptor) -> JmaResult<&'a [u8]> {
    let start = checked_usize(descriptor.offset, "section offset")?;
    let length = checked_usize(descriptor.length, "section length")?;
    let end = start.checked_add(length).ok_or(JmaError::OffsetOverflow)?;
    let payload = bytes
        .get(start..end)
        .ok_or_else(|| JmaError::CorruptSection("JMA section is truncated".to_string()))?;
    if checksum(payload) != descriptor.checksum {
        return Err(JmaError::ChecksumMismatch(format!(
            "section kind {}",
            descriptor.kind
        )));
    }
    if descriptor.uncompressed_length != descriptor.length {
        return Err(JmaError::CorruptSection(format!(
            "compressed section kind {} is unsupported in JMA v1",
            descriptor.kind
        )));
    }
    Ok(payload)
}

fn scan_sequence_blocks(
    payload: &[u8],
    descriptor: &SectionDescriptor,
) -> JmaResult<Vec<SequenceBlockIndex>> {
    if payload.len() < 8 || get_u32(payload, 0)? != 1 {
        return Err(JmaError::CorruptSection(
            "invalid sequence section header".to_string(),
        ));
    }
    let count = get_u32(payload, 4)?;
    let mut cursor = 8usize;
    let mut result = Vec::with_capacity(count as usize);
    for _ in 0..count {
        let block_start = cursor;
        let header_end = cursor
            .checked_add(SEQUENCE_BLOCK_HEADER_SIZE as usize)
            .ok_or(JmaError::OffsetOverflow)?;
        if header_end > payload.len() {
            return Err(JmaError::CorruptSection(
                "sequence block header is truncated".to_string(),
            ));
        }
        let packed_length = get_u64(payload, cursor + 20)?;
        let mask_length = get_u64(payload, cursor + 28)?;
        cursor = header_end;
        let packed = checked_usize(packed_length, "packed sequence length")?;
        let mask = checked_usize(mask_length, "unknown mask length")?;
        cursor = cursor.checked_add(packed).ok_or(JmaError::OffsetOverflow)?;
        cursor = cursor.checked_add(mask).ok_or(JmaError::OffsetOverflow)?;
        if cursor > payload.len() {
            return Err(JmaError::CorruptSection(
                "sequence block payload is truncated".to_string(),
            ));
        }
        let encoded = &payload[block_start..cursor];
        let block = decode_sequence_block(encoded)?;
        let offset = descriptor
            .offset
            .checked_add(u64::try_from(block_start).map_err(|_| JmaError::OffsetOverflow)?)
            .ok_or(JmaError::OffsetOverflow)?;
        result.push(SequenceBlockIndex {
            contig_id: block.contig_id,
            start: block.start,
            base_length: block.base_length,
            offset,
            length: u64::try_from(encoded.len()).map_err(|_| JmaError::OffsetOverflow)?,
            checksum: checksum(encoded),
        });
    }
    if cursor != payload.len() {
        return Err(JmaError::CorruptSection(format!(
            "sequence section has {} trailing bytes",
            payload.len() - cursor
        )));
    }
    Ok(result)
}

fn scan_seed_section(
    payload: &[u8],
    descriptor: &SectionDescriptor,
) -> JmaResult<(Vec<SeedLevelIndex>, Vec<SeedBucketIndex>)> {
    let section = decode_seed_section(payload)?;
    let mut cursor =
        usize::try_from(SEED_SECTION_HEADER_SIZE).map_err(|_| JmaError::OffsetOverflow)?;
    let mut levels = Vec::with_capacity(section.levels.len());
    let mut buckets = Vec::new();
    for level in &section.levels {
        let level_header_start = cursor;
        let level_header_end = cursor
            .checked_add(SEED_LEVEL_HEADER_SIZE as usize)
            .ok_or(JmaError::OffsetOverflow)?;
        if level_header_end > payload.len() {
            return Err(JmaError::CorruptSection(
                "seed level header is truncated".to_string(),
            ));
        }
        let scale = get_u64(payload, cursor)?;
        let count = get_u64(payload, cursor + 8)?;
        if scale != level.level.scale {
            return Err(JmaError::CorruptSection(
                "sidecar seed level scale does not match payload".to_string(),
            ));
        }
        cursor = level_header_end;
        let record_length = count
            .checked_mul(SEED_RECORD_SIZE)
            .ok_or(JmaError::OffsetOverflow)?;
        let record_length_usize = checked_usize(record_length, "seed record range")?;
        let record_end = cursor
            .checked_add(record_length_usize)
            .ok_or(JmaError::OffsetOverflow)?;
        if record_end > payload.len() {
            return Err(JmaError::CorruptSection(
                "seed record range is truncated".to_string(),
            ));
        }
        let records = &payload[cursor..record_end];
        let record_offset = descriptor
            .offset
            .checked_add(u64::try_from(cursor).map_err(|_| JmaError::OffsetOverflow)?)
            .ok_or(JmaError::OffsetOverflow)?;
        let section_offset = descriptor.offset;
        let level_index = SeedLevelIndex {
            k: section.k,
            scale,
            section_offset,
            section_length: descriptor.length,
            level_header_offset: descriptor
                .offset
                .checked_add(
                    u64::try_from(level_header_start).map_err(|_| JmaError::OffsetOverflow)?,
                )
                .ok_or(JmaError::OffsetOverflow)?,
            record_offset,
            record_count: count,
            record_length,
            checksum: checksum(records),
        };
        levels.push(level_index);

        let mut record_index = 0u64;
        while record_index < count {
            let first = checked_usize(
                record_index
                    .checked_mul(SEED_RECORD_SIZE)
                    .ok_or(JmaError::OffsetOverflow)?,
                "seed bucket offset",
            )?;
            let first_hash = get_u64(records, first)?;
            let prefix = (first_hash >> (64 - u32::from(HASH_PREFIX_BITS))) as u16;
            let mut end_index = record_index
                .checked_add(1)
                .ok_or(JmaError::OffsetOverflow)?;
            while end_index < count {
                let offset = checked_usize(
                    end_index
                        .checked_mul(SEED_RECORD_SIZE)
                        .ok_or(JmaError::OffsetOverflow)?,
                    "seed bucket offset",
                )?;
                let hash = get_u64(records, offset)?;
                let next_prefix = (hash >> (64 - u32::from(HASH_PREFIX_BITS))) as u16;
                if next_prefix != prefix {
                    break;
                }
                end_index += 1;
            }
            let bucket_record_count = end_index - record_index;
            let bucket_offset_in_level = record_index
                .checked_mul(SEED_RECORD_SIZE)
                .ok_or(JmaError::OffsetOverflow)?;
            let bucket_length = bucket_record_count
                .checked_mul(SEED_RECORD_SIZE)
                .ok_or(JmaError::OffsetOverflow)?;
            let bucket_offset = record_offset
                .checked_add(bucket_offset_in_level)
                .ok_or(JmaError::OffsetOverflow)?;
            let bucket_start = checked_usize(bucket_offset_in_level, "seed bucket start")?;
            let bucket_end = bucket_start
                .checked_add(checked_usize(bucket_length, "seed bucket length")?)
                .ok_or(JmaError::OffsetOverflow)?;
            buckets.push(SeedBucketIndex {
                k: section.k,
                scale,
                hash_prefix: prefix,
                offset: bucket_offset,
                length: bucket_length,
                record_count: bucket_record_count,
                checksum: checksum(&records[bucket_start..bucket_end]),
            });
            record_index = end_index;
        }
        cursor = record_end;
    }
    if cursor != payload.len() {
        return Err(JmaError::CorruptSection(format!(
            "seed section has {} trailing bytes",
            payload.len() - cursor
        )));
    }
    Ok((levels, buckets))
}

fn validate_sequence_index(
    blocks: &[SequenceBlockIndex],
    descriptor: &SectionDescriptor,
    contigs: &[ContigMetadata],
) -> JmaResult<()> {
    let first_offset = descriptor
        .offset
        .checked_add(SEQUENCE_SECTION_HEADER_SIZE)
        .ok_or(JmaError::OffsetOverflow)?;
    let section_end = checked_end(descriptor.offset, descriptor.length)?;
    let mut expected_offset = first_offset;
    let mut block_index = 0usize;
    for contig in contigs {
        let mut expected_base = 0u64;
        while block_index < blocks.len() && blocks[block_index].contig_id == contig.id {
            let block = &blocks[block_index];
            if block.offset != expected_offset || block.start != expected_base {
                return Err(JmaError::CorruptSection(format!(
                    "sequence sidecar has a gap or overlap at contig {}:{}",
                    contig.id, expected_base
                )));
            }
            if block.length < SEQUENCE_BLOCK_HEADER_SIZE {
                return Err(JmaError::CorruptSection(
                    "sequence sidecar block is shorter than its header".to_string(),
                ));
            }
            if block.base_length == 0 {
                return Err(JmaError::CorruptSection(
                    "sequence sidecar contains a zero-length block".to_string(),
                ));
            }
            if block.length > MAX_INDEXED_RANGE_BYTES {
                return Err(JmaError::CorruptSection(
                    "sequence sidecar block exceeds the indexed read limit".to_string(),
                ));
            }
            let end = checked_end(block.offset, block.length)?;
            if end > section_end {
                return Err(JmaError::CorruptSection(
                    "sequence sidecar block exceeds sequence section".to_string(),
                ));
            }
            expected_offset = end;
            expected_base = expected_base
                .checked_add(block.base_length)
                .ok_or(JmaError::OffsetOverflow)?;
            if expected_base > contig.length {
                return Err(JmaError::CorruptSection(
                    "sequence sidecar block exceeds contig length".to_string(),
                ));
            }
            block_index += 1;
        }
        if expected_base != contig.length {
            return Err(JmaError::CorruptSection(format!(
                "sequence sidecar covers {expected_base} of contig {}, expected {}",
                contig.id, contig.length
            )));
        }
    }
    if block_index != blocks.len() || expected_offset != section_end {
        return Err(JmaError::CorruptSection(
            "sequence sidecar does not cover the complete sequence section".to_string(),
        ));
    }
    Ok(())
}

fn validate_seed_index(
    levels: &[SeedLevelIndex],
    buckets: &[SeedBucketIndex],
    sections: &[SectionDescriptor],
    archive: &ArchiveHeader,
) -> JmaResult<()> {
    if levels
        .iter()
        .map(|level| SeedLevel {
            k: level.k,
            scale: level.scale,
        })
        .collect::<Vec<_>>()
        != archive.seed_levels
    {
        return Err(JmaError::CorruptSection(
            "JMA sidecar seed levels do not match archive header".to_string(),
        ));
    }
    let mut sorted_levels = levels.to_vec();
    sorted_levels.sort_by_key(|level| (level.section_offset, level.level_header_offset));
    for level in &sorted_levels {
        let descriptor = sections
            .iter()
            .find(|section| section.offset == level.section_offset)
            .ok_or_else(|| {
                JmaError::CorruptSection("sidecar seed section is missing".to_string())
            })?;
        if descriptor.kind != section_kind(level.k)? {
            return Err(JmaError::CorruptSection(
                "sidecar seed level k does not match its section".to_string(),
            ));
        }
        if level.section_length != descriptor.length {
            return Err(JmaError::CorruptSection(
                "sidecar seed section length does not match archive".to_string(),
            ));
        }
        let section_end = checked_end(level.section_offset, level.section_length)?;
        let level_header_end = checked_end(level.level_header_offset, SEED_LEVEL_HEADER_SIZE)?;
        let expected_record_length = level
            .record_count
            .checked_mul(SEED_RECORD_SIZE)
            .ok_or(JmaError::OffsetOverflow)?;
        if level_header_end > section_end
            || level.record_offset != level_header_end
            || checked_end(level.record_offset, level.record_length)? > section_end
            || level.record_length != expected_record_length
        {
            return Err(JmaError::CorruptSection(
                "sidecar seed level range is invalid".to_string(),
            ));
        }
    }
    for descriptor in sections
        .iter()
        .filter(|section| matches!(section.kind, SECTION_SEEDS_K21 | SECTION_SEEDS_K31))
    {
        let mut section_levels = sorted_levels
            .iter()
            .filter(|level| level.section_offset == descriptor.offset)
            .collect::<Vec<_>>();
        section_levels.sort_by_key(|level| level.level_header_offset);
        let mut expected = descriptor
            .offset
            .checked_add(SEED_SECTION_HEADER_SIZE)
            .ok_or(JmaError::OffsetOverflow)?;
        for level in section_levels {
            if level.level_header_offset != expected {
                return Err(JmaError::CorruptSection(
                    "sidecar seed levels are not contiguous".to_string(),
                ));
            }
            expected = checked_end(level.record_offset, level.record_length)?;
        }
        if expected != checked_end(descriptor.offset, descriptor.length)? {
            return Err(JmaError::CorruptSection(
                "sidecar seed levels do not cover their section".to_string(),
            ));
        }
    }

    for level in levels {
        let mut level_buckets = buckets
            .iter()
            .filter(|bucket| bucket.k == level.k && bucket.scale == level.scale)
            .collect::<Vec<_>>();
        level_buckets.sort_by_key(|bucket| bucket.offset);
        let level_end = checked_end(level.record_offset, level.record_length)?;
        let mut expected_offset = level.record_offset;
        let mut previous_prefix = None;
        for bucket in level_buckets {
            let expected_length = bucket
                .record_count
                .checked_mul(SEED_RECORD_SIZE)
                .ok_or(JmaError::OffsetOverflow)?;
            if bucket.length != expected_length || bucket.length == 0 {
                return Err(JmaError::CorruptSection(
                    "sidecar seed bucket length is invalid".to_string(),
                ));
            }
            if bucket.length > MAX_INDEXED_RANGE_BYTES {
                return Err(JmaError::CorruptSection(
                    "sidecar seed bucket exceeds the indexed read limit".to_string(),
                ));
            }
            if bucket.offset != expected_offset {
                return Err(JmaError::CorruptSection(
                    "sidecar seed buckets are not contiguous".to_string(),
                ));
            }
            if previous_prefix.is_some_and(|previous| bucket.hash_prefix <= previous) {
                return Err(JmaError::CorruptSection(
                    "sidecar seed bucket prefixes are not strictly increasing".to_string(),
                ));
            }
            let bucket_end = checked_end(bucket.offset, bucket.length)?;
            if bucket_end > level_end {
                return Err(JmaError::CorruptSection(
                    "sidecar seed bucket exceeds its level".to_string(),
                ));
            }
            previous_prefix = Some(bucket.hash_prefix);
            expected_offset = bucket_end;
        }
        if expected_offset != level_end
            && (level.record_count != 0 || expected_offset != level.record_offset)
        {
            return Err(JmaError::CorruptSection(
                "sidecar seed buckets do not cover their level".to_string(),
            ));
        }
    }
    Ok(())
}

fn section_kind(k: u8) -> JmaResult<u32> {
    match k {
        21 => Ok(SECTION_SEEDS_K21),
        31 => Ok(SECTION_SEEDS_K31),
        other => Err(JmaError::CorruptSection(format!(
            "unsupported JMA seed k={other}"
        ))),
    }
}

fn unique_section(sections: &[SectionDescriptor], kind: u32) -> JmaResult<&SectionDescriptor> {
    let mut matches = sections.iter().filter(|section| section.kind == kind);
    let first = matches.next().ok_or_else(|| {
        JmaError::CorruptSection(format!("required JMA section kind {kind} is missing"))
    })?;
    if matches.next().is_some() {
        return Err(JmaError::CorruptSection(format!(
            "JMA section kind {kind} occurs more than once"
        )));
    }
    Ok(first)
}

fn validate_section_kinds(sections: &[SectionDescriptor]) -> JmaResult<()> {
    for section in sections {
        if !matches!(
            section.kind,
            SECTION_CONTIGS | SECTION_SEQUENCES | SECTION_SEEDS_K21 | SECTION_SEEDS_K31
        ) {
            return Err(JmaError::CorruptSection(format!(
                "unsupported JMA section kind {}",
                section.kind
            )));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jma::sequence::pack_bases;
    use crate::jma::writer::{
        ArchiveParts, SeedLevelData, SeedRecord, SeedSection, encode_archive,
    };
    use crate::jma::{ContigMetadata, SeedLevel, SeedOccurrence, SeedQuery};
    use std::path::PathBuf;

    fn archive() -> Vec<u8> {
        let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let (packed, mask) = pack_bases(sequence);
        encode_archive(&ArchiveParts {
            flags: 0,
            source_sha256: [3; 32],
            contigs: vec![ContigMetadata {
                id: 0,
                name: "c".to_string(),
                length: sequence.len() as u64,
            }],
            sequence_blocks: vec![crate::jma::sequence::SequenceBlock {
                contig_id: 0,
                start: 0,
                base_length: sequence.len() as u64,
                packed,
                unknown_mask: mask,
            }],
            seed_sections: vec![SeedSection {
                k: 21,
                levels: vec![SeedLevelData {
                    level: SeedLevel { k: 21, scale: 1 },
                    records: vec![SeedRecord {
                        query: SeedQuery {
                            k: 21,
                            hash: u64::MAX / 2,
                            canonical_kmer: 1,
                        },
                        occurrence: SeedOccurrence {
                            contig_id: 0,
                            position: 0,
                            reverse: false,
                        },
                    }],
                }],
            }],
        })
        .unwrap()
    }

    #[test]
    fn sidecar_roundtrips_and_binds_header() {
        let bytes = archive();
        let index = build_index(&bytes).unwrap();
        let encoded = encode_index(&index).unwrap();
        assert_eq!(encoded, encode_index(&index).unwrap());
        let decoded = decode_index(&encoded).unwrap();
        assert_eq!(decoded, index);
        assert_eq!(decoded.archive_sha256, checksum(&bytes));
        let (_, parsed, sections, contigs) = parse_archive_layout(&bytes).unwrap();
        validate_against_archive(
            &decoded,
            bytes.len() as u64,
            &bytes[..HEADER_SIZE],
            &parsed.archive,
            &sections,
            &contigs,
        )
        .unwrap();
        let mut changed = bytes.clone();
        changed[16] ^= 1;
        let changed_header = parse_header(&changed[..HEADER_SIZE]).unwrap_err();
        assert!(changed_header.to_string().contains("checksum"));
    }

    #[test]
    fn sidecar_path_is_only_an_extension() {
        assert_eq!(
            sidecar_path("results/archive.jma"),
            PathBuf::from("results/archive.jma.idx.json")
        );
    }
}
