//! Explicit, checked encoding primitives for JMA version 1.
//!
//! The archive format is deliberately boring: all integer fields are encoded
//! little-endian and all variable-length values are preceded by a checked
//! length. This module contains no `unsafe` layout casts, so the on-disk
//! representation is independent of the compiler and target platform.

use crate::jma::{ArchiveHeader, JMA_FORMAT_VERSION, JMA_MAGIC, JmaError, JmaResult, SeedLevel};
use sha2::{Digest, Sha256};

/// Number of bytes in the fixed JMA v1 header.
pub const HEADER_SIZE: usize = 160;
/// Number of bytes in one section-directory entry.
pub const SECTION_ENTRY_SIZE: usize = 64;
/// Version of each variable-length section payload.
pub const SECTION_VERSION: u32 = 1;

/// Section containing the contig table and names.
pub const SECTION_CONTIGS: u32 = 1;
/// Section containing packed sequence blocks.
pub const SECTION_SEQUENCES: u32 = 2;
/// Section containing k=31 seed occurrences.
pub const SECTION_SEEDS_K31: u32 = 3;
/// Section containing k=21 seed occurrences.
pub const SECTION_SEEDS_K21: u32 = 4;

/// A checked section-directory entry.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SectionDescriptor {
    pub kind: u32,
    pub flags: u32,
    pub offset: u64,
    pub length: u64,
    pub uncompressed_length: u64,
    pub checksum: [u8; 32],
}

/// Header fields that are not part of the public shared `ArchiveHeader`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParsedHeader {
    pub archive: ArchiveHeader,
    pub section_count: u32,
    pub section_directory_offset: u64,
    pub section_directory_length: u64,
}

/// Computes the SHA-256 checksum used for headers and section payloads.
#[must_use]
pub fn checksum(bytes: &[u8]) -> [u8; 32] {
    let digest = Sha256::digest(bytes);
    let mut result = [0u8; 32];
    result.copy_from_slice(&digest);
    result
}

/// Returns an error before converting a potentially large on-disk value to
/// `usize`.
pub fn checked_usize(value: u64, what: &str) -> JmaResult<usize> {
    usize::try_from(value)
        .map_err(|_| JmaError::CorruptSection(format!("{what} does not fit in usize: {value}")))
}

/// Returns the checked end of a half-open byte range.
pub fn checked_end(offset: u64, length: u64) -> JmaResult<u64> {
    offset.checked_add(length).ok_or(JmaError::OffsetOverflow)
}

/// Encodes the fixed header. The two 16-byte seed-level slots are enough for
/// the frozen k=31 primary and k=21 rescue levels in JMA v1.
pub fn encode_header(
    archive: &ArchiveHeader,
    section_count: u32,
    section_directory_offset: u64,
    section_directory_length: u64,
) -> JmaResult<Vec<u8>> {
    if archive.format_version != JMA_FORMAT_VERSION {
        return Err(JmaError::UnsupportedVersion(archive.format_version));
    }
    if archive.seed_levels.len() > 2 {
        return Err(JmaError::CorruptSection(
            "JMA v1 supports at most two seed levels".to_string(),
        ));
    }
    checked_end(section_directory_offset, section_directory_length)?;

    let mut bytes = vec![0u8; HEADER_SIZE];
    bytes[0..4].copy_from_slice(&JMA_MAGIC);
    put_u16(&mut bytes, 4, archive.format_version);
    put_u16(
        &mut bytes,
        6,
        u16::try_from(HEADER_SIZE).expect("fixed JMA header fits u16"),
    );
    put_u32(&mut bytes, 8, archive.flags);
    put_u32(&mut bytes, 12, archive.contig_count);
    put_u64(&mut bytes, 16, archive.total_bases);
    put_u32(&mut bytes, 24, section_count);
    put_u32(
        &mut bytes,
        28,
        u32::try_from(archive.seed_levels.len()).map_err(|_| JmaError::OffsetOverflow)?,
    );
    bytes[32..64].copy_from_slice(&archive.source_sha256);
    put_u64(&mut bytes, 64, section_directory_offset);
    put_u64(&mut bytes, 72, section_directory_length);

    for (index, level) in archive.seed_levels.iter().enumerate() {
        let offset = 112 + index * 16;
        bytes[offset] = level.k;
        put_u64(&mut bytes, offset + 8, level.scale);
    }

    // Header checksum covers all fixed fields and reserved bytes. The
    // checksum field itself is zero while calculating its value.
    let digest = checksum(&bytes);
    bytes[80..112].copy_from_slice(&digest);
    Ok(bytes)
}

/// Parses and validates the fixed header, including its checksum.
pub fn decode_header(bytes: &[u8]) -> JmaResult<ParsedHeader> {
    if bytes.len() < HEADER_SIZE {
        return Err(JmaError::CorruptSection(format!(
            "header is truncated: {} bytes, expected {HEADER_SIZE}",
            bytes.len()
        )));
    }
    if bytes[0] != JMA_MAGIC[0]
        || bytes[1] != JMA_MAGIC[1]
        || bytes[2] != JMA_MAGIC[2]
        || bytes[3] != JMA_MAGIC[3]
    {
        return Err(JmaError::InvalidMagic);
    }
    let version = get_u16(bytes, 4)?;
    if version != JMA_FORMAT_VERSION {
        return Err(JmaError::UnsupportedVersion(version));
    }
    let header_length = get_u16(bytes, 6)? as usize;
    if header_length != HEADER_SIZE {
        return Err(JmaError::CorruptSection(format!(
            "unsupported header length {header_length}, expected {HEADER_SIZE}"
        )));
    }

    let expected_checksum = array_32(bytes, 80)?;
    let mut checksum_input = bytes[..HEADER_SIZE].to_vec();
    checksum_input[80..112].fill(0);
    if checksum(&checksum_input) != expected_checksum {
        return Err(JmaError::ChecksumMismatch("header".to_string()));
    }

    let contig_count = get_u32(bytes, 12)?;
    let seed_count = get_u32(bytes, 28)?;
    if seed_count > 2 {
        return Err(JmaError::CorruptSection(format!(
            "header declares unsupported seed-level count {seed_count}"
        )));
    }
    let mut seed_levels = Vec::with_capacity(seed_count as usize);
    for index in 0..seed_count as usize {
        let offset = 112 + index * 16;
        let k = bytes[offset];
        let scale = get_u64(bytes, offset + 8)?;
        if !(1..=32).contains(&k) || scale == 0 {
            return Err(JmaError::CorruptSection(format!(
                "invalid seed level k={k}, scale={scale}"
            )));
        }
        seed_levels.push(SeedLevel { k, scale });
    }

    let section_directory_offset = get_u64(bytes, 64)?;
    let section_directory_length = get_u64(bytes, 72)?;
    checked_end(section_directory_offset, section_directory_length)?;
    let archive = ArchiveHeader {
        format_version: version,
        flags: get_u32(bytes, 8)?,
        contig_count,
        total_bases: get_u64(bytes, 16)?,
        source_sha256: array_32(bytes, 32)?,
        seed_levels,
    };
    Ok(ParsedHeader {
        archive,
        section_count: get_u32(bytes, 24)?,
        section_directory_offset,
        section_directory_length,
    })
}

/// Encodes a section directory and checks every offset/length arithmetic.
pub fn encode_section_directory(entries: &[SectionDescriptor]) -> JmaResult<Vec<u8>> {
    let size = entries
        .len()
        .checked_mul(SECTION_ENTRY_SIZE)
        .ok_or(JmaError::OffsetOverflow)?;
    let mut bytes = vec![0u8; size];
    for (index, entry) in entries.iter().enumerate() {
        checked_end(entry.offset, entry.length)?;
        checked_usize(entry.length, "section length")?;
        let offset = index * SECTION_ENTRY_SIZE;
        put_u32(&mut bytes, offset, entry.kind);
        put_u32(&mut bytes, offset + 4, entry.flags);
        put_u64(&mut bytes, offset + 8, entry.offset);
        put_u64(&mut bytes, offset + 16, entry.length);
        put_u64(&mut bytes, offset + 24, entry.uncompressed_length);
        bytes[offset + 32..offset + 64].copy_from_slice(&entry.checksum);
    }
    Ok(bytes)
}

/// Parses an exact section directory. Extra or truncated entries are both
/// rejected so a corrupt count cannot shift subsequent data interpretation.
pub fn decode_section_directory(
    bytes: &[u8],
    section_count: u32,
) -> JmaResult<Vec<SectionDescriptor>> {
    let expected = (section_count as usize)
        .checked_mul(SECTION_ENTRY_SIZE)
        .ok_or(JmaError::OffsetOverflow)?;
    if bytes.len() != expected {
        return Err(JmaError::CorruptSection(format!(
            "section directory has {} bytes, expected {expected}",
            bytes.len()
        )));
    }
    let mut entries = Vec::with_capacity(section_count as usize);
    for index in 0..section_count as usize {
        let offset = index * SECTION_ENTRY_SIZE;
        let entry = SectionDescriptor {
            kind: get_u32(bytes, offset)?,
            flags: get_u32(bytes, offset + 4)?,
            offset: get_u64(bytes, offset + 8)?,
            length: get_u64(bytes, offset + 16)?,
            uncompressed_length: get_u64(bytes, offset + 24)?,
            checksum: array_32_at(bytes, offset + 32)?,
        };
        checked_end(entry.offset, entry.length)?;
        if entry.uncompressed_length < entry.length {
            return Err(JmaError::CorruptSection(format!(
                "section {index} has invalid uncompressed length {} < {}",
                entry.uncompressed_length, entry.length
            )));
        }
        entries.push(entry);
    }
    Ok(entries)
}

fn put_u16(bytes: &mut [u8], offset: usize, value: u16) {
    bytes[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
}

fn put_u32(bytes: &mut [u8], offset: usize, value: u32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

fn put_u64(bytes: &mut [u8], offset: usize, value: u64) {
    bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

pub(crate) fn get_u16(bytes: &[u8], offset: usize) -> JmaResult<u16> {
    let value = bytes
        .get(offset..offset + 2)
        .ok_or_else(|| JmaError::CorruptSection("truncated u16".to_string()))?;
    Ok(u16::from_le_bytes([value[0], value[1]]))
}

pub(crate) fn get_u32(bytes: &[u8], offset: usize) -> JmaResult<u32> {
    let value = bytes
        .get(offset..offset + 4)
        .ok_or_else(|| JmaError::CorruptSection("truncated u32".to_string()))?;
    Ok(u32::from_le_bytes([value[0], value[1], value[2], value[3]]))
}

pub(crate) fn get_u64(bytes: &[u8], offset: usize) -> JmaResult<u64> {
    let value = bytes
        .get(offset..offset + 8)
        .ok_or_else(|| JmaError::CorruptSection("truncated u64".to_string()))?;
    Ok(u64::from_le_bytes([
        value[0], value[1], value[2], value[3], value[4], value[5], value[6], value[7],
    ]))
}

pub(crate) fn array_32(bytes: &[u8], offset: usize) -> JmaResult<[u8; 32]> {
    array_32_at(bytes, offset)
}

pub(crate) fn array_32_at(bytes: &[u8], offset: usize) -> JmaResult<[u8; 32]> {
    let value = bytes
        .get(offset..offset + 32)
        .ok_or_else(|| JmaError::CorruptSection("truncated 32-byte value".to_string()))?;
    let mut result = [0u8; 32];
    result.copy_from_slice(value);
    Ok(result)
}

pub(crate) fn checked_slice<'a>(
    bytes: &'a [u8],
    offset: usize,
    length: u64,
    what: &str,
) -> JmaResult<&'a [u8]> {
    let length = checked_usize(length, what)?;
    let end = offset.checked_add(length).ok_or(JmaError::OffsetOverflow)?;
    bytes
        .get(offset..end)
        .ok_or_else(|| JmaError::CorruptSection(format!("truncated {what}")))
}
