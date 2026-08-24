//! The on-disk JMA format 1 primitives.
//!
//! JMA is intentionally a binary, self-contained container. No Rust struct
//! layout is used as a disk contract: every field is written explicitly in
//! little endian order and every range is checked before conversion to a
//! platform `usize`.

use crate::jma::{ArchiveHeader, JMA_FORMAT_VERSION, JmaError, JmaResult};
use sha2::{Digest, Sha256};

/// Eight bytes that identify the format and prevent accidental decoding by
/// unrelated binary layouts.
pub const JMA_FORMAT_MAGIC: [u8; 8] = *b"JMAF1\0\0\0";
/// Layout identifier is independent of the format version.
pub const JMA_LAYOUT_IDENTIFIER: u32 = 0x3141_4d4a;

/// Fixed superblock size. It keeps the first remote request bounded.
pub const SUPERBLOCK_SIZE: usize = 256;
/// Named alias for callers that only need the fixed prefix length.
pub const HEADER_SIZE: usize = SUPERBLOCK_SIZE;
/// Fixed encoded size of a section-directory entry.
pub const SECTION_ENTRY_SIZE: usize = 64;
/// Version of variable-length section payloads.
pub const SECTION_VERSION: u32 = 1;

/// Unknown sections carrying this bit are fatal.
pub const SECTION_FLAG_REQUIRED: u32 = 1;
/// Payload is independently decoded by its codec.
pub const SECTION_FLAG_COMPRESSED: u32 = 1 << 1;

/// Required archive metadata (builder, source commit, checksums, and IDs).
pub const SECTION_METADATA: u32 = 1;
/// Contig records followed by the contig-name string table.
pub const SECTION_CONTIGS: u32 = 2;
/// Embedded sketch hashes used by candidate screening.
pub const SECTION_SKETCH: u32 = 3;
/// Open seed-scheme descriptors, page directories, keys, and occurrences.
pub const SECTION_SEEDS: u32 = 4;
/// Fixed sequence-block directory.
pub const SECTION_SEQUENCE_DIRECTORY: u32 = 5;
/// Independently addressable sequence block payloads.
pub const SECTION_SEQUENCE_PAYLOAD: u32 = 6;
/// Optional content-dependent gear block directory.
pub const SECTION_GEAR_DIRECTORY: u32 = 7;

/// Byte range in the superblock containing the checksum of the complete
/// section-directory bytes. It is covered by the superblock checksum.
pub const SECTION_DIRECTORY_CHECKSUM_OFFSET: usize = 160;

/// A checked section-directory entry. Coordinates are absolute object bytes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SectionDescriptor {
    pub kind: u32,
    pub flags: u32,
    pub offset: u64,
    pub length: u64,
    pub uncompressed_length: u64,
    pub checksum: [u8; 32],
}

/// Archive metadata stored in the required metadata section.
#[derive(Clone, Debug, PartialEq)]
pub struct ArchiveMetadataFields {
    pub format_identifier: String,
    pub format_version: u16,
    pub layout_identifier: u32,
    pub source_assembly_sha256: [u8; 32],
    pub archive_sha256: Option<[u8; 32]>,
    pub builder_version: String,
    pub source_commit: String,
    pub hash_algorithm: String,
    pub min_entropy: Option<f64>,
}

impl ArchiveMetadataFields {
    pub fn for_source(source: [u8; 32]) -> Self {
        Self {
            format_identifier: "JMA".to_string(),
            format_version: JMA_FORMAT_VERSION,
            layout_identifier: JMA_LAYOUT_IDENTIFIER,
            source_assembly_sha256: source,
            archive_sha256: None,
            builder_version: env!("CARGO_PKG_VERSION").to_string(),
            source_commit: option_env!("JAM_RS_SOURCE_COMMIT")
                .unwrap_or("unknown")
                .to_string(),
            hash_algorithm: "jamhash_u64_v1".to_string(),
            min_entropy: None,
        }
    }
}

/// Parsed critical fields from the superblock. Variable metadata is decoded
/// from `SECTION_METADATA` by the reader after this structure is checked.
#[derive(Clone, Debug, PartialEq)]
pub struct ParsedHeader {
    pub archive: ArchiveHeader,
    pub layout_identifier: u32,
    pub section_count: u32,
    pub section_directory_offset: u64,
    pub section_directory_length: u64,
    pub object_size: u64,
    pub archive_sha256: [u8; 32],
    pub section_directory_checksum: [u8; 32],
    pub superblock_checksum: [u8; 32],
}

#[must_use]
pub fn checksum(bytes: &[u8]) -> [u8; 32] {
    let digest = Sha256::digest(bytes);
    let mut result = [0u8; 32];
    result.copy_from_slice(&digest);
    result
}

/// Computes an object checksum without recursively including the object and
/// superblock checksum fields. The directory checksum remains authenticated
/// as part of the object bytes.
#[must_use]
pub fn checksum_object(bytes: &[u8]) -> [u8; 32] {
    if bytes.len() < SUPERBLOCK_SIZE {
        return checksum(bytes);
    }
    let mut digest = Sha256::new();
    digest.update(&bytes[..96]);
    digest.update([0; 64]);
    digest.update(&bytes[160..]);
    let mut output = [0u8; 32];
    output.copy_from_slice(&digest.finalize());
    output
}

pub fn checked_usize(value: u64, what: &str) -> JmaResult<usize> {
    usize::try_from(value)
        .map_err(|_| JmaError::CorruptSection(format!("{what} does not fit in usize: {value}")))
}

pub fn checked_end(offset: u64, length: u64) -> JmaResult<u64> {
    offset.checked_add(length).ok_or(JmaError::OffsetOverflow)
}

/// Encodes a format-1 superblock. The caller supplies final object size and
/// directory coordinates, so a remote reader can reject truncation early.
pub fn encode_header(
    archive: &ArchiveHeader,
    section_count: u32,
    section_directory_offset: u64,
    section_directory_length: u64,
    object_size: u64,
    section_directory_checksum: [u8; 32],
) -> JmaResult<Vec<u8>> {
    if archive.format_version != JMA_FORMAT_VERSION {
        return Err(JmaError::UnsupportedVersion(archive.format_version));
    }
    if object_size < u64::try_from(SUPERBLOCK_SIZE).map_err(|_| JmaError::OffsetOverflow)?
        || checked_end(section_directory_offset, section_directory_length)? > object_size
    {
        return Err(JmaError::CorruptSection(
            "section directory exceeds declared object size".to_string(),
        ));
    }
    let mut bytes = vec![0u8; SUPERBLOCK_SIZE];
    bytes[0..8].copy_from_slice(&JMA_FORMAT_MAGIC);
    put_u16(&mut bytes, 8, JMA_FORMAT_VERSION);
    put_u16(
        &mut bytes,
        10,
        u16::try_from(SUPERBLOCK_SIZE).expect("JMA superblock fits u16"),
    );
    put_u32(&mut bytes, 12, JMA_LAYOUT_IDENTIFIER);
    put_u32(&mut bytes, 16, archive.flags);
    put_u32(&mut bytes, 20, section_count);
    put_u64(&mut bytes, 24, section_directory_offset);
    put_u64(&mut bytes, 32, section_directory_length);
    put_u64(&mut bytes, 40, object_size);
    put_u32(&mut bytes, 48, archive.contig_count);
    // This is a count hint only. The open descriptor list is authoritative.
    put_u32(
        &mut bytes,
        52,
        u32::try_from(archive.seed_levels.len()).map_err(|_| JmaError::OffsetOverflow)?,
    );
    put_u64(&mut bytes, 56, archive.total_bases);
    bytes[64..96].copy_from_slice(&archive.source_sha256);
    // 96..128: object checksum, filled by the finalizer.
    // 128..160: superblock checksum, filled below.
    bytes[SECTION_DIRECTORY_CHECKSUM_OFFSET..SECTION_DIRECTORY_CHECKSUM_OFFSET + 32]
        .copy_from_slice(&section_directory_checksum);
    let digest = checksum(&bytes);
    bytes[128..160].copy_from_slice(&digest);
    Ok(bytes)
}

/// Parses and validates only the fixed superblock.
pub fn decode_header(bytes: &[u8]) -> JmaResult<ParsedHeader> {
    if bytes.len() < SUPERBLOCK_SIZE {
        return Err(JmaError::CorruptSection(format!(
            "superblock is truncated: {} bytes, expected {SUPERBLOCK_SIZE}",
            bytes.len()
        )));
    }
    if bytes[0..8] != JMA_FORMAT_MAGIC {
        return Err(JmaError::InvalidMagic);
    }
    let version = get_u16(bytes, 8)?;
    if version != JMA_FORMAT_VERSION {
        return Err(JmaError::UnsupportedVersion(version));
    }
    let header_length = usize::from(get_u16(bytes, 10)?);
    if header_length != SUPERBLOCK_SIZE {
        return Err(JmaError::CorruptSection(format!(
            "unsupported superblock length {header_length}, expected {SUPERBLOCK_SIZE}"
        )));
    }
    let layout_identifier = get_u32(bytes, 12)?;
    if layout_identifier != JMA_LAYOUT_IDENTIFIER {
        return Err(JmaError::CorruptSection(format!(
            "unsupported JMA layout identifier 0x{layout_identifier:08x}"
        )));
    }
    let expected = array_32(bytes, 128)?;
    let mut input = bytes[..SUPERBLOCK_SIZE].to_vec();
    input[128..160].fill(0);
    if checksum(&input) != expected {
        return Err(JmaError::ChecksumMismatch("superblock".to_string()));
    }
    let section_directory_offset = get_u64(bytes, 24)?;
    let section_directory_length = get_u64(bytes, 32)?;
    checked_end(section_directory_offset, section_directory_length)?;
    let object_size = get_u64(bytes, 40)?;
    if object_size < u64::try_from(SUPERBLOCK_SIZE).map_err(|_| JmaError::OffsetOverflow)? {
        return Err(JmaError::CorruptSection(
            "declared JMA object size is smaller than the superblock".to_string(),
        ));
    }
    if checked_end(section_directory_offset, section_directory_length)? > object_size {
        return Err(JmaError::CorruptSection(
            "section directory exceeds declared object size".to_string(),
        ));
    }
    Ok(ParsedHeader {
        archive: ArchiveHeader {
            format_version: version,
            flags: get_u32(bytes, 16)?,
            contig_count: get_u32(bytes, 48)?,
            total_bases: get_u64(bytes, 56)?,
            source_sha256: array_32(bytes, 64)?,
            // Filled from the open seed descriptor list by the reader.
            seed_levels: Vec::new(),
            algorithm_id: Some("jam-seed-chain-align-v1".to_string()),
            algorithm_version: Some(1),
            min_entropy: None,
        },
        layout_identifier,
        section_count: get_u32(bytes, 20)?,
        section_directory_offset,
        section_directory_length,
        object_size,
        archive_sha256: array_32(bytes, 96)?,
        section_directory_checksum: array_32(bytes, SECTION_DIRECTORY_CHECKSUM_OFFSET)?,
        superblock_checksum: expected,
    })
}

/// Encodes an exact section directory. Required kinds must be unique.
pub fn encode_section_directory(entries: &[SectionDescriptor]) -> JmaResult<Vec<u8>> {
    let size = entries
        .len()
        .checked_mul(SECTION_ENTRY_SIZE)
        .ok_or(JmaError::OffsetOverflow)?;
    let mut bytes = vec![0u8; size];
    let mut seen = std::collections::BTreeSet::new();
    for (index, entry) in entries.iter().enumerate() {
        if entry.flags & !(SECTION_FLAG_REQUIRED | SECTION_FLAG_COMPRESSED) != 0 {
            return Err(JmaError::CorruptSection(format!(
                "section kind {} has unsupported flags 0x{:08x}",
                entry.kind, entry.flags
            )));
        }
        if entry.flags & SECTION_FLAG_COMPRESSED == 0 && entry.uncompressed_length != entry.length {
            return Err(JmaError::CorruptSection(
                "uncompressed section has a decoded length different from stored length"
                    .to_string(),
            ));
        }
        checked_end(entry.offset, entry.length)?;
        checked_usize(entry.length, "section length")?;
        if entry.flags & SECTION_FLAG_REQUIRED != 0 && !seen.insert(entry.kind) {
            return Err(JmaError::CorruptSection(format!(
                "duplicate required section kind {}",
                entry.kind
            )));
        }
        let offset = index
            .checked_mul(SECTION_ENTRY_SIZE)
            .ok_or(JmaError::OffsetOverflow)?;
        put_u32(&mut bytes, offset, entry.kind);
        put_u32(&mut bytes, offset + 4, entry.flags);
        put_u64(&mut bytes, offset + 8, entry.offset);
        put_u64(&mut bytes, offset + 16, entry.length);
        put_u64(&mut bytes, offset + 24, entry.uncompressed_length);
        bytes[offset + 32..offset + 64].copy_from_slice(&entry.checksum);
    }
    Ok(bytes)
}

pub fn decode_section_directory(
    bytes: &[u8],
    section_count: u32,
) -> JmaResult<Vec<SectionDescriptor>> {
    let expected = usize::try_from(section_count)
        .ok()
        .and_then(|count| count.checked_mul(SECTION_ENTRY_SIZE))
        .ok_or(JmaError::OffsetOverflow)?;
    if bytes.len() != expected {
        return Err(JmaError::CorruptSection(format!(
            "section directory has {} bytes, expected {expected}",
            bytes.len()
        )));
    }
    let mut entries = Vec::with_capacity(usize::try_from(section_count).unwrap_or(0));
    let mut required = std::collections::BTreeSet::new();
    for index in 0..usize::try_from(section_count).map_err(|_| JmaError::OffsetOverflow)? {
        let offset = index
            .checked_mul(SECTION_ENTRY_SIZE)
            .ok_or(JmaError::OffsetOverflow)?;
        let entry = SectionDescriptor {
            kind: get_u32(bytes, offset)?,
            flags: get_u32(bytes, offset + 4)?,
            offset: get_u64(bytes, offset + 8)?,
            length: get_u64(bytes, offset + 16)?,
            uncompressed_length: get_u64(bytes, offset + 24)?,
            checksum: array_32_at(bytes, offset + 32)?,
        };
        checked_end(entry.offset, entry.length)?;
        if entry.flags & SECTION_FLAG_COMPRESSED == 0 && entry.uncompressed_length != entry.length {
            return Err(JmaError::CorruptSection(format!(
                "uncompressed section {index} has decoded length {} != stored length {}",
                entry.uncompressed_length, entry.length
            )));
        }
        if entry.flags & SECTION_FLAG_REQUIRED != 0 && !required.insert(entry.kind) {
            return Err(JmaError::CorruptSection(format!(
                "duplicate required section kind {}",
                entry.kind
            )));
        }
        entries.push(entry);
    }
    Ok(entries)
}

/// Encodes metadata. Strings are UTF-8 and have checked u32 lengths.
pub fn encode_metadata(metadata: &ArchiveMetadataFields) -> JmaResult<Vec<u8>> {
    if metadata.format_version != JMA_FORMAT_VERSION {
        return Err(JmaError::UnsupportedVersion(metadata.format_version));
    }
    if metadata.layout_identifier != JMA_LAYOUT_IDENTIFIER {
        return Err(JmaError::CorruptSection(
            "metadata layout identifier does not match JMA format 1".to_string(),
        ));
    }
    let fields = [
        metadata.format_identifier.as_bytes(),
        metadata.builder_version.as_bytes(),
        metadata.source_commit.as_bytes(),
        metadata.hash_algorithm.as_bytes(),
    ];
    let mut bytes = Vec::with_capacity(108 + fields.iter().map(|field| field.len()).sum::<usize>());
    bytes.extend_from_slice(&SECTION_VERSION.to_le_bytes());
    bytes.extend_from_slice(&metadata.format_version.to_le_bytes());
    bytes.extend_from_slice(&u16::from(metadata.min_entropy.is_some()).to_le_bytes());
    bytes.extend_from_slice(&metadata.layout_identifier.to_le_bytes());
    bytes.extend_from_slice(&metadata.source_assembly_sha256);
    bytes.extend_from_slice(&metadata.archive_sha256.unwrap_or([0; 32]));
    let entropy_bits = match metadata.min_entropy {
        None => 0,
        Some(value) if value.is_finite() && (0.0..=2.0).contains(&value) => value.to_bits(),
        Some(value) => {
            return Err(JmaError::CorruptSection(format!(
                "minimum entropy {value} is not finite and between 0.0 and 2.0"
            )));
        }
    };
    bytes.extend_from_slice(&entropy_bits.to_le_bytes());
    for field in fields {
        bytes.extend_from_slice(
            &u32::try_from(field.len())
                .map_err(|_| JmaError::OffsetOverflow)?
                .to_le_bytes(),
        );
    }
    for field in fields {
        bytes.extend_from_slice(field);
    }
    Ok(bytes)
}

pub fn decode_metadata(bytes: &[u8]) -> JmaResult<ArchiveMetadataFields> {
    // 4 version + 2 version + 2 reserved + 4 layout + 32 source + 32 object
    // + entropy bits + 4 string lengths.
    if bytes.len() < 100 {
        return Err(JmaError::CorruptSection(
            "metadata section header is truncated".to_string(),
        ));
    }
    if get_u32(bytes, 0)? != SECTION_VERSION {
        return Err(JmaError::CorruptSection(
            "unsupported metadata section version".to_string(),
        ));
    }
    let format_version = get_u16(bytes, 4)?;
    if format_version != JMA_FORMAT_VERSION {
        return Err(JmaError::UnsupportedVersion(format_version));
    }
    let layout_identifier = get_u32(bytes, 8)?;
    if layout_identifier != JMA_LAYOUT_IDENTIFIER {
        return Err(JmaError::CorruptSection(
            "metadata layout identifier does not match JMA format 1".to_string(),
        ));
    }
    let source = array_32(bytes, 12)?;
    let archive = array_32(bytes, 44)?;
    let reserved = get_u16(bytes, 6)?;
    if reserved & !1 != 0 {
        return Err(JmaError::CorruptSection(
            "metadata reserved flags are non-zero".to_string(),
        ));
    }
    let entropy_bits = get_u64(bytes, 76)?;
    let min_entropy = if reserved & 1 == 0 {
        if entropy_bits != 0 {
            return Err(JmaError::CorruptSection(
                "minimum entropy is present without a metadata flag".to_string(),
            ));
        }
        None
    } else {
        let value = f64::from_bits(entropy_bits);
        if !value.is_finite() || !(0.0..=2.0).contains(&value) {
            return Err(JmaError::CorruptSection(format!(
                "invalid JMA minimum entropy {value}"
            )));
        }
        Some(value)
    };
    let mut cursor = 84usize;
    let mut lengths = [0u32; 4];
    for length in &mut lengths {
        *length = get_u32(bytes, cursor)?;
        cursor = cursor.checked_add(4).ok_or(JmaError::OffsetOverflow)?;
    }
    let mut strings = Vec::with_capacity(4);
    for length in lengths {
        let end = cursor
            .checked_add(usize::try_from(length).map_err(|_| JmaError::OffsetOverflow)?)
            .ok_or(JmaError::OffsetOverflow)?;
        let value = bytes
            .get(cursor..end)
            .ok_or_else(|| JmaError::CorruptSection("metadata string is truncated".to_string()))?;
        let value = std::str::from_utf8(value)
            .map_err(|_| JmaError::CorruptSection("metadata string is not UTF-8".to_string()))?
            .to_string();
        strings.push(value);
        cursor = end;
    }
    if cursor != bytes.len() {
        return Err(JmaError::CorruptSection(format!(
            "metadata section has {} trailing bytes",
            bytes.len() - cursor
        )));
    }
    let archive_sha256 = (archive != [0; 32]).then_some(archive);
    Ok(ArchiveMetadataFields {
        format_identifier: strings.remove(0),
        format_version,
        layout_identifier,
        source_assembly_sha256: source,
        archive_sha256,
        builder_version: strings.remove(0),
        source_commit: strings.remove(0),
        hash_algorithm: strings.remove(0),
        min_entropy,
    })
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

/// Checks a packed k-mer before it is used as exact verification material.
pub fn valid_canonical_kmer(k: u8, value: u64) -> bool {
    k == 32 || (1..=32).contains(&k) && value < (1u64 << (2 * u32::from(k)))
}

/// Returns a high-order hash prefix for an embedded seed page.
pub fn hash_prefix(hash: u64, bucket_bits: u8) -> JmaResult<u32> {
    if !(1..=32).contains(&bucket_bits) {
        return Err(JmaError::CorruptSection(format!(
            "invalid seed bucket width {bucket_bits}"
        )));
    }
    u32::try_from(hash >> (64 - u32::from(bucket_bits))).map_err(|_| JmaError::OffsetOverflow)
}
