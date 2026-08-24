//! The immutable JAM Witness Router (JWR) format 1.
//!
//! The format deliberately has no dependency on Rust layout or allocator
//! details.  All integers are little endian and every offset is checked
//! before it is converted to a platform `usize`.  Payload sections are kept
//! separate so a range reader can fetch only the directory, filter, key
//! block, or posting bytes needed for a lookup.

use crate::router::{HashAlgorithmId, WitnessScheme};
use sha2::{Digest, Sha256};
use std::collections::BTreeSet;
use std::fmt;
use thiserror::Error;

/// Eight bytes identifying a JAM Witness Router object.
pub const ROUTER_FORMAT_MAGIC: [u8; 8] = *b"JWRF1\0\0\0";
/// The immutable router format version.
pub const ROUTER_FORMAT_VERSION: u16 = 1;
/// A fixed first request for a remote reader.
pub const SUPERBLOCK_SIZE: usize = 256;
/// The section directory uses fixed-width records.
pub const SECTION_ENTRY_SIZE: usize = 64;
/// Current on-disk layout identifier.
pub const LAYOUT_IDENTIFIER: u32 = 0x3152_574a;

pub const SECTION_METAGENOMES: u32 = 1;
pub const SECTION_SCHEMES: u32 = 2;
pub const SECTION_FILTER: u32 = 3;
pub const SECTION_PREFIX: u32 = 4;
pub const SECTION_KEYS: u32 = 5;
pub const SECTION_POSTING_HEADERS: u32 = 6;
pub const SECTION_POSTING_PAYLOAD: u32 = 7;
pub const SECTION_POSITION_PAYLOAD: u32 = 8;

pub const SECTION_FLAG_REQUIRED: u32 = 1;

/// A checked half-open object range.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ObjectRange {
    pub offset: u64,
    pub length: u64,
}

impl ObjectRange {
    pub fn new(offset: u64, length: u64) -> Result<Self, RouterFormatError> {
        checked_end(offset, length)?;
        Ok(Self { offset, length })
    }

    pub fn end(self) -> Result<u64, RouterFormatError> {
        checked_end(self.offset, self.length)
    }
}

/// A section in the immutable object. Coordinates are absolute object bytes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SectionDescriptor {
    pub kind: u32,
    pub flags: u32,
    pub offset: u64,
    pub length: u64,
    pub checksum: [u8; 32],
}

impl SectionDescriptor {
    pub fn range(self) -> Result<ObjectRange, RouterFormatError> {
        ObjectRange::new(self.offset, self.length)
    }
}

/// Fixed fields validated before variable sections are read.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RouterHeader {
    pub flags: u32,
    pub section_count: u32,
    pub section_directory_offset: u64,
    pub section_directory_length: u64,
    pub object_size: u64,
    pub metagenome_count: u32,
    pub scheme_count: u32,
    pub catalog_checksum: [u8; 32],
    pub witness_scheme_checksum: [u8; 32],
    pub object_checksum: [u8; 32],
    pub section_directory_checksum: [u8; 32],
    pub superblock_checksum: [u8; 32],
}

/// The identity fields bound into the fixed superblock.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct RouterIdentity {
    pub catalog_checksum: [u8; 32],
    pub witness_scheme_checksum: [u8; 32],
    pub object_checksum: [u8; 32],
}

/// A metagenome row used by the collection router.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MetagenomeEntry {
    pub id: String,
    pub object_uri: String,
    pub checksum: [u8; 32],
}

/// The exact packed-key encoding used by the key section.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PackedKeyWidth {
    Six = 6,
    Eight = 8,
}

impl PackedKeyWidth {
    pub const fn bytes(self) -> usize {
        self as usize
    }

    pub fn from_bytes(value: u8) -> Result<Self, RouterFormatError> {
        match value {
            6 => Ok(Self::Six),
            8 => Ok(Self::Eight),
            other => Err(RouterFormatError::UnsupportedKeyWidth(other)),
        }
    }
}

/// Configuration for physical exact-key blocks.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct KeyBlockConfig {
    pub target_block_bytes: u32,
    pub packed_width: PackedKeyWidth,
    /// This is intentionally opt-in.  A hash is recomputed from each packed
    /// key by default, avoiding a redundant per-record jamhash field.
    pub store_jamhash: bool,
}

impl Default for KeyBlockConfig {
    fn default() -> Self {
        Self {
            target_block_bytes: 8 * 1024,
            packed_width: PackedKeyWidth::Six,
            store_jamhash: false,
        }
    }
}

impl KeyBlockConfig {
    pub fn validate(self) -> Result<(), RouterFormatError> {
        if !matches!(self.target_block_bytes, 4_096 | 8_192 | 16_384 | 32_768) {
            return Err(RouterFormatError::InvalidBlockTarget(
                self.target_block_bytes,
            ));
        }
        Ok(())
    }

    pub const fn encoded_record_bytes(self) -> usize {
        self.packed_width.bytes() + if self.store_jamhash { 8 } else { 0 }
    }
}

/// A sorted block in the exact-key section.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct KeyBlockDescriptor {
    pub block_id: u32,
    pub first_key_index: u64,
    pub key_count: u32,
    pub offset: u64,
    pub length: u64,
    pub first_hash: u64,
    pub last_hash: u64,
    pub checksum: [u8; 32],
}

/// One sparse prefix-directory entry. Multiple entries may share a prefix
/// when a very common prefix does not fit one target block.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PrefixEntry {
    pub prefix: u32,
    pub block_id: u32,
    pub first_key_index: u64,
    pub key_count: u32,
    pub offset: u64,
    pub length: u64,
    pub first_hash: u64,
    pub last_hash: u64,
}

/// Header for one key's raw posting and optional position payload.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PostingHeader {
    pub document_frequency: u32,
    pub posting_count: u32,
    pub flags: u32,
    pub posting_offset: u64,
    pub posting_length: u64,
    pub position_offset: u64,
    pub position_length: u64,
    pub checksum: [u8; 32],
}

impl PostingHeader {
    pub const FLAG_COMMON: u32 = 1;
    pub const FLAG_POSITION_BEARING: u32 = 1 << 1;
    pub const FLAG_SUPPRESSED: u32 = 1 << 2;

    pub const fn has_positions(self) -> bool {
        self.flags & Self::FLAG_POSITION_BEARING != 0 && self.position_length != 0
    }
}

/// Parsed section directory plus checked descriptors.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SectionDirectory {
    pub entries: Vec<SectionDescriptor>,
}

impl SectionDirectory {
    pub fn get(&self, kind: u32) -> Option<SectionDescriptor> {
        self.entries
            .iter()
            .find(|entry| entry.kind == kind)
            .copied()
    }

    pub fn required(&self, kind: u32) -> Result<SectionDescriptor, RouterFormatError> {
        self.get(kind)
            .ok_or(RouterFormatError::MissingSection(kind))
    }
}

#[derive(Debug, Error, Eq, PartialEq)]
pub enum RouterFormatError {
    #[error("router object is truncated: expected at least {expected} bytes, got {actual}")]
    Truncated { expected: usize, actual: usize },
    #[error("invalid router magic")]
    InvalidMagic,
    #[error("unsupported router format version {0}")]
    UnsupportedVersion(u16),
    #[error("unsupported router layout identifier 0x{0:08x}")]
    InvalidLayout(u32),
    #[error("router checksum mismatch in {0}")]
    ChecksumMismatch(String),
    #[error("router offset arithmetic overflow")]
    OffsetOverflow,
    #[error("router {what} does not fit in usize: {value}")]
    NotUsableAsUsize { what: &'static str, value: u64 },
    #[error(
        "router {what} range is outside object: offset={offset}, length={length}, object={object_size}"
    )]
    RangeOutOfBounds {
        what: &'static str,
        offset: u64,
        length: u64,
        object_size: u64,
    },
    #[error("router section directory is invalid: {0}")]
    InvalidDirectory(String),
    #[error("router section {0} is missing")]
    MissingSection(u32),
    #[error("router section {0} appears more than once")]
    DuplicateSection(u32),
    #[error("router section kind {0} is unknown and required")]
    UnknownRequiredSection(u32),
    #[error("router section checksum mismatch for kind {0}")]
    SectionChecksumMismatch(u32),
    #[error("router object identity checksum mismatch")]
    ObjectChecksumMismatch,
    #[error("router reserved bytes are non-zero in {0}")]
    NonZeroReserved(&'static str),
    #[error("router count is too large: {0}")]
    CountTooLarge(&'static str),
    #[error("invalid UTF-8 in router {0}")]
    InvalidUtf8(&'static str),
    #[error("duplicate metagenome identifier {0}")]
    DuplicateMetagenome(String),
    #[error("invalid packed key width {0}; expected 6 or 8")]
    UnsupportedKeyWidth(u8),
    #[error("invalid key-block target size {0}; expected 4096, 8192, 16384, or 32768")]
    InvalidBlockTarget(u32),
    #[error("key block is not sorted or has an invalid range")]
    InvalidKeyBlock,
    #[error("invalid prefix directory")]
    InvalidPrefixDirectory,
    #[error("invalid membership filter")]
    InvalidFilter,
    #[error("invalid scheme table")]
    InvalidSchemeTable,
    #[error("invalid metagenome table")]
    InvalidMetagenomeTable,
    #[error("invalid posting header table")]
    InvalidPostingHeaders,
}

#[must_use]
pub fn checksum(bytes: &[u8]) -> [u8; 32] {
    let digest = Sha256::digest(bytes);
    let mut output = [0u8; 32];
    output.copy_from_slice(&digest);
    output
}

/// The object identity digest excludes its own two recursive checksum fields
/// (object and superblock) but includes all sections and the directory.
#[must_use]
pub fn checksum_object(bytes: &[u8]) -> [u8; 32] {
    let mut copy = bytes.to_vec();
    if copy.len() >= SUPERBLOCK_SIZE {
        copy[OBJECT_CHECKSUM_OFFSET..OBJECT_CHECKSUM_OFFSET + 32].fill(0);
        copy[SUPERBLOCK_CHECKSUM_OFFSET..SUPERBLOCK_CHECKSUM_OFFSET + 32].fill(0);
    }
    checksum(&copy)
}

pub const OBJECT_CHECKSUM_OFFSET: usize = 128;
pub const SUPERBLOCK_CHECKSUM_OFFSET: usize = 160;
pub const DIRECTORY_CHECKSUM_OFFSET: usize = 192;
pub const CATALOG_CHECKSUM_OFFSET: usize = 64;
pub const SCHEME_CHECKSUM_OFFSET: usize = 96;

pub fn checked_end(offset: u64, length: u64) -> Result<u64, RouterFormatError> {
    offset
        .checked_add(length)
        .ok_or(RouterFormatError::OffsetOverflow)
}

pub fn checked_usize(value: u64, what: &'static str) -> Result<usize, RouterFormatError> {
    usize::try_from(value).map_err(|_| RouterFormatError::NotUsableAsUsize { what, value })
}

pub fn range_in_object(
    offset: u64,
    length: u64,
    object_size: u64,
    what: &'static str,
) -> Result<ObjectRange, RouterFormatError> {
    let end = checked_end(offset, length)?;
    if end > object_size {
        return Err(RouterFormatError::RangeOutOfBounds {
            what,
            offset,
            length,
            object_size,
        });
    }
    Ok(ObjectRange { offset, length })
}

fn is_known_section(kind: u32) -> bool {
    matches!(
        kind,
        SECTION_METAGENOMES
            | SECTION_SCHEMES
            | SECTION_FILTER
            | SECTION_PREFIX
            | SECTION_KEYS
            | SECTION_POSTING_HEADERS
            | SECTION_POSTING_PAYLOAD
            | SECTION_POSITION_PAYLOAD
    )
}

/// Encode the fixed section directory. The caller supplies absolute ranges.
pub fn encode_section_directory(
    entries: &[SectionDescriptor],
) -> Result<Vec<u8>, RouterFormatError> {
    let mut sorted = entries.to_vec();
    sorted.sort_by_key(|entry| entry.kind);
    let mut seen = BTreeSet::new();
    let mut bytes = Vec::with_capacity(sorted.len() * SECTION_ENTRY_SIZE);
    for entry in sorted {
        if !seen.insert(entry.kind) {
            return Err(RouterFormatError::DuplicateSection(entry.kind));
        }
        if !is_known_section(entry.kind) && entry.flags & SECTION_FLAG_REQUIRED != 0 {
            return Err(RouterFormatError::UnknownRequiredSection(entry.kind));
        }
        bytes.extend_from_slice(&entry.kind.to_le_bytes());
        bytes.extend_from_slice(&entry.flags.to_le_bytes());
        bytes.extend_from_slice(&entry.offset.to_le_bytes());
        bytes.extend_from_slice(&entry.length.to_le_bytes());
        bytes.extend_from_slice(&entry.checksum);
        bytes.extend_from_slice(&0u64.to_le_bytes());
        debug_assert_eq!(bytes.len() % SECTION_ENTRY_SIZE, 0);
    }
    Ok(bytes)
}

/// Decode and validate a section directory against the declared object.
pub fn decode_section_directory(
    bytes: &[u8],
    section_count: u32,
    object_size: u64,
    directory_offset: u64,
    directory_length: u64,
) -> Result<SectionDirectory, RouterFormatError> {
    let expected_length = u64::from(section_count)
        .checked_mul(u64::try_from(SECTION_ENTRY_SIZE).expect("constant fits u64"))
        .ok_or(RouterFormatError::OffsetOverflow)?;
    if directory_length != expected_length
        || bytes.len() != checked_usize(directory_length, "section directory")?
    {
        return Err(RouterFormatError::InvalidDirectory(
            "directory length does not match section count".to_string(),
        ));
    }
    range_in_object(
        directory_offset,
        directory_length,
        object_size,
        "section directory",
    )?;
    let mut entries = Vec::with_capacity(usize::try_from(section_count).unwrap_or(0));
    let mut seen = BTreeSet::new();
    for index in 0..usize::try_from(section_count)
        .map_err(|_| RouterFormatError::CountTooLarge("section count"))?
    {
        let start = index
            .checked_mul(SECTION_ENTRY_SIZE)
            .ok_or(RouterFormatError::OffsetOverflow)?;
        let end = start
            .checked_add(SECTION_ENTRY_SIZE)
            .ok_or(RouterFormatError::OffsetOverflow)?;
        let entry = &bytes[start..end];
        let kind = get_u32(entry, 0)?;
        let flags = get_u32(entry, 4)?;
        if !seen.insert(kind) {
            return Err(RouterFormatError::DuplicateSection(kind));
        }
        if !is_known_section(kind) && flags & SECTION_FLAG_REQUIRED != 0 {
            return Err(RouterFormatError::UnknownRequiredSection(kind));
        }
        if entry[56..64].iter().any(|byte| *byte != 0) {
            return Err(RouterFormatError::NonZeroReserved("section directory"));
        }
        let offset = get_u64(entry, 8)?;
        let length = get_u64(entry, 16)?;
        range_in_object(offset, length, object_size, "section")?;
        let checksum = array_32(entry, 24)?;
        entries.push(SectionDescriptor {
            kind,
            flags,
            offset,
            length,
            checksum,
        });
    }
    // Section payloads and the directory must be disjoint, and sections may
    // not overlap one another. Empty sections are allowed but still bounded.
    let directory_end = checked_end(directory_offset, directory_length)?;
    let mut ranges = entries
        .iter()
        .map(|entry| (entry.offset, entry.length, entry.kind))
        .collect::<Vec<_>>();
    ranges.sort_unstable_by_key(|range| range.0);
    let mut previous_end = u64::try_from(SUPERBLOCK_SIZE).expect("constant fits u64");
    for (offset, length, kind) in ranges {
        let end = checked_end(offset, length)?;
        if offset < previous_end || (offset < directory_end && end > directory_offset) {
            return Err(RouterFormatError::InvalidDirectory(format!(
                "section {kind} overlaps an earlier range or the directory"
            )));
        }
        previous_end = end;
    }
    if directory_offset < u64::try_from(SUPERBLOCK_SIZE).expect("constant fits u64") {
        return Err(RouterFormatError::InvalidDirectory(
            "directory overlaps superblock".to_string(),
        ));
    }
    Ok(SectionDirectory { entries })
}

/// Encode the fixed superblock. The object and superblock checksums are
/// filled by [`finalize_object`].
#[allow(clippy::too_many_arguments)]
pub fn encode_header(
    flags: u32,
    section_count: u32,
    section_directory_offset: u64,
    section_directory_length: u64,
    object_size: u64,
    metagenome_count: u32,
    scheme_count: u32,
    catalog_checksum: [u8; 32],
    witness_scheme_checksum: [u8; 32],
    section_directory_checksum: [u8; 32],
) -> Result<Vec<u8>, RouterFormatError> {
    if object_size < u64::try_from(SUPERBLOCK_SIZE).expect("constant fits u64")
        || checked_end(section_directory_offset, section_directory_length)? > object_size
    {
        return Err(RouterFormatError::InvalidDirectory(
            "directory exceeds declared object size".to_string(),
        ));
    }
    let mut bytes = vec![0u8; SUPERBLOCK_SIZE];
    bytes[0..8].copy_from_slice(&ROUTER_FORMAT_MAGIC);
    put_u16(&mut bytes, 8, ROUTER_FORMAT_VERSION);
    put_u16(
        &mut bytes,
        10,
        u16::try_from(SUPERBLOCK_SIZE).expect("constant fits u16"),
    );
    put_u32(&mut bytes, 12, LAYOUT_IDENTIFIER);
    put_u32(&mut bytes, 16, flags);
    put_u32(&mut bytes, 20, section_count);
    put_u64(&mut bytes, 24, section_directory_offset);
    put_u64(&mut bytes, 32, section_directory_length);
    put_u64(&mut bytes, 40, object_size);
    put_u32(&mut bytes, 48, metagenome_count);
    put_u32(&mut bytes, 52, scheme_count);
    bytes[CATALOG_CHECKSUM_OFFSET..SCHEME_CHECKSUM_OFFSET].copy_from_slice(&catalog_checksum);
    bytes[SCHEME_CHECKSUM_OFFSET..OBJECT_CHECKSUM_OFFSET].copy_from_slice(&witness_scheme_checksum);
    bytes[DIRECTORY_CHECKSUM_OFFSET..DIRECTORY_CHECKSUM_OFFSET + 32]
        .copy_from_slice(&section_directory_checksum);
    // 128..160: object identity, filled by finalize_object.
    // 160..192: superblock checksum, filled below.
    let mut check = bytes.clone();
    check[SUPERBLOCK_CHECKSUM_OFFSET..SUPERBLOCK_CHECKSUM_OFFSET + 32].fill(0);
    let digest = checksum(&check);
    bytes[SUPERBLOCK_CHECKSUM_OFFSET..SUPERBLOCK_CHECKSUM_OFFSET + 32].copy_from_slice(&digest);
    Ok(bytes)
}

/// Finish an object after all section bytes and the directory are laid out.
/// The input must contain a valid fixed superblock at byte zero.
pub fn finalize_object(
    mut object: Vec<u8>,
) -> Result<(Vec<u8>, RouterIdentity), RouterFormatError> {
    if object.len() < SUPERBLOCK_SIZE {
        return Err(RouterFormatError::Truncated {
            expected: SUPERBLOCK_SIZE,
            actual: object.len(),
        });
    }
    let object_size = get_u64(&object, 40)?;
    if object_size != u64::try_from(object.len()).map_err(|_| RouterFormatError::OffsetOverflow)? {
        return Err(RouterFormatError::InvalidDirectory(
            "object length does not match superblock".to_string(),
        ));
    }
    let identity_checksum = checksum_object(&object);
    object[OBJECT_CHECKSUM_OFFSET..OBJECT_CHECKSUM_OFFSET + 32].copy_from_slice(&identity_checksum);
    let mut check = object[..SUPERBLOCK_SIZE].to_vec();
    check[SUPERBLOCK_CHECKSUM_OFFSET..SUPERBLOCK_CHECKSUM_OFFSET + 32].fill(0);
    let superblock_checksum = checksum(&check);
    object[SUPERBLOCK_CHECKSUM_OFFSET..SUPERBLOCK_CHECKSUM_OFFSET + 32]
        .copy_from_slice(&superblock_checksum);
    let catalog_checksum = array_32(&object, CATALOG_CHECKSUM_OFFSET)?;
    let witness_scheme_checksum = array_32(&object, SCHEME_CHECKSUM_OFFSET)?;
    Ok((
        object,
        RouterIdentity {
            catalog_checksum,
            witness_scheme_checksum,
            object_checksum: identity_checksum,
        },
    ))
}

/// Parse and validate the fixed superblock. Directory and object checksums are
/// checked separately because remote readers may intentionally defer payload
/// reads.
pub fn decode_header(bytes: &[u8]) -> Result<RouterHeader, RouterFormatError> {
    if bytes.len() < SUPERBLOCK_SIZE {
        return Err(RouterFormatError::Truncated {
            expected: SUPERBLOCK_SIZE,
            actual: bytes.len(),
        });
    }
    if bytes[0..8] != ROUTER_FORMAT_MAGIC {
        return Err(RouterFormatError::InvalidMagic);
    }
    let version = get_u16(bytes, 8)?;
    if version != ROUTER_FORMAT_VERSION {
        return Err(RouterFormatError::UnsupportedVersion(version));
    }
    let header_length = usize::from(get_u16(bytes, 10)?);
    if header_length != SUPERBLOCK_SIZE {
        return Err(RouterFormatError::Truncated {
            expected: SUPERBLOCK_SIZE,
            actual: header_length,
        });
    }
    let layout = get_u32(bytes, 12)?;
    if layout != LAYOUT_IDENTIFIER {
        return Err(RouterFormatError::InvalidLayout(layout));
    }
    let expected = array_32(bytes, SUPERBLOCK_CHECKSUM_OFFSET)?;
    let mut check = bytes[..SUPERBLOCK_SIZE].to_vec();
    check[SUPERBLOCK_CHECKSUM_OFFSET..SUPERBLOCK_CHECKSUM_OFFSET + 32].fill(0);
    if checksum(&check) != expected {
        return Err(RouterFormatError::ChecksumMismatch(
            "superblock".to_string(),
        ));
    }
    let object_size = get_u64(bytes, 40)?;
    if object_size < u64::try_from(SUPERBLOCK_SIZE).expect("constant fits u64") {
        return Err(RouterFormatError::InvalidDirectory(
            "object is smaller than the superblock".to_string(),
        ));
    }
    let directory_offset = get_u64(bytes, 24)?;
    let directory_length = get_u64(bytes, 32)?;
    range_in_object(
        directory_offset,
        directory_length,
        object_size,
        "section directory",
    )?;
    Ok(RouterHeader {
        flags: get_u32(bytes, 16)?,
        section_count: get_u32(bytes, 20)?,
        section_directory_offset: directory_offset,
        section_directory_length: directory_length,
        object_size,
        metagenome_count: get_u32(bytes, 48)?,
        scheme_count: get_u32(bytes, 52)?,
        catalog_checksum: array_32(bytes, CATALOG_CHECKSUM_OFFSET)?,
        witness_scheme_checksum: array_32(bytes, SCHEME_CHECKSUM_OFFSET)?,
        object_checksum: array_32(bytes, OBJECT_CHECKSUM_OFFSET)?,
        section_directory_checksum: array_32(bytes, DIRECTORY_CHECKSUM_OFFSET)?,
        superblock_checksum: expected,
    })
}

pub fn verify_object_checksum(
    object: &[u8],
    header: &RouterHeader,
) -> Result<(), RouterFormatError> {
    if u64::try_from(object.len()).ok() != Some(header.object_size) {
        return Err(RouterFormatError::RangeOutOfBounds {
            what: "complete object",
            offset: 0,
            length: u64::try_from(object.len()).unwrap_or(u64::MAX),
            object_size: header.object_size,
        });
    }
    if checksum_object(object) != header.object_checksum {
        return Err(RouterFormatError::ObjectChecksumMismatch);
    }
    Ok(())
}

pub fn verify_section_checksum(
    bytes: &[u8],
    descriptor: SectionDescriptor,
) -> Result<(), RouterFormatError> {
    if checksum(bytes) != descriptor.checksum {
        return Err(RouterFormatError::SectionChecksumMismatch(descriptor.kind));
    }
    Ok(())
}

/// Encode a metagenome table. Offsets are section-relative.
pub fn encode_metagenomes(entries: &[MetagenomeEntry]) -> Result<Vec<u8>, RouterFormatError> {
    const HEADER: usize = 16;
    const RECORD: usize = 64;
    let count = u32::try_from(entries.len())
        .map_err(|_| RouterFormatError::CountTooLarge("metagenomes"))?;
    let records_end = HEADER
        .checked_add(
            entries
                .len()
                .checked_mul(RECORD)
                .ok_or(RouterFormatError::OffsetOverflow)?,
        )
        .ok_or(RouterFormatError::OffsetOverflow)?;
    let mut bytes = vec![0u8; records_end];
    put_u32(&mut bytes, 0, 1);
    put_u32(&mut bytes, 4, count);
    put_u64(
        &mut bytes,
        8,
        u64::try_from(records_end).map_err(|_| RouterFormatError::OffsetOverflow)?,
    );
    let mut strings = Vec::new();
    let mut ids = BTreeSet::new();
    for (index, entry) in entries.iter().enumerate() {
        if !ids.insert(entry.id.clone()) {
            return Err(RouterFormatError::DuplicateMetagenome(entry.id.clone()));
        }
        if entry.id.is_empty() || entry.object_uri.is_empty() {
            return Err(RouterFormatError::InvalidMetagenomeTable);
        }
        let id_offset = u64::try_from(records_end + strings.len())
            .map_err(|_| RouterFormatError::OffsetOverflow)?;
        let id_length = u32::try_from(entry.id.len())
            .map_err(|_| RouterFormatError::CountTooLarge("metagenome identifier"))?;
        strings.extend_from_slice(entry.id.as_bytes());
        let uri_offset = u64::try_from(records_end + strings.len())
            .map_err(|_| RouterFormatError::OffsetOverflow)?;
        let uri_length = u32::try_from(entry.object_uri.len())
            .map_err(|_| RouterFormatError::CountTooLarge("metagenome URI"))?;
        strings.extend_from_slice(entry.object_uri.as_bytes());
        let start = HEADER + index * RECORD;
        put_u64(&mut bytes, start, id_offset);
        put_u32(&mut bytes, start + 8, id_length);
        put_u64(&mut bytes, start + 16, uri_offset);
        put_u32(&mut bytes, start + 24, uri_length);
        bytes[start + 32..start + 64].copy_from_slice(&entry.checksum);
    }
    bytes.extend_from_slice(&strings);
    Ok(bytes)
}

pub fn decode_metagenomes(bytes: &[u8]) -> Result<Vec<MetagenomeEntry>, RouterFormatError> {
    const HEADER: usize = 16;
    const RECORD: usize = 64;
    if bytes.len() < HEADER || get_u32(bytes, 0)? != 1 {
        return Err(RouterFormatError::InvalidMetagenomeTable);
    }
    let count = usize::try_from(get_u32(bytes, 4)?)
        .map_err(|_| RouterFormatError::CountTooLarge("metagenomes"))?;
    let records_end = HEADER
        .checked_add(
            count
                .checked_mul(RECORD)
                .ok_or(RouterFormatError::OffsetOverflow)?,
        )
        .ok_or(RouterFormatError::OffsetOverflow)?;
    let strings_offset = checked_usize(get_u64(bytes, 8)?, "metagenome strings")?;
    if strings_offset != records_end || strings_offset > bytes.len() {
        return Err(RouterFormatError::InvalidMetagenomeTable);
    }
    let mut result = Vec::with_capacity(count);
    let mut ids = BTreeSet::new();
    for index in 0..count {
        let start = HEADER + index * RECORD;
        if start + RECORD > bytes.len()
            || bytes[start + 12..start + 16].iter().any(|byte| *byte != 0)
            || bytes[start + 28..start + 32].iter().any(|byte| *byte != 0)
        {
            return Err(RouterFormatError::InvalidMetagenomeTable);
        }
        let id_offset = get_u64(bytes, start)?;
        let uri_offset = get_u64(bytes, start + 16)?;
        if id_offset
            < u64::try_from(strings_offset).map_err(|_| RouterFormatError::OffsetOverflow)?
            || uri_offset
                < u64::try_from(strings_offset).map_err(|_| RouterFormatError::OffsetOverflow)?
        {
            return Err(RouterFormatError::InvalidMetagenomeTable);
        }
        let id = read_string(
            bytes,
            id_offset,
            get_u32(bytes, start + 8)?,
            "metagenome identifier",
        )?;
        let uri = read_string(
            bytes,
            uri_offset,
            get_u32(bytes, start + 24)?,
            "metagenome URI",
        )?;
        if !ids.insert(id.clone()) {
            return Err(RouterFormatError::DuplicateMetagenome(id));
        }
        result.push(MetagenomeEntry {
            id,
            object_uri: uri,
            checksum: array_32(bytes, start + 32)?,
        });
    }
    Ok(result)
}

/// Encode the shared nested witness scheme descriptors.
pub fn encode_schemes(schemes: &[WitnessScheme]) -> Result<Vec<u8>, RouterFormatError> {
    const HEADER: usize = 16;
    const RECORD: usize = 48;
    let count =
        u32::try_from(schemes.len()).map_err(|_| RouterFormatError::CountTooLarge("schemes"))?;
    let records_end = HEADER
        .checked_add(
            schemes
                .len()
                .checked_mul(RECORD)
                .ok_or(RouterFormatError::OffsetOverflow)?,
        )
        .ok_or(RouterFormatError::OffsetOverflow)?;
    let mut bytes = vec![0u8; records_end];
    put_u32(&mut bytes, 0, 1);
    put_u32(&mut bytes, 4, count);
    put_u64(
        &mut bytes,
        8,
        u64::try_from(records_end).map_err(|_| RouterFormatError::OffsetOverflow)?,
    );
    let mut scales = Vec::new();
    let mut ids = BTreeSet::new();
    for (index, scheme) in schemes.iter().enumerate() {
        scheme
            .validate()
            .map_err(|_| RouterFormatError::InvalidSchemeTable)?;
        if !ids.insert(scheme.scheme_id) {
            return Err(RouterFormatError::InvalidSchemeTable);
        }
        let scales_offset = u64::try_from(records_end + scales.len())
            .map_err(|_| RouterFormatError::OffsetOverflow)?;
        let scales_length = u32::try_from(scheme.available_scales.len())
            .map_err(|_| RouterFormatError::CountTooLarge("scheme scales"))?;
        for scale in &scheme.available_scales {
            scales.extend_from_slice(&scale.to_le_bytes());
        }
        let start = HEADER + index * RECORD;
        put_u32(&mut bytes, start, scheme.scheme_id);
        bytes[start + 4] = scheme.k;
        bytes[start + 5] = u8::from(scheme.zero_excluded);
        bytes[start + 6] = match scheme.hash_id {
            HashAlgorithmId::JamhashU64V1 => 1,
        };
        put_u32(&mut bytes, start + 8, scheme.base_scale);
        put_u32(&mut bytes, start + 12, scales_length);
        put_u64(&mut bytes, start + 16, scales_offset);
    }
    bytes.extend_from_slice(&scales);
    Ok(bytes)
}

pub fn decode_schemes(bytes: &[u8]) -> Result<Vec<WitnessScheme>, RouterFormatError> {
    const HEADER: usize = 16;
    const RECORD: usize = 48;
    if bytes.len() < HEADER || get_u32(bytes, 0)? != 1 {
        return Err(RouterFormatError::InvalidSchemeTable);
    }
    let count = usize::try_from(get_u32(bytes, 4)?)
        .map_err(|_| RouterFormatError::CountTooLarge("schemes"))?;
    let records_end = HEADER
        .checked_add(
            count
                .checked_mul(RECORD)
                .ok_or(RouterFormatError::OffsetOverflow)?,
        )
        .ok_or(RouterFormatError::OffsetOverflow)?;
    let scales_offset = checked_usize(get_u64(bytes, 8)?, "scheme scales")?;
    if scales_offset != records_end || scales_offset > bytes.len() {
        return Err(RouterFormatError::InvalidSchemeTable);
    }
    let mut result = Vec::with_capacity(count);
    let mut ids = BTreeSet::new();
    for index in 0..count {
        let start = HEADER + index * RECORD;
        if start + RECORD > bytes.len()
            || bytes[start + 7] != 0
            || bytes[start + 24..start + RECORD]
                .iter()
                .any(|byte| *byte != 0)
        {
            return Err(RouterFormatError::InvalidSchemeTable);
        }
        let scheme_id = get_u32(bytes, start)?;
        if !ids.insert(scheme_id) {
            return Err(RouterFormatError::InvalidSchemeTable);
        }
        let k = bytes[start + 4];
        let zero_excluded = bytes[start + 5];
        if zero_excluded > 1 || bytes[start + 6] != 1 {
            return Err(RouterFormatError::InvalidSchemeTable);
        }
        let base_scale = get_u32(bytes, start + 8)?;
        let scale_count = usize::try_from(get_u32(bytes, start + 12)?)
            .map_err(|_| RouterFormatError::CountTooLarge("scheme scales"))?;
        let scale_offset = checked_usize(get_u64(bytes, start + 16)?, "scheme scale offset")?;
        let scale_end = scale_offset
            .checked_add(
                scale_count
                    .checked_mul(4)
                    .ok_or(RouterFormatError::OffsetOverflow)?,
            )
            .ok_or(RouterFormatError::OffsetOverflow)?;
        if scale_offset < scales_offset || scale_end > bytes.len() {
            return Err(RouterFormatError::InvalidSchemeTable);
        }
        let mut available_scales = Vec::with_capacity(scale_count);
        for offset in (scale_offset..scale_end).step_by(4) {
            available_scales.push(get_u32(bytes, offset)?);
        }
        let scheme = WitnessScheme {
            scheme_id,
            k,
            base_scale,
            available_scales,
            hash_id: HashAlgorithmId::JamhashU64V1,
            zero_excluded: zero_excluded == 1,
        };
        scheme
            .validate()
            .map_err(|_| RouterFormatError::InvalidSchemeTable)?;
        result.push(scheme);
    }
    Ok(result)
}

fn read_string(
    bytes: &[u8],
    offset: u64,
    length: u32,
    what: &'static str,
) -> Result<String, RouterFormatError> {
    let offset = checked_usize(offset, what)?;
    let length = usize::try_from(length).map_err(|_| RouterFormatError::CountTooLarge(what))?;
    let end = offset
        .checked_add(length)
        .ok_or(RouterFormatError::OffsetOverflow)?;
    let value = bytes
        .get(offset..end)
        .ok_or(RouterFormatError::InvalidMetagenomeTable)?;
    String::from_utf8(value.to_vec()).map_err(|_| RouterFormatError::InvalidUtf8(what))
}

/// A deterministic blocked Bloom filter. False positives only permit an exact
/// key-page read; they never create witness evidence.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MembershipFilter {
    pub bit_count: u64,
    pub hash_count: u32,
    pub seed: u64,
    pub bits: Vec<u8>,
}

impl MembershipFilter {
    pub fn build(
        hashes: impl IntoIterator<Item = u64>,
        bits_per_key: u32,
    ) -> Result<Self, RouterFormatError> {
        if bits_per_key < 4 {
            return Err(RouterFormatError::InvalidFilter);
        }
        let hashes = hashes
            .into_iter()
            .filter(|hash| *hash != 0)
            .collect::<Vec<_>>();
        let bit_count = (u64::try_from(hashes.len())
            .map_err(|_| RouterFormatError::OffsetOverflow)?
            .saturating_mul(u64::from(bits_per_key)))
        .max(64)
        .div_ceil(64)
        .saturating_mul(64);
        let hash_count = (u64::from(bits_per_key) * 69 / 100).clamp(2, 8) as u32;
        let mut filter = Self {
            bit_count,
            hash_count,
            seed: 0x9e37_79b9_7f4a_7c15,
            bits: vec![0; checked_usize(bit_count / 8, "filter bits")?],
        };
        for hash in hashes {
            filter.insert(hash);
        }
        Ok(filter)
    }

    pub fn contains_hash(&self, hash: u64) -> bool {
        if hash == 0 || self.bit_count == 0 || self.hash_count == 0 {
            return false;
        }
        let (first, step) = filter_hashes(hash, self.seed);
        (0..self.hash_count).all(|index| {
            let bit = first.wrapping_add(u64::from(index).wrapping_mul(step)) % self.bit_count;
            self.bits[(bit / 8) as usize] & (1 << (bit % 8)) != 0
        })
    }

    fn insert(&mut self, hash: u64) {
        let (first, step) = filter_hashes(hash, self.seed);
        for index in 0..self.hash_count {
            let bit = first.wrapping_add(u64::from(index).wrapping_mul(step)) % self.bit_count;
            self.bits[(bit / 8) as usize] |= 1 << (bit % 8);
        }
    }
}

pub fn encode_filter(filter: &MembershipFilter) -> Result<Vec<u8>, RouterFormatError> {
    if filter.bit_count == 0
        || filter.hash_count == 0
        || !filter.bit_count.is_multiple_of(8)
        || checked_usize(filter.bit_count / 8, "filter bits")? != filter.bits.len()
    {
        return Err(RouterFormatError::InvalidFilter);
    }
    let mut bytes = vec![0u8; 32];
    put_u32(&mut bytes, 0, 1);
    bytes[4] = 1;
    put_u32(&mut bytes, 8, filter.hash_count);
    put_u64(&mut bytes, 12, filter.bit_count);
    put_u64(&mut bytes, 20, filter.seed);
    bytes.extend_from_slice(&filter.bits);
    Ok(bytes)
}

pub fn decode_filter(bytes: &[u8]) -> Result<MembershipFilter, RouterFormatError> {
    if bytes.len() < 32
        || get_u32(bytes, 0)? != 1
        || bytes[4] != 1
        || bytes[5..8].iter().any(|byte| *byte != 0)
    {
        return Err(RouterFormatError::InvalidFilter);
    }
    let hash_count = get_u32(bytes, 8)?;
    let bit_count = get_u64(bytes, 12)?;
    let seed = get_u64(bytes, 20)?;
    if bit_count == 0
        || bit_count % 8 != 0
        || hash_count == 0
        || bytes[28..32].iter().any(|byte| *byte != 0)
    {
        return Err(RouterFormatError::InvalidFilter);
    }
    let bit_bytes = checked_usize(bit_count / 8, "filter bits")?;
    if bytes.len() != 32 + bit_bytes {
        return Err(RouterFormatError::InvalidFilter);
    }
    Ok(MembershipFilter {
        bit_count,
        hash_count,
        seed,
        bits: bytes[32..].to_vec(),
    })
}

fn filter_hashes(hash: u64, seed: u64) -> (u64, u64) {
    let first = splitmix64(hash ^ seed);
    let step = splitmix64(first ^ 0xa076_1d64_78bd_642f) | 1;
    (first, step)
}

fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

pub fn get_u16(bytes: &[u8], offset: usize) -> Result<u16, RouterFormatError> {
    let slice = bytes
        .get(offset..offset + 2)
        .ok_or(RouterFormatError::Truncated {
            expected: offset + 2,
            actual: bytes.len(),
        })?;
    Ok(u16::from_le_bytes([slice[0], slice[1]]))
}

pub fn get_u32(bytes: &[u8], offset: usize) -> Result<u32, RouterFormatError> {
    let slice = bytes
        .get(offset..offset + 4)
        .ok_or(RouterFormatError::Truncated {
            expected: offset + 4,
            actual: bytes.len(),
        })?;
    Ok(u32::from_le_bytes([slice[0], slice[1], slice[2], slice[3]]))
}

pub fn get_u64(bytes: &[u8], offset: usize) -> Result<u64, RouterFormatError> {
    let slice = bytes
        .get(offset..offset + 8)
        .ok_or(RouterFormatError::Truncated {
            expected: offset + 8,
            actual: bytes.len(),
        })?;
    Ok(u64::from_le_bytes([
        slice[0], slice[1], slice[2], slice[3], slice[4], slice[5], slice[6], slice[7],
    ]))
}

pub fn array_32(bytes: &[u8], offset: usize) -> Result<[u8; 32], RouterFormatError> {
    let slice = bytes
        .get(offset..offset + 32)
        .ok_or(RouterFormatError::Truncated {
            expected: offset + 32,
            actual: bytes.len(),
        })?;
    let mut output = [0u8; 32];
    output.copy_from_slice(slice);
    Ok(output)
}

pub fn put_u16(bytes: &mut [u8], offset: usize, value: u16) {
    bytes[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
}

pub fn put_u32(bytes: &mut [u8], offset: usize, value: u32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

pub fn put_u64(bytes: &mut [u8], offset: usize, value: u64) {
    bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

/// Hash prefix used by the sparse exact-key directory.
pub fn hash_prefix(hash: u64, bits: u8) -> Result<u32, RouterFormatError> {
    if bits > 32 {
        return Err(RouterFormatError::InvalidPrefixDirectory);
    }
    if bits == 0 {
        Ok(0)
    } else {
        Ok((hash >> (64 - u32::from(bits))) as u32)
    }
}

impl fmt::Display for PackedKeyWidth {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}-byte", self.bytes())
    }
}
