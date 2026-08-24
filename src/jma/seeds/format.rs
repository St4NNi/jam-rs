use super::JAMHASH_ALGORITHM_ID;
use crate::jma::{JmaError, JmaResult};

pub const COMPACT_SEED_MAGIC: [u8; 8] = *b"JMS1SEED";
pub const COMPACT_SEED_VERSION: u32 = 1;
pub const COMPACT_HEADER_SIZE: usize = 320;
pub const PREFIX_RECORD_SIZE: usize = 40;
pub const KEY_PAGE_RECORD_SIZE: usize = 84;
pub const POSTING_HEADER_RECORD_SIZE: usize = 64;
pub const POSITION_PAGE_RECORD_SIZE: usize = 64;
pub const MAX_BUCKET_BITS: u8 = 24;
pub const POSTING_FLAG_POSITION_BEARING: u16 = 1;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
pub enum KeyEncoding {
    Packed8 = 1,
    Packed6 = 2,
}

impl KeyEncoding {
    pub const fn record_size(self) -> usize {
        match self {
            Self::Packed8 => 8,
            Self::Packed6 => 6,
        }
    }

    pub const fn from_byte(value: u8) -> Option<Self> {
        match value {
            1 => Some(Self::Packed8),
            2 => Some(Self::Packed6),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
pub enum OccurrenceEncoding {
    Fixed16 = 1,
    DeltaVarint = 2,
}

impl OccurrenceEncoding {
    pub const fn from_byte(value: u8) -> Option<Self> {
        match value {
            1 => Some(Self::Fixed16),
            2 => Some(Self::DeltaVarint),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
pub enum PostingClass {
    Rare = 0,
    ModeratelyCommon = 1,
    CollectionCommon = 2,
    Suppressed = 3,
}

impl PostingClass {
    pub const fn from_byte(value: u8) -> Option<Self> {
        match value {
            0 => Some(Self::Rare),
            1 => Some(Self::ModeratelyCommon),
            2 => Some(Self::CollectionCommon),
            3 => Some(Self::Suppressed),
            _ => None,
        }
    }
}

/// Explicit archive/catalog binding. A zero checksum means that the
/// corresponding binding is intentionally left open by the builder.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedBinding {
    pub scheme_id: u32,
    pub k: u8,
    pub scale: u32,
    pub catalog_checksum: [u8; 32],
    pub archive_checksum: [u8; 32],
    pub scheme_checksum: [u8; 32],
}

impl SeedBinding {
    pub fn validate(&self) -> JmaResult<()> {
        if !(21..=32).contains(&self.k) {
            return Err(JmaError::CorruptSection(format!(
                "seed k must be between 21 and 32, got {}",
                self.k
            )));
        }
        if self.scale == 0 {
            return Err(JmaError::CorruptSection(
                "seed scale must be greater than zero".to_string(),
            ));
        }
        Ok(())
    }

    pub fn matches(&self, header: &SeedHeader) -> JmaResult<()> {
        self.validate()?;
        if self.scheme_id != header.scheme_id || self.k != header.k || self.scale != header.scale {
            return Err(JmaError::CorruptSection(
                "seed scheme binding does not match the section".to_string(),
            ));
        }
        if self.catalog_checksum != [0; 32] && self.catalog_checksum != header.catalog_checksum {
            return Err(JmaError::ChecksumMismatch(
                "seed catalog binding".to_string(),
            ));
        }
        if self.archive_checksum != [0; 32] && self.archive_checksum != header.archive_checksum {
            return Err(JmaError::ChecksumMismatch(
                "seed archive binding".to_string(),
            ));
        }
        if self.scheme_checksum != [0; 32] && self.scheme_checksum != header.scheme_checksum {
            return Err(JmaError::ChecksumMismatch(
                "seed scheme binding".to_string(),
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedHeader {
    pub scheme_id: u32,
    pub k: u8,
    pub key_encoding: KeyEncoding,
    pub occurrence_encoding: OccurrenceEncoding,
    pub bucket_bits: u8,
    pub scale: u32,
    pub hash_algorithm_id: u32,
    pub key_count: u64,
    pub posting_count: u64,
    pub filter_offset: u64,
    pub filter_length: u64,
    pub prefix_offset: u64,
    pub prefix_count: u32,
    pub key_page_dir_offset: u64,
    pub key_page_count: u32,
    pub posting_header_offset: u64,
    pub posting_header_count: u32,
    pub position_page_dir_offset: u64,
    pub position_page_count: u32,
    pub key_data_offset: u64,
    pub key_data_length: u64,
    pub position_data_offset: u64,
    pub position_data_length: u64,
    pub section_length: u64,
    pub flags: u64,
    pub catalog_checksum: [u8; 32],
    pub archive_checksum: [u8; 32],
    pub scheme_checksum: [u8; 32],
    pub section_checksum: [u8; 32],
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedPrefix {
    pub hash_prefix: u32,
    pub first_key: u64,
    pub key_count: u32,
    pub first_page: u32,
    pub page_count: u32,
    pub first_hash: u64,
    pub last_hash: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedKeyPage {
    pub page_id: u32,
    pub first_key: u64,
    pub key_count: u32,
    pub first_hash: u64,
    pub last_hash: u64,
    pub data_offset: u64,
    pub data_length: u64,
    pub checksum: [u8; 32],
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedPostingHeader {
    pub key_index: u64,
    pub occurrence_count: u32,
    pub class: PostingClass,
    pub encoding: OccurrenceEncoding,
    pub flags: u16,
    pub position_page_id: u32,
    pub position_offset: u32,
    pub position_length: u32,
    pub checksum: [u8; 32],
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedPositionPage {
    pub page_id: u32,
    pub first_header: u32,
    pub header_count: u32,
    pub data_offset: u64,
    pub data_length: u64,
    pub checksum: [u8; 32],
}

pub(crate) fn put_u32(bytes: &mut [u8], offset: usize, value: u32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

pub(crate) fn put_u64(bytes: &mut [u8], offset: usize, value: u64) {
    bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

pub(crate) fn get_u8(bytes: &[u8], offset: usize) -> JmaResult<u8> {
    bytes
        .get(offset)
        .copied()
        .ok_or_else(|| JmaError::CorruptSection("compact seed value is truncated".to_string()))
}

pub(crate) fn get_u32(bytes: &[u8], offset: usize) -> JmaResult<u32> {
    let value = bytes
        .get(offset..offset + 4)
        .ok_or_else(|| JmaError::CorruptSection("compact seed u32 is truncated".to_string()))?;
    Ok(u32::from_le_bytes(
        value.try_into().expect("checked length"),
    ))
}

pub(crate) fn get_u16(bytes: &[u8], offset: usize) -> JmaResult<u16> {
    let value = bytes
        .get(offset..offset + 2)
        .ok_or_else(|| JmaError::CorruptSection("compact seed u16 is truncated".to_string()))?;
    Ok(u16::from_le_bytes(
        value.try_into().expect("checked length"),
    ))
}

pub(crate) fn get_u64(bytes: &[u8], offset: usize) -> JmaResult<u64> {
    let value = bytes
        .get(offset..offset + 8)
        .ok_or_else(|| JmaError::CorruptSection("compact seed u64 is truncated".to_string()))?;
    Ok(u64::from_le_bytes(
        value.try_into().expect("checked length"),
    ))
}

pub(crate) fn get_array(bytes: &[u8], offset: usize) -> JmaResult<[u8; 32]> {
    let value = bytes.get(offset..offset + 32).ok_or_else(|| {
        JmaError::CorruptSection("compact seed checksum is truncated".to_string())
    })?;
    Ok(value.try_into().expect("checked length"))
}

pub(crate) fn section_range(
    offset: u64,
    length: u64,
    section_length: u64,
    bytes_length: usize,
    what: &str,
) -> JmaResult<std::ops::Range<usize>> {
    let end = offset.checked_add(length).ok_or(JmaError::OffsetOverflow)?;
    if end > section_length {
        return Err(JmaError::CorruptSection(format!(
            "compact seed {what} lies outside the section"
        )));
    }
    let start = usize::try_from(offset).map_err(|_| JmaError::OffsetOverflow)?;
    let end = usize::try_from(end).map_err(|_| JmaError::OffsetOverflow)?;
    if end > bytes_length {
        return Err(JmaError::CorruptSection(format!(
            "compact seed {what} is not present in the supplied range"
        )));
    }
    Ok(start..end)
}

pub(crate) fn encode_header(header: SeedHeader) -> Vec<u8> {
    let mut bytes = vec![0u8; COMPACT_HEADER_SIZE];
    bytes[0..8].copy_from_slice(&COMPACT_SEED_MAGIC);
    put_u32(&mut bytes, 8, COMPACT_SEED_VERSION);
    put_u32(
        &mut bytes,
        12,
        u32::try_from(COMPACT_HEADER_SIZE).expect("header fits u32"),
    );
    put_u32(&mut bytes, 16, header.scheme_id);
    bytes[20] = header.k;
    bytes[21] = header.key_encoding as u8;
    bytes[22] = header.occurrence_encoding as u8;
    bytes[23] = header.bucket_bits;
    put_u32(&mut bytes, 24, header.scale);
    put_u32(&mut bytes, 28, header.hash_algorithm_id);
    put_u64(&mut bytes, 32, header.key_count);
    put_u64(&mut bytes, 40, header.posting_count);
    put_u64(&mut bytes, 48, header.filter_offset);
    put_u64(&mut bytes, 56, header.filter_length);
    put_u64(&mut bytes, 64, header.prefix_offset);
    put_u32(&mut bytes, 72, header.prefix_count);
    put_u32(&mut bytes, 76, PREFIX_RECORD_SIZE as u32);
    put_u64(&mut bytes, 80, header.key_page_dir_offset);
    put_u32(&mut bytes, 88, header.key_page_count);
    put_u32(&mut bytes, 92, KEY_PAGE_RECORD_SIZE as u32);
    put_u64(&mut bytes, 96, header.posting_header_offset);
    put_u32(&mut bytes, 104, header.posting_header_count);
    put_u32(&mut bytes, 108, POSTING_HEADER_RECORD_SIZE as u32);
    put_u64(&mut bytes, 112, header.position_page_dir_offset);
    put_u32(&mut bytes, 120, header.position_page_count);
    put_u32(&mut bytes, 124, POSITION_PAGE_RECORD_SIZE as u32);
    put_u64(&mut bytes, 128, header.key_data_offset);
    put_u64(&mut bytes, 136, header.key_data_length);
    put_u64(&mut bytes, 144, header.position_data_offset);
    put_u64(&mut bytes, 152, header.position_data_length);
    put_u64(&mut bytes, 160, header.section_length);
    put_u64(&mut bytes, 168, header.flags);
    bytes[192..224].copy_from_slice(&header.catalog_checksum);
    bytes[224..256].copy_from_slice(&header.archive_checksum);
    bytes[256..288].copy_from_slice(&header.scheme_checksum);
    bytes[288..320].copy_from_slice(&header.section_checksum);
    bytes
}

pub(crate) fn decode_header(bytes: &[u8], section_length: u64) -> JmaResult<SeedHeader> {
    if bytes.len() < COMPACT_HEADER_SIZE {
        return Err(JmaError::CorruptSection(
            "compact seed header is truncated".to_string(),
        ));
    }
    if bytes[0..8] != COMPACT_SEED_MAGIC {
        return Err(JmaError::InvalidMagic);
    }
    if get_u32(bytes, 8)? != COMPACT_SEED_VERSION
        || get_u32(bytes, 12)? != u32::try_from(COMPACT_HEADER_SIZE).unwrap_or(u32::MAX)
    {
        return Err(JmaError::UnsupportedVersion(COMPACT_SEED_VERSION as u16));
    }
    let k = get_u8(bytes, 20)?;
    if !(21..=32).contains(&k) {
        return Err(JmaError::CorruptSection(format!(
            "compact seed k must be between 21 and 32, got {k}"
        )));
    }
    let key_encoding = KeyEncoding::from_byte(get_u8(bytes, 21)?).ok_or_else(|| {
        JmaError::CorruptSection("unsupported compact seed key encoding".to_string())
    })?;
    let occurrence_encoding =
        OccurrenceEncoding::from_byte(get_u8(bytes, 22)?).ok_or_else(|| {
            JmaError::CorruptSection("unsupported compact seed occurrence encoding".to_string())
        })?;
    let bucket_bits = get_u8(bytes, 23)?;
    if !(1..=MAX_BUCKET_BITS).contains(&bucket_bits) {
        return Err(JmaError::CorruptSection(
            "invalid compact seed bucket width".to_string(),
        ));
    }
    let header = SeedHeader {
        scheme_id: get_u32(bytes, 16)?,
        k,
        key_encoding,
        occurrence_encoding,
        bucket_bits,
        scale: get_u32(bytes, 24)?,
        hash_algorithm_id: get_u32(bytes, 28)?,
        key_count: get_u64(bytes, 32)?,
        posting_count: get_u64(bytes, 40)?,
        filter_offset: get_u64(bytes, 48)?,
        filter_length: get_u64(bytes, 56)?,
        prefix_offset: get_u64(bytes, 64)?,
        prefix_count: get_u32(bytes, 72)?,
        key_page_dir_offset: get_u64(bytes, 80)?,
        key_page_count: get_u32(bytes, 88)?,
        posting_header_offset: get_u64(bytes, 96)?,
        posting_header_count: get_u32(bytes, 104)?,
        position_page_dir_offset: get_u64(bytes, 112)?,
        position_page_count: get_u32(bytes, 120)?,
        key_data_offset: get_u64(bytes, 128)?,
        key_data_length: get_u64(bytes, 136)?,
        position_data_offset: get_u64(bytes, 144)?,
        position_data_length: get_u64(bytes, 152)?,
        section_length: get_u64(bytes, 160)?,
        flags: get_u64(bytes, 168)?,
        catalog_checksum: get_array(bytes, 192)?,
        archive_checksum: get_array(bytes, 224)?,
        scheme_checksum: get_array(bytes, 256)?,
        section_checksum: get_array(bytes, 288)?,
    };
    let record_sizes_valid = get_u32(bytes, 76)? == PREFIX_RECORD_SIZE as u32
        && get_u32(bytes, 92)? == KEY_PAGE_RECORD_SIZE as u32
        && get_u32(bytes, 108)? == POSTING_HEADER_RECORD_SIZE as u32
        && get_u32(bytes, 124)? == POSITION_PAGE_RECORD_SIZE as u32
        && bytes[176..192].iter().all(|byte| *byte == 0);
    if header.section_length != section_length
        || header.hash_algorithm_id != JAMHASH_ALGORITHM_ID
        || header.scale == 0
        || header.flags != 0
        || header.key_count != u64::from(header.posting_header_count)
        || header.key_count != header.posting_count
        || header.key_count > u64::from(u32::MAX)
        || header.prefix_count > header.key_count as u32
        || header.key_page_count > header.key_count as u32
        || !record_sizes_valid
    {
        return Err(JmaError::CorruptSection(
            "compact seed header counts or identity are invalid".to_string(),
        ));
    }
    Ok(header)
}
