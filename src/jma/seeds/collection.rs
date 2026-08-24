//! Deterministic container for the multiple seed schemes in one JMA object.
//!
//! The compact scheme blobs remain independently checksummed and addressable;
//! this directory only binds each scheme identity to its relative range.

use super::format::{get_array, get_u8, get_u32, get_u64, put_u32, put_u64};
use super::{SeedBinding, SeedIndex};
use crate::jma::{JmaError, JmaResult};
use std::collections::BTreeSet;

pub const COLLECTION_MAGIC: [u8; 8] = *b"JMSC1SE1";
pub const COLLECTION_VERSION: u32 = 1;
pub const COLLECTION_HEADER_SIZE: usize = 128;
pub const COLLECTION_ENTRY_SIZE: usize = 64;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedCollectionHeader {
    pub scheme_count: u32,
    pub directory_offset: u64,
    pub directory_length: u64,
    pub collection_length: u64,
    pub flags: u64,
    pub checksum: [u8; 32],
}

/// One compact seed blob range within the complete `SECTION_SEEDS` payload.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedCollectionEntry {
    pub scheme_id: u32,
    pub k: u8,
    pub scale: u32,
    pub offset: u64,
    pub length: u64,
    pub checksum: [u8; 32],
}

#[derive(Clone, Debug)]
pub struct SeedCollection {
    bytes: Vec<u8>,
    header: SeedCollectionHeader,
    entries: Vec<SeedCollectionEntry>,
}

/// Encode all compact scheme blobs into one deterministic seed section.
pub fn encode_seed_collection(schemes: &[(SeedBinding, Vec<u8>)]) -> JmaResult<Vec<u8>> {
    let mut ordered = schemes.to_vec();
    ordered.sort_by_key(|(binding, _)| (binding.scheme_id, binding.k, binding.scale));
    let mut identities = BTreeSet::new();
    let mut blobs = Vec::with_capacity(ordered.len());
    for (binding, bytes) in ordered {
        binding.validate()?;
        let index = SeedIndex::open(bytes.clone(), Some(binding))?;
        let header = index.header();
        if !identities.insert((binding.scheme_id, binding.k, binding.scale)) {
            return Err(JmaError::CorruptSection(
                "duplicate compact seed scheme identity".to_string(),
            ));
        }
        if header.scheme_id != binding.scheme_id
            || header.k != binding.k
            || header.scale != binding.scale
        {
            return Err(JmaError::CorruptSection(
                "compact seed blob identity does not match its directory entry".to_string(),
            ));
        }
        blobs.push((binding, bytes, header.section_checksum));
    }
    let scheme_count = u32::try_from(blobs.len()).map_err(|_| JmaError::OffsetOverflow)?;
    let directory_offset =
        u64::try_from(COLLECTION_HEADER_SIZE).map_err(|_| JmaError::OffsetOverflow)?;
    let directory_length = u64::try_from(
        blobs
            .len()
            .checked_mul(COLLECTION_ENTRY_SIZE)
            .ok_or(JmaError::OffsetOverflow)?,
    )
    .map_err(|_| JmaError::OffsetOverflow)?;
    let mut next_offset = directory_offset
        .checked_add(directory_length)
        .ok_or(JmaError::OffsetOverflow)?;
    let mut entries = Vec::with_capacity(blobs.len());
    for (binding, bytes, checksum) in &blobs {
        let offset = next_offset;
        let length = u64::try_from(bytes.len()).map_err(|_| JmaError::OffsetOverflow)?;
        next_offset = next_offset
            .checked_add(length)
            .ok_or(JmaError::OffsetOverflow)?;
        entries.push(SeedCollectionEntry {
            scheme_id: binding.scheme_id,
            k: binding.k,
            scale: binding.scale,
            offset,
            length,
            checksum: *checksum,
        });
    }
    let mut header = SeedCollectionHeader {
        scheme_count,
        directory_offset,
        directory_length,
        collection_length: next_offset,
        flags: 0,
        checksum: [0; 32],
    };
    let mut output = encode_collection_header(header);
    for entry in &entries {
        encode_collection_entry(&mut output, *entry);
    }
    for (_, bytes, _) in &blobs {
        output.extend_from_slice(bytes);
    }
    if u64::try_from(output.len()).map_err(|_| JmaError::OffsetOverflow)? != next_offset {
        return Err(JmaError::CorruptSection(
            "compact seed collection length is inconsistent".to_string(),
        ));
    }
    header.checksum = collection_checksum(&output);
    output[..COLLECTION_HEADER_SIZE].copy_from_slice(&encode_collection_header(header));
    Ok(output)
}

impl SeedCollection {
    pub fn open(bytes: Vec<u8>) -> JmaResult<Self> {
        let length = u64::try_from(bytes.len()).map_err(|_| JmaError::OffsetOverflow)?;
        let header = decode_collection_header(&bytes, length)?;
        let mut checksum_input = bytes.clone();
        checksum_input[96..128].fill(0);
        if crate::jma::format::checksum(&checksum_input) != header.checksum {
            return Err(JmaError::ChecksumMismatch(
                "compact seed collection".to_string(),
            ));
        }
        let entries = decode_collection_directory(&bytes, &header)?;
        validate_collection_entries(&header, &entries)?;
        for entry in &entries {
            let range = byte_range(entry.offset, entry.length)?;
            let blob = bytes.get(range).ok_or_else(|| {
                JmaError::CorruptSection("compact seed scheme blob is truncated".to_string())
            })?;
            let index = SeedIndex::open(blob.to_vec(), None)?;
            if index.header().scheme_id != entry.scheme_id
                || index.header().k != entry.k
                || index.header().scale != entry.scale
                || index.header().section_checksum != entry.checksum
            {
                return Err(JmaError::CorruptSection(
                    "compact seed collection entry identity mismatch".to_string(),
                ));
            }
        }
        Ok(Self {
            bytes,
            header,
            entries,
        })
    }

    #[must_use]
    pub const fn header(&self) -> &SeedCollectionHeader {
        &self.header
    }

    #[must_use]
    pub fn entries(&self) -> &[SeedCollectionEntry] {
        &self.entries
    }

    pub fn entry(&self, scheme_id: u32) -> JmaResult<SeedCollectionEntry> {
        self.entries
            .binary_search_by_key(&scheme_id, |entry| entry.scheme_id)
            .map(|index| self.entries[index])
            .map_err(|_| {
                JmaError::CorruptSection(format!("compact seed scheme {scheme_id} is unavailable"))
            })
    }

    pub fn open_scheme(
        &self,
        scheme_id: u32,
        expected: Option<SeedBinding>,
    ) -> JmaResult<SeedIndex> {
        let entry = self.entry(scheme_id)?;
        if let Some(expected) = expected
            && (expected.scheme_id != entry.scheme_id
                || expected.k != entry.k
                || expected.scale != entry.scale)
        {
            return Err(JmaError::CorruptSection(
                "compact seed scheme binding does not match collection entry".to_string(),
            ));
        }
        let range = byte_range(entry.offset, entry.length)?;
        let blob = self.bytes.get(range).ok_or_else(|| {
            JmaError::CorruptSection("compact seed scheme blob is unavailable".to_string())
        })?;
        SeedIndex::open(blob.to_vec(), expected)
    }
}

pub(crate) fn encode_collection_header(header: SeedCollectionHeader) -> Vec<u8> {
    let mut bytes = vec![0u8; COLLECTION_HEADER_SIZE];
    bytes[0..8].copy_from_slice(&COLLECTION_MAGIC);
    put_u32(&mut bytes, 8, COLLECTION_VERSION);
    put_u32(
        &mut bytes,
        12,
        u32::try_from(COLLECTION_HEADER_SIZE).expect("collection header fits u32"),
    );
    put_u32(&mut bytes, 16, header.scheme_count);
    put_u32(&mut bytes, 20, COLLECTION_ENTRY_SIZE as u32);
    put_u64(&mut bytes, 24, header.directory_offset);
    put_u64(&mut bytes, 32, header.directory_length);
    put_u64(&mut bytes, 40, header.collection_length);
    put_u64(&mut bytes, 48, header.flags);
    bytes[96..128].copy_from_slice(&header.checksum);
    bytes
}

pub fn decode_collection_header(
    bytes: &[u8],
    section_length: u64,
) -> JmaResult<SeedCollectionHeader> {
    if bytes.len() < COLLECTION_HEADER_SIZE {
        return Err(JmaError::CorruptSection(
            "compact seed collection header is truncated".to_string(),
        ));
    }
    if bytes[0..8] != COLLECTION_MAGIC {
        return Err(JmaError::InvalidMagic);
    }
    if get_u32(bytes, 8)? != COLLECTION_VERSION
        || get_u32(bytes, 12)? != u32::try_from(COLLECTION_HEADER_SIZE).unwrap_or(u32::MAX)
        || get_u32(bytes, 20)? != COLLECTION_ENTRY_SIZE as u32
    {
        return Err(JmaError::UnsupportedVersion(COLLECTION_VERSION as u16));
    }
    let header = SeedCollectionHeader {
        scheme_count: get_u32(bytes, 16)?,
        directory_offset: get_u64(bytes, 24)?,
        directory_length: get_u64(bytes, 32)?,
        collection_length: get_u64(bytes, 40)?,
        flags: get_u64(bytes, 48)?,
        checksum: get_array(bytes, 96)?,
    };
    let expected_directory_length = u64::from(header.scheme_count)
        .checked_mul(COLLECTION_ENTRY_SIZE as u64)
        .ok_or(JmaError::OffsetOverflow)?;
    let directory_end = header
        .directory_offset
        .checked_add(header.directory_length)
        .ok_or(JmaError::OffsetOverflow)?;
    if header.collection_length != section_length
        || header.directory_offset != u64::try_from(COLLECTION_HEADER_SIZE).unwrap_or(u64::MAX)
        || header.directory_length != expected_directory_length
        || directory_end > section_length
        || header.flags != 0
        || bytes[56..96].iter().any(|byte| *byte != 0)
    {
        return Err(JmaError::CorruptSection(
            "compact seed collection header is invalid".to_string(),
        ));
    }
    Ok(header)
}

pub fn decode_collection_directory(
    bytes: &[u8],
    header: &SeedCollectionHeader,
) -> JmaResult<Vec<SeedCollectionEntry>> {
    let range = byte_range(header.directory_offset, header.directory_length)?;
    if range.end > bytes.len() {
        return Err(JmaError::CorruptSection(
            "compact seed collection directory is truncated".to_string(),
        ));
    }
    let mut entries = Vec::with_capacity(header.scheme_count as usize);
    for offset in (range.start..range.end).step_by(COLLECTION_ENTRY_SIZE) {
        if bytes[offset + 5..offset + 8].iter().any(|byte| *byte != 0)
            || get_u32(bytes, offset + 12)? != 0
        {
            return Err(JmaError::CorruptSection(
                "compact seed collection entry reserved bytes are non-zero".to_string(),
            ));
        }
        entries.push(SeedCollectionEntry {
            scheme_id: get_u32(bytes, offset)?,
            k: get_u8(bytes, offset + 4)?,
            scale: get_u32(bytes, offset + 8)?,
            offset: get_u64(bytes, offset + 16)?,
            length: get_u64(bytes, offset + 24)?,
            checksum: get_array(bytes, offset + 32)?,
        });
    }
    Ok(entries)
}

pub fn validate_collection_entries(
    header: &SeedCollectionHeader,
    entries: &[SeedCollectionEntry],
) -> JmaResult<()> {
    let mut scheme_ids = BTreeSet::new();
    if entries.len()
        != usize::try_from(header.scheme_count).map_err(|_| JmaError::OffsetOverflow)?
        || entries.windows(2).any(|pair| {
            (pair[0].scheme_id, pair[0].k, pair[0].scale)
                >= (pair[1].scheme_id, pair[1].k, pair[1].scale)
        })
        || entries
            .iter()
            .any(|entry| !scheme_ids.insert(entry.scheme_id))
    {
        return Err(JmaError::CorruptSection(
            "compact seed collection entries are not deterministic".to_string(),
        ));
    }
    let payload_start = header
        .directory_offset
        .checked_add(header.directory_length)
        .ok_or(JmaError::OffsetOverflow)?;
    let mut previous_end = payload_start;
    for entry in entries {
        if !(21..=32).contains(&entry.k)
            || entry.scale == 0
            || entry.length == 0
            || entry.offset < payload_start
            || entry.offset < previous_end
        {
            return Err(JmaError::CorruptSection(
                "compact seed collection entry range is invalid".to_string(),
            ));
        }
        let end = entry
            .offset
            .checked_add(entry.length)
            .ok_or(JmaError::OffsetOverflow)?;
        if end > header.collection_length {
            return Err(JmaError::CorruptSection(
                "compact seed collection entry exceeds section".to_string(),
            ));
        }
        previous_end = end;
    }
    Ok(())
}

pub(crate) fn collection_checksum(bytes: &[u8]) -> [u8; 32] {
    let mut input = bytes.to_vec();
    if input.len() >= COLLECTION_HEADER_SIZE {
        input[96..128].fill(0);
    }
    crate::jma::format::checksum(&input)
}

fn encode_collection_entry(bytes: &mut Vec<u8>, entry: SeedCollectionEntry) {
    bytes.extend_from_slice(&entry.scheme_id.to_le_bytes());
    bytes.push(entry.k);
    bytes.extend_from_slice(&[0, 0, 0]);
    bytes.extend_from_slice(&entry.scale.to_le_bytes());
    bytes.extend_from_slice(&0u32.to_le_bytes());
    bytes.extend_from_slice(&entry.offset.to_le_bytes());
    bytes.extend_from_slice(&entry.length.to_le_bytes());
    bytes.extend_from_slice(&entry.checksum);
}

fn byte_range(offset: u64, length: u64) -> JmaResult<std::ops::Range<usize>> {
    let end = offset.checked_add(length).ok_or(JmaError::OffsetOverflow)?;
    Ok(
        usize::try_from(offset).map_err(|_| JmaError::OffsetOverflow)?
            ..usize::try_from(end).map_err(|_| JmaError::OffsetOverflow)?,
    )
}
