//! Deterministic writer for JAM Witness Router format 1.

use crate::router::format::{
    KeyBlockConfig, KeyBlockDescriptor, MembershipFilter, MetagenomeEntry, PackedKeyWidth,
    PostingHeader, PrefixEntry, RouterFormatError, SECTION_FILTER, SECTION_FLAG_REQUIRED,
    SECTION_KEYS, SECTION_METAGENOMES, SECTION_POSITION_PAYLOAD, SECTION_POSTING_HEADERS,
    SECTION_POSTING_PAYLOAD, SECTION_PREFIX, SECTION_SCHEMES, SUPERBLOCK_SIZE, checksum,
    encode_filter, encode_header, encode_metagenomes, encode_schemes, encode_section_directory,
    finalize_object, hash_prefix, put_u32, put_u64,
};
use crate::router::{WitnessKey, WitnessScheme};
use std::collections::BTreeSet;
use std::fs::File;
use std::io::Write;
use std::path::Path;

const KEY_HEADER_SIZE: usize = 40;
const KEY_BLOCK_RECORD_SIZE: usize = 88;
const PREFIX_HEADER_SIZE: usize = 24;
const PREFIX_RECORD_SIZE: usize = 56;
const POSTING_HEADER_SIZE: usize = 24;
const POSTING_RECORD_SIZE: usize = 80;

/// Raw data associated with one exact packed witness key.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RouterKeyInput {
    pub key: WitnessKey,
    pub document_frequency: u32,
    pub flags: u32,
    pub posting: Vec<u8>,
    pub positions: Option<Vec<u8>>,
}

impl RouterKeyInput {
    #[must_use]
    pub fn new(key: WitnessKey, posting: Vec<u8>) -> Self {
        Self {
            key,
            document_frequency: 1,
            flags: 0,
            posting,
            positions: None,
        }
    }
}

/// All immutable source data needed to construct one router object.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RouterBuildInput {
    pub metagenomes: Vec<MetagenomeEntry>,
    pub schemes: Vec<WitnessScheme>,
    pub keys: Vec<RouterKeyInput>,
}

/// Physical options are intentionally small and measurable. They do not
/// alter the logical witness set.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct RouterWriterOptions {
    pub key_blocks: KeyBlockConfig,
    pub prefix_bits: u8,
    pub filter_bits_per_key: u32,
    pub flags: u32,
}

impl Default for RouterWriterOptions {
    fn default() -> Self {
        Self {
            key_blocks: KeyBlockConfig::default(),
            prefix_bits: 16,
            filter_bits_per_key: 10,
            flags: 0,
        }
    }
}

impl RouterWriterOptions {
    pub fn validate(self) -> Result<(), RouterFormatError> {
        self.key_blocks.validate()?;
        if self.prefix_bits > 32 || self.prefix_bits == 0 {
            return Err(RouterFormatError::InvalidPrefixDirectory);
        }
        if self.filter_bits_per_key < 4 {
            return Err(RouterFormatError::InvalidFilter);
        }
        Ok(())
    }
}

/// Build an entire deterministic JWR object in memory.
pub fn build_router(
    input: &RouterBuildInput,
    options: RouterWriterOptions,
) -> Result<Vec<u8>, RouterFormatError> {
    options.validate()?;
    validate_input(input)?;

    let metagenome_bytes = encode_metagenomes(&input.metagenomes)?;
    let scheme_bytes = encode_schemes(&input.schemes)?;

    let mut ordered = input.keys.clone();
    ordered.sort_by_key(|record| (record.key.jamhash, record.key.packed));
    let key_bytes = encode_keys(&ordered, options.key_blocks, options.prefix_bits)?;
    let prefix_bytes = encode_prefix(&key_bytes.prefix_entries, options.prefix_bits)?;
    let filter = MembershipFilter::build(
        ordered.iter().map(|record| record.key.jamhash),
        options.filter_bits_per_key,
    )?;
    let filter_bytes = encode_filter(&filter)?;
    let (header_bytes, posting_bytes, position_bytes) = encode_postings(&ordered)?;

    let mut sections = vec![
        (SECTION_METAGENOMES, SECTION_FLAG_REQUIRED, metagenome_bytes),
        (SECTION_SCHEMES, SECTION_FLAG_REQUIRED, scheme_bytes),
        (SECTION_FILTER, SECTION_FLAG_REQUIRED, filter_bytes),
        (SECTION_PREFIX, SECTION_FLAG_REQUIRED, prefix_bytes),
        (SECTION_KEYS, SECTION_FLAG_REQUIRED, key_bytes.bytes),
        (SECTION_POSTING_HEADERS, SECTION_FLAG_REQUIRED, header_bytes),
        (
            SECTION_POSTING_PAYLOAD,
            SECTION_FLAG_REQUIRED,
            posting_bytes,
        ),
    ];
    if !position_bytes.is_empty() {
        sections.push((
            SECTION_POSITION_PAYLOAD,
            SECTION_FLAG_REQUIRED,
            position_bytes,
        ));
    }

    let mut object = vec![0u8; SUPERBLOCK_SIZE];
    let mut descriptors = Vec::with_capacity(sections.len());
    for (kind, flags, payload) in sections {
        align_object(&mut object, 8);
        let offset = u64::try_from(object.len()).map_err(|_| RouterFormatError::OffsetOverflow)?;
        let length = u64::try_from(payload.len()).map_err(|_| RouterFormatError::OffsetOverflow)?;
        let section_checksum = checksum(&payload);
        object.extend_from_slice(&payload);
        descriptors.push(crate::router::format::SectionDescriptor {
            kind,
            flags,
            offset,
            length,
            checksum: section_checksum,
        });
    }
    align_object(&mut object, 8);
    let directory_offset =
        u64::try_from(object.len()).map_err(|_| RouterFormatError::OffsetOverflow)?;
    let directory = encode_section_directory(&descriptors)?;
    let directory_length =
        u64::try_from(directory.len()).map_err(|_| RouterFormatError::OffsetOverflow)?;
    object.extend_from_slice(&directory);
    let object_size = u64::try_from(object.len()).map_err(|_| RouterFormatError::OffsetOverflow)?;
    let header = encode_header(
        options.flags,
        u32::try_from(descriptors.len())
            .map_err(|_| RouterFormatError::CountTooLarge("sections"))?,
        directory_offset,
        directory_length,
        object_size,
        u32::try_from(input.metagenomes.len())
            .map_err(|_| RouterFormatError::CountTooLarge("metagenomes"))?,
        u32::try_from(input.schemes.len())
            .map_err(|_| RouterFormatError::CountTooLarge("schemes"))?,
        checksum(&encode_metagenomes(&input.metagenomes)?),
        checksum(&encode_schemes(&input.schemes)?),
        checksum(&directory),
    )?;
    object[..SUPERBLOCK_SIZE].copy_from_slice(&header);
    let (object, _) = finalize_object(object)?;
    Ok(object)
}

/// Build and write one router without creating a temporary, mutable format.
pub fn write_router(
    path: impl AsRef<Path>,
    input: &RouterBuildInput,
    options: RouterWriterOptions,
) -> Result<(), RouterFormatError> {
    let bytes = build_router(input, options)?;
    let mut file = File::create(path)
        .map_err(|error| RouterFormatError::InvalidDirectory(error.to_string()))?;
    file.write_all(&bytes)
        .map_err(|error| RouterFormatError::InvalidDirectory(error.to_string()))?;
    file.sync_all()
        .map_err(|error| RouterFormatError::InvalidDirectory(error.to_string()))
}

fn validate_input(input: &RouterBuildInput) -> Result<(), RouterFormatError> {
    if input.schemes.is_empty() {
        return Err(RouterFormatError::InvalidSchemeTable);
    }
    let mut packed = BTreeSet::new();
    for record in &input.keys {
        let checked = WitnessKey::checked(record.key.packed, record.key.jamhash)
            .map_err(|_| RouterFormatError::InvalidKeyBlock)?;
        if checked != record.key || !packed.insert(record.key.packed) {
            return Err(RouterFormatError::InvalidKeyBlock);
        }
    }
    Ok(())
}

struct EncodedKeys {
    bytes: Vec<u8>,
    prefix_entries: Vec<PrefixEntry>,
}

fn encode_keys(
    records: &[RouterKeyInput],
    config: KeyBlockConfig,
    prefix_bits: u8,
) -> Result<EncodedKeys, RouterFormatError> {
    let record_bytes = config.encoded_record_bytes();
    let per_block = (usize::try_from(config.target_block_bytes)
        .map_err(|_| RouterFormatError::OffsetOverflow)?
        .saturating_sub(KEY_HEADER_SIZE + KEY_BLOCK_RECORD_SIZE)
        / record_bytes)
        .max(1);
    let block_count = records.len().div_ceil(per_block);
    let directory_offset = KEY_HEADER_SIZE;
    let data_offset = directory_offset
        .checked_add(
            block_count
                .checked_mul(KEY_BLOCK_RECORD_SIZE)
                .ok_or(RouterFormatError::OffsetOverflow)?,
        )
        .ok_or(RouterFormatError::OffsetOverflow)?;
    let mut descriptors = Vec::with_capacity(block_count);
    let mut prefix_entries = Vec::new();
    let mut payload = Vec::new();
    for (block_id, chunk) in records.chunks(per_block).enumerate() {
        let offset = u64::try_from(data_offset + payload.len())
            .map_err(|_| RouterFormatError::OffsetOverflow)?;
        let mut block = Vec::with_capacity(chunk.len() * record_bytes);
        let mut prefix = None;
        for (local_index, record) in chunk.iter().enumerate() {
            let current_prefix = hash_prefix(record.key.jamhash, prefix_bits)?;
            if prefix != Some(current_prefix) {
                if let Some(previous_prefix) = prefix.take() {
                    let block_local_start = chunk[..local_index]
                        .iter()
                        .position(|candidate| {
                            hash_prefix(candidate.key.jamhash, prefix_bits).ok()
                                == Some(previous_prefix)
                        })
                        .unwrap_or(0);
                    let entry_offset = offset
                        .checked_add(
                            u64::try_from(block_local_start * record_bytes)
                                .map_err(|_| RouterFormatError::OffsetOverflow)?,
                        )
                        .ok_or(RouterFormatError::OffsetOverflow)?;
                    let entry_length =
                        u64::try_from((local_index - block_local_start) * record_bytes)
                            .map_err(|_| RouterFormatError::OffsetOverflow)?;
                    prefix_entries.push(PrefixEntry {
                        prefix: previous_prefix,
                        block_id: u32::try_from(block_id)
                            .map_err(|_| RouterFormatError::CountTooLarge("key blocks"))?,
                        first_key_index: u64::try_from((block_id * per_block) + block_local_start)
                            .map_err(|_| RouterFormatError::OffsetOverflow)?,
                        key_count: u32::try_from(local_index - block_local_start)
                            .map_err(|_| RouterFormatError::CountTooLarge("key block keys"))?,
                        offset: entry_offset,
                        length: entry_length,
                        first_hash: chunk[block_local_start].key.jamhash,
                        last_hash: chunk[local_index - 1].key.jamhash,
                    });
                }
                prefix = Some(current_prefix);
            }
            append_packed(&mut block, record.key.packed, config.packed_width);
            if config.store_jamhash {
                block.extend_from_slice(&record.key.jamhash.to_le_bytes());
            }
        }
        if let Some(last_prefix) = prefix {
            let start = chunk
                .iter()
                .position(|candidate| {
                    hash_prefix(candidate.key.jamhash, prefix_bits).ok() == Some(last_prefix)
                })
                .unwrap_or(0);
            prefix_entries.push(PrefixEntry {
                prefix: last_prefix,
                block_id: u32::try_from(block_id)
                    .map_err(|_| RouterFormatError::CountTooLarge("key blocks"))?,
                first_key_index: u64::try_from((block_id * per_block) + start)
                    .map_err(|_| RouterFormatError::OffsetOverflow)?,
                key_count: u32::try_from(chunk.len() - start)
                    .map_err(|_| RouterFormatError::CountTooLarge("key block keys"))?,
                offset: offset
                    + u64::try_from(start * record_bytes)
                        .map_err(|_| RouterFormatError::OffsetOverflow)?,
                length: u64::try_from((chunk.len() - start) * record_bytes)
                    .map_err(|_| RouterFormatError::OffsetOverflow)?,
                first_hash: chunk[start].key.jamhash,
                last_hash: chunk[chunk.len() - 1].key.jamhash,
            });
        }
        descriptors.push(KeyBlockDescriptor {
            block_id: u32::try_from(block_id)
                .map_err(|_| RouterFormatError::CountTooLarge("key blocks"))?,
            first_key_index: u64::try_from(block_id * per_block)
                .map_err(|_| RouterFormatError::OffsetOverflow)?,
            key_count: u32::try_from(chunk.len())
                .map_err(|_| RouterFormatError::CountTooLarge("key block keys"))?,
            offset,
            length: u64::try_from(block.len()).map_err(|_| RouterFormatError::OffsetOverflow)?,
            first_hash: chunk.first().map_or(0, |record| record.key.jamhash),
            last_hash: chunk.last().map_or(0, |record| record.key.jamhash),
            checksum: checksum(&block),
        });
        payload.extend_from_slice(&block);
    }
    let mut bytes = vec![0u8; data_offset];
    put_u32(&mut bytes, 0, 1);
    bytes[4] = config.packed_width.bytes() as u8;
    bytes[5] = u8::from(config.store_jamhash);
    bytes[6] = prefix_bits;
    put_u32(
        &mut bytes,
        8,
        u32::try_from(block_count).map_err(|_| RouterFormatError::CountTooLarge("key blocks"))?,
    );
    put_u64(
        &mut bytes,
        12,
        u64::try_from(records.len()).map_err(|_| RouterFormatError::OffsetOverflow)?,
    );
    put_u64(
        &mut bytes,
        20,
        u64::try_from(directory_offset).map_err(|_| RouterFormatError::OffsetOverflow)?,
    );
    put_u64(
        &mut bytes,
        28,
        u64::try_from(data_offset).map_err(|_| RouterFormatError::OffsetOverflow)?,
    );
    for (index, descriptor) in descriptors.iter().enumerate() {
        encode_key_block_descriptor(
            &mut bytes,
            directory_offset + index * KEY_BLOCK_RECORD_SIZE,
            *descriptor,
        )?;
    }
    bytes.extend_from_slice(&payload);
    Ok(EncodedKeys {
        bytes,
        prefix_entries,
    })
}

fn encode_key_block_descriptor(
    bytes: &mut [u8],
    offset: usize,
    descriptor: KeyBlockDescriptor,
) -> Result<(), RouterFormatError> {
    put_u32(bytes, offset, descriptor.block_id);
    put_u64(bytes, offset + 8, descriptor.first_key_index);
    put_u32(bytes, offset + 16, descriptor.key_count);
    put_u64(bytes, offset + 24, descriptor.offset);
    put_u64(bytes, offset + 32, descriptor.length);
    put_u64(bytes, offset + 40, descriptor.first_hash);
    put_u64(bytes, offset + 48, descriptor.last_hash);
    bytes[offset + 56..offset + 88].copy_from_slice(&descriptor.checksum);
    Ok(())
}

fn encode_prefix(entries: &[PrefixEntry], prefix_bits: u8) -> Result<Vec<u8>, RouterFormatError> {
    let mut entries = entries.to_vec();
    entries.sort_by_key(|entry| (entry.prefix, entry.block_id, entry.first_key_index));
    let mut bytes = vec![0u8; PREFIX_HEADER_SIZE + entries.len() * PREFIX_RECORD_SIZE];
    put_u32(&mut bytes, 0, 1);
    bytes[4] = prefix_bits;
    put_u32(
        &mut bytes,
        8,
        u32::try_from(entries.len())
            .map_err(|_| RouterFormatError::CountTooLarge("prefix entries"))?,
    );
    put_u64(
        &mut bytes,
        12,
        u64::try_from(PREFIX_HEADER_SIZE).map_err(|_| RouterFormatError::OffsetOverflow)?,
    );
    for (index, entry) in entries.iter().enumerate() {
        let offset = PREFIX_HEADER_SIZE + index * PREFIX_RECORD_SIZE;
        put_u32(&mut bytes, offset, entry.prefix);
        put_u32(&mut bytes, offset + 4, entry.block_id);
        put_u64(&mut bytes, offset + 8, entry.first_key_index);
        put_u32(&mut bytes, offset + 16, entry.key_count);
        put_u64(&mut bytes, offset + 24, entry.offset);
        put_u64(&mut bytes, offset + 32, entry.length);
        put_u64(&mut bytes, offset + 40, entry.first_hash);
        put_u64(&mut bytes, offset + 48, entry.last_hash);
    }
    Ok(bytes)
}

#[allow(clippy::type_complexity)]
fn encode_postings(
    records: &[RouterKeyInput],
) -> Result<(Vec<u8>, Vec<u8>, Vec<u8>), RouterFormatError> {
    let headers_offset = POSTING_HEADER_SIZE;
    let records_end = headers_offset
        .checked_add(
            records
                .len()
                .checked_mul(POSTING_RECORD_SIZE)
                .ok_or(RouterFormatError::OffsetOverflow)?,
        )
        .ok_or(RouterFormatError::OffsetOverflow)?;
    let mut headers = vec![0u8; records_end];
    put_u32(&mut headers, 0, 1);
    put_u32(&mut headers, 4, POSTING_RECORD_SIZE as u32);
    put_u64(
        &mut headers,
        8,
        u64::try_from(records.len()).map_err(|_| RouterFormatError::OffsetOverflow)?,
    );
    put_u64(
        &mut headers,
        16,
        u64::try_from(headers_offset).map_err(|_| RouterFormatError::OffsetOverflow)?,
    );
    let mut posting_payload = Vec::new();
    let mut position_payload = Vec::new();
    for (index, record) in records.iter().enumerate() {
        let posting_offset =
            u64::try_from(posting_payload.len()).map_err(|_| RouterFormatError::OffsetOverflow)?;
        posting_payload.extend_from_slice(&record.posting);
        let position_offset =
            u64::try_from(position_payload.len()).map_err(|_| RouterFormatError::OffsetOverflow)?;
        if let Some(positions) = &record.positions {
            position_payload.extend_from_slice(positions);
        }
        let position_length = u64::try_from(position_payload.len())
            .map_err(|_| RouterFormatError::OffsetOverflow)?
            .checked_sub(position_offset)
            .ok_or(RouterFormatError::OffsetOverflow)?;
        let mut digest_input = Vec::with_capacity(
            record.posting.len() + usize::try_from(position_length).unwrap_or(0),
        );
        digest_input.extend_from_slice(&record.posting);
        if let Some(positions) = &record.positions {
            digest_input.extend_from_slice(positions);
        }
        let mut flags = record.flags;
        if position_length != 0 {
            flags |= PostingHeader::FLAG_POSITION_BEARING;
        }
        let descriptor = PostingHeader {
            document_frequency: record.document_frequency,
            posting_count: record.document_frequency,
            flags,
            posting_offset,
            posting_length: u64::try_from(record.posting.len())
                .map_err(|_| RouterFormatError::OffsetOverflow)?,
            position_offset,
            position_length,
            checksum: checksum(&digest_input),
        };
        encode_posting_header(
            &mut headers,
            headers_offset + index * POSTING_RECORD_SIZE,
            descriptor,
        );
    }
    Ok((headers, posting_payload, position_payload))
}

fn encode_posting_header(bytes: &mut [u8], offset: usize, header: PostingHeader) {
    put_u32(bytes, offset, header.document_frequency);
    put_u32(bytes, offset + 4, header.posting_count);
    put_u32(bytes, offset + 8, header.flags);
    put_u64(bytes, offset + 16, header.posting_offset);
    put_u64(bytes, offset + 24, header.posting_length);
    put_u64(bytes, offset + 32, header.position_offset);
    put_u64(bytes, offset + 40, header.position_length);
    bytes[offset + 48..offset + 80].copy_from_slice(&header.checksum);
}

fn append_packed(bytes: &mut Vec<u8>, packed: u64, width: PackedKeyWidth) {
    let raw = packed.to_le_bytes();
    bytes.extend_from_slice(&raw[..width.bytes()]);
}

fn align_object(bytes: &mut Vec<u8>, alignment: usize) {
    let remainder = bytes.len() % alignment;
    if remainder != 0 {
        bytes.resize(bytes.len() + (alignment - remainder), 0);
    }
}
