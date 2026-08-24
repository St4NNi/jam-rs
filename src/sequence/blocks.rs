//! Independent sequence blocks and their fixed-record directory.

use super::{
    EncodedContig, SequenceError, SequenceResult, complement_base, decode_ambiguity_payload,
    decode_range, encode_ambiguity_payload, encode_contig, validate_encoded,
};
use sha2::{Digest, Sha256};
use std::ops::Range;

/// A sequence payload codec.  The decoded payload is always the canonical
/// two-bit byte stream; ambiguity records are stored independently.
#[repr(u8)]
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum BlockCodec {
    Raw2Bit = 0,
    Zstd2Bit = 1,
}

impl BlockCodec {
    #[must_use]
    pub const fn id(self) -> u8 {
        match self {
            Self::Raw2Bit => 0,
            Self::Zstd2Bit => 1,
        }
    }
}

impl TryFrom<u8> for BlockCodec {
    type Error = SequenceError;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Self::Raw2Bit),
            1 => Ok(Self::Zstd2Bit),
            other => Err(SequenceError::Codec(format!(
                "unknown sequence block codec {other}"
            ))),
        }
    }
}

/// Alias used by callers that want to distinguish block errors in signatures
/// while retaining one stable error type for the shared codec.
pub type BlockDecodeError = SequenceError;

/// A base interval selected by a fixed or content-defined block policy.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SequenceSpan {
    pub start: u64,
    pub base_count: u64,
}

impl SequenceSpan {
    #[must_use]
    pub fn end(self) -> Option<u64> {
        self.start.checked_add(self.base_count)
    }
}

/// Deterministic DNA Gear table families used by content-defined blocks.
#[repr(u8)]
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum GearTable {
    #[default]
    SingleBase = 0,
    Dinucleotide = 1,
    PackedFourBase = 2,
}

impl GearTable {
    #[must_use]
    pub const fn id(self) -> u8 {
        match self {
            Self::SingleBase => 0,
            Self::Dinucleotide => 1,
            Self::PackedFourBase => 2,
        }
    }
}

/// Sequence block split policies.  Gear boundaries are independent of the
/// storage codec and therefore can be compared fairly against fixed blocks.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SequenceBlockPolicy {
    Fixed {
        block_bases: usize,
    },
    Gear {
        min_bases: usize,
        target_bases: usize,
        max_bases: usize,
        table: GearTable,
    },
}

impl Default for SequenceBlockPolicy {
    fn default() -> Self {
        Self::Fixed {
            block_bases: 1 << 20,
        }
    }
}

/// Fixed-width binary metadata for one independently readable sequence block.
///
/// All integer fields are encoded little-endian by
/// [`SequenceBlockDirectory::encode`].  Offsets are absolute archive byte
/// offsets; the sequence module intentionally does not own the surrounding
/// archive layout.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SequenceBlockRecord {
    pub contig_id: u32,
    pub base_start: u64,
    pub base_count: u64,
    pub payload_offset: u64,
    pub stored_length: u64,
    pub decoded_length: u64,
    pub ambiguity_offset: u64,
    pub ambiguity_length: u64,
    pub codec: BlockCodec,
    pub flags: u8,
    pub checksum: [u8; 32],
}

impl SequenceBlockRecord {
    pub const ENCODED_SIZE: usize = 96;

    #[must_use]
    pub fn base_end(self) -> Option<u64> {
        self.base_start.checked_add(self.base_count)
    }

    #[must_use]
    pub fn payload_end(self) -> Option<u64> {
        self.payload_offset.checked_add(self.stored_length)
    }

    #[must_use]
    pub fn ambiguity_end(self) -> Option<u64> {
        self.ambiguity_offset.checked_add(self.ambiguity_length)
    }

    /// Returns a copy with offsets filled by an archive writer after payloads
    /// have been laid out.  The checksum remains bound to sequence content,
    /// not to placement, so deterministic layout code can set both offsets.
    #[must_use]
    pub const fn with_offsets(self, payload_offset: u64, ambiguity_offset: u64) -> Self {
        Self {
            payload_offset,
            ambiguity_offset,
            ..self
        }
    }
}

/// Encoded payloads and their checksum-bound directory record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EncodedSequenceBlock {
    pub record: SequenceBlockRecord,
    pub payload: Vec<u8>,
    pub ambiguity_payload: Vec<u8>,
}

/// Fixed directory header and records.  It is safe to decode directly from a
/// mapped or range-fetched slice without allocating the payload bytes.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SequenceBlockDirectory {
    pub records: Vec<SequenceBlockRecord>,
}

pub type SequenceBlockDirectoryError = SequenceError;

const DIRECTORY_MAGIC: [u8; 4] = *b"SQD1";
const DIRECTORY_VERSION: u16 = 1;
const DIRECTORY_HEADER_SIZE: usize = 16;
/// Maximum decoded bytes accepted by the convenience block readers.  This
/// bound is checked before zstd allocates its output buffer.
pub const DEFAULT_MAX_DECODED_BLOCK_BYTES: usize = 256 * 1024 * 1024;

impl SequenceBlockDirectory {
    pub const HEADER_SIZE: usize = DIRECTORY_HEADER_SIZE;
    pub const RECORD_SIZE: usize = SequenceBlockRecord::ENCODED_SIZE;

    pub fn new(mut records: Vec<SequenceBlockRecord>) -> SequenceResult<Self> {
        records.sort_by_key(|record| {
            (
                record.contig_id,
                record.base_start,
                record.payload_offset,
                record.ambiguity_offset,
            )
        });
        let directory = Self { records };
        directory.validate()?;
        Ok(directory)
    }

    /// Encodes the deterministic fixed-record directory.
    pub fn encode(&self) -> SequenceResult<Vec<u8>> {
        self.validate()?;
        let count = u64::try_from(self.records.len()).map_err(|_| SequenceError::LengthOverflow)?;
        let body = usize::try_from(count)
            .ok()
            .and_then(|count| count.checked_mul(Self::RECORD_SIZE))
            .ok_or(SequenceError::LengthOverflow)?;
        let length = DIRECTORY_HEADER_SIZE
            .checked_add(body)
            .ok_or(SequenceError::LengthOverflow)?;
        let mut bytes = Vec::with_capacity(length);
        bytes.extend_from_slice(&DIRECTORY_MAGIC);
        bytes.extend_from_slice(&DIRECTORY_VERSION.to_le_bytes());
        let record_size =
            u16::try_from(Self::RECORD_SIZE).map_err(|_| SequenceError::LengthOverflow)?;
        bytes.extend_from_slice(&record_size.to_le_bytes());
        bytes.extend_from_slice(&count.to_le_bytes());
        for record in &self.records {
            encode_record(&mut bytes, record);
        }
        Ok(bytes)
    }

    /// Decodes all records and validates arithmetic, ordering, and per-contig
    /// non-overlap.  Object-size checks are supplied separately because a
    /// directory may be read before the full remote object size is known.
    pub fn decode(bytes: &[u8]) -> SequenceResult<Self> {
        if bytes.len() < DIRECTORY_HEADER_SIZE || bytes[..4] != DIRECTORY_MAGIC {
            return Err(SequenceError::Directory(
                "sequence block directory header is invalid".to_string(),
            ));
        }
        let version = u16::from_le_bytes(bytes[4..6].try_into().expect("directory version"));
        if version != DIRECTORY_VERSION {
            return Err(SequenceError::Directory(format!(
                "unsupported sequence directory version {version}"
            )));
        }
        let record_size = u16::from_le_bytes(bytes[6..8].try_into().expect("record size"));
        if usize::from(record_size) != Self::RECORD_SIZE {
            return Err(SequenceError::Directory(format!(
                "unsupported sequence directory record size {record_size}"
            )));
        }
        let count = u64::from_le_bytes(bytes[8..16].try_into().expect("directory count"));
        let body = usize::try_from(count)
            .ok()
            .and_then(|count| count.checked_mul(Self::RECORD_SIZE))
            .ok_or(SequenceError::LengthOverflow)?;
        let expected = DIRECTORY_HEADER_SIZE
            .checked_add(body)
            .ok_or(SequenceError::LengthOverflow)?;
        if bytes.len() != expected {
            return Err(SequenceError::Directory(format!(
                "directory length {} does not match {count} records",
                bytes.len()
            )));
        }
        let count = usize::try_from(count).map_err(|_| SequenceError::LengthOverflow)?;
        let mut records = Vec::with_capacity(count);
        let mut cursor = DIRECTORY_HEADER_SIZE;
        for _ in 0..count {
            records.push(decode_record(&bytes[cursor..cursor + Self::RECORD_SIZE])?);
            cursor += Self::RECORD_SIZE;
        }
        let directory = Self { records };
        directory.validate()?;
        Ok(directory)
    }

    /// Validates ranges and ensures records for one contig cannot overlap.
    pub fn validate(&self) -> SequenceResult<()> {
        let mut previous: Option<SequenceBlockRecord> = None;
        for record in &self.records {
            if record.base_count == 0 || record.stored_length == 0 {
                return Err(SequenceError::Directory(
                    "empty sequence block record".to_string(),
                ));
            }
            if record.flags != 0 {
                return Err(SequenceError::UnsupportedBlockFlags(record.flags));
            }
            let expected_decoded = record
                .base_count
                .checked_add(3)
                .ok_or(SequenceError::LengthOverflow)?
                / 4;
            if record.decoded_length != expected_decoded {
                return Err(SequenceError::Directory(format!(
                    "decoded length {} does not match {} bases",
                    record.decoded_length, record.base_count
                )));
            }
            record.base_end().ok_or(SequenceError::LengthOverflow)?;
            record.payload_end().ok_or(SequenceError::LengthOverflow)?;
            record
                .ambiguity_end()
                .ok_or(SequenceError::LengthOverflow)?;
            if let Some(previous_record) = previous {
                if (record.contig_id, record.base_start)
                    < (previous_record.contig_id, previous_record.base_start)
                {
                    return Err(SequenceError::Directory(
                        "sequence block records are not sorted".to_string(),
                    ));
                }
                let previous_end = previous_record
                    .base_end()
                    .ok_or(SequenceError::LengthOverflow)?;
                if previous_record.contig_id == record.contig_id && previous_end > record.base_start
                {
                    return Err(SequenceError::Directory(
                        "sequence blocks overlap within a contig".to_string(),
                    ));
                }
            }
            previous = Some(*record);
        }
        Ok(())
    }

    /// Checks payload locations against a known archive/object size.
    pub fn validate_object_size(&self, object_size: u64) -> SequenceResult<()> {
        for record in &self.records {
            if record.payload_end().ok_or(SequenceError::LengthOverflow)? > object_size
                || record
                    .ambiguity_end()
                    .ok_or(SequenceError::LengthOverflow)?
                    > object_size
            {
                return Err(SequenceError::Directory(
                    "sequence block payload exceeds object size".to_string(),
                ));
            }
        }
        Ok(())
    }

    /// Returns the record at a sorted index without decoding unrelated
    /// payload bytes.  This is useful for mmap and paged directory readers.
    pub fn record_at(bytes: &[u8], index: usize) -> SequenceResult<SequenceBlockRecord> {
        let header = bytes
            .get(..DIRECTORY_HEADER_SIZE)
            .ok_or_else(|| SequenceError::Directory("directory header is truncated".to_string()))?;
        if header[..4] != DIRECTORY_MAGIC {
            return Err(SequenceError::Directory(
                "sequence block directory magic is invalid".to_string(),
            ));
        }
        let version = u16::from_le_bytes(header[4..6].try_into().expect("directory version"));
        let record_size = u16::from_le_bytes(header[6..8].try_into().expect("record size"));
        if version != DIRECTORY_VERSION || usize::from(record_size) != Self::RECORD_SIZE {
            return Err(SequenceError::Directory(
                "sequence block directory header is unsupported".to_string(),
            ));
        }
        let count = u64::from_le_bytes(header[8..16].try_into().expect("directory count"));
        if u64::try_from(index).map_err(|_| SequenceError::LengthOverflow)? >= count {
            return Err(SequenceError::Directory(
                "sequence block directory index is out of range".to_string(),
            ));
        }
        let offset = DIRECTORY_HEADER_SIZE
            .checked_add(
                index
                    .checked_mul(Self::RECORD_SIZE)
                    .ok_or(SequenceError::LengthOverflow)?,
            )
            .ok_or(SequenceError::LengthOverflow)?;
        let end = offset
            .checked_add(Self::RECORD_SIZE)
            .ok_or(SequenceError::LengthOverflow)?;
        decode_record(bytes.get(offset..end).ok_or_else(|| {
            SequenceError::Directory("sequence block record is truncated".to_string())
        })?)
    }
}

fn encode_record(bytes: &mut Vec<u8>, record: &SequenceBlockRecord) {
    bytes.extend_from_slice(&record.contig_id.to_le_bytes());
    bytes.push(record.codec.id());
    bytes.push(record.flags);
    bytes.extend_from_slice(&0u16.to_le_bytes());
    bytes.extend_from_slice(&record.base_start.to_le_bytes());
    bytes.extend_from_slice(&record.base_count.to_le_bytes());
    bytes.extend_from_slice(&record.payload_offset.to_le_bytes());
    bytes.extend_from_slice(&record.stored_length.to_le_bytes());
    bytes.extend_from_slice(&record.decoded_length.to_le_bytes());
    bytes.extend_from_slice(&record.ambiguity_offset.to_le_bytes());
    bytes.extend_from_slice(&record.ambiguity_length.to_le_bytes());
    bytes.extend_from_slice(&record.checksum);
}

fn decode_record(bytes: &[u8]) -> SequenceResult<SequenceBlockRecord> {
    if bytes.len() != SequenceBlockRecord::ENCODED_SIZE {
        return Err(SequenceError::Directory(
            "sequence block record has an invalid length".to_string(),
        ));
    }
    let reserved = u16::from_le_bytes(bytes[6..8].try_into().expect("reserved bytes"));
    if reserved != 0 {
        return Err(SequenceError::Directory(
            "sequence block record reserved bytes are non-zero".to_string(),
        ));
    }
    let mut checksum = [0u8; 32];
    checksum.copy_from_slice(&bytes[64..96]);
    let record = SequenceBlockRecord {
        contig_id: u32::from_le_bytes(bytes[..4].try_into().expect("contig id")),
        codec: BlockCodec::try_from(bytes[4])?,
        flags: bytes[5],
        base_start: u64::from_le_bytes(bytes[8..16].try_into().expect("base start")),
        base_count: u64::from_le_bytes(bytes[16..24].try_into().expect("base count")),
        payload_offset: u64::from_le_bytes(bytes[24..32].try_into().expect("payload offset")),
        stored_length: u64::from_le_bytes(bytes[32..40].try_into().expect("stored length")),
        decoded_length: u64::from_le_bytes(bytes[40..48].try_into().expect("decoded length")),
        ambiguity_offset: u64::from_le_bytes(bytes[48..56].try_into().expect("ambiguity offset")),
        ambiguity_length: u64::from_le_bytes(bytes[56..64].try_into().expect("ambiguity length")),
        checksum,
    };
    if record.flags != 0 {
        return Err(SequenceError::UnsupportedBlockFlags(record.flags));
    }
    Ok(record)
}

/// Splits one contig into deterministic, independently addressable spans.
pub fn split_contig(
    sequence: &[u8],
    policy: SequenceBlockPolicy,
) -> SequenceResult<Vec<SequenceSpan>> {
    match policy {
        SequenceBlockPolicy::Fixed { block_bases } => split_fixed(sequence.len(), block_bases),
        SequenceBlockPolicy::Gear {
            min_bases,
            target_bases,
            max_bases,
            table,
        } => split_gear(sequence, min_bases, target_bases, max_bases, table),
    }
}

fn split_fixed(length: usize, block_bases: usize) -> SequenceResult<Vec<SequenceSpan>> {
    if block_bases == 0 {
        return Err(SequenceError::Directory(
            "fixed sequence block size must be greater than zero".to_string(),
        ));
    }
    let mut spans = Vec::new();
    for (index, chunk) in (0..length).step_by(block_bases).enumerate() {
        let count = (length - chunk).min(block_bases);
        let start = index
            .checked_mul(block_bases)
            .and_then(|value| u64::try_from(value).ok())
            .ok_or(SequenceError::LengthOverflow)?;
        spans.push(SequenceSpan {
            start,
            base_count: u64::try_from(count).map_err(|_| SequenceError::LengthOverflow)?,
        });
    }
    Ok(spans)
}

fn split_gear(
    sequence: &[u8],
    min_bases: usize,
    target_bases: usize,
    max_bases: usize,
    table: GearTable,
) -> SequenceResult<Vec<SequenceSpan>> {
    if min_bases == 0 || target_bases < min_bases || max_bases < target_bases {
        return Err(SequenceError::Directory(
            "Gear block sizes must satisfy 0 < min <= target <= max".to_string(),
        ));
    }
    if sequence.is_empty() {
        return Ok(Vec::new());
    }
    let mask = target_bases
        .checked_next_power_of_two()
        .ok_or(SequenceError::LengthOverflow)?
        .saturating_sub(1);
    let mask = u64::try_from(mask).map_err(|_| SequenceError::LengthOverflow)?;
    let mut spans = Vec::new();
    let mut start = 0usize;
    let mut rolling = 0u64;
    for position in 0..sequence.len() {
        let token = canonical_token_hash(sequence, position, table)?;
        rolling = rolling.rotate_left(7) ^ token;
        let block_len = position + 1 - start;
        let boundary = block_len >= min_bases && (rolling & mask == 0 || block_len >= max_bases);
        if boundary {
            spans.push(SequenceSpan {
                start: u64::try_from(start).map_err(|_| SequenceError::LengthOverflow)?,
                base_count: u64::try_from(block_len).map_err(|_| SequenceError::LengthOverflow)?,
            });
            start = position + 1;
            rolling = 0;
        }
    }
    if start < sequence.len() {
        spans.push(SequenceSpan {
            start: u64::try_from(start).map_err(|_| SequenceError::LengthOverflow)?,
            base_count: u64::try_from(sequence.len() - start)
                .map_err(|_| SequenceError::LengthOverflow)?,
        });
    }
    Ok(spans)
}

fn canonical_token_hash(sequence: &[u8], position: usize, table: GearTable) -> SequenceResult<u64> {
    let width = match table {
        GearTable::SingleBase => 1,
        GearTable::Dinucleotide => 2,
        GearTable::PackedFourBase => 4,
    };
    let start = position.saturating_add(1).saturating_sub(width);
    let available = position + 1 - start;
    let mut forward = 0u64;
    for base in &sequence[start..=position] {
        forward = (forward << 4) | u64::from(base_code(*base));
    }
    let mut reverse = 0u64;
    for base in sequence[start..=position].iter().rev() {
        reverse = (reverse << 4) | u64::from(complement_code(base_code(*base)));
    }
    let available = u8::try_from(available).map_err(|_| SequenceError::LengthOverflow)?;
    let token = forward.min(reverse) ^ u64::from(available);
    Ok(gear_value(table, token))
}

fn gear_value(table: GearTable, token: u64) -> u64 {
    // SplitMix constants are deterministic table-specific Gear values.  They
    // avoid a large static table while preserving stable archive bytes.
    let salt = match table {
        GearTable::SingleBase => 0x7369_6e67_6c65_6261,
        GearTable::Dinucleotide => 0x6469_6e75_636c_656f,
        GearTable::PackedFourBase => 0x7061_636b_3462_6173,
    };
    splitmix64(token ^ salt)
}

fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

fn base_code(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 4,
    }
}

fn complement_code(code: u8) -> u8 {
    match code {
        0 => 3,
        1 => 2,
        2 => 1,
        3 => 0,
        _ => 4,
    }
}

/// Encodes one independent sequence block.
pub fn encode_sequence_block(
    contig_id: u32,
    base_start: u64,
    sequence: &[u8],
    codec: BlockCodec,
) -> SequenceResult<EncodedSequenceBlock> {
    if sequence.is_empty() {
        return Err(SequenceError::Directory(
            "cannot encode an empty sequence block".to_string(),
        ));
    }
    let encoded = encode_contig(sequence)?;
    base_start
        .checked_add(encoded.base_count)
        .ok_or(SequenceError::LengthOverflow)?;
    let ambiguity_payload = encode_ambiguity_payload(&encoded.ambiguities)?;
    let payload = match codec {
        BlockCodec::Raw2Bit => encoded.two_bit.clone(),
        BlockCodec::Zstd2Bit => zstd::bulk::compress(&encoded.two_bit, 3)
            .map_err(|error| SequenceError::Codec(error.to_string()))?,
    };
    let checksum = sequence_block_checksum(
        contig_id,
        base_start,
        encoded.base_count,
        codec,
        &encoded.two_bit,
        &ambiguity_payload,
    );
    let record = SequenceBlockRecord {
        contig_id,
        base_start,
        base_count: encoded.base_count,
        payload_offset: 0,
        stored_length: u64::try_from(payload.len()).map_err(|_| SequenceError::LengthOverflow)?,
        decoded_length: u64::try_from(encoded.two_bit.len())
            .map_err(|_| SequenceError::LengthOverflow)?,
        ambiguity_offset: 0,
        ambiguity_length: u64::try_from(ambiguity_payload.len())
            .map_err(|_| SequenceError::LengthOverflow)?,
        codec,
        flags: 0,
        checksum,
    };
    Ok(EncodedSequenceBlock {
        record,
        payload,
        ambiguity_payload,
    })
}

/// Decodes a complete independent block after checking its checksum.
pub fn decode_stored_block(
    record: &SequenceBlockRecord,
    payload: &[u8],
    ambiguity_payload: &[u8],
) -> SequenceResult<Vec<u8>> {
    decode_stored_block_bounded(
        record,
        payload,
        ambiguity_payload,
        DEFAULT_MAX_DECODED_BLOCK_BYTES,
    )
}

/// Decodes a complete block with an explicit hard-bounded output limit.  The
/// configured limit may be lowered for a query, but cannot exceed
/// [`DEFAULT_MAX_DECODED_BLOCK_BYTES`].  The limit is checked before zstd
/// allocates its output buffer.
pub fn decode_stored_block_bounded(
    record: &SequenceBlockRecord,
    payload: &[u8],
    ambiguity_payload: &[u8],
    max_decoded_bytes: usize,
) -> SequenceResult<Vec<u8>> {
    let encoded = decode_block_parts(record, payload, ambiguity_payload, max_decoded_bytes)?;
    decode_range(&encoded, 0..encoded.base_count)
}

/// Decodes only the requested global range from an independent block.
pub fn decode_block_range(
    record: &SequenceBlockRecord,
    payload: &[u8],
    ambiguity_payload: &[u8],
    range: Range<u64>,
) -> SequenceResult<Vec<u8>> {
    decode_block_range_bounded(
        record,
        payload,
        ambiguity_payload,
        range,
        DEFAULT_MAX_DECODED_BLOCK_BYTES,
    )
}

/// Decodes a block range with an explicit hard-bounded output limit.  This is
/// the remote-reader entry point: a corrupt directory cannot turn a declared
/// zstd decoded length into an unbounded allocation.
pub fn decode_block_range_bounded(
    record: &SequenceBlockRecord,
    payload: &[u8],
    ambiguity_payload: &[u8],
    range: Range<u64>,
    max_decoded_bytes: usize,
) -> SequenceResult<Vec<u8>> {
    let block_end = record.base_end().ok_or(SequenceError::LengthOverflow)?;
    if range.start > range.end || range.start < record.base_start || range.end > block_end {
        return Err(SequenceError::InvalidRange {
            start: range.start,
            end: range.end,
            length: block_end,
        });
    }
    if record.codec == BlockCodec::Raw2Bit {
        return decode_raw_range(record, payload, ambiguity_payload, range, max_decoded_bytes);
    }
    let encoded = decode_block_parts(record, payload, ambiguity_payload, max_decoded_bytes)?;
    decode_range(
        &encoded,
        (range.start - record.base_start)..(range.end - record.base_start),
    )
}

fn decode_raw_range(
    record: &SequenceBlockRecord,
    payload: &[u8],
    ambiguity_payload: &[u8],
    range: Range<u64>,
    max_decoded_bytes: usize,
) -> SequenceResult<Vec<u8>> {
    let ambiguities = validate_raw_payload(record, payload, ambiguity_payload, max_decoded_bytes)?;
    let local_start = usize::try_from(range.start - record.base_start)
        .map_err(|_| SequenceError::LengthOverflow)?;
    let local_end = usize::try_from(range.end - record.base_start)
        .map_err(|_| SequenceError::LengthOverflow)?;
    let local_start_u64 = u64::try_from(local_start).map_err(|_| SequenceError::LengthOverflow)?;
    let mut ambiguity_index =
        ambiguities.partition_point(|ambiguity| ambiguity.position < local_start_u64);
    let mut output = Vec::with_capacity(local_end - local_start);
    for position in local_start..local_end {
        let position_u64 = u64::try_from(position).map_err(|_| SequenceError::LengthOverflow)?;
        if ambiguities
            .get(ambiguity_index)
            .is_some_and(|ambiguity| ambiguity.position == position_u64)
        {
            output.push(ambiguities[ambiguity_index].base);
            ambiguity_index += 1;
        } else {
            let shift = (position % 4) * 2;
            let value = (payload[position / 4] >> shift) & 0b11;
            output.push(b"ACGT"[usize::from(value)]);
        }
    }
    Ok(output)
}

fn validate_raw_payload(
    record: &SequenceBlockRecord,
    payload: &[u8],
    ambiguity_payload: &[u8],
    max_decoded_bytes: usize,
) -> SequenceResult<Vec<super::AmbiguityBase>> {
    validate_decode_limit(record.decoded_length, max_decoded_bytes)?;
    if record.flags != 0 {
        return Err(SequenceError::UnsupportedBlockFlags(record.flags));
    }
    let stored_length =
        usize::try_from(record.stored_length).map_err(|_| SequenceError::LengthOverflow)?;
    let decoded_length =
        usize::try_from(record.decoded_length).map_err(|_| SequenceError::LengthOverflow)?;
    let ambiguity_length =
        usize::try_from(record.ambiguity_length).map_err(|_| SequenceError::LengthOverflow)?;
    if payload.len() != stored_length
        || payload.len() != decoded_length
        || ambiguity_payload.len() != ambiguity_length
    {
        return Err(SequenceError::Codec(
            "raw sequence block payload length does not match directory record".to_string(),
        ));
    }
    let expected_decoded = record
        .base_count
        .checked_add(3)
        .ok_or(SequenceError::LengthOverflow)?
        / 4;
    if record.base_count == 0 || record.decoded_length != expected_decoded {
        return Err(SequenceError::Codec(
            "raw sequence block base length is invalid".to_string(),
        ));
    }
    let ambiguities = decode_ambiguity_payload(ambiguity_payload)?;
    for ambiguity in &ambiguities {
        if ambiguity.position >= record.base_count {
            return Err(SequenceError::InvalidAmbiguityPosition(ambiguity.position));
        }
    }
    if !record.base_count.is_multiple_of(4) && !payload.is_empty() {
        let used = u32::try_from((record.base_count % 4) * 2)
            .map_err(|_| SequenceError::LengthOverflow)?;
        let mask = !((1u8 << used) - 1);
        if payload[payload.len() - 1] & mask != 0 {
            return Err(SequenceError::NonZeroPaddingBits);
        }
    }
    let expected = sequence_block_checksum(
        record.contig_id,
        record.base_start,
        record.base_count,
        record.codec,
        payload,
        ambiguity_payload,
    );
    if expected != record.checksum {
        return Err(SequenceError::ChecksumMismatch);
    }
    Ok(ambiguities)
}

/// Decodes only a requested block range and returns it in reverse-complement
/// orientation.  Coordinates remain in forward contig space.
pub fn decode_block_reverse_complement_range(
    record: &SequenceBlockRecord,
    payload: &[u8],
    ambiguity_payload: &[u8],
    range: Range<u64>,
) -> SequenceResult<Vec<u8>> {
    let mut decoded = decode_block_range(record, payload, ambiguity_payload, range)?;
    decoded.reverse();
    for base in &mut decoded {
        *base = complement_base(*base)?;
    }
    Ok(decoded)
}

fn decode_block_parts(
    record: &SequenceBlockRecord,
    payload: &[u8],
    ambiguity_payload: &[u8],
    max_decoded_bytes: usize,
) -> SequenceResult<EncodedContig> {
    if record.flags != 0 {
        return Err(SequenceError::UnsupportedBlockFlags(record.flags));
    }
    let expected_decoded = record
        .base_count
        .checked_add(3)
        .ok_or(SequenceError::LengthOverflow)?
        / 4;
    if record.base_count == 0 || record.decoded_length != expected_decoded {
        return Err(SequenceError::Codec(
            "sequence block base length is invalid".to_string(),
        ));
    }
    validate_decode_limit(record.decoded_length, max_decoded_bytes)?;
    let stored_length =
        usize::try_from(record.stored_length).map_err(|_| SequenceError::LengthOverflow)?;
    let decoded_length =
        usize::try_from(record.decoded_length).map_err(|_| SequenceError::LengthOverflow)?;
    let ambiguity_length =
        usize::try_from(record.ambiguity_length).map_err(|_| SequenceError::LengthOverflow)?;
    if payload.len() != stored_length || ambiguity_payload.len() != ambiguity_length {
        return Err(SequenceError::Codec(
            "sequence block payload length does not match directory record".to_string(),
        ));
    }
    let raw = match record.codec {
        BlockCodec::Raw2Bit => {
            if payload.len() != decoded_length {
                return Err(SequenceError::Codec(
                    "raw sequence block decoded length mismatch".to_string(),
                ));
            }
            payload.to_vec()
        }
        BlockCodec::Zstd2Bit => zstd::bulk::decompress(payload, decoded_length)
            .map_err(|error| SequenceError::Codec(error.to_string()))?,
    };
    let ambiguities = decode_ambiguity_payload(ambiguity_payload)?;
    let encoded = EncodedContig {
        base_count: record.base_count,
        two_bit: raw,
        ambiguities,
    };
    validate_encoded(&encoded)?;
    let expected = sequence_block_checksum(
        record.contig_id,
        record.base_start,
        record.base_count,
        record.codec,
        &encoded.two_bit,
        ambiguity_payload,
    );
    if expected != record.checksum {
        return Err(SequenceError::ChecksumMismatch);
    }
    Ok(encoded)
}

fn validate_decode_limit(decoded_length: u64, max_decoded_bytes: usize) -> SequenceResult<()> {
    let maximum = u64::try_from(DEFAULT_MAX_DECODED_BLOCK_BYTES)
        .map_err(|_| SequenceError::LengthOverflow)?;
    let configured = u64::try_from(max_decoded_bytes).map_err(|_| SequenceError::LengthOverflow)?;
    if configured > maximum {
        return Err(SequenceError::DecodedLengthLimit {
            requested: configured,
            maximum,
        });
    }
    if decoded_length > configured {
        return Err(SequenceError::DecodedLengthLimit {
            requested: decoded_length,
            maximum: configured,
        });
    }
    Ok(())
}

/// Hashes canonical uncompressed block content and coordinates.  Placement
/// offsets are deliberately excluded so moving a block within a deterministic
/// archive does not alter its content identity.
pub fn sequence_block_checksum(
    contig_id: u32,
    base_start: u64,
    base_count: u64,
    codec: BlockCodec,
    raw_two_bit: &[u8],
    ambiguity_payload: &[u8],
) -> [u8; 32] {
    let mut hasher = Sha256::new();
    hasher.update(b"jam.sequence.block.1\0");
    hasher.update(contig_id.to_le_bytes());
    hasher.update(base_start.to_le_bytes());
    hasher.update(base_count.to_le_bytes());
    hasher.update([codec.id()]);
    hasher.update(raw_two_bit);
    hasher.update([0]);
    hasher.update(ambiguity_payload);
    hasher.finalize().into()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::{decode_reverse_complement_range, encode_contig};

    #[test]
    fn fixed_blocks_cover_without_crossing() {
        let spans =
            split_contig(b"ACGTNACGT", SequenceBlockPolicy::Fixed { block_bases: 4 }).unwrap();
        assert_eq!(
            spans,
            vec![
                SequenceSpan {
                    start: 0,
                    base_count: 4
                },
                SequenceSpan {
                    start: 4,
                    base_count: 4
                },
                SequenceSpan {
                    start: 8,
                    base_count: 1
                },
            ]
        );
    }

    #[test]
    fn gear_policies_are_deterministic_and_bounded() {
        let sequence = b"ACGT".repeat(1000);
        for table in [
            GearTable::SingleBase,
            GearTable::Dinucleotide,
            GearTable::PackedFourBase,
        ] {
            let policy = SequenceBlockPolicy::Gear {
                min_bases: 64,
                target_bases: 128,
                max_bases: 256,
                table,
            };
            let first = split_contig(&sequence, policy).unwrap();
            assert_eq!(first, split_contig(&sequence, policy).unwrap());
            assert!(first.iter().all(|span| span.base_count <= 256));
            assert_eq!(
                first.iter().map(|span| span.base_count).sum::<u64>(),
                u64::try_from(sequence.len()).unwrap()
            );
        }
    }

    #[test]
    fn raw_and_zstd_blocks_roundtrip_ranges_and_reverse_complements() {
        let sequence = b"ACGTRYKMBVDHN-ACGT";
        for codec in [BlockCodec::Raw2Bit, BlockCodec::Zstd2Bit] {
            let block = encode_sequence_block(3, 17, sequence, codec).unwrap();
            let record = block.record.with_offsets(100, 200);
            assert_eq!(
                decode_block_range(&record, &block.payload, &block.ambiguity_payload, 20..29)
                    .unwrap(),
                b"TRYKMBVDH"
            );
            let encoded = encode_contig(sequence).unwrap();
            let expected = decode_reverse_complement_range(&encoded, 2..11).unwrap();
            let mut actual =
                decode_block_range(&record, &block.payload, &block.ambiguity_payload, 19..28)
                    .unwrap();
            actual.reverse();
            for base in &mut actual {
                *base = crate::sequence::complement_base(*base).unwrap();
            }
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn directory_is_fixed_width_deterministic_and_corruption_checked() {
        let block = encode_sequence_block(1, 0, b"ACGTN", BlockCodec::Raw2Bit).unwrap();
        let record = block.record.with_offsets(256, 512);
        let directory = SequenceBlockDirectory::new(vec![record]).unwrap();
        let bytes = directory.encode().unwrap();
        assert_eq!(bytes, directory.encode().unwrap());
        assert_eq!(SequenceBlockDirectory::decode(&bytes).unwrap(), directory);
        assert_eq!(
            SequenceBlockDirectory::record_at(&bytes, 0).unwrap(),
            record
        );
        let mut corrupt = bytes.clone();
        corrupt[22] = 0xff;
        assert!(SequenceBlockDirectory::decode(&corrupt).is_err());
    }

    #[test]
    fn checksum_collision_or_payload_corruption_never_decodes() {
        let block = encode_sequence_block(2, 0, b"ACGTACGT", BlockCodec::Raw2Bit).unwrap();
        let mut payload = block.payload.clone();
        payload[0] ^= 1;
        assert_eq!(
            decode_stored_block(&block.record, &payload, &block.ambiguity_payload),
            Err(SequenceError::ChecksumMismatch)
        );
    }
}
