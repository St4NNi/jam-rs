//! Backend-neutral sequence codecs used by native JMA and alternate archives.
//!
//! Coordinates are zero-based half-open base ranges.  The two-bit payload is
//! little-endian within each byte (the first base occupies the low two bits),
//! while ambiguity records carry their original IUPAC symbol separately.  A
//! block never owns bases from more than one contig.

mod blocks;

pub use blocks::{
    BlockCodec, BlockDecodeError, DEFAULT_MAX_DECODED_BLOCK_BYTES, EncodedSequenceBlock, GearTable,
    SequenceBlockDirectory, SequenceBlockDirectoryError, SequenceBlockPolicy, SequenceBlockRecord,
    SequenceSpan, decode_block_range, decode_block_range_bounded,
    decode_block_reverse_complement_range, decode_stored_block, decode_stored_block_bounded,
    encode_sequence_block, sequence_block_checksum, split_contig,
};

use std::ops::Range;
use thiserror::Error;

/// A non-ACGT symbol stored at a contig-relative base position.
///
/// The value is normalized to uppercase and includes every IUPAC ambiguity
/// code (`R Y S W K M B D H V N`) plus `-`, which is accepted by FASTA
/// readers as a gap/unknown symbol. A record stores the symbol rather than
/// reducing all ambiguities to `N`, so a decoded query can preserve its input.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct AmbiguityBase {
    pub position: u64,
    pub base: u8,
}

/// Canonically packed representation of one complete contig.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EncodedContig {
    pub base_count: u64,
    pub two_bit: Vec<u8>,
    pub ambiguities: Vec<AmbiguityBase>,
}

/// Fixed-width ambiguity record size in the binary payload.
pub const AMBIGUITY_RECORD_SIZE: usize = 9;

/// Encodes A/C/G/T (and U as T) into two bits and keeps all other supported
/// IUPAC symbols in a sorted side payload.
pub fn encode_contig(sequence: &[u8]) -> SequenceResult<EncodedContig> {
    let base_count = u64::try_from(sequence.len()).map_err(|_| SequenceError::LengthOverflow)?;
    let packed_len = sequence
        .len()
        .checked_add(3)
        .ok_or(SequenceError::LengthOverflow)?
        / 4;
    let mut two_bit = vec![0u8; packed_len];
    for (index, base) in sequence.iter().copied().enumerate() {
        let value = encode_base(base).unwrap_or_default();
        // Low-order bits are used deliberately: this permits direct range
        // decoding without shifting an entire block into a temporary buffer.
        let shift = (index % 4) * 2;
        two_bit[index / 4] |= value << shift;
    }
    let ambiguities = encode_ambiguity(sequence)?;
    let encoded = EncodedContig {
        base_count,
        two_bit,
        ambiguities,
    };
    validate_encoded(&encoded)?;
    Ok(encoded)
}

/// Encodes all supported non-ACGT symbols, retaining their exact IUPAC code.
pub fn encode_ambiguity(sequence: &[u8]) -> SequenceResult<Vec<AmbiguityBase>> {
    sequence
        .iter()
        .copied()
        .enumerate()
        .filter(|(_, base)| encode_base(*base).is_none())
        .map(|(position, base)| {
            Ok(AmbiguityBase {
                position: u64::try_from(position).map_err(|_| SequenceError::LengthOverflow)?,
                base: normalize_ambiguity(base)?,
            })
        })
        .collect()
}

/// Validates the complete packed representation, including sorted ambiguity
/// records and unused high bits in the final packed byte.
pub fn validate_encoded(encoded: &EncodedContig) -> SequenceResult<()> {
    let expected_packed = encoded
        .base_count
        .checked_add(3)
        .ok_or(SequenceError::LengthOverflow)?
        / 4;
    let expected_packed =
        usize::try_from(expected_packed).map_err(|_| SequenceError::LengthOverflow)?;
    if encoded.two_bit.len() != expected_packed {
        return Err(SequenceError::InvalidPackedLength {
            actual: encoded.two_bit.len(),
            expected: expected_packed,
        });
    }
    if !encoded.base_count.is_multiple_of(4) && !encoded.two_bit.is_empty() {
        let used = u32::try_from((encoded.base_count % 4) * 2)
            .map_err(|_| SequenceError::LengthOverflow)?;
        let mask = !((1u8 << used) - 1);
        if encoded.two_bit[encoded.two_bit.len() - 1] & mask != 0 {
            return Err(SequenceError::NonZeroPaddingBits);
        }
    }
    let mut previous = None;
    for ambiguity in &encoded.ambiguities {
        if ambiguity.position >= encoded.base_count
            || previous.is_some_and(|position| ambiguity.position <= position)
        {
            return Err(SequenceError::InvalidAmbiguityPosition(ambiguity.position));
        }
        normalize_ambiguity(ambiguity.base)?;
        previous = Some(ambiguity.position);
    }
    Ok(())
}

/// Validates a range and returns the representable `usize` coordinates used
/// by an in-memory decoder.
pub fn validate_encoded_range(
    encoded: &EncodedContig,
    range: Range<u64>,
) -> SequenceResult<Range<usize>> {
    validate_encoded(encoded)?;
    if range.start > range.end || range.end > encoded.base_count {
        return Err(SequenceError::InvalidRange {
            start: range.start,
            end: range.end,
            length: encoded.base_count,
        });
    }
    Ok(
        usize::try_from(range.start).map_err(|_| SequenceError::LengthOverflow)?
            ..usize::try_from(range.end).map_err(|_| SequenceError::LengthOverflow)?,
    )
}

/// Decodes only the requested forward range.
pub fn decode_range(encoded: &EncodedContig, range: Range<u64>) -> SequenceResult<Vec<u8>> {
    let range = validate_encoded_range(encoded, range)?;
    let range_start = u64::try_from(range.start).map_err(|_| SequenceError::LengthOverflow)?;
    let mut decoded = Vec::with_capacity(range.len());
    let mut ambiguity_index = encoded
        .ambiguities
        .partition_point(|ambiguity| ambiguity.position < range_start);
    for position in range {
        let position_u64 = u64::try_from(position).map_err(|_| SequenceError::LengthOverflow)?;
        let base = if encoded
            .ambiguities
            .get(ambiguity_index)
            .is_some_and(|ambiguity| ambiguity.position == position_u64)
        {
            let base = encoded.ambiguities[ambiguity_index].base;
            ambiguity_index += 1;
            base
        } else {
            let shift = (position % 4) * 2;
            decode_base((encoded.two_bit[position / 4] >> shift) & 0b11)
        };
        decoded.push(base);
    }
    Ok(decoded)
}

/// Decodes a requested range in reverse-complement orientation. The range is
/// interpreted in forward contig coordinates before reverse-complementing.
pub fn decode_reverse_complement_range(
    encoded: &EncodedContig,
    range: Range<u64>,
) -> SequenceResult<Vec<u8>> {
    let mut decoded = decode_range(encoded, range)?;
    decoded.reverse();
    for base in &mut decoded {
        *base = complement_base(*base)?;
    }
    Ok(decoded)
}

/// Encodes ambiguity records as a deterministic fixed-width little-endian
/// payload. The payload starts with a u32 record count followed by
/// `(position:u64, symbol:u8)` records.
pub fn encode_ambiguity_payload(ambiguities: &[AmbiguityBase]) -> SequenceResult<Vec<u8>> {
    let count = u32::try_from(ambiguities.len()).map_err(|_| SequenceError::LengthOverflow)?;
    let body_len = usize::try_from(count)
        .ok()
        .and_then(|count| count.checked_mul(AMBIGUITY_RECORD_SIZE))
        .ok_or(SequenceError::LengthOverflow)?;
    let mut payload = Vec::with_capacity(4 + body_len);
    payload.extend_from_slice(&count.to_le_bytes());
    let mut previous = None;
    for ambiguity in ambiguities {
        if previous.is_some_and(|position| ambiguity.position <= position) {
            return Err(SequenceError::InvalidAmbiguityPosition(ambiguity.position));
        }
        normalize_ambiguity(ambiguity.base)?;
        payload.extend_from_slice(&ambiguity.position.to_le_bytes());
        payload.push(ambiguity.base.to_ascii_uppercase());
        previous = Some(ambiguity.position);
    }
    Ok(payload)
}

/// Decodes and validates a fixed-width ambiguity payload.
pub fn decode_ambiguity_payload(payload: &[u8]) -> SequenceResult<Vec<AmbiguityBase>> {
    if payload.len() < 4 {
        return Err(SequenceError::InvalidAmbiguityPayload);
    }
    let count = u32::from_le_bytes(payload[..4].try_into().expect("four-byte count"));
    let body_len = usize::try_from(count)
        .ok()
        .and_then(|count| count.checked_mul(AMBIGUITY_RECORD_SIZE))
        .ok_or(SequenceError::LengthOverflow)?;
    let expected = 4usize
        .checked_add(body_len)
        .ok_or(SequenceError::LengthOverflow)?;
    if payload.len() != expected {
        return Err(SequenceError::InvalidAmbiguityPayload);
    }
    let count = usize::try_from(count).map_err(|_| SequenceError::LengthOverflow)?;
    let mut result = Vec::with_capacity(count);
    let mut cursor = 4usize;
    let mut previous = None;
    for _ in 0..count {
        let position = u64::from_le_bytes(
            payload[cursor..cursor + 8]
                .try_into()
                .expect("eight-byte ambiguity position"),
        );
        let base = normalize_ambiguity(payload[cursor + 8])?;
        if previous.is_some_and(|old| position <= old) {
            return Err(SequenceError::InvalidAmbiguityPosition(position));
        }
        result.push(AmbiguityBase { position, base });
        previous = Some(position);
        cursor += AMBIGUITY_RECORD_SIZE;
    }
    Ok(result)
}

/// Returns the canonical two-bit value for an unambiguous DNA base.
pub fn encode_base(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

fn decode_base(value: u8) -> u8 {
    b"ACGT"[usize::from(value)]
}

/// Normalizes one ambiguity symbol and rejects non-IUPAC input.
pub fn normalize_ambiguity(base: u8) -> SequenceResult<u8> {
    let base = base.to_ascii_uppercase();
    match base {
        b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' | b'N' | b'-' => {
            Ok(base)
        }
        _ => Err(SequenceError::UnsupportedBase(base)),
    }
}

/// Complements an IUPAC base, preserving ambiguity semantics.
pub fn complement_base(base: u8) -> SequenceResult<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Ok(b'T'),
        b'C' => Ok(b'G'),
        b'G' => Ok(b'C'),
        b'T' | b'U' => Ok(b'A'),
        b'R' => Ok(b'Y'),
        b'Y' => Ok(b'R'),
        b'S' => Ok(b'S'),
        b'W' => Ok(b'W'),
        b'K' => Ok(b'M'),
        b'M' => Ok(b'K'),
        b'B' => Ok(b'V'),
        b'D' => Ok(b'H'),
        b'H' => Ok(b'D'),
        b'V' => Ok(b'B'),
        b'N' => Ok(b'N'),
        b'-' => Ok(b'-'),
        other => Err(SequenceError::UnsupportedBase(other)),
    }
}

pub type SequenceResult<T> = Result<T, SequenceError>;

#[derive(Debug, Error, Eq, PartialEq)]
pub enum SequenceError {
    #[error("sequence length does not fit the codec coordinate space")]
    LengthOverflow,
    #[error("invalid sequence range [{start}, {end}) for contig length {length}")]
    InvalidRange { start: u64, end: u64, length: u64 },
    #[error("invalid two-bit payload length {actual}; expected {expected}")]
    InvalidPackedLength { actual: usize, expected: usize },
    #[error("invalid ambiguity position {0}")]
    InvalidAmbiguityPosition(u64),
    #[error("unsupported sequence base 0x{0:02x}")]
    UnsupportedBase(u8),
    #[error("non-zero padding bits in the final two-bit byte")]
    NonZeroPaddingBits,
    #[error("invalid ambiguity payload")]
    InvalidAmbiguityPayload,
    #[error("unsupported sequence block flags 0x{0:02x}")]
    UnsupportedBlockFlags(u8),
    #[error("decoded sequence block length {requested} exceeds limit {maximum}")]
    DecodedLengthLimit { requested: u64, maximum: u64 },
    #[error("sequence block codec failure: {0}")]
    Codec(String),
    #[error("sequence block checksum mismatch")]
    ChecksumMismatch,
    #[error("sequence block directory is corrupt: {0}")]
    Directory(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn range_and_reverse_complement_roundtrip_ambiguities() {
        let encoded = encode_contig(b"acgtRYKMBVDHN-").unwrap();
        assert_eq!(
            decode_range(&encoded, 0..encoded.base_count).unwrap(),
            b"ACGTRYKMBVDHN-"
        );
        assert_eq!(
            decode_reverse_complement_range(&encoded, 0..encoded.base_count).unwrap(),
            b"-NDHBVKMRYACGT"
        );
    }

    #[test]
    fn invalid_ranges_and_payloads_fail_closed() {
        let mut encoded = encode_contig(b"ACGTN").unwrap();
        assert!(decode_range(&encoded, 4..6).is_err());
        encoded.two_bit.clear();
        assert!(decode_range(&encoded, 0..1).is_err());
    }

    #[test]
    fn every_iupac_symbol_roundtrips_in_binary_payload() {
        let input = b"RYSWKMBDHVN-";
        let encoded = encode_contig(input).unwrap();
        let payload = encode_ambiguity_payload(&encoded.ambiguities).unwrap();
        assert_eq!(
            decode_ambiguity_payload(&payload).unwrap(),
            encoded.ambiguities
        );
        assert_eq!(
            decode_range(&encoded, 0..u64::try_from(input.len()).unwrap()).unwrap(),
            input
        );
    }
}
