//! Shared two-bit sequence codec used by archive backends.

use std::ops::Range;
use thiserror::Error;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AmbiguityBase {
    pub position: u64,
    pub base: u8,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EncodedContig {
    pub base_count: u64,
    pub two_bit: Vec<u8>,
    pub ambiguities: Vec<AmbiguityBase>,
}

pub fn encode_contig(sequence: &[u8]) -> SequenceResult<EncodedContig> {
    let base_count = u64::try_from(sequence.len()).map_err(|_| SequenceError::LengthOverflow)?;
    let packed_len = sequence
        .len()
        .checked_add(3)
        .ok_or(SequenceError::LengthOverflow)?
        / 4;
    let mut two_bit = vec![0u8; packed_len];
    for (index, base) in sequence.iter().copied().enumerate() {
        let value = encode_base(base).unwrap_or(0);
        let shift = (index % 4) * 2;
        two_bit[index / 4] |= value << shift;
    }
    Ok(EncodedContig {
        base_count,
        two_bit,
        ambiguities: encode_ambiguity(sequence)?,
    })
}

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

pub fn validate_encoded_range(
    encoded: &EncodedContig,
    range: Range<u64>,
) -> SequenceResult<Range<usize>> {
    if range.start > range.end || range.end > encoded.base_count {
        return Err(SequenceError::InvalidRange {
            start: range.start,
            end: range.end,
            length: encoded.base_count,
        });
    }
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
    Ok(
        usize::try_from(range.start).map_err(|_| SequenceError::LengthOverflow)?
            ..usize::try_from(range.end).map_err(|_| SequenceError::LengthOverflow)?,
    )
}

pub fn decode_range(encoded: &EncodedContig, range: Range<u64>) -> SequenceResult<Vec<u8>> {
    let range = validate_encoded_range(encoded, range)?;
    let mut decoded = Vec::with_capacity(range.len());
    for position in range {
        let shift = (position % 4) * 2;
        let value = (encoded.two_bit[position / 4] >> shift) & 0b11;
        let position_u64 = u64::try_from(position).map_err(|_| SequenceError::LengthOverflow)?;
        let base = encoded
            .ambiguities
            .binary_search_by_key(&position_u64, |ambiguity| ambiguity.position)
            .map_or_else(
                |_| decode_base(value),
                |index| encoded.ambiguities[index].base,
            );
        decoded.push(base);
    }
    Ok(decoded)
}

pub fn decode_reverse_complement_range(
    encoded: &EncodedContig,
    range: Range<u64>,
) -> SequenceResult<Vec<u8>> {
    let mut decoded = decode_range(encoded, range)?;
    decoded.reverse();
    for base in &mut decoded {
        *base = complement(*base)?;
    }
    Ok(decoded)
}

fn encode_base(base: u8) -> Option<u8> {
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

fn normalize_ambiguity(base: u8) -> SequenceResult<u8> {
    let base = base.to_ascii_uppercase();
    match base {
        b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' | b'N' | b'-' => {
            Ok(base)
        }
        _ => Err(SequenceError::UnsupportedBase(base)),
    }
}

fn complement(base: u8) -> SequenceResult<u8> {
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

#[derive(Debug, Error)]
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
}
