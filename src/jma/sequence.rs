//! Packed sequence block encoding and checked range decoding.

use crate::jma::contigs::find_contig;
use crate::jma::format::{SECTION_VERSION, checked_slice, get_u32, get_u64};
use crate::jma::{ContigId, ContigMetadata, JmaError, JmaResult, SequenceRange};

const BLOCK_HEADER_SIZE: usize = 36;

/// A packed sequence block. `unknown_mask` has one bit per base; set bits
/// decode as `N`, while unset bits decode from the two-bit `packed` payload.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SequenceBlock {
    pub contig_id: ContigId,
    pub start: u64,
    pub base_length: u64,
    pub packed: Vec<u8>,
    pub unknown_mask: Vec<u8>,
}

impl SequenceBlock {
    #[must_use]
    pub fn end(&self) -> Option<u64> {
        self.start.checked_add(self.base_length)
    }
}

/// Encodes all blocks in a sequence section.
pub fn encode_sequence_blocks(blocks: &[SequenceBlock]) -> JmaResult<Vec<u8>> {
    let count = u32::try_from(blocks.len()).map_err(|_| JmaError::OffsetOverflow)?;
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&SECTION_VERSION.to_le_bytes());
    bytes.extend_from_slice(&count.to_le_bytes());
    let mut previous = None;
    for block in blocks {
        validate_block_shape(block)?;
        let key = (block.contig_id, block.start);
        if previous.is_some_and(|previous_key| key < previous_key) {
            return Err(JmaError::CorruptSection(
                "sequence blocks must be sorted by contig and start".to_string(),
            ));
        }
        previous = Some(key);
        bytes.extend_from_slice(&block.contig_id.to_le_bytes());
        bytes.extend_from_slice(&block.start.to_le_bytes());
        bytes.extend_from_slice(&block.base_length.to_le_bytes());
        bytes.extend_from_slice(
            &u64::try_from(block.packed.len())
                .map_err(|_| JmaError::OffsetOverflow)?
                .to_le_bytes(),
        );
        bytes.extend_from_slice(
            &u64::try_from(block.unknown_mask.len())
                .map_err(|_| JmaError::OffsetOverflow)?
                .to_le_bytes(),
        );
        bytes.extend_from_slice(&block.packed);
        bytes.extend_from_slice(&block.unknown_mask);
    }
    Ok(bytes)
}

/// Decodes every block in a sequence section and checks its byte boundaries.
pub fn decode_sequence_blocks(bytes: &[u8]) -> JmaResult<Vec<SequenceBlock>> {
    if bytes.len() < 8 {
        return Err(JmaError::CorruptSection(
            "sequence section header is truncated".to_string(),
        ));
    }
    if get_u32(bytes, 0)? != SECTION_VERSION {
        return Err(JmaError::CorruptSection(
            "unsupported sequence section version".to_string(),
        ));
    }
    let count = get_u32(bytes, 4)?;
    let minimum_block_bytes = (count as usize)
        .checked_mul(BLOCK_HEADER_SIZE)
        .ok_or(JmaError::OffsetOverflow)?;
    let minimum_length = 8usize
        .checked_add(minimum_block_bytes)
        .ok_or(JmaError::OffsetOverflow)?;
    if bytes.len() < minimum_length {
        return Err(JmaError::CorruptSection(format!(
            "sequence section is too short for {count} blocks"
        )));
    }
    let mut cursor = 8usize;
    let mut previous = None;
    let mut blocks = Vec::with_capacity(count as usize);
    for _ in 0..count {
        let end = cursor
            .checked_add(BLOCK_HEADER_SIZE)
            .ok_or(JmaError::OffsetOverflow)?;
        if end > bytes.len() {
            return Err(JmaError::CorruptSection(
                "sequence block header is truncated".to_string(),
            ));
        }
        let contig_id = get_u32(bytes, cursor)?;
        let start = get_u64(bytes, cursor + 4)?;
        let base_length = get_u64(bytes, cursor + 12)?;
        let packed_length = get_u64(bytes, cursor + 20)?;
        let mask_length = get_u64(bytes, cursor + 28)?;
        cursor = end;
        let packed = checked_slice(bytes, cursor, packed_length, "packed sequence")?.to_vec();
        cursor = cursor
            .checked_add(usize::try_from(packed_length).map_err(|_| JmaError::OffsetOverflow)?)
            .ok_or(JmaError::OffsetOverflow)?;
        let unknown_mask = checked_slice(bytes, cursor, mask_length, "unknown mask")?.to_vec();
        cursor = cursor
            .checked_add(usize::try_from(mask_length).map_err(|_| JmaError::OffsetOverflow)?)
            .ok_or(JmaError::OffsetOverflow)?;
        let block = SequenceBlock {
            contig_id,
            start,
            base_length,
            packed,
            unknown_mask,
        };
        validate_block_shape(&block)?;
        let key = (contig_id, start);
        if previous.is_some_and(|previous_key| key < previous_key) {
            return Err(JmaError::CorruptSection(
                "sequence blocks are not sorted by contig and start".to_string(),
            ));
        }
        previous = Some(key);
        blocks.push(block);
    }
    if cursor != bytes.len() {
        return Err(JmaError::CorruptSection(format!(
            "sequence section has {} trailing bytes",
            bytes.len() - cursor
        )));
    }
    Ok(blocks)
}

/// Verifies block ordering, packed lengths, and optional contig bounds.
pub fn validate_blocks(blocks: &[SequenceBlock], contigs: &[ContigMetadata]) -> JmaResult<()> {
    let mut previous = None;
    for block in blocks {
        validate_block_shape(block)?;
        let metadata = find_contig(contigs, block.contig_id)?;
        let end = block.end().ok_or(JmaError::OffsetOverflow)?;
        if end > metadata.length {
            return Err(JmaError::CorruptSection(format!(
                "sequence block {}:{} exceeds contig length {}",
                block.contig_id, block.start, metadata.length
            )));
        }
        let key = (block.contig_id, block.start);
        if previous.is_some_and(|previous_key| key < previous_key) {
            return Err(JmaError::CorruptSection(
                "sequence blocks are not sorted".to_string(),
            ));
        }
        previous = Some(key);
    }
    Ok(())
}

/// Encodes bases using two bits per base and an ambiguity bitmap. U is
/// accepted as T; all other non-ACGT symbols are represented as N.
pub fn pack_bases(sequence: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let packed_len = sequence.len().div_ceil(4);
    let mask_len = sequence.len().div_ceil(8);
    let mut packed = vec![0u8; packed_len];
    let mut unknown_mask = vec![0u8; mask_len];
    for (index, base) in sequence.iter().copied().enumerate() {
        let (code, unknown) = match base.to_ascii_uppercase() {
            b'A' => (0, false),
            b'C' => (1, false),
            b'G' => (2, false),
            b'T' | b'U' => (3, false),
            _ => (0, true),
        };
        let shift = 6 - ((index % 4) * 2);
        packed[index / 4] |= code << shift;
        if unknown {
            unknown_mask[index / 8] |= 1 << (index % 8);
        }
    }
    (packed, unknown_mask)
}

/// Decodes a packed block to its canonical uppercase representation.
pub fn unpack_bases(packed: &[u8], unknown_mask: &[u8], base_length: u64) -> JmaResult<Vec<u8>> {
    let base_length = usize::try_from(base_length).map_err(|_| JmaError::OffsetOverflow)?;
    let expected_packed = base_length.div_ceil(4);
    let expected_mask = base_length.div_ceil(8);
    if packed.len() != expected_packed || unknown_mask.len() != expected_mask {
        return Err(JmaError::CorruptSection(format!(
            "packed sequence lengths {}, {} do not match base length {base_length}",
            packed.len(),
            unknown_mask.len()
        )));
    }
    let mut sequence = Vec::with_capacity(base_length);
    for index in 0..base_length {
        if unknown_mask[index / 8] & (1 << (index % 8)) != 0 {
            sequence.push(b'N');
            continue;
        }
        let shift = 6 - ((index % 4) * 2);
        let code = (packed[index / 4] >> shift) & 0b11;
        sequence.push(match code {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        });
    }
    Ok(sequence)
}

/// Reads a half-open range from blocks without allocating data outside the
/// requested range. Blocks may be adjacent or overlapping; each base is
/// copied at most once and gaps are rejected as corrupt archive data.
pub fn decode_range(
    blocks: &[SequenceBlock],
    contigs: &[ContigMetadata],
    contig_id: ContigId,
    range: SequenceRange,
) -> JmaResult<Vec<u8>> {
    let metadata = find_contig(contigs, contig_id)?;
    if range.start > range.end || range.end > metadata.length {
        return Err(JmaError::InvalidSequenceRange {
            start: range.start,
            end: range.end,
        });
    }
    let wanted = usize::try_from(range.len()).map_err(|_| JmaError::OffsetOverflow)?;
    let mut output = vec![b'N'; wanted];
    let mut covered = vec![false; wanted];
    for block in blocks.iter().filter(|block| block.contig_id == contig_id) {
        let block_end = block.end().ok_or(JmaError::OffsetOverflow)?;
        let overlap_start = block.start.max(range.start);
        let overlap_end = block_end.min(range.end);
        if overlap_start >= overlap_end {
            continue;
        }
        let sequence = unpack_bases(&block.packed, &block.unknown_mask, block.base_length)?;
        for position in overlap_start..overlap_end {
            let source =
                usize::try_from(position - block.start).map_err(|_| JmaError::OffsetOverflow)?;
            let target =
                usize::try_from(position - range.start).map_err(|_| JmaError::OffsetOverflow)?;
            if covered[target] && output[target] != sequence[source] {
                return Err(JmaError::CorruptSection(format!(
                    "overlapping sequence blocks disagree at contig {contig_id}:{position}"
                )));
            }
            output[target] = sequence[source];
            covered[target] = true;
        }
    }
    if covered.iter().any(|covered| !covered) {
        return Err(JmaError::CorruptSection(format!(
            "sequence range {}:{} is not covered by archive blocks",
            contig_id, range.start
        )));
    }
    Ok(output)
}

fn validate_block_shape(block: &SequenceBlock) -> JmaResult<()> {
    let base_length = usize::try_from(block.base_length).map_err(|_| JmaError::OffsetOverflow)?;
    if block.packed.len() != base_length.div_ceil(4)
        || block.unknown_mask.len() != base_length.div_ceil(8)
    {
        return Err(JmaError::CorruptSection(format!(
            "sequence block packed lengths do not match {} bases",
            block.base_length
        )));
    }
    block.end().ok_or(JmaError::OffsetOverflow)?;
    Ok(())
}
