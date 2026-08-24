//! Final JMA sequence-block wrappers.
//!
//! Sequence payloads and directory records are defined by the shared
//! backend-neutral codec in [`crate::sequence`]. This module only exposes the
//! JMA-facing names used by the format builder and reader; it contains no
//! separate sequence layout or compatibility decoder.

use crate::jma::{ContigId, JmaError, JmaResult};

pub use crate::sequence::{
    BlockCodec, BlockDecodeError, DEFAULT_MAX_DECODED_BLOCK_BYTES, EncodedSequenceBlock, GearTable,
    SequenceBlockDirectory, SequenceBlockDirectoryError, SequenceBlockPolicy, SequenceBlockRecord,
    SequenceSpan, decode_block_range, decode_block_range_bounded,
    decode_block_reverse_complement_range, decode_stored_block, decode_stored_block_bounded,
    encode_sequence_block, sequence_block_checksum, split_contig,
};

/// Builds checksum-bound independent blocks for one contig. Every returned
/// block is local to `contig_id`; the archive writer assigns absolute payload
/// offsets before encoding the final directory.
pub fn build_indexed_sequence_blocks(
    contig_id: ContigId,
    sequence: &[u8],
    policy: SequenceBlockPolicy,
    codec: BlockCodec,
) -> JmaResult<Vec<EncodedSequenceBlock>> {
    let spans = crate::sequence::split_contig(sequence, policy)
        .map_err(|error| JmaError::CorruptSection(error.to_string()))?;
    spans
        .into_iter()
        .map(|span| {
            let start = usize::try_from(span.start).map_err(|_| JmaError::OffsetOverflow)?;
            let count = usize::try_from(span.base_count).map_err(|_| JmaError::OffsetOverflow)?;
            let end = start.checked_add(count).ok_or(JmaError::OffsetOverflow)?;
            let bases = sequence
                .get(start..end)
                .ok_or_else(|| JmaError::CorruptSection("sequence span exceeds contig".into()))?;
            crate::sequence::encode_sequence_block(contig_id, span.start, bases, codec)
                .map_err(|error| JmaError::CorruptSection(error.to_string()))
        })
        .collect()
}

/// Encodes the fixed-width sequence-block directory after the archive writer
/// has assigned absolute payload and ambiguity offsets.
pub fn encode_indexed_sequence_directory(records: &[SequenceBlockRecord]) -> JmaResult<Vec<u8>> {
    SequenceBlockDirectory::new(records.to_vec())
        .and_then(|directory| directory.encode())
        .map_err(|error| JmaError::CorruptSection(error.to_string()))
}

/// Decodes and validates a self-contained sequence-block directory.
pub fn decode_indexed_sequence_directory(bytes: &[u8]) -> JmaResult<SequenceBlockDirectory> {
    SequenceBlockDirectory::decode(bytes)
        .map_err(|error| JmaError::CorruptSection(error.to_string()))
}
