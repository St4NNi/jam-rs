//! JMA sequence-block construction over the shared sequence codec.

use crate::jma::{ContigId, JmaError, JmaResult};
use crate::sequence::{BlockCodec, EncodedSequenceBlock, SequenceBlockPolicy, SequenceSpan};

/// Default fixed sequence block size used by indexed archive construction.
pub const DEFAULT_BLOCK_BASES: usize = 1 << 20;

/// Configuration for final independent sequence blocks.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct IndexedSequenceBuildConfig {
    pub policy: SequenceBlockPolicy,
    pub codec: BlockCodec,
}

impl Default for IndexedSequenceBuildConfig {
    fn default() -> Self {
        Self {
            policy: SequenceBlockPolicy::Fixed {
                block_bases: DEFAULT_BLOCK_BASES,
            },
            codec: BlockCodec::Raw2Bit,
        }
    }
}

/// Plans deterministic block spans without encoding sequence payloads.
pub fn plan_indexed_sequence_blocks(
    sequence: &[u8],
    policy: SequenceBlockPolicy,
) -> JmaResult<Vec<SequenceSpan>> {
    crate::sequence::split_contig(sequence, policy)
        .map_err(|error| JmaError::CorruptSection(error.to_string()))
}

/// Builds complete, independently readable blocks for one contig. No block
/// crosses a contig boundary because this function accepts exactly one contig
/// and binds the contig ID into every checksum.
pub fn build_indexed_sequence_blocks(
    contig_id: ContigId,
    sequence: &[u8],
    config: IndexedSequenceBuildConfig,
) -> JmaResult<Vec<EncodedSequenceBlock>> {
    let spans = plan_indexed_sequence_blocks(sequence, config.policy)?;
    spans
        .into_iter()
        .map(|span| {
            let start = usize::try_from(span.start).map_err(|_| JmaError::OffsetOverflow)?;
            let end = span
                .start
                .checked_add(span.base_count)
                .ok_or(JmaError::OffsetOverflow)?;
            let end = usize::try_from(end).map_err(|_| JmaError::OffsetOverflow)?;
            let bases = sequence
                .get(start..end)
                .ok_or_else(|| JmaError::CorruptSection("sequence block exceeds contig".into()))?;
            crate::sequence::encode_sequence_block(contig_id, span.start, bases, config.codec)
                .map_err(|error| JmaError::CorruptSection(error.to_string()))
        })
        .collect()
}
