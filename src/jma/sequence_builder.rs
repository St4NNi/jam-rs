//! Sequence block construction for JMA v1 archives.
//!
//! The builder keeps the format-facing sequence code deliberately small: a
//! contig is divided into deterministic, sorted blocks and each block is
//! packed with the checked encoder used by the JMA reader.  Ambiguous bases
//! are retained in the block's mask rather than being silently removed.

use crate::jma::sequence::{SequenceBlock, pack_bases};
use crate::jma::{ContigId, JmaError, JmaResult};

/// Default maximum number of bases in one packed sequence block.
pub const DEFAULT_BLOCK_BASES: usize = 1 << 20;

/// Configuration for packed sequence construction.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SequenceBuildConfig {
    /// Maximum number of bases in each block.  A value of zero is invalid.
    pub block_bases: usize,
}

impl Default for SequenceBuildConfig {
    fn default() -> Self {
        Self {
            block_bases: DEFAULT_BLOCK_BASES,
        }
    }
}

/// Builds sorted packed blocks covering the complete input sequence.
pub fn build_sequence_blocks(
    contig_id: ContigId,
    sequence: &[u8],
    config: SequenceBuildConfig,
) -> JmaResult<Vec<SequenceBlock>> {
    if config.block_bases == 0 {
        return Err(JmaError::CorruptSection(
            "sequence block size must be greater than zero".to_string(),
        ));
    }

    let mut blocks = Vec::new();
    for (block_index, chunk) in sequence.chunks(config.block_bases).enumerate() {
        let start = block_index
            .checked_mul(config.block_bases)
            .and_then(|value| u64::try_from(value).ok())
            .ok_or(JmaError::OffsetOverflow)?;
        let base_length = u64::try_from(chunk.len()).map_err(|_| JmaError::OffsetOverflow)?;
        let (packed, unknown_mask) = pack_bases(chunk);
        blocks.push(SequenceBlock {
            contig_id,
            start,
            base_length,
            packed,
            unknown_mask,
        });
    }
    Ok(blocks)
}

/// Alias used by callers that refer to packed blocks simply as blocks.
pub fn build_blocks(
    contig_id: ContigId,
    sequence: &[u8],
    block_bases: usize,
) -> JmaResult<Vec<SequenceBlock>> {
    build_sequence_blocks(contig_id, sequence, SequenceBuildConfig { block_bases })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jma::sequence::unpack_bases;

    #[test]
    fn blocks_cover_sequence_and_preserve_ambiguity() {
        let sequence = b"acgtnuXG";
        let blocks =
            build_sequence_blocks(4, sequence, SequenceBuildConfig { block_bases: 3 }).unwrap();

        assert_eq!(
            blocks
                .iter()
                .map(|block| block.base_length)
                .collect::<Vec<_>>(),
            vec![3, 3, 2]
        );
        assert_eq!(blocks[0].start, 0);
        assert_eq!(blocks[1].start, 3);
        assert_eq!(blocks[2].start, 6);

        let decoded = blocks
            .iter()
            .flat_map(|block| {
                unpack_bases(&block.packed, &block.unknown_mask, block.base_length).unwrap()
            })
            .collect::<Vec<_>>();
        assert_eq!(decoded, b"ACGTNTNG");
    }

    #[test]
    fn zero_block_size_is_rejected() {
        let error =
            build_sequence_blocks(0, b"ACGT", SequenceBuildConfig { block_bases: 0 }).unwrap_err();
        assert!(error.to_string().contains("block size"));
    }

    #[test]
    fn empty_sequence_has_no_blocks() {
        let blocks = build_sequence_blocks(0, b"", SequenceBuildConfig::default()).unwrap();
        assert!(blocks.is_empty());
    }
}
