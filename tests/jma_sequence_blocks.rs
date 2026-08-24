use jam_rs::jma::sequence::{
    BlockCodec, SequenceBlockDirectory, SequenceBlockPolicy, build_indexed_sequence_blocks,
    decode_indexed_sequence_directory, encode_indexed_sequence_directory,
};
use jam_rs::jma::sequence_builder::{
    IndexedSequenceBuildConfig, build_indexed_sequence_blocks as build_from_config,
};

#[test]
fn jma_wrapper_builds_contig_local_blocks_and_directory() {
    let sequence = b"ACGTRYKMBVDHN-".repeat(3000);
    let blocks = build_indexed_sequence_blocks(
        11,
        &sequence,
        SequenceBlockPolicy::Fixed { block_bases: 4096 },
        BlockCodec::Raw2Bit,
    )
    .unwrap();
    assert!(blocks.iter().all(|block| block.record.contig_id == 11));
    assert_eq!(
        blocks
            .iter()
            .map(|block| block.record.base_count)
            .sum::<u64>(),
        sequence.len() as u64
    );

    let records = blocks.iter().map(|block| block.record).collect::<Vec<_>>();
    let directory_bytes = encode_indexed_sequence_directory(&records).unwrap();
    let decoded = decode_indexed_sequence_directory(&directory_bytes).unwrap();
    assert_eq!(decoded, SequenceBlockDirectory::new(records).unwrap());
}

#[test]
fn jma_builder_supports_gear_and_zstd_policies() {
    let sequence = b"ACGT".repeat(50_000);
    let blocks = build_from_config(
        2,
        &sequence,
        IndexedSequenceBuildConfig {
            policy: SequenceBlockPolicy::Gear {
                min_bases: 1024,
                target_bases: 2048,
                max_bases: 4096,
                table: jam_rs::sequence::GearTable::Dinucleotide,
            },
            codec: BlockCodec::Zstd2Bit,
        },
    )
    .unwrap();
    assert!(!blocks.is_empty());
    assert!(blocks.iter().all(|block| block.record.contig_id == 2));
    assert!(
        blocks
            .iter()
            .all(|block| block.record.codec == BlockCodec::Zstd2Bit)
    );
    assert!(blocks.iter().all(|block| block.record.base_count <= 4096));
}
