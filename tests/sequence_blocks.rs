use jam_rs::sequence::{
    BlockCodec, DEFAULT_MAX_DECODED_BLOCK_BYTES, GearTable, SequenceBlockDirectory,
    SequenceBlockPolicy, SequenceError, decode_ambiguity_payload, decode_block_range,
    decode_block_reverse_complement_range, decode_range, decode_reverse_complement_range,
    decode_stored_block, decode_stored_block_bounded, encode_ambiguity_payload, encode_contig,
    encode_sequence_block, split_contig,
};

#[test]
fn full_iupac_payload_and_direct_ranges_are_lossless() {
    let sequence = b"AACGTRYKMBVDHN-tt";
    let encoded = encode_contig(sequence).unwrap();
    assert_eq!(decode_range(&encoded, 2..15).unwrap(), b"CGTRYKMBVDHN-");
    assert_eq!(
        decode_reverse_complement_range(&encoded, 2..15).unwrap(),
        b"-NDHBVKMRYACG"
    );
    let payload = encode_ambiguity_payload(&encoded.ambiguities).unwrap();
    assert_eq!(
        decode_ambiguity_payload(&payload).unwrap(),
        encoded.ambiguities
    );
}

#[test]
fn fixed_block_sizes_cover_exactly_and_never_cross_a_contig() {
    let sequence = b"ACGTN".repeat(300_000);
    for block_bases in [16 * 1024, 64 * 1024, 256 * 1024, 1024 * 1024] {
        let spans = split_contig(&sequence, SequenceBlockPolicy::Fixed { block_bases }).unwrap();
        assert_eq!(
            spans.iter().map(|span| span.base_count).sum::<u64>(),
            sequence.len() as u64
        );
        assert!(spans.windows(2).all(|window| {
            window[0].start + window[0].base_count == window[1].start
                && window[0].base_count <= block_bases as u64
        }));
    }
}

#[test]
fn each_gear_table_is_deterministic_and_respects_maximum() {
    let sequence = b"ACGT".repeat(100_000);
    for (min_bases, target_bases, max_bases) in [
        (1024, 2048, 4096),
        (4096, 8192, 16384),
        (16384, 32768, 65536),
    ] {
        for table in [
            GearTable::SingleBase,
            GearTable::Dinucleotide,
            GearTable::PackedFourBase,
        ] {
            let policy = SequenceBlockPolicy::Gear {
                min_bases,
                target_bases,
                max_bases,
                table,
            };
            let first = split_contig(&sequence, policy).unwrap();
            assert_eq!(first, split_contig(&sequence, policy).unwrap());
            assert!(first.iter().all(|span| span.base_count <= max_bases as u64));
            assert_eq!(
                first.iter().map(|span| span.base_count).sum::<u64>(),
                sequence.len() as u64
            );
        }
    }
}

#[test]
fn raw_and_zstd_blocks_decode_only_requested_ranges_and_bind_checksum() {
    let sequence = b"ACGTRYKMBVDHN-ACGT";
    for codec in [BlockCodec::Raw2Bit, BlockCodec::Zstd2Bit] {
        let block = encode_sequence_block(7, 100, sequence, codec).unwrap();
        let record = block.record.with_offsets(4096, 8192);
        assert_eq!(
            decode_block_range(&record, &block.payload, &block.ambiguity_payload, 103..111)
                .unwrap(),
            b"TRYKMBVD"
        );
        assert_eq!(
            decode_block_reverse_complement_range(
                &record,
                &block.payload,
                &block.ambiguity_payload,
                103..111,
            )
            .unwrap(),
            b"HBVKMRYA"
        );
        assert_eq!(
            decode_stored_block(&record, &block.payload, &block.ambiguity_payload).unwrap(),
            sequence
                .iter()
                .map(|base| base.to_ascii_uppercase())
                .collect::<Vec<_>>()
        );

        let mut corrupt = block.payload.clone();
        corrupt[0] ^= 1;
        assert!(decode_stored_block(&record, &corrupt, &block.ambiguity_payload).is_err());

        let mut flagged = record;
        flagged.flags = 1;
        assert_eq!(
            decode_stored_block(&flagged, &block.payload, &block.ambiguity_payload),
            Err(SequenceError::UnsupportedBlockFlags(1))
        );
    }
}

#[test]
fn zstd_decode_rejects_declared_lengths_before_allocation() {
    let block = encode_sequence_block(8, 0, b"ACGT", BlockCodec::Zstd2Bit).unwrap();
    assert!(matches!(
        decode_stored_block_bounded(&block.record, &block.payload, &block.ambiguity_payload, 0),
        Err(SequenceError::DecodedLengthLimit {
            requested: 1,
            maximum: 0
        })
    ));

    let mut corrupt = block.record;
    let default_limit = u64::try_from(DEFAULT_MAX_DECODED_BLOCK_BYTES).unwrap();
    corrupt.base_count = default_limit
        .checked_add(1)
        .unwrap()
        .checked_mul(4)
        .unwrap();
    corrupt.decoded_length = default_limit + 1;
    assert!(matches!(
        decode_stored_block(&corrupt, &block.payload, &block.ambiguity_payload),
        Err(SequenceError::DecodedLengthLimit { .. })
    ));
}

#[test]
fn directory_records_are_fixed_width_sorted_and_range_checked() {
    let first = encode_sequence_block(1, 0, b"ACGT", BlockCodec::Raw2Bit)
        .unwrap()
        .record
        .with_offsets(256, 512);
    let second = encode_sequence_block(1, 4, b"NNNN", BlockCodec::Raw2Bit)
        .unwrap()
        .record
        .with_offsets(260, 516);
    let directory = SequenceBlockDirectory::new(vec![second, first]).unwrap();
    let bytes = directory.encode().unwrap();
    assert_eq!(
        bytes.len(),
        SequenceBlockDirectory::HEADER_SIZE + 2 * SequenceBlockDirectory::RECORD_SIZE
    );
    assert_eq!(SequenceBlockDirectory::decode(&bytes).unwrap(), directory);
    let object_size = directory
        .records
        .iter()
        .map(|record| record.ambiguity_end().unwrap())
        .max()
        .unwrap();
    assert!(directory.validate_object_size(object_size).is_ok());
    assert!(directory.validate_object_size(object_size - 1).is_err());

    let mut corrupt = bytes;
    corrupt[6] = 0;
    assert!(SequenceBlockDirectory::decode(&corrupt).is_err());

    let mut flagged = directory.encode().unwrap();
    flagged[SequenceBlockDirectory::HEADER_SIZE + 5] = 1;
    assert_eq!(
        SequenceBlockDirectory::decode(&flagged),
        Err(SequenceError::UnsupportedBlockFlags(1))
    );
}

#[test]
fn invalid_policy_and_ranges_fail_closed() {
    assert!(split_contig(b"ACGT", SequenceBlockPolicy::Fixed { block_bases: 0 }).is_err());
    assert!(
        split_contig(
            b"ACGT",
            SequenceBlockPolicy::Gear {
                min_bases: 8,
                target_bases: 4,
                max_bases: 16,
                table: GearTable::SingleBase,
            }
        )
        .is_err()
    );
    let encoded = encode_contig(b"ACGT").unwrap();
    assert!(decode_range(&encoded, 3..5).is_err());
}
