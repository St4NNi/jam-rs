//! Deterministic malformed-input and boundary properties for the trace data
//! formats.  These tests deliberately mutate encoded bytes rather than
//! deriving expected answers from a current trace result.

use jam_rs::jma::contigs::decode_contigs;
use jam_rs::jma::format::{
    HEADER_SIZE, SECTION_ENTRY_SIZE, checksum, decode_header, decode_section_directory,
};
use jam_rs::jma::header::{parse_header, parse_section_directory};
use jam_rs::jma::sequence::{SequenceBlock, decode_sequence_blocks, pack_bases, unpack_bases};
use jam_rs::jma::writer::{
    ArchiveParts, SeedLevelData, SeedRecord, SeedSection, decode_seed_section, encode_archive,
    encode_seed_section,
};
use jam_rs::jma::{ContigMetadata, SeedLevel, SeedOccurrence, SeedQuery};
use jam_rs::resource::ResourceMetrics;
use jam_rs::trace::coverage::{
    cigar_from_edit_script, parse_cigar, project_cigar, project_edit_script,
};
use jam_rs::trace::intervals::{
    circular_gap_complement, circular_union, covered_length, split_circular, union,
};
use jam_rs::trace::model::{BaseInterval, EditOperation, EditRun, TraceRecord, TraceRunFooter};
use jam_rs::trace::output::SCHEMA_VERSION;
use std::panic::{AssertUnwindSafe, catch_unwind};

fn no_panic<T>(label: &str, f: impl FnOnce() -> T) -> T {
    catch_unwind(AssertUnwindSafe(f))
        .unwrap_or_else(|_| panic!("{label} panicked on malformed or generated input"))
}

fn put_u32(bytes: &mut [u8], offset: usize, value: u32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

fn put_u64(bytes: &mut [u8], offset: usize, value: u64) {
    bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

fn get_u64(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(bytes[offset..offset + 8].try_into().unwrap())
}

fn refresh_header_checksum(bytes: &mut [u8]) {
    let mut input = bytes[..HEADER_SIZE].to_vec();
    input[80..112].fill(0);
    bytes[80..112].copy_from_slice(&checksum(&input));
}

fn archive_fixture() -> Vec<u8> {
    let sequence = b"ACGTNNACGTACGTACGTACGTACGTACGTACGT";
    let (packed, unknown_mask) = pack_bases(sequence);
    let contig = ContigMetadata {
        id: 0,
        name: "fixture-contig".to_string(),
        length: sequence.len() as u64,
    };
    encode_archive(&ArchiveParts {
        flags: 0,
        source_sha256: [7; 32],
        contigs: vec![contig],
        sequence_blocks: vec![SequenceBlock {
            contig_id: 0,
            start: 0,
            base_length: sequence.len() as u64,
            packed,
            unknown_mask,
        }],
        seed_sections: vec![
            SeedSection {
                k: 21,
                levels: vec![SeedLevelData {
                    level: SeedLevel { k: 21, scale: 500 },
                    records: vec![SeedRecord {
                        query: SeedQuery {
                            k: 21,
                            hash: 11,
                            canonical_kmer: 3,
                        },
                        occurrence: SeedOccurrence {
                            contig_id: 0,
                            position: 2,
                            reverse: false,
                        },
                    }],
                }],
            },
            SeedSection {
                k: 31,
                levels: vec![SeedLevelData {
                    level: SeedLevel { k: 31, scale: 200 },
                    records: vec![SeedRecord {
                        query: SeedQuery {
                            k: 31,
                            hash: 17,
                            canonical_kmer: 5,
                        },
                        occurrence: SeedOccurrence {
                            contig_id: 0,
                            position: 4,
                            reverse: true,
                        },
                    }],
                }],
            },
        ],
    })
    .unwrap()
}

#[test]
fn header_mutations_are_rejected_without_panics() {
    let archive = archive_fixture();
    let valid = &archive[..HEADER_SIZE];
    assert!(decode_header(valid).is_ok());
    assert!(parse_header(valid).is_ok());

    let mut fixed_header_mutations = Vec::new();
    for (offset, value) in [
        (0, 0x01),
        (4, 0xff),
        (6, 0x01),
        (8, 0x80),
        (80, 0x01),
        (112, 0x80),
    ] {
        let mut bytes = valid.to_vec();
        bytes[offset] ^= value;
        fixed_header_mutations.push(bytes);
    }

    for (index, bytes) in fixed_header_mutations.into_iter().enumerate() {
        let decoded = no_panic(&format!("fixed header mutation {index}"), || {
            decode_header(&bytes)
        });
        assert!(
            decoded.is_err(),
            "fixed header mutation {index} was accepted"
        );
        let parsed = no_panic(&format!("parsed fixed header mutation {index}"), || {
            parse_header(&bytes)
        });
        assert!(
            parsed.is_err(),
            "parsed fixed header mutation {index} was accepted"
        );
    }

    let mut structural_headers = Vec::new();
    let mut wrong_directory_length = valid.to_vec();
    put_u64(&mut wrong_directory_length, 72, 1);
    refresh_header_checksum(&mut wrong_directory_length);
    structural_headers.push(wrong_directory_length);

    let mut overflowing_directory = valid.to_vec();
    put_u64(&mut overflowing_directory, 64, u64::MAX);
    refresh_header_checksum(&mut overflowing_directory);
    structural_headers.push(overflowing_directory);

    let mut before_header = valid.to_vec();
    put_u64(&mut before_header, 64, 1);
    refresh_header_checksum(&mut before_header);
    structural_headers.push(before_header);

    let mut impossible_count = valid.to_vec();
    put_u32(&mut impossible_count, 24, u32::MAX);
    refresh_header_checksum(&mut impossible_count);
    structural_headers.push(impossible_count);

    for (index, bytes) in structural_headers.into_iter().enumerate() {
        let _decoded = no_panic(&format!("structural header mutation {index}"), || {
            decode_header(&bytes)
        });
        let parsed = no_panic(
            &format!("parsed structural header mutation {index}"),
            || parse_header(&bytes),
        );
        assert!(
            parsed.is_err(),
            "structural header mutation {index} was accepted"
        );
    }

    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    for index in 0..96 {
        let mut bytes = vec![0u8; HEADER_SIZE];
        for byte in &mut bytes {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            *byte = (state >> 56) as u8;
        }
        let result = no_panic(&format!("random header {index}"), || decode_header(&bytes));
        assert!(
            result.is_err(),
            "random malformed header {index} was accepted"
        );
    }
}

#[test]
fn section_directory_mutations_are_checked_before_allocation() {
    let archive = archive_fixture();
    let directory_offset = get_u64(&archive, 64) as usize;
    let directory_length = get_u64(&archive, 72) as usize;
    let directory = &archive[directory_offset..directory_offset + directory_length];
    let section_count = (directory_length / SECTION_ENTRY_SIZE) as u32;
    assert_eq!(
        decode_section_directory(directory, section_count)
            .unwrap()
            .len(),
        4
    );
    assert_eq!(
        parse_section_directory(
            directory,
            section_count,
            directory_offset as u64,
            archive.len() as u64,
        )
        .unwrap()
        .len(),
        4
    );

    let malformed = [
        directory[..directory.len() - 1].to_vec(),
        [directory, &[0]].concat(),
        {
            let mut bytes = directory.to_vec();
            put_u64(&mut bytes, 8, u64::MAX);
            bytes
        },
        {
            let mut bytes = directory.to_vec();
            put_u64(&mut bytes, 16, u64::MAX);
            bytes
        },
        {
            let mut bytes = directory.to_vec();
            put_u64(&mut bytes, 24, 0);
            bytes
        },
    ];

    for (index, bytes) in malformed.into_iter().enumerate() {
        let decoded = no_panic(&format!("directory mutation {index}"), || {
            decode_section_directory(&bytes, section_count)
        });
        assert!(decoded.is_err(), "directory mutation {index} was accepted");
        let parsed = no_panic(&format!("parsed directory mutation {index}"), || {
            parse_section_directory(
                &bytes,
                section_count,
                directory_offset as u64,
                archive.len() as u64,
            )
        });
        assert!(
            parsed.is_err(),
            "parsed directory mutation {index} was accepted"
        );
    }

    let mut oversized_count = directory.to_vec();
    let result = no_panic("oversized section count", || {
        decode_section_directory(&oversized_count, u32::MAX)
    });
    assert!(result.is_err());
    oversized_count.clear();
}

#[test]
fn seed_payload_roundtrip_and_malformed_lengths_are_safe() {
    let section = SeedSection {
        k: 31,
        levels: vec![SeedLevelData {
            level: SeedLevel { k: 31, scale: 100 },
            records: vec![
                SeedRecord {
                    query: SeedQuery {
                        k: 31,
                        hash: 0,
                        canonical_kmer: 1,
                    },
                    occurrence: SeedOccurrence {
                        contig_id: 0,
                        position: 4,
                        reverse: false,
                    },
                },
                SeedRecord {
                    query: SeedQuery {
                        k: 31,
                        hash: 99,
                        canonical_kmer: 2,
                    },
                    occurrence: SeedOccurrence {
                        contig_id: 3,
                        position: 8,
                        reverse: true,
                    },
                },
            ],
        }],
    };
    let encoded = encode_seed_section(&section).unwrap();
    let decoded = decode_seed_section(&encoded).unwrap();
    assert_eq!(encode_seed_section(&decoded).unwrap(), encoded);

    let mut mutations = Vec::new();
    let mut version = encoded.clone();
    put_u32(&mut version, 0, 2);
    mutations.push(version);
    let mut invalid_k = encoded.clone();
    invalid_k[4] = 0;
    mutations.push(invalid_k);
    let mut reserved = encoded.clone();
    reserved[5] = 1;
    mutations.push(reserved);
    let mut levels = encoded.clone();
    put_u32(&mut levels, 8, u32::MAX);
    mutations.push(levels);
    let mut scale = encoded.clone();
    put_u64(&mut scale, 12, 0);
    mutations.push(scale);
    let mut records = encoded.clone();
    put_u64(&mut records, 20, u64::MAX);
    mutations.push(records);
    let mut orientation = encoded.clone();
    orientation[12 + 16 + 20] = 2;
    mutations.push(orientation);
    let mut record_reserved = encoded.clone();
    record_reserved[12 + 16 + 21] = 1;
    mutations.push(record_reserved);
    let mut trailing = encoded.clone();
    trailing.push(0);
    mutations.push(trailing);

    for (index, bytes) in mutations.into_iter().enumerate() {
        let result = no_panic(&format!("seed mutation {index}"), || {
            decode_seed_section(&bytes)
        });
        assert!(result.is_err(), "seed mutation {index} was accepted");
    }

    let mut huge_record_count = vec![0u8; 28];
    put_u32(&mut huge_record_count, 0, 1);
    huge_record_count[4] = 31;
    put_u64(&mut huge_record_count, 12, 100);
    put_u64(&mut huge_record_count, 20, u64::MAX);
    let result = no_panic("huge seed record count", || {
        decode_seed_section(&huge_record_count)
    });
    assert!(result.is_err());
}

#[test]
fn packed_sequences_roundtrip_and_declared_lengths_stay_bounded() {
    for length in 0..128usize {
        let input: Vec<u8> = (0..length)
            .map(|index| match index % 7 {
                0 => b'a',
                1 => b'C',
                2 => b'G',
                3 => b't',
                4 => b'N',
                5 => b'U',
                _ => b'?',
            })
            .collect();
        let (packed, unknown_mask) = pack_bases(&input);
        let decoded = unpack_bases(&packed, &unknown_mask, length as u64).unwrap();
        let expected: Vec<u8> = input
            .into_iter()
            .map(|base| match base.to_ascii_uppercase() {
                b'A' | b'C' | b'G' => base.to_ascii_uppercase(),
                b'T' | b'U' => b'T',
                _ => b'N',
            })
            .collect();
        assert_eq!(decoded, expected);
        assert!(decoded.len() <= 128);
    }

    let malformed = [
        (vec![0], vec![0], u64::MAX),
        (vec![0; 8], vec![0; 8], 1),
        (vec![], vec![], 1),
        (vec![0; 1], vec![], 4),
        (vec![0; 1], vec![0; 1], 5),
    ];
    for (index, (packed, mask, length)) in malformed.into_iter().enumerate() {
        let result = no_panic(&format!("packed mutation {index}"), || {
            unpack_bases(&packed, &mask, length)
        });
        assert!(result.is_err(), "packed mutation {index} was accepted");
    }

    let mut huge_block = vec![0u8; 44];
    put_u32(&mut huge_block, 0, 1);
    put_u32(&mut huge_block, 8, 0);
    put_u64(&mut huge_block, 12, 0);
    put_u64(&mut huge_block, 20, u64::MAX);
    put_u64(&mut huge_block, 28, 0);
    put_u64(&mut huge_block, 36, 0);
    let result = no_panic("huge sequence block", || {
        decode_sequence_blocks(&huge_block)
    });
    assert!(result.is_err());

    let mut state = 0x243f_6a88_85a3_08d3u64;
    for index in 0..128 {
        let length = (index * 13) % 384;
        let mut bytes = vec![0u8; length];
        for byte in &mut bytes {
            state = state
                .wrapping_mul(2_862_933_555_777_941_757)
                .wrapping_add(3_037_000_493);
            *byte = (state >> 32) as u8;
        }
        let result = no_panic(&format!("random sequence section {index}"), || {
            decode_sequence_blocks(&bytes)
        });
        if let Ok(blocks) = result {
            assert!(blocks.len() <= bytes.len());
            for block in blocks {
                assert!(block.packed.len() <= bytes.len());
                assert!(block.unknown_mask.len() <= bytes.len());
            }
        }
    }

    let malformed_contigs = {
        let mut bytes = vec![0u8; 24];
        put_u32(&mut bytes, 0, 1);
        put_u32(&mut bytes, 4, 1);
        put_u32(&mut bytes, 8, 0);
        put_u64(&mut bytes, 12, 1);
        put_u32(&mut bytes, 20, u32::MAX);
        bytes
    };
    let result = no_panic("huge contig name length", || {
        decode_contigs(&malformed_contigs, 1)
    });
    assert!(result.is_err());
}

#[test]
fn cigar_and_edit_script_properties_are_coordinate_safe() {
    let valid = ["3=2X1I4D5S", "10M", "1=1X1I1D1S", "4294967295="];
    for cigar in valid {
        let runs = no_panic(cigar, || parse_cigar(cigar)).unwrap();
        let canonical = cigar_from_edit_script(&runs);
        assert_eq!(parse_cigar(&canonical).unwrap(), runs);
    }

    let malformed = [
        "",
        "3",
        "0=",
        "=3",
        "1Q",
        "4294967296=",
        "18446744073709551616=",
        "1=2",
        "1==",
    ];
    for cigar in malformed {
        let result = no_panic(cigar, || parse_cigar(cigar));
        assert!(result.is_err(), "malformed CIGAR {cigar:?} was accepted");
    }

    let segments = [
        BaseInterval { start: 8, end: 10 },
        BaseInterval { start: 0, end: 7 },
    ];
    let projection = project_cigar(&segments, "2=1I2D5=").unwrap();
    assert_eq!(projection.query_consumed, 9);
    assert_eq!(projection.target_consumed, 8);
    for span in &projection.spans {
        for interval in &span.query_segments {
            assert!(interval.end <= 10);
            assert!(interval.start <= interval.end);
        }
        assert!(span.target_interval.start <= span.target_interval.end);
    }

    let operations = [
        EditOperation::Equal,
        EditOperation::Substitution,
        EditOperation::Insertion,
        EditOperation::Deletion,
        EditOperation::SoftClip,
    ];
    for seed in 0..128u64 {
        let mut state = seed.wrapping_add(1);
        let mut runs = Vec::new();
        let mut query_length = 0u64;
        for _ in 0..8 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let operation = operations[(state as usize) % operations.len()];
            let length = (state % 5) + 1;
            if matches!(
                operation,
                EditOperation::Equal
                    | EditOperation::Substitution
                    | EditOperation::Deletion
                    | EditOperation::SoftClip
            ) {
                query_length += length;
            }
            runs.push(EditRun {
                operation,
                length: length as u32,
            });
        }
        let projection = project_edit_script(
            &[BaseInterval {
                start: 0,
                end: query_length,
            }],
            &runs,
        )
        .unwrap();
        assert_eq!(projection.query_consumed, query_length);
        let canonical_runs = parse_cigar(&cigar_from_edit_script(&runs)).unwrap();
        assert_eq!(parse_cigar(&projection.cigar).unwrap(), canonical_runs);
        for span in projection.spans {
            for interval in span.query_segments {
                assert!(interval.start <= interval.end);
                assert!(interval.end <= query_length);
            }
            assert!(span.target_interval.start <= span.target_interval.end);
        }
    }
}

#[test]
fn circular_union_and_gap_complement_preserve_bounds_and_total_length() {
    let mut state = 0x517c_c1b7_2722_0a95u64;
    for case in 0..256 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let length = 1 + state % 256;
        let mut segments = Vec::new();
        for _ in 0..(state as usize % 12) {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let start = state % (length + 1);
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let end = start + (state % (length - start + 1));
            segments.push(BaseInterval { start, end });
        }
        let covered = no_panic(&format!("circular union {case}"), || {
            circular_union(&segments, length)
        })
        .unwrap();
        assert_eq!(covered, union(segments.clone()));
        assert!(covered.windows(2).all(|pair| pair[0].end < pair[1].start));
        assert!(covered.iter().all(|interval| interval.end <= length));
        assert!(covered_length(&covered) <= length);

        let gaps = circular_gap_complement(&covered, length).unwrap();
        assert!(gaps.iter().all(|interval| interval.end <= length));
        assert_eq!(covered_length(&covered) + covered_length(&gaps), length);
        let mut all = covered.clone();
        all.extend(gaps);
        assert_eq!(
            union(all),
            vec![BaseInterval {
                start: 0,
                end: length
            }]
        );

        let wrap_start = length.saturating_sub(2);
        let wrap_end = if length > 2 { 1 } else { 0 };
        let split = split_circular(wrap_start, wrap_end, length).unwrap();
        assert!(split.iter().all(|interval| interval.end <= length));
    }

    assert!(circular_union(&[BaseInterval { start: 0, end: 2 }], 0).is_err());
    assert!(circular_union(&[BaseInterval { start: 2, end: 1 }], 3).is_err());
    assert!(circular_union(&[BaseInterval { start: 0, end: 4 }], 3).is_err());
    assert!(circular_gap_complement(&[BaseInterval { start: 0, end: 1 }], 0).is_err());
}

#[test]
fn json_serialization_is_stable_and_malformed_values_do_not_panic() {
    let footer = TraceRunFooter {
        schema_version: SCHEMA_VERSION.to_string(),
        run_id: "property-run".to_string(),
        completed_at_utc: "2026-01-01T00:00:00Z".to_string(),
        metagenomes_total: 3,
        metagenomes_with_candidates: 2,
        metagenomes_aligned: 1,
        metagenomes_failed: 1,
        alignments_total: 4,
        resource_metrics: ResourceMetrics::default(),
    };
    let record = TraceRecord::RunFooter(footer);
    let first = serde_json::to_vec(&record).unwrap();
    for _ in 0..32 {
        assert_eq!(serde_json::to_vec(&record).unwrap(), first);
    }
    let decoded: TraceRecord = serde_json::from_slice(&first).unwrap();
    assert_eq!(decoded, record);

    let malformed = [
        b"".as_slice(),
        b"{".as_slice(),
        b"null".as_slice(),
        br#"{"record_type":"unknown"}"#.as_slice(),
        br#"{"record_type":"run_footer","schema_version":null}"#.as_slice(),
        &[0xff, 0xfe, 0xfd],
    ];
    for (index, bytes) in malformed.into_iter().enumerate() {
        let result = no_panic(&format!("JSON mutation {index}"), || {
            serde_json::from_slice::<TraceRecord>(bytes)
        });
        assert!(
            result.is_err(),
            "malformed JSON mutation {index} was accepted"
        );
    }
}
