use jam_rs::jam_index::{JamIndexPartReader, MetagenomeSource, ScreenSelectionPolicy, write_part};
use std::fs;
use tempfile::Builder;

fn directory() -> tempfile::TempDir {
    Builder::new()
        .prefix("jam-index-part-")
        .tempdir_in("target")
        .unwrap()
}

fn fixture() -> (tempfile::TempDir, Vec<MetagenomeSource>) {
    let directory = directory();
    let first = directory.path().join("first.fasta");
    let second = directory.path().join("second.fasta");
    fs::write(
        &first,
        ">short\nACGTTGCATGTCAGTACGATCGTACGTTAGCTAGCTGACTGATCGTAGCTAGTCGATCGTACGTAGCTAGTCGATGCTAGCTGATCGTAGCTAGTCGATGCTAGCTGATCGTAGCTAGTCGATGCTAGCTGATCGTAGCTAGTCGATGCTA\n>ambiguous\nACGTNNRYACGT\n",
    )
    .unwrap();
    fs::write(&second, ">reverse\nTAGCATCGACTAGCTACGATCAGCTAGCATCGACTAGCTACGATCAGCTAGCATCGACTAGCTACGATCAGCTAGCATCGACTAGCTACGATCAGCTAGCATCGACTAGCTACGATCAGCTAGCATCGACTAGCTACGATCAGCTAGCATCGT\n")
        .unwrap();
    let sources = vec![
        MetagenomeSource {
            metagenome_id: "mg-a".to_string(),
            sequence_path: first,
        },
        MetagenomeSource {
            metagenome_id: "mg-b".to_string(),
            sequence_path: second,
        },
    ];
    (directory, sources)
}

#[test]
fn external_part_roundtrip() {
    let (directory, sources) = fixture();
    let output = directory.path().join("part.bin");
    let result = write_part(
        &output,
        &sources,
        &ScreenSelectionPolicy::default_signatures(),
    )
    .unwrap();
    assert_eq!(result.metagenome_count, 2);
    assert_eq!(result.contig_count, 3);
    assert_eq!(result.source_reference_bytes, 2 * 160);
    assert_eq!(result.metagenome_directory_bytes, 2 * 48);
    assert_eq!(result.contig_length_bytes, 3 * 4);
    assert_eq!(result.exceptional_length_bytes, 0);
    assert_eq!(result.screen_samples.len(), 2);
    let reader = JamIndexPartReader::open(&output).unwrap();
    assert_eq!(reader.metagenomes()[0].metagenome_id, "mg-a");
    assert_eq!(
        reader.metagenomes()[0].source_path,
        sources[0].sequence_path
    );
    let loaded = reader.read_contigs(0, &[1]).unwrap();
    assert_eq!(loaded.contigs[&1].name, "ambiguous");
    assert_eq!(loaded.contigs[&1].bases.as_ref(), b"ACGTNNRYACGT");
    assert_eq!(reader.contig_length(0, 1).unwrap(), 12);
    assert_eq!(reader.metagenome_contigs(0).unwrap(), 0..2);
    assert_eq!(reader.posting_count(), result.posting_count);
    let bytes = fs::read(&output).unwrap();
    assert!(
        !bytes
            .windows(b"short".len())
            .any(|window| window == b"short")
    );
    assert!(
        !bytes
            .windows(b"ambiguous".len())
            .any(|window| window == b"ambiguous")
    );
    assert_eq!(
        result.single_contig_mappings + result.overflow_mappings,
        result.posting_count
    );
    assert!((0..reader.posting_count()).all(|ordinal| {
        let kind = reader.posting_kind(ordinal).unwrap();
        kind.spatial || kind.whole_sample
    }));
    assert!(
        (0..reader.posting_count()).all(|ordinal| !reader.posting(ordinal).unwrap().is_empty())
    );
}

#[test]
fn duplicate_contigs_use_compact_overflow_mappings() {
    let directory = directory();
    let source = directory.path().join("duplicates.fasta");
    let sequence = "ACGTTGCA".repeat(40);
    fs::write(
        &source,
        format!(">first\n{sequence}\n>second\n{sequence}\n"),
    )
    .unwrap();
    let output = directory.path().join("compact-part.bin");
    let result = write_part(
        &output,
        &[MetagenomeSource {
            metagenome_id: "duplicates".to_string(),
            sequence_path: source,
        }],
        &ScreenSelectionPolicy::spatial_256(512),
    )
    .unwrap();
    assert!(result.overflow_mappings > 0);
    assert!(result.maximum_overflow_count >= 2);
    assert!(result.contig_posting_bytes < result.posting_count.saturating_mul(16));
    let reader = JamIndexPartReader::open(output).unwrap();
    assert!((0..reader.posting_count()).any(|ordinal| {
        reader
            .posting(ordinal)
            .is_ok_and(|posting| posting.len() == 2)
    }));
}

#[test]
fn external_corruption_fails() {
    let (directory, sources) = fixture();
    let output = directory.path().join("part.bin");
    write_part(
        &output,
        &sources,
        &ScreenSelectionPolicy::default_signatures(),
    )
    .unwrap();
    let original = fs::read(&output).unwrap();

    let mut posting_corrupt = original;
    let posting_offset = u64::from_le_bytes(posting_corrupt[64..72].try_into().unwrap());
    posting_corrupt[posting_offset as usize] ^= 1;
    let posting_path = directory.path().join("posting-corrupt.bin");
    fs::write(&posting_path, posting_corrupt).unwrap();
    assert!(JamIndexPartReader::open(posting_path).is_err());

    let reader = JamIndexPartReader::open(output).unwrap();
    let changed = fs::read_to_string(&sources[0].sequence_path)
        .unwrap()
        .replace("ACGTNNRYACGT", "ACGTNNRAAACGT");
    fs::write(&sources[0].sequence_path, changed).unwrap();
    assert!(reader.read_contigs(0, &[1]).is_err());
}
