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
fn independent_part_roundtrips_directories_signatures_and_complete_sequences() {
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
    assert!(result.packed_sequence_bytes > 0);
    assert_eq!(result.screen_samples.len(), 2);
    let reader = JamIndexPartReader::open(&output).unwrap();
    reader.verify_sequence_section().unwrap();
    assert_eq!(reader.metagenomes()[0].metagenome_id, "mg-a");
    assert_eq!(reader.contigs()[0].name, "short");
    assert_eq!(reader.read_contig(1).unwrap(), b"ACGTNNRYACGT");
    assert_eq!(reader.contig_ids_for_metagenome(0).unwrap(), 0..2);
    let shared = result.screen_samples[0]
        .hashes
        .iter()
        .find_map(|hash| {
            let hits = reader.signature_hits(*hash).unwrap();
            (!hits.is_empty()).then_some(hits)
        })
        .expect("screen hash maps to a contig");
    assert!(shared.iter().any(|hit| hit.contig_id < 2));
}

#[test]
fn signature_corruption_fails_open_and_sequence_corruption_fails_selected_read() {
    let (directory, sources) = fixture();
    let output = directory.path().join("part.bin");
    write_part(
        &output,
        &sources,
        &ScreenSelectionPolicy::default_signatures(),
    )
    .unwrap();
    let original = fs::read(&output).unwrap();

    let mut signature_corrupt = original.clone();
    let signature_offset = u64::from_le_bytes(signature_corrupt[64..72].try_into().unwrap());
    signature_corrupt[signature_offset as usize] ^= 1;
    let signature_path = directory.path().join("signature-corrupt.bin");
    fs::write(&signature_path, signature_corrupt).unwrap();
    assert!(JamIndexPartReader::open(signature_path).is_err());

    let mut sequence_corrupt = original;
    sequence_corrupt[512] ^= 1;
    let sequence_path = directory.path().join("sequence-corrupt.bin");
    fs::write(&sequence_path, sequence_corrupt).unwrap();
    let reader = JamIndexPartReader::open(sequence_path).unwrap();
    assert!(reader.read_contig(0).is_err());
}
