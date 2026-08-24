use jam_rs::archive::{SeedKey, SeedSchemeId, SequenceRequest, TraceArchive};
use jam_rs::jam_index::archive::{K21_SCHEME, K31_SCHEME};
use jam_rs::jam_index::{
    ExternalPartReader, JamIndexArchive, MetagenomeSource, ScreenSelectionPolicy,
    write_external_part,
};
use needletail::Sequence;
use std::fs;
use tempfile::Builder;

fn fixture(sequence: &[u8]) -> (tempfile::TempDir, ExternalPartReader) {
    let directory = Builder::new()
        .prefix("jam-index-archive-")
        .tempdir_in("target")
        .unwrap();
    let source = directory.path().join("source.fasta");
    fs::write(
        &source,
        format!(">selected\n{}\n", String::from_utf8_lossy(sequence)),
    )
    .unwrap();
    let part = directory.path().join("part.bin");
    write_external_part(
        &part,
        &[MetagenomeSource {
            metagenome_id: "target".to_string(),
            sequence_path: source,
        }],
        &ScreenSelectionPolicy::default_signatures(),
    )
    .unwrap();
    let reader = ExternalPartReader::open(part).unwrap();
    (directory, reader)
}

fn seed_key(sequence: &[u8], k: u8) -> SeedKey {
    let normalized = sequence.normalize(false);
    let (_, kmer, _) = normalized.bit_kmers(k, true).next().unwrap();
    let width = usize::from(k).div_ceil(4);
    let bytes = kmer.0.to_be_bytes();
    SeedKey {
        digest: jam_rs::jamhash_u64_v1(kmer.0),
        verification: bytes[bytes.len() - width..].to_vec(),
    }
}

#[test]
fn dense_lookup_exact() {
    let sequence = b"ACGTTGCATGTCAGTACGATCGTACGTTAGCTAGCTGACTGATCGTAGCTAGTCGATCGTACGTAGCTAGTCGATGCTAGCTGATCGTAGCTAGTCGATGCTAGCTGATCGTAGCTAGTCGATGCTA";
    let (_directory, reader) = fixture(sequence);
    let archive = JamIndexArchive::load(&reader, 0, [0]).unwrap();
    for (scheme, k) in [(K21_SCHEME, 21), (K31_SCHEME, 31)] {
        let result = archive
            .lookup_seeds(SeedSchemeId(scheme), &[seed_key(sequence, k)])
            .unwrap();
        assert_eq!(result.matches.len(), 1);
        assert_eq!(result.matches[0].occurrences[0].position, 0);
        assert_eq!(result.matches[0].occurrences[0].span, u16::from(k));
    }
    let mut collision = seed_key(sequence, 21);
    collision.verification[0] ^= 1;
    assert!(
        archive
            .lookup_seeds(SeedSchemeId(K21_SCHEME), &[collision])
            .unwrap()
            .matches
            .is_empty()
    );
}

#[test]
fn bounds_apply_early() {
    let sequence = vec![b'A'; 500];
    let (_directory, reader) = fixture(&sequence);
    let archive = JamIndexArchive::load(&reader, 0, [0]).unwrap();
    let result = archive
        .lookup_seeds_bounded(
            SeedSchemeId(K21_SCHEME),
            &[seed_key(&sequence, 21)],
            Some(2),
        )
        .unwrap();
    assert!(result.matches.is_empty());
    assert!(result.metrics.occurrences_before_limits > 2);
    assert_eq!(result.metrics.occurrences_after_limits, 0);
}

#[test]
fn complete_sequence_reads() {
    let sequence = b"ACGTNNRYACGT";
    let (_directory, reader) = fixture(sequence);
    let archive = JamIndexArchive::load(&reader, 0, [0]).unwrap();
    let forward = archive
        .read_sequences(&[SequenceRequest::new(0, 0, sequence.len() as u64, false).unwrap()])
        .unwrap();
    assert_eq!(forward[0].bases, sequence);
    let reverse = archive
        .read_sequences(&[SequenceRequest::new(0, 0, sequence.len() as u64, true).unwrap()])
        .unwrap();
    assert_eq!(reverse[0].bases, b"ACGTRYNNACGT");
}
