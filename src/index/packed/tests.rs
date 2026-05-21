use super::*;
use crate::index::MappedFile;
use std::path::Path;

fn random(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

/// Contigs of random bases with runs of N, lower case bases and other bytes.
fn contigs(seed: u64) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut state = seed;
    [0usize, 1, 70_000, 5, 33_000]
        .iter()
        .enumerate()
        .map(|(index, &length)| {
            let mut sequence: Vec<u8> = (0..length)
                .map(|_| b"ACGTacgt"[(random(&mut state) % 8) as usize])
                .collect();
            if length > 1000 {
                sequence[100..400].fill(b'N');
                sequence[401] = b'r';
                sequence[length - 1] = b'n';
            }
            (format!("contig_{index}").into_bytes(), sequence)
        })
        .collect()
}

fn write(path: &Path, contigs: &[(Vec<u8>, Vec<u8>)]) {
    let mut writer = PackedWriter::new();
    for (name, sequence) in contigs {
        writer.add_contig(name, sequence).unwrap();
    }
    let (head, packed) = writer.finish();
    std::fs::write(path, [head, packed].concat()).unwrap();
}

fn open(path: &Path) -> Result<PackedSequence, IndexError> {
    PackedSequence::open(Arc::new(MappedFile::open(path, false)?))
}

#[test]
fn packed_round_trip() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("entry.packed");
    let contigs = contigs(3);
    write(&path, &contigs);
    let packed = open(&path).unwrap();
    assert_eq!(packed.contig_count(), contigs.len());
    let total: usize = contigs.iter().map(|(_, sequence)| sequence.len()).sum();
    assert_eq!(packed.bases(), total as u64);
    let mut out = Vec::new();
    let mut start = 0u64;
    for (index, (name, sequence)) in contigs.iter().enumerate() {
        let expected = sequence.to_ascii_uppercase();
        assert_eq!(packed.contig_name(index).unwrap(), *name);
        let length = sequence.len() as u64;
        assert_eq!(packed.contig_range(index).unwrap(), start..start + length);
        packed.read(index, 0..length, &mut out).unwrap();
        assert_eq!(out, expected);
        if length > 20_000 {
            packed.read(index, 16_383..16_390, &mut out).unwrap();
            assert_eq!(out, expected[16_383..16_390]);
            packed.read(index, 350..402, &mut out).unwrap();
            assert_eq!(out, expected[350..402]);
            assert_eq!(packed.contig_at(start + 17), Some(index));
        }
        start += length;
    }
    assert_eq!(packed.contig_at(total as u64), None);
    assert!(packed.read(2, 0..70_001, &mut out).is_err());
    assert!(packed.contig_name(5).is_err());
}

#[test]
fn empty_entries_and_bad_names() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("empty.packed");
    write(&path, &[]);
    let packed = open(&path).unwrap();
    assert_eq!((packed.contig_count(), packed.bases()), (0, 0));
    let mut writer = PackedWriter::new();
    assert!(writer.add_contig(b"", b"ACGT").is_err());
    assert!(writer.add_contig(b"two words", b"ACGT").is_err());
}

#[test]
fn damage_is_detected() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("entry.packed");
    write(&path, &contigs(5));
    let original = std::fs::read(&path).unwrap();
    let damaged = directory.path().join("damaged.packed");
    for at in [0, 30, 70, 200] {
        let mut bytes = original.clone();
        bytes[at] ^= 1;
        std::fs::remove_file(&damaged).ok();
        std::fs::write(&damaged, &bytes).unwrap();
        assert!(open(&damaged).is_err(), "byte {at}");
    }
    let mut bytes = original.clone();
    let last = bytes.len() - 3;
    bytes[last] ^= 1;
    std::fs::remove_file(&damaged).ok();
    std::fs::write(&damaged, &bytes).unwrap();
    let packed = open(&damaged).unwrap();
    let mut out = Vec::new();
    assert!(packed.read_bases(0..10, &mut out).is_ok());
    assert!(packed.read_bases(0..packed.bases(), &mut out).is_err());
    std::fs::write(&damaged, &original[..original.len() - 1]).unwrap();
    assert!(open(&damaged).is_err());
}
