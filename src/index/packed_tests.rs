//! Round-trip and damage tests of the packed sequence format.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::storage::MappedFile;
use crate::test_support::rng;
use rand::RngExt;
use std::path::Path;

/// Contigs of random bases with runs of N, lower case bases and other bytes.
fn contigs(seed: u64) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut rng = rng(seed);
    [0usize, 1, 70_000, 5, 33_000]
        .iter()
        .enumerate()
        .map(|(index, &length)| {
            let mut sequence: Vec<u8> = (0..length)
                .map(|_| b"ACGTacgt"[rng.random_range(0..8)])
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

/// Writes the packed file of the contigs and returns their contig table.
fn write(path: &Path, contigs: &[(Vec<u8>, Vec<u8>)]) -> ContigTable {
    let mut writer = PackedWriter::new();
    for (name, sequence) in contigs {
        writer.add_contig(name, sequence).unwrap();
    }
    let table = writer.contig_table();
    std::fs::write(path, writer.finish().unwrap().parts().concat()).unwrap();
    table
}

fn open(path: &Path, table: &ContigTable) -> PackedSequence {
    PackedSequence::open(
        Arc::new(MappedFile::open(path, false).unwrap()),
        table.clone(),
    )
}

#[test]
fn packed_round_trip() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("entry.packed");
    let contigs = contigs(3);
    let table = write(&path, &contigs);
    let packed = open(&path, &table);
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
        }
        start += length;
    }
    let table = packed.contigs([1, 0]).unwrap();
    assert_eq!(table.len(), 2);
    assert_eq!(table[&1].length, contigs[1].1.len() as u64);
    assert_eq!(table[&1].name.as_bytes(), contigs[1].0);
    assert!(packed.contigs([5]).is_err());
    packed.will_need([(1, &(0..10)), (9, &(0..10))]);
    assert!(packed.read(2, 0..70_001, &mut out).is_err());
    assert!(packed.contig_name(5).is_err());
}

#[test]
fn empty_bad_names() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("empty.packed");
    let table = write(&path, &[]);
    let packed = open(&path, &table);
    assert_eq!((packed.contig_count(), packed.bases()), (0, 0));
    assert!(packed.contigs([]).unwrap().is_empty());
    let mut writer = PackedWriter::new();
    assert!(writer.add_contig(b"", b"ACGT").is_err());
    assert!(writer.add_contig(b"two words", b"ACGT").is_err());
}

/// Many contigs fill several name blocks; a read of some names reads only their blocks.
#[test]
fn names_in_blocks() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("many.packed");
    let contigs: Vec<(Vec<u8>, Vec<u8>)> = (0..600)
        .map(|index| (format!("c{index}").into_bytes(), vec![b'A'; index % 7 + 1]))
        .collect();
    let table = write(&path, &contigs);
    let packed = open(&path, &table);
    for index in [0, 255, 256, 511, 599] {
        assert_eq!(packed.contig_name(index).unwrap(), contigs[index].0);
    }
    let found = packed.contigs([599, 3, 300]).unwrap();
    assert_eq!(found.keys().copied().collect::<Vec<_>>(), [3, 300, 599]);
    assert_eq!(found[&300].length, 300 % 7 + 1);
    assert!(packed.contig_name(600).is_err());
}

/// An entry with more contigs than the header page indexes in blocks of 256 uses larger blocks.
#[test]
fn larger_name_blocks() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("larger.packed");
    let contigs: Vec<(Vec<u8>, Vec<u8>)> = (0..200_000)
        .map(|index| (format!("n{index}").into_bytes(), b"C".to_vec()))
        .collect();
    let table = write(&path, &contigs);
    assert_eq!(name_block(200_000), 512);
    let packed = open(&path, &table);
    for index in [0, 511, 512, 199_999] {
        assert_eq!(packed.contig_name(index).unwrap(), contigs[index].0);
    }
}

#[test]
fn damage_is_detected() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("entry.packed");
    let table = write(&path, &contigs(5));
    let original = std::fs::read(&path).unwrap();
    let damaged = directory.path().join("damaged.packed");
    let mut out = Vec::new();
    for at in [0, 30, 50] {
        let mut bytes = original.clone();
        bytes[at] ^= 1;
        std::fs::write(&damaged, &bytes).unwrap();
        let packed = open(&damaged, &table);
        assert!(packed.read_bases(0..10, &mut out).is_err(), "byte {at}");
        assert!(packed.contigs([0]).is_err(), "byte {at}");
    }
    // Damage in the names is found when their block is read.
    let mut bytes = original.clone();
    let last = bytes.len() - 6;
    bytes[last] ^= 1;
    std::fs::write(&damaged, &bytes).unwrap();
    let packed = open(&damaged, &table);
    assert!(packed.contig_name(0).is_err());
    assert!(packed.read_bases(0..10, &mut out).is_ok());
    // Damage in a block of bases is found when that block is read.
    let mut bytes = original.clone();
    bytes[BASES_AT as usize + BLOCK_BYTES as usize + 10] ^= 1;
    std::fs::write(&damaged, &bytes).unwrap();
    let packed = open(&damaged, &table);
    assert!(packed.read_bases(0..10, &mut out).is_ok());
    assert!(packed.read_bases(0..packed.bases(), &mut out).is_err());
    // A file of another size, or a table of another entry, fails at the first read.
    std::fs::write(&damaged, &original[..original.len() - 1]).unwrap();
    assert!(open(&damaged, &table).read_bases(0..10, &mut out).is_err());
    let other = ContigTable::new([1, 2], Vec::new());
    assert!(open(&path, &other).read_bases(0..1, &mut out).is_err());
}

/// Packing four bases at once gives the bytes and runs that packing base by base gives, for
/// contigs that start at every offset within a byte and hold other bytes anywhere.
#[test]
fn quads_match_bases() {
    let mut rng = rng(0x9a11);
    let (mut quads, mut bases) = (PackedWriter::new(), PackedWriter::new());
    for index in 0..400 {
        let length = rng.random_range(0..90);
        let sequence: Vec<u8> = (0..length)
            .map(|_| b"ACGTacgtACGTNnRy-"[rng.random_range(0..17)])
            .collect();
        let name = format!("contig_{index}").into_bytes();
        quads.add_contig(&name, &sequence).unwrap();
        bases.starts.push(bases.bases);
        sequence.iter().for_each(|&byte| bases.push_base(byte));
    }
    assert_eq!(quads.bases, bases.bases);
    assert_eq!(quads.packed, bases.packed);
    assert_eq!(quads.runs, bases.runs);
    assert!(quads.runs.len() > 100);
}
