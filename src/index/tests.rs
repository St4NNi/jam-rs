//! End-to-end tests: build shards from generated entries, find every seed, reshard.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::seeds::Seeder;
use crate::test_support::{
    build_options, fasta_text, random_bases, rng, runtime, test_shard, write_entries,
};
use rand::RngExt;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Entries of a few contigs each. They share stretches of one common sequence, so many seeds
/// belong to several entries, and hold runs of N, lower case bases and a low-complexity run.
fn generated_entries(count: usize) -> Vec<Vec<Vec<u8>>> {
    let mut rng = rng(0x1dea);
    let shared = random_bases(&mut rng, 30_000);
    (0..count)
        .map(|entry| {
            (0..1 + entry % 3)
                .map(|contig| {
                    let mut sequence = random_bases(&mut rng, 2_000 + 3_000 * contig);
                    let from = rng.random_range(0..20_000);
                    sequence.extend_from_slice(&shared[from..from + 5_000]);
                    sequence.extend(std::iter::repeat_n(b'A', 300));
                    sequence[100..140].fill(b'N');
                    sequence[500..600].make_ascii_lowercase();
                    sequence
                })
                .collect()
        })
        .collect()
}

fn options(directory: &Path, threads: usize) -> BuildOptions {
    let mut options = build_options(directory, threads);
    options.parameters.bin_bits = 6;
    options
}

/// Checks that every seed of every entry is found through the sketch and the positions, and
/// that the packed sequences give back the contigs.
fn check_every_seed(shard: &Shard, entries: &[Vec<Vec<u8>>]) {
    let parameters = shard.manifest().parameters;
    let seeder = Seeder::new(parameters.scheme).unwrap();
    let mut locations = Vec::new();
    let mut bases = Vec::new();
    for (entry, contigs) in entries.iter().enumerate() {
        let entry = entry as u32;
        let packed = shard.packed(entry).unwrap();
        let mut start = 0u64;
        for (index, sequence) in contigs.iter().enumerate() {
            packed
                .read(index, 0..sequence.len() as u64, &mut bases)
                .unwrap();
            assert_eq!(bases, sequence.to_ascii_uppercase());
            let mut seeds = 0;
            seeder.seeds(sequence, |seed| {
                seeds += 1;
                let id = shard.sketch().color_set(seed.hash).expect("seed in sketch");
                let set = shard.sketch().color_sets().entries(id).unwrap();
                assert!(set.collect::<Vec<_>>().contains(&entry));
                locations.clear();
                shard
                    .positions()
                    .locations(entry, seed.hash, &mut locations)
                    .unwrap();
                let expected = SeedLocation {
                    bin: ((start + seed.position) >> parameters.bin_bits) as u32,
                    reverse: seed.reverse,
                };
                assert!(locations.contains(&expected));
            });
            assert!(seeds > 100);
            start += sequence.len() as u64;
        }
        let lengths = contigs.iter().map(|sequence| sequence.len() as u64);
        let bounds: Vec<u64> = std::iter::once(0)
            .chain(lengths.scan(0, |end, length| {
                *end += length;
                Some(*end)
            }))
            .collect();
        assert_eq!(shard.positions().contig_starts(entry).unwrap(), bounds);
    }
}

#[test]
fn every_seed_found() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(7);
    let inputs = write_entries(directory.path(), &entries);
    let output = directory.path().join("shard");
    let manifest = build_shard(&inputs, &output, &options(directory.path(), 3)).unwrap();
    assert_eq!(manifest.entries[2].name, "entry2");
    assert_eq!(manifest.entries[2].contigs, 3);
    let shard = Shard::open(&output).unwrap();
    let packed = shard.packed(2).unwrap();
    assert_eq!(packed.contig_name(1).unwrap(), b"entry2_contig1");
    shard.verify().unwrap();
    check_every_seed(&shard, &entries);
    let stats = shard.stats();
    assert_eq!(stats.entries, 7);
    assert!(stats.color_sets > 7 && stats.seeds > 1000 && stats.seed_locations >= stats.seeds);
    // Seeds of shared stretches have color sets of several entries.
    let sets = shard.sketch().color_sets();
    let largest = (0..stats.color_sets as u32)
        .map(|id| sets.entries(ColorSetId(id)).unwrap().len())
        .max();
    assert!(largest > Some(2));
}

#[test]
fn threads_keep_output() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(5);
    let inputs = write_entries(directory.path(), &entries);
    let one = directory.path().join("one");
    let four = directory.path().join("four");
    build_shard(&inputs, &one, &options(directory.path(), 1)).unwrap();
    build_shard(&inputs, &four, &options(directory.path(), 4)).unwrap();
    for file in [
        manifest::MANIFEST_FILE,
        sketch::SKETCH_FILE,
        positions::POSITIONS_FILE,
        "packed/entry3.packed",
    ] {
        let read = |shard: &PathBuf| std::fs::read(shard.join(file)).unwrap();
        assert!(read(&one) == read(&four), "{file}");
    }
}

fn reshard_options(directory: &Path, split: Split, packed: PackedOutput) -> ReshardOptions {
    ReshardOptions {
        split,
        packed,
        threads: 3,
        temp_dir: directory.to_path_buf(),
    }
}

/// Checks that two shards hold the same manifest, sketch and positions.
fn assert_same_shard(one: &Path, other: &Path) {
    assert_eq!(Manifest::read(one).unwrap(), Manifest::read(other).unwrap());
    for file in [sketch::SKETCH_FILE, positions::POSITIONS_FILE] {
        let read = |shard: &Path| std::fs::read(shard.join(file)).unwrap();
        assert!(read(one) == read(other), "{file}");
    }
}

#[test]
fn merged_equals_whole() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(6);
    let inputs = write_entries(directory.path(), &entries);
    let options = options(directory.path(), 2);
    let whole = directory.path().join("whole");
    build_shard(&inputs, &whole, &options).unwrap();
    let parts: Vec<PathBuf> = [0..2, 2..3, 3..6]
        .into_iter()
        .enumerate()
        .map(|(index, range)| {
            let part = directory.path().join(format!("part{index}"));
            build_shard(&inputs[range], &part, &options).unwrap();
            part
        })
        .collect();
    for packed in [PackedOutput::Link, PackedOutput::Copy] {
        let merged = directory.path().join(format!("merged-{packed:?}"));
        let split = Split::Shards(1);
        let index = reshard(
            &parts,
            &merged,
            &reshard_options(directory.path(), split, packed),
        );
        let index = index.unwrap();
        assert_eq!(index, IndexFile::read(&merged).unwrap());
        assert_eq!(index.shards.len(), 1);
        assert_eq!(index.shards[0].entries, 6);
        let shard_directory = merged.join(&index.shards[0].directory);
        assert_same_shard(&whole, &shard_directory);
        let link = shard_directory.join("packed/entry0.packed");
        let is_link = std::fs::symlink_metadata(link).unwrap().is_symlink();
        assert_eq!(is_link, packed == PackedOutput::Link);
        let shard = Shard::open(&shard_directory).unwrap();
        shard.verify().unwrap();
        check_every_seed(&shard, &entries);
    }
    // A shard can read its packed sequences from another directory.
    let shard = Shard::open_with_packed(&parts[0], &whole.join("packed")).unwrap();
    check_every_seed(&shard, &entries[..2]);
}

#[test]
fn splits_equal_direct() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(7);
    let inputs = write_entries(directory.path(), &entries);
    let options = options(directory.path(), 2);
    let whole = directory.path().join("whole");
    build_shard(&inputs, &whole, &options).unwrap();
    let link = PackedOutput::Link;
    let merged = directory.path().join("merged");
    reshard(
        std::slice::from_ref(&whole),
        &merged,
        &reshard_options(directory.path(), Split::Shards(1), link),
    )
    .unwrap();
    let total = Manifest::read(&whole).unwrap().seed_locations;
    for split in [Split::Shards(3), Split::MaxSeedLocations(total / 4)] {
        let output = directory.path().join(format!("split-{split:?}"));
        // The merged index is given as an index directory, so its shards are read in order.
        let options_of = reshard_options(directory.path(), split, link);
        let index = reshard(std::slice::from_ref(&merged), &output, &options_of).unwrap();
        let mut names = Vec::new();
        for (number, record) in index.shards.iter().enumerate() {
            let range = record.first_entry as usize..(record.first_entry + record.entries) as usize;
            let direct = directory.path().join(format!("direct-{split:?}-{number}"));
            build_shard(&inputs[range.clone()], &direct, &options).unwrap();
            let shard_directory = output.join(&record.directory);
            assert_same_shard(&direct, &shard_directory);
            let shard = Shard::open(&shard_directory).unwrap();
            check_every_seed(&shard, &entries[range]);
            let manifest = shard.manifest();
            assert_eq!(manifest.seed_locations, record.seed_locations);
            names.extend(manifest.entries.iter().map(|entry| entry.name.clone()));
            if let Split::MaxSeedLocations(limit) = split {
                assert!(record.seed_locations <= limit || record.entries == 1);
            }
        }
        let expected: Vec<_> = (0..7).map(|index| format!("entry{index}")).collect();
        assert_eq!(names, expected);
        if split == Split::Shards(3) {
            assert_eq!(index.shards.len(), 3);
            let most = index.shards.iter().map(|shard| shard.seed_locations).max();
            assert!(most.unwrap() < total / 2, "groups are balanced");
        }
    }
    let error = reshard(
        &[whole.clone(), whole.clone()],
        &directory.path().join("twice"),
        &reshard_options(directory.path(), Split::Shards(1), link),
    );
    assert!(error.unwrap_err().to_string().contains("used twice"));
}

#[test]
fn forced_bucket_bits() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(3);
    let inputs = write_entries(directory.path(), &entries);
    let output = directory.path().join("shard");
    let options = BuildOptions {
        bucket_bits: Some(12),
        ..options(directory.path(), 2)
    };
    let manifest = build_shard(&inputs, &output, &options).unwrap();
    assert_eq!(manifest.parameters.bucket_bits, 12);
    let shard = Shard::open(&output).unwrap();
    assert_eq!(shard.sketch().bucket_bits(), 12);
    check_every_seed(&shard, &entries);
    let options = BuildOptions {
        bucket_bits: Some(29),
        ..options
    };
    assert!(build_shard(&inputs, &directory.path().join("bad"), &options).is_err());
}

#[test]
fn records_as_entries() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(4);
    let records: Vec<Vec<Vec<u8>>> = entries
        .iter()
        .flatten()
        .map(|contig| vec![contig.clone()])
        .collect();
    let path = directory.path().join("all.fa");
    let text = entries
        .iter()
        .enumerate()
        .map(|(entry, contigs)| fasta_text(entry, contigs));
    std::fs::write(&path, text.collect::<Vec<_>>().concat()).unwrap();
    let output = directory.path().join("shard");
    let options = BuildOptions {
        entry_per_record: true,
        ..options(directory.path(), 3)
    };
    let manifest = build_shard(&[EntryInput::from_path(path)], &output, &options).unwrap();
    assert_eq!(manifest.entries.len(), records.len());
    assert_eq!(manifest.entries[1].name, "entry1_contig0");
    check_every_seed(&Shard::open(&output).unwrap(), &records);
    let twice = directory.path().join("twice.fa");
    std::fs::write(&twice, fasta_text(0, &[records[0][0].clone()]).repeat(2)).unwrap();
    let error = build_shard(
        &[EntryInput::from_path(twice)],
        &output.with_file_name("x"),
        &options,
    );
    assert!(error.unwrap_err().to_string().contains("used twice"));
}

#[test]
fn bad_builds_cleaned() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(2);
    let mut inputs = write_entries(directory.path(), &entries);
    let output = directory.path().join("shard");
    inputs[1].name = inputs[0].name.clone();
    assert!(build_shard(&inputs, &output, &options(directory.path(), 2)).is_err());
    inputs[1].name = "other".into();
    inputs[1].path = directory.path().join("missing.fa");
    assert!(build_shard(&inputs, &output, &options(directory.path(), 2)).is_err());
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 2);
    assert_eq!(
        EntryInput::from_path("x/GCF_1.2.fna.gz".into()).name,
        "GCF_1.2"
    );
}

#[test]
fn packed_elsewhere_cached() {
    use crate::storage::{BlockCache, CacheConfig};
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(4);
    test_shard(directory.path(), &entries);
    let output = directory.path().join("shard");
    let elsewhere = directory.path().join("elsewhere");
    std::fs::rename(output.join("packed"), &elsewhere).unwrap();
    assert!(Shard::open(&output).unwrap().packed(0).is_err());
    let location = PackedLocation::Directory(elsewhere.clone());
    let options = PackedOptions {
        location,
        cache: None,
    };
    let shard = Shard::open_with(&output, &options).unwrap();
    shard.verify().unwrap();
    check_every_seed(&shard, &entries);
    assert_eq!(
        shard.packed_path(1).unwrap(),
        elsewhere.join("entry1.packed")
    );
    let cache = BlockCache::open(&CacheConfig::new(1 << 20, 4096), runtime()).unwrap();
    let options = PackedOptions {
        cache: Some(Arc::clone(&cache)),
        ..options
    };
    let shard = Shard::open_with(&output, &options).unwrap();
    check_every_seed(&shard, &entries);
    check_every_seed(&shard, &entries);
}

#[test]
fn packed_on_s3() {
    use crate::storage::{S3Config, S3Limits};
    let directory = tempfile::tempdir().unwrap();
    let inputs = write_entries(directory.path(), &generated_entries(1));
    let output = directory.path().join("shard");
    build_shard(&inputs, &output, &options(directory.path(), 1)).unwrap();
    let limits = S3Limits {
        tries: 1,
        ..S3Limits::default()
    };
    // Port 9 (discard) is closed on test machines, so the store cannot be reached.
    let endpoint = Some("http://127.0.0.1:9");
    let config = S3Config::new("local", endpoint, true, limits, runtime()).unwrap();
    let config = config.unsigned();
    let open = |prefix: &str| {
        let location = PackedLocation::S3 {
            prefix: prefix.into(),
            config: Box::new(config.clone()),
        };
        let options = PackedOptions {
            location,
            cache: None,
        };
        Shard::open_with(&output, &options)
    };
    assert!(open("s3:///packed/").is_err());
    assert!(open("/local/packed/").is_err());
    let shard = open("s3://bucket/packed").unwrap();
    assert!(shard.packed_path(0).is_err());
    assert!(shard.packed(0).is_err());
}
