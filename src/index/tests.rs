//! End-to-end tests: build shards from generated entries, find every seed, reshard.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::plan::{Part, Place, Plan, part_pins};
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

/// Checks that every seed of every entry is found in the seeds file at its seed location, and
/// that the packed sequences give back the contigs.
fn check_every_seed(shard: &Shard, entries: &[Vec<Vec<u8>>]) {
    let parameters = shard.manifest().parameters;
    let seeder = Seeder::new(parameters.scheme).unwrap();
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
                let expected = SeedLocation {
                    bin: ((start + seed.position) >> parameters.bin_bits) as u32,
                    reverse: seed.reverse,
                };
                check_seed_file(shard.seed_file(), seed.hash, entry, expected);
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
        assert_eq!(shard.contig_table(entry).unwrap().starts(), bounds);
    }
}

/// Checks that the seeds file holds a seed in an entry at a seed location.
fn check_seed_file(file: &SeedFile, hash: u64, entry: u32, expected: SeedLocation) {
    let stored = file.find_seeds(&[seed_key(hash)]).unwrap()[0].expect("seed in the seeds file");
    let sets = file.read_color_sets([stored.color_set]).unwrap();
    let set = sets.entries(stored.color_set).unwrap();
    let rank = set.binary_search(&entry).expect("entry in the color set") as u32;
    let mut values = Vec::new();
    file.read_locator_buckets(&[stored.bucket], |_, bucket| {
        bucket.locations(stored.place, rank, &mut values)
    })
    .unwrap();
    assert!(values.contains(&expected.value()));
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
    let ids = (0..stats.color_sets as u32).map(ColorSetId);
    let sets = shard.seed_file().read_color_sets(ids.clone()).unwrap();
    let largest = ids.map(|id| sets.entries(id).unwrap().len()).max();
    assert!(largest > Some(2));
}

#[test]
fn threads_keep_output() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(5);
    let inputs = write_entries(directory.path(), &entries);
    let [one, four] = ["one", "four"].map(|name| directory.path().join(name));
    for (output, threads) in [(&one, 1), (&four, 4)] {
        std::fs::create_dir(output).unwrap();
        let options = options(directory.path(), threads);
        build_shard(&inputs, &output.join("shard"), &options).unwrap();
    }
    for file in ["shard/seeds", "shard/packed/entry3.packed"] {
        let read = |output: &PathBuf| std::fs::read(output.join(file)).unwrap();
        assert!(read(&one) == read(&four), "{file}");
    }
}

fn reshard_options(directory: &Path, split: Split) -> ReshardOptions {
    ReshardOptions {
        split,
        threads: 3,
        temporary_directory: directory.to_path_buf(),
    }
}

/// Checks that two shards hold the same seeds file, with the same manifest.
fn assert_same_shard(one: &Path, other: &Path) {
    let read = |shard: &Path| std::fs::read(shard.join(SEEDS_FILE)).unwrap();
    assert!(read(one) == read(other));
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
    let merged = directory.path().join("merged");
    let options = reshard_options(directory.path(), Split::Shards(1));
    let index = reshard(&parts, &merged, &options).unwrap();
    assert_eq!(index, IndexFile::read(&merged).unwrap());
    assert_eq!(index.shards.len(), 1);
    assert_eq!(index.shards[0].entries, 6);
    let shard_directory = merged.join(&index.shards[0].directory);
    assert_same_shard(&whole, &shard_directory);
    let shard = Shard::open(&shard_directory).unwrap();
    shard.verify().unwrap();
    check_every_seed(&shard, &entries);
}

#[test]
fn splits_equal_direct() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(7);
    let inputs = write_entries(directory.path(), &entries);
    let options = options(directory.path(), 2);
    let whole = directory.path().join("whole");
    build_shard(&inputs, &whole, &options).unwrap();
    let merged = directory.path().join("merged");
    reshard(
        std::slice::from_ref(&whole),
        &merged,
        &reshard_options(directory.path(), Split::Shards(1)),
    )
    .unwrap();
    let total = Shard::open(&whole).unwrap().manifest().seed_locations;
    for split in [Split::Shards(3), Split::MaxSeedLocations(total / 4)] {
        let output = directory.path().join(format!("split-{split:?}"));
        // The merged index is given as an index directory, so its shards are read in order.
        let options_of = reshard_options(directory.path(), split);
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
        &reshard_options(directory.path(), Split::Shards(1)),
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
    assert_eq!(shard.seed_file().bucket_bits(), 12);
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
        singleton: true,
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
fn packed_read_cached() {
    use crate::storage::{BlockCache, CacheConfig};
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(4);
    test_shard(directory.path(), &entries);
    let local = Location::Directory(directory.path().join("shard"));
    let cache = BlockCache::open(&CacheConfig::new(1 << 20, 4096), runtime()).unwrap();
    let options = ShardOptions {
        cache: Some(Arc::clone(&cache)),
        ..ShardOptions::default()
    };
    let shard = Shard::open_with(&local, &options).unwrap();
    shard.verify().unwrap();
    check_every_seed(&shard, &entries);
    check_every_seed(&shard, &entries);
}

#[test]
fn shard_on_s3() {
    use crate::storage::{BlockCache, CacheConfig, S3Config, S3Limits};
    use crate::test_support::serve_directory;
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(3);
    test_shard(directory.path(), &entries);
    let (server, _) = serve_directory(directory.path());
    let endpoint = Some(server.address.as_str());
    let limits = S3Limits::default();
    let config = S3Config::new("local", endpoint, true, limits, runtime()).unwrap();
    let config = config.unsigned();
    assert!(Location::s3("s3:///shard", &config).is_err());
    assert!(Location::s3("/local/shard", &config).is_err());
    let location = Location::s3("s3://bucket/shard", &config).unwrap();
    let cache = BlockCache::open(&CacheConfig::new(1 << 20, 4096), runtime()).unwrap();
    let shard = Shard::open_with(&location, &ShardOptions::default()).unwrap();
    assert!(shard.packed_path(0).is_err());
    shard.verify().unwrap();
    let local = Shard::open(&directory.path().join("shard")).unwrap();
    assert_eq!(shard.digest(), local.digest());
    // Each seed lookup is a small read, so the check reads through the cache.
    let options = ShardOptions {
        cache: Some(cache),
        ..ShardOptions::default()
    };
    let shard = Shard::open_with(&location, &options).unwrap();
    check_every_seed(&shard, &entries);
    let missing = Location::s3("s3://bucket/missing", &config).unwrap();
    assert!(Shard::open_with(&missing, &ShardOptions::default()).is_err());
}

#[test]
fn pins_read_locally() {
    use crate::storage::{S3Config, S3Limits};
    use crate::test_support::serve_directory;
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(3);
    let local = test_shard(directory.path(), &entries);
    let (server, asked) = serve_directory(directory.path());
    let endpoint = Some(server.address.as_str());
    let config = S3Config::new("local", endpoint, true, S3Limits::default(), runtime());
    let location = Location::s3("s3://bucket/shard", &config.unwrap().unsigned()).unwrap();
    let head = local.head();
    let mut pins = part_pins(head, &Part::Sketch, Place::Ram).unwrap();
    pins.extend(part_pins(head, &Part::Locators, Place::Ssd).unwrap());
    let named = Part::Packed(Some(vec!["entry1".into()]));
    pins.extend(part_pins(head, &named, Place::Ram).unwrap());
    let plan = Plan {
        shard: local.digest().to_string(),
        memory: 0,
        disk: 0,
        pins,
    };
    let cache = directory.path().join("cache");
    plan.write(&cache).unwrap();
    for pin in &plan.pins {
        let path = plan.pin_path(&cache, pin);
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        let bytes = std::fs::read(directory.path().join("shard").join(&pin.file)).unwrap();
        std::fs::write(&path, &bytes[pin.start as usize..pin.end as usize]).unwrap();
        let etag = crate::test_support::etag_of(&bytes);
        std::fs::write(super::plan::etag_file(&path), etag).unwrap();
    }
    let options = ShardOptions {
        pins: Some(cache),
        ..ShardOptions::default()
    };
    // The head read to find the plan also opens the shard: two requests, which also check the
    // seeds file for this process, so its pins need no check read.
    let head = location.head().unwrap();
    let shard = Shard::open_with_head(&location, &options, head).unwrap();
    check_every_seed(&shard, &entries);
    let asked = asked.lock().unwrap().clone();
    assert_eq!(asked.iter().filter(|key| *key == "shard/seeds").count(), 2);
    assert!(asked.contains(&"shard/packed/entry0.packed".to_string()));
    // A pinned file gets nothing but the check read of its first use.
    let pinned = asked
        .iter()
        .filter(|key| *key == "shard/packed/entry1.packed");
    assert!(pinned.count() <= 1);
}

/// A shard with more buckets than partitions holds several bucket blocks per partition, and its
/// many color sets give rare ids that blocks store as escapes. Every seed must still find its
/// entry.
#[test]
fn escaped_sets_found() {
    let mut rng = rng(0xe5c);
    let shared = random_bases(&mut rng, 40_000);
    let entries: Vec<Vec<Vec<u8>>> = (0..60)
        .map(|_| {
            let mut sequence = random_bases(&mut rng, 40_000);
            for _ in 0..4 {
                let from = rng.random_range(0..38_000);
                sequence.extend_from_slice(&shared[from..from + 2_000]);
            }
            vec![sequence]
        })
        .collect();
    let directory = tempfile::tempdir().unwrap();
    let inputs = write_entries(directory.path(), &entries);
    let output = directory.path().join("shard");
    let options = BuildOptions {
        bucket_bits: Some(13),
        ..options(directory.path(), 2)
    };
    build_shard(&inputs, &output, &options).unwrap();
    let shard = Shard::open(&output).unwrap();
    check_every_seed(&shard, &entries);
}

#[test]
fn sketch_only_shards() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(5);
    let inputs = write_entries(directory.path(), &entries);
    let sketch_only = |name: &str, range: std::ops::Range<usize>| {
        let output = directory.path().join(name);
        let options = BuildOptions {
            sketch_only: true,
            ..options(directory.path(), 2)
        };
        build_shard(&inputs[range], &output, &options).unwrap();
        output
    };
    let whole = sketch_only("whole", 0..5);
    assert!(!whole.join("packed").exists());
    let shard = Shard::open(&whole).unwrap();
    shard.verify().unwrap();
    assert!(!shard.seed_file().has_locators());
    assert_eq!(shard.stats().locator_bytes, 0);
    assert_eq!(shard.stats().packed_bytes, 0);
    assert!(shard.packed(0).is_err() && shard.contig_table(0).is_err());
    let full = directory.path().join("full");
    build_shard(&inputs, &full, &options(directory.path(), 2)).unwrap();
    let full_sketch = Shard::open(&full).unwrap().manifest().sketch_part.clone();
    let full_bytes = std::fs::read(full.join(SEEDS_FILE)).unwrap();
    let whole_bytes = std::fs::read(whole.join(SEEDS_FILE)).unwrap();
    // The sketch part is the same as in a full shard; only the header and table differ.
    let sketch = shard.manifest().sketch_part.clone();
    assert_eq!(
        whole_bytes[sketch.start as usize..sketch.end as usize],
        full_bytes[full_sketch.start as usize..full_sketch.end as usize]
    );
    let output = directory.path().join("split");
    let options = reshard_options(directory.path(), Split::Shards(2));
    let index = reshard(std::slice::from_ref(&whole), &output, &options).unwrap();
    for record in &index.shards {
        let range = record.first_entry as usize..(record.first_entry + record.entries) as usize;
        let direct = sketch_only(&format!("direct-{}", record.directory), range);
        assert_same_shard(&direct, &output.join(&record.directory));
    }
    let mixed = [whole, full];
    assert!(reshard(&mixed, &directory.path().join("mixed"), &options).is_err());
}
