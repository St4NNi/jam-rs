//! End-to-end tests: build shards from generated entries, find every seed, join shards.

use super::*;
use crate::seeds::Seeder;
use std::path::{Path, PathBuf};

fn random(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

fn random_bases(length: usize, state: &mut u64) -> Vec<u8> {
    (0..length)
        .map(|_| b"ACGT"[(random(state) % 4) as usize])
        .collect()
}

/// Entries of a few contigs each. They share stretches of one common sequence, so many seeds
/// belong to several entries, and hold runs of N, lower case bases and a low-complexity run.
fn generated_entries(count: usize) -> Vec<Vec<(String, Vec<u8>)>> {
    let mut state = 0x1dea_u64;
    let shared = random_bases(30_000, &mut state);
    (0..count)
        .map(|entry| {
            (0..1 + entry % 3)
                .map(|contig| {
                    let mut sequence = random_bases(2_000 + 3_000 * contig, &mut state);
                    let from = (random(&mut state) % 20_000) as usize;
                    sequence.extend_from_slice(&shared[from..from + 5_000]);
                    sequence.extend(std::iter::repeat_n(b'A', 300));
                    sequence[100..140].fill(b'N');
                    sequence[500..600].make_ascii_lowercase();
                    (format!("entry{entry}_contig{contig}"), sequence)
                })
                .collect()
        })
        .collect()
}

fn write_fasta(path: &Path, contigs: &[(String, Vec<u8>)]) {
    let mut text = Vec::new();
    for (name, sequence) in contigs {
        text.extend_from_slice(format!(">{name} some description\n").as_bytes());
        for line in sequence.chunks(70) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    std::fs::write(path, text).unwrap();
}

fn parameters() -> Parameters {
    Parameters {
        k: 21,
        window: 20,
        bin_bits: 6,
        ..Parameters::default()
    }
}

fn options(directory: &Path, threads: usize) -> BuildOptions {
    BuildOptions {
        parameters: parameters(),
        bucket_bits: None,
        threads,
        temp_dir: directory.to_path_buf(),
        entry_per_record: false,
    }
}

/// Writes one FASTA file per entry and returns the inputs.
fn inputs(directory: &Path, entries: &[Vec<(String, Vec<u8>)>]) -> Vec<EntryInput> {
    entries
        .iter()
        .enumerate()
        .map(|(index, contigs)| {
            let path = directory.join(format!("genome{index}.fa"));
            write_fasta(&path, contigs);
            EntryInput::from_path(path)
        })
        .collect()
}

/// Checks that every seed of every entry is found through the sketch and the positions, and
/// that the packed sequences give back the contigs.
fn check_every_seed(shard: &Shard, entries: &[Vec<(String, Vec<u8>)>]) {
    let parameters = shard.manifest().parameters;
    let seeder = Seeder::new(parameters.scheme()).unwrap();
    let mut locations = Vec::new();
    let mut bases = Vec::new();
    for (entry, contigs) in entries.iter().enumerate() {
        let entry = entry as u32;
        let packed = shard.packed(entry).unwrap();
        let mut start = 0u64;
        for (index, (name, sequence)) in contigs.iter().enumerate() {
            assert_eq!(packed.contig_name(index).unwrap(), name.as_bytes());
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
    }
}

#[test]
fn every_seed_is_found_through_the_sketch_and_the_positions() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(7);
    let inputs = inputs(directory.path(), &entries);
    let output = directory.path().join("shard");
    let manifest = build_shard(&inputs, &output, &options(directory.path(), 3)).unwrap();
    assert_eq!(manifest.entries[2].name, "genome2");
    assert_eq!(manifest.entries[2].contigs, 3);
    let shard = Shard::open(&output).unwrap();
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
fn output_does_not_depend_on_the_thread_count() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(5);
    let inputs = inputs(directory.path(), &entries);
    let one = directory.path().join("one");
    let four = directory.path().join("four");
    build_shard(&inputs, &one, &options(directory.path(), 1)).unwrap();
    build_shard(&inputs, &four, &options(directory.path(), 4)).unwrap();
    for file in [
        manifest::MANIFEST_FILE,
        sketch::SKETCH_FILE,
        positions::POSITIONS_FILE,
        "packed/genome3.packed",
    ] {
        let read = |shard: &PathBuf| std::fs::read(shard.join(file)).unwrap();
        assert!(read(&one) == read(&four), "{file}");
    }
}

#[test]
fn joined_shards_equal_one_shard_built_at_once() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(6);
    let inputs = inputs(directory.path(), &entries);
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
    let joined = directory.path().join("joined");
    let join_options = JoinOptions {
        threads: 3,
        temp_dir: directory.path().to_path_buf(),
        bucket_bits: None,
        link_packed: true,
    };
    let manifest = join_shards(&parts, &joined, &join_options).unwrap();
    let chosen = manifest::choose_bucket_bits(manifest.seed_locations);
    assert_eq!(manifest.parameters.bucket_bits, chosen);
    assert_eq!(manifest, Manifest::read(&whole).unwrap());
    for file in [sketch::SKETCH_FILE, positions::POSITIONS_FILE] {
        let read = |shard: &PathBuf| std::fs::read(shard.join(file)).unwrap();
        assert!(read(&whole) == read(&joined), "{file}");
    }
    let shard = Shard::open(&joined).unwrap();
    shard.verify().unwrap();
    check_every_seed(&shard, &entries);
    // Without links, the joined shard reads the packed sequences from a given directory.
    let unlinked = directory.path().join("unlinked");
    let join_options = JoinOptions {
        link_packed: false,
        ..join_options
    };
    join_shards(&parts, &unlinked, &join_options).unwrap();
    assert!(Shard::open(&unlinked).unwrap().packed(0).is_err());
    let shard = Shard::open_with_packed(&unlinked, &whole.join("packed")).unwrap();
    check_every_seed(&shard, &entries);
}

#[test]
fn bucket_bits_can_be_forced() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(3);
    let inputs = inputs(directory.path(), &entries);
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
fn each_record_can_be_an_entry() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(4);
    let records: Vec<Vec<(String, Vec<u8>)>> = entries
        .iter()
        .flatten()
        .map(|contig| vec![contig.clone()])
        .collect();
    let path = directory.path().join("all.fa");
    write_fasta(&path, &records.concat());
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
    write_fasta(&twice, &[records[0][0].clone(), records[0][0].clone()]);
    let error = build_shard(
        &[EntryInput::from_path(twice)],
        &output.with_file_name("x"),
        &options,
    );
    assert!(error.unwrap_err().to_string().contains("used twice"));
}

#[test]
fn bad_builds_leave_nothing_behind() {
    let directory = tempfile::tempdir().unwrap();
    let entries = generated_entries(2);
    let mut inputs = inputs(directory.path(), &entries);
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
