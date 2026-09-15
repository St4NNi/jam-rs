//! Reshards an index: input shards into output shards of contiguous entry groups, balanced by their
//! seed locations.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

//! Entries keep their packed sequences; each output shard writes its seeds file from the seed
//! records of its entries, read from the seeds files of the inputs. The groups depend only on the
//! inputs and the split, never on the thread count.

use super::build::{complete_shard, in_partial_directory};
use super::index_file::{INDEX_FILE, INDEX_FORMAT_VERSION, IndexFile, ShardRecord};
use super::manifest::{
    MANIFEST_FILE, Manifest, PACKED_DIRECTORY, Parameters, check_names, sha256_file,
};
use super::seed_file::{ColorSetId, SEED_KEY_BITS, SEEDS_FILE, SeedFile};
use super::spill::{SeedSpill, seed_record};
use super::{AtPath, IndexError};
use crate::storage::LocalFile;
use rayon::prelude::*;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// How the entries are split into output shards.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Split {
    /// This many output shards, as even in seed locations as contiguous groups allow.
    Shards(usize),
    /// As few output shards as hold at most this many seed locations each; an entry with more
    /// seed locations gets a shard of its own.
    MaxSeedLocations(u64),
}

/// How the packed sequences reach the output shards, so each output shard opens on its own.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PackedOutput {
    /// A symbolic link to the original file of each packed sequence.
    Link,
    /// A copy of each packed sequence.
    Copy,
}

/// How an index is resharded.
#[derive(Clone, Debug)]
pub struct ReshardOptions {
    /// How the entries are split into output shards.
    pub split: Split,
    /// How the packed sequences reach the output shards.
    pub packed: PackedOutput,
    /// Threads that write each seeds file.
    pub threads: usize,
    /// Directory for temporary files.
    pub temporary_directory: PathBuf,
}

/// Reshards `inputs`, in order, into the new index directory `output`. An input is a shard
/// directory or an index directory (with an index file), whose shards are taken in order.
pub fn reshard(
    inputs: &[PathBuf],
    output: &Path,
    options: &ReshardOptions,
) -> Result<IndexFile, IndexError> {
    let shards = open_inputs(inputs)?;
    let mut sources = Vec::new();
    for (shard, input) in shards.iter().enumerate() {
        let entries = 0..input.manifest.entries.len() as u32;
        sources.extend(entries.map(|entry| EntrySource { shard, entry }));
    }
    let weights: Vec<u64> = (sources.iter())
        .map(|source| shards[source.shard].manifest.entries[source.entry as usize].seed_locations)
        .collect();
    let groups = split_groups(&weights, options.split)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.threads.max(1))
        .build()?;
    in_partial_directory(output, |directory| {
        let records = (groups.iter().enumerate())
            .map(|(number, group)| {
                let part = Part {
                    number,
                    first_entry: group.start,
                    sources: &sources[group.clone()],
                };
                pool.install(|| output_shard(&shards, &part, directory, options))
            })
            .collect::<Result<Vec<_>, IndexError>>()?;
        let index = IndexFile {
            format_version: INDEX_FORMAT_VERSION,
            shards: records,
        };
        index.write(directory)?;
        Ok(index)
    })
}

/// One input shard, opened for copying its entries.
struct InputShard {
    directory: PathBuf,
    manifest: Manifest,
    seeds: SeedFile,
}

/// The entries of one output shard, in index order.
struct Part<'a> {
    number: usize,
    first_entry: usize,
    sources: &'a [EntrySource],
}

/// An entry of an input shard.
#[derive(Clone, Copy)]
struct EntrySource {
    /// The input shard, by its place among the inputs.
    shard: usize,
    /// The entry in that shard.
    entry: u32,
}

/// Opens every input shard, expanding index directories, and checks that they fit together.
fn open_inputs(inputs: &[PathBuf]) -> Result<Vec<InputShard>, IndexError> {
    let mut directories = Vec::new();
    for input in inputs {
        if input.join(INDEX_FILE).exists() {
            let index = IndexFile::read(input)?;
            index.verify(input)?;
            directories.extend(index.shard_directories(input));
        } else {
            directories.push(input.clone());
        }
    }
    let shards = directories
        .into_iter()
        .map(|directory| {
            let manifest = Manifest::read(&directory)?;
            let file = LocalFile::open(directory.join(SEEDS_FILE))?;
            let seeds = SeedFile::open(Arc::new(file))?;
            if seeds.entry_count() != manifest.entries.len() as u64 || !seeds.has_locators() {
                return Err(IndexError::format(format!(
                    "the seeds file of {} differs from its manifest or has no locator part",
                    directory.display()
                )));
            }
            Ok(InputShard {
                directory,
                manifest,
                seeds,
            })
        })
        .collect::<Result<Vec<_>, IndexError>>()?;
    let first = shards
        .first()
        .ok_or_else(|| IndexError::Input("no shard to reshard".into()))?;
    let parameters = first.manifest.parameters;
    // Bucket bits follow each shard's size, and every output shard chooses its own.
    let same = |other: &Parameters| {
        Parameters {
            bucket_bits: parameters.bucket_bits,
            ..*other
        } == parameters
    };
    if shards.iter().any(|shard| !same(&shard.manifest.parameters)) {
        return Err(IndexError::Input(
            "shards with different parameters cannot be resharded together".into(),
        ));
    }
    check_names(shards.iter().flat_map(|shard| {
        shard
            .manifest
            .entries
            .iter()
            .map(|entry| entry.name.as_str())
    }))?;
    Ok(shards)
}

/// Writes output shard `part` into its directory under `directory` and returns its record.
fn output_shard(
    shards: &[InputShard],
    part: &Part<'_>,
    directory: &Path,
    options: &ReshardOptions,
) -> Result<ShardRecord, IndexError> {
    let name = format!("shard-{:04}", part.number);
    let shard_directory = directory.join(&name);
    std::fs::create_dir(&shard_directory).at(&shard_directory)?;
    let manifest = write_shard(shards, part.sources, &shard_directory, options)?;
    let manifest_digest = sha256_file(&shard_directory.join(MANIFEST_FILE))?;
    Ok(ShardRecord {
        directory: name,
        first_entry: part.first_entry as u64,
        entries: part.sources.len() as u64,
        seed_locations: manifest.seed_locations,
        manifest_sha256: manifest_digest.sha256,
    })
}

/// Writes one output shard of the entries `part` into `directory`.
fn write_shard(
    shards: &[InputShard],
    part: &[EntrySource],
    directory: &Path,
    options: &ReshardOptions,
) -> Result<Manifest, IndexError> {
    let packed = directory.join(PACKED_DIRECTORY);
    std::fs::create_dir(&packed).at(&packed)?;
    let mut entries = Vec::with_capacity(part.len());
    for &EntrySource { shard, entry } in part {
        let input = &shards[shard];
        let record = input.manifest.entries[entry as usize].clone();
        let source = input
            .directory
            .join(PACKED_DIRECTORY)
            .join(&record.packed.path);
        let target = packed.join(&record.packed.path);
        match options.packed {
            // Links point at the original file, so they stay valid when the inputs are removed
            // after an earlier reshard that linked them.
            PackedOutput::Link => {
                let original = std::fs::canonicalize(&source).at(&source)?;
                std::os::unix::fs::symlink(original, &target).at(&target)?
            }
            PackedOutput::Copy => drop(std::fs::copy(&source, &target).at(&target)?),
        }
        entries.push(record);
    }
    let spill = SeedSpill::create(&options.temporary_directory)?;
    for (shard, input) in shards.iter().enumerate() {
        // The new entry of each input entry that this output shard holds.
        let mut new_entries = vec![u32::MAX; input.manifest.entries.len()];
        for (new, source) in part.iter().enumerate() {
            if source.shard == shard {
                new_entries[source.entry as usize] = new as u32;
            }
        }
        if new_entries.iter().any(|&new| new != u32::MAX) {
            spill_input(&input.seeds, &new_entries, &spill)?;
        }
    }
    let parameters = shards[0].manifest.parameters;
    complete_shard(directory, parameters, entries, spill, None)
}

/// Buckets of an input that one parallel step reads.
const CHUNK_BUCKETS: u32 = 1 << 12;

/// Spills the seed records of the entries of an input seeds file that `new_entries` maps to an
/// entry of the output shard (others are `u32::MAX`).
fn spill_input(file: &SeedFile, new_entries: &[u32], spill: &SeedSpill) -> Result<(), IndexError> {
    let buckets = 1u32 << file.bucket_bits();
    let check_bits = SEED_KEY_BITS - file.bucket_bits();
    (0..buckets.div_ceil(CHUNK_BUCKETS))
        .into_par_iter()
        .try_for_each(|chunk| {
            let first = chunk * CHUNK_BUCKETS;
            let list: Vec<u32> = (first..buckets.min(first + CHUNK_BUCKETS)).collect();
            // Every stored seed of the chunk: its bucket, key and color set.
            let mut seeds: Vec<(u32, u64, ColorSetId)> = Vec::new();
            file.read_sketch_buckets(&list, |bucket, part| {
                for place in 0..part.len() {
                    let key = (u64::from(bucket) << check_bits) | part.check(place);
                    seeds.push((bucket, key, ColorSetId(part.color_set(place))));
                }
                Ok(())
            })?;
            let sets = file.read_color_sets(seeds.iter().map(|seed| seed.2))?;
            let mut records = Vec::new();
            file.read_locator_buckets(&list, |bucket, part| {
                let first = seeds.partition_point(|seed| seed.0 < bucket);
                part.for_each(|place, rank, location| {
                    let damaged = || IndexError::format("seeds locator bucket");
                    let (_, key, id) = *seeds.get(first + place as usize).ok_or_else(damaged)?;
                    let entries = sets.entries(id).ok_or_else(damaged)?;
                    let entry = *entries.get(rank as usize).ok_or_else(damaged)?;
                    match new_entries.get(entry as usize) {
                        Some(&new) if new != u32::MAX => {
                            records.push(seed_record(key, new, location));
                        }
                        _ => {}
                    }
                    Ok(())
                })
            })?;
            spill.push(&mut records)
        })
}

/// Cuts entries with these weights, in order, into contiguous non-empty groups by `split`.
pub(crate) fn split_groups(weights: &[u64], split: Split) -> Result<Vec<Range<usize>>, IndexError> {
    if weights.is_empty() {
        return Err(IndexError::Input("the inputs hold no entries".into()));
    }
    match split {
        Split::MaxSeedLocations(0) => Err(IndexError::Input(
            "the largest seed locations per shard must be at least 1".into(),
        )),
        Split::MaxSeedLocations(limit) => Ok(cut(weights, limit, usize::MAX)),
        Split::Shards(count) if count == 0 || count > weights.len() => {
            Err(IndexError::Input(format!(
                "the shard count must be 1 to {}, the entries of the inputs",
                weights.len()
            )))
        }
        Split::Shards(count) => {
            // The smallest limit whose greedy cut needs at most `count` groups.
            let (mut low, mut high) = (*weights.iter().max().unwrap_or(&0), weights.iter().sum());
            while low < high {
                let middle = low + (high - low) / 2;
                if cut(weights, middle, usize::MAX).len() <= count {
                    high = middle;
                } else {
                    low = middle + 1;
                }
            }
            Ok(cut(weights, low, count))
        }
    }
}

/// Fills each group up to `limit` weight. A group also ends when the entries left are just
/// enough to give each of the remaining groups one, so a count of groups is met exactly.
fn cut(weights: &[u64], limit: u64, groups: usize) -> Vec<Range<usize>> {
    let mut ranges = Vec::new();
    let (mut start, mut weight) = (0, 0u64);
    for (index, &next) in weights.iter().enumerate().skip(1) {
        weight += weights[index - 1];
        let groups_after = groups.saturating_sub(ranges.len() + 1);
        if weight + next > limit || weights.len() - index == groups_after {
            ranges.push(start..index);
            (start, weight) = (index, 0);
        }
    }
    ranges.push(start..weights.len());
    ranges
}

#[cfg(test)]
mod tests {
    use super::*;

    fn weights_of(weights: &[u64], groups: &[Range<usize>]) -> Vec<u64> {
        groups
            .iter()
            .map(|g| weights[g.clone()].iter().sum())
            .collect()
    }

    #[test]
    fn groups_stay_balanced() {
        let weights = [5, 1, 1, 1, 1, 1, 9, 2, 2, 2, 3, 1];
        for count in 1..=weights.len() {
            let groups = split_groups(&weights, Split::Shards(count)).unwrap();
            assert_eq!(groups.len(), count);
            assert_eq!(groups[0].start, 0);
            assert_eq!(groups.last().unwrap().end, weights.len());
            assert!(groups.windows(2).all(|pair| pair[0].end == pair[1].start));
            assert!(groups.iter().all(|group| !group.is_empty()));
        }
        let groups = split_groups(&weights, Split::Shards(3)).unwrap();
        assert_eq!(weights_of(&weights, &groups), [10, 9, 10]);
        let groups = split_groups(&[1, 1, 1, 1], Split::Shards(3)).unwrap();
        assert_eq!(groups, [0..2, 2..3, 3..4]);
        let groups = split_groups(&weights, Split::MaxSeedLocations(8)).unwrap();
        assert_eq!(weights_of(&weights, &groups), [8, 2, 9, 6, 4]);
        assert!(split_groups(&weights, Split::Shards(13)).is_err());
        assert!(split_groups(&weights, Split::Shards(0)).is_err());
        assert!(split_groups(&weights, Split::MaxSeedLocations(0)).is_err());
        assert!(split_groups(&[], Split::Shards(1)).is_err());
    }
}
