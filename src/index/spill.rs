//! Builds the seeds file of a shard from seed records spilled to temporary files by key range:
//! a first pass counts and ranks the color sets, a second one encodes the buckets.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::buckets::{encode_locator_bucket, encode_sketch_bucket};
use super::seed_file::{SEED_KEY_BITS, encode_color_set};
use super::seed_writer::{SeedFileShape, SeedWriter, WrittenSeeds};
use super::{AtPath, IndexError};
use crate::storage::lock;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Mutex, PoisonError};

/// The spill files split the key space into 2^PARTITION_BITS ranges.
const PARTITION_BITS: u32 = 12;
/// Lock stripes of the color set counts, so threads rarely wait for each other.
const COUNT_STRIPES: usize = 64;
/// Units encoded per parallel round, so encoded buckets wait in memory only briefly.
const WRITE_ROUND: usize = 64;

/// Keys per distinct color set (encoded), in lock stripes.
type SetCounts = Vec<FxHashMap<Box<[u8]>, u64>>;

/// One seed location of one entry: the seed key in bits 64 and up, then the entry, then the
/// seed location value. Sorting records sorts by key, entry and seed location.
pub(crate) fn seed_record(key: u64, entry: u32, location: u32) -> u128 {
    (u128::from(key) << 64) | (u128::from(entry) << 32) | u128::from(location)
}

/// Totals of a written seeds file.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct SeedTotals {
    /// Distinct seed keys.
    pub(crate) seeds: u64,
    pub(crate) color_sets: u64,
}

/// Seed records of a shard being built, in temporary files by key range.
pub(crate) struct SeedSpill {
    directory: tempfile::TempDir,
    paths: Vec<PathBuf>,
    writers: Vec<Mutex<BufWriter<File>>>,
}

impl SeedSpill {
    /// Creates the spill files in `temporary_directory`.
    pub(crate) fn create(temporary_directory: &Path) -> Result<Self, IndexError> {
        let directory = tempfile::Builder::new()
            .prefix("jam-spill")
            .tempdir_in(temporary_directory)
            .at(temporary_directory)?;
        let paths: Vec<PathBuf> = (0..1usize << PARTITION_BITS)
            .map(|partition| directory.path().join(format!("partition-{partition:04}")))
            .collect();
        let writers = (paths.iter())
            .map(|path| {
                let file = File::create_new(path).at(path)?;
                Ok(Mutex::new(BufWriter::with_capacity(1 << 16, file)))
            })
            .collect::<Result<_, IndexError>>()?;
        Ok(Self {
            directory,
            paths,
            writers,
        })
    }

    /// Adds seed records made by `seed_record`, from any thread; they are sorted here.
    pub(crate) fn push(&self, records: &mut [u128]) -> Result<(), IndexError> {
        records.sort_unstable();
        let partition = |record: &u128| (record >> (128 - 20 - PARTITION_BITS)) as usize;
        let mut bytes = Vec::new();
        for group in records.chunk_by(|a, b| partition(a) == partition(b)) {
            bytes.clear();
            bytes.extend(group.iter().flat_map(|record| record.to_le_bytes()));
            let index = partition(&group[0]);
            let path = &self.paths[index];
            lock(&self.writers[index]).write_all(&bytes).at(path)?;
        }
        Ok(())
    }

    /// Writes the seeds of the index file at `path` from the records, up to the manifest, and
    /// returns their totals.
    pub(crate) fn write(
        self,
        path: &Path,
        shape: SeedFileShape,
        contig_tables: &[Vec<u8>],
    ) -> Result<(SeedTotals, WrittenSeeds), IndexError> {
        for (writer, path) in self.writers.into_iter().zip(&self.paths) {
            let writer = writer.into_inner().unwrap_or_else(PoisonError::into_inner);
            writer
                .into_inner()
                .map_err(|error| error.into_error())
                .at(path)?;
        }
        let units = Units {
            paths: self.paths,
            bucket_bits: shape.bucket_bits,
        };
        let ranked = rank(count_color_sets(&units)?);
        let ids: FxHashMap<&[u8], u32> = (ranked.iter().enumerate())
            .map(|(id, set)| (&set[..], id as u32))
            .collect();
        let mut writer = SeedWriter::create(path, self.directory.path(), shape)?;
        let mut seeds = 0u64;
        let all: Vec<usize> = (0..units.count()).collect();
        for round in all.chunks(WRITE_ROUND) {
            let encoded = (round.par_iter())
                .map(|&unit| encode_unit(&units, unit, &ids, shape.locators))
                .collect::<Result<Vec<_>, IndexError>>()?;
            for unit in encoded {
                for bucket in unit {
                    writer.push_bucket(bucket.bucket, &bucket.sketch, &bucket.locators)?;
                    seeds += bucket.seeds;
                }
            }
        }
        let written = writer.finish(&ranked, contig_tables)?;
        let totals = SeedTotals {
            seeds,
            color_sets: ranked.len() as u64,
        };
        Ok((totals, written))
    }
}

/// The spill files grouped into units: a unit is one partition, or all partitions of one bucket
/// when buckets are fewer than partitions.
struct Units {
    paths: Vec<PathBuf>,
    bucket_bits: u32,
}

impl Units {
    fn count(&self) -> usize {
        1 << self.bucket_bits.min(PARTITION_BITS)
    }

    /// The records of one unit, sorted.
    fn read(&self, unit: usize) -> Result<Vec<u128>, IndexError> {
        let per_unit = self.paths.len() / self.count();
        let mut bytes = Vec::new();
        for path in &self.paths[unit * per_unit..(unit + 1) * per_unit] {
            File::open(path)
                .and_then(|mut file| file.read_to_end(&mut bytes))
                .at(path)?;
        }
        let mut records: Vec<u128> = (bytes.as_chunks::<16>().0.iter())
            .map(|chunk| u128::from_le_bytes(*chunk))
            .collect();
        records.sort_unstable();
        records.dedup();
        Ok(records)
    }
}

/// Calls `visit(key, entries, records)` for every key of sorted records: its distinct entries
/// and all its records.
fn for_each_key(
    records: &[u128],
    mut visit: impl FnMut(u64, &[u32], &[u128]) -> Result<(), IndexError>,
) -> Result<(), IndexError> {
    let mut entries = Vec::new();
    for group in records.chunk_by(|a, b| a >> 64 == b >> 64) {
        entries.clear();
        entries.extend(group.iter().map(|record| (record >> 32) as u32));
        entries.dedup();
        visit((group[0] >> 64) as u64, &entries, group)?;
    }
    Ok(())
}

/// Counts the keys of every distinct color set.
fn count_color_sets(units: &Units) -> Result<SetCounts, IndexError> {
    let stripes: Vec<Mutex<FxHashMap<Box<[u8]>, u64>>> =
        (0..COUNT_STRIPES).map(|_| Mutex::default()).collect();
    (0..units.count())
        .into_par_iter()
        .with_max_len(1)
        .try_for_each(|unit| {
            let records = units.read(unit)?;
            let mut local: FxHashMap<Box<[u8]>, u64> = FxHashMap::default();
            let mut encoded = Vec::new();
            for_each_key(&records, |_, entries, _| {
                encoded.clear();
                encode_color_set(entries, &mut encoded);
                match local.get_mut(&encoded[..]) {
                    Some(count) => *count += 1,
                    None => drop(local.insert(encoded[..].into(), 1)),
                }
                Ok(())
            })?;
            for (set, count) in local {
                let stripe = jamhash::jamhash_bytes(&set) as usize % COUNT_STRIPES;
                *lock(&stripes[stripe]).entry(set).or_default() += count;
            }
            Ok::<_, IndexError>(())
        })?;
    Ok(stripes
        .into_iter()
        .map(|stripe| stripe.into_inner().unwrap_or_else(PoisonError::into_inner))
        .collect())
}

/// The distinct color sets, most frequent first; ties are ordered by their encoding, so the
/// order never depends on the thread count.
fn rank(counts: SetCounts) -> Vec<Box<[u8]>> {
    let mut sets: Vec<(u64, Box<[u8]>)> = counts
        .into_iter()
        .flat_map(|stripe| stripe.into_iter().map(|(set, count)| (count, set)))
        .collect();
    sets.par_sort_unstable_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1)));
    sets.into_iter().map(|(_, set)| set).collect()
}

/// The encoded buckets of one bucket.
struct EncodedBucket {
    bucket: usize,
    seeds: u64,
    sketch: Vec<u8>,
    locators: Vec<u8>,
}

/// Encodes the buckets of one unit.
fn encode_unit(
    units: &Units,
    unit: usize,
    ids: &FxHashMap<&[u8], u32>,
    locators: bool,
) -> Result<Vec<EncodedBucket>, IndexError> {
    let records = units.read(unit)?;
    let check_bits = SEED_KEY_BITS - units.bucket_bits;
    let bucket_of = |record: &u128| (record >> (64 + check_bits)) as usize;
    let mut encoded = Vec::new();
    let mut set = Vec::new();
    for group in records.chunk_by(|a, b| bucket_of(a) == bucket_of(b)) {
        let mut seeds = Vec::new();
        let (mut seed_pairs, mut pair_locations, mut values) = (Vec::new(), Vec::new(), Vec::new());
        for_each_key(group, |key, entries, records| {
            set.clear();
            encode_color_set(entries, &mut set);
            let id = ids
                .get(&set[..])
                .ok_or_else(|| IndexError::format("a color set changed between two reads"))?;
            seeds.push((key & ((1 << check_bits) - 1), *id));
            seed_pairs.push(entries.len() as u32);
            for pair in records.chunk_by(|a, b| a >> 32 == b >> 32) {
                pair_locations.push(pair.len() as u32);
                values.extend(pair.iter().map(|record| *record as u32));
            }
            Ok(())
        })?;
        let mut bucket = EncodedBucket {
            bucket: bucket_of(&group[0]),
            seeds: seeds.len() as u64,
            sketch: Vec::new(),
            locators: Vec::new(),
        };
        encode_sketch_bucket(&seeds, check_bits, &mut bucket.sketch);
        if locators {
            encode_locator_bucket(&seed_pairs, &pair_locations, &values, &mut bucket.locators);
        }
        encoded.push(bucket);
    }
    Ok(encoded)
}
