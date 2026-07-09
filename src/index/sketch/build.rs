//! Builds the sketch of a shard from its positions file, so building and joining shards share
//! one path. The seed keys of every entry are spilled into partitions of the key space; each
//! partition is then sorted, which gives the color set of every key. Color sets are counted and
//! ranked by frequency, and a second pass over the partitions writes the bucket blocks.

use super::{SketchWriter, encode_bucket, encode_color_set};
use crate::index::IndexError;
use crate::index::positions::{POSITION_KEY_BITS, Positions, decode_section};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

/// Partitions of the key space are 2^PARTITION_BITS (fewer when there are fewer buckets).
const PARTITION_BITS: u32 = 12;
/// Shards of the color set counts, so threads rarely wait for each other.
const COUNT_SHARDS: usize = 64;
/// Partitions written per parallel round, so encoded blocks wait in memory only briefly.
const WRITE_ROUND: usize = 64;

/// Keys per distinct color set (encoded), in shards.
type SetCounts = Vec<HashMap<Box<[u8]>, u64>>;

/// The bucket blocks of one partition: the block bytes, and for each non-empty bucket its index,
/// block length and seed count.
type EncodedPartition = (Vec<u8>, Vec<(usize, usize, u64)>);

/// Totals of a written sketch.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct SketchTotals {
    pub(crate) seeds: u64,
    pub(crate) color_sets: u64,
}

/// The spill files: one per partition of the key space. A record is a u64: the key bits below
/// the partition bits in the high half, the entry in the low half.
struct Spill {
    paths: Vec<PathBuf>,
    key_bits: u32,
    partition_bits: u32,
}

impl Spill {
    fn low_bits(&self) -> u32 {
        self.key_bits - self.partition_bits
    }

    fn read(&self, partition: usize) -> Result<Vec<u64>, IndexError> {
        let mut bytes = Vec::new();
        File::open(&self.paths[partition])?.read_to_end(&mut bytes)?;
        let mut records: Vec<u64> = bytes
            .as_chunks::<8>()
            .0
            .iter()
            .map(|chunk| u64::from_le_bytes(*chunk))
            .collect();
        records.sort_unstable();
        Ok(records)
    }
}

/// Writes the sketch at `path` from `positions`, with temporary files in `temp_dir`.
pub(crate) fn build_sketch(
    positions: &Positions,
    bucket_bits: u32,
    path: &Path,
    temp_dir: &Path,
) -> Result<SketchTotals, IndexError> {
    let directory = tempfile::Builder::new()
        .prefix("jam-sketch")
        .tempdir_in(temp_dir)?;
    let partition_bits = PARTITION_BITS.min(bucket_bits);
    let spill = Spill {
        paths: (0..1usize << partition_bits)
            .map(|partition| directory.path().join(format!("partition-{partition:04}")))
            .collect(),
        key_bits: bucket_bits + u32::from(crate::index::FINGERPRINT_BITS),
        partition_bits,
    };
    spill_keys(positions, &spill)?;
    let counts = count_color_sets(&spill)?;
    let ranked = rank(counts);
    let ids: HashMap<&[u8], u32> = ranked
        .iter()
        .enumerate()
        .map(|(id, set)| (&set[..], id as u32))
        .collect();
    let mut writer = SketchWriter::create(path, bucket_bits, u64::from(positions.entries()))?;
    let mut seeds = 0u64;
    let partitions: Vec<usize> = (0..spill.paths.len()).collect();
    for round in partitions.chunks(WRITE_ROUND) {
        let encoded = round
            .par_iter()
            .map(|&partition| encode_partition(&spill, partition, &ids))
            .collect::<Result<Vec<_>, IndexError>>()?;
        for (bytes, blocks) in encoded {
            let mut at = 0;
            for (bucket, length, count) in blocks {
                writer.push_bucket(bucket, &bytes[at..at + length])?;
                at += length;
                seeds += count;
            }
        }
    }
    writer.finish(&ranked)?;
    Ok(SketchTotals {
        seeds,
        color_sets: ranked.len() as u64,
    })
}

/// Writes the distinct keys of every entry into the spill partitions.
fn spill_keys(positions: &Positions, spill: &Spill) -> Result<(), IndexError> {
    let writers = spill
        .paths
        .iter()
        .map(|path| {
            Ok(Mutex::new(BufWriter::with_capacity(
                1 << 16,
                File::create_new(path)?,
            )))
        })
        .collect::<Result<Vec<_>, IndexError>>()?;
    let low_bits = spill.low_bits();
    let low_mask = (1u64 << low_bits) - 1;
    (0..positions.entries())
        .into_par_iter()
        .with_max_len(1)
        .try_for_each(|entry| {
            let bytes = positions.section_bytes(entry)?;
            let mut keys = Vec::new();
            let shift = POSITION_KEY_BITS - spill.key_bits;
            decode_section(positions.section(entry)?, &bytes, |key, _| {
                let key = key >> shift;
                if keys.last() != Some(&key) {
                    keys.push(key);
                }
            })?;
            drop(bytes);
            let mut records = Vec::with_capacity(8 * 4096);
            for group in keys.chunk_by(|a, b| a >> low_bits == b >> low_bits) {
                records.clear();
                for key in group {
                    let record = ((key & low_mask) << 32) | u64::from(entry);
                    records.extend_from_slice(&record.to_le_bytes());
                }
                let partition = (group[0] >> low_bits) as usize;
                let mut writer = writers[partition].lock().unwrap_or_else(|e| e.into_inner());
                writer.write_all(&records)?;
            }
            Ok::<_, IndexError>(())
        })?;
    for writer in writers {
        let writer = writer.into_inner().unwrap_or_else(|e| e.into_inner());
        writer.into_inner().map_err(|error| error.into_error())?;
    }
    Ok(())
}

/// Calls `visit(key low bits, encoded color set)` for every key of a sorted partition.
fn for_each_key(records: &[u64], mut visit: impl FnMut(u64, &[u8])) {
    let mut entries = Vec::new();
    let mut encoded = Vec::new();
    for group in records.chunk_by(|a, b| a >> 32 == b >> 32) {
        entries.clear();
        entries.extend(group.iter().map(|record| *record as u32));
        encoded.clear();
        encode_color_set(&entries, &mut encoded);
        visit(group[0] >> 32, &encoded);
    }
}

/// Counts the keys of every distinct color set.
fn count_color_sets(spill: &Spill) -> Result<SetCounts, IndexError> {
    let shards: Vec<Mutex<HashMap<Box<[u8]>, u64>>> =
        (0..COUNT_SHARDS).map(|_| Mutex::default()).collect();
    (0..spill.paths.len())
        .into_par_iter()
        .with_max_len(1)
        .try_for_each(|partition| {
            let records = spill.read(partition)?;
            let mut local: HashMap<Box<[u8]>, u64> = HashMap::new();
            for_each_key(&records, |_, set| match local.get_mut(set) {
                Some(count) => *count += 1,
                None => {
                    local.insert(set.into(), 1);
                }
            });
            for (set, count) in local {
                let shard = jamhash::jamhash_bytes(&set) as usize % COUNT_SHARDS;
                let mut shard = shards[shard].lock().unwrap_or_else(|e| e.into_inner());
                *shard.entry(set).or_default() += count;
            }
            Ok::<_, IndexError>(())
        })?;
    Ok(shards
        .into_iter()
        .map(|shard| shard.into_inner().unwrap_or_else(|e| e.into_inner()))
        .collect())
}

/// The distinct color sets, most frequent first; ties are ordered by their encoding, so the
/// order never depends on the thread count.
fn rank(counts: SetCounts) -> Vec<Box<[u8]>> {
    let mut sets: Vec<(u64, Box<[u8]>)> = counts
        .into_iter()
        .flat_map(|shard| shard.into_iter().map(|(set, count)| (count, set)))
        .collect();
    sets.par_sort_unstable_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1)));
    sets.into_iter().map(|(_, set)| set).collect()
}

/// Encodes the bucket blocks of one partition.
fn encode_partition(
    spill: &Spill,
    partition: usize,
    ids: &HashMap<&[u8], u32>,
) -> Result<EncodedPartition, IndexError> {
    let records = spill.read(partition)?;
    let mut seeds: Vec<(u64, u32)> = Vec::new();
    let mut missing = false;
    for_each_key(&records, |low, set| match ids.get(set) {
        Some(id) => seeds.push((low, *id)),
        None => missing = true,
    });
    if missing {
        return Err(IndexError::format("a color set changed between two reads"));
    }
    let first_key = (partition as u64) << spill.low_bits();
    let fingerprint_bits = u32::from(crate::index::FINGERPRINT_BITS);
    let mut bytes = Vec::new();
    let mut blocks = Vec::new();
    let mut bucket_seeds = Vec::new();
    for group in seeds.chunk_by(|a, b| a.0 >> fingerprint_bits == b.0 >> fingerprint_bits) {
        let start = bytes.len();
        bucket_seeds.clear();
        bucket_seeds.extend(group.iter().map(|(low, id)| (*low as u16, *id)));
        encode_bucket(&bucket_seeds, &mut bytes);
        let bucket = ((first_key + group[0].0) >> fingerprint_bits) as usize;
        blocks.push((bucket, bytes.len() - start, group.len() as u64));
    }
    Ok((bytes, blocks))
}
