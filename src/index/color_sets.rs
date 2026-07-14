//! Builds the sketch of a shard from its positions file, so building and joining shards share
//! one path. The seed keys of every entry are spilled into partitions of the key space; each
//! partition is then sorted, which gives the color set of every key. Color sets are counted and
//! ranked by frequency, and a second pass over the partitions writes the bucket blocks.

use super::sketch::{SketchWriter, encode_bucket, encode_color_set};
use crate::index::positions::{POSITION_KEY_BITS, Positions, decode_section};
use crate::index::{AtPath, IndexError};
use crate::storage::lock;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Mutex, PoisonError};

/// Partitions of the key space are 2^PARTITION_BITS (fewer when there are fewer buckets).
const PARTITION_BITS: u32 = 12;
/// Lock stripes of the color set counts, so threads rarely wait for each other.
const COUNT_STRIPES: usize = 64;
/// Partitions written per parallel round, so encoded blocks wait in memory only briefly.
const WRITE_ROUND: usize = 64;

/// Keys per distinct color set (encoded), in lock stripes.
type SetCounts = Vec<HashMap<Box<[u8]>, u64>>;

/// The bucket blocks of one partition, one after the other.
struct EncodedPartition {
    bytes: Vec<u8>,
    blocks: Vec<EncodedBucket>,
}

/// One non-empty bucket of an encoded partition.
struct EncodedBucket {
    bucket: usize,
    length: usize,
    seeds: u64,
}

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
        let path = &self.paths[partition];
        File::open(path)
            .and_then(|mut file| file.read_to_end(&mut bytes))
            .at(path)?;
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

/// Writes the sketch at `path` from `positions`, with temporary files in `temporary_directory`.
pub(crate) fn build_sketch(
    positions: &Positions,
    bucket_bits: u32,
    path: &Path,
    temporary_directory: &Path,
) -> Result<SketchTotals, IndexError> {
    let directory = tempfile::Builder::new()
        .prefix("jam-sketch")
        .tempdir_in(temporary_directory)
        .at(temporary_directory)?;
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
    let mut writer = SketchWriter::create(path, bucket_bits, u64::from(positions.entry_count()))?;
    let mut seeds = 0u64;
    let partitions: Vec<usize> = (0..spill.paths.len()).collect();
    for round in partitions.chunks(WRITE_ROUND) {
        let encoded = round
            .par_iter()
            .map(|&partition| encode_partition(&spill, partition, &ids))
            .collect::<Result<Vec<_>, IndexError>>()?;
        for partition in encoded {
            let mut at = 0;
            for block in partition.blocks {
                let end = at + block.length;
                writer.push_bucket(block.bucket, &partition.bytes[at..end])?;
                at = end;
                seeds += block.seeds;
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
                File::create_new(path).at(path)?,
            )))
        })
        .collect::<Result<Vec<_>, IndexError>>()?;
    let low_bits = spill.low_bits();
    let low_mask = (1u64 << low_bits) - 1;
    (0..positions.entry_count())
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
                let mut writer = lock(&writers[partition]);
                writer.write_all(&records).at(&spill.paths[partition])?;
            }
            Ok::<_, IndexError>(())
        })?;
    for (writer, path) in writers.into_iter().zip(&spill.paths) {
        let writer = writer.into_inner().unwrap_or_else(PoisonError::into_inner);
        writer
            .into_inner()
            .map_err(|error| error.into_error())
            .at(path)?;
    }
    Ok(())
}

/// Calls `visit(key low bits, encoded color set)` for every key of a sorted partition.
fn for_each_key(
    records: &[u64],
    mut visit: impl FnMut(u64, &[u8]) -> Result<(), IndexError>,
) -> Result<(), IndexError> {
    let mut entries = Vec::new();
    let mut encoded = Vec::new();
    for group in records.chunk_by(|a, b| a >> 32 == b >> 32) {
        entries.clear();
        entries.extend(group.iter().map(|record| *record as u32));
        encoded.clear();
        encode_color_set(&entries, &mut encoded);
        visit(group[0] >> 32, &encoded)?;
    }
    Ok(())
}

/// Counts the keys of every distinct color set.
fn count_color_sets(spill: &Spill) -> Result<SetCounts, IndexError> {
    let stripes: Vec<Mutex<HashMap<Box<[u8]>, u64>>> =
        (0..COUNT_STRIPES).map(|_| Mutex::default()).collect();
    (0..spill.paths.len())
        .into_par_iter()
        .with_max_len(1)
        .try_for_each(|partition| {
            let records = spill.read(partition)?;
            let mut local: HashMap<Box<[u8]>, u64> = HashMap::new();
            for_each_key(&records, |_, set| {
                match local.get_mut(set) {
                    Some(count) => *count += 1,
                    None => drop(local.insert(set.into(), 1)),
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

/// Encodes the bucket blocks of one partition.
fn encode_partition(
    spill: &Spill,
    partition: usize,
    ids: &HashMap<&[u8], u32>,
) -> Result<EncodedPartition, IndexError> {
    let records = spill.read(partition)?;
    let mut seeds: Vec<(u64, u32)> = Vec::new();
    for_each_key(&records, |low, set| {
        let id = ids
            .get(set)
            .ok_or_else(|| IndexError::format("a color set changed between two reads"))?;
        seeds.push((low, *id));
        Ok(())
    })?;
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
        blocks.push(EncodedBucket {
            bucket,
            length: bytes.len() - start,
            seeds: group.len() as u64,
        });
    }
    Ok(EncodedPartition { bytes, blocks })
}
