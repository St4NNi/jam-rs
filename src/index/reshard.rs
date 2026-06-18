//! Reshards an index: any number of input shards into any number of output shards. Entries keep
//! the input order and are cut into contiguous groups, balanced by their seed locations (a
//! shard's sketch RAM follows them). Every entry keeps its positions section and its packed
//! sequence unchanged; each output shard builds only its sketch, with bucket bits chosen from its
//! own seed locations. One output shard joins the inputs; one input into several splits it. The
//! groups depend only on the inputs and the split, never on the thread count.

use super::build::{complete_shard, in_partial_directory};
use super::index_file::{INDEX_FILE, INDEX_FORMAT_VERSION, IndexFile, ShardRecord};
use super::manifest::{MANIFEST_FILE, Manifest, PACKED_DIRECTORY, Parameters, check_names};
use super::positions::{POSITIONS_FILE, Positions};
use super::{IndexError, MappedFile};
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
    /// Threads that build each sketch.
    pub threads: usize,
    /// Directory for temporary files.
    pub temp_dir: PathBuf,
    /// Bucket bits to use instead of the ones chosen from each output shard's size.
    pub bucket_bits: Option<u8>,
}

/// One input shard, opened for copying its entries.
struct InputShard {
    directory: PathBuf,
    manifest: Manifest,
    positions: Positions,
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
        sources.extend((0..input.positions.entries()).map(|entry| (shard, entry)));
    }
    let weights = sources
        .iter()
        .map(|&(shard, entry)| shards[shard].positions.entry_locations(entry))
        .collect::<Result<Vec<_>, IndexError>>()?;
    let groups = split_groups(&weights, options.split)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.threads.max(1))
        .build()
        .map_err(|error| IndexError::Input(format!("thread pool: {error}")))?;
    in_partial_directory(output, |directory| {
        let mut records = Vec::with_capacity(groups.len());
        for (number, group) in groups.iter().enumerate() {
            let name = format!("shard-{number:04}");
            let part = &sources[group.clone()];
            let shard_directory = directory.join(&name);
            std::fs::create_dir(&shard_directory)?;
            let manifest =
                pool.install(|| write_shard(&shards, part, &shard_directory, options))?;
            records.push(ShardRecord {
                directory: name,
                first_entry: group.start as u64,
                entries: group.len() as u64,
                seed_locations: manifest.seed_locations,
                manifest_sha256: super::manifest::sha256_file(
                    &shard_directory.join(MANIFEST_FILE),
                )?
                .1,
            });
        }
        let index = IndexFile {
            format_version: INDEX_FORMAT_VERSION,
            shards: records,
        };
        index.write(directory)?;
        Ok(index)
    })
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
            let file = MappedFile::open(&directory.join(POSITIONS_FILE), false)?;
            let positions = Positions::open(Arc::new(file))?;
            if positions.entries() as usize != manifest.entries.len() {
                return Err(IndexError::format(format!(
                    "the positions of {} differ from its manifest",
                    directory.display()
                )));
            }
            Ok(InputShard {
                directory,
                manifest,
                positions,
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

/// Writes one output shard from the entries `part` (input shard, entry) into `directory`.
fn write_shard(
    shards: &[InputShard],
    part: &[(usize, u32)],
    directory: &Path,
    options: &ReshardOptions,
) -> Result<Manifest, IndexError> {
    let packed = directory.join(PACKED_DIRECTORY);
    std::fs::create_dir(&packed)?;
    let mut entries = Vec::with_capacity(part.len());
    for &(shard, entry) in part {
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
                std::os::unix::fs::symlink(std::fs::canonicalize(&source)?, &target)?
            }
            PackedOutput::Copy => drop(std::fs::copy(&source, &target)?),
        }
        entries.push(record);
    }
    let parameters = shards[0].manifest.parameters;
    complete_shard(
        directory,
        parameters,
        entries,
        &options.temp_dir,
        options.bucket_bits,
        |writer| {
            for &(shard, entry) in part {
                let positions = &shards[shard].positions;
                // The section's CRC-32 is checked as it is read, so damage is not copied.
                let bytes = positions.section_bytes(entry)?;
                writer.push_section(positions.section(entry)?, &bytes)?;
            }
            Ok(())
        },
    )
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
    fn groups_are_contiguous_and_balanced() {
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
