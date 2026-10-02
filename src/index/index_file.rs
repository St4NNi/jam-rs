//! The list `index.json` of a resharded index, which names its index files in order.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::sources::Location;
use super::{AtPath, IndexError};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// The file name of the list in an index directory.
pub const INDEX_FILE: &str = "index.json";
/// The format version of the list that this code reads and writes.
pub const INDEX_FORMAT_VERSION: u32 = 1;

/// One shard of an index.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ShardRecord {
    /// The index file of the shard, relative to the index directory.
    pub file: String,
    /// Index entry id of the shard's first entry.
    pub first_entry: u64,
    /// Entries of the shard.
    pub entries: u64,
    /// Seed locations of the shard, which the size of its seeds file follows.
    pub seed_locations: u64,
    /// SHA-256 digest of the shard's manifest, lower case hex.
    pub manifest_sha256: String,
}

/// The list of the index files of an index, in order.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct IndexFile {
    /// Format version of the index file.
    pub format_version: u32,
    /// Shards in index order.
    pub shards: Vec<ShardRecord>,
}

impl IndexFile {
    /// Reads the list of the index in `directory` and checks that its entry ranges follow each
    /// other without gaps.
    pub fn read(directory: &Path) -> Result<Self, IndexError> {
        let path = directory.join(INDEX_FILE);
        let bytes = std::fs::read(&path).at(&path)?;
        let index: Self = serde_json::from_slice(&bytes)
            .map_err(|error| IndexError::format(format!("index file: {error}")))?;
        if index.format_version != INDEX_FORMAT_VERSION {
            return Err(IndexError::format(format!(
                "index file format version {} is not {INDEX_FORMAT_VERSION}",
                index.format_version
            )));
        }
        let mut next = 0u64;
        for shard in &index.shards {
            let plain = Path::new(&shard.file)
                .components()
                .all(|part| matches!(part, std::path::Component::Normal(_)));
            if shard.first_entry != next || shard.entries == 0 || !plain {
                return Err(IndexError::format("index file shard list"));
            }
            next += shard.entries;
        }
        Ok(index)
    }

    /// Writes the list into `directory`.
    pub fn write(&self, directory: &Path) -> Result<(), IndexError> {
        let bytes = serde_json::to_vec_pretty(self)
            .map_err(|error| IndexError::format(format!("index file: {error}")))?;
        let path = directory.join(INDEX_FILE);
        let mut file = std::fs::File::create_new(&path).at(&path)?;
        std::io::Write::write_all(&mut file, &bytes).at(&path)?;
        file.sync_all().at(&path)?;
        Ok(())
    }

    /// The index files of the index in `directory`, in order.
    pub fn shard_files(&self, directory: &Path) -> Vec<PathBuf> {
        let join = |shard: &ShardRecord| directory.join(&shard.file);
        self.shards.iter().map(join).collect()
    }

    /// Checks that the manifest of every index file still has its recorded digest.
    pub fn verify(&self, directory: &Path) -> Result<(), IndexError> {
        for shard in &self.shards {
            let location = Location::File(directory.join(&shard.file));
            if location.head()?.digest != shard.manifest_sha256 {
                return Err(IndexError::Checksum(format!(
                    "the manifest of {} differs from the list",
                    shard.file
                )));
            }
        }
        Ok(())
    }

    /// Entries of the whole index.
    pub fn entries(&self) -> u64 {
        self.shards.iter().map(|shard| shard.entries).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(file: &str, first_entry: u64, entries: u64) -> ShardRecord {
        ShardRecord {
            file: file.into(),
            first_entry,
            entries,
            seed_locations: 10 * entries,
            manifest_sha256: "00".into(),
        }
    }

    #[test]
    fn index_round_trip() {
        let directory = tempfile::tempdir().unwrap();
        let index = IndexFile {
            format_version: INDEX_FORMAT_VERSION,
            shards: vec![
                record("shard-0000.jam", 0, 3),
                record("shard-0001.jam", 3, 2),
            ],
        };
        index.write(directory.path()).unwrap();
        assert_eq!(IndexFile::read(directory.path()).unwrap(), index);
        assert_eq!(index.entries(), 5);
        assert_eq!(
            index.shard_files(Path::new("/x"))[1],
            Path::new("/x/shard-0001.jam")
        );
        assert!(index.verify(directory.path()).is_err());
    }

    #[test]
    fn bad_shard_lists() {
        for shards in [
            vec![record("a", 0, 3), record("b", 4, 2)],
            vec![record("a", 0, 0)],
            vec![record("../a", 0, 1)],
            vec![record("/a", 0, 1)],
        ] {
            let directory = tempfile::tempdir().unwrap();
            let index = IndexFile {
                format_version: INDEX_FORMAT_VERSION,
                shards,
            };
            index.write(directory.path()).unwrap();
            assert!(IndexFile::read(directory.path()).is_err());
        }
    }
}
