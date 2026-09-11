//! The plan of a warmed cache directory: which byte ranges of a shard's files are pinned (held
//! as local files, never evicted), and whether each is read into RAM or from the local disk.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::manifest::{Manifest, PACKED_DIRECTORY};
use super::{AtPath, IndexError};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// The file name of the plan in the directory of a shard in the cache directory.
const PLAN_FILE: &str = "plan.json";

/// A part of a shard that can be pinned.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Part {
    /// The seeds and color sets of stage 1.
    Sketch,
    /// The seed locations of stage 2.
    Locators,
    /// Packed sequences of stage 3: all entries, or the named ones.
    Packed(Option<Vec<String>>),
}

/// Where a pinned range is read from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Place {
    /// Read into RAM at open.
    Ram,
    /// Read from the local disk.
    Ssd,
}

/// A byte range of one shard file, held locally.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Pin {
    /// The file, relative to the shard: `sketch`, `positions` or `packed/NAME.packed`.
    pub file: String,
    /// The first byte of the range.
    pub start: u64,
    /// The end of the range (exclusive).
    pub end: u64,
    /// The SHA-256 digest of the whole file, checked when the range covers it.
    pub sha256: String,
    /// Where the range is read from.
    pub place: Place,
}

impl Pin {
    /// Bytes of the range.
    pub fn bytes(&self) -> u64 {
        self.end - self.start
    }
}

/// The plan of one shard in a cache directory. It lies in the directory named by the shard's
/// manifest digest, next to the pinned files.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Plan {
    /// The SHA-256 digest of the shard's manifest file.
    pub shard: String,
    /// The RAM budget in bytes, for pins and the block cache.
    pub memory: u64,
    /// The disk budget in bytes, for pins and the block cache.
    pub disk: u64,
    /// The pinned ranges.
    pub pins: Vec<Pin>,
}

impl Plan {
    /// The plan of the shard with manifest digest `shard` in `cache`, or `None` before a warm.
    pub fn read(cache: &Path, shard: &str) -> Result<Option<Self>, IndexError> {
        let path = cache.join(shard).join(PLAN_FILE);
        let bytes = match std::fs::read(&path) {
            Ok(bytes) => bytes,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
            Err(error) => return Err(error).at(&path),
        };
        let plan: Self = serde_json::from_slice(&bytes)
            .map_err(|error| IndexError::format(format!("{}: {error}", path.display())))?;
        if plan.shard != shard {
            return Err(IndexError::format(format!(
                "{} belongs to another shard",
                path.display()
            )));
        }
        Ok(Some(plan))
    }

    /// Writes the plan into `cache`, replacing an older one at once.
    pub fn write(&self, cache: &Path) -> Result<(), IndexError> {
        let directory = cache.join(&self.shard);
        std::fs::create_dir_all(&directory).at(&directory)?;
        let bytes = serde_json::to_vec_pretty(self)
            .map_err(|error| IndexError::format(format!("plan: {error}")))?;
        let path = directory.join(PLAN_FILE);
        let partial = directory.join(format!("{PLAN_FILE}.partial"));
        std::fs::write(&partial, bytes).at(&partial)?;
        std::fs::rename(&partial, &path).at(&path)
    }

    /// The local file of a pin in `cache`.
    pub fn pin_path(&self, cache: &Path, pin: &Pin) -> PathBuf {
        cache.join(&self.shard).join(&pin.file)
    }
}

/// The byte ranges that hold `part`, as pins in `place`. This is the one place that maps part
/// names to files.
pub fn part_pins(manifest: &Manifest, part: &Part, place: Place) -> Result<Vec<Pin>, IndexError> {
    let whole = |file: String, record: &super::FileRecord| Pin {
        file,
        start: 0,
        end: record.bytes,
        sha256: record.sha256.clone(),
        place,
    };
    let packed = |record: &super::EntryRecord| {
        let file = format!("{PACKED_DIRECTORY}/{}", record.packed.path);
        whole(file, &record.packed)
    };
    Ok(match part {
        Part::Sketch => vec![whole(manifest.sketch.path.clone(), &manifest.sketch)],
        Part::Locators => vec![whole(manifest.positions.path.clone(), &manifest.positions)],
        Part::Packed(None) => manifest.entries.iter().map(packed).collect(),
        Part::Packed(Some(names)) => {
            let mut pins = Vec::new();
            for name in names {
                let record = (manifest.entries.iter().find(|entry| &entry.name == name))
                    .ok_or_else(|| IndexError::Input(format!("no entry is named {name}")))?;
                pins.push(packed(record));
            }
            pins
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::test_shard;

    #[test]
    fn plan_round_trip() {
        let directory = tempfile::tempdir().unwrap();
        let entries = vec![vec![vec![b'A'; 2_000]], vec![vec![b'C'; 3_000]]];
        let shard = test_shard(directory.path(), &entries);
        let manifest = shard.manifest();
        let mut pins = part_pins(manifest, &Part::Sketch, Place::Ram).unwrap();
        pins.extend(part_pins(manifest, &Part::Locators, Place::Ssd).unwrap());
        assert_eq!(
            part_pins(manifest, &Part::Packed(None), Place::Ssd)
                .unwrap()
                .len(),
            2
        );
        let named = Part::Packed(Some(vec!["entry1".into()]));
        pins.extend(part_pins(manifest, &named, Place::Ssd).unwrap());
        assert_eq!(pins[2].file, "packed/entry1.packed");
        assert_eq!(pins[0].bytes(), manifest.sketch.bytes);
        let missing = Part::Packed(Some(vec!["other".into()]));
        assert!(part_pins(manifest, &missing, Place::Ssd).is_err());
        let cache = directory.path().join("cache");
        assert_eq!(Plan::read(&cache, shard.digest()).unwrap(), None);
        let plan = Plan {
            shard: shard.digest().to_string(),
            memory: 1 << 30,
            disk: 1 << 40,
            pins,
        };
        plan.write(&cache).unwrap();
        plan.write(&cache).unwrap();
        assert_eq!(
            Plan::read(&cache, shard.digest()).unwrap(),
            Some(plan.clone())
        );
        let path = plan.pin_path(&cache, &plan.pins[2]);
        assert!(path.ends_with("packed/entry1.packed"));
        let other = "0".repeat(64);
        std::fs::rename(cache.join(shard.digest()), cache.join(&other)).unwrap();
        assert!(Plan::read(&cache, &other).is_err());
    }
}
