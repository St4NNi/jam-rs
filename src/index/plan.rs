//! The plan of a warmed cache directory: which byte ranges of a shard's files are pinned (held
//! as local files, never evicted), and whether each is read into RAM or from the local disk.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::packed::BLOCK_BYTES;
use super::sources::IndexHead;
use super::{AtPath, IndexError};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// The file name of the plan in the directory of a shard in the cache directory.
pub(crate) const PLAN_FILE: &str = "plan.json";

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

impl std::fmt::Display for Part {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(match self {
            Self::Sketch => "sketch",
            Self::Locators => "locators",
            Self::Packed(_) => "packed",
        })
    }
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
    /// The part it belongs to: `sketch`, `locators` or `packed`.
    pub part: String,
    /// The file, relative to the shard: `seeds` or a sequence file.
    pub file: String,
    /// The first byte of the range.
    pub start: u64,
    /// The end of the range (exclusive).
    pub end: u64,
    /// The SHA-256 digest of the whole file, checked when the range covers it; for the seeds
    /// file the digest of its manifest.
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
    /// The SHA-256 digest of the shard's manifest.
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

    /// The plan for the shard of `head` within `memory` and `disk` bytes: `parts`, or else the
    /// sketch and the locators while the budgets last, each in RAM if it fits, else on disk.
    /// Packed sequences go on disk; a pin in RAM also takes disk bytes.
    pub fn choose(
        head: &IndexHead,
        memory: u64,
        disk: u64,
        parts: &[Part],
    ) -> Result<Self, IndexError> {
        let mut plan = Self {
            shard: head.digest.clone(),
            memory,
            disk,
            pins: Vec::new(),
        };
        let chosen = if parts.is_empty() {
            &[Part::Sketch, Part::Locators][..]
        } else {
            parts
        };
        for part in chosen {
            let pins = part_pins(head, part, Place::Ssd)?;
            let bytes: u64 = pins.iter().map(Pin::bytes).sum();
            let (memory_left, disk_left) = plan.cache_budget();
            let place = match part {
                _ if bytes > disk_left && parts.is_empty() => continue,
                _ if bytes > disk_left => {
                    return Err(IndexError::Input(format!(
                        "pinning {part} needs {bytes} bytes of disk, {disk_left} are left"
                    )));
                }
                Part::Packed(_) => Place::Ssd,
                _ if bytes <= memory_left => Place::Ram,
                _ => Place::Ssd,
            };
            plan.pins
                .extend(pins.into_iter().map(|pin| Pin { place, ..pin }));
        }
        Ok(plan)
    }

    /// The RAM and disk bytes that the pins leave of the budgets, for the block cache.
    pub fn cache_budget(&self) -> (u64, u64) {
        let ram: u64 = (self.pins.iter())
            .filter(|pin| pin.place == Place::Ram)
            .map(Pin::bytes)
            .sum();
        let disk: u64 = self.pins.iter().map(Pin::bytes).sum();
        (
            self.memory.saturating_sub(ram),
            self.disk.saturating_sub(disk),
        )
    }

    /// The local file of a pin in `cache`, by part and range, as one file may hold several.
    pub fn pin_path(&self, cache: &Path, pin: &Pin) -> PathBuf {
        let file = pin.file.replace('%', "%25").replace('/', "%2F");
        let name = format!("{file}.{}-{}", pin.start, pin.end);
        cache.join(&self.shard).join(&pin.part).join(name)
    }

    /// The file next to a pin that holds the ETag its bytes were downloaded under.
    #[cfg(test)]
    pub fn etag_path(&self, cache: &Path, pin: &Pin) -> PathBuf {
        etag_file(&self.pin_path(cache, pin))
    }
}

/// The file next to the pin file `path` that holds the ETag of its bytes.
pub(crate) fn etag_file(path: &Path) -> PathBuf {
    PathBuf::from(format!("{}.etag", path.display()))
}

/// The byte ranges that hold `part`, as pins in `place`. This is the one place that maps part
/// names to files.
pub fn part_pins(head: &IndexHead, part: &Part, place: Place) -> Result<Vec<Pin>, IndexError> {
    let manifest = &head.manifest;
    let seeds = |range: &std::ops::Range<u64>| Pin {
        part: part.to_string(),
        file: head.file.clone(),
        start: range.start,
        end: range.end,
        sha256: head.digest.clone(),
        place,
    };
    // The packed sequences of the entries in each sequence file, joined where they follow each
    // other; a shard built with the sketch part only has none.
    let packed = |entries: &mut dyn Iterator<Item = &super::EntryRecord>| {
        let mut sequences: Vec<_> = entries
            .filter_map(|entry| entry.sequence.as_ref())
            .collect();
        sequences.sort_by_key(|sequence| (sequence.file, sequence.offset));
        let mut pins: Vec<Pin> = Vec::new();
        for sequence in sequences {
            let file = &manifest.sequence_files[sequence.file];
            match pins.last_mut() {
                Some(last)
                    if last.file == file.path
                        && last.end.next_multiple_of(BLOCK_BYTES) == sequence.offset =>
                {
                    last.end = sequence.range().end;
                }
                _ => pins.push(Pin {
                    file: file.path.clone(),
                    sha256: file.sha256.clone(),
                    ..seeds(&sequence.range())
                }),
            }
        }
        pins
    };
    Ok(match part {
        // The header and the bucket table of the seeds file belong to the sketch part.
        Part::Sketch => vec![seeds(&(0..manifest.sketch_part.end))],
        Part::Locators => manifest.locator_part.iter().map(seeds).collect(),
        Part::Packed(None) => packed(&mut manifest.entries.iter()),
        Part::Packed(Some(names)) => {
            let mut pins = Vec::new();
            for name in names {
                let record = (manifest.entries.iter().find(|entry| &entry.name == name))
                    .ok_or_else(|| IndexError::Input(format!("no entry is named {name}")))?;
                pins.extend(packed(&mut std::iter::once(record)));
            }
            pins
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::{Shard, build_shard};
    use crate::test_support::{build_options, test_shard, write_entries};

    #[test]
    fn budgets_place_parts() {
        let directory = tempfile::tempdir().unwrap();
        let entries = vec![vec![vec![b'A'; 2_000]], vec![vec![b'C'; 3_000]]];
        let shard = test_shard(directory.path(), &entries);
        let manifest = shard.manifest();
        let locator_part = manifest.locator_part.clone().unwrap();
        let (sketch, locators) = (
            manifest.sketch_part.end,
            locator_part.end - locator_part.start,
        );
        let places = |memory, disk, parts: &[Part]| {
            let plan = Plan::choose(shard.head(), memory, disk, parts).unwrap();
            let places: Vec<(String, Place)> = (plan.pins.iter())
                .map(|pin| (pin.part.clone(), pin.place))
                .collect();
            (places, plan.cache_budget())
        };
        let total = sketch + locators;
        let (pinned, budget) = places(total, total + 100, &[]);
        let both = vec![
            ("sketch".to_string(), Place::Ram),
            ("locators".to_string(), Place::Ram),
        ];
        assert_eq!((pinned, budget), (both, (0, 100)));
        let (pinned, budget) = places(sketch, total, &[]);
        assert_eq!(pinned[1], ("locators".to_string(), Place::Ssd));
        assert_eq!(budget, (0, 0));
        let (pinned, budget) = places(10, sketch, &[]);
        assert_eq!(pinned, vec![("sketch".to_string(), Place::Ssd)]);
        assert_eq!(budget, (10, 0));
        assert!(places(0, 0, &[]).0.is_empty());
        let (pinned, _) = places(1 << 30, 1 << 30, &[Part::Packed(None)]);
        assert!(pinned.iter().all(|(_, place)| *place == Place::Ssd));
        assert_eq!(pinned.len(), 1);
        let error = Plan::choose(shard.head(), 1 << 30, 10, &[Part::Sketch]).unwrap_err();
        assert!(error.to_string().starts_with("pinning sketch needs"));
    }

    /// Pins are byte ranges: the seeds parts of the index file, a whole chunk for all packed
    /// sequences, and the exact range of a named entry.
    #[test]
    fn pins_are_ranges() {
        let directory = tempfile::tempdir().unwrap();
        let entries: Vec<Vec<Vec<u8>>> = (0..3).map(|e| vec![vec![b"ACG"[e]; 9_000]]).collect();
        let inputs = write_entries(directory.path(), &entries);
        let options = crate::index::BuildOptions {
            chunks: 2,
            ..build_options(directory.path(), 2)
        };
        build_shard(&inputs, &directory.path().join("shard"), &options).unwrap();
        let shard = Shard::open(&directory.path().join("shard.jam")).unwrap();
        let (head, manifest) = (shard.head(), shard.manifest());
        let ranges = |part: &Part| {
            let pins = part_pins(head, part, Place::Ssd).unwrap();
            let ranges = pins
                .iter()
                .map(|pin| (pin.file.clone(), pin.start..pin.end));
            ranges.collect::<Vec<_>>()
        };
        let index = "shard.jam".to_string();
        assert_eq!(
            ranges(&Part::Sketch),
            [(index.clone(), 0..manifest.sketch_part.end)]
        );
        let locators = manifest.locator_part.clone().unwrap();
        assert_eq!(ranges(&Part::Locators), [(index, locators)]);
        let files = &manifest.sequence_files;
        let whole: Vec<_> = files.iter().map(|f| (f.path.clone(), 0..f.bytes)).collect();
        assert_eq!(ranges(&Part::Packed(None)), whole);
        let sequence = manifest.entries[2].sequence.clone().unwrap();
        let file = files[sequence.file].path.clone();
        let named = Part::Packed(Some(vec!["entry2".into()]));
        assert_eq!(ranges(&named), [(file, sequence.range())]);
        let missing = Part::Packed(Some(vec!["other".into()]));
        assert!(part_pins(head, &missing, Place::Ssd).is_err());
    }

    #[test]
    fn plan_round_trip() {
        let directory = tempfile::tempdir().unwrap();
        let entries = vec![vec![vec![b'A'; 2_000]], vec![vec![b'C'; 3_000]]];
        let shard = test_shard(directory.path(), &entries);
        let manifest = shard.manifest();
        let mut pins = part_pins(shard.head(), &Part::Sketch, Place::Ram).unwrap();
        pins.extend(part_pins(shard.head(), &Part::Locators, Place::Ssd).unwrap());
        assert_eq!(
            part_pins(shard.head(), &Part::Packed(None), Place::Ssd)
                .unwrap()
                .len(),
            1
        );
        let named = Part::Packed(Some(vec!["entry1".into()]));
        pins.extend(part_pins(shard.head(), &named, Place::Ssd).unwrap());
        assert_eq!(pins[2].file, "shard.jamseq");
        assert_eq!(pins[0].bytes(), manifest.sketch_part.end);
        let missing = Part::Packed(Some(vec!["other".into()]));
        assert!(part_pins(shard.head(), &missing, Place::Ssd).is_err());
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
        let pin = &plan.pins[2];
        let path = plan.pin_path(&cache, pin);
        let name = format!("packed/shard.jamseq.{}-{}", pin.start, pin.end);
        assert!(path.ends_with(name));
        let other = "0".repeat(64);
        std::fs::rename(cache.join(shard.digest()), cache.join(&other)).unwrap();
        assert!(Plan::read(&cache, &other).is_err());
    }
}
