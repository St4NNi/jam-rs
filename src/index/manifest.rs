//! The manifest at the end of the seeds file of a shard (JSON): parameters, entries, and the size
//! and SHA-256 digest of every packed sequence file.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::IndexError;
use super::seed_file::{HEADER_BYTES, PAGE_BITS, header_manifest, header_sections};
use super::{MAX_BUCKET_BITS, MIN_BUCKET_BITS};
use crate::seeds::SeedScheme;
use crate::storage::ReadAt;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::ops::Range;
/// The directory of the packed sequences in a shard directory.
pub(crate) const PACKED_DIRECTORY: &str = "packed";
/// The format version that this code reads and writes.
pub(crate) const FORMAT_VERSION: u32 = 2;

/// The parameters of a shard.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Parameters {
    /// How seeds are chosen: k, w and the entropy limit.
    #[serde(flatten)]
    pub scheme: SeedScheme,
    /// Bases per bin, as a power of two.
    pub bin_bits: u8,
    /// Sketch buckets, as a power of two; the build chooses them from the shard size.
    pub bucket_bits: u8,
    /// Bytes of a page of the seeds file, as a power of two; no bucket crosses a page.
    pub page_bits: u8,
}

impl Default for Parameters {
    fn default() -> Self {
        Self {
            scheme: SeedScheme::default(),
            bin_bits: 10,
            bucket_bits: 24,
            page_bits: 16,
        }
    }
}

impl Parameters {
    /// Checks that every value lies in its range.
    pub fn validate(&self) -> Result<(), IndexError> {
        self.scheme.validate()?;
        if !(4..=20).contains(&self.bin_bits) {
            return Err(IndexError::Input("bin bits must be 4 to 20".into()));
        }
        if !(MIN_BUCKET_BITS..=MAX_BUCKET_BITS).contains(&self.bucket_bits) {
            return Err(IndexError::Input(format!(
                "bucket bits must be {MIN_BUCKET_BITS} to {MAX_BUCKET_BITS}"
            )));
        }
        if !PAGE_BITS.contains(&u32::from(self.page_bits)) {
            return Err(IndexError::Input("pages must be 64 bytes to 64 KiB".into()));
        }
        Ok(())
    }
}

/// One file of a shard with its size and digest.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FileRecord {
    /// Path relative to the shard directory, or to the packed directory for packed files.
    pub path: String,
    /// Size in bytes.
    pub bytes: u64,
    /// SHA-256 digest, lower case hex.
    pub sha256: String,
}

/// One entry of a shard.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EntryRecord {
    /// Entry name, unique in the shard.
    pub name: String,
    /// Contigs of the entry.
    pub contigs: u64,
    /// Bases of all contigs.
    pub bases: u64,
    /// Seed locations of the entry: its distinct pairs of seed and seed location.
    pub seed_locations: u64,
    /// The packed sequence file; none in a shard built with the sketch part only.
    pub packed: Option<FileRecord>,
}

/// The manifest of a shard.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Manifest {
    /// Format version of the shard files.
    pub format_version: u32,
    /// Parameters of the shard.
    pub parameters: Parameters,
    /// Entries in entry id order.
    pub entries: Vec<EntryRecord>,
    /// SHA-256 digest of the seeds: the bytes of the seeds file from its header to the manifest.
    pub seeds_sha256: String,
    /// The bytes of the sketch part in the seeds file.
    pub sketch_part: Range<u64>,
    /// The bytes of the locator part in the seeds file, if it has one.
    pub locator_part: Option<Range<u64>>,
    /// Stored seeds of the seeds file: distinct seed keys.
    pub stored_seeds: u64,
    /// Distinct color sets of the seeds file.
    pub color_sets: u64,
    /// Seed locations of all entries.
    pub seed_locations: u64,
}

impl Manifest {
    /// Parses and checks the bytes of a manifest.
    pub fn parse(bytes: &[u8]) -> Result<Self, IndexError> {
        let manifest: Self = serde_json::from_slice(bytes)
            .map_err(|error| IndexError::format(format!("manifest: {error}")))?;
        if manifest.format_version != FORMAT_VERSION {
            return Err(IndexError::format(format!(
                "manifest format version {} is not {FORMAT_VERSION}",
                manifest.format_version
            )));
        }
        manifest.parameters.validate()?;
        check_names(manifest.entries.iter().map(|entry| entry.name.as_str()))?;
        Ok(manifest)
    }

    /// The bytes of the manifest.
    pub fn to_bytes(&self) -> Result<Vec<u8>, IndexError> {
        serde_json::to_vec_pretty(self)
            .map_err(|error| IndexError::format(format!("manifest: {error}")))
    }

    /// Reads the manifest at the end of the seeds file `source`, checks it, and gives the
    /// SHA-256 digest of its bytes, which names the shard, and the size of the seeds file.
    pub fn read_from(source: &dyn ReadAt) -> Result<(Self, String, u64), IndexError> {
        let mut header = [0u8; HEADER_BYTES as usize];
        source.read_exact_at(0, &mut header)?;
        header_sections(&header)?;
        let (range, crc) = header_manifest(&header);
        let length = usize::try_from(range.end - range.start)
            .map_err(|_| IndexError::format("seeds file manifest"))?;
        let mut bytes = vec![0u8; length];
        source.read_exact_at(range.start, &mut bytes)?;
        if length == 0 || crc32fast::hash(&bytes) != crc {
            return Err(IndexError::Checksum("seeds file manifest".into()));
        }
        Ok((Self::parse(&bytes)?, sha256_hex(&[&bytes]), range.end))
    }
}

/// The fewest bucket bits (8 to 28) that hold at most 256 seeds per bucket on average, given an
/// upper bound of the distinct seeds. A query k-mer then matches a stored seed by chance at most
/// about 256 / 65,536 = 0.4% of the time, and the directory costs at most 8 / 256 bytes per seed.
pub(crate) fn choose_bucket_bits(seeds: u64) -> u8 {
    (MIN_BUCKET_BITS..MAX_BUCKET_BITS)
        .find(|bits| seeds <= 256 << bits)
        .unwrap_or(MAX_BUCKET_BITS)
}

/// The packed file path of an entry, relative to the packed directory.
pub(crate) fn packed_file_name(entry: &str) -> String {
    format!("{entry}.packed")
}

/// Checks that entry names are valid file names and distinct.
pub(crate) fn check_names<'a>(names: impl Iterator<Item = &'a str>) -> Result<(), IndexError> {
    let mut seen = std::collections::HashSet::new();
    for name in names {
        let valid = !name.is_empty()
            && name != "."
            && name != ".."
            && !name
                .bytes()
                .any(|byte| matches!(byte, b'/' | b'\\' | 0) || byte.is_ascii_whitespace());
        if !valid {
            return Err(IndexError::Input(format!(
                "{name:?} is not a valid entry name: it must not be empty or hold a slash or white space"
            )));
        }
        if !seen.insert(name) {
            return Err(IndexError::Input(format!(
                "the entry name {name} is used twice"
            )));
        }
    }
    Ok(())
}

/// The SHA-256 digest of `parts` one after the other, lower case hex.
pub(crate) fn sha256_hex(parts: &[&[u8]]) -> String {
    let mut digest = Sha256::new();
    for part in parts {
        digest.update(part);
    }
    hex::encode(digest.finalize())
}

/// Checks a file read from `source` against its record, for files that may not be local.
pub(crate) fn verify_source(source: &dyn ReadAt, record: &FileRecord) -> Result<(), IndexError> {
    let size = source.size();
    if size != record.bytes || sha256_range(source, 0..size)? != record.sha256 {
        return Err(IndexError::Checksum(format!(
            "{} differs from its manifest record",
            record.path
        )));
    }
    Ok(())
}

/// The SHA-256 digest of the bytes `range` of `source`, lower case hex.
pub(crate) fn sha256_range(source: &dyn ReadAt, range: Range<u64>) -> Result<String, IndexError> {
    let mut digest = Sha256::new();
    let mut buffer = vec![0u8; 1 << 20];
    let mut offset = range.start;
    while offset < range.end {
        let length = (range.end - offset).min(buffer.len() as u64) as usize;
        source.read_exact_at(offset, &mut buffer[..length])?;
        digest.update(&buffer[..length]);
        offset += length as u64;
    }
    Ok(hex::encode(digest.finalize()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn file(path: &str) -> FileRecord {
        FileRecord {
            path: path.into(),
            bytes: 3,
            sha256: sha256_hex(&[b"abc"]),
        }
    }

    #[test]
    fn manifest_round_trip() {
        let manifest = Manifest {
            format_version: FORMAT_VERSION,
            parameters: Parameters {
                scheme: SeedScheme {
                    min_entropy: 1.1,
                    ..SeedScheme::default()
                },
                ..Parameters::default()
            },
            entries: vec![EntryRecord {
                name: "entry_1".into(),
                contigs: 2,
                bases: 1000,
                seed_locations: 90,
                packed: Some(file("entry_1.packed")),
            }],
            seeds_sha256: sha256_hex(&[b""]),
            sketch_part: 4096..8192,
            locator_part: None,
            stored_seeds: 10,
            color_sets: 1,
            seed_locations: 11,
        };
        let bytes = manifest.to_bytes().unwrap();
        assert_eq!(Manifest::parse(&bytes).unwrap(), manifest);
        let record = file("abc");
        assert_eq!(
            record.sha256,
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        verify_source(&b"abc".to_vec(), &record).unwrap();
        assert!(verify_source(&b"abd".to_vec(), &record).is_err());
    }

    #[test]
    fn names_parameters_checked() {
        assert!(check_names(["a", "b"].into_iter()).is_ok());
        for bad in [&["a", "a"][..], &[""], &["a b"], &["a/b"], &[".."]] {
            assert!(check_names(bad.iter().copied()).is_err(), "{bad:?}");
        }
        let parameters = Parameters::default();
        assert!(parameters.validate().is_ok());
        assert_eq!(choose_bucket_bits(0), 8);
        assert_eq!(choose_bucket_bits(256 << 8), 8);
        assert_eq!(choose_bucket_bits((256 << 8) + 1), 9);
        assert_eq!(choose_bucket_bits(17_000_000_000), 26);
        assert_eq!(choose_bucket_bits(u64::MAX), 28);
        for bad in [
            Parameters {
                bin_bits: 3,
                ..parameters
            },
            Parameters {
                bucket_bits: 29,
                ..parameters
            },
            Parameters {
                page_bits: 17,
                ..parameters
            },
            Parameters {
                scheme: SeedScheme {
                    k: 33,
                    ..parameters.scheme
                },
                ..parameters
            },
        ] {
            assert!(bad.validate().is_err());
        }
    }
}
