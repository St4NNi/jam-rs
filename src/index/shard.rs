//! An open shard: its manifest, the sketch in RAM, the positions, and the packed sequences on
//! request, read from a local directory or an S3 prefix.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::IndexError;
use super::manifest::{EntryRecord, MANIFEST_FILE, Manifest, sha256_hex, verify_source};
use super::packed::PackedSequence;
use super::positions::Positions;
use super::sketch::Sketch;
use super::sources::{Location, ShardFiles, ShardOptions};
use serde::Serialize;
use std::path::{Path, PathBuf};
use std::time::Instant;

/// Sizes and counts of a shard, for `jam stats`.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize)]
pub struct ShardStats {
    /// Entries of the shard.
    pub entries: u64,
    /// Contigs of all entries.
    pub contigs: u64,
    /// Bases of all entries.
    pub bases: u64,
    /// Stored seeds of the sketch.
    pub seeds: u64,
    /// Distinct color sets of the sketch.
    pub color_sets: u64,
    /// Seed locations of the positions file.
    pub seed_locations: u64,
    /// Bytes of the sketch file.
    pub sketch_bytes: u64,
    /// Bytes of the positions file.
    pub positions_bytes: u64,
    /// Bytes of all packed sequence files.
    pub packed_bytes: u64,
}

/// A shard opened from a local directory or an S3 prefix.
pub struct Shard {
    files: ShardFiles,
    manifest: Manifest,
    /// The SHA-256 digest of the manifest file, lower case hex.
    digest: String,
    sketch: Sketch,
    positions: Positions,
    /// Wall seconds that opening took.
    open_seconds: f64,
}

impl Shard {
    /// Opens the shard in `directory`, with its packed sequences in `directory/packed`.
    pub fn open(directory: &Path) -> Result<Self, IndexError> {
        let location = Location::Directory(directory.to_path_buf());
        Self::open_with(&location, &ShardOptions::default())
    }

    /// Opens the shard in `directory`, with its packed sequences in `packed_directory`, which
    /// several shards may share.
    pub fn open_with_packed(directory: &Path, packed_directory: &Path) -> Result<Self, IndexError> {
        let options = ShardOptions {
            packed: Some(Location::Directory(packed_directory.to_path_buf())),
            ..ShardOptions::default()
        };
        Self::open_with(&Location::Directory(directory.to_path_buf()), &options)
    }

    /// Opens the shard at `location` as `options` say. Only headers and sizes are checked;
    /// `verify` reads everything.
    pub fn open_with(location: &Location, options: &ShardOptions) -> Result<Self, IndexError> {
        let started = Instant::now();
        let manifest_bytes = location.read_whole(MANIFEST_FILE)?;
        let manifest = Manifest::parse(&manifest_bytes)?;
        let files = ShardFiles::new(location, options, manifest.entries.len());
        let sketch_file = files.open_sketch(&manifest.sketch)?;
        let sketch_bytes = sketch_file.bytes().len() as u64;
        let sketch = Sketch::open(sketch_file)?;
        let positions_file = files.open(&manifest.positions)?;
        let positions_bytes = positions_file.size();
        let positions = Positions::open(positions_file)?;
        let parameters = manifest.parameters;
        let entries = manifest.entries.len() as u64;
        if sketch.entry_count() != entries
            || sketch.bucket_bits() != u32::from(parameters.bucket_bits)
            || sketch.seed_count() != manifest.seeds
            || u64::from(positions.entry_count()) != entries
            || positions.bin_bits() != u32::from(parameters.bin_bits)
            || sketch_bytes != manifest.sketch.bytes
            || positions_bytes != manifest.positions.bytes
            || (manifest.entries.iter().enumerate()).any(|(entry, record)| {
                positions
                    .section(entry as u32)
                    .map(|section| section.contigs)
                    .ok()
                    != Some(record.contigs)
            })
        {
            return Err(IndexError::format(
                "the sketch or positions differ from the manifest",
            ));
        }
        Ok(Self {
            files,
            manifest,
            digest: sha256_hex(&[&manifest_bytes]),
            sketch,
            positions,
            open_seconds: started.elapsed().as_secs_f64(),
        })
    }

    /// The SHA-256 digest of the manifest file, which names this version of the shard.
    pub fn digest(&self) -> &str {
        &self.digest
    }

    /// Wall seconds that opening the shard took.
    pub fn open_seconds(&self) -> f64 {
        self.open_seconds
    }

    /// The manifest of the shard.
    pub fn manifest(&self) -> &Manifest {
        &self.manifest
    }

    /// The sketch (stage 1).
    pub fn sketch(&self) -> &Sketch {
        &self.sketch
    }

    /// The positions (stage 2).
    pub fn positions(&self) -> &Positions {
        &self.positions
    }

    /// The local path of an entry's packed sequence; an error when it lies on S3.
    pub fn packed_path(&self, entry: u32) -> Result<PathBuf, IndexError> {
        let record = self.entry_record(entry)?;
        self.files.packed_path(&record.packed).ok_or_else(|| {
            IndexError::Input(format!("the packed sequence of {} is on S3", record.name))
        })
    }

    /// Opens the packed sequence of an entry (stage 3) with the contig starts of the positions,
    /// and checks it against the manifest. Reading its names checks the starts against the file.
    pub fn packed(&self, entry: u32) -> Result<PackedSequence, IndexError> {
        let record = self.entry_record(entry)?;
        let source = self.files.open_packed(entry as usize, &record.packed)?;
        let size = source.size();
        let starts = self.positions.contig_starts(entry)?;
        let packed = PackedSequence::open_with_starts(source, starts)?;
        if size != record.packed.bytes
            || packed.bases() != record.bases
            || packed.contig_count() as u64 != record.contigs
        {
            return Err(IndexError::format(format!(
                "the packed sequence of {} differs from the manifest",
                record.name
            )));
        }
        Ok(packed)
    }

    /// Checks the SHA-256 digest of every file and the checksum of every section.
    pub fn verify(&self) -> Result<(), IndexError> {
        for record in [&self.manifest.sketch, &self.manifest.positions] {
            verify_source(self.files.open_uncached(record)?.as_ref(), record)?;
        }
        self.sketch.verify_checksums()?;
        self.positions.verify_checksums()?;
        for record in &self.manifest.entries {
            let source = self.files.open_packed_uncached(&record.packed)?;
            verify_source(source.as_ref(), &record.packed)?;
        }
        Ok(())
    }

    /// Sizes and counts of the shard.
    pub fn stats(&self) -> ShardStats {
        let entries = &self.manifest.entries;
        ShardStats {
            entries: entries.len() as u64,
            contigs: entries.iter().map(|entry| entry.contigs).sum(),
            bases: entries.iter().map(|entry| entry.bases).sum(),
            seeds: self.manifest.seeds,
            color_sets: self.manifest.color_sets,
            seed_locations: self.manifest.seed_locations,
            sketch_bytes: self.manifest.sketch.bytes,
            positions_bytes: self.manifest.positions.bytes,
            packed_bytes: entries.iter().map(|entry| entry.packed.bytes).sum(),
        }
    }

    /// The manifest record of an entry.
    fn entry_record(&self, entry: u32) -> Result<&EntryRecord, IndexError> {
        self.manifest
            .entries
            .get(entry as usize)
            .ok_or_else(|| IndexError::Input(format!("entry {entry} is not in this shard")))
    }
}
