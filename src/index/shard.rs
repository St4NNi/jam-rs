//! An open shard: its manifest, the seeds file, and the packed sequences on request, read from
//! a local directory or an S3 prefix.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::IndexError;
use super::manifest::{EntryRecord, FileRecord, Manifest, verify_source};
use super::packed::{PackedSequence, contig_table_bytes, read_contig_starts};
use super::seed_file::SeedFile;
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
    /// Seed locations of all entries.
    pub seed_locations: u64,
    /// Bytes of the sketch part of the seeds file.
    pub sketch_bytes: u64,
    /// Bytes of the locator part of the seeds file, 0 without it.
    pub locator_bytes: u64,
    /// Bytes of all packed sequence files.
    pub packed_bytes: u64,
}

/// A shard opened from a local directory or an S3 prefix.
pub struct Shard {
    files: ShardFiles,
    manifest: Manifest,
    /// The SHA-256 digest of the manifest file, lower case hex.
    digest: String,
    seeds: SeedFile,
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
        let (manifest, digest) = location.manifest()?;
        let files = ShardFiles::new(location, options, manifest.entries.len(), &digest)?;
        let parameters = manifest.parameters;
        let entries = manifest.entries.len() as u64;
        let seeds_file = files.open(&manifest.seeds)?;
        let seeds_bytes = seeds_file.size();
        let seeds = SeedFile::open(seeds_file)?;
        if seeds.entry_count() != entries
            || seeds.bucket_bits() != u32::from(parameters.bucket_bits)
            || seeds.seed_count() != manifest.stored_seeds
            || seeds_bytes != manifest.seeds.bytes
            || seeds.sketch_part() != manifest.sketch_part
            || seeds.locator_part() != manifest.locator_part
        {
            return Err(IndexError::format(
                "the seeds file differs from the manifest",
            ));
        }
        Ok(Self {
            files,
            manifest,
            digest,
            seeds,
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

    /// The seeds file (stages 1 and 2).
    pub fn seed_file(&self) -> &SeedFile {
        &self.seeds
    }

    /// The local path of an entry's packed sequence; an error when it lies on S3.
    pub fn packed_path(&self, entry: u32) -> Result<PathBuf, IndexError> {
        let (record, file) = self.packed_record(entry)?;
        self.files.packed_path(file).ok_or_else(|| {
            IndexError::Input(format!("the packed sequence of {} is on S3", record.name))
        })
    }

    /// The start of each contig of an entry in entry coordinates, followed by the entry's bases,
    /// so contig `i` spans `starts[i]..starts[i + 1]`. One read at the front of its packed file.
    pub fn contig_starts(&self, entry: u32) -> Result<Vec<u64>, IndexError> {
        let (record, file) = self.packed_record(entry)?;
        let source = self.files.open_packed(entry as usize, file)?;
        read_contig_starts(source.as_ref(), record.contigs)
    }

    /// Tells the source of an entry's packed file that its contig starts will be read soon.
    pub fn prefetch_contig_starts(&self, entry: u32) {
        let Ok((record, file)) = self.packed_record(entry) else {
            return;
        };
        if let Ok(source) = self.files.open_packed(entry as usize, file) {
            source.will_need(0, contig_table_bytes(record.contigs));
        }
    }

    /// Opens the packed sequence of an entry (stage 3) and checks it against the manifest.
    pub fn packed(&self, entry: u32) -> Result<PackedSequence, IndexError> {
        let (record, file) = self.packed_record(entry)?;
        let source = self.files.open_packed(entry as usize, file)?;
        let size = source.size();
        let packed = PackedSequence::open(source)?;
        if size != file.bytes
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
        let seeds = &self.manifest.seeds;
        verify_source(self.files.open_uncached(seeds)?.as_ref(), seeds)?;
        self.seeds.verify_checksums()?;
        for file in self
            .manifest
            .entries
            .iter()
            .filter_map(|e| e.packed.as_ref())
        {
            let source = self.files.open_packed_uncached(file)?;
            verify_source(source.as_ref(), file)?;
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
            seeds: self.manifest.stored_seeds,
            color_sets: self.manifest.color_sets,
            seed_locations: self.manifest.seed_locations,
            sketch_bytes: part_bytes(Some(&self.manifest.sketch_part)),
            locator_bytes: part_bytes(self.manifest.locator_part.as_ref()),
            packed_bytes: (entries.iter().filter_map(|entry| entry.packed.as_ref()))
                .map(|file| file.bytes)
                .sum(),
        }
    }

    /// The manifest record of an entry and of its packed file; an error for a shard built with
    /// the sketch part only.
    fn packed_record(&self, entry: u32) -> Result<(&EntryRecord, &FileRecord), IndexError> {
        let record = (self.manifest.entries.get(entry as usize))
            .ok_or_else(|| IndexError::Input(format!("entry {entry} is not in this shard")))?;
        let file = record.packed.as_ref().ok_or_else(|| {
            IndexError::Input("the shard was built with the sketch part only".into())
        })?;
        Ok((record, file))
    }
}

/// Bytes of a part of the seeds file, 0 without it.
fn part_bytes(part: Option<&std::ops::Range<u64>>) -> u64 {
    part.map_or(0, |part| part.end - part.start)
}
