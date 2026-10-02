//! An open shard: its manifest, the seeds file, and the packed sequences on request, read from
//! a local directory or an S3 prefix.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::IndexError;
use super::contig_table::ContigTable;
use super::manifest::{EntryRecord, Manifest, SequenceRecord, sha256_range, verify_source};
use super::packed::PackedSequence;
use super::seed_file::{HEADER_BYTES, SEEDS_FILE, SeedFile};
use super::sources::{IndexHead, Location, ShardFiles, ShardOptions};
use crate::storage::SliceSource;
use serde::Serialize;
use std::path::{Path, PathBuf};
use std::sync::Arc;
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
    /// Bytes of the packed sequences of all entries.
    pub packed_bytes: u64,
}

/// A shard opened from a local directory or an S3 prefix.
pub struct Shard {
    files: ShardFiles,
    head: IndexHead,
    seeds: SeedFile,
    /// Wall seconds that opening took.
    open_seconds: f64,
}

impl Shard {
    /// Whether an S3 object was found replaced since the shard was opened. Its reads then fail
    /// until the shard is opened again, with a new manifest.
    pub fn is_stale(&self) -> bool {
        self.files.is_stale()
    }

    /// Opens the shard in `directory`, with its packed sequences in `directory/packed`.
    pub fn open(directory: &Path) -> Result<Self, IndexError> {
        let location = Location::File(directory.join(SEEDS_FILE));
        Self::open_with(&location, &ShardOptions::default())
    }

    /// Opens the shard at `location` as `options` say. Only headers and sizes are checked;
    /// `verify` reads everything.
    pub fn open_with(location: &Location, options: &ShardOptions) -> Result<Self, IndexError> {
        Self::open_with_head(location, options, location.head()?)
    }

    /// Opens the shard at `location` as `options` say, with the head of its seeds file already
    /// read.
    pub fn open_with_head(
        location: &Location,
        options: &ShardOptions,
        head: IndexHead,
    ) -> Result<Self, IndexError> {
        let started = Instant::now();
        let files = ShardFiles::new(location, options, &head)?;
        let manifest = &head.manifest;
        let parameters = manifest.parameters;
        let entries = manifest.entries.len() as u64;
        let seeds = SeedFile::open(files.open_seeds()?)?;
        if seeds.entry_count() != entries
            || seeds.bucket_bits() != u32::from(parameters.bucket_bits)
            || seeds.seed_count() != manifest.stored_seeds
            || seeds.sketch_part() != manifest.sketch_part
            || seeds.locator_part() != manifest.locator_part
            || seeds.manifest_range().end != head.bytes
        {
            return Err(IndexError::format(
                "the seeds file differs from the manifest",
            ));
        }
        Ok(Self {
            files,
            head,
            seeds,
            open_seconds: started.elapsed().as_secs_f64(),
        })
    }

    /// The SHA-256 digest of the manifest, which names this version of the shard.
    pub fn digest(&self) -> &str {
        &self.head.digest
    }

    /// The head of the seeds file: manifest, digest, size and name.
    pub fn head(&self) -> &IndexHead {
        &self.head
    }

    /// Wall seconds that opening the shard took.
    pub fn open_seconds(&self) -> f64 {
        self.open_seconds
    }

    /// The manifest of the shard.
    pub fn manifest(&self) -> &Manifest {
        &self.head.manifest
    }

    /// The seeds file (stages 1 and 2).
    pub fn seed_file(&self) -> &SeedFile {
        &self.seeds
    }

    /// The local path of the sequence file of an entry; an error when it lies on S3.
    pub fn sequence_path(&self, entry: u32) -> Result<PathBuf, IndexError> {
        let (record, sequence) = self.packed_record(entry)?;
        self.files.sequence_path(sequence.file).ok_or_else(|| {
            IndexError::Input(format!("the packed sequence of {} is on S3", record.name))
        })
    }

    /// The contig table of an entry from the locator part of the seeds file, checked against
    /// the manifest.
    pub fn contig_table(&self, entry: u32) -> Result<ContigTable, IndexError> {
        let (record, _) = self.packed_record(entry)?;
        let table = self.seeds.contig_table(entry)?;
        if table.contig_count() as u64 != record.contigs || table.bases() != record.bases {
            return Err(IndexError::format(format!(
                "the contig table of {} differs from the manifest",
                record.name
            )));
        }
        Ok(table)
    }

    /// Opens the packed sequence of an entry (stage 3) with its contig table. Its header is
    /// checked by the first read.
    pub fn packed(&self, entry: u32) -> Result<PackedSequence, IndexError> {
        let table = self.contig_table(entry)?;
        let (_, sequence) = self.packed_record(entry)?;
        let file = self.files.sequence(sequence.file)?;
        let source = SliceSource::new(file, sequence.range())?;
        Ok(PackedSequence::open(Arc::new(source), table))
    }

    /// Checks the SHA-256 digest of the seeds and of every sequence file, and the checksums of
    /// every section and of the manifest.
    pub fn verify(&self) -> Result<(), IndexError> {
        let file = self.files.open_seeds_uncached()?;
        let seeds = HEADER_BYTES..self.seeds.manifest_range().start;
        if sha256_range(file.as_ref(), seeds)? != self.manifest().seeds_sha256 {
            return Err(IndexError::Checksum(
                "the seeds differ from the manifest".into(),
            ));
        }
        self.seeds.verify_checksums()?;
        for file in 0..self.manifest().sequence_files.len() {
            let (source, record) = self.files.open_uncached(file)?;
            verify_source(source.as_ref(), record)?;
        }
        Ok(())
    }

    /// Sizes and counts of the shard.
    pub fn stats(&self) -> ShardStats {
        let entries = &self.manifest().entries;
        ShardStats {
            entries: entries.len() as u64,
            contigs: entries.iter().map(|entry| entry.contigs).sum(),
            bases: entries.iter().map(|entry| entry.bases).sum(),
            seeds: self.manifest().stored_seeds,
            color_sets: self.manifest().color_sets,
            seed_locations: self.manifest().seed_locations,
            sketch_bytes: part_bytes(Some(&self.manifest().sketch_part)),
            locator_bytes: part_bytes(self.manifest().locator_part.as_ref()),
            packed_bytes: (entries.iter().filter_map(|entry| entry.sequence.as_ref()))
                .map(|sequence| sequence.bytes)
                .sum(),
        }
    }

    /// The manifest record of an entry and where its packed sequence lies; an error for a shard
    /// built with the sketch part only.
    fn packed_record(&self, entry: u32) -> Result<(&EntryRecord, &SequenceRecord), IndexError> {
        let record = (self.manifest().entries.get(entry as usize))
            .ok_or_else(|| IndexError::Input(format!("entry {entry} is not in this shard")))?;
        let sequence = record.sequence.as_ref().ok_or_else(|| {
            IndexError::Input("the shard was built with the sketch part only".into())
        })?;
        Ok((record, sequence))
    }
}

/// Bytes of a part of the seeds file, 0 without it.
fn part_bytes(part: Option<&std::ops::Range<u64>>) -> u64 {
    part.map_or(0, |part| part.end - part.start)
}
