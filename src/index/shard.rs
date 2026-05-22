//! An open shard: its manifest, the sketch in RAM, the positions, and the packed sequences on
//! request.

use super::manifest::{Manifest, PACKED_DIRECTORY, verify_file};
use super::packed::PackedSequence;
use super::positions::{POSITIONS_FILE, Positions};
use super::sketch::{SKETCH_FILE, Sketch};
use super::{IndexError, MappedFile};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Sizes and counts of a shard, for `jam stats`.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
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
    /// Postings of the positions file.
    pub postings: u64,
    /// Bytes of the sketch file.
    pub sketch_bytes: u64,
    /// Bytes of the positions file.
    pub positions_bytes: u64,
    /// Bytes of all packed sequence files.
    pub packed_bytes: u64,
}

/// A shard opened from local files.
pub struct Shard {
    directory: PathBuf,
    packed_directory: PathBuf,
    manifest: Manifest,
    sketch: Sketch,
    positions: Positions,
}

impl Shard {
    /// Opens the shard in `directory`, with its packed sequences in `directory/packed`.
    pub fn open(directory: &Path) -> Result<Self, IndexError> {
        Self::open_with_packed(directory, &directory.join(PACKED_DIRECTORY))
    }

    /// Opens the shard in `directory`, with its packed sequences in `packed_directory`, which
    /// joined shards may share. Only headers and sizes are checked; `verify` reads everything.
    pub fn open_with_packed(directory: &Path, packed_directory: &Path) -> Result<Self, IndexError> {
        let manifest = Manifest::read(directory)?;
        let parameters = manifest.parameters;
        let entries = manifest.entries.len() as u64;
        let sketch = Sketch::open(&directory.join(SKETCH_FILE))?;
        let positions = Positions::open(Arc::new(MappedFile::open(
            &directory.join(POSITIONS_FILE),
            false,
        )?))?;
        let file_size = |name: &str| std::fs::metadata(directory.join(name)).map(|m| m.len());
        if sketch.entries() != entries
            || sketch.bucket_bits() != u32::from(parameters.bucket_bits)
            || sketch.seeds() != manifest.seeds
            || u64::from(positions.entries()) != entries
            || positions.key_bits() != parameters.key_bits()
            || positions.bin_bits() != u32::from(parameters.bin_bits)
            || file_size(SKETCH_FILE)? != manifest.sketch.bytes
            || file_size(POSITIONS_FILE)? != manifest.positions.bytes
        {
            return Err(IndexError::format(
                "the sketch or positions differ from the manifest",
            ));
        }
        Ok(Self {
            directory: directory.to_path_buf(),
            packed_directory: packed_directory.to_path_buf(),
            manifest,
            sketch,
            positions,
        })
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

    /// The local path of an entry's packed sequence.
    pub fn packed_path(&self, entry: u32) -> Result<PathBuf, IndexError> {
        let record = self
            .manifest
            .entries
            .get(entry as usize)
            .ok_or_else(|| IndexError::Input(format!("entry {entry} is not in this shard")))?;
        Ok(self.packed_directory.join(&record.packed.path))
    }

    /// Opens the packed sequence of an entry (stage 3) and checks it against the manifest.
    pub fn packed(&self, entry: u32) -> Result<PackedSequence, IndexError> {
        let file = MappedFile::open(&self.packed_path(entry)?, false)?;
        let record = &self.manifest.entries[entry as usize];
        if file.bytes().len() as u64 != record.packed.bytes {
            return Err(IndexError::format(format!(
                "the packed sequence of {} differs from the manifest",
                record.name
            )));
        }
        let packed = PackedSequence::open(Arc::new(file))?;
        if packed.bases() != record.bases || packed.contig_count() as u64 != record.contigs {
            return Err(IndexError::format(format!(
                "the packed sequence of {} differs from the manifest",
                record.name
            )));
        }
        Ok(packed)
    }

    /// Checks the SHA-256 digest of every file and the checksum of every section.
    pub fn verify(&self) -> Result<(), IndexError> {
        verify_file(&self.directory.join(SKETCH_FILE), &self.manifest.sketch)?;
        verify_file(
            &self.directory.join(POSITIONS_FILE),
            &self.manifest.positions,
        )?;
        self.sketch.verify_checksums()?;
        self.positions.verify_checksums()?;
        for (entry, record) in self.manifest.entries.iter().enumerate() {
            verify_file(&self.packed_path(entry as u32)?, &record.packed)?;
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
            postings: self.manifest.postings,
            sketch_bytes: self.manifest.sketch.bytes,
            positions_bytes: self.manifest.positions.bytes,
            packed_bytes: entries.iter().map(|entry| entry.packed.bytes).sum(),
        }
    }
}
