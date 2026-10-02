//! Builds an index file and its sequence files from FASTA files, one entry per file or per
//! record, with entries indexed in parallel.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

//! The output never depends on the thread count, and entries keep the input order.

use super::manifest::{
    EntryRecord, FORMAT_VERSION, FileRecord, Manifest, Parameters, SequenceRecord, check_names,
    choose_bucket_bits, sha256_hex,
};
use super::packed::PackedWriter;
use super::seed_file::{SeedLocation, seed_key};
use super::seed_writer::{SeedFileShape, WrittenSeeds};
use super::sequence_file::{
    INDEX_EXTENSION, Piece, Place, place_entries, sequence_file_names, write_sequence_files,
};
use super::spill::{SeedSpill, SeedTotals, seed_record};
use super::{AtPath, IndexError};
use crate::seeds::Seeder;
use rayon::prelude::*;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Sequence bytes collected before a parallel round, in singleton mode.
const RECORD_ROUND_BYTES: usize = 256 << 20;

/// One input file and the name of its entry.
#[derive(Clone, Debug)]
pub struct EntryInput {
    /// Entry name.
    pub name: String,
    /// FASTA or FASTQ file, plain or compressed.
    pub path: PathBuf,
}

impl EntryInput {
    /// An input named after its file, without the usual FASTA and compression extensions.
    pub fn from_path(path: PathBuf) -> Self {
        let mut name = path
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_default();
        for extension in [".gz", ".bgz", ".xz", ".bz2", ".zst"] {
            name = name.strip_suffix(extension).unwrap_or(&name).to_string();
        }
        for extension in [".fasta", ".fna", ".fa", ".fastq", ".fq"] {
            name = name.strip_suffix(extension).unwrap_or(&name).to_string();
        }
        Self { name, path }
    }
}

/// How a shard is built.
#[derive(Clone, Debug)]
pub struct BuildOptions {
    /// Parameters of the shard; the build chooses its bucket bits.
    pub parameters: Parameters,
    /// Bucket bits to use instead of the ones chosen from the shard size, so that tests can give a
    /// small shard the buckets of a large one.
    pub bucket_bits: Option<u8>,
    /// Threads that index entries and write the seeds file.
    pub threads: usize,
    /// Directory for temporary files.
    pub temporary_directory: PathBuf,
    /// Make each FASTA record an entry, named by its first word, instead of each file; the names
    /// of the inputs are then not used.
    pub singleton: bool,
    /// Write only the sketch part of the seeds file and no packed sequences, so the shard
    /// answers stage 1 only.
    pub sketch_only: bool,
    /// Sequence files to split the packed sequences into, of about equal size.
    pub chunks: usize,
}

/// The index file of the output `NAME` or `NAME.jam`: `NAME.jam`.
pub fn index_file_path(output: &Path) -> PathBuf {
    let extension = output.extension().and_then(|e| e.to_str());
    match extension {
        Some(INDEX_EXTENSION) => output.to_path_buf(),
        _ => PathBuf::from(format!("{}.{INDEX_EXTENSION}", output.display())),
    }
}

/// Builds the index file `NAME.jam` of the output `NAME` (or `NAME.jam`) and its sequence files
/// from `inputs`, and returns its manifest.
pub fn build_shard(
    inputs: &[EntryInput],
    output: &Path,
    options: &BuildOptions,
) -> Result<Manifest, IndexError> {
    options.parameters.validate()?;
    if options.threads == 0 {
        return Err(IndexError::Input("at least one thread is needed".into()));
    }
    if !options.singleton {
        check_names(inputs.iter().map(|input| input.name.as_str()))?;
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.threads)
        .build()?;
    let temporary = &options.temporary_directory;
    let pieces = (!options.sketch_only)
        .then(|| {
            tempfile::Builder::new()
                .prefix("jam-packed")
                .tempdir_in(temporary)
        })
        .transpose()
        .at(temporary)?;
    let index = index_file_path(output);
    let name = (index.file_stem().and_then(|name| name.to_str()))
        .filter(|name| !name.is_empty())
        .ok_or_else(|| IndexError::Input(format!("{} needs a file name", index.display())))?;
    let chunks = if options.sketch_only {
        0
    } else {
        options.chunks
    };
    let names = sequence_file_names(name, chunks.max(1));
    let mut outputs: Vec<PathBuf> = (names.iter().take(chunks))
        .map(|name| index.with_file_name(name))
        .collect();
    outputs.push(index.clone());
    in_partial_files(&outputs, |partial| {
        let context = Context {
            seeder: Seeder::new(options.parameters.scheme)?,
            parameters: options.parameters,
            packed: pieces.as_ref().map(|pieces| pieces.path().to_path_buf()),
            spill: SeedSpill::create(temporary)?,
        };
        pool.install(|| {
            let built = if options.singleton {
                index_records(inputs, &context)?
            } else {
                index_files(inputs, &context)?
            };
            let (mut entries, tables): (Vec<_>, Vec<_>) = built.into_iter().unzip();
            check_names(entries.iter().map(|entry| entry.name.as_str()))?;
            let (index, chunk_paths) = partial.split_last().unwrap_or((&partial[0], &[]));
            let pieces = context.packed.as_deref().unwrap_or(Path::new(""));
            let (pieces, places) = place_sequences(&mut entries, pieces, chunks)?;
            let tables = (!options.sketch_only).then_some(tables);
            let names = &names[..chunk_paths.len()];
            // The chunks are written while the seeds are, as both take long.
            std::thread::scope(|scope| {
                let files = scope
                    .spawn(|| write_sequence_files(chunk_paths, names, &pieces, &places, true));
                let parameters = options.parameters;
                let bucket_bits = options.bucket_bits;
                let seeds = write_seeds(
                    index,
                    parameters,
                    &entries,
                    context.spill,
                    bucket_bits,
                    tables,
                );
                let files = files
                    .join()
                    .unwrap_or_else(|panic| std::panic::resume_unwind(panic))?;
                finish_index(seeds?, entries, files)
            })
        })
    })
}

/// Gives each entry its place in `chunks` chunks, and returns the packed sequence of each entry,
/// in the file named by its entry id in `pieces`, with the places.
fn place_sequences(
    entries: &mut [EntryRecord],
    pieces: &Path,
    chunks: usize,
) -> Result<(Vec<Piece>, Vec<Place>), IndexError> {
    let pieces: Vec<Piece> = (entries.iter().enumerate())
        .filter_map(|(entry, record)| Some((entry, record.sequence.as_ref()?.bytes)))
        .map(|(entry, bytes)| Piece {
            path: pieces.join(entry.to_string()),
            offset: 0,
            bytes,
        })
        .collect();
    if chunks == 0 {
        return Ok((pieces, Vec::new()));
    }
    let sizes: Vec<u64> = pieces.iter().map(|piece| piece.bytes).collect();
    let places = place_entries(&sizes, chunks)?.0;
    let sequences = entries
        .iter_mut()
        .filter_map(|entry| entry.sequence.as_mut());
    for (sequence, place) in sequences.zip(&places) {
        (sequence.file, sequence.offset) = (place.chunk, place.offset);
    }
    Ok((pieces, places))
}

/// Runs `build` on fresh partial files next to `outputs`, and moves them to `outputs` in order
/// when it succeeds, so the index file, given last, comes last; a failed build leaves nothing
/// behind.
fn in_partial_files<T>(
    outputs: &[PathBuf],
    build: impl FnOnce(&[PathBuf]) -> Result<T, IndexError>,
) -> Result<T, IndexError> {
    if let Some(output) = outputs.iter().find(|output| output.exists()) {
        return Err(IndexError::Input(format!(
            "{} already exists: remove it or choose another output",
            output.display()
        )));
    }
    let partial: Vec<PathBuf> = (outputs.iter())
        .map(|output| {
            PathBuf::from(format!(
                "{}.partial-{}",
                output.display(),
                std::process::id()
            ))
        })
        .collect();
    let built = build(&partial).and_then(|built| {
        for (from, to) in partial.iter().zip(outputs) {
            std::fs::rename(from, to).at(to)?;
        }
        Ok(built)
    });
    if built.is_err() {
        partial
            .iter()
            .for_each(|path| drop(std::fs::remove_file(path)));
    }
    built
}

/// Writes the index file `path` from the seed records in `spill` (with the locator part when there
/// are contig tables) and the manifest of these entries and sequence files; the build chooses
/// the bucket bits unless `bucket_bits` is given.
pub(crate) fn complete_shard(
    path: &Path,
    parameters: Parameters,
    entries: Vec<EntryRecord>,
    sequence_files: Vec<FileRecord>,
    spill: SeedSpill,
    bucket_bits: Option<u8>,
    contig_tables: Option<Vec<Vec<u8>>>,
) -> Result<Manifest, IndexError> {
    let seeds = write_seeds(
        path,
        parameters,
        &entries,
        spill,
        bucket_bits,
        contig_tables,
    )?;
    finish_index(seeds, entries, sequence_files)
}

/// The seeds of an index file, written up to the manifest.
struct Seeds {
    written: WrittenSeeds,
    totals: SeedTotals,
    parameters: Parameters,
    seed_locations: u64,
}

/// Writes the seeds of the index file `path` from the seed records in `spill`, up to the manifest
/// that `finish_index` adds.
fn write_seeds(
    path: &Path,
    parameters: Parameters,
    entries: &[EntryRecord],
    spill: SeedSpill,
    bucket_bits: Option<u8>,
    contig_tables: Option<Vec<Vec<u8>>>,
) -> Result<Seeds, IndexError> {
    if u32::try_from(entries.len()).is_err() {
        return Err(IndexError::Input(
            "a shard holds at most 2^32 - 1 entries".into(),
        ));
    }
    let locations = entries.iter().map(|entry| entry.seed_locations).sum();
    let parameters = Parameters {
        bucket_bits: bucket_bits.unwrap_or_else(|| choose_bucket_bits(locations)),
        ..parameters
    };
    parameters.validate()?;
    let shape = SeedFileShape {
        bucket_bits: u32::from(parameters.bucket_bits),
        page_bits: u32::from(parameters.page_bits),
        entries: entries.len() as u64,
        locators: contig_tables.is_some(),
    };
    let tables = contig_tables.unwrap_or_default();
    let (totals, written) = spill.write(path, shape, &tables)?;
    Ok(Seeds {
        written,
        totals,
        parameters,
        seed_locations: locations,
    })
}

/// Appends the manifest of these entries and sequence files to the seeds of an index file.
fn finish_index(
    seeds: Seeds,
    entries: Vec<EntryRecord>,
    sequence_files: Vec<FileRecord>,
) -> Result<Manifest, IndexError> {
    let written = seeds.written;
    let manifest = Manifest {
        format_version: FORMAT_VERSION,
        parameters: seeds.parameters,
        entries,
        sequence_files,
        seeds_sha256: written.sha256.clone(),
        sketch_part: written.sketch_part.clone(),
        locator_part: written.locator_part.clone(),
        stored_seeds: seeds.totals.seeds,
        color_sets: seeds.totals.color_sets,
        seed_locations: seeds.seed_locations,
    };
    written.append_manifest(&manifest.to_bytes()?)?;
    Ok(manifest)
}

/// What every entry job needs.
struct Context {
    seeder: Seeder,
    parameters: Parameters,
    /// The directory of each entry's packed sequence until the sequence file is written; none
    /// when they are not written.
    packed: Option<PathBuf>,
    spill: SeedSpill,
}

/// An entry's record and its encoded contig table.
type BuiltEntry = (EntryRecord, Vec<u8>);

/// One entry per file: files run in parallel, the largest first.
fn index_files(inputs: &[EntryInput], context: &Context) -> Result<Vec<BuiltEntry>, IndexError> {
    let mut order: Vec<(u64, usize)> = inputs
        .iter()
        .enumerate()
        .map(|(index, input)| Ok((std::fs::metadata(&input.path).at(&input.path)?.len(), index)))
        .collect::<Result<_, IndexError>>()?;
    order.sort_unstable_by(|a, b| b.cmp(a));
    let mut built = order
        .par_iter()
        .with_max_len(1)
        .map(|&(_, index)| {
            let input = &inputs[index];
            let mut entry = EntryBuilder::new(context, index);
            for_each_record(&input.path, |name, sequence| {
                entry.add_contig(name, sequence)
            })?;
            Ok((index, entry.finish(input.name.clone())?))
        })
        .collect::<Result<Vec<_>, IndexError>>()?;
    built.sort_unstable_by_key(|(index, _)| *index);
    Ok(built.into_iter().map(|(_, entry)| entry).collect())
}

/// One entry per record: records are read in order and indexed in parallel rounds.
fn index_records(inputs: &[EntryInput], context: &Context) -> Result<Vec<BuiltEntry>, IndexError> {
    let mut built = Vec::new();
    let mut round: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    let mut round_bytes = 0usize;
    let run_round = |round: &mut Vec<(Vec<u8>, Vec<u8>)>, built: &mut Vec<BuiltEntry>| {
        let first = built.len();
        let entries = round
            .par_iter()
            .enumerate()
            .with_max_len(1)
            .map(|(offset, (name, sequence))| {
                let mut entry = EntryBuilder::new(context, first + offset);
                entry.add_contig(name, sequence)?;
                entry.finish(String::from_utf8_lossy(name).into_owned())
            })
            .collect::<Result<Vec<_>, IndexError>>()?;
        built.extend(entries);
        round.clear();
        Ok::<_, IndexError>(())
    };
    for input in inputs {
        for_each_record(&input.path, |name, sequence| {
            round.push((name.to_vec(), sequence.to_vec()));
            round_bytes += sequence.len();
            if round_bytes >= RECORD_ROUND_BYTES {
                round_bytes = 0;
                run_round(&mut round, &mut built)?;
            }
            Ok(())
        })?;
    }
    run_round(&mut round, &mut built)?;
    Ok(built)
}

/// Collects the contigs of one entry: its packed sequence and its seed records.
struct EntryBuilder<'a> {
    context: &'a Context,
    entry: u32,
    packed: PackedWriter,
    records: Vec<u128>,
}

impl<'a> EntryBuilder<'a> {
    fn new(context: &'a Context, entry: usize) -> Self {
        Self {
            context,
            entry: entry as u32,
            packed: PackedWriter::new(),
            records: Vec::new(),
        }
    }

    fn add_contig(&mut self, name: &[u8], sequence: &[u8]) -> Result<(), IndexError> {
        let start = self.packed.bases();
        self.packed.add_contig(name, sequence)?;
        let bin_bits = self.context.parameters.bin_bits;
        if self.packed.bases() >> bin_bits >= 1 << 31 {
            return Err(IndexError::Input(
                "an entry has too many bases for its bins".into(),
            ));
        }
        let (records, entry) = (&mut self.records, self.entry);
        self.context.seeder.seeds(sequence, |seed| {
            let bin = ((start + seed.position) >> bin_bits) as u32;
            let location = SeedLocation {
                bin,
                reverse: seed.reverse,
            };
            records.push(seed_record(seed_key(seed.hash), entry, location.value()));
        });
        Ok(())
    }

    /// Writes the packed sequence of the entry into a file named by its entry id and spills its
    /// seed records; returns its record and its encoded contig table.
    fn finish(mut self, name: String) -> Result<BuiltEntry, IndexError> {
        self.records.sort_unstable();
        self.records.dedup();
        let seed_locations = self.records.len() as u64;
        self.context.spill.push(&mut self.records)?;
        self.records = Vec::new();
        let (contigs, bases) = (self.packed.contig_count(), self.packed.bases());
        let mut table = Vec::new();
        self.packed.contig_table().encode(&mut table);
        let sequence = match &self.context.packed {
            Some(pieces) => {
                let path = pieces.join(self.entry.to_string());
                Some(write_packed(self.packed, &path)?)
            }
            None => None,
        };
        let record = EntryRecord {
            name,
            contigs,
            bases,
            seed_locations,
            sequence,
        };
        Ok((record, table))
    }
}

/// Writes the packed sequence of an entry to `path`; its record has its size and digest, and
/// gets its place when the chunks are written.
fn write_packed(writer: PackedWriter, path: &Path) -> Result<SequenceRecord, IndexError> {
    let parts = writer.finish()?;
    let file = std::fs::File::create_new(path).at(path)?;
    let parts = parts.parts();
    let mut out = std::io::BufWriter::with_capacity(1 << 20, file);
    for part in &parts {
        out.write_all(part).at(path)?;
    }
    out.flush().at(path)?;
    Ok(SequenceRecord {
        file: 0,
        offset: 0,
        bytes: parts.iter().map(|part| part.len() as u64).sum(),
        sha256: sha256_hex(&parts),
    })
}

/// Calls `visit(name, sequence)` for every record of a FASTA or FASTQ file.
fn for_each_record(
    path: &Path,
    mut visit: impl FnMut(&[u8], &[u8]) -> Result<(), IndexError>,
) -> Result<(), IndexError> {
    let failed =
        |error: &dyn std::fmt::Display| IndexError::Input(format!("{}: {error}", path.display()));
    let mut reader = match needletail::parse_fastx_file(path) {
        Ok(reader) => reader,
        Err(error) if error.kind == needletail::errors::ParseErrorKind::EmptyFile => {
            return Ok(());
        }
        Err(error) => return Err(failed(&error)),
    };
    while let Some(record) = reader.next() {
        let record = record.map_err(|error| failed(&error))?;
        let name = record.id().split(u8::is_ascii_whitespace).next();
        visit(name.unwrap_or_default(), &record.seq())?;
    }
    Ok(())
}
