//! Builds a shard from FASTA files, one entry per file or per record, with entries indexed in
//! parallel.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

//! The output never depends on the thread count, and entries keep the input order.

use super::manifest::{
    EntryRecord, FORMAT_VERSION, FileRecord, Manifest, PACKED_DIRECTORY, Parameters, check_names,
    choose_bucket_bits, packed_file_name, sha256_file, sha256_hex,
};
use super::packed::PackedWriter;
use super::seed_file::{SEEDS_FILE, SeedFile, SeedLocation, seed_key};
use super::seed_writer::SeedFileShape;
use super::spill::{SeedSpill, seed_record};
use super::{AtPath, IndexError};
use crate::seeds::Seeder;
use crate::storage::LocalFile;
use rayon::prelude::*;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Sequence bytes collected before a parallel round, in one-entry-per-record mode.
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
    pub entry_per_record: bool,
    /// Write only the sketch part of the seeds file and no packed sequences, so the shard
    /// answers stage 1 only.
    pub sketch_only: bool,
}

/// Builds the shard directory `output` from `inputs` and returns its manifest.
pub fn build_shard(
    inputs: &[EntryInput],
    output: &Path,
    options: &BuildOptions,
) -> Result<Manifest, IndexError> {
    options.parameters.validate()?;
    if options.threads == 0 {
        return Err(IndexError::Input("at least one thread is needed".into()));
    }
    if !options.entry_per_record {
        check_names(inputs.iter().map(|input| input.name.as_str()))?;
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.threads)
        .build()?;
    in_partial_directory(output, |directory| {
        let packed = (!options.sketch_only).then(|| directory.join(PACKED_DIRECTORY));
        if let Some(packed) = &packed {
            std::fs::create_dir(packed).at(packed)?;
        }
        let context = Context {
            seeder: Seeder::new(options.parameters.scheme)?,
            parameters: options.parameters,
            packed,
            spill: SeedSpill::create(&options.temporary_directory)?,
        };
        pool.install(|| {
            let built = if options.entry_per_record {
                index_records(inputs, &context)?
            } else {
                index_files(inputs, &context)?
            };
            let (entries, tables): (Vec<_>, Vec<_>) = built.into_iter().unzip();
            check_names(entries.iter().map(|entry| entry.name.as_str()))?;
            let parameters = options.parameters;
            let spill = context.spill;
            let tables = (!options.sketch_only).then_some(tables);
            complete_shard(
                directory,
                parameters,
                entries,
                spill,
                options.bucket_bits,
                tables,
            )
        })
    })
}

/// Runs `build` on a fresh partial directory next to `output`, and moves it to `output` when
/// it succeeds; a failed build leaves nothing behind.
pub(crate) fn in_partial_directory<T>(
    output: &Path,
    build: impl FnOnce(&Path) -> Result<T, IndexError>,
) -> Result<T, IndexError> {
    if output.exists() {
        return Err(IndexError::Input(format!(
            "{} already exists: remove it or choose another output",
            output.display()
        )));
    }
    let name = output
        .file_name()
        .ok_or_else(|| IndexError::Input("the output needs a directory name".into()))?
        .to_string_lossy()
        .into_owned();
    let partial = output.with_file_name(format!("{name}.partial-{}", std::process::id()));
    std::fs::create_dir(&partial).at(&partial)?;
    match build(&partial) {
        Ok(built) => {
            std::fs::rename(&partial, output).at(output)?;
            Ok(built)
        }
        Err(error) => {
            let _ = std::fs::remove_dir_all(&partial);
            Err(error)
        }
    }
}

/// Writes the seeds file from the seed records in `spill` (with the locator part when `locators`)
/// and the manifest of a shard with these entries; the build chooses its bucket bits unless
/// `bucket_bits` is given.
pub(crate) fn complete_shard(
    directory: &Path,
    parameters: Parameters,
    entries: Vec<EntryRecord>,
    spill: SeedSpill,
    bucket_bits: Option<u8>,
    contig_tables: Option<Vec<Vec<u8>>>,
) -> Result<Manifest, IndexError> {
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
    let seeds_path = directory.join(SEEDS_FILE);
    let shape = SeedFileShape {
        bucket_bits: u32::from(parameters.bucket_bits),
        page_bits: u32::from(parameters.page_bits),
        entries: entries.len() as u64,
        locators: contig_tables.is_some(),
    };
    let tables = contig_tables.unwrap_or_default();
    let totals = spill.write(&seeds_path, shape, &tables)?;
    let seed_file = SeedFile::open(Arc::new(LocalFile::open(&seeds_path)?))?;
    let record = |path: &Path, name: &str| -> Result<FileRecord, IndexError> {
        let digest = sha256_file(path)?;
        Ok(FileRecord {
            path: name.to_string(),
            bytes: digest.bytes,
            sha256: digest.sha256,
        })
    };
    let manifest = Manifest {
        format_version: FORMAT_VERSION,
        parameters,
        entries,
        seeds: record(&seeds_path, SEEDS_FILE)?,
        sketch_part: seed_file.sketch_part(),
        locator_part: seed_file.locator_part(),
        stored_seeds: totals.seeds,
        color_sets: totals.color_sets,
        seed_locations: locations,
    };
    manifest.write(directory)?;
    Ok(manifest)
}

/// What every entry job needs.
struct Context {
    seeder: Seeder,
    parameters: Parameters,
    /// The directory of the packed sequences; none when they are not written.
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
                // The name becomes a file name, so it is checked before the file is made.
                let name = String::from_utf8_lossy(name).into_owned();
                check_names(std::iter::once(name.as_str()))?;
                entry.finish(name)
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

    /// Writes the packed sequence of the entry and spills its seed records; returns its record
    /// and its encoded contig table.
    fn finish(mut self, name: String) -> Result<BuiltEntry, IndexError> {
        self.records.sort_unstable();
        self.records.dedup();
        let seed_locations = self.records.len() as u64;
        self.context.spill.push(&mut self.records)?;
        self.records = Vec::new();
        let (contigs, bases) = (self.packed.contig_count(), self.packed.bases());
        let mut table = Vec::new();
        self.packed.contig_table().encode(&mut table);
        let packed = match &self.context.packed {
            Some(directory) => Some(write_packed(self.packed, directory, &name)?),
            None => None,
        };
        let record = EntryRecord {
            name,
            contigs,
            bases,
            seed_locations,
            packed,
        };
        Ok((record, table))
    }
}

/// Writes the packed sequence of the entry `name` into `directory` and returns its record.
fn write_packed(
    writer: PackedWriter,
    directory: &Path,
    name: &str,
) -> Result<FileRecord, IndexError> {
    let parts = writer.finish();
    let file_name = packed_file_name(name);
    let path = directory.join(&file_name);
    let mut file = std::fs::File::create_new(&path).map_err(|error| match error.kind() {
        std::io::ErrorKind::AlreadyExists => {
            IndexError::Input(format!("the entry name {name} is used twice"))
        }
        _ => IndexError::File {
            path: path.clone(),
            source: error,
        },
    })?;
    let (head, packed) = (&parts.header_and_tables, &parts.packed_bases);
    file.write_all(head).at(&path)?;
    file.write_all(packed).at(&path)?;
    Ok(FileRecord {
        path: file_name,
        bytes: (head.len() + packed.len()) as u64,
        sha256: sha256_hex(&[head, packed]),
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
