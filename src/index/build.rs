//! Builds a shard from FASTA files: one entry per file, or one entry per record. Entries are
//! indexed in parallel; each one writes its packed sequence and its positions section at once.
//! The positions file then takes the sections in entry order, and the sketch is built from it.
//! The output never depends on the thread count, and entries keep the input order.

use super::manifest::{
    EntryRecord, FORMAT_VERSION, FileRecord, Manifest, PACKED_DIRECTORY, Parameters, check_names,
    packed_file_name, sha256_file, sha256_hex,
};
use super::packed::PackedWriter;
use super::positions::{POSITIONS_FILE, Positions, PositionsWriter, Section, encode_section};
use super::sketch::SKETCH_FILE;
use super::sketch::build::build_sketch;
use super::{IndexError, MappedFile, seed_key};
use crate::seeds::Seeder;
use rayon::prelude::*;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Sequence bytes collected before a parallel round, in one-entry-per-record mode.
const RECORD_ROUND_BYTES: usize = 256 << 20;

/// One input file and the name of its entry.
#[derive(Clone, Debug)]
pub struct EntryInput {
    /// Entry name (unused when each record is an entry).
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
    /// Parameters of the shard.
    pub parameters: Parameters,
    /// Threads that index entries and build the sketch.
    pub threads: usize,
    /// Directory for temporary files.
    pub temp_dir: PathBuf,
    /// Make each FASTA record an entry, named by its first word, instead of each file.
    pub entry_per_record: bool,
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
        .build()
        .map_err(|error| IndexError::Input(format!("thread pool: {error}")))?;
    in_partial_directory(output, |directory| {
        std::fs::create_dir(directory.join(PACKED_DIRECTORY))?;
        let temp = tempfile::Builder::new()
            .prefix("jam-build")
            .tempdir_in(&options.temp_dir)?;
        let context = Context {
            seeder: Seeder::new(options.parameters.scheme())?,
            parameters: options.parameters,
            packed: directory.join(PACKED_DIRECTORY),
            temp: temp.path().to_path_buf(),
        };
        pool.install(|| {
            let built = if options.entry_per_record {
                index_records(inputs, &context)?
            } else {
                index_files(inputs, &context)?
            };
            check_names(built.iter().map(|entry| entry.record.name.as_str()))?;
            let entries = built.iter().map(|entry| entry.record.clone()).collect();
            complete_shard(
                directory,
                options.parameters,
                entries,
                &options.temp_dir,
                |writer| {
                    for entry in &built {
                        let bytes = std::fs::read(&entry.section_path)?;
                        writer.push_section(&entry.section, &bytes)?;
                        std::fs::remove_file(&entry.section_path)?;
                    }
                    Ok(())
                },
            )
        })
    })
}

/// Runs `build` on a fresh partial directory next to `output`, and moves it to `output` when
/// it succeeds; a failed build leaves nothing behind.
pub(crate) fn in_partial_directory(
    output: &Path,
    build: impl FnOnce(&Path) -> Result<Manifest, IndexError>,
) -> Result<Manifest, IndexError> {
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
    std::fs::create_dir(&partial)?;
    match build(&partial) {
        Ok(manifest) => {
            std::fs::rename(&partial, output)?;
            Ok(manifest)
        }
        Err(error) => {
            let _ = std::fs::remove_dir_all(&partial);
            Err(error)
        }
    }
}

/// Writes the positions (sections given by `write_sections`, in entry order), builds the sketch
/// from them, and writes the manifest.
pub(crate) fn complete_shard(
    directory: &Path,
    parameters: Parameters,
    entries: Vec<EntryRecord>,
    temp_dir: &Path,
    write_sections: impl FnOnce(&mut PositionsWriter) -> Result<(), IndexError>,
) -> Result<Manifest, IndexError> {
    if u32::try_from(entries.len()).is_err() {
        return Err(IndexError::Input(
            "a shard holds fewer than 2^32 entries".into(),
        ));
    }
    let positions_path = directory.join(POSITIONS_FILE);
    let key_bits = parameters.key_bits();
    let mut writer = PositionsWriter::create(
        &positions_path,
        key_bits,
        u32::from(parameters.bin_bits),
        entries.len(),
    )?;
    write_sections(&mut writer)?;
    writer.finish()?;
    let positions = Positions::open(Arc::new(MappedFile::open(&positions_path, false)?))?;
    let sketch_path = directory.join(SKETCH_FILE);
    let totals = build_sketch(
        &positions,
        u32::from(parameters.bucket_bits),
        &sketch_path,
        temp_dir,
    )?;
    let postings = (0..positions.entries())
        .map(|entry| positions.entry_postings(entry))
        .sum::<Result<u64, IndexError>>()?;
    let record = |path: &Path, name: &str| -> Result<FileRecord, IndexError> {
        let (bytes, sha256) = sha256_file(path)?;
        Ok(FileRecord {
            path: name.to_string(),
            bytes,
            sha256,
        })
    };
    let manifest = Manifest {
        format_version: FORMAT_VERSION,
        parameters,
        entries,
        sketch: record(&sketch_path, SKETCH_FILE)?,
        positions: record(&positions_path, POSITIONS_FILE)?,
        seeds: totals.seeds,
        color_sets: totals.color_sets,
        postings,
    };
    manifest.write(directory)?;
    Ok(manifest)
}

/// What every entry job needs.
struct Context {
    seeder: Seeder,
    parameters: Parameters,
    packed: PathBuf,
    temp: PathBuf,
}

/// An indexed entry whose positions section waits in a temporary file.
struct BuiltEntry {
    record: EntryRecord,
    section: Section,
    section_path: PathBuf,
}

/// Collects the contigs of one entry: its packed sequence and its seeds as (key, posting).
struct EntryBuilder<'a> {
    context: &'a Context,
    packed: PackedWriter,
    seeds: Vec<(u64, u32)>,
}

impl<'a> EntryBuilder<'a> {
    fn new(context: &'a Context) -> Self {
        Self {
            context,
            packed: PackedWriter::new(),
            seeds: Vec::new(),
        }
    }

    fn add_contig(&mut self, name: &[u8], sequence: &[u8]) -> Result<(), IndexError> {
        let start = self.packed.bases();
        self.packed.add_contig(name, sequence)?;
        let bin_bits = self.context.parameters.bin_bits;
        let key_bits = self.context.parameters.key_bits();
        if self.packed.bases() >> bin_bits >= 1 << 31 {
            return Err(IndexError::Input(
                "an entry has too many bases for its bins".into(),
            ));
        }
        let seeds = &mut self.seeds;
        self.context.seeder.seeds(sequence, |seed| {
            let bin = ((start + seed.position) >> bin_bits) as u32;
            seeds.push((
                seed_key(seed.hash, key_bits),
                (bin << 1) | u32::from(seed.reverse),
            ));
        });
        Ok(())
    }

    /// Writes the packed sequence and the section of entry `index`.
    fn finish(mut self, index: usize, name: String) -> Result<BuiltEntry, IndexError> {
        self.seeds.sort_unstable();
        self.seeds.dedup();
        let (section, bytes) = encode_section(&self.seeds, self.context.parameters.key_bits())?;
        drop(self.seeds);
        let section_path = self.context.temp.join(format!("section-{index:010}"));
        std::fs::write(&section_path, &bytes)?;
        let (contigs, bases) = (self.packed.contigs(), self.packed.bases());
        let (head, packed) = self.packed.finish();
        let file_name = packed_file_name(&name);
        let path = self.context.packed.join(&file_name);
        let mut file = std::fs::File::create_new(path).map_err(|error| match error.kind() {
            std::io::ErrorKind::AlreadyExists => {
                IndexError::Input(format!("the entry name {name} is used twice"))
            }
            _ => error.into(),
        })?;
        file.write_all(&head)?;
        file.write_all(&packed)?;
        let record = EntryRecord {
            name,
            contigs,
            bases,
            packed: FileRecord {
                path: file_name,
                bytes: (head.len() + packed.len()) as u64,
                sha256: sha256_hex(&[&head, &packed]),
            },
        };
        Ok(BuiltEntry {
            record,
            section,
            section_path,
        })
    }
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

/// One entry per file: files run in parallel, the largest first.
fn index_files(inputs: &[EntryInput], context: &Context) -> Result<Vec<BuiltEntry>, IndexError> {
    let mut order: Vec<(u64, usize)> = inputs
        .iter()
        .enumerate()
        .map(|(index, input)| Ok((std::fs::metadata(&input.path)?.len(), index)))
        .collect::<Result<_, IndexError>>()?;
    order.sort_unstable_by(|a, b| b.cmp(a));
    let mut built = order
        .par_iter()
        .with_max_len(1)
        .map(|&(_, index)| {
            let input = &inputs[index];
            let mut entry = EntryBuilder::new(context);
            for_each_record(&input.path, |name, sequence| {
                entry.add_contig(name, sequence)
            })?;
            Ok((index, entry.finish(index, input.name.clone())?))
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
                let mut entry = EntryBuilder::new(context);
                entry.add_contig(name, sequence)?;
                let name = String::from_utf8_lossy(name).into_owned();
                check_names(std::iter::once(name.as_str()))?;
                entry.finish(first + offset, name)
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
