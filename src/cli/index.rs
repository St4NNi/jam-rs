//! `jam index`: build a shard from FASTA files, reshard shards (merge and split them), or warm
//! a cache directory with pinned parts of a shard.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::storage::ShardStorage;
use super::{CliError, read_list};
use crate::index::{
    BuildOptions, EntryInput, Manifest, PackedOutput, Parameters, Part, Place, Plan,
    ReshardOptions, Split, build_shard, record_etags, reshard, warm,
};
use bytesize::ByteSize;
use clap::{Args, Subcommand};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::Instant;

/// The options of `jam index`.
#[derive(Debug, Args)]
pub struct IndexArgs {
    #[command(subcommand)]
    command: IndexCommand,
}

#[derive(Debug, Subcommand)]
enum IndexCommand {
    /// Build one shard from FASTA files
    Build(BuildArgs),
    /// Merge and split shards into an index directory of new shards
    Reshard(ReshardArgs),
    /// Pin parts of a shard in a local cache directory, within RAM and disk budgets
    Warm(WarmArgs),
    /// Record the ETags of an uploaded shard in its local manifest, which binds every S3 read
    /// to them; then upload the manifest again
    Etags(EtagsArgs),
}

/// The options of `jam index etags`.
#[derive(Debug, Args)]
struct EtagsArgs {
    /// The local shard directory whose files were just uploaded
    shard: PathBuf,
    /// Where they were uploaded: s3://bucket/prefix
    uploaded: PathBuf,
}

/// The options of `jam index build`.
#[derive(Debug, Args)]
struct BuildArgs {
    /// The shard directory to create; it must not exist yet
    #[arg(short, long)]
    output: PathBuf,
    /// FASTA files (plain or compressed), one entry per file
    inputs: Vec<PathBuf>,
    /// A file with one FASTA path per line, added to the paths
    #[arg(short, long)]
    list: Option<PathBuf>,
    /// Make each FASTA record an entry instead of each file
    #[arg(long)]
    entry_per_record: bool,
    /// Write only the sketch part and no packed sequences: the shard answers stage 1 only
    #[arg(long)]
    sketch_only: bool,
    #[command(flatten)]
    common: CommonArgs,
    /// Bases per k-mer, 1 to 32 [default: 25]
    #[arg(short, long)]
    k: Option<u8>,
    /// K-mers per minimizer window [default: 100]
    #[arg(short, long)]
    window: Option<u16>,
}

/// The options of `jam index reshard`.
#[derive(Debug, Args)]
#[command(group = clap::ArgGroup::new("split").required(true))]
struct ReshardArgs {
    /// The index directory to create; it must not exist yet
    #[arg(short, long)]
    output: PathBuf,
    /// Shard directories, or index directories of earlier reshards; entries keep this order
    #[arg(required = true)]
    inputs: Vec<PathBuf>,
    /// Number of output shards, balanced by seed locations
    #[arg(long, group = "split")]
    shards: Option<usize>,
    /// Most seed locations per output shard, which bounds the size of each seeds file
    #[arg(long, group = "split")]
    max_seed_locations: Option<u64>,
    /// Copy the packed sequences instead of linking them
    #[arg(long)]
    copy_packed: bool,
    #[command(flatten)]
    common: CommonArgs,
}

/// The options of `jam index warm`.
#[derive(Debug, Args)]
struct WarmArgs {
    /// The shard: a directory or s3://bucket/prefix
    shard: PathBuf,
    /// The cache directory to fill; search and serve read the pins from there
    #[arg(long)]
    cache: PathBuf,
    /// RAM for pins and the block cache, for example 16GiB
    #[arg(long)]
    memory: ByteSize,
    /// Disk for pins and the block cache, for example 200GiB
    #[arg(long)]
    disk: ByteSize,
    /// Parts to pin instead of the default plan, as a comma list of sketch, locators, packed,
    /// or packed=NAME,NAME for some entries: every item after packed= is an entry name, so it
    /// comes last. The flag may repeat
    #[arg(long, value_parser = parse_parts)]
    pin: Vec<Parts>,
}

/// The parts of one `--pin` value.
#[derive(Clone, Debug, PartialEq)]
struct Parts(Vec<Part>);

/// Options of the build and reshard subcommands.
#[derive(Debug, Args)]
struct CommonArgs {
    /// Threads to use [default: all cores]
    #[arg(short, long)]
    threads: Option<usize>,
    /// Directory for temporary files [default: the directory of the output]
    #[arg(long)]
    temp_dir: Option<PathBuf>,
}

impl CommonArgs {
    /// The thread count and the directory for temporary files.
    fn resolve(&self, output: &Path) -> Result<(usize, PathBuf), CliError> {
        let threads = match self.threads {
            Some(0) => return Err(CliError::Usage("--threads must be at least 1".into())),
            Some(threads) => threads,
            None => std::thread::available_parallelism().map_or(1, |n| n.get()),
        };
        let temp_dir = self.temp_dir.clone().unwrap_or_else(|| parent_of(output));
        Ok((threads, temp_dir))
    }
}

/// Runs `jam index build` or `jam index reshard`.
pub fn run(args: &IndexArgs) -> Result<(), CliError> {
    match &args.command {
        IndexCommand::Build(args) => run_build(args),
        IndexCommand::Reshard(args) => run_reshard(args),
        IndexCommand::Warm(args) => run_warm(args),
        IndexCommand::Etags(args) => run_etags(args),
    }
}

fn run_etags(args: &EtagsArgs) -> Result<(), CliError> {
    let storage = ShardStorage::new()?;
    let uploaded = storage.location(&args.uploaded)?;
    let manifest = record_etags(&args.shard, &uploaded).map_err(|source| CliError::Shard {
        path: args.shard.clone(),
        source,
    })?;
    let files = 1 + manifest
        .entries
        .iter()
        .filter(|e| e.packed.is_some())
        .count();
    println!(
        "recorded {files} ETags in {}; upload its manifest.json again",
        args.shard.display()
    );
    Ok(())
}

fn run_warm(args: &WarmArgs) -> Result<(), CliError> {
    let started = Instant::now();
    let storage = ShardStorage::new()?;
    let location = storage.location(&args.shard)?;
    let failed = |source| CliError::Shard {
        path: args.shard.clone(),
        source,
    };
    let (manifest, digest) = location.manifest().map_err(failed)?;
    let (memory, disk) = (args.memory.as_u64(), args.disk.as_u64());
    let parts: Vec<Part> = args.pin.iter().flat_map(|parts| parts.0.clone()).collect();
    let plan = Plan::choose(&manifest, &digest, memory, disk, &parts)?;
    for line in plan_lines(&plan) {
        println!("{line}");
    }
    warm(&location, &manifest, &plan, &args.cache).map_err(failed)?;
    let seconds = started.elapsed().as_secs_f64();
    println!("warmed {} in {seconds:.1} s", args.cache.display());
    Ok(())
}

/// What the plan pins where, one line per file kind, and what is left for the block cache.
fn plan_lines(plan: &Plan) -> Vec<String> {
    let mut lines = vec![format!("shard {}", plan.shard)];
    for place in [Place::Ram, Place::Ssd] {
        let mut sizes: BTreeMap<&str, (usize, u64)> = BTreeMap::new();
        for pin in plan.pins.iter().filter(|pin| pin.place == place) {
            let size = sizes.entry(&pin.part).or_default();
            *size = (size.0 + 1, size.1 + pin.bytes());
        }
        for (part, (files, bytes)) in sizes {
            let place = if place == Place::Ram {
                "in RAM"
            } else {
                "on disk"
            };
            let bytes = ByteSize(bytes).display().iec();
            lines.push(format!("pinned {part} {place}: {bytes} in {files} files"));
        }
    }
    let (ram, disk) = plan.cache_budget();
    let (ram, disk) = (
        ByteSize(ram).display().iec(),
        ByteSize(disk).display().iec(),
    );
    lines.push(format!("block cache: {ram} RAM, {disk} disk"));
    lines
}

/// The parts of a `--pin` value: a comma list of `sketch`, `locators` and `packed`, which may
/// end with `packed=NAME,NAME`; every item after `packed=` is an entry name.
fn parse_parts(text: &str) -> Result<Parts, String> {
    let mut parts = Vec::new();
    let mut items = text.split(',');
    while let Some(item) = items.next() {
        parts.push(match item.split_once('=') {
            None if item == "sketch" => Part::Sketch,
            None if item == "locators" => Part::Locators,
            None if item == "packed" => Part::Packed(None),
            Some(("packed", first)) => {
                let names: Vec<String> = (std::iter::once(first).chain(items.by_ref()))
                    .map(str::to_string)
                    .collect();
                if names.iter().any(String::is_empty) {
                    return Err(format!("{text:?} has an empty entry name"));
                }
                Part::Packed(Some(names))
            }
            _ => {
                return Err(format!(
                    "{item:?} is not sketch, locators, packed or packed=NAME,NAME"
                ));
            }
        });
    }
    Ok(Parts(parts))
}

fn run_build(args: &BuildArgs) -> Result<(), CliError> {
    raise_open_files()?;
    let mut inputs = args.inputs.clone();
    if let Some(list) = &args.list {
        inputs.extend(read_list(list)?);
    }
    if inputs.is_empty() {
        return Err(CliError::Usage(
            "give at least one input path or a --list file".into(),
        ));
    }
    let (threads, temp_dir) = args.common.resolve(&args.output)?;
    let options = build_options(args, threads, temp_dir);
    let manifest = build_shard(&inputs_of(&inputs), &args.output, &options)?;
    println!("{}", summary_line(&args.output, &manifest));
    Ok(())
}

fn run_reshard(args: &ReshardArgs) -> Result<(), CliError> {
    raise_open_files()?;
    let (threads, temp_dir) = args.common.resolve(&args.output)?;
    let split = match (args.shards, args.max_seed_locations) {
        (Some(shards), _) => Split::Shards(shards),
        (None, Some(limit)) => Split::MaxSeedLocations(limit),
        (None, None) => {
            return Err(CliError::Usage(
                "give --shards or --max-seed-locations".into(),
            ));
        }
    };
    let options = ReshardOptions {
        split,
        packed: if args.copy_packed {
            PackedOutput::Copy
        } else {
            PackedOutput::Link
        },
        threads,
        temporary_directory: temp_dir,
    };
    let index = reshard(&args.inputs, &args.output, &options)?;
    println!(
        "wrote {}: {} entries in {} shards",
        args.output.display(),
        index.entries(),
        index.shards.len()
    );
    Ok(())
}

/// Open files a build or reshard needs: its 4,096 spill files, plus the files that its threads
/// read and write at the same time.
const MIN_OPEN_FILES: u64 = 8192;

/// Raises the soft limit of open files to the hard limit, and fails when it stays below
/// `MIN_OPEN_FILES`.
fn raise_open_files() -> Result<(), CliError> {
    let limit = rlimit::increase_nofile_limit(u64::MAX).map_err(CliError::FileLimit)?;
    if limit < MIN_OPEN_FILES {
        return Err(CliError::OpenFiles {
            limit,
            needed: MIN_OPEN_FILES,
        });
    }
    Ok(())
}

fn build_options(args: &BuildArgs, threads: usize, temp_dir: PathBuf) -> BuildOptions {
    let mut parameters = Parameters::default();
    parameters.scheme.k = args.k.unwrap_or(parameters.scheme.k);
    parameters.scheme.w = args.window.unwrap_or(parameters.scheme.w);
    BuildOptions {
        parameters,
        bucket_bits: None,
        threads,
        temporary_directory: temp_dir,
        entry_per_record: args.entry_per_record,
        sketch_only: args.sketch_only,
    }
}

fn inputs_of(paths: &[PathBuf]) -> Vec<EntryInput> {
    paths.iter().cloned().map(EntryInput::from_path).collect()
}

fn parent_of(output: &Path) -> PathBuf {
    match output.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent.to_path_buf(),
        _ => PathBuf::from("."),
    }
}

fn summary_line(output: &Path, manifest: &Manifest) -> String {
    let bases: u64 = manifest.entries.iter().map(|entry| entry.bases).sum();
    format!(
        "wrote {}: {} entries, {bases} bases",
        output.display(),
        manifest.entries.len()
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Command};
    use clap::Parser;

    #[test]
    fn pin_parts_parsed() {
        let parts = |text| parse_parts(text).map(|parts| parts.0);
        assert_eq!(parts("sketch"), Ok(vec![Part::Sketch]));
        assert_eq!(parts("locators"), Ok(vec![Part::Locators]));
        assert_eq!(parts("packed"), Ok(vec![Part::Packed(None)]));
        let all = vec![Part::Sketch, Part::Locators, Part::Packed(None)];
        assert_eq!(parts("sketch,locators,packed"), Ok(all));
        let named = Part::Packed(Some(vec!["a".into(), "b".into()]));
        assert_eq!(parts("packed=a,b"), Ok(vec![named.clone()]));
        assert_eq!(parts("sketch,packed=a,b"), Ok(vec![Part::Sketch, named]));
        // Items after packed= are entry names, even one named like a part.
        let named = Part::Packed(Some(vec!["a".into(), "sketch".into()]));
        assert_eq!(parts("packed=a,sketch"), Ok(vec![named]));
        for bad in ["packed=", "packed=a,,b", "seeds", "sketch=a", "sketch,", ""] {
            assert!(parse_parts(bad).is_err(), "{bad}");
        }
        // The flag repeats, and each value may list several parts.
        let parsed = Cli::try_parse_from([
            "jam",
            "index",
            "warm",
            "s3://b/s",
            "--cache",
            "c",
            "--memory",
            "1GiB",
            "--disk",
            "1GiB",
            "--pin",
            "sketch,locators",
            "--pin",
            "packed=a",
        ]);
        let Ok(Cli {
            command: Command::Index(args),
        }) = parsed
        else {
            panic!("--pin did not parse");
        };
        let IndexCommand::Warm(warm) = args.command else {
            panic!("not warm");
        };
        assert_eq!(warm.pin.len(), 2);
        assert_eq!(warm.pin[0].0, [Part::Sketch, Part::Locators]);
    }
}
