//! `jam index`: build a shard from FASTA files, or reshard shards (merge and split them).
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{CliError, read_list};
use crate::index::{
    BuildOptions, EntryInput, Manifest, PackedOutput, Parameters, ReshardOptions, Split,
    build_shard, reshard,
};
use clap::{Args, Subcommand};
use std::path::{Path, PathBuf};

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
    /// Most seed locations per output shard, which bounds each shard's sketch RAM
    #[arg(long, group = "split")]
    max_seed_locations: Option<u64>,
    /// Copy the packed sequences instead of linking them
    #[arg(long)]
    copy_packed: bool,
    #[command(flatten)]
    common: CommonArgs,
}

/// Options of both subcommands.
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
    }
}

fn run_build(args: &BuildArgs) -> Result<(), CliError> {
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
