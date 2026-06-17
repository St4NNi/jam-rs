//! `jam index`: build a shard from FASTA files, or join shards.

use super::CliError;
use crate::index::{
    BuildOptions, EntryInput, JoinOptions, Manifest, Parameters, build_shard, join_shards,
};
use clap::Args;
use std::path::{Path, PathBuf};

/// The options of `jam index`.
#[derive(Debug, Args)]
pub struct IndexArgs {
    /// The shard directory to create; it must not exist yet
    #[arg(short, long)]
    output: PathBuf,
    /// FASTA files (plain or compressed), one entry per file; with --join, the shards to join
    inputs: Vec<PathBuf>,
    /// A file with one FASTA path per line, added to the paths
    #[arg(short, long)]
    list: Option<PathBuf>,
    /// Make each FASTA record an entry instead of each file
    #[arg(long)]
    entry_per_record: bool,
    /// Join the given shards into one instead of reading FASTA files
    #[arg(long)]
    join: bool,
    /// Threads to use [default: all cores]
    #[arg(short, long)]
    threads: Option<usize>,
    /// Directory for temporary files [default: the directory of the output]
    #[arg(long)]
    temp_dir: Option<PathBuf>,
    /// Bases per k-mer, 1 to 32 [default: 25]
    #[arg(short, long)]
    k: Option<u8>,
    /// K-mers per minimizer window [default: 100]
    #[arg(short, long)]
    window: Option<u16>,
}

/// Builds or joins the shard.
pub fn run(args: &IndexArgs) -> Result<(), CliError> {
    let mut inputs = args.inputs.clone();
    if let Some(list) = &args.list {
        inputs.extend(read_list(list)?);
    }
    if inputs.is_empty() {
        return Err(CliError::Usage(
            "give at least one input path or a --list file".into(),
        ));
    }
    let threads = match args.threads {
        Some(0) => return Err(CliError::Usage("--threads must be at least 1".into())),
        Some(threads) => threads,
        None => std::thread::available_parallelism().map_or(1, |n| n.get()),
    };
    let temp_dir = args
        .temp_dir
        .clone()
        .unwrap_or_else(|| parent_of(&args.output));
    let manifest = if args.join {
        let options = JoinOptions {
            threads,
            temp_dir,
            bucket_bits: None,
            link_packed: true,
        };
        join_shards(&inputs, &args.output, &options)?
    } else {
        build_shard(
            &inputs_of(&inputs),
            &args.output,
            &build_options(args, threads, temp_dir),
        )?
    };
    println!("{}", summary_line(&args.output, &manifest));
    Ok(())
}

fn build_options(args: &IndexArgs, threads: usize, temp_dir: PathBuf) -> BuildOptions {
    let mut parameters = Parameters::default();
    parameters.k = args.k.unwrap_or(parameters.k);
    parameters.window = args.window.unwrap_or(parameters.window);
    BuildOptions {
        parameters,
        bucket_bits: None,
        threads,
        temp_dir,
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

fn read_list(list: &Path) -> Result<Vec<PathBuf>, CliError> {
    let text = std::fs::read_to_string(list).map_err(|source| CliError::File {
        path: list.to_path_buf(),
        source,
    })?;
    Ok(text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(PathBuf::from)
        .collect())
}

fn summary_line(output: &Path, manifest: &Manifest) -> String {
    let bases: u64 = manifest.entries.iter().map(|entry| entry.bases).sum();
    format!(
        "wrote {}: {} entries, {bases} bases",
        output.display(),
        manifest.entries.len()
    )
}
