//! The commands: `index`, `search`, `stats`, `screen` and `serve`.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

mod error;
mod index;
mod rule;
mod screen;
mod search;
mod serve;
mod stats;

use clap::{Parser, Subcommand};
pub use error::CliError;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// jam (just another minimizer) finds small sequences in very large sets of assemblies and aligns
/// them with CIGAR strings.
#[derive(Debug, Parser)]
#[command(name = "jam", version, about, arg_required_else_help = true)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Build a shard from FASTA files, or merge and split shards
    Index(index::IndexArgs),
    /// Search query sequences in a shard and align them
    Search(search::SearchArgs),
    /// Show what a shard holds
    Stats(stats::StatsArgs),
    /// Stream large assemblies past a shard of small sequences and report what each holds
    Screen(screen::ScreenArgs),
    /// Serve a shard over iroh, as the coordinator of its peers when it has some
    Serve(serve::ServeArgs),
}

/// Parses the command line and runs the chosen command.
pub fn run() -> Result<(), CliError> {
    match Cli::parse().command {
        Command::Index(args) => index::run(&args),
        Command::Search(args) => search::run(&args),
        Command::Stats(args) => stats::run(&args),
        Command::Screen(args) => screen::run(&args),
        Command::Serve(args) => serve::run(&args),
    }
}

/// The output of a command: the file at `path`, buffered, or the standard output.
fn output_writer(path: Option<&Path>) -> Result<Box<dyn Write>, CliError> {
    let Some(path) = path else {
        return Ok(Box::new(std::io::stdout()));
    };
    let file = File::create(path).map_err(|source| CliError::File {
        path: path.to_path_buf(),
        source,
    })?;
    Ok(Box::new(BufWriter::new(file)))
}

/// The paths of a list file: one per line, blank lines left out.
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
