//! The commands: `index`, `search`, `stats`, `screen` and `serve`.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

mod error;
mod index;
mod rule;
mod search;
mod serve;
mod stats;

use clap::{Parser, Subcommand};
pub use error::CliError;

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
    /// Stream large assemblies past small query sequences (not built yet)
    Screen,
    /// Serve a shard over iroh, as the coordinator of its peers when it has some
    Serve(serve::ServeArgs),
}

/// Parses the command line and runs the chosen command.
pub fn run() -> Result<(), CliError> {
    match Cli::parse().command {
        Command::Index(args) => index::run(&args),
        Command::Search(args) => search::run(&args),
        Command::Stats(args) => stats::run(&args),
        Command::Screen => Err(CliError::Usage("jam screen is not built yet".into())),
        Command::Serve(args) => serve::run(&args),
    }
}
