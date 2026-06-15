//! The commands: `index`, `search`, `stats`, `screen` and `serve`.

mod index;
mod stats;

use anyhow::{Result, bail};
use clap::{Parser, Subcommand};

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
    /// Build a shard from FASTA files, or join shards into one
    Index(index::IndexArgs),
    /// Show what a shard holds
    Stats(stats::StatsArgs),
    /// Stream large assemblies past small query sequences (not built yet)
    Screen,
    /// Serve a shard over TCP (not built yet)
    Serve,
}

/// Parses the command line and runs the chosen command.
pub fn run() -> Result<()> {
    match Cli::parse().command {
        Command::Index(args) => index::run(&args),
        Command::Stats(args) => stats::run(&args),
        Command::Screen => bail!("jam screen is not built yet"),
        Command::Serve => bail!("jam serve is not built yet"),
    }
}
