pub mod cli;
pub mod core;
pub mod db;
pub mod distance;
pub mod io;
pub mod sketch;
pub mod stats;
pub mod writer;

pub use core as core_utils;
pub use io::{expand_input_paths, is_sequence_file};
pub use cli::handlers::{handle_sketch_command, handle_distance_command, handle_stats_command};
pub use jamhash::jamhash_u64;

use anyhow::Result;
use clap::Parser;
use cli::{Cli, Commands};

pub fn run() -> Result<()> {
    let cli = Cli::parse();

    if let Some(threads) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()?;
    }

    match cli.command {
        Commands::Sketch {
            input,
            output,
            kmer_size,
            fscale,
            nmax,
            complexity,
            singleton,
            temp_dir,
        } => {
            let expanded_inputs = expand_input_paths(&input)?;
            handle_sketch_command(
                expanded_inputs,
                output,
                kmer_size,
                fscale,
                nmax,
                singleton,
                cli.threads.unwrap_or(1),
                cli.memory.unwrap_or(2),
                cli.force,
                cli.silent,
                complexity,
                temp_dir,
            )
        }

        Commands::Dist {
            input,
            database,
            output,
            cutoff,
            singleton,
        } => handle_distance_command(input, database, output, cutoff, singleton, cli.silent),

        Commands::Stats { input, short, full } => {
            handle_stats_command(input, short, full, cli.silent)
        }
    }
}
