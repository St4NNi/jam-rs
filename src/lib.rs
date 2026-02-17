pub mod bias;
pub mod cli;
pub mod core_utils;
pub mod format;
pub mod io;
pub mod query;
pub mod reader;
pub mod sketch;
pub mod writer;
pub use cli::handlers::{
    handle_bias_create_command, handle_bias_stats_command, handle_distance_command,
    handle_sketch_command, handle_stats_command,
};
pub use io::{expand_input_paths, is_sequence_file};
pub use jamhash::jamhash_u64;

use anyhow::Result;
use clap::Parser;
use cli::{BiasCommands, Cli, Commands};

pub fn run() -> Result<()> {
    let cli = Cli::parse();

    if let Some(threads) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .stack_size(8 * 1024 * 1024)
            .build_global()?;
    }

    match cli.command {
        Commands::Sketch {
            input,
            output,
            kmer_size,
            fscale,
            complexity,
            singleton,
            temp_dir,
            bias_table,
        } => {
            let expanded_inputs = expand_input_paths(&input)?;
            handle_sketch_command(
                expanded_inputs,
                output,
                kmer_size,
                fscale,
                singleton,
                cli.threads.unwrap_or(1),
                cli.memory.unwrap_or(2),
                cli.force,
                cli.silent,
                complexity,
                temp_dir,
                bias_table,
            )
        }

        Commands::Bias { command } => match command {
            BiasCommands::Create {
                positive,
                negative,
                output,
                kmer_size,
                fscale,
                cms_width,
                cms_depth,
                alpha,
                positive_retention,
                negative_retention,
                positive_fscale,
                negative_fscale,
                unbiased_fscale,
                threads,
            } => handle_bias_create_command(
                positive,
                negative,
                output,
                kmer_size,
                fscale,
                cms_width,
                cms_depth,
                alpha,
                positive_retention,
                negative_retention,
                positive_fscale,
                negative_fscale,
                unbiased_fscale,
                threads.or(cli.threads),
                cli.force,
                cli.silent,
            ),
            BiasCommands::Stats { input, output } => {
                handle_bias_stats_command(input, output, cli.silent)
            }
        },

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
