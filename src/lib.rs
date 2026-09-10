pub mod alignment;
pub mod bgzf;
pub mod bgzf_cache;
pub mod bias;
pub mod cli;
pub mod collection;
pub mod core_utils;
pub mod format;
pub mod io;
pub mod jidx;
pub mod jidx_builder;
mod jidx_filters;
mod jidx_postings;
pub mod jidx_reader;
mod jidx_runs;
pub mod jidx_writer;
pub mod mosaic;
#[cfg(test)]
mod owner_export_tests;
mod owner_file;
mod owner_format;
#[cfg(test)]
mod owner_integrity_tests;
pub mod owner_observer;
pub mod owner_postings;
mod owner_reader;
#[cfg(test)]
mod owner_tests;
pub mod owner_writer;
pub mod query;
pub mod range_source;
pub mod reader;
mod shared_file;
pub mod shared_format;
pub mod shared_seed;
pub mod shared_writer;
pub mod sketch;
pub mod trace;
mod trace_index;
pub mod writer;
pub use cli::handlers::{
    handle_bias_create_command, handle_bias_stats_command, handle_distance_command,
    handle_sketch_command, handle_stats_command,
};
pub use io::{expand_input_paths, is_sequence_file};
pub use jamhash::jamhash_u64;

use anyhow::Result;
use clap::Parser;
use cli::handlers::{TraceArgs, TraceInput, handle_jidx_build_command, handle_trace_command};
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
                target_fscale,
                max_fscale,
                unseen_fscale,
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
                target_fscale,
                max_fscale,
                unseen_fscale,
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
        } => handle_distance_command(
            input,
            database,
            output,
            cutoff,
            singleton,
            cli.force,
            cli.silent,
            cli.memory.unwrap_or(2),
        ),

        Commands::Jidx {
            database,
            manifest,
            output,
            kmer_size,
            minimizer_window,
            rescue_k15,
        } => handle_jidx_build_command(
            database,
            manifest,
            output,
            jidx_builder::JidxBuildConfig {
                k: kmer_size,
                minimizer_window,
                rescue_k15,
            },
            cli.force,
        ),

        Commands::SharedIndex {
            reference_index,
            output,
            minimizer_window,
        } => {
            let stats =
                shared_writer::build_shared_index(reference_index, output, minimizer_window)?;
            if !cli.silent {
                println!("{}", serde_json::to_string(&stats)?);
            }
            Ok(())
        }

        Commands::Trace {
            query,
            database,
            index,
            manifest,
            collection,
            owner_index,
            audit_index,
            output,
            query_id,
            linear,
            no_sketch,
            min_containment,
            max_metagenomes,
            min_seed_hits,
            verify_resources,
            s3_region,
            s3_endpoint,
            s3_path_style,
        } => {
            let s3 = if let Some(region) = s3_region {
                Some(range_source::S3Config::new(
                    &region,
                    s3_endpoint.as_deref(),
                    s3_path_style,
                    s3::creds::Credentials::default()?,
                )?)
            } else {
                None
            };
            handle_trace_command(TraceArgs {
                query,
                input: match (owner_index, collection, database, index, manifest) {
                    (Some(root), None, None, None, None) => TraceInput::Owner(root),
                    (None, Some(root), None, None, None) => TraceInput::Collection(root),
                    (None, None, Some(database), Some(index), Some(manifest)) => {
                        TraceInput::Shard {
                            database,
                            index,
                            manifest,
                        }
                    }
                    _ => {
                        return Err(anyhow::anyhow!(
                            "Use --owner-index, --collection, or --database, --index, and --manifest"
                        ));
                    }
                },
                audit_index,
                output,
                query_id,
                config: trace::TraceConfig {
                    min_containment,
                    max_metagenomes: if max_metagenomes == 0 {
                        usize::MAX
                    } else {
                        max_metagenomes
                    },
                    use_sketch: !no_sketch,
                    min_seed_hits,
                    circular: !linear,
                    verify_resources,
                    ..trace::TraceConfig::default()
                },
                s3,
                force: cli.force,
            })
        }

        Commands::Stats { input, short, full } => {
            handle_stats_command(input, short, full, cli.silent)
        }
    }
}
