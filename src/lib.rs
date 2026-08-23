pub mod bias;
pub mod cli;
pub mod core_utils;
pub mod format;
pub mod io;
pub mod jma;
pub mod provenance;
pub mod query;
pub mod reader;
pub mod resource;
pub mod screen;
pub mod sketch;
pub mod trace;
pub mod writer;
pub use cli::handlers::{
    handle_archive_command, handle_bias_create_command, handle_bias_stats_command,
    handle_distance_command, handle_screen_command, handle_sketch_command, handle_stats_command,
    handle_trace_command,
};
pub use io::{expand_input_paths, is_sequence_file};
pub use jamhash::jamhash_u64;
pub use jamhash::jamhash_u64 as jamhash_u64_v1;

use anyhow::Result;
use clap::Parser;
use cli::{BiasCommands, Cli, Commands, QueryKindArg, TopologyArg, TraceSensitivityArg};

pub fn run() -> Result<()> {
    let cli = Cli::parse();
    let threads = cli.threads.unwrap_or(1);
    let memory_target = cli.memory_target.unwrap_or(2);
    if threads == 0 {
        anyhow::bail!("--threads must be > 0");
    }
    if memory_target == 0 {
        anyhow::bail!("--memory-target must be > 0");
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(8 * 1024 * 1024)
        .build_global()?;

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
                threads,
                memory_target,
                cli.force,
                cli.silent,
                complexity,
                temp_dir,
                bias_table,
            )
        }

        Commands::Screen {
            input,
            database,
            output,
            summary,
            metadata,
            assembly_name,
            min_shared,
            min_query_containment,
            min_reference_containment,
            top_per_contig,
            top_references,
        } => handle_screen_command(
            crate::screen::ScreenConfig {
                input,
                database,
                output,
                summary,
                metadata,
                assembly_name,
                min_shared,
                min_query_containment,
                min_reference_containment,
                top_per_contig,
                top_references,
                threads,
                memory_gb: memory_target,
                force: cli.force,
            },
            cli.silent,
        ),

        Commands::Archive {
            input,
            output,
            block_bases,
            primary_scale,
            rescue_scale,
            no_rescue,
            complexity,
        } => handle_archive_command(
            input,
            output,
            block_bases,
            primary_scale,
            (!no_rescue).then_some(rescue_scale),
            complexity,
            cli.force,
            cli.silent,
        ),

        Commands::Trace {
            query,
            plasmid,
            query_kind,
            topology,
            database,
            metagenomes,
            output,
            upload_to,
            query_id,
            sensitivity,
            min_shared,
            min_query_containment,
            min_metagenome_containment,
            top_candidates,
            max_alignments,
            io_concurrency,
            topology_margin_bases,
            cache_dir,
            cache_block_bytes,
            request_timeout_seconds,
            max_retries,
            no_full_download_fallback,
        } => {
            let used_plasmid_alias = plasmid.is_some();
            let query = query
                .or(plasmid)
                .ok_or_else(|| anyhow::anyhow!("--query is required"))?;
            let query_kind = if used_plasmid_alias {
                if !matches!(query_kind, QueryKindArg::Unknown | QueryKindArg::Plasmid) {
                    anyhow::bail!("--plasmid implies --query-kind plasmid");
                }
                eprintln!("warning: --plasmid is deprecated; use --query and --query-kind plasmid");
                trace::model::QueryKind::Plasmid
            } else {
                match query_kind {
                    QueryKindArg::Plasmid => trace::model::QueryKind::Plasmid,
                    QueryKindArg::Phage => trace::model::QueryKind::Phage,
                    QueryKindArg::Other => trace::model::QueryKind::Other,
                    QueryKindArg::Unknown => trace::model::QueryKind::Unknown,
                }
            };
            let topology = match topology {
                TopologyArg::Linear => trace::model::TopologyRequested::Linear,
                TopologyArg::Circular => trace::model::TopologyRequested::Circular,
                TopologyArg::Auto => trace::model::TopologyRequested::Auto,
                TopologyArg::Unknown => trace::model::TopologyRequested::Unknown,
            };
            handle_trace_command(
                query,
                query_kind,
                topology,
                database,
                metagenomes,
                output,
                upload_to,
                query_id,
                match sensitivity {
                    TraceSensitivityArg::Fast => trace::config::SensitivityProfile::Fast,
                    TraceSensitivityArg::Balanced => trace::config::SensitivityProfile::Balanced,
                    TraceSensitivityArg::Sensitive => trace::config::SensitivityProfile::Sensitive,
                },
                min_shared,
                min_query_containment,
                min_metagenome_containment,
                top_candidates,
                max_alignments,
                threads,
                io_concurrency.unwrap_or(threads),
                topology_margin_bases,
                memory_target,
                cache_dir,
                cache_block_bytes,
                request_timeout_seconds,
                max_retries,
                !no_full_download_fallback,
                cli.force,
                cli.silent,
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
                min_positive_retention,
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
                min_positive_retention,
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
            memory_target,
        ),

        Commands::Stats {
            input,
            short,
            full,
            json,
        } => handle_stats_command(input, short, full, json, cli.silent),
    }
}
