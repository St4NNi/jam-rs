pub mod archive;
pub mod bias;
pub mod cli;
pub mod core_utils;
pub mod format;
pub mod io;
pub mod jam_index;
pub mod jma;
pub mod profiling;
pub mod provenance;
pub mod query;
pub mod reader;
pub mod resource;
pub mod router;
pub mod screen;
pub mod sequence;
pub mod sketch;
pub mod trace;
pub mod writer;
pub use cli::handlers::{
    IndexBatchTraceArgs, IndexBuildArgs, IndexTraceArgs, handle_archive_command,
    handle_bias_create_command, handle_bias_stats_command, handle_distance_command,
    handle_index_append, handle_index_batch_trace, handle_index_build, handle_index_build_fragment,
    handle_index_diagnose_spatial, handle_index_finalize, handle_index_merge_part,
    handle_index_plan, handle_index_trace, handle_screen_command, handle_sketch_command,
    handle_stats_command, handle_trace_command,
};
pub use io::{expand_input_paths, is_sequence_file};
pub use jamhash::jamhash_u64;
pub use jamhash::jamhash_u64 as jamhash_u64_v1;

use anyhow::Result;
use clap::Parser;
use cli::{
    ArchiveBlockCodecArg, ArchiveBlockPolicyArg, ArchiveGearTableArg, BiasCommands, Cli, Commands,
    IndexCommands, IndexScreenPolicyArg, QueryKindArg, TopologyArg, TraceSensitivityArg,
};

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
            sequence_block_policy,
            sequence_block_codec,
            gear_min_bases,
            gear_target_bases,
            gear_max_bases,
            gear_table,
            primary_scale,
            rescue_scale,
            no_rescue,
            complexity,
        } => handle_archive_command(
            input,
            output,
            match sequence_block_policy {
                ArchiveBlockPolicyArg::Fixed => {
                    crate::sequence::SequenceBlockPolicy::Fixed { block_bases }
                }
                ArchiveBlockPolicyArg::Gear => crate::sequence::SequenceBlockPolicy::Gear {
                    min_bases: gear_min_bases,
                    target_bases: gear_target_bases,
                    max_bases: gear_max_bases,
                    table: match gear_table {
                        ArchiveGearTableArg::SingleBase => crate::sequence::GearTable::SingleBase,
                        ArchiveGearTableArg::Dinucleotide => {
                            crate::sequence::GearTable::Dinucleotide
                        }
                        ArchiveGearTableArg::PackedFourBase => {
                            crate::sequence::GearTable::PackedFourBase
                        }
                    },
                },
            },
            match sequence_block_codec {
                ArchiveBlockCodecArg::Raw2bit => crate::sequence::BlockCodec::Raw2Bit,
                ArchiveBlockCodecArg::Zstd2bit => crate::sequence::BlockCodec::Zstd2Bit,
            },
            primary_scale,
            (!no_rescue).then_some(rescue_scale),
            complexity,
            cli.force,
            cli.silent,
        ),

        Commands::Index { command } => match command {
            IndexCommands::Plan {
                metagenomes,
                output,
                parts,
                fragments_per_part,
                estimated_expansion,
                screen_policy,
                adaptive_second_minimum_threshold,
                whole_metagenome_hashes,
            } => handle_index_plan(
                metagenomes,
                output,
                parts,
                fragments_per_part,
                estimated_expansion,
                match (screen_policy, adaptive_second_minimum_threshold) {
                    (IndexScreenPolicyArg::Baseline, None) => {
                        crate::jam_index::ScreenSelectionPolicy::default_signatures()
                    }
                    (IndexScreenPolicyArg::Spatial256One, None) => {
                        crate::jam_index::ScreenSelectionPolicy::spatial_256(
                            whole_metagenome_hashes,
                        )
                    }
                    (IndexScreenPolicyArg::Spatial256One, Some(threshold)) => {
                        crate::jam_index::ScreenSelectionPolicy::spatial_256_adaptive(
                            threshold,
                            whole_metagenome_hashes,
                        )?
                    }
                    (IndexScreenPolicyArg::Spatial256Two, None) => {
                        crate::jam_index::ScreenSelectionPolicy::spatial_256_two(
                            whole_metagenome_hashes,
                        )
                    }
                    (_, Some(_)) => {
                        return Err(anyhow::anyhow!(
                            "--adaptive-second-minimum-threshold requires --screen-policy spatial256-one"
                        ));
                    }
                },
                cli.force,
                cli.silent,
            ),
            IndexCommands::BuildFragment {
                plan,
                fragment_id,
                staged_metagenomes,
                output,
            } => handle_index_build_fragment(
                plan,
                fragment_id,
                staged_metagenomes,
                output,
                cli.force,
                cli.silent,
            ),
            IndexCommands::MergePart {
                plan,
                part_id,
                fragments_root,
                output,
            } => handle_index_merge_part(
                plan,
                part_id,
                fragments_root,
                output,
                cli.force,
                cli.silent,
            ),
            IndexCommands::Finalize { plan, output } => {
                handle_index_finalize(plan, output, cli.force, cli.silent)
            }
            IndexCommands::DiagnoseSpatial {
                index,
                source_catalog,
                queries,
                query_id,
                metagenome_id,
                contig_header,
                query_start,
                rare_rescue_df,
                output,
            } => handle_index_diagnose_spatial(
                index,
                source_catalog,
                queries,
                query_id,
                metagenome_id,
                contig_header,
                query_start,
                rare_rescue_df,
                output,
                cli.force,
            ),
            IndexCommands::Build {
                metagenomes,
                output,
                max_part_bases,
                max_part_signatures,
                parts,
                parallel_parts,
            } => handle_index_build(IndexBuildArgs {
                metagenomes,
                output,
                policy: crate::jam_index::ScreenSelectionPolicy::default_signatures(),
                max_part_bases,
                max_part_signatures,
                target_parts: parts,
                parallel_parts: parallel_parts.unwrap_or(threads),
                force: cli.force,
                silent: cli.silent,
            }),
            IndexCommands::Append {
                metagenomes,
                output,
                max_part_bases,
                max_part_signatures,
                parts,
                parallel_parts,
            } => handle_index_append(IndexBuildArgs {
                metagenomes,
                output,
                policy: crate::jam_index::ScreenSelectionPolicy::default_signatures(),
                max_part_bases,
                max_part_signatures,
                target_parts: parts,
                parallel_parts: parallel_parts.unwrap_or(threads),
                force: cli.force,
                silent: cli.silent,
            }),
        },

        Commands::Trace {
            query,
            plasmid,
            query_kind,
            topology,
            database,
            index,
            metagenomes,
            output,
            upload_to,
            query_id,
            sensitivity,
            min_shared,
            min_query_windows,
            rare_rescue_df,
            hamming1_rescue,
            whole_sample_min_shared,
            screen_only,
            min_query_containment,
            min_metagenome_containment,
            top_candidates,
            initial_contigs,
            max_contigs,
            max_contig_bases,
            expansion_batch,
            max_alignments,
            max_group_contig_bases,
            fallback_contigs_per_chunk,
            io_concurrency,
            topology_margin_bases,
            cache_dir,
            cache_block_bytes,
            request_timeout_seconds,
            max_retries,
        } => {
            let used_plasmid_alias = plasmid.is_some();
            let query = if let Some(plasmid) = plasmid {
                vec![plasmid]
            } else {
                query
            };
            if query.is_empty() {
                anyhow::bail!("--query is required");
            }
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
            let profile = match sensitivity {
                TraceSensitivityArg::Fast => trace::config::SensitivityProfile::Fast,
                TraceSensitivityArg::Balanced => trace::config::SensitivityProfile::Balanced,
                TraceSensitivityArg::Sensitive => trace::config::SensitivityProfile::Sensitive,
            };
            if let Some(index) = index {
                if upload_to.is_some() {
                    anyhow::bail!("Jam Index trace is local-only and cannot use --upload-to");
                }
                handle_index_batch_trace(IndexBatchTraceArgs {
                    queries: query,
                    query_id,
                    query_kind,
                    topology,
                    index,
                    source_catalog: metagenomes,
                    output,
                    candidates: None,
                    work: None,
                    status: None,
                    profile,
                    min_shared,
                    min_query_windows,
                    rare_rescue_df,
                    hamming1_rescue,
                    whole_sample_min_shared,
                    screen_only,
                    top_candidates,
                    initial_contigs,
                    max_contigs,
                    max_contig_bases,
                    expansion_batch,
                    max_alignments,
                    max_group_contig_bases,
                    fallback_contigs_per_chunk,
                    threads,
                    topology_margin: topology_margin_bases,
                    memory_gb: memory_target,
                    force: cli.force,
                    silent: cli.silent,
                })
            } else {
                if query.len() != 1 {
                    anyhow::bail!("multi-file trace requires --index");
                }
                handle_trace_command(
                    query
                        .into_iter()
                        .next()
                        .expect("one query path was checked"),
                    query_kind,
                    topology,
                    database,
                    None,
                    metagenomes.ok_or_else(|| anyhow::anyhow!("--metagenomes is required"))?,
                    output,
                    upload_to,
                    query_id,
                    profile,
                    min_shared,
                    min_query_containment,
                    min_metagenome_containment,
                    top_candidates,
                    160,
                    0.99,
                    0.01,
                    false,
                    256,
                    crate::router::WitnessHandoffMode::SampleOnly,
                    max_alignments,
                    threads,
                    io_concurrency.unwrap_or(threads),
                    topology_margin_bases,
                    memory_target,
                    cache_dir,
                    cache_block_bytes,
                    request_timeout_seconds,
                    max_retries,
                    cli.force,
                    cli.silent,
                )
            }
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
