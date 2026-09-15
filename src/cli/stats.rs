//! `jam stats`: what a shard holds.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::CliError;
use crate::index::{Parameters, Shard, ShardStats};
use clap::Args;
use serde::Serialize;
use std::path::PathBuf;

/// The options of `jam stats`.
#[derive(Debug, Args)]
pub struct StatsArgs {
    /// The shard directory
    shard: PathBuf,
    /// Print JSON instead of a short table
    #[arg(long)]
    json: bool,
}

/// The parameters and the counts of a shard, as one JSON object.
#[derive(Serialize)]
struct StatsRecord {
    parameters: Parameters,
    #[serde(flatten)]
    stats: ShardStats,
}

/// Prints the stats of the shard.
pub fn run(args: &StatsArgs) -> Result<(), CliError> {
    let shard = Shard::open(&args.shard).map_err(|source| CliError::Shard {
        path: args.shard.clone(),
        source,
    })?;
    let record = StatsRecord {
        parameters: shard.manifest().parameters,
        stats: shard.stats(),
    };
    if args.json {
        println!("{}", serde_json::to_string(&record)?);
    } else {
        print!("{}", table(&record));
    }
    Ok(())
}

fn table(record: &StatsRecord) -> String {
    let parameters = &record.parameters;
    let rows = [
        ("k", parameters.scheme.k.to_string()),
        ("minimizer window", parameters.scheme.w.to_string()),
        ("bin bases", (1u64 << parameters.bin_bits).to_string()),
        ("entries", record.stats.entries.to_string()),
        ("contigs", record.stats.contigs.to_string()),
        ("bases", record.stats.bases.to_string()),
        ("seeds", record.stats.seeds.to_string()),
        ("color sets", record.stats.color_sets.to_string()),
        ("seed locations", record.stats.seed_locations.to_string()),
        ("sketch bytes", record.stats.sketch_bytes.to_string()),
        ("locator bytes", record.stats.locator_bytes.to_string()),
        ("packed bytes", record.stats.packed_bytes.to_string()),
    ];
    rows.iter()
        .map(|(name, value)| format!("{name:<18}{value}\n"))
        .collect()
}
