//! `jam screen`: stream assemblies past a shard of small sequences and write what each holds.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::rule::RuleArgs;
use super::storage::{CacheArgs, ShardStorage};
use super::{CliError, output_writer, read_list};
use crate::screen::{AssemblyFailure, AssemblySummary, ScreenConfig, ScreenEntry, Screener};
use crate::search::{SearchControl, Stopped};
use clap::Args;
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::AtomicBool;
use std::sync::mpsc::{SyncSender, sync_channel};

/// Records held between the screening threads and the writer.
const RECORD_QUEUE: usize = 64;

/// The options of `jam screen`.
#[derive(Debug, Args)]
pub struct ScreenArgs {
    /// The shard of the small sequences to look for: NAME.jam or s3://bucket/prefix/NAME.jam
    index: PathBuf,
    /// A file with one assembly path per line (plain or gzip FASTA); the file name names the
    /// assembly
    #[arg(long)]
    assemblies: PathBuf,
    #[command(flatten)]
    rule: RuleArgs,
    /// A contig counts only for the entries within this margin of its best one (score: aligned
    /// share of the contig times identity), 0 to 1; without it, for every entry it qualifies for
    #[arg(long)]
    closest_margin: Option<f64>,
    /// Bases of contigs searched together
    #[arg(long, default_value_t = ScreenConfig::default().batch_bases)]
    batch_bases: usize,
    /// Write the contigs that share enough with an entry to this folder, one gzip FASTA for each
    /// assembly
    #[arg(long)]
    contig_dir: Option<PathBuf>,
    /// Stop after the survey: align nothing and report no entries, only write --contig-dir
    #[arg(long, requires = "contig_dir")]
    survey_only: bool,
    /// Threads to use, each screening one assembly at a time [default: all cores]
    #[arg(short, long)]
    threads: Option<usize>,
    /// Write to this file instead of the standard output
    #[arg(short, long)]
    output: Option<PathBuf>,
}

/// One output line of a screen.
enum Line {
    Entry(ScreenEntry),
    Summary(AssemblySummary),
    Failure(AssemblyFailure),
}

/// Screens every assembly of the list on the threads and writes the records as they come, on
/// this thread. Fails after all records are written when an assembly could not be screened.
pub fn run(args: &ScreenArgs) -> Result<(), CliError> {
    if let Some(note) = args.rule.rule().identity_note() {
        eprintln!("jam: {note}");
    }
    let config = ScreenConfig {
        closest_margin: args.closest_margin,
        batch_bases: args.batch_bases,
        contig_dir: args.contig_dir.clone(),
        survey_only: args.survey_only,
        ..ScreenConfig::new(args.rule.rule())
    };
    // A shard on S3 is read through a block cache in RAM.
    let mut storage = ShardStorage::new()?;
    let paths = std::slice::from_ref(&args.index);
    let mut shards = storage.open_all(paths, &CacheArgs::default(), false)?;
    let shard = shards
        .pop()
        .ok_or_else(|| CliError::Usage("no shard".into()))?;
    let screener = Screener::new(&shard, config)?;
    let assemblies = read_list(&args.assemblies)?;
    check_names(&assemblies, args.contig_dir.is_some())?;
    let threads = match args.threads {
        Some(0) => return Err(CliError::Usage("--threads must be at least 1".into())),
        Some(threads) => threads,
        None => 0,
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(crate::search::SearchError::from)?;
    let mut writer = output_writer(args.output.as_deref())?;
    let (sender, receiver) = sync_channel(RECORD_QUEUE);
    let failed = std::thread::scope(|scope| {
        scope.spawn(move || {
            pool.install(|| {
                assemblies
                    .par_iter()
                    .for_each_with(sender, |sender, path| screen_one(&screener, path, sender));
            });
        });
        // A write error returns at once; the dropped receiver then stops the screens.
        let mut failed = 0;
        for line in receiver {
            let text = match &line {
                Line::Entry(entry) => serde_json::to_string(entry)?,
                Line::Summary(summary) => serde_json::to_string(summary)?,
                Line::Failure(failure) => serde_json::to_string(failure)?,
            };
            writeln!(writer, "{text}")?;
            failed += usize::from(matches!(line, Line::Failure(_)));
        }
        Ok::<_, CliError>(failed)
    })?;
    writer.flush()?;
    match failed {
        0 => Ok(()),
        count => Err(CliError::FailedAssemblies(count)),
    }
}

/// The name of an assembly in the records: its file name.
fn assembly_name(path: &Path) -> String {
    path.file_name().map_or_else(
        || path.display().to_string(),
        |name| name.to_string_lossy().into_owned(),
    )
}

/// Fails when two assemblies have the same name, or, with a contig folder, the same contig
/// file there (`x.fa` and `x.fa.gz` both write `x.fa.gz`), as one would replace the other.
fn check_names(assemblies: &[PathBuf], contig_dir: bool) -> Result<(), CliError> {
    let mut seen: HashMap<String, &Path> = HashMap::new();
    for path in assemblies {
        let name = assembly_name(path);
        let key = match name.strip_suffix(".gz") {
            Some(stem) if contig_dir => stem.to_string(),
            _ => name,
        };
        if let Some(other) = seen.insert(key, path) {
            return Err(CliError::Usage(format!(
                "assemblies {} and {} have the same name or contig file",
                other.display(),
                path.display()
            )));
        }
    }
    Ok(())
}

/// Screens one assembly and sends its entry records and its summary, or its failure record.
fn screen_one(screener: &Screener<'_>, path: &Path, lines: &SyncSender<Line>) {
    let name = assembly_name(path);
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let mut emit = |entry| (lines.send(Line::Entry(entry))).map_err(|_| Stopped::Asked.into());
    let line = match screener.screen_assembly(path, &name, control, &mut emit) {
        Ok(summary) => Line::Summary(summary),
        Err(error) => Line::Failure(AssemblyFailure::new(&name, &error)),
    };
    // A dropped receiver means the writer failed and returns its own error.
    let _ = lines.send(line);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn colliding_names_refused() {
        let paths = |names: &[&str]| names.iter().map(PathBuf::from).collect::<Vec<_>>();
        let same = paths(&["/sample-a/assembly.fa", "/sample-b/assembly.fa"]);
        assert!(check_names(&same, false).is_err());
        let gzip = paths(&["/a/x.fa", "/b/x.fa.gz"]);
        assert!(check_names(&gzip, false).is_ok());
        assert!(check_names(&gzip, true).is_err());
        assert!(check_names(&paths(&["/a/x.fa", "/a/y.fa"]), true).is_ok());
    }
}
