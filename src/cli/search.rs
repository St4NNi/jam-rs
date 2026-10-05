//! `jam search`: search every record of a query file in one or more shards, or on a server, and
//! write the records.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::rule::RuleArgs;
use super::storage::{CacheArgs, ShardStorage, shard_name};
use super::{CliError, output_writer};
use crate::index::Shard;
use crate::loss::{LossBudget, MissModel};
use crate::report::{
    QueryTotals, Record, ReportError, ShardRecord, Status, TSV_HEADER, paf_lines, result_line,
    tsv_line,
};
use crate::search::{
    EntryEvent, MissBound, SearchConfig, SearchControl, SearchStage, SearchSummary, Searcher,
    SeedMismatches,
};
use crate::service::{SearchRequest, SearchSettings, ServiceError, bind, parse_address, search};
use clap::{Args, ValueEnum};
use iroh::{Endpoint, EndpointAddr, SecretKey};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::AtomicBool;
use std::sync::mpsc::sync_channel;
use std::time::{Duration, Instant};

/// Events held between the search and the writer.
const EVENT_QUEUE: usize = 64;

/// The long help of `--seed-mismatches`.
const SEED_MISMATCHES_HELP: &str = "Mismatches a query k-mer may have against a seed. 1 also \
looks up every k-mer one substitution away (3 k lookups per k-mer), which finds entries at lower \
identity. auto takes 1 for a query when both hold: (a) the proven upper bound on the chance that \
exact seeds miss a present piece of the required query coverage, at the --identity value, is \
above the loss budget, and one mismatch lowers it; (b) the extra CPU, estimated as 1 ms per query \
k-mer, stays within 10 CPU seconds. With the default coverage this takes 1 at 80% identity for \
queries up to about 10 kb, and at 95% for queries below about 7.2 kb. 2 also looks up every k-mer \
two substitutions away (2,700 more lookups per k-mer). auto takes 2 instead of 1 when also: the \
identity is below 0.9; the query has at most 1,800 k-mers; the proven bound with one mismatch is \
above the loss budget; and two mismatches lower it. There is no stored estimate table for two \
mismatches, so with 2 both --miss-bound modes use the proven bound.";

/// The output format.
#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum Format {
    /// One JSON line per reported entry, then one summary line per query
    Json,
    /// A header line, then one line per confirmed entry: query, assembly, coverage, identity,
    /// covered bases and fragments
    Tsv,
    /// One PAF line per fragment of a confirmed entry
    Paf,
}

/// The options of `jam search`.
#[derive(Debug, Args)]
pub struct SearchArgs {
    /// One or more index files (NAME.jam or s3://bucket/prefix/NAME.jam), then a FASTA file
    /// with the query records; each record is searched on its own, in every shard. With --server
    /// only the query file
    #[arg(required = true, num_args = 1.., value_name = "SHARD... QUERY")]
    paths: Vec<PathBuf>,
    /// Search on the server at NODE_ID@ADDRESS, as `jam serve` prints it, instead of local shards
    #[arg(long, value_parser = parse_address)]
    server: Option<EndpointAddr>,
    /// Treat the queries as circular instead of linear
    #[arg(long)]
    circular: bool,
    #[command(flatten)]
    rule: RuleArgs,
    #[command(flatten)]
    cache: CacheArgs,
    /// Chance of missing a qualifying entry that the coverage bound may spend
    #[arg(long, default_value_t = LossBudget::default().chance())]
    loss_budget: f64,
    /// How the coverage bound allows for stretches without a seed hit: estimate or proven
    #[arg(long, value_enum, default_value = "estimate")]
    miss_bound: MissBound,
    /// Mismatches a query k-mer may have against a seed: auto, 0, 1 or 2
    #[arg(long, value_enum, default_value = "auto", long_help = SEED_MISMATCHES_HELP)]
    seed_mismatches: SeedMismatches,
    /// The last stage that runs
    #[arg(long, value_enum, default_value = "align")]
    level: SearchStage,
    /// Threads to use [default: all cores]
    #[arg(short, long)]
    threads: Option<usize>,
    /// Seconds after which the search of each query stops
    #[arg(long)]
    deadline: Option<f64>,
    /// The output format
    #[arg(long, value_enum, default_value = "json")]
    format: Format,
    /// Write to this file instead of the standard output
    #[arg(short, long)]
    output: Option<PathBuf>,
    /// Write a progress line per query to the standard error; twice (-vv) also the events and
    /// full summaries as JSON lines
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,
}

/// The shards and the query file of the command line.
struct Inputs<'a> {
    shards: Vec<PathBuf>,
    query: &'a Path,
}

impl SearchArgs {
    fn inputs(&self) -> Result<Inputs<'_>, CliError> {
        let (query, places) = self
            .paths
            .split_last()
            .filter(|(_, places)| !places.is_empty())
            .ok_or_else(|| CliError::Usage("give a shard and a query file".into()))?;
        let mut shards = Vec::new();
        for place in places {
            shards.extend(shard_directories(place)?);
        }
        Ok(Inputs { shards, query })
    }
}

/// The shard directories that one path of the command line stands for. Today a path is one
/// shard. An index directory that lists several shards belongs here: it will return them all.
fn shard_directories(path: &Path) -> Result<Vec<PathBuf>, CliError> {
    Ok(vec![path.to_path_buf()])
}

/// A shard opened for searching, with the name its records carry.
struct NamedShard {
    name: String,
    shard: Shard,
}

fn open_shards(
    directories: &[PathBuf],
    storage: &mut ShardStorage,
    cache: &CacheArgs,
) -> Result<Vec<NamedShard>, CliError> {
    let mut names: Vec<String> = Vec::new();
    for directory in directories {
        let name = shard_name(directory);
        if names.contains(&name) {
            return Err(CliError::Usage(format!(
                "two shards are named {name}; records could not tell them apart"
            )));
        }
        names.push(name);
    }
    let shards = storage.open_all(directories, cache, false)?;
    let named = names.into_iter().zip(shards);
    Ok(named
        .map(|(name, shard)| NamedShard { name, shard })
        .collect())
}

/// Searches every query record in every shard and writes the records. Fails when a search did
/// not complete, after all records are written.
pub fn run(args: &SearchArgs) -> Result<(), CliError> {
    if let Some(server) = &args.server {
        return run_remote(args, server);
    }
    let base = SearchConfig {
        threads: args.threads.unwrap_or(0),
        ..SearchConfig::default()
    };
    let config = settings_of(args)?.apply(&base)?;
    let inputs = args.inputs()?;
    let mut storage = ShardStorage::new()?;
    let shards = open_shards(&inputs.shards, &mut storage, &args.cache)?;
    let parameters = shards[0].shard.manifest().parameters;
    if shards.iter().any(|open| {
        let other = open.shard.manifest().parameters;
        (other.scheme.k, other.scheme.w) != (parameters.scheme.k, parameters.scheme.w)
    }) {
        return Err(CliError::Usage(
            "the shards differ in k or minimizer window".into(),
        ));
    }
    let exact = MissModel::new(parameters.scheme)?;
    let one_mismatch = MissModel::one_mismatch(parameters.scheme)?;
    let searchers = shards
        .iter()
        .map(|open| Searcher::new(&open.shard, &exact, &one_mismatch))
        .collect::<Result<Vec<_>, _>>()?;
    let mut sink = Sink::new(args)?;
    let incomplete = each_query(inputs.query, |name, sequence| {
        sink.totals = QueryTotals::new(name);
        for (open, searcher) in shards.iter().zip(&searchers) {
            let query = Query { name, sequence };
            let summary = search_one(args, searcher, open, &config, &mut sink, &query)?;
            sink.record(&open.name, Record::summary(name, &summary))?;
        }
        sink.finish_query()
    })?;
    match incomplete {
        0 => Ok(()),
        count => Err(CliError::Incomplete(count)),
    }
}

/// Searches every query record on the server and writes its records as a local search does.
/// Fails when a shard did not complete, after all records are written.
fn run_remote(args: &SearchArgs, server: &EndpointAddr) -> Result<(), CliError> {
    let [query] = args.paths.as_slice() else {
        return Err(CliError::Usage(
            "with --server give only the query file".into(),
        ));
    };
    let settings = settings_of(args)?;
    let deadline_ms = args
        .deadline
        .map(|seconds| (seconds * 1_000.0).ceil() as u64);
    let mut sink = Sink::new(args)?;
    let runtime = tokio::runtime::Runtime::new().map_err(ServiceError::from)?;
    let endpoint = runtime.block_on(bind(SecretKey::generate(), &[]))?;
    let incomplete = each_query(query, |name, sequence| {
        let request = SearchRequest {
            query_name: name.to_string(),
            bases: sequence.to_vec(),
            circular: args.circular,
            settings: settings.clone(),
            deadline_ms,
            local_only: false,
        };
        runtime.block_on(search_remote(&endpoint, server, &request, &mut sink))
    })?;
    match incomplete {
        0 => Ok(()),
        count => Err(CliError::Incomplete(count)),
    }
}

/// Writes the records of one remote search as they arrive; returns 1 when a shard did not
/// complete.
async fn search_remote(
    endpoint: &Endpoint,
    server: &EndpointAddr,
    request: &SearchRequest,
    sink: &mut Sink,
) -> Result<usize, CliError> {
    let mut answer = search(endpoint, server.clone(), request).await?;
    sink.totals = QueryTotals::new(&request.query_name);
    while let Some(record) = answer.next().await?.and_then(|reply| reply.into_record()) {
        sink.record(&record.shard, record.record)?;
    }
    sink.finish_query()
}

/// Calls `search` with the name (first word of the id) and the bases of every record
/// of the query file, and sums the incomplete queries it returns. Fails on a file without
/// records.
fn each_query(
    path: &Path,
    mut search: impl FnMut(&str, &[u8]) -> Result<usize, CliError>,
) -> Result<usize, CliError> {
    let query_error = |source| CliError::Query {
        path: path.to_path_buf(),
        source,
    };
    let mut reader = needletail::parse_fastx_file(path).map_err(query_error)?;
    let (mut records, mut incomplete) = (0, 0);
    while let Some(record) = reader.next() {
        let record = record.map_err(query_error)?;
        records += 1;
        incomplete += search(&first_word(record.id()), &record.seq())?;
    }
    match records {
        0 => Err(CliError::NoQueries(path.to_path_buf())),
        _ => Ok(incomplete),
    }
}

/// One query record.
struct Query<'a> {
    name: &'a str,
    sequence: &'a [u8],
}

fn search_one(
    args: &SearchArgs,
    searcher: &Searcher<'_>,
    open: &NamedShard,
    config: &SearchConfig,
    sink: &mut Sink,
    query: &Query<'_>,
) -> Result<SearchSummary, CliError> {
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: args
            .deadline
            .map(|seconds| Instant::now() + Duration::from_secs_f64(seconds)),
        stop: &stop,
    };
    let (sender, receiver) = sync_channel(EVENT_QUEUE);
    // A write error returns at once; the dropped receiver then stops the search.
    std::thread::scope(|scope| {
        let search = scope.spawn(|| {
            let summary = searcher.search(query.sequence, args.circular, config, control, &sender);
            drop(sender);
            summary
        });
        for event in receiver {
            sink.event(open, query.name, &event)?;
        }
        search
            .join()
            .map_err(|_| CliError::Panic)?
            .map_err(|source| CliError::Search {
                query: query.name.to_string(),
                source,
            })
    })
}

/// The search settings of the command line, checked; the threads and the deadline stay local.
fn settings_of(args: &SearchArgs) -> Result<SearchSettings, CliError> {
    if args.threads == Some(0) {
        return Err(CliError::Usage("--threads must be at least 1".into()));
    }
    if args
        .deadline
        .is_some_and(|seconds| !(seconds > 0.0 && seconds.is_finite()))
    {
        return Err(CliError::Usage(
            "--deadline must be a number of seconds above 0".into(),
        ));
    }
    if args.format != Format::Json && args.level != SearchStage::Align {
        return Err(CliError::Usage(
            "--format tsv and paf need --level align".into(),
        ));
    }
    let mut settings = SearchSettings::of(&SearchConfig::default());
    settings.last_stage = args.level;
    settings.rule = args.rule.rule();
    settings.loss_budget = args.loss_budget;
    settings.seed_mismatches = args.seed_mismatches;
    settings.miss_bound = args.miss_bound;
    settings.apply(&SearchConfig::default())?;
    Ok(settings)
}

fn first_word(id: &[u8]) -> String {
    let id = String::from_utf8_lossy(id);
    id.split_whitespace().next().unwrap_or_default().to_string()
}

/// The output of a search: results on the writer, progress and debug lines on the standard
/// error.
struct Sink {
    writer: Box<dyn Write>,
    format: Format,
    level: SearchStage,
    verbose: u8,
    /// The summary of the query being searched.
    totals: QueryTotals,
}

impl Sink {
    fn new(args: &SearchArgs) -> Result<Self, CliError> {
        let mut writer = output_writer(args.output.as_deref())?;
        if args.format == Format::Tsv {
            writeln!(writer, "{TSV_HEADER}")?;
        }
        Ok(Self {
            writer,
            format: args.format,
            level: args.level,
            verbose: args.verbose,
            totals: QueryTotals::new(""),
        })
    }

    fn event(
        &mut self,
        open: &NamedShard,
        query: &str,
        event: &EntryEvent,
    ) -> Result<(), CliError> {
        let index = event.entry();
        let entry = (open.shard.manifest().entries.get(index as usize))
            .ok_or(ReportError::UnknownEntry(index))?;
        self.record(&open.name, Record::from_event(query, &entry.name, event)?)
    }

    /// Writes the result lines of a record of `shard` and its debug line with -vv, and adds a
    /// summary to the totals of the query.
    fn record(&mut self, shard: &str, record: Record) -> Result<(), CliError> {
        if let Record::Summary(summary) = &record {
            self.totals.add(summary);
        }
        let lines = self.result_lines(&record)?;
        if self.verbose >= 2 && !matches!(record, Record::Confirmed(_)) {
            let shard = shard.to_string();
            let debug = serde_json::to_string(&ShardRecord { shard, record })?;
            writeln!(std::io::stderr(), "{debug}")?;
        }
        for line in lines {
            writeln!(self.writer, "{line}")?;
        }
        Ok(())
    }

    /// The output lines of one record: its compact JSON line, or the PAF lines of a confirmed
    /// entry.
    fn result_lines(&self, record: &Record) -> Result<Vec<String>, CliError> {
        Ok(match (self.format, record) {
            (Format::Json, _) => result_line(record, self.level)?.into_iter().collect(),
            (Format::Paf, Record::Confirmed(confirmed)) => {
                paf_lines(&confirmed.query, &confirmed.decision)
            }
            (Format::Tsv, Record::Confirmed(confirmed)) => vec![tsv_line(confirmed)],
            (Format::Tsv | Format::Paf, _) => Vec::new(),
        })
    }

    /// Writes the summary line of a query (JSON only) and its progress line with -v, and
    /// flushes. Returns 1 when the query did not complete, else 0.
    fn finish_query(&mut self) -> Result<usize, CliError> {
        if self.format == Format::Json {
            writeln!(self.writer, "{}", self.totals.line()?)?;
        }
        self.writer.flush()?;
        if self.verbose >= 1 {
            writeln!(std::io::stderr(), "{}", self.totals.progress())?;
        }
        Ok(usize::from(self.totals.status() != Status::Complete))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Command};
    use clap::Parser;

    /// The seed mismatches of a search command line with `extra` flags, or `None` if it fails.
    fn parsed(extra: &[&str]) -> Option<SeedMismatches> {
        let line = ["jam", "search", "shard", "query.fa"].iter().chain(extra);
        match Cli::try_parse_from(line).ok()?.command {
            Command::Search(args) => Some(settings_of(&args).ok()?.seed_mismatches),
            _ => None,
        }
    }

    #[test]
    fn seed_mismatches_parsed() {
        assert_eq!(parsed(&[]), Some(SeedMismatches::Auto));
        assert_eq!(
            parsed(&["--seed-mismatches", "auto"]),
            Some(SeedMismatches::Auto)
        );
        assert_eq!(
            parsed(&["--seed-mismatches", "0"]),
            Some(SeedMismatches::Zero)
        );
        assert_eq!(
            parsed(&["--seed-mismatches", "1"]),
            Some(SeedMismatches::One)
        );
        assert_eq!(
            parsed(&["--seed-mismatches", "2"]),
            Some(SeedMismatches::Two)
        );
        for bad in ["3", "one", ""] {
            assert_eq!(parsed(&["--seed-mismatches", bad]), None, "{bad}");
        }
    }
}
