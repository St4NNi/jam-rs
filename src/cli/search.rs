//! `jam search`: search every record of a query file in one or more shards, or on a server, and
//! write the matches.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::rule::RuleArgs;
use super::storage::{CacheArgs, ShardStorage, shard_name};
use super::{CliError, output_writer};
use crate::index::Shard;
use crate::loss::MaxMiss;
use crate::report::{
    QueryTotals, Record, ReportError, ShardRecord, Status, TSV_HEADER, paf_lines, result_line,
    tsv_line,
};
use crate::search::{
    EntryEvent, SearchConfig, SearchControl, SearchStage, SearchSummary, Searcher,
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

/// The output format.
#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum Format {
    /// One JSON line per reported entry
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
    /// The miss target: the total chance, per query and matching copy, that jam may lose a copy
    /// through its seed lookups and the pruning of candidates. jam can prove it only for some
    /// query lengths and identities; a query that cannot reach it is searched without pruning
    #[arg(long, default_value_t = MaxMiss::default().chance())]
    max_miss: f64,
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
    /// Write a progress line per query to the standard error, with its status, counts and miss
    /// chance; twice (-vv) also the events and full summaries as JSON lines
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
    let searchers = shards
        .iter()
        .map(|open| Searcher::new(&open.shard))
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
    let from = searchers
        .first()
        .and_then(|searcher| searcher.target_length(&config));
    sink.hint(&config, from)?;
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
    settings.max_miss = args.max_miss;
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
    /// The queries searched, those over the target, and the largest miss chance among them.
    seen: MissCount,
}

/// How many queries missed `--max-miss`.
#[derive(Default)]
struct MissCount {
    queries: u64,
    over: u64,
    worst: f64,
}

impl MissCount {
    /// The hint line for stderr, if any query was over the target. `from` is the length from
    /// which queries reach it, when known.
    fn hint(&self, max_miss: f64, identity: f64, from: Option<u64>) -> Option<String> {
        if self.over == 0 {
            return None;
        }
        let (over, queries) = (self.over, self.queries);
        let worst = if self.worst >= 0.01 {
            format!("{:.2}", self.worst)
        } else {
            format!("{:.1e}", self.worst)
        };
        let mut line = format!(
            "jam: {over} of {queries} queries are above the miss target --max-miss {max_miss} at -i {identity} \
             (worst-case bound up to {worst})"
        );
        if let Some(length) = from {
            line += &format!(
                "; at this identity the bound reaches it from about {}",
                rounded(length)
            );
        }
        Some(line + ". -v shows each query.")
    }
}

/// A length in bases rounded up to three digits, as bases or kb.
fn rounded(length: u64) -> String {
    let digits = length.to_string().len() as u32;
    let step = 10u64.pow(digits.saturating_sub(if length < 10_000 { 3 } else { 2 }));
    let up = length.div_ceil(step) * step;
    if up < 10_000 {
        let text = up.to_string();
        let (head, tail) = text.split_at(text.len() - 3);
        format!(
            "{head}{}{tail} bases",
            if head.is_empty() { "" } else { "," }
        )
    } else {
        format!("{} kb", up / 1_000)
    }
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
            seen: MissCount::default(),
        })
    }

    /// Writes the hint on the standard error when some queries missed `--max-miss`.
    fn hint(&self, config: &SearchConfig, from: Option<u64>) -> Result<(), CliError> {
        let (miss, identity) = (config.max_miss.chance(), config.rule.min_identity);
        if let Some(line) = self.seen.hint(miss, identity, from) {
            writeln!(std::io::stderr(), "{line}")?;
        }
        Ok(())
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

    /// Flushes, then writes the progress line with -v, else the problem line of a query that did
    /// not complete, on the standard error. Returns 1 when the query did not complete, else 0.
    fn finish_query(&mut self) -> Result<usize, CliError> {
        self.writer.flush()?;
        self.seen.queries += 1;
        if self.totals.over_target() {
            self.seen.over += 1;
            self.seen.worst = self.seen.worst.max(self.totals.miss());
        }
        if self.verbose >= 1 {
            writeln!(std::io::stderr(), "{}", self.totals.progress())?;
        } else if let Some(problem) = self.totals.problem() {
            writeln!(std::io::stderr(), "{problem}")?;
        }
        Ok(usize::from(self.totals.status() != Status::Complete))
    }
}

#[cfg(test)]
mod tests {
    use super::MissCount;
    use crate::cli::Cli;
    use clap::Parser;

    #[test]
    fn hint_when_over() {
        let mut seen = MissCount {
            queries: 10,
            ..MissCount::default()
        };
        assert_eq!(seen.hint(0.0001, 0.8, Some(600_000)), None);
        seen.over = 3;
        seen.worst = 0.39;
        let line = seen.hint(0.0001, 0.8, Some(598_123)).unwrap();
        assert!(line.starts_with(
            "jam: 3 of 10 queries are above the miss target --max-miss 0.0001 at -i 0.8"
        ));
        assert!(line.contains("up to 0.39"));
        assert!(line.contains("from about 600 kb."), "{line}");
        assert!(
            seen.hint(0.0001, 0.95, Some(1_248))
                .unwrap()
                .contains("about 1,250 bases.")
        );
    }

    #[test]
    fn seed_mismatches_rejected() {
        let line = [
            "jam",
            "search",
            "shard",
            "query.fa",
            "--seed-mismatches",
            "1",
        ];
        assert!(Cli::try_parse_from(line).is_err());
    }
}
