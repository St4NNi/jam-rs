//! `jam search`: search every record of a query file in one or more shards, or on a server, and
//! write the records.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::CliError;
use crate::decide::QueryCoverage;
use crate::index::Shard;
use crate::loss::{LossBudget, MissModel};
use crate::report::{ContigInfo, Ending, Record, ReportError, ShardRecord, paf_lines};
use crate::search::{
    Completion, EntryEvent, SearchConfig, SearchControl, SearchStage, SearchSummary, Searcher,
};
use crate::service::{SearchRequest, SearchSettings, ServiceError, bind, parse_address, search};
use clap::{Args, ValueEnum};
use iroh::{Endpoint, EndpointAddr, SecretKey};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::AtomicBool;
use std::sync::mpsc::sync_channel;
use std::time::{Duration, Instant};

/// Events held between the search and the writer.
const EVENT_QUEUE: usize = 64;

/// The last stage that runs.
#[derive(Clone, Copy, Debug, ValueEnum)]
enum Level {
    /// Stage 1 only: the entries that may hold the query
    Sketch,
    /// Stages 1 and 2: where in those entries the query may lie
    Locate,
    /// All three stages: aligned and decided entries
    Align,
}

/// The output format.
#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum Format {
    /// One JSON line per entry event, then one summary line per query
    Json,
    /// One PAF line per fragment of a confirmed entry
    Paf,
}

/// The options of `jam search`.
#[derive(Debug, Args)]
pub struct SearchArgs {
    /// One or more shard directories, then a FASTA file with the query records; each record is
    /// searched on its own, in every shard. With --server only the query file
    #[arg(required = true, num_args = 1.., value_name = "SHARD... QUERY")]
    paths: Vec<PathBuf>,
    /// Search on the server at NODE_ID@ADDRESS, as `jam serve` prints it, instead of local shards
    #[arg(long, value_parser = parse_address)]
    server: Option<EndpointAddr>,
    /// Treat the queries as circular instead of linear
    #[arg(long)]
    circular: bool,
    /// Least query coverage: a share such as 0.5, or a base count such as 800bp [default: 0.5]
    #[arg(short, long)]
    coverage: Option<QueryCoverage>,
    /// Least identity over all counted fragments together [default: 0.8]
    #[arg(short, long)]
    identity: Option<f64>,
    /// Least identity of one fragment to count [default: 0]
    #[arg(long)]
    fragment_identity: Option<f64>,
    /// Least alignment columns of one fragment to count [default: 100]
    #[arg(long)]
    fragment_columns: Option<u64>,
    /// Contig rule: least share of a contig's bases that must align [default: 0, rule off]
    #[arg(long)]
    contig_coverage: Option<f64>,
    /// Contig rule: least identity of the fragments of a contig [default: 0, rule off]
    #[arg(long)]
    contig_identity: Option<f64>,
    /// Largest chance of missing a qualifying entry [default: 0.0001]
    #[arg(long)]
    loss_budget: Option<f64>,
    /// The last stage that runs
    #[arg(long, value_enum, default_value = "align")]
    level: Level,
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

fn open_shards(directories: &[PathBuf]) -> Result<Vec<NamedShard>, CliError> {
    let mut shards: Vec<NamedShard> = Vec::new();
    for directory in directories {
        let name = directory.file_name().map_or_else(
            || directory.display().to_string(),
            |n| n.to_string_lossy().into_owned(),
        );
        if shards.iter().any(|open| open.name == name) {
            return Err(CliError::Usage(format!(
                "two shards are named {name}; records could not tell them apart"
            )));
        }
        let shard = Shard::open(directory).map_err(|source| CliError::Shard {
            path: directory.clone(),
            source,
        })?;
        shards.push(NamedShard { name, shard });
    }
    Ok(shards)
}

/// Searches every query record in every shard and writes the records. Fails when a search did
/// not complete, after all records are written.
pub fn run(args: &SearchArgs) -> Result<(), CliError> {
    if let Some(server) = &args.server {
        return run_remote(args, server);
    }
    let config = config_of(args)?;
    let inputs = args.inputs()?;
    let shards = open_shards(&inputs.shards)?;
    let parameters = shards[0].shard.manifest().parameters;
    if shards.iter().any(|open| {
        let other = open.shard.manifest().parameters;
        (other.k, other.window) != (parameters.k, parameters.window)
    }) {
        return Err(CliError::Usage(
            "the shards differ in k or minimizer window".into(),
        ));
    }
    let model = MissModel::new(usize::from(parameters.k), usize::from(parameters.window))?;
    let searchers = shards
        .iter()
        .map(|open| Searcher::new(&open.shard, &model))
        .collect::<Result<Vec<_>, _>>()?;
    let mut sink = Sink::new(args)?;
    let incomplete = each_query(inputs.query, |name, sequence| {
        let mut incomplete = 0;
        for (open, searcher) in shards.iter().zip(&searchers) {
            let query = Query { name, sequence };
            let summary = search_one(args, searcher, open, &config, &mut sink, &query)?;
            sink.finish_query(args.format, &open.name, name, &summary)?;
            incomplete += usize::from(summary.completion != Completion::Complete);
        }
        Ok(incomplete)
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
    let config = config_of(args)?;
    let settings = SearchSettings::of(&config);
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
        runtime.block_on(search_remote(
            &endpoint,
            server,
            &request,
            args.format,
            &mut sink,
        ))
    })?;
    match incomplete {
        0 => Ok(()),
        count => Err(CliError::Incomplete(count)),
    }
}

/// Writes the records of one remote search as they arrive; returns how many shards did not
/// complete.
async fn search_remote(
    endpoint: &Endpoint,
    server: &EndpointAddr,
    request: &SearchRequest,
    format: Format,
    sink: &mut Sink,
) -> Result<usize, CliError> {
    let mut answer = search(endpoint, server.clone(), request).await?;
    let mut incomplete = 0;
    while let Some(record) = answer.next().await?.and_then(|reply| reply.into_record()) {
        if let Record::Summary(summary) = &record.record {
            incomplete += usize::from(summary.ending != Ending::Complete);
        }
        sink.write_lines(&record_lines(format, &record)?)?;
    }
    sink.writer.flush()?;
    Ok(incomplete)
}

/// Calls `search` with the name (first word of the id) and the upper-case bases of every record
/// of the query file, and sums the incomplete searches it returns. Fails on a file without
/// records.
fn each_query(
    path: &Path,
    mut search: impl FnMut(&str, &[u8]) -> Result<usize, CliError>,
) -> Result<usize, CliError> {
    let query_error = |message: String| CliError::Query {
        path: path.to_path_buf(),
        message,
    };
    let mut reader =
        needletail::parse_fastx_file(path).map_err(|error| query_error(error.to_string()))?;
    let (mut records, mut incomplete) = (0, 0);
    while let Some(record) = reader.next() {
        let record = record.map_err(|error| query_error(error.to_string()))?;
        records += 1;
        incomplete += search(&first_word(record.id()), &record.seq().to_ascii_uppercase())?;
    }
    match records {
        0 => Err(query_error("it holds no records".into())),
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
            sink.event(args.format, open, query.name, &event)?;
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

fn config_of(args: &SearchArgs) -> Result<SearchConfig, CliError> {
    let mut config = SearchConfig {
        level: match args.level {
            Level::Sketch => SearchStage::Sketch,
            Level::Locate => SearchStage::Locate,
            Level::Align => SearchStage::Align,
        },
        ..SearchConfig::default()
    };
    let rule = &mut config.rule;
    rule.min_query_coverage = args.coverage.unwrap_or(rule.min_query_coverage);
    rule.min_identity = args.identity.unwrap_or(rule.min_identity);
    rule.min_fragment_identity = args.fragment_identity.unwrap_or(rule.min_fragment_identity);
    rule.min_fragment_columns = args.fragment_columns.unwrap_or(rule.min_fragment_columns);
    rule.min_contig_coverage = args.contig_coverage.unwrap_or(rule.min_contig_coverage);
    rule.min_contig_identity = args.contig_identity.unwrap_or(rule.min_contig_identity);
    if let Some(chance) = args.loss_budget {
        config.loss_budget = LossBudget::new(chance)?;
    }
    if args.threads == Some(0) {
        return Err(CliError::Usage("--threads must be at least 1".into()));
    }
    config.threads = args.threads.unwrap_or(0);
    if args
        .deadline
        .is_some_and(|seconds| !(seconds > 0.0 && seconds.is_finite()))
    {
        return Err(CliError::Usage(
            "--deadline must be a number of seconds above 0".into(),
        ));
    }
    config.validate()?;
    Ok(config)
}

fn first_word(id: &[u8]) -> String {
    let id = String::from_utf8_lossy(id);
    id.split_whitespace().next().unwrap_or_default().to_string()
}

/// The output of a search.
struct Sink {
    writer: Box<dyn Write>,
}

impl Sink {
    fn new(args: &SearchArgs) -> Result<Self, CliError> {
        let writer: Box<dyn Write> = match &args.output {
            Some(path) => Box::new(BufWriter::new(File::create(path).map_err(|source| {
                CliError::File {
                    path: path.clone(),
                    source,
                }
            })?)),
            None => Box::new(std::io::stdout()),
        };
        Ok(Self { writer })
    }

    fn event(
        &mut self,
        format: Format,
        open: &NamedShard,
        query: &str,
        event: &EntryEvent,
    ) -> Result<(), CliError> {
        let index = event.entry();
        let shard = &open.shard;
        let entry = shard
            .manifest()
            .entries
            .get(index as usize)
            .ok_or(ReportError::UnknownEntry(index))?;
        let record = ShardRecord {
            shard: open.name.clone(),
            record: Record::from_event(query, &entry.name, event, || contigs_of(shard, index))?,
        };
        self.write_lines(&record_lines(format, &record)?)
    }

    fn write_lines(&mut self, lines: &[String]) -> Result<(), CliError> {
        for line in lines {
            writeln!(self.writer, "{line}")?;
        }
        Ok(())
    }

    /// Writes the summary line of a query in a shard (JSON only) and flushes.
    fn finish_query(
        &mut self,
        format: Format,
        shard: &str,
        query: &str,
        summary: &SearchSummary,
    ) -> Result<(), CliError> {
        if format == Format::Json {
            let line = json_line(shard, Record::summary(query, summary))?;
            self.write_lines(&[line])?;
        }
        self.writer.flush()?;
        Ok(())
    }
}

/// The output lines of one record: its JSON line, or the PAF lines of a confirmed entry.
fn record_lines(format: Format, record: &ShardRecord) -> Result<Vec<String>, CliError> {
    Ok(match (format, &record.record) {
        (Format::Json, _) => vec![serde_json::to_string(record)?],
        (Format::Paf, Record::Confirmed(confirmed)) => {
            paf_lines(&confirmed.query, &confirmed.decision)
        }
        (Format::Paf, _) => Vec::new(),
    })
}

fn json_line(shard: &str, record: Record) -> Result<String, CliError> {
    let record = ShardRecord {
        shard: shard.to_string(),
        record,
    };
    Ok(serde_json::to_string(&record)?)
}

/// The names and lengths of the contigs of an entry.
fn contigs_of(shard: &Shard, entry: u32) -> Result<Vec<ContigInfo>, CliError> {
    let packed = shard.packed(entry)?;
    (0..packed.contig_count())
        .map(|contig| {
            let range = packed.contig_range(contig)?;
            Ok(ContigInfo {
                name: String::from_utf8_lossy(&packed.contig_name(contig)?).into_owned(),
                length: range.end - range.start,
            })
        })
        .collect()
}
