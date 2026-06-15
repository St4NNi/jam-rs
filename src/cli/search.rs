//! `jam search`: search every record of a query file in a shard and write the records.

use crate::decide::QueryCoverage;
use crate::index::Shard;
use crate::loss::{LossBudget, MissModel};
use crate::report::{ContigInfo, Record, paf_lines};
use crate::search::{
    EntryEvent, SearchConfig, SearchControl, SearchStage, SearchSummary, Searcher,
};
use anyhow::{Context, Result, anyhow, bail};
use clap::{Args, ValueEnum};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::Mutex;
use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};

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
    /// The shard directory
    shard: PathBuf,
    /// A FASTA file with the query records; each record is searched on its own
    query: PathBuf,
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

/// Searches every query record and writes the records.
pub fn run(args: &SearchArgs) -> Result<()> {
    let config = config_of(args)?;
    let shard = Shard::open(&args.shard)
        .with_context(|| format!("cannot open the shard {}", args.shard.display()))?;
    let parameters = shard.manifest().parameters;
    let model = MissModel::new(usize::from(parameters.k), usize::from(parameters.window))?;
    let searcher = Searcher::new(&shard, &model)?;
    let sink = Sink::new(args)?;
    let mut reader = needletail::parse_fastx_file(&args.query)
        .with_context(|| format!("cannot read the query file {}", args.query.display()))?;
    let mut records = 0;
    while let Some(record) = reader.next() {
        let record = record.context("bad record in the query file")?;
        let name = first_word(record.id());
        let sequence = record.seq().to_ascii_uppercase();
        records += 1;
        let summary = search_one(args, &searcher, &shard, &config, &sink, &name, &sequence)
            .with_context(|| format!("searching {name} failed"))?;
        sink.finish_query(args.format, &name, &summary)?;
    }
    if records == 0 {
        bail!("the query file {} holds no records", args.query.display());
    }
    sink.flush()
}

fn search_one(
    args: &SearchArgs,
    searcher: &Searcher<'_>,
    shard: &Shard,
    config: &SearchConfig,
    sink: &Sink,
    name: &str,
    sequence: &[u8],
) -> Result<SearchSummary> {
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: args
            .deadline
            .map(|seconds| Instant::now() + Duration::from_secs_f64(seconds)),
        stop: &stop,
    };
    let on_event = |event: EntryEvent| sink.event(args.format, shard, name, &event);
    let summary = searcher.search(sequence, args.circular, config, control, &on_event)?;
    sink.take_error()?;
    Ok(summary)
}

fn config_of(args: &SearchArgs) -> Result<SearchConfig> {
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
        bail!("--threads must be at least 1");
    }
    config.threads = args.threads.unwrap_or(0);
    if args
        .deadline
        .is_some_and(|seconds| !(seconds > 0.0 && seconds.is_finite()))
    {
        bail!("--deadline must be a number of seconds above 0");
    }
    config.validate()?;
    Ok(config)
}

fn first_word(id: &[u8]) -> String {
    let id = String::from_utf8_lossy(id);
    id.split_whitespace().next().unwrap_or_default().to_string()
}

/// The output and the first error that writing or building a record met. Events arrive from many
/// threads, so each line is written under a lock.
struct Sink {
    writer: Mutex<Box<dyn Write + Send>>,
    error: Mutex<Option<anyhow::Error>>,
}

impl Sink {
    fn new(args: &SearchArgs) -> Result<Self> {
        let writer: Box<dyn Write + Send> = match &args.output {
            Some(path) => Box::new(BufWriter::new(
                File::create(path).with_context(|| format!("cannot create {}", path.display()))?,
            )),
            None => Box::new(std::io::stdout()),
        };
        Ok(Self {
            writer: Mutex::new(writer),
            error: Mutex::new(None),
        })
    }

    fn event(&self, format: Format, shard: &Shard, query: &str, event: &EntryEvent) {
        if self.error.lock().is_ok_and(|error| error.is_some()) {
            return;
        }
        if let Err(error) = self.write_event(format, shard, query, event)
            && let Ok(mut slot) = self.error.lock()
        {
            slot.get_or_insert(error);
        }
    }

    fn write_event(
        &self,
        format: Format,
        shard: &Shard,
        query: &str,
        event: &EntryEvent,
    ) -> Result<()> {
        let index = event.entry();
        let entry = shard
            .manifest()
            .entries
            .get(index as usize)
            .ok_or_else(|| anyhow!("entry {index} is not in the shard"))?;
        let record = Record::from_event(query, &entry.name, event, || contigs_of(shard, index))?;
        let lines = match (format, &record) {
            (Format::Json, record) => vec![serde_json::to_string(record)?],
            (Format::Paf, Record::Confirmed { decision, .. }) => paf_lines(query, decision),
            (Format::Paf, _) => return Ok(()),
        };
        self.write_lines(&lines)
    }

    fn write_lines(&self, lines: &[String]) -> Result<()> {
        let mut writer = self.writer.lock().map_err(|_| anyhow!("output lock"))?;
        for line in lines {
            writeln!(writer, "{line}")?;
        }
        Ok(())
    }

    fn take_error(&self) -> Result<()> {
        match self
            .error
            .lock()
            .map_err(|_| anyhow!("output lock"))?
            .take()
        {
            Some(error) => Err(error),
            None => Ok(()),
        }
    }

    /// Writes the summary line of a query (JSON only) and flushes.
    fn finish_query(&self, format: Format, query: &str, summary: &SearchSummary) -> Result<()> {
        if format == Format::Json {
            let line = serde_json::to_string(&Record::summary(query, summary))?;
            self.write_lines(&[line])?;
        }
        self.flush()
    }

    fn flush(&self) -> Result<()> {
        let mut writer = self.writer.lock().map_err(|_| anyhow!("output lock"))?;
        writer.flush().context("cannot write the output")
    }
}

/// The names and lengths of the contigs of an entry.
fn contigs_of(shard: &Shard, entry: u32) -> Result<Vec<ContigInfo>> {
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
