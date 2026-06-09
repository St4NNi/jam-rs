//! The driver of one search against one shard: it connects the stages, streams the entry events,
//! and sums the work.
//!
//! Stage 1 runs in parallel chunks. The candidate entries then go to a work queue, highest
//! coverage bound first, and each worker takes one entry through stages 2 and 3, so the events of
//! one entry come from one thread, in order. Memory per request stays bounded: stage 1 holds at
//! most `max_seed_hits` seed hits, and stages 2 and 3 hold one entry per worker; results leave as
//! events.

use super::align_stage::{
    AlignStageConfig, AlignStageError, AlignStageWorkspace, AlignWork, AlignedEntry, PreparedQuery,
    align_packed_entry,
};
use super::bound::CoverageBound;
use super::control::SearchControl;
use super::error::SearchError;
use super::events::{Completion, DropReason, EntryEvent, SearchStage};
use super::locate_stage::{
    LocateOutcome, LocateRequest, LocateWork, LocateWorkspace, locate_entry,
};
use super::located::LocatedEntry;
use super::query::QueryKmers;
use super::sketch_stage::{PossibleEntry, SketchLimits, SketchWork, sketch_stage};
use crate::decide::DecisionRule;
use crate::index::Shard;
use crate::loss::{LossBudget, MissModel};
use crate::seeds::Seeder;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::{Duration, Instant};

/// The settings of one search.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SearchConfig {
    /// The last stage that runs.
    pub level: SearchStage,
    /// The rule that reports an entry; its query coverage and identity also drive the pruning.
    pub rule: DecisionRule,
    /// The chance of missing a qualifying entry that stages 1 and 2 may spend, half each.
    pub loss_budget: LossBudget,
    /// A seed held by more entries is not expanded in stage 1. `None` takes 1% of the shard's
    /// entries, at least `MIN_FREQUENT_SEED_ENTRIES`.
    pub frequent_seed_entries: Option<u64>,
    /// The most seed hits stage 1 holds; seeds over it count as frequent.
    pub max_seed_hits: u64,
    /// Bases added on each side of a candidate bin to make an alignment window.
    pub flank: u64,
    /// Worker threads; 0 uses the current thread pool.
    pub threads: usize,
    /// The settings of stage 3.
    pub align: AlignStageConfig,
}

/// The least frequent-seed limit that `SearchConfig::frequent_seed_entries` takes by default.
pub const MIN_FREQUENT_SEED_ENTRIES: u64 = 32;

impl Default for SearchConfig {
    fn default() -> Self {
        Self {
            level: SearchStage::Align,
            rule: DecisionRule::default(),
            loss_budget: LossBudget::default(),
            frequent_seed_entries: None,
            // 32 Mi pairs of entry and query k-mer, 256 MiB.
            max_seed_hits: 1 << 25,
            flank: 8_192,
            threads: 0,
            align: AlignStageConfig::default(),
        }
    }
}

impl SearchConfig {
    /// Checks that the settings make sense.
    pub fn validate(&self) -> Result<(), SearchError> {
        self.rule
            .validate()
            .map_err(|_| SearchError::InvalidConfig("the decision rule is invalid"))?;
        if self.flank == 0 {
            return Err(SearchError::InvalidConfig("the flank must be at least 1"));
        }
        self.align.validate()?;
        Ok(())
    }
}

/// The work of one search and how it ended.
#[derive(Clone, Debug, PartialEq)]
pub struct SearchSummary {
    pub query_length: u64,
    pub completion: Completion,
    pub sketch: SketchWork,
    pub locate: LocateWork,
    pub align: AlignStageWork,
}

/// The work of stage 3, summed over entries.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AlignStageWork {
    /// Located entries that stage 3 aligned.
    pub entries: u64,
    /// Entries that qualify.
    pub confirmed_entries: u64,
    /// Windows, bases read, anchors, chains and alignments.
    pub counts: AlignWork,
    /// Time spent in the stage, summed over threads.
    pub time: Duration,
}

impl AlignStageWork {
    /// Adds the work of another part of the search.
    pub fn add(&mut self, other: &Self) {
        self.entries += other.entries;
        self.confirmed_entries += other.confirmed_entries;
        self.counts.add(&other.counts);
        self.time += other.time;
    }
}

/// Searches one shard. The miss model depends only on k and w, so shards share it.
pub struct Searcher<'a> {
    shard: &'a Shard,
    model: &'a MissModel,
    seeder: Seeder,
}

impl<'a> Searcher<'a> {
    /// A searcher of `shard`, with a miss model of the shard's k and minimizer window.
    pub fn new(shard: &'a Shard, model: &'a MissModel) -> Result<Self, SearchError> {
        let parameters = shard.manifest().parameters;
        if model.k() != usize::from(parameters.k)
            || model.window() != usize::from(parameters.window)
        {
            return Err(SearchError::InvalidConfig(
                "the miss model has another k or minimizer window than the shard",
            ));
        }
        let seeder = Seeder::new(parameters.scheme())
            .map_err(|_| SearchError::InvalidConfig("the seed scheme of the shard"))?;
        Ok(Self {
            shard,
            model,
            seeder,
        })
    }

    /// Searches `query` and calls `on_event` each time an entry changes its state, from any
    /// worker thread. Bad input fails at once; a stop or a read failure ends the search early,
    /// which the summary's completion tells.
    pub fn search(
        &self,
        query: &[u8],
        circular: bool,
        config: &SearchConfig,
        control: SearchControl<'_>,
        on_event: &(dyn Fn(EntryEvent) + Sync),
    ) -> Result<SearchSummary, SearchError> {
        config.validate()?;
        if query.len() < usize::from(self.seeder.scheme().k) || query.len() >= u32::MAX as usize {
            return Err(SearchError::InvalidQuery(
                "the query must hold at least k bases and fewer than 2^32",
            ));
        }
        let budget = LossBudget::new(config.loss_budget.chance() / 2.0)?;
        let prepared = match config.level {
            SearchStage::Align => Some(PreparedQuery::new(query, circular, &config.align)?),
            _ => None,
        };
        let request = Request {
            config,
            prepared: prepared.as_ref(),
            control,
            on_event,
        };
        if config.threads == 0 {
            return Ok(self.run(query, circular, budget, &request));
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .build()
            .map_err(|_| SearchError::InvalidConfig("the thread pool cannot start"))?;
        Ok(pool.install(|| self.run(query, circular, budget, &request)))
    }

    fn run(
        &self,
        query: &[u8],
        circular: bool,
        budget: LossBudget,
        request: &Request<'_>,
    ) -> SearchSummary {
        let Request {
            config,
            control,
            on_event,
            ..
        } = *request;
        let query_kmers = QueryKmers::new(&self.seeder, query, circular);
        let length = query.len() as u64;
        let rule = &config.rule;
        let bound = CoverageBound::new(
            self.model,
            budget,
            rule.min_identity,
            length,
            circular,
            rule.min_query_coverage.required_bases(length),
        );
        let mut summary = SearchSummary {
            query_length: length,
            completion: Completion::Complete,
            sketch: SketchWork::default(),
            locate: LocateWork::default(),
            align: AlignStageWork::default(),
        };
        let sketch = match sketch_stage(
            self.shard.sketch(),
            &query_kmers,
            &bound,
            self.sketch_limits(config),
            control,
        ) {
            Ok(sketch) => sketch,
            Err(stopped) => {
                summary.completion = Completion::Stopped(stopped);
                return summary;
            }
        };
        summary.sketch = sketch.work;
        for candidate in &sketch.candidates {
            on_event(EntryEvent::Possible {
                entry: candidate.entry,
                coverage_bound: candidate.coverage_bound,
                containment_estimate: candidate.containment_estimate,
            });
        }
        if config.level == SearchStage::Sketch {
            return summary;
        }
        let locate = LocateRequest {
            query: &query_kmers,
            bound: &bound,
            frequent_seeds: &sketch.frequent_seeds,
            flank: config.flank,
        };
        let queue = EntryQueue::new(&sketch.candidates, control, on_event);
        queue.run(|possible, worker| self.finish_entry(possible, &locate, request, worker));
        (summary.completion, summary.locate, summary.align) = queue.finish();
        summary
    }

    fn sketch_limits(&self, config: &SearchConfig) -> SketchLimits {
        let entries = self.shard.manifest().entries.len() as u64;
        SketchLimits {
            frequent_seed_entries: config
                .frequent_seed_entries
                .unwrap_or((entries / 100).max(MIN_FREQUENT_SEED_ENTRIES)),
            max_seed_hits: config.max_seed_hits,
        }
    }

    /// Takes one candidate entry through the stages after stage 1 and reports its events.
    fn finish_entry(
        &self,
        possible: &PossibleEntry,
        locate: &LocateRequest<'_>,
        request: &Request<'_>,
        worker: &mut Worker,
    ) -> Result<(), (SearchStage, SearchError)> {
        let entry = possible.entry;
        let on_event = request.on_event;
        let outcome = locate_entry(
            self.shard,
            locate,
            possible,
            &mut worker.locate_workspace,
            &mut worker.locate,
            request.control,
        )
        .map_err(|error| (SearchStage::Locate, error))?;
        let located = match outcome {
            LocateOutcome::Located(located) => located,
            LocateOutcome::Dropped(reason) => {
                on_event(EntryEvent::Dropped {
                    entry,
                    stage: SearchStage::Locate,
                    reason,
                });
                return Ok(());
            }
        };
        let Some(prepared) = request.prepared else {
            on_event(EntryEvent::Located(located));
            return Ok(());
        };
        on_event(EntryEvent::Located(located.clone()));
        let started = Instant::now();
        let aligned = self.align(prepared, &located, request, worker);
        worker.align.entries += 1;
        worker.align.time += started.elapsed();
        let aligned = aligned.map_err(|error| (SearchStage::Align, error.into()))?;
        worker.align.counts.add(&aligned.work);
        if aligned.decision.qualifies {
            worker.align.confirmed_entries += 1;
            on_event(EntryEvent::Confirmed(aligned));
        } else {
            on_event(EntryEvent::Dropped {
                entry,
                stage: SearchStage::Align,
                reason: DropReason::Decision,
            });
        }
        Ok(())
    }

    /// Stage 3 for one located entry. The packed sequence is opened here, so it can later come
    /// from another source, such as a block cache.
    fn align(
        &self,
        prepared: &PreparedQuery,
        located: &LocatedEntry,
        request: &Request<'_>,
        worker: &mut Worker,
    ) -> Result<AlignedEntry, AlignStageError> {
        let packed = self.shard.packed(located.entry)?;
        align_packed_entry(
            &packed,
            prepared,
            located,
            &request.config.align,
            &request.config.rule,
            &mut worker.align_workspace,
            request.control,
        )
    }
}

/// What every entry of one search shares besides the stage 2 request.
#[derive(Clone, Copy)]
struct Request<'a> {
    config: &'a SearchConfig,
    /// The query prepared for stage 3, when the search goes that far.
    prepared: Option<&'a PreparedQuery>,
    control: SearchControl<'a>,
    on_event: &'a (dyn Fn(EntryEvent) + Sync),
}

/// What one worker thread keeps across entries: its buffers and its share of the work.
#[derive(Default)]
struct Worker {
    locate_workspace: LocateWorkspace,
    align_workspace: AlignStageWorkspace,
    locate: LocateWork,
    align: AlignStageWork,
}

/// The candidate entries of one search as a work queue: workers take the next entry in order
/// until the queue is empty. After a stop or a failure the rest are reported as unfinished.
struct EntryQueue<'a> {
    candidates: &'a [PossibleEntry],
    control: SearchControl<'a>,
    on_event: &'a (dyn Fn(EntryEvent) + Sync),
    next: AtomicUsize,
    halted: AtomicBool,
    completion: Mutex<Completion>,
    work: Mutex<(LocateWork, AlignStageWork)>,
}

impl<'a> EntryQueue<'a> {
    fn new(
        candidates: &'a [PossibleEntry],
        control: SearchControl<'a>,
        on_event: &'a (dyn Fn(EntryEvent) + Sync),
    ) -> Self {
        Self {
            candidates,
            control,
            on_event,
            next: AtomicUsize::new(0),
            halted: AtomicBool::new(false),
            completion: Mutex::new(Completion::Complete),
            work: Mutex::default(),
        }
    }

    /// Runs `finish` for every entry on the workers of the current thread pool.
    fn run<F>(&self, finish: F)
    where
        F: Fn(&PossibleEntry, &mut Worker) -> Result<(), (SearchStage, SearchError)> + Sync,
    {
        let workers = rayon::current_num_threads().min(self.candidates.len());
        rayon::scope(|scope| {
            for _ in 0..workers {
                scope.spawn(|_| self.work_on(&finish));
            }
        });
    }

    fn work_on<F>(&self, finish: &F)
    where
        F: Fn(&PossibleEntry, &mut Worker) -> Result<(), (SearchStage, SearchError)>,
    {
        let mut worker = Worker::default();
        while let Some(possible) = self
            .candidates
            .get(self.next.fetch_add(1, Ordering::Relaxed))
        {
            if self.halted.load(Ordering::Relaxed) {
                self.unfinished(possible.entry, SearchStage::Locate);
                continue;
            }
            let result = match self.control.check() {
                Err(stopped) => Err((SearchStage::Locate, SearchError::Stopped(stopped))),
                Ok(()) => finish(possible, &mut worker),
            };
            if let Err((stage, error)) = result {
                self.halt(error);
                self.unfinished(possible.entry, stage);
            }
        }
        let mut work = lock(&self.work);
        work.0.add(&worker.locate);
        work.1.add(&worker.align);
    }

    /// Stops the queue and keeps the first reason.
    fn halt(&self, error: SearchError) {
        self.halted.store(true, Ordering::Relaxed);
        let mut completion = lock(&self.completion);
        if *completion == Completion::Complete {
            *completion = match error {
                SearchError::Stopped(stopped)
                | SearchError::Align(AlignStageError::Stopped(stopped)) => {
                    Completion::Stopped(stopped)
                }
                error => Completion::Failed(error.to_string()),
            };
        }
    }

    fn unfinished(&self, entry: u32, stage: SearchStage) {
        (self.on_event)(EntryEvent::Unfinished { entry, stage });
    }

    fn finish(self) -> (Completion, LocateWork, AlignStageWork) {
        let completion = self
            .completion
            .into_inner()
            .unwrap_or_else(|e| e.into_inner());
        let (locate, align) = self.work.into_inner().unwrap_or_else(|e| e.into_inner());
        (completion, locate, align)
    }
}

fn lock<T>(mutex: &Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    mutex
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}
