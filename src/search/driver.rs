//! The driver of one search against one shard: it connects the stages, streams the entry events,
//! and sums the work.
//!
//! Stage 1 runs in parallel chunks. The candidate entries then go to a work queue, highest
//! coverage bound first, and each worker takes one entry through stages 2 and 3, so the events of
//! one entry come from one thread, in order. Memory per request stays bounded: stage 1 holds at
//! most `max_seed_hits` seed hits, and stages 2 and 3 hold one entry per worker; results leave as
//! events.

use super::bound::CoverageBound;
use super::control::SearchControl;
use super::error::SearchError;
use super::events::{Completion, EntryEvent, SearchStage};
use super::locate_stage::{
    LocateOutcome, LocateRequest, LocateWork, LocateWorkspace, locate_entry,
};
use super::query::QueryKmers;
use super::sketch_stage::{PossibleEntry, SketchLimits, SketchWork, sketch_stage};
use crate::decide::DecisionRule;
use crate::index::Shard;
use crate::loss::{LossBudget, MissModel};
use crate::seeds::Seeder;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

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
        if self.level == SearchStage::Align {
            return Err(SearchError::InvalidConfig("stage 3 is not built yet"));
        }
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
        if config.threads == 0 {
            return Ok(self.run(query, circular, config, budget, control, on_event));
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .build()
            .map_err(|_| SearchError::InvalidConfig("the thread pool cannot start"))?;
        Ok(pool.install(|| self.run(query, circular, config, budget, control, on_event)))
    }

    fn run(
        &self,
        query: &[u8],
        circular: bool,
        config: &SearchConfig,
        budget: LossBudget,
        control: SearchControl<'_>,
        on_event: &(dyn Fn(EntryEvent) + Sync),
    ) -> SearchSummary {
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
        let request = LocateRequest {
            query: &query_kmers,
            bound: &bound,
            frequent_seeds: &sketch.frequent_seeds,
            flank: config.flank,
        };
        let queue = EntryQueue::new(&sketch.candidates, control, on_event);
        queue.run(|possible, workspace, work| {
            self.finish_entry(possible, &request, workspace, work, control, on_event)
        });
        (summary.completion, summary.locate) = queue.finish();
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
        request: &LocateRequest<'_>,
        workspace: &mut LocateWorkspace,
        work: &mut LocateWork,
        control: SearchControl<'_>,
        on_event: &(dyn Fn(EntryEvent) + Sync),
    ) -> Result<(), (SearchStage, SearchError)> {
        let entry = possible.entry;
        let outcome = locate_entry(self.shard, request, possible, workspace, work, control)
            .map_err(|error| (SearchStage::Locate, error))?;
        on_event(match outcome {
            LocateOutcome::Located(located) => EntryEvent::Located(located),
            LocateOutcome::Dropped(reason) => EntryEvent::Dropped {
                entry,
                stage: SearchStage::Locate,
                reason,
            },
        });
        Ok(())
    }
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
    work: Mutex<LocateWork>,
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
            work: Mutex::new(LocateWork::default()),
        }
    }

    /// Runs `finish` for every entry on the workers of the current thread pool.
    fn run<F>(&self, finish: F)
    where
        F: Fn(
                &PossibleEntry,
                &mut LocateWorkspace,
                &mut LocateWork,
            ) -> Result<(), (SearchStage, SearchError)>
            + Sync,
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
        F: Fn(
            &PossibleEntry,
            &mut LocateWorkspace,
            &mut LocateWork,
        ) -> Result<(), (SearchStage, SearchError)>,
    {
        let mut workspace = LocateWorkspace::default();
        let mut work = LocateWork::default();
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
                Ok(()) => finish(possible, &mut workspace, &mut work),
            };
            if let Err((stage, error)) = result {
                self.halt(error);
                self.unfinished(possible.entry, stage);
            }
        }
        lock(&self.work).add(&work);
    }

    /// Stops the queue and keeps the first reason.
    fn halt(&self, error: SearchError) {
        self.halted.store(true, Ordering::Relaxed);
        let mut completion = lock(&self.completion);
        if *completion == Completion::Complete {
            *completion = match error {
                SearchError::Stopped(stopped) => Completion::Stopped(stopped),
                error => Completion::Failed(error.to_string()),
            };
        }
    }

    fn unfinished(&self, entry: u32, stage: SearchStage) {
        (self.on_event)(EntryEvent::Unfinished { entry, stage });
    }

    fn finish(self) -> (Completion, LocateWork) {
        let completion = self
            .completion
            .into_inner()
            .unwrap_or_else(|e| e.into_inner());
        let work = self.work.into_inner().unwrap_or_else(|e| e.into_inner());
        (completion, work)
    }
}

fn lock<T>(mutex: &Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    mutex
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}
