//! The driver of one search against one shard: it connects the stages, streams the entry events,
//! and sums the work.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::align_stage::{
    AlignStageConfig, AlignStageWorkspace, AlignWork, PreparedQuery, align_packed_entry,
};
use super::bound::CoverageBound;
use super::control::{SearchControl, Stopped};
use super::error::SearchError;
use super::events::{Completion, DropReason, EntryEvent, SearchStage};
use super::locate_stage::{
    LocateOutcome, LocateRequest, LocateWork, LocateWorkspace, locate_entry,
};
use super::located::LocatedEntry;
use super::query::QueryKmers;
use super::sketch_stage::{PossibleEntry, SketchLimits, SketchWork, sketch_stage};
use super::times::{PartStart, SearchTimes};
use crate::decide::DecisionRule;
use crate::index::Shard;
use crate::loss::{LossBudget, MissModel};
use crate::seeds::Seeder;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::mpsc::SyncSender;
use std::time::Instant;

/// The settings of one search.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SearchConfig {
    /// The last stage that runs.
    pub last_stage: SearchStage,
    /// The rule that reports an entry; its query coverage and identity also drive the pruning.
    pub rule: DecisionRule,
    /// The chance of missing a qualifying entry that stages 1 and 2 may spend, half each.
    pub loss_budget: LossBudget,
    /// A seed held by more entries is not expanded in stage 1. `None` takes 1% of the shard's
    /// entries, at least `MIN_FREQUENT_LIMIT`.
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
pub const MIN_FREQUENT_LIMIT: u64 = 32;

impl Default for SearchConfig {
    fn default() -> Self {
        Self {
            last_stage: SearchStage::Align,
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
        self.rule.validate()?;
        if self.flank == 0 {
            return Err(SearchError::InvalidConfig("the flank must be at least 1"));
        }
        self.align.validate()?;
        Ok(())
    }
}

/// The work of one search and how it ended.
#[derive(Debug)]
pub struct SearchSummary {
    pub query_length: u64,
    pub completion: Completion,
    pub sketch: SketchWork,
    pub locate: LocateWork,
    pub align: AlignStageWork,
    /// Wall seconds, CPU seconds and major page faults of the parts of the search.
    pub times: SearchTimes,
}

/// The work of stage 3, summed over entries.
#[derive(Clone, Copy, Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct AlignStageWork {
    /// Located entries that stage 3 aligned.
    pub entries: u64,
    /// Entries that qualify.
    pub confirmed_entries: u64,
    /// Windows, bases read, anchors, chains and alignments.
    pub counts: AlignWork,
    /// Seconds spent in the stage, summed over threads.
    pub seconds: f64,
}

impl AlignStageWork {
    /// Adds the work of another part of the search.
    pub fn add(&mut self, other: &Self) {
        self.entries += other.entries;
        self.confirmed_entries += other.confirmed_entries;
        self.counts.add(&other.counts);
        self.seconds += other.seconds;
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
        if model.k() != usize::from(parameters.scheme.k)
            || model.w() != usize::from(parameters.scheme.w)
        {
            return Err(SearchError::InvalidConfig(
                "the miss model has another k or minimizer window than the shard",
            ));
        }
        let seeder = Seeder::new(parameters.scheme)?;
        Ok(Self {
            shard,
            model,
            seeder,
        })
    }

    /// Searches `query` and sends an event into `events` each time an entry changes its state;
    /// another thread must read the receiver, and dropping it stops the search as if asked.
    /// Bad input fails at once; a later stop or failure ends it early, as the summary tells.
    pub fn search(
        &self,
        query: &[u8],
        circular: bool,
        config: &SearchConfig,
        control: SearchControl<'_>,
        events: &SyncSender<EntryEvent>,
    ) -> Result<SearchSummary, SearchError> {
        config.validate()?;
        if query.len() < usize::from(self.seeder.scheme().k) || query.len() >= u32::MAX as usize {
            return Err(SearchError::InvalidQuery(
                "the query must hold at least k bases and fewer than 2^32",
            ));
        }
        let budget = LossBudget::new(config.loss_budget.chance() / 2.0)?;
        // The query tables are built on the threads of the search too.
        let search = || {
            let start = PartStart::now();
            let prepared = match config.last_stage {
                SearchStage::Align => Some(PreparedQuery::new(query, circular, &config.align)?),
                _ => None,
            };
            let workspace = AlignStageWorkspace::default();
            let request = Request {
                config,
                prepared: prepared.as_ref(),
                workspace: &workspace,
                control,
                events,
            };
            Ok(self.run(query, circular, budget, &request, start))
        };
        if config.threads == 0 {
            return search();
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .build()?;
        pool.install(search)
    }

    /// Stage 1 runs in parallel chunks; then each worker takes the next candidate entry, highest
    /// coverage bound first, through stages 2 and 3, so the events of an entry come in order.
    fn run(
        &self,
        query: &[u8],
        circular: bool,
        budget: LossBudget,
        request: &Request<'_>,
        start: PartStart,
    ) -> SearchSummary {
        let Request {
            config,
            control,
            events,
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
            times: SearchTimes {
                open_shard_seconds: self.shard.open_seconds(),
                ..SearchTimes::default()
            },
        };
        summary.times.prepare_query = start.finish();
        let start = PartStart::now();
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
        summary.times.sketch = start.finish();
        for candidate in &sketch.candidates {
            let possible = EntryEvent::Possible {
                entry: candidate.entry,
                coverage_bound: candidate.coverage_bound,
                containment_estimate: candidate.containment_estimate,
            };
            if events.send(possible).is_err() {
                summary.completion = Completion::Stopped(Stopped::Asked);
                return summary;
            }
        }
        if config.last_stage == SearchStage::Sketch {
            return summary;
        }
        let locate = LocateRequest {
            query: &query_kmers,
            bound: &bound,
            frequent_seeds: &sketch.frequent_seeds,
            flank: config.flank,
        };
        let queue = EntryQueue::new(&sketch.candidates, control, events);
        let finish = |possible: &PossibleEntry, worker: &mut Worker| {
            self.finish_entry(possible, &locate, request, worker)
        };
        let start = PartStart::now();
        queue.run(finish, &mut summary);
        summary.times.entries = start.finish();
        summary
    }

    fn sketch_limits(&self, config: &SearchConfig) -> SketchLimits {
        let entries = self.shard.manifest().entries.len() as u64;
        SketchLimits {
            frequent_seed_entries: config
                .frequent_seed_entries
                .unwrap_or((entries / 100).max(MIN_FREQUENT_LIMIT)),
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
    ) -> Result<(), Halt> {
        let entry = possible.entry;
        let events = request.events;
        let outcome = locate_entry(
            self.shard,
            locate,
            possible,
            &mut worker.locate_workspace,
            &mut worker.locate,
            request.control,
        )
        .map_err(|error| Halt {
            stage: SearchStage::Locate,
            error,
        })?;
        let located = match outcome {
            LocateOutcome::Located(located) => located,
            LocateOutcome::Dropped(reason) => {
                let dropped = EntryEvent::Dropped {
                    entry,
                    stage: SearchStage::Locate,
                    reason,
                };
                return emit(events, dropped, SearchStage::Locate);
            }
        };
        let Some(prepared) = request.prepared else {
            return emit(events, EntryEvent::Located(located), SearchStage::Locate);
        };
        emit(
            events,
            EntryEvent::Located(located.clone()),
            SearchStage::Locate,
        )?;
        let started = Instant::now();
        let last = self.align(prepared, &located, request, worker);
        worker.align.entries += 1;
        worker.align.seconds += started.elapsed().as_secs_f64();
        let last = last.map_err(|error| Halt {
            stage: SearchStage::Align,
            error,
        })?;
        emit(events, last, SearchStage::Align)
    }

    /// Stage 3 for one located entry: opens its packed sequence once, aligns and decides, and
    /// returns the last event of the entry, confirmed (with the contig table) or dropped.
    fn align(
        &self,
        prepared: &PreparedQuery,
        located: &LocatedEntry,
        request: &Request<'_>,
        worker: &mut Worker,
    ) -> Result<EntryEvent, SearchError> {
        let packed = self.shard.packed(located.entry)?;
        for windows in &located.alignment_windows {
            for range in &windows.ranges {
                packed.will_need(windows.contig as usize, range.clone());
            }
        }
        let aligned = align_packed_entry(
            &packed,
            prepared,
            located,
            &request.config.align,
            &request.config.rule,
            request.workspace,
            request.control,
        )?;
        worker.align.counts.add(&aligned.work);
        if !aligned.decision.qualifies {
            return Ok(EntryEvent::Dropped {
                entry: located.entry,
                stage: SearchStage::Align,
                reason: DropReason::Decision,
            });
        }
        worker.align.confirmed_entries += 1;
        let mosaic = &aligned.decision.mosaic;
        let primary = mosaic.primary.iter().map(|selected| &selected.fragment);
        let fragments = primary.chain(&mosaic.alternatives);
        let contigs = packed.contigs(fragments.map(|fragment| fragment.contig_index))?;
        Ok(EntryEvent::Confirmed { aligned, contigs })
    }
}

/// What every entry of one search shares besides the stage 2 request.
#[derive(Clone, Copy)]
struct Request<'a> {
    config: &'a SearchConfig,
    /// The query prepared for stage 3, when the search goes that far.
    prepared: Option<&'a PreparedQuery>,
    /// The buffers of stage 3, lent to the threads that align windows.
    workspace: &'a AlignStageWorkspace,
    control: SearchControl<'a>,
    events: &'a SyncSender<EntryEvent>,
}

/// What one worker thread keeps across entries: its buffers and its share of the work.
#[derive(Default)]
struct Worker {
    locate_workspace: LocateWorkspace,
    locate: LocateWork,
    align: AlignStageWork,
}

/// Why a worker stopped taking entries: the error, and the stage the entry was in.
struct Halt {
    stage: SearchStage,
    error: SearchError,
}

/// What a worker hands back once the queue is empty: its work, and the error that halted it.
struct WorkerEnd {
    worker: Worker,
    result: Result<(), SearchError>,
}

/// The candidate entries of one search as a work queue: workers take the next entry in order
/// until the queue is empty. After a stop or a failure the rest are reported as unfinished.
struct EntryQueue<'a> {
    candidates: &'a [PossibleEntry],
    control: SearchControl<'a>,
    events: &'a SyncSender<EntryEvent>,
    next: AtomicUsize,
    /// Set by the first worker that stops or fails, so the others take no more entries.
    halted: AtomicBool,
}

impl<'a> EntryQueue<'a> {
    fn new(
        candidates: &'a [PossibleEntry],
        control: SearchControl<'a>,
        events: &'a SyncSender<EntryEvent>,
    ) -> Self {
        Self {
            candidates,
            control,
            events,
            next: AtomicUsize::new(0),
            halted: AtomicBool::new(false),
        }
    }

    /// Runs `finish` for every entry on the workers of the current thread pool, then adds their
    /// work to `summary`. The first error, by worker, decides how the search ended.
    fn run<F>(&self, finish: F, summary: &mut SearchSummary)
    where
        F: Fn(&PossibleEntry, &mut Worker) -> Result<(), Halt> + Sync,
    {
        let workers = rayon::current_num_threads().min(self.candidates.len());
        let ends: Vec<WorkerEnd> = (0..workers)
            .into_par_iter()
            .map(|_| self.work_on(&finish))
            .collect();
        for end in ends {
            summary.locate.add(&end.worker.locate);
            summary.align.add(&end.worker.align);
            if let (Completion::Complete, Err(error)) = (&summary.completion, end.result) {
                summary.completion = error.into();
            }
        }
    }

    fn work_on<F>(&self, finish: &F) -> WorkerEnd
    where
        F: Fn(&PossibleEntry, &mut Worker) -> Result<(), Halt>,
    {
        let mut worker = Worker::default();
        while let Some(possible) = self.take() {
            if self.halted.load(Ordering::Relaxed) {
                self.unfinished(possible.entry, SearchStage::Locate);
                continue;
            }
            let result = match self.control.check() {
                Err(stopped) => Err(Halt {
                    stage: SearchStage::Locate,
                    error: SearchError::Stopped(stopped),
                }),
                Ok(()) => finish(possible, &mut worker),
            };
            if let Err(Halt { stage, error }) = result {
                self.halted.store(true, Ordering::Relaxed);
                self.unfinished(possible.entry, stage);
                while let Some(rest) = self.take() {
                    self.unfinished(rest.entry, SearchStage::Locate);
                }
                return WorkerEnd {
                    worker,
                    result: Err(error),
                };
            }
        }
        WorkerEnd {
            worker,
            result: Ok(()),
        }
    }

    fn take(&self) -> Option<&'a PossibleEntry> {
        self.candidates
            .get(self.next.fetch_add(1, Ordering::Relaxed))
    }

    fn unfinished(&self, entry: u32, stage: SearchStage) {
        // A dropped receiver means nobody listens any more.
        let _ = self.events.send(EntryEvent::Unfinished { entry, stage });
    }
}

/// Sends one event; a dropped receiver stops the search as if asked.
fn emit(
    events: &SyncSender<EntryEvent>,
    event: EntryEvent,
    stage: SearchStage,
) -> Result<(), Halt> {
    events.send(event).map_err(|_| Halt {
        stage,
        error: SearchError::Stopped(Stopped::Asked),
    })
}
