//! The driver of one search against one shard: it connects the stages, streams the entry events,
//! and sums the work.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::align_stage::{
    AlignStageConfig, AlignStageError, AlignStageWorkspace, AlignWork, PreparedQuery,
    align_packed_entry,
};
use super::bound::CoverageBound;
use super::control::{SearchControl, Stopped};
use super::error::SearchError;
use super::events::{Completion, DropReason, EntryEvent, SearchStage};
use super::locate_stage::{
    LocateOutcome, LocateRequest, LocateWork, LocateWorkspace, locate_entry, read_entry_locations,
    will_need_entry,
};
use super::located::LocatedEntry;
use super::query::{QueryKmers, lookup_bound};
use super::seed_mismatches::{LookupPlan, lookup_plan, target_length};
use super::sketch_stage::{
    PossibleEntry, SketchLimits, SketchWork, held_keys, keep_found, sketch_stage,
};
use super::times::{PartStart, PartTime, SearchTimes};
use crate::align::{AlignmentError, Interval};
use crate::decide::DecisionRule;
use crate::index::{SeedLocation, Shard, seed_key};
use crate::loss::{MaxMiss, MissModel};
use crate::seeds::{Seeder, low_complexity_mask};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::Arc;
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
    /// The largest total chance of missing a qualifying entry: the seed miss bound, plus what the
    /// coverage bound may spend. Stages 1 and 2 share one fixed allowance and stage 2 sees a subset
    /// of the hits, so one share covers both.
    pub max_miss: MaxMiss,
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

/// The sketch bytes per lookup above which a shared lookup pass pays. A query's own pass reads up
/// to its whole sketch part, at about 0.3 ns per byte on CephFS, and the shared pass makes its
/// lookups a second time, at about 70 ns each (measured on worm): 230 bytes, rounded up.
const SKETCH_BYTES_PER_LOOKUP: u64 = 256;
/// Candidate entries ahead of the workers whose reads are asked for, so on S3 many requests are
/// in flight; their blocks take a few MiB of the block cache.
const ENTRIES_AHEAD: usize = 128;

/// The least frequent-seed limit that `SearchConfig::frequent_seed_entries` takes by default.
pub const MIN_FREQUENT_LIMIT: u64 = 32;

impl Default for SearchConfig {
    fn default() -> Self {
        Self {
            last_stage: SearchStage::Align,
            rule: DecisionRule::default(),
            max_miss: MaxMiss::default(),
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
    /// The lookup level: 0 exact, 1 one mismatch, 2 two mismatches, 3 also one indel.
    pub seed_mismatches: u8,
    pub completion: Completion,
    pub sketch: SketchWork,
    pub locate: LocateWork,
    pub align: AlignStageWork,
    /// Wall seconds, CPU seconds and major page faults of the parts of the search.
    pub times: SearchTimes,
}

/// The seed keys held by the sketch, found by `Searcher::prepare_kmers` for a batch of queries,
/// and one query's share of the time of the pass.
#[derive(Debug)]
pub struct PreparedKmers {
    held: Arc<Vec<u64>>,
    time: PartTime,
}

/// The work of stage 3, summed over entries.
#[derive(Clone, Copy, Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct AlignStageWork {
    /// Located entries that stage 3 aligned.
    pub entries: u64,
    /// Entries that qualify.
    pub confirmed_entries: u64,
    /// Entries dropped for a resource limit, and so not decided.
    pub limited_entries: u64,
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
        self.limited_entries += other.limited_entries;
        self.counts.add(&other.counts);
        self.seconds += other.seconds;
    }
}

/// Searches one shard.
pub struct Searcher<'a> {
    shard: &'a Shard,
    /// The miss models of the shard's k and w, for seeds hit with 0, 1 and 2 mismatches.
    models: [MissModel; 3],
    seeder: Seeder,
}

impl<'a> Searcher<'a> {
    /// A searcher of `shard`.
    pub fn new(shard: &'a Shard) -> Result<Self, SearchError> {
        let scheme = shard.manifest().parameters.scheme;
        let seeder = Seeder::new(scheme)?;
        let model = |mismatches| {
            MissModel::new(scheme, mismatches)
                .map_err(|_| SearchError::InvalidConfig("the seed scheme of the shard is invalid"))
        };
        Ok(Self {
            shard,
            models: [model(0)?, model(1)?, model(2)?],
            seeder,
        })
    }

    /// The smallest query length from which every longer query reaches `--max-miss` with the
    /// lookups the rule chooses, or `None` if none does.
    pub fn target_length(&self, config: &SearchConfig) -> Option<u64> {
        target_length(&self.models, &config.rule, config.max_miss)
    }

    /// The most lookups that the search of a query of `length` bases makes, by its plan.
    pub fn lookup_bound(&self, length: u64, circular: bool, config: &SearchConfig) -> u64 {
        let plan = lookup_plan(
            &self.models,
            &config.rule,
            config.max_miss,
            length,
            circular,
        );
        lookup_bound(&plan, length)
    }

    /// Whether a query of `length` bases gains from a shared lookup pass: when its own pass could
    /// read more sketch bytes than making its lookups again costs, at `SKETCH_BYTES_PER_LOOKUP`.
    pub fn shares_lookups(&self, length: u64, circular: bool, config: &SearchConfig) -> bool {
        let sketch = self.shard.seed_file().sketch_part();
        let bound = self.lookup_bound(length, circular, config);
        sketch.end - sketch.start >= bound.saturating_mul(SKETCH_BYTES_PER_LOOKUP)
    }

    /// The seed keys that the sketch holds among the lookups of several queries, found in one
    /// pass: the keys of all of them are checked together, so close keys of different queries
    /// share bucket reads. Only the distinct keys of each query are held during the pass. Each
    /// search then makes its lookups again, by its own plan, and keeps those whose keys are held,
    /// so it is the same as without the pass. Each query gets the share of the time of the pass
    /// that its lookups make up. A query the search would reject gets `None`, and so do all when
    /// the pass fails; their searches then look up their own keys.
    pub fn prepare_kmers(
        &self,
        queries: &[(&[u8], bool)],
        config: &SearchConfig,
    ) -> Vec<Option<PreparedKmers>> {
        let pass = || {
            let start = PartStart::now();
            let k = usize::from(self.seeder.scheme().k);
            let (mut keys, mut counts) = (Vec::new(), Vec::with_capacity(queries.len()));
            for &(query, circular) in queries {
                if query.len() < k || query.len() >= u32::MAX as usize {
                    counts.push(None);
                    continue;
                }
                let length = query.len() as u64;
                let rule = &config.rule;
                let plan = lookup_plan(&self.models, rule, config.max_miss, length, circular);
                let (count, mut own) = QueryKmers::made_keys(&self.seeder, query, circular, &plan);
                counts.push(Some(count));
                keys.append(&mut own);
            }
            let held = Arc::new(held_keys(self.shard.seed_file(), keys)?);
            let time = start.finish();
            let total: u64 = counts.iter().flatten().sum();
            Ok::<_, SearchError>(
                (counts.into_iter())
                    .map(|count| {
                        count.map(|count| PreparedKmers {
                            held: Arc::clone(&held),
                            time: time.share(count as f64 / total.max(1) as f64),
                        })
                    })
                    .collect(),
            )
        };
        let prepared = match config.threads {
            0 => pass(),
            threads => rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(SearchError::from)
                .and_then(|pool| pool.install(pass)),
        };
        prepared.unwrap_or_else(|_| queries.iter().map(|_| None).collect())
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
        self.search_prepared(query, circular, config, control, events, None)
    }

    /// As `search`, with the held seed keys from `prepare_kmers` when given. The query's share
    /// of the time of the pass counts toward its preparation.
    pub fn search_prepared(
        &self,
        query: &[u8],
        circular: bool,
        config: &SearchConfig,
        control: SearchControl<'_>,
        events: &SyncSender<EntryEvent>,
        prepared: Option<PreparedKmers>,
    ) -> Result<SearchSummary, SearchError> {
        config.validate()?;
        if config.last_stage != SearchStage::Sketch && !self.shard.seed_file().has_locators() {
            return Err(SearchError::InvalidConfig(
                "the shard was built with the sketch part only: search it with --level sketch",
            ));
        }
        if query.len() < usize::from(self.seeder.scheme().k) || query.len() >= u32::MAX as usize {
            return Err(SearchError::InvalidQuery(
                "the query must hold at least k bases and fewer than 2^32",
            ));
        }
        let budget = config.max_miss;
        let plan = lookup_plan(
            &self.models,
            &config.rule,
            budget,
            query.len() as u64,
            circular,
        );
        let masked: Vec<Interval> = (low_complexity_mask(query).into_iter())
            .map(|stretch| Interval {
                start: stretch.start as u64,
                end: stretch.end as u64,
            })
            .collect();
        let low_complexity = Interval::total_length(&masked);
        // A search that waited past its deadline, or was stopped meanwhile, builds no tables.
        if let Err(stopped) = control.check() {
            let mut summary = self.empty_summary(query.len() as u64, plan.level());
            summary.completion = Completion::Stopped(stopped);
            return Ok(summary);
        }
        // The query tables are built on the threads of the search too.
        let (held, prepared_time) = match prepared {
            Some(prepared) => (Some(prepared.held), prepared.time),
            None => (None, PartTime::default()),
        };
        let search = || {
            let start = PartStart::now();
            let (prepared, query_kmers) = rayon::join(
                || match config.last_stage {
                    SearchStage::Align => {
                        let prepared = PreparedQuery::new(query, circular, &config.align);
                        prepared.map(|mut prepared| {
                            prepared.low_complexity = masked.clone();
                            Some(prepared)
                        })
                    }
                    _ => Ok(None),
                },
                || {
                    let seeds = self.shard.seed_file();
                    QueryKmers::filtered(&self.seeder, query, circular, &plan, |lookups| {
                        control.check()?;
                        match &held {
                            Some(held) => {
                                let held = |hash| held.binary_search(&seed_key(hash)).is_ok();
                                lookups.retain(|&(hash, _)| held(hash));
                                Ok(())
                            }
                            None => keep_found(seeds, lookups),
                        }
                    })
                },
            );
            let prepared = prepared?;
            let query_kmers = match query_kmers {
                Ok(query_kmers) => query_kmers,
                Err(error) => {
                    let mut summary = self.empty_summary(query.len() as u64, plan.level());
                    summary.completion = error.into();
                    return Ok(summary);
                }
            };
            let workspace = AlignStageWorkspace::default();
            let request = Request {
                config,
                prepared: prepared.as_ref(),
                workspace: &workspace,
                control,
                events,
                low_complexity,
            };
            let mut summary = self.run(&query_kmers, circular, &plan, budget, &request, start);
            summary.seed_mismatches = plan.level() as u8;
            summary.times.prepare_query = summary.times.prepare_query.plus(&prepared_time);
            Ok(summary)
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
        query_kmers: &QueryKmers,
        circular: bool,
        plan: &LookupPlan,
        budget: MaxMiss,
        request: &Request<'_>,
        start: PartStart,
    ) -> SearchSummary {
        let Request {
            config,
            control,
            events,
            low_complexity,
            ..
        } = *request;
        let length = query_kmers.length();
        let rule = &config.rule;
        let model = &self.models[plan.full];
        let bound = CoverageBound::new(
            model,
            budget,
            rule.min_identity,
            length,
            circular,
            rule.min_query_coverage.required_bases(length),
        );
        let mut summary = self.empty_summary(length, model.mismatches());
        summary.times.prepare_query = start.finish();
        let (miss_upper_bound, total_miss) = (bound.seed_miss_bound(), bound.total_miss());
        let miss_over_budget = !bound.prunes();
        summary.sketch = SketchWork {
            low_complexity_bases: low_complexity,
            sampled_kmers: plan.sampled,
            miss_upper_bound,
            total_miss,
            miss_over_budget,
            ..SketchWork::default()
        };
        // Half of the required coverage must be normal bases, so a query with fewer normal bases
        // has no answer, and nothing is looked up.
        if length - low_complexity < rule.min_query_coverage.required_normal_bases(length) {
            return summary;
        }
        let start = PartStart::now();
        let sketch = match sketch_stage(
            self.shard.seed_file(),
            query_kmers,
            &bound,
            self.sketch_limits(config),
            control,
        ) {
            Ok(sketch) => sketch,
            Err(error) => {
                summary.completion = error.into();
                return summary;
            }
        };
        summary.sketch = SketchWork {
            low_complexity_bases: low_complexity,
            sampled_kmers: plan.sampled,
            miss_upper_bound,
            total_miss,
            miss_over_budget,
            ..sketch.work
        };
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
            query: query_kmers,
            bound: &bound,
            frequent_seeds: &sketch.frequent_seeds,
            flank: config.flank,
            verify: Some(&self.seeder),
        };
        let start = PartStart::now();
        let locations = read_entry_locations(
            self.shard.seed_file(),
            &locate,
            &sketch.found,
            &sketch.candidates,
            control,
            &mut summary.locate,
        );
        let queue = EntryQueue::new(&sketch.candidates, control, events);
        let locations = match locations {
            Ok(locations) => locations,
            Err(error) => {
                queue.give_up(SearchStage::Locate);
                summary.completion = error.into();
                summary.times.entries = start.finish();
                return summary;
            }
        };
        let finish = |index: usize, possible: &PossibleEntry, worker: &mut Worker| {
            self.finish_entry(possible, &locations[index], &locate, request, worker)
        };
        let reach = config.flank.max(u64::from(self.seeder.scheme().k));
        let prefetch = |index: usize, possible: &PossibleEntry| {
            will_need_entry(self.shard, possible.entry, &locations[index], reach);
        };
        queue.run(finish, prefetch, &mut summary);
        summary.times.entries = start.finish();
        summary
    }

    /// The summary of a complete search of a query of `length` bases that did no work yet.
    fn empty_summary(&self, length: u64, mismatches: usize) -> SearchSummary {
        SearchSummary {
            query_length: length,
            seed_mismatches: mismatches as u8,
            completion: Completion::Complete,
            sketch: SketchWork::default(),
            locate: LocateWork::default(),
            align: AlignStageWork::default(),
            times: SearchTimes {
                open_shard_seconds: self.shard.take_open_seconds(),
                ..SearchTimes::default()
            },
        }
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
        locations: &[(u32, SeedLocation)],
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
            locations,
            &mut worker.locate_workspace,
            &mut worker.locate,
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
            let windows = located.alignment_windows.iter().map(|w| w.contig);
            let contigs = (self.shard.packed(entry))
                .and_then(|packed| packed.contigs(windows))
                .map_err(|error| Halt {
                    stage: SearchStage::Locate,
                    error: error.into(),
                })?;
            let event = EntryEvent::Located(located, contigs);
            return emit(events, event, SearchStage::Locate);
        };
        let event = EntryEvent::Located(located.clone(), BTreeMap::new());
        emit(events, event, SearchStage::Locate)?;
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
        let windows = located.alignment_windows.iter();
        packed.will_need(windows.flat_map(|w| w.ranges.iter().map(|r| (w.contig as usize, r))));
        let aligned = align_packed_entry(
            &packed,
            prepared,
            located,
            &request.config.align,
            &request.config.rule,
            request.workspace,
            request.control,
        );
        let dropped = |reason| EntryEvent::Dropped {
            entry: located.entry,
            stage: SearchStage::Align,
            reason,
        };
        let aligned = match aligned {
            Ok(aligned) => aligned,
            Err(error) if is_resource_limit(&error) => {
                worker.align.limited_entries += 1;
                return Ok(dropped(DropReason::ResourceLimit));
            }
            Err(error) => return Err(error.into()),
        };
        worker.align.counts.add(&aligned.work);
        if !aligned.decision.qualifies {
            return Ok(dropped(DropReason::Decision));
        }
        worker.align.confirmed_entries += 1;
        let mosaic = &aligned.decision.mosaic;
        let primary = mosaic.primary.iter().map(|selected| &selected.fragment);
        let fragments = primary.chain(&mosaic.alternatives);
        let contigs = packed.contigs(fragments.map(|fragment| fragment.contig_index))?;
        Ok(EntryEvent::Confirmed {
            aligned: Box::new(aligned),
            contigs,
        })
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
    /// Query bases in low-complexity stretches.
    low_complexity: u64,
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
    /// work to `summary`. The first error, by worker, decides how the search ended. Taking an
    /// entry calls `prefetch` (with its place) for the entry `ENTRIES_AHEAD` later, so many
    /// reads are in flight while the workers compute.
    fn run<F, P>(&self, finish: F, prefetch: P, summary: &mut SearchSummary)
    where
        F: Fn(usize, &PossibleEntry, &mut Worker) -> Result<(), Halt> + Sync,
        P: Fn(usize, &PossibleEntry) + Sync,
    {
        let workers = rayon::current_num_threads().min(self.candidates.len());
        for (index, possible) in self.candidates.iter().enumerate().take(ENTRIES_AHEAD) {
            prefetch(index, possible);
        }
        let ends: Vec<WorkerEnd> = (0..workers)
            .into_par_iter()
            .map(|_| self.work_on(&finish, &prefetch, ENTRIES_AHEAD))
            .collect();
        for end in ends {
            summary.locate.add(&end.worker.locate);
            summary.align.add(&end.worker.align);
            if let (Completion::Complete, Err(error)) = (&summary.completion, end.result) {
                summary.completion = error.into();
            }
        }
        // A stop or a failure says more than a resource drop, so it takes precedence.
        if matches!(summary.completion, Completion::Complete) && summary.align.limited_entries > 0 {
            summary.completion = Completion::ResourceIncomplete;
        }
    }

    fn work_on<F, P>(&self, finish: &F, prefetch: &P, ahead: usize) -> WorkerEnd
    where
        F: Fn(usize, &PossibleEntry, &mut Worker) -> Result<(), Halt>,
        P: Fn(usize, &PossibleEntry),
    {
        let mut worker = Worker::default();
        while let Some((index, possible)) = self.take_indexed() {
            if self.halted.load(Ordering::Relaxed) {
                self.unfinished(possible.entry, SearchStage::Locate);
                continue;
            }
            if let Some(later) = self.candidates.get(index + ahead) {
                prefetch(index + ahead, later);
            }
            let result = match self.control.check() {
                Err(stopped) => Err(Halt {
                    stage: SearchStage::Locate,
                    error: SearchError::Stopped(stopped),
                }),
                Ok(()) => finish(index, possible, &mut worker),
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

    /// Reports every entry of the queue as unfinished in `stage`, when the search ends before
    /// any entry is taken.
    fn give_up(&self, stage: SearchStage) {
        while let Some(rest) = self.take() {
            self.unfinished(rest.entry, stage);
        }
    }

    fn take(&self) -> Option<&'a PossibleEntry> {
        self.take_indexed().map(|(_, possible)| possible)
    }

    /// The next entry and its place in the queue.
    fn take_indexed(&self) -> Option<(usize, &'a PossibleEntry)> {
        let index = self.next.fetch_add(1, Ordering::Relaxed);
        self.candidates.get(index).map(|possible| (index, possible))
    }

    fn unfinished(&self, entry: u32, stage: SearchStage) {
        // A dropped receiver means nobody listens any more.
        let _ = self.events.send(EntryEvent::Unfinished { entry, stage });
    }
}

/// Whether stage 3 gave up on an entry for a memory or work limit rather than for an error;
/// such an entry is dropped, the other entries go on, and the search ends `ResourceIncomplete`.
fn is_resource_limit(error: &AlignStageError) -> bool {
    matches!(
        error,
        AlignStageError::TooLarge { .. }
            | AlignStageError::Alignment(
                AlignmentError::MatrixTooLarge { .. } | AlignmentError::ChainingBudget(_)
            )
    )
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
