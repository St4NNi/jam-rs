//! The driver of one search against one shard: it connects the stages, streams the entry events,
//! and sums the work.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::align_stage::{
    AlignStageConfig, AlignStageError, AlignStageWorkspace, AlignWork, Limited, PreparedQuery,
    align_packed_entry,
};
use super::control::{SearchControl, Stopped};
use super::error::SearchError;
use super::events::{Completion, DropReason, EntryEvent, SearchStage};
use super::locate_stage::{
    LocateOutcome, LocateRequest, LocateWork, LocateWorkspace, locate_entry, read_entry_locations,
    will_need_entry,
};
use super::located::LocatedEntry;
use super::query::{QueryKmers, lookup_bound};
use super::seed_mismatches::{LookupPlan, lookup_plan, search_bound, target_length};
use super::sketch_stage::{
    PossibleEntry, SketchLimits, SketchWork, held_seeds, keep_found_seeds, sketch_stage_known,
};
use super::times::{PartStart, PartTime, SearchTimes};
use crate::align::Interval;
use crate::decide::DecisionRule;
use crate::index::{SeedLocation, Shard, StoredSeed, seed_key};
use crate::loss::{MaxMiss, MissModel};
use crate::seeds::{Seeder, low_complexity_mask};
use crate::storage::lock;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::mpsc::SyncSender;
use std::sync::{Arc, Mutex};
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
    /// An entry that is left out after stage 1 and gets no events, such as the indexed record
    /// of the query itself.
    pub skip_entry: Option<u32>,
    /// The settings of stage 3.
    pub align: AlignStageConfig,
}

/// The sketch bytes per lookup above which a shared lookup pass pays. A query's own pass reads up
/// to its whole sketch part, at about 0.3 ns per byte on CephFS, and the shared pass makes its
/// lookups a second time, at about 70 ns each (measured on worm): 230 bytes, rounded up.
const SKETCH_BYTES_PER_LOOKUP: u64 = 256;
/// Candidate entries ahead of the workers whose reads are asked for, so on S3 as many requests
/// are in flight as the store allows (512); their blocks take a share of the block cache.
pub(super) const ENTRIES_AHEAD: usize = 512;

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
            skip_entry: None,
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
    /// One line for each entry that passed a limit of stage 3 and was decided once the limit
    /// was raised.
    pub raised: Vec<String>,
    /// One line for each entry that passed a limit of stage 3 that the memory budget could not
    /// raise, so the entry was decided within its limits (or, should even that fail, not
    /// decided).
    pub limited: Vec<String>,
}

/// The seed keys held by the sketch, found by `Searcher::prepare_kmers` for a batch of queries,
/// and one query's share of the time of the pass.
#[derive(Debug)]
pub struct PreparedKmers {
    held: Arc<Vec<(u64, StoredSeed)>>,
    time: PartTime,
}

/// The work of stage 3, summed over entries.
#[derive(Clone, Copy, Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct AlignStageWork {
    /// Located entries that stage 3 aligned.
    pub entries: u64,
    /// Entries that qualify.
    pub confirmed_entries: u64,
    /// Entries that passed a limit that the memory budget could not raise: decided within their
    /// limits, with part of their fragments or contigs (or, should even that fail, not decided).
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
            let held = Arc::new(held_seeds(self.shard.seed_file(), keys)?);
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
                    // The stored seeds of the held keys, so stage 1 reads no bucket again.
                    let mut known = Vec::new();
                    let kmers =
                        QueryKmers::filtered(&self.seeder, query, circular, &plan, |lookups| {
                            control.check()?;
                            match &held {
                                Some(held) => {
                                    let held = |hash| {
                                        held.binary_search_by_key(&seed_key(hash), |&(k, _)| k)
                                    };
                                    lookups.retain(|&(hash, _)| held(hash).is_ok());
                                }
                                None => known.extend(keep_found_seeds(seeds, lookups)?),
                            }
                            Ok::<(), SearchError>(())
                        });
                    known.sort_unstable_by_key(|&(key, _)| key);
                    known.dedup_by_key(|&mut (key, _)| key);
                    kmers.map(|kmers| (kmers, known))
                },
            );
            let prepared = prepared?;
            let (query_kmers, own_known) = match query_kmers {
                Ok(found) => found,
                Err(error) => {
                    let mut summary = self.empty_summary(query.len() as u64, plan.level());
                    summary.completion = error.into();
                    return Ok(summary);
                }
            };
            let workspace = AlignStageWorkspace::default();
            let passed = Mutex::new(Vec::new());
            let request = Request {
                config,
                prepared: prepared.as_ref(),
                workspace: &workspace,
                control,
                events,
                low_complexity,
                known: held.as_deref().map_or(&own_known, |held| held),
                passed: &passed,
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
        let bound = search_bound(&self.models, rule, budget, length, circular, plan);
        let mut summary = self.empty_summary(length, plan.full);
        summary.times.prepare_query = start.finish();
        let (miss_upper_bound, total_miss) = (bound.seed_miss_bound(), bound.total_miss());
        let miss_over_budget = !bound.prunes();
        summary.sketch = SketchWork {
            low_complexity_bases: low_complexity,
            sampled_kmers: plan.sampled,
            miss_upper_bound,
            total_miss,
            miss_over_budget,
            min_hit_starts: bound.min_hits(),
            ..SketchWork::default()
        };
        // Half of the required coverage must be normal bases, so a query with fewer normal bases
        // has no answer, and nothing is looked up.
        if length - low_complexity < rule.min_query_coverage.required_normal_bases(length) {
            return summary;
        }
        let start = PartStart::now();
        let mut sketch = match sketch_stage_known(
            self.shard.seed_file(),
            query_kmers,
            request.known,
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
        if let Some(skipped) = config.skip_entry {
            sketch
                .candidates
                .retain(|candidate| candidate.entry != skipped);
        }
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
        let passed = std::mem::take(&mut *lock(request.passed));
        self.align_again(passed, request, &mut summary);
        summary.times.entries = start.finish();
        summary
    }

    /// Aligns the entries that passed a limit of stage 3 again, one at a time on all threads of
    /// the search, and reports their last events. Each time an entry passes a limit, that limit
    /// is raised within the memory budget. An entry whose limit cannot be raised is decided
    /// within its limits, as the last resort. After a stop or a failure, the entries left are
    /// unfinished.
    fn align_again(
        &self,
        passed: Vec<(LocatedEntry, AlignStageError)>,
        request: &Request<'_>,
        summary: &mut SearchSummary,
    ) {
        let Some(prepared) = request.prepared else {
            return;
        };
        let threads = rayon::current_num_threads();
        let mut worker = Worker::default();
        let mut undecided = false;
        for (located, first) in passed {
            let entry = located.entry;
            if !matches!(summary.completion, Completion::Complete) {
                let unfinished = EntryEvent::Unfinished {
                    entry,
                    stage: SearchStage::Align,
                };
                // A dropped receiver means nobody listens any more.
                let _ = request.events.send(unfinished);
                continue;
            }
            let started = Instant::now();
            let (mut config, mut error) = (*request.config, first);
            let first = error.to_string();
            let last = loop {
                let passed = error.passed_limit();
                let raised =
                    passed.and_then(|(limit, needed)| config.align.raised(limit, needed, threads));
                let Some(raised) = raised else {
                    if config.align.keep_within_limits {
                        break Ok(Err(error));
                    }
                    // The last resort: the entry is decided within the limits it has.
                    config.align.keep_within_limits = true;
                    let within = Request {
                        config: &config,
                        ..*request
                    };
                    break self.align(prepared, &located, &within, &mut worker);
                };
                config.align = raised;
                let again = Request {
                    config: &config,
                    ..*request
                };
                match self.align(prepared, &located, &again, &mut worker) {
                    Ok(Err(passed)) => error = passed,
                    last => break last,
                }
            };
            worker.align.seconds += started.elapsed().as_secs_f64();
            let name = self.entry_name(entry);
            let budget = config.align.memory_budget;
            let last = match last {
                Ok(Ok((last, None))) => {
                    summary.raised.push(format!(
                        "entry {name} passed a limit ({first}) and was aligned again alone \
                         with it raised"
                    ));
                    last
                }
                Ok(Ok((last, Some(limited)))) => {
                    worker.align.limited_entries += 1;
                    let left_out = match limited.left_out_contigs {
                        0 => String::new(),
                        contigs => format!(", {contigs} contigs left out"),
                    };
                    summary.limited.push(format!(
                        "entry {name} was decided within its limits, as the memory budget of \
                         {budget} bytes cannot raise them: {}; {:.1}% of its fragment bytes \
                         kept{left_out}",
                        limited.limit,
                        100.0 * limited.kept
                    ));
                    last
                }
                Ok(Err(error)) => {
                    worker.align.limited_entries += 1;
                    undecided = true;
                    summary.limited.push(format!(
                        "entry {name} was not decided: {error}; it cannot be raised further \
                         within the memory budget of {budget} bytes"
                    ));
                    EntryEvent::Dropped {
                        entry,
                        stage: SearchStage::Align,
                        reason: DropReason::ResourceLimit,
                    }
                }
                Err(error) => {
                    summary.completion = error.into();
                    EntryEvent::Unfinished {
                        entry,
                        stage: SearchStage::Align,
                    }
                }
            };
            if request.events.send(last).is_err()
                && matches!(summary.completion, Completion::Complete)
            {
                summary.completion = Completion::Stopped(Stopped::Asked);
            }
        }
        summary.align.add(&worker.align);
        // A stop or a failure says more than an entry left undecided, so it takes precedence.
        if undecided && matches!(summary.completion, Completion::Complete) {
            summary.completion = Completion::ResourceIncomplete;
        }
    }

    /// The name of an entry of the shard, or its index when the manifest does not hold it.
    fn entry_name(&self, entry: u32) -> String {
        let entries = &self.shard.manifest().entries;
        entries
            .get(entry as usize)
            .map_or_else(|| entry.to_string(), |held| held.name.clone())
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
            raised: Vec::new(),
            limited: Vec::new(),
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
        match last {
            Ok((last, _)) => emit(events, last, SearchStage::Align),
            // Aligned again alone, once the other entries are done.
            Err(passed) => {
                lock(request.passed).push((located, passed));
                Ok(())
            }
        }
    }

    /// Stage 3 for one located entry: opens its packed sequence once, aligns and decides, and
    /// returns the last event of the entry, confirmed (with the contig table) or dropped, and how
    /// it was decided within its limits, or the error of a limit of stage 3 that it passed.
    fn align(
        &self,
        prepared: &PreparedQuery,
        located: &LocatedEntry,
        request: &Request<'_>,
        worker: &mut Worker,
    ) -> Result<Result<(EntryEvent, Option<Limited>), AlignStageError>, SearchError> {
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
            Err(error) if error.passed_limit().is_some() => return Ok(Err(error)),
            Err(error) => return Err(error.into()),
        };
        worker.align.counts.add(&aligned.work);
        let limited = aligned.limited.clone();
        if !aligned.decision.qualifies {
            return Ok(Ok((dropped(DropReason::Decision), limited)));
        }
        worker.align.confirmed_entries += 1;
        let mosaic = &aligned.decision.mosaic;
        let primary = mosaic.primary.iter().map(|selected| &selected.fragment);
        let fragments = primary.chain(&mosaic.alternatives);
        let contigs = packed.contigs(fragments.map(|fragment| fragment.contig_index))?;
        let confirmed = EntryEvent::Confirmed {
            aligned: Box::new(aligned),
            contigs,
        };
        Ok(Ok((confirmed, limited)))
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
    /// The stored seeds of seed keys that the lookups found already, ascending by key.
    known: &'a [(u64, StoredSeed)],
    /// The entries that passed a limit of stage 3, with the error of the limit.
    passed: &'a Mutex<Vec<(LocatedEntry, AlignStageError)>>,
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
