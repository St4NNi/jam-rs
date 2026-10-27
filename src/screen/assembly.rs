//! The screen of one assembly: the survey of all its batches, the alignments of the batches that
//! reached a kept entry, the decisions, and the contig output.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::batch::{AssemblyContigs, Batch, read_batches};
use super::decision::decide_entries;
use super::fragments::{BatchAlignment, EntryContig, FragmentWork, ScreenFragment};
use super::key_filter::KeyFilter;
use super::records::{AssemblySummary, ScreenEntry};
use super::survey::{
    PieceSeed, SEEDS_PART, SurveyLimits, Surveyed, expected_seeds, found_seeds, kept_seed_hits,
    query_seeds, survey,
};
use super::{ScreenConfig, ScreenError, write_contigs};
use crate::decide::Fragment;
use crate::index::{ColorSets, Shard};
use crate::loss::MissModel;
use crate::search::QueryKmers;
use crate::search::{SearchControl, Stopped};
use crate::seeds::Seeder;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::Mutex;
use std::sync::mpsc::sync_channel;
use std::time::Instant;

/// Screens assemblies against the entries of one shard, several at once from different threads.
/// The contigs of an assembly are joined into batches, surveyed (stage 1), and the batches that
/// reached a kept entry are located and aligned (stages 2 and 3), with the batch as the query.
pub struct Screener<'a> {
    shard: &'a Shard,
    config: ScreenConfig,
    seeder: Seeder,
    /// The miss model of stage 2, whose bound keeps every entry in the screen.
    model: MissModel,
    /// The stored seed keys of the shard, so that most k-mers skip the sketch buckets.
    filter: KeyFilter,
    /// The shortest contig of each entry with more than one contig, once read.
    shortest: Mutex<HashMap<u32, u64>>,
    /// The threads that survey and align the batches, when the caller is not one of them.
    pool: Option<&'a rayon::ThreadPool>,
}

/// The searched contigs of an assembly and, for each batch, the entries its pieces count for.
struct SurveyedAssembly {
    contigs: AssemblyContigs,
    surveyed: Vec<Vec<Surveyed>>,
    /// The found seeds of each batch in the pieces that count for some entry, and their color
    /// sets, so that stage 2 does not look them up again.
    found: Vec<BatchSeeds>,
}

/// Found seeds of the pieces of a batch, sorted by piece, and their color sets.
type BatchSeeds = (Vec<PieceSeed>, ColorSets);

/// Batches handled at once by the threads of one assembly, at most.
const MAX_GROUP: usize = 32;

/// The kept pieces of a batch joined into a smaller batch, as preparing the query of stage 3
/// costs time in proportion to its length, with the kept pieces of each entry and the found seeds
/// renumbered to it.
fn kept_subset(
    batch: &Batch,
    kept: &HashMap<u32, Vec<u32>>,
    (seeds, color_sets): BatchSeeds,
) -> (Batch, HashMap<u32, Vec<u32>>, BatchSeeds) {
    let mut pieces: Vec<u32> = kept.values().flatten().copied().collect();
    pieces.sort_unstable();
    pieces.dedup();
    let renumber = |piece: &u32| pieces.binary_search(piece).ok().map(|piece| piece as u32);
    let kept = (kept.iter())
        .map(|(&entry, old)| (entry, old.iter().filter_map(renumber).collect()))
        .collect();
    let seeds = (seeds.into_iter())
        .filter_map(|seed| {
            let piece = renumber(&seed.piece)?;
            Some(PieceSeed { piece, ..seed })
        })
        .collect();
    (batch.subset(&pieces), kept, (seeds, color_sets))
}

/// Receives every entry contig found in an assembly.
pub type EntrySink<'a> = dyn FnMut(ScreenEntry) -> Result<(), ScreenError> + 'a;

impl<'a> Screener<'a> {
    /// A screener of `shard` with `config`.
    pub fn new(shard: &'a Shard, config: ScreenConfig) -> Result<Self, ScreenError> {
        config.validate()?;
        if !shard.seed_file().has_locators() {
            let message = "the shard was built with the sketch part only and cannot be screened";
            return Err(crate::index::IndexError::Input(message.into()).into());
        }
        let scheme = shard.manifest().parameters.scheme;
        let seeder = Seeder::new(scheme).map_err(crate::index::IndexError::from)?;
        let model = MissModel::new(scheme, 0).map_err(crate::search::SearchError::from)?;
        let filter = KeyFilter::new(shard.seed_file())?;
        Ok(Self {
            shard,
            config,
            seeder,
            model,
            filter,
            shortest: Mutex::new(HashMap::new()),
            pool: None,
        })
    }

    /// Hands the batches to the threads of `pool`. A caller outside the pool then waits for its
    /// batches instead of taking up the work of other assemblies while it waits.
    pub fn with_pool(mut self, pool: &'a rayon::ThreadPool) -> Self {
        self.pool = Some(pool);
        self
    }

    /// Reads the assembly at `path` (plain or gzip FASTA), screens it and passes every entry
    /// contig found to `emit`. `name` names the assembly in the records and the contig output.
    pub fn screen_assembly(
        &self,
        path: &Path,
        name: &str,
        control: SearchControl<'_>,
        emit: &mut EntrySink<'_>,
    ) -> Result<AssemblySummary, ScreenError> {
        let mut summary = AssemblySummary {
            kind: "assembly",
            assembly: name.to_string(),
            ..AssemblySummary::default()
        };
        let started = Instant::now();
        let Some(SurveyedAssembly {
            contigs,
            surveyed,
            found,
        }) = self.survey_assembly(path, control, &mut summary)?
        else {
            return Ok(summary);
        };
        let kept = self.kept_entries(&surveyed)?;
        summary.surveyed_entries = kept.len() as u64;
        let matched: HashSet<u32> = surveyed
            .iter()
            .flatten()
            .filter(|found| kept.contains(&found.entry))
            .map(|found| found.contig)
            .collect();
        summary.matched_contigs = matched.len() as u64;
        summary.matched_bases = matched
            .iter()
            .map(|&contig| contigs.list[contig as usize].length)
            .sum();
        summary.survey_seconds = started.elapsed().as_secs_f64();
        if self.config.survey_only {
            if let Some(directory) = &self.config.contig_dir {
                write_contigs(path, name, &self.config, directory, &matched)?;
            }
            return Ok(summary);
        }
        let started = Instant::now();
        let found = self.align_assembly(path, &surveyed, found, &kept, control, &mut summary)?;
        let entry_contig = |key: EntryContig| self.entry_contig(key);
        let reported = decide_entries(
            found,
            &contigs.list,
            &self.config.rule,
            self.config.closest_margin,
            name,
            &entry_contig,
        )?;
        let mut reported_contigs = HashSet::new();
        for found in reported {
            reported_contigs.extend(found.contigs.iter().map(|contig| contig.id));
            summary.reported_entries += 1;
            emit(found)?;
        }
        summary.reported_contigs = reported_contigs.len() as u64;
        summary.reported_bases = reported_contigs
            .iter()
            .map(|&contig| contigs.list[contig as usize].length)
            .sum();
        if let Some(directory) = &self.config.contig_dir {
            write_contigs(path, name, &self.config, directory, &reported_contigs)?;
        }
        summary.align_seconds = started.elapsed().as_secs_f64();
        Ok(summary)
    }

    /// Reads the batches of the assembly and surveys each. `None` for an empty file.
    fn survey_assembly(
        &self,
        path: &Path,
        control: SearchControl<'_>,
        summary: &mut AssemblySummary,
    ) -> Result<Option<SurveyedAssembly>, ScreenError> {
        let rule = &self.config.rule;
        let limits = SurveyLimits {
            identity: rule.min_identity.max(rule.min_contig_identity),
            min_contig_coverage: rule.min_contig_coverage,
            closest_margin: self.config.closest_margin,
        };
        let seed_file = self.shard.seed_file();
        let survey_batch = |batch: &Batch| {
            let found = found_seeds(seed_file, &self.filter, &self.seeder, batch, control)?;
            let counted = survey(&self.seeder, batch, &found, limits, control)?;
            // Only the pieces that count for some entry can be aligned.
            let pieces: Vec<u32> = counted.iter().map(|counted| counted.piece).collect();
            let (mut seeds, color_sets) = found;
            seeds.retain(|seed| pieces.binary_search(&seed.piece).is_ok());
            let seeds = match seeds.is_empty() {
                true => BatchSeeds::default(),
                false => (seeds, color_sets),
            };
            Ok::<_, ScreenError>((counted, seeds))
        };
        let (mut surveyed, mut kept_seeds) = (Vec::new(), Vec::new());
        let contigs = self.grouped_batches(
            path,
            |_, batch| Some(batch.clone()),
            |group| {
                let done = (group.par_iter())
                    .map(survey_batch)
                    .collect::<Result<Vec<_>, _>>()?;
                for (counted, seeds) in done {
                    surveyed.push(counted);
                    kept_seeds.push(seeds);
                }
                Ok(())
            },
        )?;
        let Some(contigs) = contigs else {
            return Ok(None);
        };
        summary.contigs = contigs.list.len() as u64;
        summary.bases = contigs.bases;
        summary.skipped_contigs = contigs.skipped;
        summary.batches = surveyed.len() as u64;
        Ok(Some(SurveyedAssembly {
            contigs,
            surveyed,
            found: kept_seeds,
        }))
    }

    /// The entries that the surveyed contigs of the assembly share enough seeds with to cover the
    /// query coverage that is asked for.
    fn kept_entries(&self, surveyed: &[Vec<Surveyed>]) -> Result<HashSet<u32>, ScreenError> {
        let mut shared: HashMap<u32, u64> = HashMap::new();
        for found in surveyed.iter().flatten() {
            *shared.entry(found.entry).or_default() += found.seeds;
        }
        let rule = &self.config.rule;
        let identity = rule.min_identity.max(rule.min_contig_identity);
        let mut kept = HashSet::new();
        for (entry, seeds) in shared {
            let length = self.query_length(entry)?;
            let required = rule.min_query_coverage.required_bases(length) as f64;
            // Rounded down, so that a short entry is kept with one seed fewer than its share.
            let needed = (SEEDS_PART * expected_seeds(&self.seeder, required, identity)).floor();
            if seeds as f64 >= needed.max(1.0) {
                kept.insert(entry);
            }
        }
        Ok(kept)
    }

    /// Reads the batches again and aligns those that reached a kept entry. Returns the fragments
    /// of each entry contig, cut at the contig borders and turned around, for the decision.
    fn align_assembly(
        &self,
        path: &Path,
        surveyed: &[Vec<Surveyed>],
        mut seeds: Vec<BatchSeeds>,
        kept: &HashSet<u32>,
        control: SearchControl<'_>,
        summary: &mut AssemblySummary,
    ) -> Result<HashMap<EntryContig, Vec<Fragment>>, ScreenError> {
        let pieces: Vec<HashMap<u32, Vec<u32>>> = surveyed
            .iter()
            .map(|batch| {
                let mut pieces: HashMap<u32, Vec<u32>> = HashMap::new();
                for found in batch.iter().filter(|found| kept.contains(&found.entry)) {
                    pieces.entry(found.entry).or_default().push(found.piece);
                }
                pieces
            })
            .collect();
        let mut found: HashMap<EntryContig, Vec<Fragment>> = HashMap::new();
        if pieces.iter().all(HashMap::is_empty) {
            return Ok(found);
        }
        let aligner = BatchAlignment {
            shard: self.shard,
            model: &self.model,
            rule: &self.config.rule,
            align: &self.config.align,
            flank: self.config.flank,
            control,
        };
        let mut work = FragmentWork::default();
        self.grouped_batches(
            path,
            |index, batch| {
                let kept = pieces.get(index).filter(|kept| !kept.is_empty())?;
                let seeds = std::mem::take(&mut seeds[index]);
                Some(kept_subset(batch, kept, seeds))
            },
            |group| {
                let done = (group.into_par_iter())
                    .map(|(batch, kept, seeds)| self.align_batch(&aligner, &batch, &kept, seeds))
                    .collect::<Result<Vec<_>, _>>()?;
                // In batch order, so that the fragments of an entry contig keep their order.
                for (turned, part) in done {
                    work.add(&part);
                    for turned in turned {
                        found
                            .entry((turned.entry, turned.entry_contig))
                            .or_default()
                            .push(turned.fragment);
                    }
                }
                Ok(())
            },
        )?;
        summary.aligned_entries = work.aligned_entries;
        summary.alignments = work.align.alignments;
        summary.fragments = found.values().map(|list| list.len() as u64).sum();
        Ok(found)
    }

    /// Aligns the pieces of a batch that are kept for some entry, joined into a smaller batch by
    /// `kept_subset`, with their found seeds. Returns the turned fragments and the work.
    fn align_batch(
        &self,
        aligner: &BatchAlignment<'_>,
        batch: &Batch,
        kept: &HashMap<u32, Vec<u32>>,
        (seeds, color_sets): BatchSeeds,
    ) -> Result<(Vec<ScreenFragment>, FragmentWork), ScreenError> {
        let kmers = QueryKmers::new(&self.seeder, &batch.bases, false);
        let seed_hits = kept_seed_hits(&kmers, &seeds, &color_sets, kept);
        let query_seeds = query_seeds(&kmers, &seeds, color_sets);
        aligner.fragments(batch, &kmers, seed_hits, &query_seeds)
    }

    /// Reads the batches of the assembly on another thread, which turns each into an item with
    /// `item`, and passes the items in order to `group`, a few at a time, so that the threads
    /// share the work of several batches. `None` for an empty file.
    fn grouped_batches<T: Send>(
        &self,
        path: &Path,
        mut item: impl FnMut(usize, &Batch) -> Option<T> + Send,
        mut group: impl FnMut(Vec<T>) -> Result<(), ScreenError> + Send,
    ) -> Result<Option<AssemblyContigs>, ScreenError> {
        let threads = self.pool.map_or_else(rayon::current_num_threads, |pool| {
            pool.current_num_threads()
        });
        let size = threads.clamp(2, MAX_GROUP);
        let mut group = |items: Vec<T>| match self.pool {
            Some(pool) => pool.install(|| group(items)),
            None => group(items),
        };
        let (sender, receiver) = sync_channel::<T>(size);
        std::thread::scope(|scope| {
            let reader = scope.spawn(move || {
                read_batches(
                    path,
                    self.config.batch_bases,
                    self.config.overlap(),
                    self.config.rule.min_fragment_columns as usize,
                    &mut |index, batch| match item(index, batch) {
                        Some(item) => sender.send(item).map_err(|_| Stopped::Asked.into()),
                        None => Ok(()),
                    },
                    None,
                )
            });
            let mut result = Ok(());
            let mut pending = Vec::with_capacity(size);
            for item in receiver.iter() {
                pending.push(item);
                if pending.len() == size {
                    result = group(std::mem::take(&mut pending));
                    if result.is_err() {
                        break;
                    }
                }
            }
            if result.is_ok() && !pending.is_empty() {
                result = group(pending);
            }
            // A failed group drops the receiver, which stops the reader.
            drop(receiver);
            let read = reader
                .join()
                .unwrap_or_else(|panic| std::panic::resume_unwind(panic));
            result.and(read)
        })
    }

    /// The length the survey compares an entry with: its bases, or for an entry with more than
    /// one contig its shortest contig, as each contig is a query of its own.
    fn query_length(&self, entry: u32) -> Result<u64, ScreenError> {
        let record = &self.shard.manifest().entries[entry as usize];
        if record.contigs <= 1 {
            return Ok(record.bases);
        }
        let mut shortest = self
            .shortest
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(&length) = shortest.get(&entry) {
            return Ok(length);
        }
        let packed = self.shard.packed(entry)?;
        let mut length = u64::MAX;
        for contig in 0..packed.contig_count() {
            let range = packed.contig_range(contig)?;
            length = length.min(range.end - range.start);
        }
        shortest.insert(entry, length);
        Ok(length)
    }

    /// The name and length of an entry contig: the entry name, with the contig name when the
    /// entry has more than one contig.
    fn entry_contig(&self, (entry, contig): EntryContig) -> Result<(String, u64), ScreenError> {
        let record = &self.shard.manifest().entries[entry as usize];
        let packed = self.shard.packed(entry)?;
        let range = packed.contig_range(contig as usize)?;
        let length = range.end - range.start;
        if record.contigs <= 1 {
            return Ok((record.name.clone(), length));
        }
        let contig_name = packed.contig_name(contig as usize)?;
        let contig_name = String::from_utf8_lossy(&contig_name);
        Ok((format!("{}/{contig_name}", record.name), length))
    }
}
