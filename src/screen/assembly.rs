//! The screen of one assembly: the survey of all its batches, the alignments of the batches that
//! reached a kept entry, the decisions, and the contig output.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::batch::{AssemblyContigs, Batch, read_batches};
use super::decision::decide_entries;
use super::fragments::{BatchAlignment, EntryContig, FragmentWork};
use super::records::{AssemblySummary, ScreenEntry};
use super::survey::{
    SEEDS_PART, SurveyLimits, Surveyed, expected_seeds, found_seeds, kept_seed_hits, survey,
};
use super::{ScreenConfig, ScreenError, write_contigs};
use crate::decide::Fragment;
use crate::index::Shard;
use crate::loss::MissModel;
use crate::search::QueryKmers;
use crate::search::SearchControl;
use crate::seeds::Seeder;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::Mutex;
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
    /// The shortest contig of each entry with more than one contig, once read.
    shortest: Mutex<HashMap<u32, u64>>,
}

/// The searched contigs of an assembly and, for each batch, the entries its pieces count for.
struct SurveyedAssembly {
    contigs: AssemblyContigs,
    surveyed: Vec<Vec<Surveyed>>,
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
        Ok(Self {
            shard,
            config,
            seeder,
            model,
            shortest: Mutex::new(HashMap::new()),
        })
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
        let Some(SurveyedAssembly { contigs, surveyed }) =
            self.survey_assembly(path, control, &mut summary)?
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
        let found = self.align_assembly(path, &surveyed, &kept, control, &mut summary)?;
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
        let mut surveyed = Vec::new();
        let contigs = read_batches(
            path,
            self.config.batch_bases,
            self.config.overlap(),
            rule.min_fragment_columns as usize,
            &mut |_, batch| {
                let found = found_seeds(seed_file, &self.seeder, batch, control)?;
                surveyed.push(survey(&self.seeder, batch, &found, limits, control)?);
                Ok(())
            },
            None,
        )?;
        let Some(contigs) = contigs else {
            return Ok(None);
        };
        summary.contigs = contigs.list.len() as u64;
        summary.bases = contigs.bases;
        summary.skipped_contigs = contigs.skipped;
        summary.batches = surveyed.len() as u64;
        Ok(Some(SurveyedAssembly { contigs, surveyed }))
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
            if seeds as f64 >= SEEDS_PART * expected_seeds(&self.seeder, required, identity) {
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
        read_batches(
            path,
            self.config.batch_bases,
            self.config.overlap(),
            self.config.rule.min_fragment_columns as usize,
            &mut |index, batch| {
                let Some(kept) = pieces.get(index).filter(|kept| !kept.is_empty()) else {
                    return Ok(());
                };
                self.align_batch(&aligner, batch, kept, &mut found, &mut work)
            },
            None,
        )?;
        summary.aligned_entries = work.aligned_entries;
        summary.alignments = work.align.alignments;
        summary.fragments = found.values().map(|list| list.len() as u64).sum();
        Ok(found)
    }

    /// Aligns the pieces of a batch that are kept for some entry, joined into a smaller batch, as
    /// preparing the query of stage 3 costs time in proportion to its length.
    fn align_batch(
        &self,
        aligner: &BatchAlignment<'_>,
        batch: &Batch,
        kept: &HashMap<u32, Vec<u32>>,
        found: &mut HashMap<EntryContig, Vec<Fragment>>,
        work: &mut FragmentWork,
    ) -> Result<(), ScreenError> {
        let mut pieces: Vec<u32> = kept.values().flatten().copied().collect();
        pieces.sort_unstable();
        pieces.dedup();
        let batch = batch.subset(&pieces);
        let kept: HashMap<u32, Vec<u32>> = kept
            .iter()
            .map(|(&entry, old)| {
                let new = old
                    .iter()
                    .filter_map(|piece| pieces.binary_search(piece).ok())
                    .map(|piece| piece as u32)
                    .collect();
                (entry, new)
            })
            .collect();
        let seed_file = self.shard.seed_file();
        let kmers = QueryKmers::new(&self.seeder, &batch.bases, false);
        let seeds = found_seeds(seed_file, &self.seeder, &batch, aligner.control)?;
        let seed_hits = kept_seed_hits(&kmers, &seeds, &kept);
        let (turned, part) = aligner.fragments(&batch, &kmers, seed_hits)?;
        work.add(&part);
        for turned in turned {
            found
                .entry((turned.entry, turned.entry_contig))
                .or_default()
                .push(turned.fragment);
        }
        Ok(())
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
