//! The screen of one assembly: the survey of all its batches and the entries it keeps.

use super::batch::{AssemblyContigs, read_batches};
use super::records::AssemblySummary;
use super::survey::{SEEDS_PART, SurveyLimits, Surveyed, expected_seeds, found_seeds, survey};
use super::{ScreenConfig, ScreenError, write_contigs};
use crate::index::Shard;
use crate::search::SearchControl;
use crate::seeds::Seeder;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::Mutex;
use std::time::Instant;

/// Screens assemblies against the entries of one shard. It can screen several assemblies at
/// once from different threads.
pub struct Screener<'a> {
    shard: &'a Shard,
    config: ScreenConfig,
    seeder: Seeder,
    /// The shortest contig of each entry with more than one contig, once read.
    shortest: Mutex<HashMap<u32, u64>>,
}

/// The searched contigs of an assembly and, for each batch, the entries its pieces count for.
struct SurveyedAssembly {
    contigs: AssemblyContigs,
    surveyed: Vec<Vec<Surveyed>>,
}

impl<'a> Screener<'a> {
    /// A screener of `shard` with `config`.
    pub fn new(shard: &'a Shard, config: ScreenConfig) -> Result<Self, ScreenError> {
        config.validate()?;
        let scheme = shard.manifest().parameters.scheme();
        let seeder = Seeder::new(scheme).map_err(crate::index::IndexError::from)?;
        Ok(Self {
            shard,
            config,
            seeder,
            shortest: Mutex::new(HashMap::new()),
        })
    }

    /// Reads the assembly at `path` (plain or gzip FASTA) and surveys it. `name` names the
    /// assembly in the records.
    pub fn screen_assembly(
        &self,
        path: &Path,
        name: &str,
        control: SearchControl<'_>,
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
        let sketch = self.shard.sketch();
        let mut surveyed = Vec::new();
        let contigs = read_batches(
            path,
            self.config.batch_bases,
            rule.min_fragment_columns as usize,
            &mut |_, batch| {
                let found = found_seeds(sketch, &self.seeder, batch, control)?;
                surveyed.push(survey(
                    self.shard,
                    &self.seeder,
                    batch,
                    &found,
                    limits,
                    control,
                )?);
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
}
