//! The screen of one assembly.

use super::batch::read_batches;
use super::records::AssemblySummary;
use super::{ScreenConfig, ScreenError};
use std::path::Path;

/// Screens assemblies against the entries of one shard. It can screen several assemblies at
/// once from different threads.
pub struct Screener {
    config: ScreenConfig,
}

impl Screener {
    /// A screener with `config`.
    pub fn new(config: ScreenConfig) -> Result<Self, ScreenError> {
        config.validate()?;
        Ok(Self { config })
    }

    /// Reads the assembly at `path` (plain or gzip FASTA) into batches. `name` names the assembly
    /// in the records.
    pub fn screen_assembly(&self, path: &Path, name: &str) -> Result<AssemblySummary, ScreenError> {
        let mut summary = AssemblySummary {
            kind: "assembly",
            assembly: name.to_string(),
            ..AssemblySummary::default()
        };
        let mut batches = 0;
        let contigs = read_batches(
            path,
            self.config.batch_bases,
            self.config.rule.min_fragment_columns as usize,
            &mut |_, _| {
                batches += 1;
                Ok(())
            },
            None,
        )?;
        let Some(contigs) = contigs else {
            return Ok(summary);
        };
        summary.contigs = contigs.list.len() as u64;
        summary.bases = contigs.bases;
        summary.skipped_contigs = contigs.skipped;
        summary.batches = batches;
        Ok(summary)
    }
}
