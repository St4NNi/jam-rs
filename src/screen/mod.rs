//! The reverse direction: small sequences form the index and large assemblies stream past it once.
//!
//! The contigs of an assembly are joined into batches that are searched like queries against the
//! index of small sequences. A survey (stage 1) counts the seeds each contig shares with each
//! entry and keeps the entries that enough contigs share enough seeds with. The batches that
//! reached a kept entry are then located and aligned (stages 2 and 3). The fragments are cut at
//! the contig borders and turned around, so that they align an entry contig (the query of the
//! screen) to contigs of the assembly. The fragments of one entry contig over all contigs of the
//! assembly then get the decision, with the assembly contigs for the contig rule.

mod assembly;
mod batch;
mod records;
mod survey;

pub use assembly::Screener;
pub use records::AssemblySummary;

use crate::align::AlignmentError;
use crate::decide::{DecideError, DecisionRule};
use crate::index::IndexError;
use crate::search::align_stage::{AlignStageConfig, AlignStageError};
use crate::search::{SearchError, Stopped};
use std::path::PathBuf;
use thiserror::Error;

/// Why an assembly could not be screened.
#[derive(Debug, Error)]
pub enum ScreenError {
    /// An assembly cannot be read.
    #[error("bad assembly: {0}")]
    Input(String),
    /// A setting of the screen is out of range.
    #[error("invalid screen setting: {0}")]
    InvalidConfig(&'static str),
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Index(#[from] IndexError),
    #[error(transparent)]
    Search(#[from] SearchError),
    #[error(transparent)]
    AlignStage(#[from] AlignStageError),
    #[error(transparent)]
    Decide(#[from] DecideError),
    #[error(transparent)]
    Alignment(#[from] AlignmentError),
    #[error(transparent)]
    Stopped(#[from] Stopped),
}

/// The settings of a screen.
#[derive(Clone, Debug, PartialEq)]
pub struct ScreenConfig {
    /// The decision of an entry contig over the contigs of an assembly; its contig limits apply
    /// to the assembly contigs.
    pub rule: DecisionRule,
    /// Stage 3. Batches are linear queries.
    pub align: AlignStageConfig,
    /// Bases added on each side of a candidate bin to make an alignment window.
    pub flank: u64,
    /// Bases of contigs searched together, and the longest piece of a contig.
    pub batch_bases: usize,
    /// An assembly contig counts only for the entry contigs closest to it: those whose score
    /// (aligned share of the contig times identity) is within this margin of its best one.
    /// `None` counts it for every entry contig that it qualifies for.
    pub closest_margin: Option<f64>,
    /// A folder for the contigs of the reported entries (or, with `survey_only`, of the kept
    /// entries): one gzip FASTA for each assembly.
    pub contig_dir: Option<PathBuf>,
    /// Stop after the survey: no alignments and no reported entries, only the contig output.
    pub survey_only: bool,
}

impl ScreenConfig {
    /// The settings for `rule`. Chains with anchors too sparse for the identity asked for, with
    /// room for the noise of the estimate, are not aligned.
    pub fn new(rule: DecisionRule) -> Self {
        let align = AlignStageConfig {
            min_chain_identity: (rule.min_identity.max(rule.min_contig_identity) - 0.04).max(0.0),
            ..AlignStageConfig::default()
        };
        Self {
            rule,
            align,
            flank: 8_192,
            batch_bases: 500_000,
            closest_margin: None,
            contig_dir: None,
            survey_only: false,
        }
    }

    /// Checks that the settings make sense.
    pub fn validate(&self) -> Result<(), ScreenError> {
        self.rule.validate()?;
        self.align.validate()?;
        if self.flank == 0 {
            return Err(ScreenError::InvalidConfig("the flank must be at least 1"));
        }
        if !(1_000..u32::MAX as usize).contains(&self.batch_bases) {
            return Err(ScreenError::InvalidConfig(
                "batch bases must be from 1,000 to below 2^32",
            ));
        }
        if self
            .closest_margin
            .is_some_and(|margin| !(0.0..=1.0).contains(&margin))
        {
            return Err(ScreenError::InvalidConfig(
                "the closest margin must be from 0 to 1",
            ));
        }
        Ok(())
    }
}

impl Default for ScreenConfig {
    fn default() -> Self {
        Self::new(DecisionRule::default())
    }
}
