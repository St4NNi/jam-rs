//! Stages 2 and 3 of the screen for one batch: the kept entries are located and aligned with the
//! batch as the query.

use super::ScreenError;
use super::batch::Batch;
use crate::decide::{DecisionRule, Fragment};
use crate::index::Shard;
use crate::loss::{LossBudget, MissModel};
use crate::search::SearchControl;
use crate::search::align_stage::{
    AlignStageConfig, AlignStageWorkspace, AlignWork, PreparedQuery, entry_fragments,
};
use crate::search::bound::CoverageBound;
use crate::search::locate_stage::{
    LocateOutcome, LocateRequest, LocateWork, LocateWorkspace, locate_entry,
};
use crate::search::query::QueryKmers;
use crate::search::sketch_stage::PossibleEntry;
use rayon::prelude::*;

/// An entry contig, the query of the screen: the entry and the contig within it.
pub(super) type EntryContig = (u32, u32);

/// A fragment of an entry contig, found with a batch as the query.
#[derive(Clone, Debug)]
pub(super) struct ScreenFragment {
    pub entry: u32,
    /// The contig of the entry.
    pub entry_contig: u32,
    /// Its query is the batch.
    pub fragment: Fragment,
}

/// What the alignment of a batch needs besides the batch.
pub(super) struct BatchAlignment<'a> {
    pub shard: &'a Shard,
    pub model: &'a MissModel,
    pub rule: &'a DecisionRule,
    pub align: &'a AlignStageConfig,
    pub flank: u64,
    pub control: SearchControl<'a>,
}

/// The work of stages 2 and 3 over the batches of an assembly.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct FragmentWork {
    pub locate: LocateWork,
    pub align: AlignWork,
    /// Entries aligned, counted once per batch.
    pub aligned_entries: u64,
}

impl FragmentWork {
    pub fn add(&mut self, other: &Self) {
        self.locate.add(&other.locate);
        self.align.add(&other.align);
        self.aligned_entries += other.aligned_entries;
    }
}

impl BatchAlignment<'_> {
    /// Locates and aligns each entry with its seed hits in the batch.
    pub fn fragments(
        &self,
        batch: &Batch,
        kmers: &QueryKmers,
        seed_hits: Vec<(u32, Vec<u32>)>,
    ) -> Result<(Vec<ScreenFragment>, FragmentWork), ScreenError> {
        // The screen decides an entry by its own coverage, not by the coverage of the batch, so
        // the bound of stage 2 keeps every entry.
        let bound = CoverageBound::new(
            self.model,
            LossBudget::default(),
            self.rule.min_identity,
            kmers.length(),
            false,
            0,
        );
        let request = LocateRequest {
            query: kmers,
            bound: &bound,
            frequent_seeds: &[],
            flank: self.flank,
        };
        let query = PreparedQuery::new(&batch.bases, false, self.align)?;
        let results = seed_hits
            .into_par_iter()
            .map_init(
                || (LocateWorkspace::default(), AlignStageWorkspace::default()),
                |(locate, align), (entry, seed_hits)| {
                    let possible = PossibleEntry {
                        entry,
                        coverage_bound: 0,
                        containment_estimate: 0.0,
                        seed_hits,
                    };
                    self.entry(&request, &query, &possible, locate, align)
                },
            )
            .collect::<Result<Vec<_>, ScreenError>>()?;
        let mut work = FragmentWork::default();
        let mut turned = Vec::new();
        for (fragments, part) in results {
            turned.extend(fragments);
            work.add(&part);
        }
        Ok((turned, work))
    }

    fn entry(
        &self,
        request: &LocateRequest<'_>,
        query: &PreparedQuery,
        possible: &PossibleEntry,
        locate: &mut LocateWorkspace,
        align: &mut AlignStageWorkspace,
    ) -> Result<(Vec<ScreenFragment>, FragmentWork), ScreenError> {
        let mut work = FragmentWork::default();
        let outcome = locate_entry(
            self.shard,
            request,
            possible,
            locate,
            &mut work.locate,
            self.control,
        )?;
        let LocateOutcome::Located(located) = outcome else {
            return Ok((Vec::new(), work));
        };
        let packed = self.shard.packed(possible.entry)?;
        let (fragments, align_work) = entry_fragments(
            &packed,
            query,
            &located,
            self.align,
            self.rule,
            align,
            self.control,
        )?;
        work.align = align_work;
        work.aligned_entries = 1;
        let turned = fragments
            .into_iter()
            .map(|fragment| ScreenFragment {
                entry: possible.entry,
                entry_contig: fragment.contig_index,
                fragment,
            })
            .collect();
        Ok((turned, work))
    }
}
