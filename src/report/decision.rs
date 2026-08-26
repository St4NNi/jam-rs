//! The record of a confirmed entry's decision: query coverage, identity, counts and the
//! fragments of its mosaic.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::ReportError;
use crate::align::{Alignment, EditRun, Interval, Strand, format_cigar};
use crate::decide::{Decision, Fragment};
use crate::index::Contig;
use serde::{Deserialize, Serialize};

/// One fragment of a mosaic, placed on the query and on a named contig.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct FragmentRecord {
    pub contig: String,
    pub contig_length: u64,
    /// The query bases the fragment spans; two segments when it runs over the end of a circular
    /// query.
    pub query_segments: Vec<Interval>,
    pub contig_interval: Interval,
    pub strand: Strand,
    pub score: i32,
    /// The query bases that only this fragment adds to the mosaic; 0 for an alternative.
    pub new_bases: u64,
    pub matches: u64,
    pub substitutions: u64,
    pub insertions: u64,
    pub deletions: u64,
    pub identity: f64,
    /// The CIGAR in SAM and PAF order: along the contig, from its lowest position on.
    pub cigar: String,
}

/// The decision about one entry: the totals of its mosaic and the fragments.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct DecisionRecord {
    pub query_length: u64,
    /// The query bases covered by the mosaic.
    pub covered_bases: u64,
    /// `covered_bases` divided by `query_length`.
    pub query_coverage: f64,
    /// Matches over alignment columns of the counted fragments together.
    pub identity: f64,
    pub matches: u64,
    pub substitutions: u64,
    pub insertions: u64,
    pub deletions: u64,
    /// Fragments left out because they fall below the fragment limits.
    pub dropped_fragments: u64,
    /// The primary fragments of the mosaic, best first.
    pub fragments: Vec<FragmentRecord>,
    /// Fragments that add no new query base, such as other copies of the same region.
    pub alternatives: Vec<FragmentRecord>,
}

impl FragmentRecord {
    /// The record of `fragment`; `contigs` is the contig table of its entry.
    pub fn new(
        fragment: &Fragment,
        new_bases: u64,
        contigs: &[Contig],
    ) -> Result<Self, ReportError> {
        let contig = contigs
            .get(fragment.contig_index as usize)
            .ok_or(ReportError::UnknownContig(fragment.contig_index))?;
        let alignment = &fragment.alignment;
        let counts = alignment.counts();
        Ok(Self {
            contig: contig.name.clone(),
            contig_length: contig.length,
            query_segments: fragment.query_segments.clone(),
            contig_interval: alignment.contig_interval,
            strand: alignment.strand,
            score: alignment.score,
            new_bases,
            matches: counts.matches,
            substitutions: counts.substitutions,
            insertions: counts.insertions,
            deletions: counts.deletions,
            identity: counts.identity(),
            cigar: contig_order_cigar(alignment),
        })
    }

    /// The alignment columns: matches, substitutions, insertions and deletions.
    pub fn columns(&self) -> u64 {
        self.matches + self.substitutions + self.insertions + self.deletions
    }

    /// The query bases the alignment spans, including the ones only in the query.
    pub fn query_bases(&self) -> u64 {
        self.matches + self.substitutions + self.insertions
    }
}

/// The CIGAR along the contig. An alignment on the reverse strand stores its runs in query order,
/// which runs backwards on the contig, so they are turned around.
fn contig_order_cigar(alignment: &Alignment) -> String {
    match alignment.strand {
        Strand::Forward => format_cigar(&alignment.runs),
        Strand::Reverse => {
            let runs: Vec<EditRun> = alignment.runs.iter().rev().copied().collect();
            format_cigar(&runs)
        }
    }
}

impl DecisionRecord {
    /// The record of `decision` for a query of `query_length` bases; `contigs` is the contig
    /// table of the entry.
    pub fn new(
        query_length: u64,
        decision: &Decision,
        contigs: &[Contig],
    ) -> Result<Self, ReportError> {
        let mosaic = &decision.mosaic;
        let counts = decision.counts;
        let fragments = mosaic
            .primary
            .iter()
            .map(|selected| FragmentRecord::new(&selected.fragment, selected.new_bases, contigs))
            .collect::<Result<_, _>>()?;
        let alternatives = mosaic
            .alternatives
            .iter()
            .map(|fragment| FragmentRecord::new(fragment, 0, contigs))
            .collect::<Result<_, _>>()?;
        Ok(Self {
            query_length,
            covered_bases: mosaic.covered_bases,
            query_coverage: mosaic.covered_bases as f64 / query_length.max(1) as f64,
            identity: decision.identity(),
            matches: counts.matches,
            substitutions: counts.substitutions,
            insertions: counts.insertions,
            deletions: counts.deletions,
            dropped_fragments: decision.dropped_fragments as u64,
            fragments,
            alternatives,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{Interval, parse_cigar};
    use crate::decide::DecisionRule;

    fn contigs() -> Vec<Contig> {
        vec![Contig {
            name: "plasmid_1".into(),
            length: 5_000,
        }]
    }

    fn fragment(strand: Strand, cigar: &str) -> Fragment {
        let runs = parse_cigar(cigar).unwrap();
        let query_end = 8;
        Fragment {
            contig_index: 0,
            query_segments: vec![Interval::new(0, query_end).unwrap()],
            alignment: Alignment {
                score: 20,
                strand,
                query_interval: Interval::new(0, query_end).unwrap(),
                contig_interval: Interval::new(100, 109).unwrap(),
                runs,
            },
        }
    }

    fn decision(fragment: Fragment) -> Decision {
        let rule = DecisionRule {
            min_fragment_columns: 1,
            ..DecisionRule::default()
        };
        rule.decide(8, vec![fragment], |_| 5_000).unwrap()
    }

    #[test]
    fn fragment_record_fields() {
        let forward = fragment(Strand::Forward, "3=1X3=1D1=");
        let record = FragmentRecord::new(&forward, 8, &contigs()).unwrap();
        assert_eq!(record.contig, "plasmid_1");
        assert_eq!(record.contig_length, 5_000);
        assert_eq!(record.cigar, "3=1X3=1D1=");
        assert_eq!((record.matches, record.columns()), (7, 9));
        let reverse = fragment(Strand::Reverse, "3=1X3=1D1=");
        let record = FragmentRecord::new(&reverse, 8, &contigs()).unwrap();
        assert_eq!(record.cigar, "1=1D3=1X3=");
        assert_eq!(record.strand, Strand::Reverse);
        assert!(FragmentRecord::new(&reverse, 0, &[]).is_err());
    }

    #[test]
    fn decision_round_trip() {
        let decision = decision(fragment(Strand::Forward, "3=1X3=1D1="));
        let record = DecisionRecord::new(8, &decision, &contigs()).unwrap();
        assert_eq!(record.fragments.len(), 1);
        assert_eq!(record.query_coverage, 1.0);
        let text = serde_json::to_string(&record).unwrap();
        let back: DecisionRecord = serde_json::from_str(&text).unwrap();
        assert_eq!(back, record);
    }
}
