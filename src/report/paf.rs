//! PAF lines: one per fragment of a confirmed entry, so standard tools read jam's alignments.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::decision::{DecisionRecord, FragmentRecord};
use crate::align::Strand;
use std::fmt::Write;

/// The PAF lines of a decision: first the primary fragments, then the alternatives. Mapping
/// quality is 255 (not available); `tp:A:P` marks a primary fragment and `tp:A:S` an alternative.
pub fn paf_lines(query: &str, decision: &DecisionRecord) -> Vec<String> {
    let primary = decision.fragments.iter().map(|fragment| (fragment, 'P'));
    let alternatives = decision.alternatives.iter().map(|fragment| (fragment, 'S'));
    primary
        .chain(alternatives)
        .map(|(fragment, kind)| paf_line(query, decision.query_length, fragment, kind))
        .collect()
}

fn paf_line(query: &str, query_length: u64, fragment: &FragmentRecord, kind: char) -> String {
    // The contig is the PAF target, and the CIGAR follows SAM with the contig as the reference.
    // A fragment over the end of a circular query ends past the query length.
    let query_start = fragment.query_segments.first().map_or(0, |span| span.start);
    let strand = match fragment.strand {
        Strand::Forward => '+',
        Strand::Reverse => '-',
    };
    let mut line = String::new();
    // Writing to a string cannot fail.
    let _ = write!(
        line,
        "{query}\t{query_length}\t{query_start}\t{}\t{strand}\t{}\t{}\t{}\t{}\t{}\t{}\t255",
        query_start + fragment.query_bases(),
        fragment.contig,
        fragment.contig_length,
        fragment.contig_interval.start,
        fragment.contig_interval.end,
        fragment.matches,
        fragment.columns(),
    );
    let edits = fragment.substitutions + fragment.insertions + fragment.deletions;
    let _ = write!(line, "\ttp:A:{kind}\tNM:i:{edits}\tcg:Z:{}", fragment.cigar);
    line
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::Interval;

    fn fragment(strand: Strand, segments: Vec<Interval>) -> FragmentRecord {
        FragmentRecord {
            contig: "plasmid_1".into(),
            contig_length: 5_000,
            query_segments: segments,
            contig_interval: Interval {
                start: 100,
                end: 109,
            },
            strand,
            score: 20,
            new_bases: 8,
            matches: 7,
            substitutions: 1,
            insertions: 0,
            deletions: 1,
            identity: 7.0 / 9.0,
            cigar: "3=1X3=1D1=".into(),
        }
    }

    fn decision(
        fragments: Vec<FragmentRecord>,
        alternatives: Vec<FragmentRecord>,
    ) -> DecisionRecord {
        DecisionRecord {
            query_length: 10,
            covered_bases: 8,
            query_coverage: 0.8,
            identity: 7.0 / 9.0,
            matches: 7,
            substitutions: 1,
            insertions: 0,
            deletions: 1,
            dropped_fragments: 0,
            fragments,
            alternatives,
        }
    }

    #[test]
    fn fragment_standard_columns() {
        let span = Interval { start: 2, end: 10 };
        let lines = paf_lines(
            "q1",
            &decision(vec![fragment(Strand::Reverse, vec![span])], vec![]),
        );
        assert_eq!(
            lines,
            [
                "q1\t10\t2\t10\t-\tplasmid_1\t5000\t100\t109\t7\t9\t255\ttp:A:P\tNM:i:2\tcg:Z:3=1X3=1D1="
            ]
        );
    }

    #[test]
    fn secondary_and_wrapped() {
        let wrapped = vec![
            Interval { start: 7, end: 10 },
            Interval { start: 0, end: 5 },
        ];
        let single = vec![Interval { start: 0, end: 8 }];
        let lines = paf_lines(
            "q1",
            &decision(
                vec![fragment(Strand::Forward, wrapped)],
                vec![fragment(Strand::Forward, single)],
            ),
        );
        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("q1\t10\t7\t15\t+\t"));
        assert!(lines[1].contains("\ttp:A:S\t"));
    }
}
