//! Tests of the decision and fragment records: names, counts and the contig order CIGAR.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::align::{Interval, parse_cigar};
use crate::decide::DecisionRule;

fn contigs() -> BTreeMap<u32, Contig> {
    let contig = Contig {
        name: "plasmid_1".into(),
        length: 5_000,
    };
    BTreeMap::from([(0, contig)])
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
        min_identity: 0.0,
        min_fragment_identity: 0.0,
        ..DecisionRule::default()
    };
    rule.decide(8, vec![fragment], &[], |_| 5_000).unwrap()
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
    assert!(FragmentRecord::new(&reverse, 0, &BTreeMap::new()).is_err());
}

#[test]
fn decision_round_trip() {
    let decision = decision(fragment(Strand::Forward, "3=1X3=1D1="));
    let record = DecisionRecord::new(&decision, &contigs()).unwrap();
    assert_eq!(record.fragments.len(), 1);
    assert_eq!(record.query_coverage, 1.0);
    let text = serde_json::to_string(&record).unwrap();
    let back: DecisionRecord = serde_json::from_str(&text).unwrap();
    assert_eq!(back, record);
}
