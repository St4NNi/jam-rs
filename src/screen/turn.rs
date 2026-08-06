//! Fragments of a batch, cut at the contig borders and turned around: the alignment of a batch
//! (the query) to an index entry becomes an alignment of the entry (now the query) to one
//! assembly contig.

use crate::align::{
    Alignment, AlignmentConfig, AlignmentError, EditOperation, EditRun, Interval, Strand,
};

/// One run cut out of an alignment: its operation and length, the query position of its start,
/// and the contig bases of the alignment before it.
#[derive(Clone, Copy)]
struct Cut {
    operation: EditOperation,
    length: u32,
    query: u64,
    contig_before: u64,
}

impl Cut {
    fn is_gap(&self) -> bool {
        matches!(
            self.operation,
            EditOperation::Insertion | EditOperation::Deletion
        )
    }
}

/// The score of runs: a gap costs its opening and every base.
fn score_runs(runs: &[EditRun], scores: &AlignmentConfig) -> Result<i32, AlignmentError> {
    let mut score = 0i64;
    for run in runs {
        let length = i64::from(run.length);
        score += match run.operation {
            EditOperation::Match => length * i64::from(scores.match_score),
            EditOperation::Substitution => length * i64::from(scores.mismatch_score),
            EditOperation::Insertion | EditOperation::Deletion => {
                i64::from(scores.gap_open_score) + length * i64::from(scores.gap_extend_score)
            }
        };
    }
    i32::try_from(score).map_err(|_| AlignmentError::LengthOverflow)
}

/// The runs of `alignment` that lie in the query positions `from..to`.
fn cut_runs(alignment: &Alignment, from: u64, to: u64) -> Vec<Cut> {
    let (low, high) = (from, to);
    let mut cuts = Vec::new();
    let (mut query, mut contig) = (alignment.query_interval.start, 0u64);
    for run in &alignment.runs {
        let length = u64::from(run.length);
        let (lo, hi) = (query.max(low), (query + length).min(high));
        match run.operation {
            EditOperation::Match | EditOperation::Substitution => {
                if lo < hi {
                    cuts.push(Cut {
                        operation: run.operation,
                        length: (hi - lo) as u32,
                        query: lo,
                        contig_before: contig + (lo - query),
                    });
                }
                query += length;
                contig += length;
            }
            // Bases only in the query.
            EditOperation::Insertion => {
                if lo < hi {
                    cuts.push(Cut {
                        operation: run.operation,
                        length: (hi - lo) as u32,
                        query: lo,
                        contig_before: contig,
                    });
                }
                query += length;
            }
            // Bases only in the contig: kept when they lie strictly inside the cut.
            EditOperation::Deletion => {
                if query > low && query < high {
                    cuts.push(Cut {
                        operation: run.operation,
                        length: run.length,
                        query,
                        contig_before: contig,
                    });
                }
                contig += length;
            }
        }
    }
    cuts
}

/// The part of `alignment` in the query positions `from..to`, trimmed so that it starts and ends
/// with a paired base, and scored again. `None` when no paired base is left.
pub(super) fn slice_alignment(
    alignment: &Alignment,
    from: u64,
    to: u64,
    scores: &AlignmentConfig,
) -> Result<Option<Alignment>, AlignmentError> {
    let cuts = cut_runs(alignment, from, to);
    let first = cuts.iter().position(|cut| !cut.is_gap());
    let last = cuts.iter().rposition(|cut| !cut.is_gap());
    let (Some(first), Some(last)) = (first, last) else {
        return Ok(None);
    };
    let cuts = &cuts[first..=last];
    let runs: Vec<EditRun> = cuts
        .iter()
        .map(|cut| EditRun {
            operation: cut.operation,
            length: cut.length,
        })
        .collect();
    let mut sliced = Alignment {
        score: score_runs(&runs, scores)?,
        strand: alignment.strand,
        query_interval: Interval::new(0, 0)?,
        contig_interval: Interval::new(0, 0)?,
        runs,
    };
    let counts = sliced.counts();
    let query_bases = counts.matches + counts.substitutions + counts.insertions;
    let contig_bases = counts.matches + counts.substitutions + counts.deletions;
    let (query_start, before) = (cuts[0].query, cuts[0].contig_before);
    sliced.query_interval = Interval::new(query_start, query_start + query_bases)?;
    // A forward alignment reads the contig from its start, a reverse one from its end.
    sliced.contig_interval = match alignment.strand {
        Strand::Forward => {
            let bottom = alignment.contig_interval.start + before;
            Interval::new(bottom, bottom + contig_bases)?
        }
        Strand::Reverse => {
            let top = alignment.contig_interval.end - before;
            Interval::new(top - contig_bases, top)?
        }
    };
    Ok(Some(sliced))
}

/// The same alignment with query and contig exchanged: bases only in the query become bases
/// only in the contig. A reverse alignment reads the old contig backwards, so its runs are
/// reversed to read the new query forwards.
pub(super) fn swap_roles(alignment: &Alignment) -> Alignment {
    let mut runs: Vec<EditRun> = alignment
        .runs
        .iter()
        .map(|run| EditRun {
            operation: match run.operation {
                EditOperation::Insertion => EditOperation::Deletion,
                EditOperation::Deletion => EditOperation::Insertion,
                other => other,
            },
            length: run.length,
        })
        .collect();
    if alignment.strand == Strand::Reverse {
        runs.reverse();
    }
    Alignment {
        score: alignment.score,
        strand: alignment.strand,
        query_interval: alignment.contig_interval,
        contig_interval: alignment.query_interval,
        runs,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{format_cigar, parse_cigar};
    use crate::seeds::complement;

    fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|&base| complement(base))
            .collect()
    }

    /// Checks that every paired base of `alignment` pairs equal bases (for `=`) or different ones
    /// (for `X`) between `query` and `contig`.
    fn check_bases(alignment: &Alignment, query: &[u8], contig: &[u8]) {
        alignment.validate().unwrap();
        let oriented = match alignment.strand {
            Strand::Forward => contig
                [alignment.contig_interval.start as usize..alignment.contig_interval.end as usize]
                .to_vec(),
            Strand::Reverse => reverse_complement(
                &contig[alignment.contig_interval.start as usize
                    ..alignment.contig_interval.end as usize],
            ),
        };
        let (mut q, mut c) = (alignment.query_interval.start as usize, 0usize);
        for run in &alignment.runs {
            for _ in 0..run.length {
                match run.operation {
                    EditOperation::Match => assert_eq!(query[q], oriented[c]),
                    EditOperation::Substitution => assert_ne!(query[q], oriented[c]),
                    _ => {}
                }
                match run.operation {
                    EditOperation::Match | EditOperation::Substitution => (q, c) = (q + 1, c + 1),
                    EditOperation::Insertion => q += 1,
                    EditOperation::Deletion => c += 1,
                }
            }
        }
    }

    /// A query and a contig with an alignment between them: 10 matching bases, a base only in the
    /// query, 10 more, a substitution, 10 more, two bases only in the contig and 10 more.
    fn example(strand: Strand) -> (Vec<u8>, Vec<u8>, Alignment) {
        let core = b"ACGTTGCAAGGCTAGCTAGGATCCATGCATGCAGGTACCAGTCAGT".to_vec();
        // The query holds an extra base after 10 bases, the contig a changed base after 20 and
        // two extra bases after 30.
        let mut query = b"TTT".to_vec();
        query.extend_from_slice(&core[..10]);
        query.push(b'G');
        query.extend_from_slice(&core[10..40]);
        let mut contig = b"CC".to_vec();
        contig.extend_from_slice(&core[..20]);
        contig.push(if core[20] == b'A' { b'C' } else { b'A' });
        contig.extend_from_slice(&core[21..30]);
        contig.extend_from_slice(b"AA");
        contig.extend_from_slice(&core[30..40]);
        let runs = parse_cigar("10=1I10=1X9=2D10=").unwrap();
        let mut alignment = Alignment {
            score: 0,
            strand,
            query_interval: Interval::new(3, 44).unwrap(),
            contig_interval: Interval::new(2, 44).unwrap(),
            runs,
        };
        if strand == Strand::Reverse {
            contig = reverse_complement(&contig);
            let length = contig.len() as u64;
            alignment.contig_interval = Interval::new(length - 44, length - 2).unwrap();
        }
        (query, contig, alignment)
    }

    #[test]
    fn slices_keep_a_valid_alignment_with_its_own_score() {
        let scores = AlignmentConfig::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            let (query, contig, alignment) = example(strand);
            check_bases(&alignment, &query, &contig);
            let whole = slice_alignment(&alignment, 0, 100, &scores)
                .unwrap()
                .unwrap();
            assert_eq!(whole.runs, alignment.runs);
            assert_eq!(whole.score, 2 * 39 - 3 - 6 - 7);
            for (from, to, cigar) in [
                (5, 30, "8=1I10=1X5="),
                (13, 14, "1I"),
                (14, 40, "10=1X9=2D6="),
                (33, 44, "1=2D10="),
                (34, 44, "10="),
            ] {
                let sliced = slice_alignment(&alignment, from, to, &scores).unwrap();
                if cigar == "1I" {
                    assert!(sliced.is_none());
                    continue;
                }
                let sliced = sliced.unwrap();
                assert_eq!(format_cigar(&sliced.runs), cigar, "{strand:?} {from}..{to}");
                check_bases(&sliced, &query, &contig);
            }
        }
    }

    #[test]
    fn turning_around_moves_the_gaps_to_the_other_side() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let (query, contig, alignment) = example(strand);
            let turned = swap_roles(&alignment);
            let (counts, back) = (alignment.counts(), turned.counts());
            assert_eq!(
                (back.insertions, back.deletions),
                (counts.deletions, counts.insertions)
            );
            // The old contig is the query now: on the reverse strand it is read as the reverse
            // complement of the old query.
            match strand {
                Strand::Forward => check_bases(&turned, &contig, &query),
                Strand::Reverse => {
                    assert_eq!(format_cigar(&turned.runs), "10=2I9=1X10=1D10=");
                    check_bases(&turned, &contig, &query);
                }
            }
            assert_eq!(swap_roles(&turned), alignment);
        }
    }
}
