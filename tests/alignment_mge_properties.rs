use jam_rs::trace::alignment::{
    AlignmentError, AlignmentOptions, AlignmentResult, align, align_circular, align_oriented,
};
use jam_rs::trace::config::AlignmentScoring;
use jam_rs::trace::model::{EditOperation, Strand};

fn scoring() -> AlignmentScoring {
    AlignmentScoring {
        match_score: 2,
        mismatch_score: -3,
        gap_open_score: -5,
        gap_extend_score: -1,
        band_width: 64,
    }
}

#[test]
fn generated_queries_match_full_affine_oracle_and_reconstruct() {
    let mut rng = Rng::new(0x5eed_3121_2026);
    for case in 0..32 {
        let query = generated_query(case, &mut rng);
        let target_core = transformed_query(&query, case, &mut rng);
        let target = with_flanks(&target_core, case, &mut rng);
        let result = align(&query, &target, scoring()).unwrap();

        assert_eq!(
            result.score,
            oracle_local_affine_score(&query, &target, scoring()),
            "oracle mismatch for generated case {case}"
        );
        assert_alignment_invariants(&query, &target, &result, scoring());
    }
}

#[test]
fn generated_reverse_targets_preserve_forward_coordinates() {
    let mut rng = Rng::new(0xface_cafe_3121);
    for case in 0..24 {
        let query = generated_query(case + 7, &mut rng);
        let oriented_core = transformed_query(&query, case + 7, &mut rng);
        let oriented_target = with_flanks(&oriented_core, case + 7, &mut rng);
        let target_forward = reverse_complement(&oriented_target);
        let target_start = 37 + u64::from(case as u32);

        let direct = align(&query, &oriented_target, scoring()).unwrap();
        let reverse = align_oriented(
            &query,
            &target_forward,
            target_start,
            Strand::Reverse,
            scoring(),
        )
        .unwrap();

        let mapped_start = target_forward.len() as u64 - direct.target_interval.end + target_start;
        let mapped_end = target_forward.len() as u64 - direct.target_interval.start + target_start;
        assert_eq!(
            reverse.score, direct.score,
            "score mismatch for reverse case {case}"
        );
        assert_eq!(
            reverse.cigar, direct.cigar,
            "CIGAR mismatch for reverse case {case}"
        );
        assert_eq!(
            reverse.target_interval.start, mapped_start,
            "reverse start mismatch for case {case}"
        );
        assert_eq!(
            reverse.target_interval.end, mapped_end,
            "reverse end mismatch for case {case}"
        );
        assert_alignment_invariants(&query, &oriented_target, &direct, scoring());
        assert!(reverse.query_interval.end <= query.len() as u64);
    }
}

#[test]
fn ambiguous_and_empty_inputs_are_explicit_and_bounded() {
    let query = b"acNRYswkmbdhvACGT";
    let target = b"TTacNRYswkmbdhvACGTAA";
    let result = align(query, target, scoring()).unwrap();
    assert_eq!(result.matches, query.len() as u64);
    assert_alignment_invariants(query, target, &result, scoring());

    assert!(matches!(
        align(b"", b"ACGT", scoring()),
        Err(AlignmentError::EmptyQuery)
    ));
    assert!(matches!(
        align(b"ACGT", b"", scoring()),
        Err(AlignmentError::EmptyTarget)
    ));

    let mut narrow_scoring = scoring();
    narrow_scoring.band_width = 0;
    let narrow = AlignmentOptions::new(narrow_scoring).with_diagonal_offset(100);
    assert!(matches!(
        align(b"ACGTACGT", b"ACGTACGT", narrow),
        Err(AlignmentError::BandExcludesInput)
    ));
}

#[test]
fn circular_full_span_is_ordered_and_does_not_exceed_query() {
    let query = b"ACGTACGTAC";
    let linearized = b"ACACGTACGT";
    let result = align_circular(
        query,
        8,
        query.len() as u64,
        linearized,
        0,
        Strand::Forward,
        scoring(),
    )
    .unwrap();

    assert!(result.origin_crossing);
    assert_eq!(result.query_segments.len(), 2);
    assert_eq!(
        result
            .query_segments
            .iter()
            .map(|segment| segment.len())
            .sum::<u64>(),
        query.len() as u64
    );
    assert_alignment_invariants(linearized, linearized, &result, scoring());
}

fn assert_alignment_invariants(
    query: &[u8],
    target: &[u8],
    result: &AlignmentResult,
    scoring: AlignmentScoring,
) {
    assert!(result.query_interval.start <= result.query_interval.end);
    assert!(result.target_interval.start <= result.target_interval.end);
    assert!(result.query_interval.end <= query.len() as u64);
    assert!(result.target_interval.end <= target.len() as u64);
    assert_eq!(result.query_length, result.query_interval.len());
    assert_eq!(result.target_length, result.target_interval.len());

    let (aligned_query, aligned_target) = result.reconstruct(query, target).unwrap();
    assert_eq!(aligned_query.len(), aligned_target.len());

    let mut query_consumed = 0u64;
    let mut target_consumed = 0u64;
    let mut matches = 0u64;
    let mut substitutions = 0u64;
    let mut insertions = 0u64;
    let mut deletions = 0u64;
    let mut score = 0i32;
    let mut cigar = String::new();
    let mut column = 0usize;

    for run in &result.edit_script {
        assert!(run.length > 0);
        let (code, consumes_query, consumes_target) = match run.operation {
            EditOperation::Equal => {
                matches += u64::from(run.length);
                score += scoring.match_score * run.length as i32;
                ('=', true, true)
            }
            EditOperation::Substitution => {
                substitutions += u64::from(run.length);
                score += scoring.mismatch_score * run.length as i32;
                ('X', true, true)
            }
            EditOperation::Insertion => {
                insertions += u64::from(run.length);
                score += scoring.gap_open_score + scoring.gap_extend_score * run.length as i32;
                ('I', false, true)
            }
            EditOperation::Deletion => {
                deletions += u64::from(run.length);
                score += scoring.gap_open_score + scoring.gap_extend_score * run.length as i32;
                ('D', true, false)
            }
            EditOperation::SoftClip => panic!("alignment kernel must not emit soft clips"),
        };

        if consumes_query {
            query_consumed += u64::from(run.length);
        }
        if consumes_target {
            target_consumed += u64::from(run.length);
        }
        let _ = std::fmt::Write::write_fmt(&mut cigar, format_args!("{}{}", run.length, code));

        for _ in 0..run.length {
            let query_base = aligned_query[column];
            let target_base = aligned_target[column];
            match run.operation {
                EditOperation::Equal => {
                    assert_ne!(query_base, b'-');
                    assert_ne!(target_base, b'-');
                    assert!(query_base.eq_ignore_ascii_case(&target_base));
                }
                EditOperation::Substitution => {
                    assert_ne!(query_base, b'-');
                    assert_ne!(target_base, b'-');
                    assert!(!query_base.eq_ignore_ascii_case(&target_base));
                }
                EditOperation::Insertion => {
                    assert_eq!(query_base, b'-');
                    assert_ne!(target_base, b'-');
                }
                EditOperation::Deletion => {
                    assert_ne!(query_base, b'-');
                    assert_eq!(target_base, b'-');
                }
                EditOperation::SoftClip => unreachable!(),
            }
            column += 1;
        }
    }

    assert_eq!(column, aligned_query.len());
    assert_eq!(query_consumed, result.query_interval.len());
    assert_eq!(target_consumed, result.target_interval.len());
    assert_eq!(matches, result.matches);
    assert_eq!(substitutions, result.substitutions);
    assert_eq!(insertions, result.insertions);
    assert_eq!(deletions, result.deletions);
    assert_eq!(cigar, result.cigar);
    assert_eq!(score, result.score);
}

fn generated_query(case: usize, rng: &mut Rng) -> Vec<u8> {
    let length = 18 + (case % 19);
    match case % 4 {
        0 => vec![b'A'; length],
        1 => (0..length).map(|index| b"ACGT"[index % 4]).collect(),
        2 => (0..length).map(|_| rng.base()).collect(),
        _ => (0..length)
            .map(|index| b"ACNRYGT"[(index + case) % 7])
            .collect(),
    }
}

fn transformed_query(query: &[u8], case: usize, rng: &mut Rng) -> Vec<u8> {
    let mut transformed = Vec::with_capacity(query.len() + 8);
    for (index, &base) in query.iter().enumerate() {
        if case % 5 == 2 && index % 11 == 4 {
            continue;
        }
        transformed.push(if case % 5 == 1 && index % 9 == 3 {
            mutate_base(base)
        } else {
            base
        });
        if case % 5 == 3 && index % 13 == 6 {
            transformed.push(rng.base());
        }
    }
    transformed
}

fn with_flanks(core: &[u8], case: usize, rng: &mut Rng) -> Vec<u8> {
    let prefix_len = case % 5;
    let suffix_len = (case / 3) % 5;
    let mut target = Vec::with_capacity(prefix_len + core.len() + suffix_len);
    for _ in 0..prefix_len {
        target.push(rng.base());
    }
    target.extend_from_slice(core);
    for _ in 0..suffix_len {
        target.push(rng.base());
    }
    target
}

fn mutate_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    }
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            _ => b'N',
        })
        .collect()
}

fn oracle_local_affine_score(query: &[u8], target: &[u8], scoring: AlignmentScoring) -> i32 {
    let mut matching = vec![vec![0i32; target.len() + 1]; query.len() + 1];
    let mut insertion = matching.clone();
    let mut deletion = matching.clone();
    let mut best = 0i32;

    for i in 1..=query.len() {
        for j in 1..=target.len() {
            let substitution = if query[i - 1].eq_ignore_ascii_case(&target[j - 1]) {
                scoring.match_score
            } else {
                scoring.mismatch_score
            };
            let previous = matching[i - 1][j - 1]
                .max(insertion[i - 1][j - 1])
                .max(deletion[i - 1][j - 1]);
            matching[i][j] = (previous + substitution).max(0);
            insertion[i][j] = (insertion[i][j - 1] + scoring.gap_extend_score)
                .max(matching[i][j - 1] + scoring.gap_open_score + scoring.gap_extend_score)
                .max(deletion[i][j - 1] + scoring.gap_open_score + scoring.gap_extend_score)
                .max(0);
            deletion[i][j] = (deletion[i - 1][j] + scoring.gap_extend_score)
                .max(matching[i - 1][j] + scoring.gap_open_score + scoring.gap_extend_score)
                .max(insertion[i - 1][j] + scoring.gap_open_score + scoring.gap_extend_score)
                .max(0);
            best = best
                .max(matching[i][j])
                .max(insertion[i][j])
                .max(deletion[i][j]);
        }
    }
    best
}

struct Rng(u64);

impl Rng {
    const fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    fn base(&mut self) -> u8 {
        b"ACGT"[(self.next() & 3) as usize]
    }
}
