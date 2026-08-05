//! Tests of the kernel: CIGAR, coordinates, endpoint completion, band geometry and scores.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::geometry::{WaveRange, wave_bounds, wave_range};
use super::*;
use crate::align::{format_cigar, parse_cigar};

fn config() -> AlignmentConfig {
    AlignmentConfig {
        band_width: 8,
        ..AlignmentConfig::default()
    }
}

fn input<'a>(query: &'a [u8], contig: &'a [u8], contig_offset: u64) -> AlignmentInput<'a> {
    AlignmentInput {
        query,
        query_offset: 0,
        contig,
        contig_offset,
    }
}

fn forward(
    workspace: &mut AlignmentWorkspace,
    query: &[u8],
    contig: &[u8],
    config: AlignmentConfig,
) -> Result<Alignment, AlignmentError> {
    workspace.align(input(query, contig, 0), Strand::Forward, config)
}

pub(super) fn test_dna(state: &mut u64, length: usize) -> Vec<u8> {
    (0..length)
        .map(|_| {
            *state ^= *state << 13;
            *state ^= *state >> 7;
            *state ^= *state << 17;
            b"ACGT"[(*state >> 33) as usize % 4]
        })
        .collect()
}

/// A copy of `source` with substitutions at `per_thousand`, and as many indels again when set.
pub(super) fn test_mutate(
    state: &mut u64,
    source: &[u8],
    per_thousand: u64,
    indels: bool,
) -> Vec<u8> {
    let mut next = |bound: u64| {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        (*state >> 20) % bound
    };
    let mut output = Vec::with_capacity(source.len() + 16);
    for &base in source {
        let roll = next(1000);
        if roll < per_thousand {
            output.push(b"ACGT"[next(4) as usize]);
        } else if indels && roll < per_thousand + per_thousand / 8 {
            // The base is left out.
        } else if indels && roll < per_thousand + per_thousand / 4 {
            output.push(base);
            output.push(b"ACGT"[next(4) as usize]);
        } else {
            output.push(base);
        }
    }
    output
}

#[test]
fn contig_base_deleted() {
    let mut workspace = AlignmentWorkspace::default();
    let alignment = forward(&mut workspace, b"ACGTACGT", b"ACGTGACGT", config()).unwrap();
    assert_eq!(format_cigar(&alignment.runs), "4=1D4=");
    assert_eq!(alignment.counts().deletions, 1);
    assert_eq!(alignment.identity(), 8.0 / 9.0);
    let alignment = forward(&mut workspace, b"ACGTGACGT", b"ACGTACGT", config()).unwrap();
    assert_eq!(format_cigar(&alignment.runs), "4=1I4=");
}

#[test]
fn reverse_interval_coordinates() {
    let mut workspace = AlignmentWorkspace::default();
    let alignment = workspace
        .align(
            input(b"AGGACTT", b"AAGTCCT", 100),
            Strand::Reverse,
            config(),
        )
        .unwrap();
    assert_eq!(alignment.contig_interval, Interval::new(100, 107).unwrap());
    assert_eq!(alignment.strand, Strand::Reverse);
}

#[test]
fn query_offset_shifts() {
    let mut workspace = AlignmentWorkspace::default();
    let shifted = AlignmentInput {
        query_offset: 50,
        ..input(b"TTACGTACGT", b"ACGTACGT", 7)
    };
    let alignment = workspace.align(shifted, Strand::Forward, config()).unwrap();
    assert_eq!(alignment.query_interval, Interval::new(52, 60).unwrap());
    assert_eq!(alignment.contig_interval, Interval::new(7, 15).unwrap());
}

fn core(strand: Strand, contig_interval: Interval) -> Alignment {
    Alignment {
        score: 8,
        strand,
        query_interval: Interval::new(2, 6).unwrap(),
        contig_interval,
        runs: parse_cigar("4=").unwrap(),
    }
}

#[test]
fn completes_both_endpoints() {
    let mut workspace = AlignmentWorkspace::default();
    let completed = workspace
        .complete_endpoints(
            core(Strand::Forward, Interval::new(102, 106).unwrap()),
            input(b"TTACGTAA", b"TTACGTAA", 100),
            2,
            config(),
        )
        .unwrap();
    assert_eq!(completed.query_interval, Interval::new(0, 8).unwrap());
    assert_eq!(completed.contig_interval, Interval::new(100, 108).unwrap());
    assert_eq!(format_cigar(&completed.runs), "8=");
    assert_eq!(completed.score, 16);
}

#[test]
fn reverse_completion_coordinates() {
    let mut workspace = AlignmentWorkspace::default();
    let completed = workspace
        .complete_endpoints(
            core(Strand::Reverse, Interval::new(102, 106).unwrap()),
            input(b"TTACGCAA", b"TTGCGTAA", 100),
            2,
            config(),
        )
        .unwrap();
    assert_eq!(completed.contig_interval, Interval::new(100, 108).unwrap());
    assert_eq!(format_cigar(&completed.runs), "8=");
}

#[test]
fn overhang_stays_free() {
    let mut workspace = AlignmentWorkspace::default();
    let core = Alignment {
        score: 2,
        strand: Strand::Forward,
        query_interval: Interval::new(0, 1).unwrap(),
        contig_interval: Interval::new(0, 1).unwrap(),
        runs: parse_cigar("1=").unwrap(),
    };
    let completed = workspace
        .complete_endpoints(core, input(b"AACGTACGT", b"AACGTACGTCCCC", 0), 64, config())
        .unwrap();
    assert_eq!(format_cigar(&completed.runs), "9=");
    assert_eq!(completed.contig_interval, Interval::new(0, 9).unwrap());
}

#[test]
fn wave_geometry_exhaustive() {
    let mut workspace = AlignmentWorkspace::default();
    for m in 0..=8 {
        for n in 0..=8 {
            for diagonal_offset in -12..=12 {
                for band_width in 0..=3 {
                    let config = AlignmentConfig {
                        diagonal_offset,
                        band_width,
                        ..AlignmentConfig::default()
                    };
                    let mut counts = vec![0usize; m + n + 1];
                    for i in 0..=m {
                        for j in 0..=n {
                            if (j as i64 - i as i64 - diagonal_offset).abs()
                                <= i64::from(band_width)
                            {
                                counts[i + j] += 1;
                            }
                        }
                    }
                    let expected = counts
                        .iter()
                        .position(|&count| count > 0)
                        .zip(counts.iter().rposition(|&count| count > 0));
                    assert_eq!(wave_bounds(m, n, config), expected, "{m} {n} {config:?}");
                    for (wave, &count) in counts.iter().enumerate() {
                        let width = wave_range(m, n, config, wave).map_or(0, WaveRange::width);
                        assert_eq!(width, count, "{m} {n} {wave} {config:?}");
                    }
                    for (first, last) in [(0, m + n), (0, (m + n) / 2), ((m + n) / 2, m + n)] {
                        let mut total = 0;
                        let mut offsets = Vec::new();
                        for &count in &counts[first..=last] {
                            offsets.push(total);
                            total += count;
                        }
                        offsets.push(total);
                        assert_eq!(
                            workspace.layout_chunk(m, n, config, first, last).unwrap(),
                            total
                        );
                        assert_eq!(workspace.wave_offsets, offsets);
                    }
                }
            }
        }
    }
    for diagonal_offset in [i64::MIN, i64::MAX] {
        let config = AlignmentConfig {
            diagonal_offset,
            band_width: u32::MAX,
            ..config()
        };
        assert_eq!(wave_bounds(8, 8, config), None);
    }
    let config = AlignmentConfig {
        diagonal_offset: 2,
        band_width: 0,
        ..config()
    };
    assert_eq!(wave_bounds(5, 9, config), Some((2, 12)));
    for wave in 2..=12 {
        assert_eq!(wave_range(5, 9, config, wave).is_some(), wave % 2 == 0);
    }
}

#[test]
fn skip_preserves_carries() {
    let range = Some(WaveRange { start: 2, end: 7 });
    let ranges = |carry: WaveCarry| {
        (
            carry.older.map(|range| (range.start, range.end)),
            carry.previous.map(|range| (range.start, range.end)),
        )
    };
    for older in [None, range] {
        for previous in [None, range] {
            for count in 0..=12 {
                let mut workspace = AlignmentWorkspace {
                    waves: [11, 22, 33].map(|value| vec![WaveChunk([value; 16])]),
                    ..AlignmentWorkspace::default()
                };
                let mut expected_waves = workspace.waves.clone();
                let mut carry = WaveCarry {
                    older,
                    previous,
                    base: 0,
                };
                let mut expected_carry = carry;
                for _ in 0..count {
                    expected_waves.rotate_left(1);
                    expected_carry.older = expected_carry.previous;
                    expected_carry.previous = None;
                }
                workspace.skip_empty_waves(count, &mut carry);
                assert_eq!(workspace.waves, expected_waves);
                assert_eq!(ranges(carry), ranges(expected_carry));
            }
        }
    }
}

/// The best local score by a plain full-matrix recurrence, without band.
fn reference_score(query: &[u8], contig: &[u8], config: AlignmentConfig) -> i32 {
    let (open, extend) = (config.gap_open(), config.gap_extend_score);
    let columns = contig.len() + 1;
    let mut previous = vec![[0i32; 3]; columns];
    let mut best = 0;
    for row in 1..=query.len() {
        let mut current = vec![[0i32; 3]; columns];
        for column in 1..columns {
            let diagonal = previous[column - 1].into_iter().max().unwrap_or(0);
            let same = query[row - 1].eq_ignore_ascii_case(&contig[column - 1]);
            let pair = if same {
                config.match_score
            } else {
                config.mismatch_score
            };
            let left = current[column - 1];
            let above = previous[column];
            current[column] = [
                (diagonal + pair).max(0),
                (left[0] + open)
                    .max(left[1] + extend)
                    .max(left[2] + open)
                    .max(0),
                (above[0] + open)
                    .max(above[1] + open)
                    .max(above[2] + extend)
                    .max(0),
            ];
            best = best.max(current[column].into_iter().max().unwrap_or(0));
        }
        previous = current;
    }
    best
}

#[test]
fn finds_reference_score() {
    let mut state = 0x2545_f491_4f6c_dd1du64;
    for round in 0..40 {
        let query = test_dna(&mut state, 20 + round * 3);
        let contig = test_mutate(&mut state, &query, 150, true);
        let config = AlignmentConfig {
            band_width: 400,
            ..AlignmentConfig::default()
        };
        let expected = reference_score(&query, &contig, config);
        for mut workspace in [AlignmentWorkspace::default(), AlignmentWorkspace::scalar()] {
            let found = forward(&mut workspace, &query, &contig, config);
            assert_eq!(found.map_or(0, |alignment| alignment.score), expected);
        }
    }
}
