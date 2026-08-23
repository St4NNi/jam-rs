use jam_rs::trace::SeedSensitivity;
use jam_rs::trace::model::BaseInterval;
use jam_rs::trace::seeds::{
    SeedError, expanded_seed_intervals, extract_seed_level, extract_seed_level_in_intervals,
    extract_seed_levels_in_intervals,
};

fn sequence(length: usize) -> Vec<u8> {
    b"ACGT".iter().copied().cycle().take(length).collect()
}

#[test]
fn rescue_intervals_keep_global_positions_and_only_hash_the_window() {
    let sequence = sequence(256);
    let config = SeedSensitivity {
        k: 31,
        scale: 1,
        max_occurrences: 32,
    };
    let interval = BaseInterval {
        start: 80,
        end: 128,
    };
    let flank = 7;
    let rescued = extract_seed_level_in_intervals(&sequence, config, &[interval], flank).unwrap();
    let expanded = expanded_seed_intervals(&[interval], sequence.len() as u64, flank).unwrap();
    assert_eq!(
        expanded,
        vec![BaseInterval {
            start: 73,
            end: 135
        }]
    );
    assert!(!rescued.seeds.is_empty());
    assert!(
        rescued
            .seeds
            .iter()
            .all(|seed| { seed.position >= 73 && seed.position + u64::from(config.k) <= 135 })
    );

    let whole = extract_seed_level(&sequence, config).unwrap();
    let expected = whole
        .seeds
        .into_iter()
        .filter(|seed| seed.position >= 73 && seed.position + u64::from(config.k) <= 135)
        .collect::<Vec<_>>();
    assert_eq!(rescued.seeds, expected);
}

#[test]
fn overlapping_gap_windows_are_merged_without_duplicate_seed_keys() {
    let sequence = sequence(256);
    let config = SeedSensitivity {
        k: 21,
        scale: 1,
        max_occurrences: 32,
    };
    let intervals = [
        BaseInterval {
            start: 60,
            end: 106,
        },
        BaseInterval {
            start: 96,
            end: 150,
        },
    ];
    let merged = expanded_seed_intervals(&intervals, sequence.len() as u64, 10).unwrap();
    assert_eq!(
        merged,
        vec![BaseInterval {
            start: 50,
            end: 160
        }]
    );
    let overlapping = extract_seed_level_in_intervals(&sequence, config, &intervals, 10).unwrap();
    let once = extract_seed_level_in_intervals(
        &sequence,
        config,
        &[BaseInterval {
            start: 50,
            end: 160,
        }],
        0,
    )
    .unwrap();
    assert_eq!(overlapping.seeds, once.seeds);
    assert!(
        overlapping
            .seeds
            .windows(2)
            .all(|seeds| seeds[0].position < seeds[1].position)
    );
}

#[test]
fn rescue_rounds_are_monotonic_and_density_levels_remain_nested() {
    let sequence = sequence(512);
    let primary = SeedSensitivity {
        k: 31,
        scale: 100,
        max_occurrences: 32,
    };
    let dense = SeedSensitivity {
        k: 31,
        scale: 1,
        max_occurrences: 32,
    };
    let rescue = SeedSensitivity {
        k: 21,
        scale: 100,
        max_occurrences: 32,
    };
    let first_round = extract_seed_level_in_intervals(
        &sequence,
        primary,
        &[BaseInterval {
            start: 32,
            end: 180,
        }],
        8,
    )
    .unwrap();
    let second_round = extract_seed_level_in_intervals(
        &sequence,
        dense,
        &[
            BaseInterval {
                start: 32,
                end: 180,
            },
            BaseInterval {
                start: 250,
                end: 360,
            },
        ],
        8,
    )
    .unwrap();
    assert!(first_round.seeds.iter().all(|seed| {
        second_round.seeds.iter().any(|candidate| {
            candidate.position == seed.position
                && candidate.hash == seed.hash
                && candidate.canonical_kmer == seed.canonical_kmer
        })
    }));

    let levels = extract_seed_levels_in_intervals(
        &sequence,
        &[dense, primary, rescue],
        &[BaseInterval {
            start: 250,
            end: 360,
        }],
        0,
    )
    .unwrap();
    let dense_level = &levels.levels[0];
    let primary_level = &levels.levels[1];
    assert!(
        primary_level
            .seeds
            .iter()
            .all(|seed| dense_level.seeds.contains(seed))
    );
    assert_eq!(levels.levels[2].k, 21);
}

#[test]
fn invalid_or_overflowing_rescue_intervals_fail_closed() {
    let sequence = sequence(64);
    let config = SeedSensitivity {
        k: 21,
        scale: 1,
        max_occurrences: 4,
    };
    assert_eq!(
        extract_seed_level_in_intervals(
            &sequence,
            config,
            &[BaseInterval { start: 40, end: 65 }],
            0,
        ),
        Err(SeedError::IntervalOutOfBounds {
            start: 40,
            end: 65,
            length: 64,
        })
    );
    assert_eq!(
        expanded_seed_intervals(&[BaseInterval { start: 1, end: 2 }], 64, u64::MAX,),
        Err(SeedError::CoordinateOverflow)
    );
    let empty = extract_seed_level_in_intervals(
        &sequence,
        config,
        &[BaseInterval { start: 12, end: 12 }],
        4,
    )
    .unwrap();
    assert!(empty.seeds.is_empty());
}
