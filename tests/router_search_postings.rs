use jam_rs::router::postings::{
    InMemoryPostingSource, PositionPosting, PostingCodec, PostingDecodeStats, PostingInput,
    PostingSource,
};
use jam_rs::router::search::{RouterSearchConfig, TieredQueryWitness, search};
use jam_rs::router::{QueryWitness, WitnessHandoffMode, WitnessKey};

fn key(value: u64) -> WitnessKey {
    WitnessKey::from_packed(value).expect("test packed kmer must be a valid witness")
}

fn query(key: WitnessKey, position: u64, windows: &[u32]) -> QueryWitness {
    QueryWitness {
        key,
        query_position: position,
        query_reverse: false,
        query_window_ids: windows.to_vec(),
    }
}

fn input(tier: u32, key: WitnessKey, documents: Vec<u32>) -> PostingInput {
    PostingInput::new(tier, key, documents)
}

fn default_config() -> RouterSearchConfig {
    RouterSearchConfig {
        top_k: 10,
        accumulator_capacity: 10,
        ..RouterSearchConfig::default()
    }
}

#[test]
fn inline_and_delta_postings_round_trip_sorted_ids() {
    let mut source = InMemoryPostingSource::new((0..300).map(|id| format!("m{id:03}")).collect());
    let inline_key = key(1);
    let block_key = key(2);
    source.insert(input(20, inline_key, vec![7, 1, 7])).unwrap();
    source
        .insert(input(20, block_key, (0..300).rev().collect()))
        .unwrap();

    assert_eq!(
        source.header(20, &inline_key).unwrap().codec,
        PostingCodec::Inline
    );
    assert_eq!(
        source.decode_document_ids(20, &inline_key).unwrap(),
        vec![1, 7]
    );
    let header = source.header(20, &block_key).unwrap();
    assert_eq!(header.codec, PostingCodec::DeltaVarByte);
    assert!(header.block_count >= 2);
    assert_eq!(
        source.decode_document_ids(20, &block_key).unwrap(),
        (0..300).collect::<Vec<_>>()
    );
}

#[test]
fn packed_key_equality_is_required_after_digest_lookup() {
    let actual = key(11);
    let forged = WitnessKey {
        packed: actual.packed + 1,
        jamhash: actual.jamhash,
    };
    let mut source = InMemoryPostingSource::new(vec!["target".to_string()]);
    source.insert(input(20, actual, vec![0])).unwrap();

    assert!(source.header(20, &forged).is_none());
    let results = search(
        &source,
        &[TieredQueryWitness::new(20, query(forged, 0, &[0]))],
        &default_config(),
    )
    .unwrap();
    assert!(results.is_empty());
    assert_eq!(
        source.stats(),
        PostingDecodeStats {
            document_payload_decodes: 0,
            position_payload_decodes: 0,
        }
    );
}

#[test]
fn suppressed_common_header_skips_payload_decode() {
    let witness = key(21);
    let mut source = InMemoryPostingSource::new(vec!["common".to_string()]);
    source
        .insert(input(20, witness, vec![0]).mark_suppressed())
        .unwrap();
    let results = search(
        &source,
        &[TieredQueryWitness::new(20, query(witness, 0, &[0]))],
        &default_config(),
    )
    .unwrap();
    assert!(results.is_empty());
    assert_eq!(source.stats().document_payload_decodes, 0);
}

#[test]
fn top_k_order_is_deterministic_for_equal_evidence() {
    let first = key(31);
    let second = key(32);
    let mut source = InMemoryPostingSource::new(vec!["zeta".to_string(), "alpha".to_string()]);
    source.insert(input(20, first, vec![0])).unwrap();
    source.insert(input(20, second, vec![1])).unwrap();
    let witnesses = vec![
        TieredQueryWitness::new(20, query(second, 1, &[0])),
        TieredQueryWitness::new(20, query(first, 0, &[0])),
    ];
    let mut config = default_config();
    config.top_k = 2;
    config.accumulator_capacity = 2;
    let left = search(&source, &witnesses, &config).unwrap();
    let right = search(&source, &witnesses, &config).unwrap();
    assert_eq!(left, right);
    assert_eq!(
        left.iter()
            .map(|candidate| candidate.metagenome_id.as_str())
            .collect::<Vec<_>>(),
        vec!["alpha", "zeta"]
    );
}

#[test]
fn later_tiers_union_without_removing_earlier_candidates() {
    let dense = key(41);
    let sparse = key(42);
    let mut source =
        InMemoryPostingSource::new(vec!["dense-hit".to_string(), "sparse-hit".to_string()]);
    source.insert(input(20, dense, vec![0])).unwrap();
    source.insert(input(100, sparse, vec![1])).unwrap();
    let mut config = default_config();
    config.top_k = 2;
    config.accumulator_capacity = 2;
    let dense_only = search(
        &source,
        &[TieredQueryWitness::new(20, query(dense, 0, &[0]))],
        &config,
    )
    .unwrap();
    let union = search(
        &source,
        &[
            TieredQueryWitness::new(20, query(dense, 0, &[0])),
            TieredQueryWitness::new(100, query(sparse, 100, &[1])),
        ],
        &config,
    )
    .unwrap();
    assert!(
        union
            .iter()
            .any(|candidate| candidate.metagenome_id == dense_only[0].metagenome_id)
    );
    assert!(
        union
            .iter()
            .any(|candidate| candidate.metagenome_id == "sparse-hit")
    );
}

#[test]
fn nested_tiers_deduplicate_one_physical_query_occurrence() {
    let witness = key(43);
    let mut source = InMemoryPostingSource::new(vec!["target".to_string()]);
    source.insert(input(20, witness, vec![0])).unwrap();
    source.insert(input(100, witness, vec![0])).unwrap();
    let results = search(
        &source,
        &[
            TieredQueryWitness::new(20, query(witness, 8, &[0])),
            TieredQueryWitness::new(100, query(witness, 8, &[1])),
        ],
        &default_config(),
    )
    .unwrap();
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].window_evidence.total_shared_witness_count, 1);
    assert_eq!(results[0].supported_query_windows, 2);
    assert_eq!(results[0].candidate_tier, 100);
}

#[test]
fn window_evidence_and_handoff_details_are_bounded() {
    let witness = key(51);
    let mut source = InMemoryPostingSource::new(vec!["positioned".to_string()]);
    source
        .insert(
            input(20, witness, vec![0])
                .with_positions(vec![
                    PositionPosting {
                        metagenome_id: 0,
                        contig_id: 2,
                        position: 200,
                        reverse: false,
                    },
                    PositionPosting {
                        metagenome_id: 0,
                        contig_id: 2,
                        position: 100,
                        reverse: true,
                    },
                ])
                .with_position_cap(8),
        )
        .unwrap();
    let mut config = default_config();
    config.handoff_mode = WitnessHandoffMode::PositionBearing;
    config.max_shared_witnesses_per_candidate = 1;
    config.max_positional_witnesses_per_candidate = 1;
    let results = search(
        &source,
        &[TieredQueryWitness::new(20, query(witness, 4, &[0, 1, 3]))],
        &config,
    )
    .unwrap();
    let candidate = &results[0];
    assert_eq!(candidate.supported_query_windows, 3);
    assert_eq!(candidate.window_evidence.longest_supported_window_run, 2);
    assert_eq!(candidate.shared_witnesses.len(), 1);
    assert_eq!(candidate.positional_witnesses.len(), 1);

    config.handoff_mode = WitnessHandoffMode::SampleOnly;
    let sample_only = search(
        &source,
        &[TieredQueryWitness::new(20, query(witness, 4, &[0]))],
        &config,
    )
    .unwrap();
    assert!(sample_only[0].positional_witnesses.is_empty());
}

#[test]
fn common_evidence_requires_multiple_windows_when_configured() {
    let witness = key(61);
    let mut source = InMemoryPostingSource::new(vec!["common".to_string(); 3]);
    source
        .insert(input(20, witness, vec![0, 1, 2]).mark_common())
        .unwrap();
    let mut config = default_config();
    config.common_query_window_requirement = 2;
    config.rare_max_document_frequency = 1;
    config.moderately_common_max_document_frequency = 2;
    let one_window = search(
        &source,
        &[TieredQueryWitness::new(20, query(witness, 0, &[0]))],
        &config,
    )
    .unwrap();
    assert!(one_window.is_empty());
    let two_windows = search(
        &source,
        &[TieredQueryWitness::new(20, query(witness, 0, &[0, 1]))],
        &config,
    )
    .unwrap();
    assert_eq!(two_windows.len(), 3);
    assert!(
        two_windows
            .iter()
            .all(|candidate| candidate.rare_shared_witnesses == 0)
    );
}

#[test]
fn exact_two_pass_top_k_keeps_a_late_gaining_candidate() {
    let first = key(71);
    let second = key(72);
    let mut source = InMemoryPostingSource::new(vec!["early".to_string(), "winner".to_string()]);
    source.insert(input(20, first, vec![0, 1])).unwrap();
    source.insert(input(20, second, vec![1])).unwrap();
    let mut config = default_config();
    config.top_k = 1;
    config.accumulator_capacity = 1;
    let results = search(
        &source,
        &[
            TieredQueryWitness::new(20, query(first, 0, &[0])),
            TieredQueryWitness::new(20, query(second, 1, &[1])),
        ],
        &config,
    )
    .unwrap();
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].metagenome_id, "winner");
    assert_eq!(results[0].window_evidence.total_shared_witness_count, 2);
}

#[test]
fn overlapping_windows_count_one_exact_witness() {
    let witness = key(73);
    let mut source = InMemoryPostingSource::new(vec![
        "target".to_string(),
        "other-a".to_string(),
        "other-b".to_string(),
    ]);
    source.insert(input(20, witness, vec![0])).unwrap();
    let single = search(
        &source,
        &[TieredQueryWitness::new(20, query(witness, 12, &[3]))],
        &default_config(),
    )
    .unwrap();
    let results = search(
        &source,
        &[TieredQueryWitness::new(20, query(witness, 12, &[3, 4, 5]))],
        &default_config(),
    )
    .unwrap();
    let candidate = &results[0];
    assert_eq!(candidate.rare_shared_witnesses, 1);
    assert_eq!(candidate.common_shared_witnesses, 0);
    assert_eq!(candidate.window_evidence.total_shared_witness_count, 1);
    assert_eq!(candidate.supported_query_windows, 3);
    assert_eq!(candidate.estimated_query_containment, 1.0);
    assert_eq!(
        candidate.weighted_witness_sum,
        single[0].weighted_witness_sum
    );
}

#[test]
fn positions_are_decoded_only_for_pass_two_selected_candidates() {
    let witness = key(74);
    let mut source = InMemoryPostingSource::new(vec!["zeta".to_string(), "alpha".to_string()]);
    let positions = vec![
        PositionPosting {
            metagenome_id: 0,
            contig_id: 1,
            position: 10,
            reverse: false,
        },
        PositionPosting {
            metagenome_id: 1,
            contig_id: 1,
            position: 20,
            reverse: false,
        },
    ];
    source
        .insert(input(20, witness, vec![0, 1]).with_positions(positions))
        .unwrap();
    let mut config = default_config();
    config.top_k = 1;
    config.accumulator_capacity = 1;
    config.handoff_mode = WitnessHandoffMode::PositionBearing;
    let results = search(
        &source,
        &[TieredQueryWitness::new(20, query(witness, 2, &[0]))],
        &config,
    )
    .unwrap();
    assert_eq!(results[0].metagenome_id, "alpha");
    assert_eq!(results[0].positional_witnesses.len(), 1);
    assert_eq!(source.stats().position_payload_decodes, 1);
}

#[test]
fn positional_handoff_is_created_from_the_exact_shared_detail() {
    let witness = key(75);
    let mut source = InMemoryPostingSource::new(vec!["target".to_string()]);
    source
        .insert(
            input(20, witness, vec![0]).with_positions(vec![PositionPosting {
                metagenome_id: 0,
                contig_id: 4,
                position: 40,
                reverse: true,
            }]),
        )
        .unwrap();
    let mut config = default_config();
    config.handoff_mode = WitnessHandoffMode::PositionBearing;
    let results = search(
        &source,
        &[TieredQueryWitness::new(20, query(witness, 99, &[7]))],
        &config,
    )
    .unwrap();
    let candidate = &results[0];
    assert_eq!(candidate.shared_witnesses.len(), 1);
    assert_eq!(candidate.positional_witnesses.len(), 1);
    assert_eq!(candidate.shared_witnesses[0].key, witness);
    assert_eq!(candidate.positional_witnesses[0].witness.key, witness);
    assert_eq!(candidate.positional_witnesses[0].witness.query_position, 99);
    assert_eq!(candidate.positional_witnesses[0].position, 40);
}
