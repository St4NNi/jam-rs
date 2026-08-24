use jam_rs::router::postings::{InMemoryPostingSource, PostingInput};
use jam_rs::router::search::{RouterSearchConfig, TieredQueryWitness, search};
use jam_rs::router::witness::extract_nested_witnesses;
use jam_rs::router::{HashAlgorithmId, WITNESS_K, WitnessScheme};
use jam_rs::trace::contracts::WitnessPlanningRequest;
use jam_rs::trace::query_plan::plan_witness_tier;

fn exact_sequence(length: usize) -> Vec<u8> {
    let bases = *b"ACGT";
    let mut state = 0x6a09_e667_f3bc_c909_u64;
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            bases[(state as usize) & 3]
        })
        .collect()
}

#[test]
fn exact_160_base_trace_is_routed_by_the_selected_microtrace_tier() {
    let scheme = WitnessScheme {
        scheme_id: 1,
        k: WITNESS_K,
        base_scale: 20,
        available_scales: vec![20, 50, 100, 200, 500],
        hash_id: HashAlgorithmId::JamhashU64V1,
        zero_excluded: true,
    };
    let plan = plan_witness_tier(
        WitnessPlanningRequest {
            min_trace_length: 160,
            target_identity: 1.0,
            max_zero_witness_risk: 0.01,
            strict: true,
        },
        &scheme,
    )
    .unwrap();
    assert_eq!(plan.selected_witness_tier, Some(20));

    let witnesses = extract_nested_witnesses(&exact_sequence(160), &scheme, 128)
        .unwrap()
        .at_scale(20)
        .unwrap();
    let shared = witnesses
        .first()
        .cloned()
        .expect("160 bp trace has a witness");
    let mut source = InMemoryPostingSource::new(vec!["short-trace-target".to_string()]);
    source
        .insert(PostingInput::new(20, shared.key, vec![0]))
        .unwrap();
    let candidates = search(
        &source,
        &witnesses
            .into_iter()
            .map(|witness| TieredQueryWitness::new(20, witness))
            .collect::<Vec<_>>(),
        &RouterSearchConfig::default(),
    )
    .unwrap();
    assert_eq!(candidates.len(), 1);
    assert_eq!(candidates[0].metagenome_id, "short-trace-target");
    assert_eq!(candidates[0].candidate_tier, 20);
}
