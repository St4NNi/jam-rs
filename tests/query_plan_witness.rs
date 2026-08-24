use jam_rs::router::{HashAlgorithmId, WITNESS_K, WitnessScheme};
use jam_rs::trace::contracts::WitnessPlanningRequest;
use jam_rs::trace::query_plan::{QueryPlanError, estimate_zero_witness_risk, plan_witness_tier};

fn scheme() -> WitnessScheme {
    WitnessScheme {
        scheme_id: 1,
        k: WITNESS_K,
        base_scale: 20,
        available_scales: vec![20, 50, 100, 200, 500],
        hash_id: HashAlgorithmId::JamhashU64V1,
        zero_excluded: true,
    }
}

#[test]
fn planner_uses_the_requested_trace_length_and_identity() {
    let request = WitnessPlanningRequest {
        min_trace_length: 1_000,
        target_identity: 0.99,
        max_zero_witness_risk: 0.01,
        strict: true,
    };
    let plan = plan_witness_tier(request, &scheme()).unwrap();
    let selected = plan.selected_witness_tier.unwrap();
    assert!(estimate_zero_witness_risk(1_000, 0.99, selected) <= 0.01);
    assert!(
        scheme()
            .available_scales
            .iter()
            .filter(|scale| **scale > selected)
            .all(|scale| estimate_zero_witness_risk(1_000, 0.99, *scale) > 0.01)
    );
}

#[test]
fn unavailable_short_trace_risk_warns_or_fails_strictly() {
    let request = WitnessPlanningRequest {
        min_trace_length: 80,
        target_identity: 0.90,
        max_zero_witness_risk: 0.01,
        strict: false,
    };
    let normal = plan_witness_tier(request, &scheme()).unwrap();
    assert!(normal.selected_witness_tier.is_none());
    assert!(normal.warning.is_some());

    let strict = WitnessPlanningRequest {
        strict: true,
        ..request
    };
    assert!(matches!(
        plan_witness_tier(strict, &scheme()),
        Err(QueryPlanError::NoTierMeetsRisk { .. })
    ));
}
