//! Query-length-aware planning for nested k=21 witness tiers.
//!
//! The estimate here is deliberately an engineering model.  It describes the
//! tier selected from the archive's available scales; it is not a biological
//! sensitivity guarantee.  The model follows the frozen contract:
//!
//! ```text
//! n = max(0, min_trace_length - k + 1)
//! survival = target_identity ^ k
//! lambda = n * survival / scale
//! risk = exp(-lambda)
//! ```

use thiserror::Error;

use crate::router::{RouterContractError, WITNESS_K, WitnessScheme};
use crate::trace::contracts::{WitnessPlanningRequest, WitnessTierPlan};

/// Human-readable assumptions persisted with every tier plan.
pub const MODEL_ASSUMPTIONS: [&str; 4] = [
    "eligible k-mer starts are approximated as independent witness trials",
    "per-k-mer survival is target_identity raised to the witness k",
    "n is max(0, min_trace_length - k + 1)",
    "zero-witness risk is estimated as exp(-n * survival / scale)",
];

/// Planning failures that should be surfaced before candidate routing begins.
#[derive(Debug, Error, PartialEq)]
pub enum QueryPlanError {
    #[error(transparent)]
    Contract(#[from] RouterContractError),
    #[error("minimum trace identity must be finite and between 0 and 1")]
    InvalidIdentity,
    #[error("maximum zero-witness risk must be finite and between 0 and 1")]
    InvalidRisk,
    #[error(
        "no witness tier meets requested zero-witness risk {requested_risk}; available scales: {available_scales:?}"
    )]
    NoTierMeetsRisk {
        requested_risk: f64,
        available_scales: Vec<u32>,
    },
}

/// Estimate the zero-witness risk for one scale using the fixed k=21 model.
#[must_use]
pub fn estimate_zero_witness_risk(min_trace_length: u64, target_identity: f64, scale: u32) -> f64 {
    estimate_zero_witness_risk_for_k(min_trace_length, target_identity, WITNESS_K, scale)
}

/// General form of [`estimate_zero_witness_risk`] used by deterministic tests
/// and by callers validating a future scheme before construction.
#[must_use]
pub fn estimate_zero_witness_risk_for_k(
    min_trace_length: u64,
    target_identity: f64,
    k: u8,
    scale: u32,
) -> f64 {
    if scale == 0 {
        return f64::NAN;
    }
    let n = min_trace_length.saturating_sub(u64::from(k).saturating_sub(1));
    let survival = target_identity.powf(f64::from(k));
    let lambda = (n as f64) * survival / f64::from(scale);
    (-lambda).exp()
}

/// Choose the largest (coarsest) available scale whose modeled risk is within
/// the requested bound.  `available_scales` must be strictly increasing, as
/// required by [`WitnessScheme::validate`].
pub fn select_largest_acceptable_scale(
    min_trace_length: u64,
    target_identity: f64,
    max_zero_witness_risk: f64,
    k: u8,
    available_scales: &[u32],
) -> Result<Option<u32>, QueryPlanError> {
    validate_request_values(target_identity, max_zero_witness_risk)?;
    if k == 0 || available_scales.is_empty() {
        return Ok(None);
    }
    if available_scales.windows(2).any(|pair| pair[0] >= pair[1]) {
        return Err(QueryPlanError::Contract(RouterContractError::InvalidScales));
    }

    let mut selected = None;
    for &scale in available_scales {
        if scale == 0 {
            return Err(QueryPlanError::Contract(RouterContractError::InvalidScales));
        }
        let risk = estimate_zero_witness_risk_for_k(min_trace_length, target_identity, k, scale);
        if risk <= max_zero_witness_risk {
            selected = Some(scale);
        }
    }
    Ok(selected)
}

/// Build a plan from the archive's validated witness scheme.
pub fn plan_witness_tier(
    request: WitnessPlanningRequest,
    scheme: &WitnessScheme,
) -> Result<WitnessTierPlan, QueryPlanError> {
    scheme.validate()?;
    validate_request_values(request.target_identity, request.max_zero_witness_risk)?;

    let selected = select_largest_acceptable_scale(
        request.min_trace_length,
        request.target_identity,
        request.max_zero_witness_risk,
        scheme.k,
        &scheme.available_scales,
    )?;

    let warning = selected.is_none().then(|| {
        format!(
            "available witness tiers cannot meet requested zero-witness risk {:.6} for a {} bp trace at identity {:.6}",
            request.max_zero_witness_risk, request.min_trace_length, request.target_identity
        )
    });
    if request.strict && warning.is_some() {
        return Err(QueryPlanError::NoTierMeetsRisk {
            requested_risk: request.max_zero_witness_risk,
            available_scales: scheme.available_scales.clone(),
        });
    }

    let estimated_zero_witness_risk = selected.map_or_else(
        || {
            // Keep the plan numerically useful in normal mode when no tier
            // meets the requested bound.  The densest available tier is the
            // best archive-supported estimate, and the warning records that
            // it is outside the requested bound.
            scheme.available_scales.first().map_or(1.0, |scale| {
                estimate_zero_witness_risk(
                    request.min_trace_length,
                    request.target_identity,
                    *scale,
                )
            })
        },
        |scale| {
            estimate_zero_witness_risk(request.min_trace_length, request.target_identity, scale)
        },
    );

    let available_denser_tiers = match selected {
        Some(scale) => scheme
            .available_scales
            .iter()
            .copied()
            .filter(|candidate| *candidate < scale)
            .collect(),
        None => scheme.available_scales.clone(),
    };

    Ok(WitnessTierPlan {
        model_assumptions: MODEL_ASSUMPTIONS
            .iter()
            .map(|assumption| (*assumption).to_string())
            .collect(),
        requested_min_trace_length: request.min_trace_length,
        requested_identity: request.target_identity,
        selected_witness_tier: selected,
        estimated_zero_witness_risk,
        available_denser_tiers,
        warning,
    })
}

/// Alias emphasizing that the planner is scheme-driven.
pub fn plan_for_scheme(
    scheme: &WitnessScheme,
    request: WitnessPlanningRequest,
) -> Result<WitnessTierPlan, QueryPlanError> {
    plan_witness_tier(request, scheme)
}

/// Return the fixed query-window candidates used by the short-trace matrix.
#[must_use]
pub fn standard_query_window_sizes() -> &'static [u32; 4] {
    &crate::router::witness::DEFAULT_QUERY_WINDOW_SIZES
}

fn validate_request_values(
    target_identity: f64,
    max_zero_witness_risk: f64,
) -> Result<(), QueryPlanError> {
    if !target_identity.is_finite() || !(0.0..=1.0).contains(&target_identity) {
        return Err(QueryPlanError::InvalidIdentity);
    }
    if !max_zero_witness_risk.is_finite() || !(0.0..=1.0).contains(&max_zero_witness_risk) {
        return Err(QueryPlanError::InvalidRisk);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::router::HashAlgorithmId;

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
    fn risk_formula_matches_contract() {
        let risk = estimate_zero_witness_risk(160, 0.99, 20);
        let n = (160_u64 - 21 + 1) as f64;
        let expected = (-(n * 0.99_f64.powf(21.0) / 20.0)).exp();
        assert!((risk - expected).abs() < 1e-15);
    }

    #[test]
    fn planner_selects_largest_acceptable_scale() {
        let request = WitnessPlanningRequest {
            min_trace_length: 1_000,
            target_identity: 0.99,
            max_zero_witness_risk: 0.01,
            strict: true,
        };
        let plan = plan_witness_tier(request, &scheme()).unwrap();
        assert!(plan.selected_witness_tier.is_some());
        let selected = plan.selected_witness_tier.unwrap();
        assert!(
            scheme()
                .available_scales
                .iter()
                .filter(|scale| **scale <= selected)
                .all(|scale| { estimate_zero_witness_risk(1_000, 0.99, *scale) <= 0.01 })
        );
        assert!(
            plan.available_denser_tiers
                .iter()
                .all(|scale| *scale < selected)
        );
    }

    #[test]
    fn short_trace_can_fail_strictly_or_warn_normally() {
        let strict = WitnessPlanningRequest {
            min_trace_length: 80,
            target_identity: 0.90,
            max_zero_witness_risk: 0.01,
            strict: true,
        };
        assert!(matches!(
            plan_witness_tier(strict, &scheme()),
            Err(QueryPlanError::NoTierMeetsRisk { .. })
        ));

        let normal = WitnessPlanningRequest {
            strict: false,
            ..strict
        };
        let plan = plan_witness_tier(normal, &scheme()).unwrap();
        assert!(plan.selected_witness_tier.is_none());
        assert!(plan.warning.is_some());
        assert_eq!(plan.available_denser_tiers, scheme().available_scales);
    }

    #[test]
    fn zero_length_and_invalid_values_are_deterministic() {
        assert_eq!(estimate_zero_witness_risk(0, 0.99, 20), 1.0);
        let invalid_identity = WitnessPlanningRequest {
            min_trace_length: 128,
            target_identity: f64::NAN,
            max_zero_witness_risk: 0.01,
            strict: false,
        };
        assert_eq!(
            plan_witness_tier(invalid_identity, &scheme()),
            Err(QueryPlanError::InvalidIdentity)
        );
    }
}
