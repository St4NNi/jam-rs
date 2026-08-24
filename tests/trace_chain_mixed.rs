use jam_rs::trace::anchors::Anchor;
use jam_rs::trace::chain::{
    AnchorClass, ChainConfig, WeightedAnchor, chain_anchors, chain_weighted_anchors,
};
use jam_rs::trace::model::{CoordinateModel, Strand};

fn anchor(query_position: u64, target_position: u64, k: u8) -> Anchor {
    Anchor {
        query_position,
        target_position,
        contig_id: 0,
        strand: Strand::Forward,
        k,
        hash: query_position + 11,
        canonical_kmer: query_position + 101,
        query_reverse: false,
        target_reverse: false,
    }
}

fn config() -> ChainConfig {
    ChainConfig {
        max_chains: 4,
        min_anchors: 3,
        max_predecessors: 32,
        max_query_gap: 100,
        max_target_gap: 100,
        gap_penalty: 1,
        coordinate_model: CoordinateModel::Linear,
    }
}

#[test]
fn legacy_k31_and_k21_anchors_join_one_mixed_chain() {
    let anchors = [
        anchor(10, 110, 31),
        anchor(50, 150, 21),
        anchor(80, 180, 21),
    ];
    let chains = chain_anchors(&anchors, 200, config()).unwrap();
    assert_eq!(chains.len(), 1);
    assert_eq!(chains[0].anchors.len(), 3);
}

#[test]
fn explicit_weights_accept_multiple_lower_specificity_anchors() {
    let anchors = [
        WeightedAnchor::new(anchor(10, 110, 21), AnchorClass::SpecificK21),
        WeightedAnchor::new(anchor(50, 150, 21), AnchorClass::SpecificK21),
        WeightedAnchor::new(anchor(80, 180, 21), AnchorClass::SpecificSpaced),
    ];
    let chains = chain_weighted_anchors(&anchors, 200, config()).unwrap();
    assert_eq!(chains.len(), 1);
    assert!(!chains[0].has_high_specificity_anchor());
    assert_eq!(chains[0].anchor_count(), 3);
}

#[test]
fn one_repetitive_rescue_hit_cannot_create_a_chain() {
    let anchors = [WeightedAnchor::new(
        anchor(10, 110, 21),
        AnchorClass::Repetitive,
    )];
    let mut one_anchor = config();
    one_anchor.min_anchors = 1;
    assert!(
        chain_weighted_anchors(&anchors, 200, one_anchor)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn one_high_specificity_anchor_is_eligible_when_the_profile_allows_it() {
    let anchors = [WeightedAnchor::new(
        anchor(10, 110, 31),
        AnchorClass::SpecificK31,
    )];
    let mut one_anchor = config();
    one_anchor.min_anchors = 3;
    let chains = chain_weighted_anchors(&anchors, 200, one_anchor).unwrap();
    assert_eq!(chains.len(), 1);
    assert!(chains[0].has_high_specificity_anchor());
}

#[test]
fn weighted_score_is_class_aware_and_input_order_independent() {
    let anchors = [
        WeightedAnchor::new(anchor(10, 110, 21), AnchorClass::SpecificK21),
        WeightedAnchor::new(anchor(50, 150, 21), AnchorClass::SpecificK21),
        WeightedAnchor::new(anchor(80, 180, 21), AnchorClass::SpecificPaired),
    ];
    let mut reversed = anchors;
    reversed.reverse();
    let left = chain_weighted_anchors(&anchors, 200, config()).unwrap();
    let right = chain_weighted_anchors(&reversed, 200, config()).unwrap();
    assert_eq!(left, right);
    assert_eq!(left[0].score, 80 + 80 + 115);
}

#[test]
fn exact_gear_run_span_is_wider_than_anchor_k() {
    let gear_run =
        WeightedAnchor::new(anchor(10, 110, 31), AnchorClass::ExactGearRun).with_span(300);
    let chains = chain_weighted_anchors(&[gear_run], 1_000, config()).unwrap();
    assert_eq!(chains.len(), 1);
    assert_eq!(chains[0].query_interval.start, 10);
    assert_eq!(chains[0].query_interval.end, 310);
    assert_eq!(chains[0].target_interval.start, 110);
    assert_eq!(chains[0].target_interval.end, 410);
    assert_eq!(chains[0].anchors[0].span, 300);
}
