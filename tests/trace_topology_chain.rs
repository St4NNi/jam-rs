use std::collections::BTreeSet;

use jam_rs::trace::anchors::Anchor;
use jam_rs::trace::chain::{ChainConfig, chain_anchors};
use jam_rs::trace::model::{CoordinateModel, Strand, TopologyRequested};

const QUERY_LENGTH: u64 = 100;

fn anchor(query_position: u64, target_position: u64) -> Anchor {
    Anchor {
        query_position,
        target_position,
        contig_id: 0,
        strand: Strand::Forward,
        k: 3,
        hash: query_position + 1,
        canonical_kmer: query_position + 101,
        query_reverse: false,
        target_reverse: false,
    }
}

fn config(coordinate_model: CoordinateModel, min_anchors: u32) -> ChainConfig {
    ChainConfig {
        max_chains: 1,
        min_anchors,
        max_predecessors: 32,
        max_query_gap: QUERY_LENGTH,
        max_target_gap: QUERY_LENGTH,
        gap_penalty: 1,
        coordinate_model,
    }
}

fn crossing_anchors() -> [Anchor; 3] {
    [anchor(90, 190), anchor(2, 202), anchor(12, 212)]
}

fn covered_bases(chains: &[jam_rs::trace::chain::AnchorChain]) -> BTreeSet<u64> {
    chains
        .iter()
        .flat_map(|chain| {
            chain
                .query_segments
                .iter()
                .flat_map(|interval| interval.start..interval.end)
        })
        .collect()
}

fn rotate_position(position: u64, rotation: u64) -> u64 {
    (position + QUERY_LENGTH - rotation) % QUERY_LENGTH
}

#[test]
fn linear_model_never_duplicates_anchors_or_wraps() {
    let chains = chain_anchors(
        &crossing_anchors(),
        QUERY_LENGTH,
        config(CoordinateModel::Linear, 2),
    )
    .unwrap();

    assert!(!chains.is_empty());
    assert!(chains.iter().all(|chain| !chain.origin_crossing));
    assert!(chains.iter().all(|chain| chain.query_segments.len() == 1));
    assert!(chains.iter().all(|chain| {
        chain
            .anchors
            .windows(2)
            .all(|window| window[0].query_position < window[1].query_position)
    }));
    assert_eq!(covered_bases(&chains).len(), 13);
}

#[test]
fn wrap_model_can_chain_across_the_origin_with_bounded_span() {
    let chains = chain_anchors(
        &crossing_anchors(),
        QUERY_LENGTH,
        config(CoordinateModel::Wrap, 3),
    )
    .unwrap();
    assert_eq!(chains.len(), 1);
    let chain = &chains[0];
    assert!(chain.origin_crossing);
    assert_eq!(chain.anchors.len(), 3);
    assert_eq!(chain.linear_query_start, 90);
    assert_eq!(chain.linear_query_end, 115);
    assert_eq!(chain.linear_query_end - chain.linear_query_start, 25);
    assert!(chain.linear_query_end - chain.linear_query_start <= QUERY_LENGTH);
    assert_eq!(chain.query_segments.len(), 2);
}

#[test]
fn wrapped_chain_support_is_invariant_to_query_origin_rotation() {
    let original = crossing_anchors();
    let original_chains =
        chain_anchors(&original, QUERY_LENGTH, config(CoordinateModel::Wrap, 3)).unwrap();
    let original_support = covered_bases(&original_chains);

    for rotation in [10, 20, 55] {
        let rotated = original.map(|mut value| {
            value.query_position = rotate_position(value.query_position, rotation);
            value
        });
        let chains =
            chain_anchors(&rotated, QUERY_LENGTH, config(CoordinateModel::Wrap, 3)).unwrap();
        let expected = original_support
            .iter()
            .map(|&position| rotate_position(position, rotation))
            .collect::<BTreeSet<_>>();
        assert_eq!(covered_bases(&chains), expected, "rotation={rotation}");
    }
}

#[test]
fn auto_and_unknown_callers_can_evaluate_a_wrap_alternative_explicitly() {
    let requested = [TopologyRequested::Auto, TopologyRequested::Unknown];
    for topology in requested {
        // The primary display model for both requests is linear. A caller
        // evaluating the alternative model must opt into duplication
        // explicitly; chaining does not infer biological topology.
        let linear = chain_anchors(
            &crossing_anchors(),
            QUERY_LENGTH,
            config(CoordinateModel::Undetermined, 3),
        )
        .unwrap();
        assert!(linear.iter().all(|chain| !chain.origin_crossing));

        let wrapped = chain_anchors(
            &crossing_anchors(),
            QUERY_LENGTH,
            config(CoordinateModel::Wrap, 3),
        )
        .unwrap();
        assert!(wrapped.iter().any(|chain| chain.origin_crossing));
        assert!(matches!(
            topology,
            TopologyRequested::Auto | TopologyRequested::Unknown
        ));
    }
}

#[test]
fn terminal_repeat_anchors_do_not_force_a_wrap() {
    let anchors = [
        // Head and tail occurrences are ordered independently on the target;
        // joining them across the query origin would contradict the target.
        anchor(2, 100),
        anchor(8, 106),
        anchor(88, 500),
        anchor(94, 506),
    ];
    let chains = chain_anchors(&anchors, QUERY_LENGTH, config(CoordinateModel::Wrap, 2)).unwrap();
    assert!(!chains.is_empty());
    assert!(chains.iter().all(|chain| !chain.origin_crossing));
    assert!(chains.iter().all(|chain| chain.query_segments.len() == 1));
}
