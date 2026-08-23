use jam_rs::trace::anchors::Anchor;
use jam_rs::trace::chain::{ChainConfig, ChainError, chain_anchors};
use jam_rs::trace::model::{BaseInterval, Strand};

fn anchor(query_position: u64, target_position: u64, strand: Strand) -> Anchor {
    Anchor {
        query_position,
        target_position,
        contig_id: 0,
        strand,
        k: 3,
        hash: query_position + 1,
        canonical_kmer: query_position + 9,
        query_reverse: false,
        target_reverse: false,
    }
}

fn config() -> ChainConfig {
    ChainConfig {
        max_chains: 4,
        min_anchors: 2,
        max_predecessors: 32,
        max_query_gap: 100,
        max_target_gap: 100,
        gap_penalty: 1,
    }
}

#[test]
fn forward_and_reverse_paths_are_orientation_aware() {
    let forward = [
        anchor(5, 105, Strand::Forward),
        anchor(15, 115, Strand::Forward),
        anchor(25, 125, Strand::Forward),
    ];
    let chains = chain_anchors(&forward, 100, config()).unwrap();
    assert_eq!(chains.len(), 1);
    assert_eq!(chains[0].anchors.len(), 3);

    let reverse = [
        anchor(5, 125, Strand::Reverse),
        anchor(15, 115, Strand::Reverse),
        anchor(25, 105, Strand::Reverse),
    ];
    let chains = chain_anchors(&reverse, 100, config()).unwrap();
    assert_eq!(chains.len(), 1);
    assert_eq!(chains[0].strand, Strand::Reverse);
}

#[test]
fn origin_crossing_is_reported_as_two_non_wrapping_segments() {
    let anchors = [
        anchor(90, 190, Strand::Forward),
        anchor(2, 202, Strand::Forward),
        anchor(12, 212, Strand::Forward),
    ];
    let chains = chain_anchors(&anchors, 100, config()).unwrap();
    assert_eq!(chains.len(), 1);
    assert!(chains[0].origin_crossing);
    assert_eq!(
        chains[0].query_segments,
        vec![
            BaseInterval::new(90, 100).unwrap(),
            BaseInterval::new(0, 15).unwrap()
        ]
    );
}

#[test]
fn bounded_chaining_rejects_invalid_limits_and_anchor_coordinates() {
    let mut invalid = config();
    invalid.max_predecessors = 0;
    assert!(matches!(
        chain_anchors(&[], 100, invalid),
        Err(ChainError::ZeroLimit {
            field: "max_predecessors"
        })
    ));
    assert!(matches!(
        chain_anchors(&[anchor(99, 1, Strand::Forward)], 100, config()),
        Err(ChainError::AnchorOutsideQuery { .. })
    ));
}
