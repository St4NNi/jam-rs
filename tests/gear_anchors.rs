use jam_rs::trace::seeds::gear::{
    FragmentationMode, GearAnchorConfig, GearAnchorIndex, GearAnchorScheme, GearConfig,
    GearTableKind, canonical_anchor_at, fragment_bytes, fragment_sequence, gear_selected_anchors,
    verify_anchor_candidate, verify_anchor_sequences,
};

fn dna(length: usize) -> Vec<u8> {
    (0..length)
        .map(|index| b"ACGT"[(index.wrapping_mul(31) ^ (index / 7)) & 3])
        .collect()
}

#[test]
fn gear_selected_anchors_have_only_allowed_spans() {
    let sequence = dna(8_000);
    let gear = GearConfig::fine(GearTableKind::PackedFourBase, 0x5eed);
    for config in [
        GearAnchorConfig::k31(),
        GearAnchorConfig::k21_rescue(),
        GearAnchorConfig::spaced31_weight21(),
    ] {
        let anchors =
            gear_selected_anchors(&sequence, gear, FragmentationMode::Forward, config).unwrap();
        assert!(anchors.iter().all(|anchor| anchor.span >= 21));
        assert!(anchors.iter().all(|anchor| anchor.informative_bases >= 21));
        assert!(
            anchors
                .windows(2)
                .all(|window| window[0].position <= window[1].position)
        );
        assert!(anchors.iter().all(|anchor| anchor.hash != 0));
    }
}

#[test]
fn anchor_hash_collision_is_filtered_by_exact_packed_key() {
    let sequence = dna(6_000);
    let config = GearConfig::fine(GearTableKind::SingleBase, 17);
    let anchors = gear_selected_anchors(
        &sequence,
        config,
        FragmentationMode::Forward,
        GearAnchorConfig::k31(),
    )
    .unwrap();
    let query = anchors.iter().copied().next().expect("anchor fixture");
    let mut collision = query;
    collision.packed_selected ^= 1;
    assert!(!verify_anchor_candidate(query, collision));
    let mut index = GearAnchorIndex::default();
    index.insert(collision);
    assert!(index.candidates(query.hash).len() == 1);
    assert!(index.verified_candidates(query).is_empty());
}

#[test]
fn canonical_anchor_verification_accepts_forward_and_reverse_targets() {
    let sequence = dna(2_000);
    let reverse = sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect::<Vec<_>>();
    let anchors = gear_selected_anchors(
        &sequence,
        GearConfig::fine(GearTableKind::Dinucleotide, 2),
        FragmentationMode::Forward,
        GearAnchorConfig {
            scheme: GearAnchorScheme::CanonicalK31,
            max_anchors: 1,
        },
    )
    .unwrap();
    let anchor = anchors[0];
    assert!(verify_anchor_sequences(
        &sequence,
        anchor,
        &sequence,
        anchor.position,
        false,
    ));
    let reverse_position = sequence.len() as u64 - anchor.position - u64::from(anchor.span);
    assert!(verify_anchor_sequences(
        &sequence,
        anchor,
        &reverse,
        reverse_position,
        true,
    ));
}

#[test]
fn reverse_fragment_anchor_uses_forward_source_orientation_and_coordinates() {
    let sequence = dna(8_000);
    let gear = GearConfig::fine(GearTableKind::SingleBase, 0x5eed);
    let records = fragment_sequence(&sequence, 0, gear, FragmentationMode::BothStrands).unwrap();
    let anchors = gear_selected_anchors(
        &sequence,
        gear,
        FragmentationMode::BothStrands,
        GearAnchorConfig::k31(),
    )
    .unwrap();
    let mut reverse_checked = false;
    for anchor in anchors {
        let record = records[usize::try_from(anchor.fragment_index).unwrap()];
        if record.orientation != jam_rs::trace::seeds::gear::FragmentOrientation::Reverse {
            continue;
        }
        let expected = canonical_anchor_at(&sequence, anchor.position, anchor.scheme).unwrap();
        assert_eq!(
            anchor.position,
            record.start + u64::from(record.length) - u64::from(anchor.span)
        );
        assert_eq!(anchor.reverse, expected.reverse);
        let oriented = fragment_bytes(&sequence, record).unwrap();
        let local = canonical_anchor_at(&oriented, 0, anchor.scheme).unwrap();
        assert_eq!(anchor.reverse, local.reverse ^ true);
        assert!(verify_anchor_sequences(
            &oriented,
            local,
            &sequence,
            anchor.position,
            true,
        ));
        reverse_checked = true;
        break;
    }
    assert!(reverse_checked, "expected a reverse Gear fragment anchor");
}
