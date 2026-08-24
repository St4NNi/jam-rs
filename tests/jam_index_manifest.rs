use jam_rs::jam_index::{JamIndexManifest, JamIndexPart, ScreenSelectionPolicy};

fn checksum(byte: char) -> String {
    std::iter::repeat_n(byte, 64).collect()
}

fn part(part_id: u32) -> JamIndexPart {
    JamIndexPart {
        part_id,
        directory: format!("parts/part-{part_id:06}"),
        screen_file: "screen.jam".to_string(),
        data_file: "part.bin".to_string(),
        metagenome_count: 2,
        contig_count: 7,
        total_bases: 1_000,
        estimated_signature_count: 50,
        screen_jam_bytes: 800,
        contig_signature_bytes: 600,
        packed_sequence_bytes: 250,
        screen_sha256: checksum('a'),
        data_sha256: checksum('b'),
    }
}

#[test]
fn append_adds_only_new_independent_parts_and_recomputes_totals() {
    let mut manifest = JamIndexManifest::empty(ScreenSelectionPolicy::default_signatures());
    manifest
        .append_parts(checksum('c'), [part(0), part(1)])
        .unwrap();
    assert_eq!(manifest.next_part_id(), 2);
    assert_eq!(manifest.total_metagenomes, 4);
    assert_eq!(manifest.total_contigs, 14);
    assert_eq!(manifest.total_bases, 2_000);
    assert_eq!(manifest.estimated_signature_count, 100);
    manifest.validate().unwrap();

    manifest.append_parts(checksum('d'), [part(2)]).unwrap();
    assert_eq!(manifest.parts.len(), 3);
    assert_eq!(manifest.parts[0].directory, "parts/part-000000");
}

#[test]
fn manifest_rejects_path_escape_policy_drift_and_nonappend_parts() {
    let mut invalid = JamIndexManifest::empty(ScreenSelectionPolicy::default_signatures());
    let mut escaped = part(0);
    escaped.directory = "../outside".to_string();
    assert!(invalid.append_parts(checksum('c'), [escaped]).is_err());

    let mut wrong_id = JamIndexManifest::empty(ScreenSelectionPolicy::default_signatures());
    assert!(wrong_id.append_parts(checksum('c'), [part(1)]).is_err());

    let mut bad_policy = ScreenSelectionPolicy::default_signatures();
    bad_policy.k = 20;
    assert!(bad_policy.validate().is_err());
}

#[test]
fn default_and_smaller_policies_are_explicit_and_length_dependent() {
    let default = ScreenSelectionPolicy::default_signatures();
    let smaller = ScreenSelectionPolicy::smaller_signatures();
    default.validate().unwrap();
    smaller.validate().unwrap();
    assert_eq!(default.contig_budget.budget_for_bases(160), 8);
    assert_eq!(smaller.contig_budget.budget_for_bases(160), 4);
    assert!(
        default.estimated_signature_count(&[160, 1_000_000])
            > smaller.estimated_signature_count(&[160, 1_000_000])
    );
}
