//! Deterministic malformed-input and coordinate properties for trace output.

use jam_rs::resource::ResourceMetrics;
use jam_rs::trace::coverage::{
    cigar_from_edit_script, parse_cigar, project_cigar, project_edit_script,
};
use jam_rs::trace::intervals::{
    circular_gap_complement, circular_union, covered_length, split_circular, union,
};
use jam_rs::trace::model::{BaseInterval, EditOperation, EditRun, TraceRecord, TraceRunFooter};
use jam_rs::trace::output::SCHEMA_VERSION;
use std::panic::{AssertUnwindSafe, catch_unwind};

fn no_panic<T>(label: &str, f: impl FnOnce() -> T) -> T {
    catch_unwind(AssertUnwindSafe(f))
        .unwrap_or_else(|_| panic!("{label} panicked on malformed or generated input"))
}

#[test]
fn cigar_and_edit_script_properties_are_coordinate_safe() {
    let valid = ["3=2X1I4D5S", "10M", "1=1X1I1D1S", "4294967295="];
    for cigar in valid {
        let runs = no_panic(cigar, || parse_cigar(cigar)).unwrap();
        let canonical = cigar_from_edit_script(&runs);
        assert_eq!(parse_cigar(&canonical).unwrap(), runs);
    }

    for cigar in [
        "",
        "3",
        "0=",
        "=3",
        "1Q",
        "4294967296=",
        "18446744073709551616=",
        "1=2",
        "1==",
    ] {
        let result = no_panic(cigar, || parse_cigar(cigar));
        assert!(result.is_err(), "malformed CIGAR {cigar:?} was accepted");
    }

    let segments = [
        BaseInterval { start: 8, end: 10 },
        BaseInterval { start: 0, end: 7 },
    ];
    let projection = project_cigar(&segments, "2=1I2D5=").unwrap();
    assert_eq!(projection.query_consumed, 9);
    assert_eq!(projection.target_consumed, 8);
    for span in &projection.spans {
        for interval in &span.query_segments {
            assert!(interval.end <= 10);
            assert!(interval.start <= interval.end);
        }
        assert!(span.target_interval.start <= span.target_interval.end);
    }

    let operations = [
        EditOperation::Equal,
        EditOperation::Substitution,
        EditOperation::Insertion,
        EditOperation::Deletion,
        EditOperation::SoftClip,
    ];
    for seed in 0..128u64 {
        let mut state = seed.wrapping_add(1);
        let mut runs = Vec::new();
        let mut query_length = 0u64;
        for _ in 0..8 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let operation = operations[(state as usize) % operations.len()];
            let length = (state % 5) + 1;
            if matches!(
                operation,
                EditOperation::Equal
                    | EditOperation::Substitution
                    | EditOperation::Deletion
                    | EditOperation::SoftClip
            ) {
                query_length += length;
            }
            runs.push(EditRun {
                operation,
                length: length as u32,
            });
        }
        let projection = project_edit_script(
            &[BaseInterval {
                start: 0,
                end: query_length,
            }],
            &runs,
        )
        .unwrap();
        assert_eq!(projection.query_consumed, query_length);
        let canonical_runs = parse_cigar(&cigar_from_edit_script(&runs)).unwrap();
        assert_eq!(parse_cigar(&projection.cigar).unwrap(), canonical_runs);
        for span in projection.spans {
            for interval in span.query_segments {
                assert!(interval.start <= interval.end);
                assert!(interval.end <= query_length);
            }
            assert!(span.target_interval.start <= span.target_interval.end);
        }
    }
}

#[test]
fn circular_union_and_gap_complement_preserve_bounds_and_total_length() {
    let mut state = 0x517c_c1b7_2722_0a95u64;
    for case in 0..256 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let length = 1 + state % 256;
        let mut segments = Vec::new();
        for _ in 0..(state as usize % 12) {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let start = state % (length + 1);
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let end = start + (state % (length - start + 1));
            segments.push(BaseInterval { start, end });
        }
        let covered = no_panic(&format!("circular union {case}"), || {
            circular_union(&segments, length)
        })
        .unwrap();
        assert_eq!(covered, union(segments.clone()));
        assert!(covered.windows(2).all(|pair| pair[0].end < pair[1].start));
        assert!(covered.iter().all(|interval| interval.end <= length));
        assert!(covered_length(&covered) <= length);

        let gaps = circular_gap_complement(&covered, length).unwrap();
        assert!(gaps.iter().all(|interval| interval.end <= length));
        assert_eq!(covered_length(&covered) + covered_length(&gaps), length);
        let mut all = covered.clone();
        all.extend(gaps);
        assert_eq!(
            union(all),
            vec![BaseInterval {
                start: 0,
                end: length
            }]
        );

        let wrap_start = length.saturating_sub(2);
        let wrap_end = if length > 2 { 1 } else { 0 };
        let split = split_circular(wrap_start, wrap_end, length).unwrap();
        assert!(split.iter().all(|interval| interval.end <= length));
    }

    assert!(circular_union(&[BaseInterval { start: 0, end: 2 }], 0).is_err());
    assert!(circular_union(&[BaseInterval { start: 2, end: 1 }], 3).is_err());
    assert!(circular_union(&[BaseInterval { start: 0, end: 4 }], 3).is_err());
    assert!(circular_gap_complement(&[BaseInterval { start: 0, end: 1 }], 0).is_err());
}

#[test]
fn json_serialization_is_stable_and_malformed_values_do_not_panic() {
    let record = TraceRecord::RunFooter(TraceRunFooter {
        schema_version: SCHEMA_VERSION.to_string(),
        run_id: "property-run".to_string(),
        completed_at_utc: "2026-01-01T00:00:00Z".to_string(),
        metagenomes_total: 3,
        metagenomes_with_candidates: 2,
        metagenomes_aligned: 1,
        metagenomes_failed: 1,
        alignments_total: 4,
        resource_metrics: ResourceMetrics::default(),
    });
    let first = serde_json::to_vec(&record).unwrap();
    for _ in 0..32 {
        assert_eq!(serde_json::to_vec(&record).unwrap(), first);
    }
    assert_eq!(
        serde_json::from_slice::<TraceRecord>(&first).unwrap(),
        record
    );

    let malformed = [
        b"".as_slice(),
        b"{".as_slice(),
        b"null".as_slice(),
        br#"{"record_type":"unknown"}"#.as_slice(),
        br#"{"record_type":"run_footer","schema_version":null}"#.as_slice(),
        &[0xff, 0xfe, 0xfd],
    ];
    for (index, bytes) in malformed.into_iter().enumerate() {
        let result = no_panic(&format!("JSON mutation {index}"), || {
            serde_json::from_slice::<TraceRecord>(bytes)
        });
        assert!(
            result.is_err(),
            "malformed JSON mutation {index} was accepted"
        );
    }
}
