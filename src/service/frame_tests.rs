use super::*;
use crate::search::SearchConfig;
use crate::service::records::{
    CompletionRecord, ConfirmedRecord, DropName, EntryState, FragmentRecord, LocatedRecord,
    StageName,
};

fn event(state: EntryState) -> Frame {
    Frame::Event {
        shard: 2,
        record: EventRecord {
            entry: 7,
            name: "entry7".into(),
            state,
        },
    }
}

fn summary() -> SummaryRecord {
    SummaryRecord {
        query_length: 5_000,
        completion: CompletionRecord::Stopped { deadline: true },
        possible_entries: 3,
        located_entries: 2,
        aligned_entries: 2,
        confirmed_entries: 1,
        seed_hits: 400,
        seed_locations: 380,
        candidate_bins: 12,
        bases_read: 40_000,
        alignments: 3,
        sketch_seconds: 0.25,
        locate_seconds: 0.5,
        align_seconds: 1.5,
    }
}

/// One frame of every kind, and an event of every state.
fn every_frame() -> Vec<Frame> {
    let confirmed = ConfirmedRecord {
        query_length: 5_000,
        covered_bases: 4_990,
        matches: 4_980,
        substitutions: 10,
        insertions: 2,
        deletions: 1,
        identity: 4_980.0 / 4_993.0,
        dropped_fragments: 1,
        fragments: vec![FragmentRecord {
            primary: true,
            contig: 1,
            query_segments: vec![(4_000, 5_000), (0, 3_990)],
            contig_start: 100,
            contig_end: 5_091,
            reverse: true,
            score: 9_000,
            cigar: "2000=1X1I2987=1D".into(),
        }],
    };
    vec![
        Frame::Hello {
            protocol_version: PROTOCOL_VERSION,
            shards: vec![ShardIdentity {
                name: "shard0".into(),
                entries: 12,
                parameters: Parameters::default(),
                manifest_sha256: "ab".repeat(32),
            }],
        },
        Frame::Search {
            protocol_version: PROTOCOL_VERSION,
            request: SearchRequest {
                query_name: "plasmid".into(),
                query: "ACGT".repeat(10),
                circular: true,
                settings: SearchSettings::of(&SearchConfig::default()),
                deadline_ms: Some(1_500),
            },
        },
        event(EntryState::Possible {
            coverage_bound: 4_000,
            containment_estimate: 0.75,
        }),
        event(EntryState::Located(LocatedRecord {
            coverage_bound: 3_000,
            containment_estimate: 0.5,
            windows: vec![(0, vec![(0, 2_048), (8_192, 9_216)]), (3, vec![(5, 6)])],
        })),
        event(EntryState::Confirmed(confirmed)),
        event(EntryState::Dropped {
            stage: StageName::Locate,
            reason: DropName::CoverageBound,
        }),
        event(EntryState::Unfinished {
            stage: StageName::Align,
        }),
        Frame::Summary {
            shard: 1,
            summary: summary(),
        },
        Frame::ShardFailed {
            shard: 3,
            message: "peer lost".into(),
        },
        Frame::Refused {
            message: "busy".into(),
        },
        Frame::Complete { complete: false },
    ]
}

#[test]
fn every_frame_round_trips() {
    let frames = every_frame();
    let mut bytes = Vec::new();
    for frame in &frames {
        write_frame(&mut bytes, frame).unwrap();
    }
    let mut reader = bytes.as_slice();
    for frame in &frames {
        assert_eq!(read_frame(&mut reader).unwrap().as_ref(), Some(frame));
    }
    assert!(read_frame(&mut reader).unwrap().is_none());
}

#[test]
fn oversized_and_cut_frames_fail() {
    let mut bytes = (MAX_FRAME_BYTES + 1).to_le_bytes().to_vec();
    bytes.extend_from_slice(b"{}");
    assert!(matches!(
        read_frame(&mut bytes.as_slice()),
        Err(ServiceError::Frame(_))
    ));

    let mut bytes = Vec::new();
    write_frame(&mut bytes, &Frame::Complete { complete: true }).unwrap();
    bytes.pop();
    assert!(read_frame(&mut bytes.as_slice()).is_err());
    assert!(read_frame(&mut &bytes[..2]).is_err());

    let mut bytes = 2u32.to_le_bytes().to_vec();
    bytes.extend_from_slice(b"{}");
    assert!(matches!(
        read_frame(&mut bytes.as_slice()),
        Err(ServiceError::Frame(_))
    ));
}

#[test]
fn another_protocol_version_is_named() {
    assert!(check_version(PROTOCOL_VERSION).is_ok());
    let error = check_version(PROTOCOL_VERSION + 1).unwrap_err();
    assert!(matches!(error, ServiceError::Version { .. }));
    assert!(error.to_string().contains("protocol version 2"));
}

#[test]
fn shards_with_other_bucket_bits_share_seeds_but_not_with_another_k() {
    let identity = |parameters| ShardIdentity {
        name: "shard".into(),
        entries: 1,
        parameters,
        manifest_sha256: String::new(),
    };
    let ours = identity(Parameters::default());
    let more_buckets = identity(Parameters {
        bucket_bits: 20,
        ..Parameters::default()
    });
    let other_k = identity(Parameters {
        k: 21,
        ..Parameters::default()
    });
    assert!(ours.same_seeds(&more_buckets));
    assert!(!ours.same_seeds(&other_k));
}
