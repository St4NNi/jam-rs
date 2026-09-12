use crate::jidx::sha256;
use crate::jidx_reader::JidxReader;
use crate::jidx_writer::{ContigInput, JidxInput, JidxWriter, MetagenomeInput};
use crate::owner_format::checksum_layout;
use crate::shared_format::{HEADER_BYTES, PAGE_BYTES, Section, SharedError, SharedHeader};
use crate::shared_reader::{CoreKind, CoreRow, SharedGroup, SharedReadStats, SharedReader};
use crate::shared_seed::{SharedKey, SharedSeed};
use crate::shared_writer::{IndexedSeed, SharedBuildStats, write_shared_index};
const TARGET_CORE: u32 = 1_000_000;
const TARGET_CONTEXT: u32 = 0xabc0_0123;

fn fixture(preceding: u32) -> (tempfile::TempDir, SharedReader, SharedBuildStats) {
    let directory = tempfile::tempdir().unwrap();
    let jidx = directory.path().join("metadata.jidx");
    let mut writer = JidxWriter::new(
        &jidx,
        &JidxInput {
            k: 15,
            rescue_k15: false,
            minimizer_window: 64,
            jam_sha256: [1; 32],
            manifest_sha256: [2; 32],
        },
    )
    .unwrap();
    for (id, name) in ["a", "b", "c"].into_iter().enumerate() {
        writer
            .begin_metagenome(MetagenomeInput {
                name: name.into(),
                bgzf_uri: format!("{name}.bgz"),
                bgzf_bytes: 100,
                bgzf_sha256: [id as u8 + 3; 32],
                gzi: vec![0; 8],
            })
            .unwrap();
        writer
            .begin_contig(ContigInput {
                name: format!("{name}-contig"),
                length: 100_000,
                fasta_offset: 4,
                line_bases: 80,
                line_width: 81,
            })
            .unwrap();
    }
    writer.finish().unwrap();
    let reference = JidxReader::open(&jidx).unwrap();
    let mut seeds = (0..preceding)
        .flat_map(|core| {
            [0, 1].map(|offset| IndexedSeed {
                member: 0,
                contig: 0,
                seed: SharedSeed {
                    core,
                    context: 0,
                    flags: 0,
                    position: u64::from(core % 90_000 + offset),
                },
            })
        })
        .collect::<Vec<_>>();
    seeds.extend([
        indexed_target(1, 1, 100, false),
        indexed_target(1, 1, 200, true),
        indexed_target(2, 2, 300, false),
        IndexedSeed {
            member: 0,
            contig: 0,
            seed: SharedSeed {
                core: TARGET_CORE + 1,
                context: 0,
                flags: 0,
                position: 500,
            },
        },
    ]);
    let shared = directory.path().join("fixture.shared");
    let stats = write_shared_index(&reference, &shared, 64, &mut seeds).unwrap();
    drop(reference);
    (
        directory,
        SharedReader::open_observed(shared).unwrap(),
        stats,
    )
}

fn indexed_target(member: u32, contig: u32, position: u64, reverse: bool) -> IndexedSeed {
    IndexedSeed {
        member,
        contig,
        seed: SharedSeed {
            core: TARGET_CORE,
            context: TARGET_CONTEXT,
            flags: u8::from(reverse) | 6,
            position,
        },
    }
}

fn delta(after: SharedReadStats, before: SharedReadStats) -> [u64; 5] {
    [
        after.core_descriptor_inspections - before.core_descriptor_inspections,
        after.group_descriptor_inspections - before.group_descriptor_inspections,
        after.member_descriptor_inspections - before.member_descriptor_inspections,
        after.references_decoded - before.references_decoded,
        after.physical_positions_decoded - before.physical_positions_decoded,
    ]
}

#[test]
fn absent_heavy_batch_preserves_successful_associations_and_reports_capacity() {
    use crate::trace_index::{TraceDocument, TraceIndex};
    use std::mem::size_of;

    let (_directory, reader, _) = fixture(16);
    let index = TraceIndex::Shared(Box::new(reader));
    let present = [
        SharedKey::core(TARGET_CORE),
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT >> 20,
            length: 21,
        },
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT,
            length: 31,
        },
    ];
    let mut requests = (2_000_000..2_032_768)
        .map(|key| (key, 0))
        .collect::<Vec<_>>();
    requests.extend(
        present
            .iter()
            .flat_map(|key| [0, 1].map(|query| (key.packed().unwrap(), query))),
    );
    let requests_capacity = requests.capacity();
    let lookup = crate::trace_batch::prepare_lookup(&index, requests, 2, true)
        .unwrap()
        .unwrap();
    let mut expected = present.map(|key| key.packed().unwrap());
    expected.sort_unstable();
    for range in &lookup.query_ranges {
        assert_eq!(
            lookup.query_entries[range.clone()]
                .iter()
                .map(|&ordinal| lookup.entries[ordinal as usize].0)
                .collect::<Vec<_>>(),
            expected
        );
    }
    assert_eq!(
        lookup
            .entries
            .iter()
            .filter(|entry| entry.1.is_some())
            .count(),
        3
    );
    assert_eq!(lookup.postings.len(), 3);
    assert!(lookup.postings_complete);
    assert_eq!(lookup.capacity_bytes, lookup._reservation.bytes);
    assert!(lookup.peak_capacity_bound > lookup.capacity_bytes);
    assert!(lookup.capacity_bytes < 2_000_000);
    assert_eq!(lookup.entries.len(), 3);
    assert_eq!(lookup.entries.capacity(), 3);
    assert_eq!(lookup.query_entries.capacity(), 6);
    for posting in lookup.postings.iter().flatten() {
        assert_eq!(posting.documents.len(), 2);
        assert_eq!(
            posting
                .documents
                .iter()
                .map(|doc| doc.occurrence_count())
                .sum::<u64>(),
            3
        );
        assert_eq!(
            posting
                .occurrences
                .as_ref()
                .unwrap()
                .iter()
                .map(Vec::len)
                .sum::<usize>(),
            3
        );
    }
    let first_key = lookup.entries[0].0;
    let direct = lookup.posting(first_key, Some(0)).unwrap().unwrap();
    let fallback = lookup.posting(first_key, None).unwrap().unwrap();
    assert!(std::ptr::eq(direct, fallback));
    assert!(lookup.posting(u64::MAX, None).unwrap().is_none());
    assert!(matches!(
        lookup.posting(first_key, Some(lookup.entries.len())),
        Err(crate::trace::TraceError::Invalid("batch posting ordinal"))
    ));
    assert!(matches!(
        lookup.posting(lookup.entries[1].0, Some(0)),
        Err(crate::trace::TraceError::Invalid("batch posting ordinal"))
    ));
    println!(
        "requests_capacity={requests_capacity} entries_length={} entries_capacity={} associations_capacity={} retained_bytes={} reserved_bytes={} raw_document_bytes={} trace_document_bytes={} shared_group_bytes={} optional_group_bytes={} shared_member_bytes={}",
        lookup.entries.len(),
        lookup.entries.capacity(),
        lookup.query_entries.capacity(),
        lookup.capacity_bytes,
        lookup._reservation.bytes,
        size_of::<crate::jidx_reader::SeedDocument>(),
        size_of::<TraceDocument>(),
        size_of::<SharedGroup>(),
        size_of::<Option<SharedGroup>>(),
        size_of::<crate::shared_reader::SharedMember>(),
    );
}

#[test]
fn full_lookup_work_is_stable_across_worker_counts() {
    use crate::trace_index::TraceIndex;

    let (directory, reader, _) = fixture(0);
    drop(reader);
    let path = directory.path().join("fixture.shared");
    let present = [
        SharedKey::core(TARGET_CORE),
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT >> 20,
            length: 21,
        },
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT,
            length: 31,
        },
    ];
    let mut requests = (0..=65_536)
        .map(|context| {
            (
                SharedKey {
                    core: TARGET_CORE,
                    context,
                    length: 31,
                }
                .packed()
                .unwrap(),
                0,
            )
        })
        .collect::<Vec<_>>();
    requests.extend(
        present
            .iter()
            .flat_map(|key| [0, 1].map(|query| (key.packed().unwrap(), query))),
    );
    let mut expected = None;
    for threads in [1, 4, 8, 16] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let actual = pool.install(|| {
            let index = TraceIndex::Shared(Box::new(SharedReader::open_observed(&path).unwrap()));
            let lookup = crate::trace_batch::prepare_lookup(&index, requests.clone(), 2, true)
                .unwrap()
                .unwrap();
            let entries = lookup
                .entries
                .iter()
                .map(|(key, seed)| (*key, seed.unwrap().document_frequency()))
                .collect::<Vec<_>>();
            let associations = lookup
                .query_ranges
                .iter()
                .map(|range| {
                    lookup.query_entries[range.clone()]
                        .iter()
                        .map(|&ordinal| lookup.entries[ordinal as usize].0)
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let postings = lookup
                .postings
                .iter()
                .enumerate()
                .filter_map(|(ordinal, posting)| posting.as_ref().map(|posting| (ordinal, posting)))
                .map(|(ordinal, posting)| {
                    let key = lookup.entries[ordinal].0;
                    let documents = posting
                        .documents
                        .iter()
                        .map(|document| (document.metagenome_id(), document.occurrence_count()))
                        .collect::<Vec<_>>();
                    let occurrences = posting
                        .occurrences
                        .as_ref()
                        .unwrap()
                        .iter()
                        .map(|values| {
                            values
                                .iter()
                                .map(|value| {
                                    (value.contig_id, value.position, value.canonical_orientation)
                                })
                                .collect::<Vec<_>>()
                        })
                        .collect::<Vec<_>>();
                    (key, documents, occurrences)
                })
                .collect::<Vec<_>>();
            let reads = match &index {
                TraceIndex::Shared(reader) => reader.stats(),
                _ => unreachable!(),
            };
            (
                lookup.lookup_tasks,
                lookup.lookup_plan_hash,
                lookup.split_core_resolutions,
                lookup.attempted_keys,
                entries,
                associations,
                postings,
                (
                    reads.core_key_inspections,
                    reads.core_descriptor_inspections,
                    reads.group_descriptor_inspections,
                    reads.member_descriptor_inspections,
                    reads.core_resolutions_present,
                    reads.core_resolutions_absent,
                    reads.grouped_core_rows,
                    reads.references_decoded,
                    reads.physical_positions_decoded,
                ),
            )
        });
        assert_eq!(actual.0, 3);
        assert_ne!(actual.1, 0);
        assert_eq!(actual.2, 2);
        assert_eq!(actual.3, 65_540);
        assert_eq!(actual.4.len(), 3);
        assert_eq!(actual.5.iter().map(Vec::len).collect::<Vec<_>>(), [3, 3]);
        if let Some(expected) = &expected {
            assert_eq!(&actual, expected, "thread count {threads}");
        } else {
            expected = Some(actual);
        }
    }
}

#[test]
fn exact_counts_and_absence_do_not_decode_payloads_as_prefix_grows() {
    for preceding in [16u32, 4096, 16_384] {
        let (_directory, reader, build) = fixture(preceding);
        assert!(reader.stats().observed);
        assert_eq!(build.repeated_core_positions, u64::from(preceding) * 2 + 3);
        assert_eq!(build.occurrence_references, u64::from(preceding) * 2 + 9);
        let before = reader.stats();
        let group = reader.find(SharedKey::core(TARGET_CORE)).unwrap().unwrap();
        assert_eq!((group.member_count(), group.occurrence_count()), (2, 3));
        let after = reader.stats();
        let work = delta(after, before);
        let logarithmic_bound = u64::from(u32::BITS - (preceding + 1).leading_zeros()) + 1;
        assert_eq!(work[0], 1);
        assert!(after.core_key_inspections - before.core_key_inspections <= logarithmic_bound);
        assert!(work[1] <= 3);
        assert_eq!(&work[2..], &[0, 0, 0]);

        let before = reader.stats();
        assert!(
            reader
                .find(SharedKey::core(CORE_LIMIT - 1))
                .unwrap()
                .is_none()
        );
        let after = reader.stats();
        let work = delta(after, before);
        assert_eq!(work[0], 0);
        assert!(after.core_key_inspections > before.core_key_inspections);
        assert!(after.core_key_inspections - before.core_key_inspections <= logarithmic_bound);
        assert_eq!(&work[1..], &[0, 0, 0, 0]);

        let before = reader.stats();
        let context = SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT >> 20,
            length: 21,
        };
        let group = reader.find(context).unwrap().unwrap();
        assert_eq!((group.member_count(), group.occurrence_count()), (2, 3));
        let after = reader.stats();
        let work = delta(after, before);
        assert_eq!(work[0], 1);
        assert!(after.core_key_inspections - before.core_key_inspections <= logarithmic_bound);
        assert!(work[1] <= 3);
        assert_eq!(&work[2..], &[0, 0, 0]);

        let cached = reader
            .group_at(group.core_ordinal(), context)
            .unwrap()
            .unwrap();
        assert_eq!(cached, group);
        let before = reader.stats();
        let selected = reader.member(group, 2).unwrap().unwrap();
        let work = delta(reader.stats(), before);
        assert_eq!(selected.occurrence_count(), 1);
        assert!(work[2] <= 2);
        assert_eq!(&work[3..], &[0, 0]);

        let before = reader.stats();
        let occurrences = reader.member_occurrences(group, selected).unwrap();
        assert_eq!(
            occurrences
                .iter()
                .map(|occurrence| (occurrence.contig_id, occurrence.position))
                .collect::<Vec<_>>(),
            [(2, 300)]
        );
        let work = delta(reader.stats(), before);
        assert_eq!(work, [0, 0, 0, 1, 1]);

        let before = reader.stats();
        assert!(reader.member(group, 0).unwrap().is_none());
        let work = delta(reader.stats(), before);
        assert!(work[2] <= 2);
        assert_eq!(&work[3..], &[0, 0]);

        let members = reader.members(group).unwrap();
        assert_eq!(
            members
                .iter()
                .map(|member| (member.metagenome_id, member.occurrence_count()))
                .collect::<Vec<_>>(),
            [(1, 2), (2, 1)]
        );
        assert_ne!(
            members[0].occurrence_storage_identity(),
            members[1].occurrence_storage_identity()
        );

        let singleton = reader
            .find(SharedKey::core(TARGET_CORE + 1))
            .unwrap()
            .unwrap();
        assert_eq!(
            (singleton.member_count(), singleton.occurrence_count()),
            (1, 1)
        );
        let singleton_member = reader.member(singleton, 0).unwrap().unwrap();
        let before = reader.stats();
        let occurrence = reader
            .member_occurrences(singleton, singleton_member)
            .unwrap();
        assert_eq!((occurrence[0].contig_id, occurrence[0].position), (0, 500));
        assert_eq!(delta(reader.stats(), before)[3..], [0, 1]);
        assert!(
            reader
                .occurrence_block(singleton, singleton_member, 0, 0)
                .is_err()
        );
        assert_eq!(reader.metagenome(0).unwrap().unwrap().gzi, [0; 8]);
        let contig = reader.contig(0).unwrap().unwrap();
        assert_eq!((contig.metagenome_id, contig.length), (0, 100_000));

        for absent_context in [
            SharedKey {
                core: TARGET_CORE,
                context: (TARGET_CONTEXT >> 20) ^ 1,
                length: 21,
            },
            SharedKey {
                core: TARGET_CORE,
                context: TARGET_CONTEXT ^ 1,
                length: 31,
            },
        ] {
            let before = reader.stats();
            assert!(reader.find(absent_context).unwrap().is_none());
            let work = delta(reader.stats(), before);
            assert!(work[0] <= logarithmic_bound);
            assert!(work[1] <= 3);
            assert_eq!(&work[2..], &[0, 0, 0]);
        }
    }
}

#[test]
fn key_only_probe_handles_empty_singleton_and_dictionary_edges() {
    let (directory, reader, _) = fixture(16);
    let metadata = directory.path().join("metadata.jidx");
    let reference = JidxReader::open(&metadata).unwrap();

    let empty = directory.path().join("empty.shared");
    write_shared_index(&reference, &empty, 64, &mut []).unwrap();
    let empty = SharedReader::open_observed(&empty).unwrap();
    assert_eq!(empty.core_count(), 0);
    assert!(empty.find(SharedKey::core(0)).unwrap().is_none());
    let stats = empty.stats();
    assert_eq!(stats.core_key_inspections, 0);
    assert_eq!(stats.core_descriptor_inspections, 0);

    let empty_packed = directory.path().join("empty-packed.shared");
    let empty_compact = directory.path().join("empty-compact.shared");
    crate::shared_pack::repack_shared_index(directory.path().join("empty.shared"), &empty_packed)
        .unwrap();
    crate::shared_pack::repack_shared_cores(&empty_packed, &empty_compact).unwrap();
    let empty = SharedReader::open_observed(&empty_compact).unwrap();
    assert!(empty.has_core_prefixes());
    assert!(empty.find(SharedKey::core(0)).unwrap().is_none());
    assert!(
        empty
            .find(SharedKey::core(CORE_LIMIT - 1))
            .unwrap()
            .is_none()
    );
    let stats = empty.stats();
    assert_eq!(stats.core_key_inspections, 0);
    assert_eq!(stats.core_descriptor_inspections, 0);
    let mut sparse = Vec::new();
    empty
        .resolve_sorted_cores_into(&[0, CORE_LIMIT - 1], &mut sparse)
        .unwrap();
    assert!(sparse.is_empty());
    assert_eq!(sparse.capacity(), 0);

    let singleton = directory.path().join("singleton.shared");
    write_shared_index(
        &reference,
        &singleton,
        64,
        &mut [IndexedSeed {
            member: 0,
            contig: 0,
            seed: SharedSeed {
                core: 7,
                context: 0,
                flags: 0,
                position: 500,
            },
        }],
    )
    .unwrap();
    let singleton = SharedReader::open_observed(&singleton).unwrap();
    assert!(singleton.find(SharedKey::core(6)).unwrap().is_none());
    assert!(singleton.find(SharedKey::core(7)).unwrap().is_some());
    assert!(singleton.find(SharedKey::core(8)).unwrap().is_none());
    let stats = singleton.stats();
    assert_eq!(stats.core_key_inspections, 3);
    assert_eq!(stats.core_descriptor_inspections, 1);

    let singleton_packed = directory.path().join("singleton-packed.shared");
    let singleton_compact = directory.path().join("singleton-compact.shared");
    crate::shared_pack::repack_shared_index(
        directory.path().join("singleton.shared"),
        &singleton_packed,
    )
    .unwrap();
    crate::shared_pack::repack_shared_cores(&singleton_packed, &singleton_compact).unwrap();
    let singleton = SharedReader::open_observed(&singleton_compact).unwrap();
    let mut sparse = Vec::new();
    singleton
        .resolve_sorted_cores_into(&[6, 7, 8], &mut sparse)
        .unwrap();
    assert_eq!(
        sparse
            .iter()
            .map(|group| group.key().core)
            .collect::<Vec<_>>(),
        [7]
    );
    assert_eq!(sparse.capacity(), 1);

    let before = reader.stats();
    let groups = reader
        .find_many(&[SharedKey::core(0), SharedKey::core(TARGET_CORE + 1)])
        .unwrap();
    assert!(groups.iter().all(Option::is_some));
    let after = reader.stats();
    assert_eq!(
        after.core_resolutions_present - before.core_resolutions_present,
        2
    );
    assert_eq!(
        after.core_descriptor_inspections - before.core_descriptor_inspections,
        2
    );
    assert!(after.core_key_inspections > before.core_key_inspections);
}

#[test]
fn sorted_core_sparse_results_match_dense_lookup_and_bound_output() {
    let (directory, reader, _) = fixture(8192);
    drop(reader);
    let source = directory.path().join("fixture.shared");
    let packed = directory.path().join("packed.shared");
    let compact = directory.path().join("compact.shared");
    crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
    crate::shared_pack::repack_shared_cores(&packed, &compact).unwrap();

    let workload = |present: usize| {
        let mut cores = (0..present as u32).collect::<Vec<_>>();
        cores.extend(8192..8192 + (4096 - present) as u32);
        cores
    };
    let mut spread = (1..=64u32)
        .map(|prefix| (prefix << 14) | 12_345)
        .filter(|&core| core != TARGET_CORE && core != TARGET_CORE + 1)
        .collect::<Vec<_>>();
    spread.extend([0, TARGET_CORE, TARGET_CORE + 1]);
    spread.sort_unstable();
    let workloads = [
        ("absent", workload(0), 0),
        ("10-percent", workload(410), 410),
        ("50-percent", workload(2048), 2048),
        ("all-present", workload(4096), 4096),
        ("spread-prefixes", spread, 3),
    ];
    let evidence =
        |group: &SharedGroup| (group.key(), group.member_count(), group.occurrence_count());
    for (name, cores, expected_present) in workloads {
        let reference = SharedReader::open(&compact).unwrap();
        let keys = cores
            .iter()
            .copied()
            .map(SharedKey::core)
            .collect::<Vec<_>>();
        let expected = reference
            .find_many(&keys)
            .unwrap()
            .into_iter()
            .flatten()
            .map(|group| evidence(&group))
            .collect::<Vec<_>>();
        assert_eq!(expected.len(), expected_present, "{name}");

        let reader = SharedReader::open_observed(&compact).unwrap();
        let before = reader.stats();
        let mut actual = Vec::new();
        reader
            .resolve_sorted_cores_into(&cores, &mut actual)
            .unwrap();
        let after = reader.stats();
        assert_eq!(
            actual.iter().map(evidence).collect::<Vec<_>>(),
            expected,
            "{name}"
        );
        assert!(actual.capacity() <= cores.len(), "{name}");
        assert!(actual.capacity() >= actual.len(), "{name}");
        assert_eq!(
            after.core_descriptor_inspections - before.core_descriptor_inspections,
            actual.len() as u64,
            "{name}"
        );
        assert_eq!(
            after.member_descriptor_inspections, before.member_descriptor_inspections,
            "{name}"
        );
        assert_eq!(
            after.references_decoded, before.references_decoded,
            "{name}"
        );
        assert_eq!(
            after.physical_positions_decoded, before.physical_positions_decoded,
            "{name}"
        );
        if expected_present == 0 {
            assert_eq!(actual.capacity(), 0);
            assert!(after.core_view_creations > before.core_view_creations);
            assert!(after.core_view_comparisons > before.core_view_comparisons);
        }
    }

    let reader = SharedReader::open(&compact).unwrap();
    for invalid in [vec![1, 0], vec![1, 1], vec![CORE_LIMIT]] {
        let mut output = Vec::new();
        assert!(matches!(
            reader.resolve_sorted_cores_into(&invalid, &mut output),
            Err(SharedError::Invalid("sorted cores"))
        ));
        assert!(output.is_empty());
    }
    let mut oversized = Vec::with_capacity(2);
    assert!(matches!(
        reader.resolve_sorted_cores_into(&[0], &mut oversized),
        Err(SharedError::ResourceLimit)
    ));

    let first = SharedReader::open(&compact).unwrap();
    let second = SharedReader::open(&compact).unwrap();
    let mut left = Vec::new();
    let mut right = Vec::new();
    first
        .resolve_sorted_cores_into(&[0, TARGET_CORE], &mut left)
        .unwrap();
    second
        .resolve_sorted_cores_into(&[0, TARGET_CORE], &mut right)
        .unwrap();
    assert_eq!(
        left.iter().map(evidence).collect::<Vec<_>>(),
        right.iter().map(evidence).collect::<Vec<_>>()
    );
    assert_ne!(first.reader_token(), second.reader_token());
    assert!(matches!(
        second.members(left[0]),
        Err(SharedError::Invalid("group handle"))
    ));
}

#[test]
fn compact_core_rows_decode_narrow_and_wide_values_exactly() {
    let hot = 7u32.to_le_bytes();
    let mut narrow = [0u8; 13];
    narrow[..4].copy_from_slice(&0xabc0_0000u32.to_le_bytes());
    narrow[4..8].copy_from_slice(&3u32.to_le_bytes());
    narrow[8..12].copy_from_slice(&u32::MAX.to_le_bytes());
    narrow[12] = 3;
    let row = CoreRow::decode_compact(&hot, &narrow, 4, 0).unwrap();
    assert_eq!(row.core, 7);
    assert!(matches!(
        row.kind,
        CoreKind::Singleton {
            context: 0xabc0_0000,
            contig_id: 3,
            flags: 3,
            position,
        } if position == u64::from(u32::MAX)
    ));

    let mut wide = [0u8; 21];
    wide[..4].copy_from_slice(&0xabc0_0123u32.to_le_bytes());
    wide[4..8].copy_from_slice(&3u32.to_le_bytes());
    wide[8..16].copy_from_slice(&(u64::from(u32::MAX) + 1).to_le_bytes());
    wide[20] = 7;
    let row = CoreRow::decode_compact(&hot, &wide, 4, 0).unwrap();
    assert!(matches!(
        row.kind,
        CoreKind::Singleton {
            context: 0xabc0_0123,
            contig_id: 3,
            flags: 7,
            position,
        } if position == u64::from(u32::MAX) + 1
    ));
    wide[16] = 1;
    assert!(matches!(
        CoreRow::decode_compact(&hot, &wide, 4, 0),
        Err(SharedError::Invalid("singleton core"))
    ));

    let hot = (9u32 | (1 << 31)).to_le_bytes();
    let mut repeated = [0u8; 21];
    let first_group = u64::from(u32::MAX) + 1;
    let occurrence_count = u64::from(u32::MAX) + 2;
    repeated[..8].copy_from_slice(&first_group.to_le_bytes());
    repeated[8..12].copy_from_slice(&2u32.to_le_bytes());
    repeated[12..20].copy_from_slice(&occurrence_count.to_le_bytes());
    let row = CoreRow::decode_compact(&hot, &repeated, 4, first_group + 2).unwrap();
    assert!(matches!(
        row.kind,
        CoreKind::Repeated {
            first_group: first,
            group_count: 2,
            occurrence_count: count,
        } if first == first_group && count == occurrence_count
    ));
    repeated[20] = 1;
    assert!(matches!(
        CoreRow::decode_compact(&hot, &repeated, 4, first_group + 2),
        Err(SharedError::Invalid("repeated core"))
    ));
}

#[test]
fn compact_handles_are_reader_bound_and_keep_resolved_locations() {
    use std::mem::size_of;

    assert_eq!(size_of::<SharedGroup>(), 64);
    assert!(size_of::<Option<SharedGroup>>() <= 64);
    assert_eq!(size_of::<crate::shared_reader::SharedMember>(), 56);
    assert!(size_of::<Option<crate::shared_reader::SharedMember>>() <= 56);
    assert_eq!(
        size_of::<crate::shared_reader::SharedOccurrenceStorage>(),
        32
    );
    println!(
        "shared_group_bytes={} optional_group_bytes={} shared_member_bytes={} optional_member_bytes={} occurrence_storage_bytes={}",
        size_of::<SharedGroup>(),
        size_of::<Option<SharedGroup>>(),
        size_of::<crate::shared_reader::SharedMember>(),
        size_of::<Option<crate::shared_reader::SharedMember>>(),
        size_of::<crate::shared_reader::SharedOccurrenceStorage>(),
    );

    let (directory, reader, _) = fixture(16);
    let path = directory.path().join("fixture.shared");
    let other = SharedReader::open_observed(&path).unwrap();
    assert_ne!(reader.reader_token(), other.reader_token());
    let group = reader.find(SharedKey::core(TARGET_CORE)).unwrap().unwrap();
    let member = reader.member(group, 2).unwrap().unwrap();
    assert!(matches!(
        other.members(group),
        Err(SharedError::Invalid("group handle"))
    ));
    assert!(matches!(
        other.member_occurrences(group, member),
        Err(SharedError::Invalid("group handle"))
    ));

    let singleton = reader
        .find(SharedKey::core(TARGET_CORE + 1))
        .unwrap()
        .unwrap();
    let singleton_member = reader.member(singleton, 0).unwrap().unwrap();
    assert!(matches!(
        reader.member_occurrences(group, singleton_member),
        Err(SharedError::Invalid("member handle"))
    ));

    let before = reader.stats();
    let occurrences = reader.member_occurrences(group, member).unwrap();
    assert_eq!(
        occurrences
            .iter()
            .map(|occurrence| (occurrence.contig_id, occurrence.position))
            .collect::<Vec<_>>(),
        [(2, 300)]
    );
    assert_eq!(delta(reader.stats(), before), [0, 0, 0, 1, 1]);

    let (inline_directory, source) = grouped_lookup_fixture(0);
    let packed = inline_directory.path().join("inline.shared");
    crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
    let inline_reader = SharedReader::open_observed(&packed).unwrap();
    let inline_group = inline_reader
        .find(SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT ^ 1,
            length: 31,
        })
        .unwrap()
        .unwrap();
    let inline_member = inline_reader.members(inline_group).unwrap().remove(0);
    let before = inline_reader.stats();
    let occurrences = inline_reader
        .member_occurrences(inline_group, inline_member)
        .unwrap();
    assert_eq!(
        occurrences
            .iter()
            .map(|occurrence| (occurrence.contig_id, occurrence.position))
            .collect::<Vec<_>>(),
        [(2, 300)]
    );
    assert_eq!(delta(inline_reader.stats(), before), [0, 0, 0, 0, 1]);
}

#[test]
fn numeric_contig_reads_only_checked_numeric_metadata() {
    let (_directory, reader, _) = fixture(0);
    let before = reader.stats();
    let numeric = reader.numeric_contig(2).unwrap().unwrap();
    let after = reader.stats();
    assert_eq!(
        (numeric.id, numeric.metagenome_id, numeric.length),
        (2, 2, 100_000)
    );
    assert_eq!(
        after.numeric_contig_resolutions - before.numeric_contig_resolutions,
        1
    );
    assert_eq!(
        after.file.requested_bytes - before.file.requested_bytes,
        u64::from(crate::jidx::CONTIG_RECORD_SIZE)
    );
    assert_eq!(after.file.requested_pages - before.file.requested_pages, 1);
    let full = reader.contig(2).unwrap().unwrap();
    assert_eq!(
        (numeric.id, numeric.metagenome_id, numeric.length),
        (full.id, full.metagenome_id, full.length)
    );

    let before = reader.stats();
    assert!(
        reader
            .numeric_contig(reader.contig_count())
            .unwrap()
            .is_none()
    );
    let after = reader.stats();
    assert_eq!(
        after.numeric_contig_resolutions,
        before.numeric_contig_resolutions
    );
    assert_eq!(after.file.requested_bytes, before.file.requested_bytes);
}

#[test]
fn numeric_contig_ignores_names_and_rejects_invalid_numeric_fields() {
    let (directory, reader, _) = fixture(0);
    drop(reader);
    let path = directory.path().join("fixture.shared");
    mutate_and_resign(&path, |bytes, header| {
        let contig = header.section(Section::Contigs).offset as usize;
        bytes[contig + 4..contig + 8].copy_from_slice(&u32::MAX.to_le_bytes());
    });
    let reader = SharedReader::open_observed(&path).unwrap();
    assert_eq!(reader.numeric_contig(0).unwrap().unwrap().length, 100_000);
    assert!(reader.contig(0).is_err());

    for invalid_document in [false, true] {
        let (directory, reader, _) = fixture(0);
        drop(reader);
        let path = directory.path().join("fixture.shared");
        mutate_and_resign(&path, |bytes, header| {
            let contig = header.section(Section::Contigs).offset as usize;
            if invalid_document {
                bytes[contig..contig + 4].copy_from_slice(&header.document_count.to_le_bytes());
            } else {
                bytes[contig + 16..contig + 24].copy_from_slice(&0u64.to_le_bytes());
            }
        });
        let reader = SharedReader::open_observed(&path).unwrap();
        assert!(matches!(
            reader.numeric_contig(0),
            Err(SharedError::Invalid("contig metadata"))
        ));
        assert_eq!(reader.stats().numeric_contig_resolutions, 0);
    }
}

#[test]
fn trace_adapter_reuses_resolved_group_and_member_handles() {
    use crate::trace_index::TraceIndex;

    let (_directory, reader, _) = fixture(4096);
    let index = TraceIndex::Shared(Box::new(reader));
    let seed = index
        .find_seed(SharedKey::core(TARGET_CORE).packed().unwrap())
        .unwrap()
        .unwrap();
    let documents = index.seed_documents(seed).unwrap();
    let document = documents
        .into_iter()
        .find(|document| document.metagenome_id() == 2)
        .unwrap();
    let before = match &index {
        TraceIndex::Shared(reader) => reader.stats(),
        _ => unreachable!(),
    };
    let occurrences = index.seed_document_occurrences(seed, document).unwrap();
    let after = match &index {
        TraceIndex::Shared(reader) => reader.stats(),
        _ => unreachable!(),
    };
    assert_eq!(
        occurrences
            .iter()
            .map(|occurrence| (occurrence.contig_id, occurrence.position))
            .collect::<Vec<_>>(),
        [(2, 300)]
    );
    assert_eq!(
        after.group_descriptor_inspections,
        before.group_descriptor_inspections
    );
    assert_eq!(
        after.member_descriptor_inspections,
        before.member_descriptor_inspections
    );
}

#[test]
fn retained_core_resolves_contexts_without_another_core_search() {
    let (_directory, path) = grouped_lookup_fixture(4096);
    let reader = SharedReader::open_observed(&path).unwrap();
    let core = SharedKey::core(TARGET_CORE);
    let core_group = reader.find(core).unwrap().unwrap();
    let requests = [
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT,
            length: 31,
        },
        core,
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT ^ 2,
            length: 31,
        },
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT >> 20,
            length: 21,
        },
    ];
    let before = reader.stats();
    let groups = reader.find_in_core(core_group, &requests).unwrap();
    let after = reader.stats();
    assert_eq!(
        after.core_resolutions_present,
        before.core_resolutions_present
    );
    assert_eq!(
        after.core_resolutions_absent,
        before.core_resolutions_absent
    );
    assert_eq!(after.grouped_core_rows, before.grouped_core_rows);
    assert_eq!(
        after.core_descriptor_inspections,
        before.core_descriptor_inspections + 1
    );

    let scalar = SharedReader::open(&path).unwrap();
    let expected = requests
        .iter()
        .map(|&key| group_evidence(&scalar, scalar.find(key).unwrap()))
        .collect::<Vec<_>>();
    let actual = groups
        .into_iter()
        .map(|group| group_evidence(&reader, group))
        .collect::<Vec<_>>();
    assert_eq!(actual, expected);
    assert!(
        reader
            .find_in_core(core_group, &[SharedKey::core(TARGET_CORE + 1)])
            .is_err()
    );
}

#[test]
fn core_first_absence_preserves_complete_linear_and_circular_accounting() {
    use crate::trace::{SearchCompletion, TraceConfig, TraceEngine};

    let (directory, reader, _) = fixture(0);
    drop(reader);
    let shared = directory.path().join("fixture.shared");
    let query = vec![b'A'; 2000];
    let config = TraceConfig {
        use_sketch: false,
        ..TraceConfig::default()
    };
    for (circular, associations, cores) in [(false, 5936, 1986), (true, 6000, 2000)] {
        let engine = TraceEngine::open_shared_observed(&shared, None, true).unwrap();
        let result = engine
            .search(
                if circular { "circular" } else { "linear" },
                &query,
                TraceConfig { circular, ..config },
            )
            .unwrap();
        assert_eq!(result.completion, SearchCompletion::Complete);
        assert!(result.metagenomes.is_empty());
        let stats = engine.batch_stats();
        assert_eq!(stats.query_context_associations, associations);
        assert_eq!(stats.query_core_occurrences, cores);
        assert_eq!(stats.query_executed_context_associations, 0);
        assert_eq!(stats.nested_context_calls, 0);
        assert_eq!(stats.query_distinct_cores, 1);
        assert_eq!(stats.core_lookup_tasks, 1);
        let reads = engine.shared_read_stats().unwrap();
        assert_eq!(reads.core_resolutions_present, 0);
        assert_eq!(reads.core_resolutions_absent, 1);
        assert!(reads.core_key_inspections > 0);
        assert_eq!(reads.core_descriptor_inspections, 0);
        assert_eq!(reads.group_descriptor_inspections, 0);
        assert_eq!(reads.member_descriptor_inspections, 0);
        assert_eq!(reads.physical_positions_decoded, 0);
    }
}

#[test]
fn absent_core_batch_work_is_stable_across_worker_counts() {
    use crate::trace::{SearchCompletion, TraceConfig, TraceEngine};

    const ISOLATED: &str = "JAM_SHARED_THREAD_WORK_ISOLATED";
    if std::env::var_os(ISOLATED).is_none() {
        let status = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "shared_tests::absent_core_batch_work_is_stable_across_worker_counts",
            ])
            .env(ISOLATED, "1")
            .status()
            .unwrap();
        assert!(status.success());
        return;
    }

    for compact in [false, true] {
        let (directory, reader, _) = fixture(0);
        drop(reader);
        let source = directory.path().join("fixture.shared");
        let shared = if compact {
            let packed = directory.path().join("packed.shared");
            let compact = directory.path().join("compact.shared");
            crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
            crate::shared_pack::repack_shared_cores(&packed, &compact).unwrap();
            compact
        } else {
            source
        };
        let queries = (0..4)
            .map(|ordinal| (format!("absent-{ordinal}"), vec![b'A'; 2000]))
            .collect::<Vec<_>>();
        let config = TraceConfig {
            use_sketch: false,
            ..TraceConfig::default()
        };
        let mut expected = None;
        for threads in [1, 4, 8, 16] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let (results, work) = pool.install(|| {
                let engine = TraceEngine::open_shared_observed(&shared, None, true).unwrap();
                let results = engine.search_batch(&queries, config).unwrap();
                let stats = engine.batch_stats();
                let reads = engine.shared_read_stats().unwrap();
                (
                    results,
                    (
                        stats.core_lookup_tasks,
                        stats.query_distinct_cores,
                        stats.query_core_occurrences,
                        stats.query_context_associations,
                        stats.query_executed_context_associations,
                        stats.nested_context_calls,
                        reads.core_resolutions_present,
                        reads.core_resolutions_absent,
                        reads.grouped_core_rows,
                    ),
                )
            });
            assert!(results.iter().all(|result| {
                result.completion == SearchCompletion::Complete && result.metagenomes.is_empty()
            }));
            if let Some((expected_results, expected_work)) = &expected {
                assert_eq!(&results, expected_results, "thread count {threads}");
                assert_eq!(&work, expected_work, "thread count {threads}");
            } else {
                expected = Some((results, work));
            }
        }
    }
}

#[test]
fn resolved_core_scope_and_short_evidence_are_checked() {
    use crate::trace_batch::prepare_cores;
    use crate::trace_index::{TraceIndex, TraceSeed};

    let (directory, reader, _) = fixture(0);
    drop(reader);
    let path = directory.path().join("fixture.shared");
    let first = TraceIndex::Shared(Box::new(SharedReader::open_observed(&path).unwrap()));
    let second = TraceIndex::Shared(Box::new(SharedReader::open_observed(&path).unwrap()));
    let unadmitted = std::iter::once(TARGET_CORE).inspect(|_| {
        panic!("unadmitted core requests were consumed");
    });
    assert!(
        prepare_cores(
            &first,
            unadmitted,
            crate::trace::LOOKUP_CACHE_BYTES / std::mem::size_of::<u32>() + 1,
            true,
        )
        .unwrap()
        .is_none()
    );
    let cores = prepare_cores(&first, [TARGET_CORE], 1, true)
        .unwrap()
        .unwrap();
    let core = SharedKey::core(TARGET_CORE).packed().unwrap();
    let absent_long = SharedKey {
        core: TARGET_CORE,
        context: TARGET_CONTEXT ^ 2,
        length: 31,
    }
    .packed()
    .unwrap();
    assert!(
        second
            .find_seeds_in_cores(&[core, absent_long], &cores)
            .is_err()
    );

    let before = match &first {
        TraceIndex::Shared(reader) => reader.stats(),
        _ => unreachable!(),
    };
    let seeds = first
        .find_seeds_in_cores(&[absent_long, core], &cores)
        .unwrap();
    assert!(seeds[0].is_none());
    assert!(matches!(seeds[1], Some(TraceSeed::Shared(_))));
    let after = match &first {
        TraceIndex::Shared(reader) => reader.stats(),
        _ => unreachable!(),
    };
    assert_eq!(
        after.core_resolutions_present,
        before.core_resolutions_present
    );
    assert_eq!(
        after.core_resolutions_absent,
        before.core_resolutions_absent
    );
    assert_eq!(
        after.core_descriptor_inspections,
        before.core_descriptor_inspections + 1
    );
}

fn grouped_lookup_fixture(preceding: u32) -> (tempfile::TempDir, std::path::PathBuf) {
    let directory = tempfile::tempdir().unwrap();
    let jidx = directory.path().join("grouped-metadata.jidx");
    let mut writer = JidxWriter::new(
        &jidx,
        &JidxInput {
            k: 15,
            rescue_k15: false,
            minimizer_window: 64,
            jam_sha256: [21; 32],
            manifest_sha256: [22; 32],
        },
    )
    .unwrap();
    for id in 0..3 {
        writer
            .begin_metagenome(MetagenomeInput {
                name: format!("grouped-{id}"),
                bgzf_uri: format!("grouped-{id}.bgz"),
                bgzf_bytes: 100,
                bgzf_sha256: [id as u8 + 23; 32],
                gzi: vec![0; 8],
            })
            .unwrap();
        writer
            .begin_contig(ContigInput {
                name: format!("grouped-{id}-contig"),
                length: 100_000,
                fasta_offset: 4,
                line_bases: 80,
                line_width: 81,
            })
            .unwrap();
    }
    writer.finish().unwrap();
    let reference = JidxReader::open(&jidx).unwrap();
    let mut seeds = (0..preceding)
        .map(|core| IndexedSeed {
            member: 0,
            contig: 0,
            seed: SharedSeed {
                core,
                context: 0,
                flags: 0,
                position: u64::from(core % 90_000),
            },
        })
        .collect::<Vec<_>>();
    seeds.extend([
        indexed_target(0, 0, 100, false),
        indexed_target(0, 0, 200, true),
        IndexedSeed {
            member: 2,
            contig: 2,
            seed: SharedSeed {
                core: TARGET_CORE,
                context: TARGET_CONTEXT ^ 1,
                flags: 6,
                position: 300,
            },
        },
        IndexedSeed {
            member: 1,
            contig: 1,
            seed: SharedSeed {
                core: TARGET_CORE,
                context: TARGET_CONTEXT & !((1 << 20) - 1),
                flags: 2,
                position: 400,
            },
        },
        IndexedSeed {
            member: 2,
            contig: 2,
            seed: SharedSeed {
                core: TARGET_CORE,
                context: 0,
                flags: 0,
                position: 5,
            },
        },
    ]);
    let shared = directory.path().join("grouped.shared");
    write_shared_index(&reference, &shared, 64, &mut seeds).unwrap();
    (directory, shared)
}

#[allow(clippy::type_complexity)]
fn group_evidence(
    reader: &SharedReader,
    group: Option<SharedGroup>,
) -> Option<(u32, u64, Vec<(u32, Vec<(u32, u64, bool)>)>)> {
    let group = group?;
    let members = reader
        .members(group)
        .unwrap()
        .into_iter()
        .map(|member| {
            let occurrences = reader
                .member_occurrences(group, member)
                .unwrap()
                .into_iter()
                .map(|occurrence| {
                    (
                        occurrence.contig_id,
                        occurrence.position,
                        occurrence.canonical_orientation,
                    )
                })
                .collect();
            (member.metagenome_id, occurrences)
        })
        .collect();
    Some((group.member_count(), group.occurrence_count(), members))
}

#[test]
fn grouped_lookup_matches_scalar_with_contexts_absence_and_chunking() {
    for preceding in [0, 4096] {
        let (_directory, path) = grouped_lookup_fixture(preceding);
        let core = SharedKey::core(TARGET_CORE);
        let context_21 = SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT >> 20,
            length: 21,
        };
        let context_31 = SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT,
            length: 31,
        };
        let different_outer = SharedKey {
            context: TARGET_CONTEXT ^ 1,
            ..context_31
        };
        let absent_long = SharedKey {
            context: TARGET_CONTEXT ^ 2,
            ..context_31
        };
        let absent_core = SharedKey::core(CORE_LIMIT - 1);
        let requests = [
            different_outer,
            core,
            absent_core,
            context_21,
            context_31,
            absent_long,
            core,
            context_31,
            different_outer,
            absent_core,
            context_21,
            core,
            absent_long,
            context_31,
            core,
            different_outer,
            context_21,
        ];
        let scalar = SharedReader::open(&path).unwrap();
        let expected = requests
            .iter()
            .map(|&key| {
                let group = scalar.find(key).unwrap();
                group_evidence(&scalar, group)
            })
            .collect::<Vec<_>>();
        assert_eq!(expected[1].as_ref().unwrap().1, 5);
        assert_eq!(expected[3].as_ref().unwrap().1, 4);
        assert_eq!(expected[4].as_ref().unwrap().1, 2);
        assert_eq!(expected[0].as_ref().unwrap().1, 1);
        assert!(expected[2].is_none());
        assert!(expected[5].is_none());

        for chunk_size in [1, 3, 17] {
            let reader = SharedReader::open_observed(&path).unwrap();
            let mut groups = Vec::new();
            for chunk in requests.chunks(chunk_size) {
                groups.extend(reader.find_many(chunk).unwrap());
            }
            let stats = reader.stats();
            assert_eq!(
                stats.grouped_core_rows,
                stats.core_resolutions_present + stats.grouped_core_rows_without_match
            );
            let expected_present = requests
                .chunks(chunk_size)
                .filter(|chunk| chunk.iter().any(|key| key.core == TARGET_CORE))
                .count() as u64;
            let expected_absent = requests
                .chunks(chunk_size)
                .filter(|chunk| chunk.iter().any(|key| key.core == CORE_LIMIT - 1))
                .count() as u64;
            assert_eq!(stats.core_resolutions_present, expected_present);
            assert_eq!(stats.core_resolutions_absent, expected_absent);
            let actual = groups
                .into_iter()
                .map(|group| group_evidence(&reader, group))
                .collect::<Vec<_>>();
            assert_eq!(actual, expected, "chunk size {chunk_size}");

            if chunk_size == requests.len() {
                let unique = [
                    core,
                    context_21,
                    context_31,
                    different_outer,
                    absent_long,
                    absent_core,
                ];
                let distinct_reader = SharedReader::open_observed(&path).unwrap();
                distinct_reader.find_many(&unique).unwrap();
                let distinct = distinct_reader.stats();
                assert_eq!(stats.core_resolutions_present, 1);
                assert_eq!(stats.core_resolutions_absent, 1);
                assert_eq!(
                    stats.core_descriptor_inspections,
                    distinct.core_descriptor_inspections
                );
                assert_eq!(
                    stats.group_descriptor_inspections,
                    distinct.group_descriptor_inspections
                );
                assert!(stats.context_comparisons >= distinct.context_comparisons);
            }
        }
    }
}

const CORE_LIMIT: u32 = 1 << 30;

fn member_prefix_fixture(preceding: u32) -> (tempfile::TempDir, SharedReader, u32, u64) {
    let directory = tempfile::tempdir().unwrap();
    let jidx = directory.path().join("members.jidx");
    let mut writer = JidxWriter::new(
        &jidx,
        &JidxInput {
            k: 15,
            rescue_k15: false,
            minimizer_window: 64,
            jam_sha256: [11; 32],
            manifest_sha256: [12; 32],
        },
    )
    .unwrap();
    for id in 0..=preceding {
        let name = if id == preceding {
            "target".to_owned()
        } else {
            format!("member-{id:05}")
        };
        writer
            .begin_metagenome(MetagenomeInput {
                name: name.clone(),
                bgzf_uri: format!("{name}.bgz"),
                bgzf_bytes: 100,
                bgzf_sha256: [13; 32],
                gzi: vec![0; 8],
            })
            .unwrap();
        writer
            .begin_contig(ContigInput {
                name: format!("{name}-contig"),
                length: 10_000,
                fasta_offset: 4,
                line_bases: 80,
                line_width: 81,
            })
            .unwrap();
    }
    writer.finish().unwrap();
    let reference = JidxReader::open(&jidx).unwrap();
    let common_core = 500_000;
    let mut seeds = (0..preceding)
        .map(|id| IndexedSeed {
            member: id,
            contig: id,
            seed: SharedSeed {
                core: common_core,
                context: 0,
                flags: 0,
                position: 20,
            },
        })
        .collect::<Vec<_>>();
    seeds.extend((0..5000).map(|offset| IndexedSeed {
        member: preceding,
        contig: preceding,
        seed: SharedSeed {
            core: common_core,
            context: 0,
            flags: 0,
            position: 1000 + offset,
        },
    }));
    let shared = directory.path().join("members.shared");
    write_shared_index(&reference, &shared, 64, &mut seeds).unwrap();
    (
        directory,
        SharedReader::open_observed(shared).unwrap(),
        common_core,
        5999,
    )
}

fn compact_member_prefix_fixture(preceding: u32) -> (tempfile::TempDir, SharedReader, u32, u64) {
    let (directory, reader, core, last_position) = member_prefix_fixture(preceding);
    drop(reader);
    let source = directory.path().join("members.shared");
    let packed = directory.path().join("members-packed.shared");
    let compact = directory.path().join("members-compact.shared");
    crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
    crate::shared_pack::repack_shared_cores(&packed, &compact).unwrap();
    (
        directory,
        SharedReader::open_observed(compact).unwrap(),
        core,
        last_position,
    )
}

fn singleton_posting_fixture(count: u32) -> (tempfile::TempDir, std::path::PathBuf, Vec<u32>) {
    let directory = tempfile::tempdir().unwrap();
    let metadata = directory.path().join("singletons.jidx");
    let mut writer = JidxWriter::new(
        &metadata,
        &JidxInput {
            k: 15,
            rescue_k15: false,
            minimizer_window: 64,
            jam_sha256: [51; 32],
            manifest_sha256: [52; 32],
        },
    )
    .unwrap();
    writer
        .begin_metagenome(MetagenomeInput {
            name: "singletons".into(),
            bgzf_uri: "singletons.bgz".into(),
            bgzf_bytes: 100,
            bgzf_sha256: [53; 32],
            gzi: vec![0; 8],
        })
        .unwrap();
    writer
        .begin_contig(ContigInput {
            name: "singletons".into(),
            length: u64::from(count) + 100,
            fasta_offset: 4,
            line_bases: 80,
            line_width: 81,
        })
        .unwrap();
    writer.finish().unwrap();
    let reference = JidxReader::open(&metadata).unwrap();
    let cores = (0..count).collect::<Vec<_>>();
    let mut seeds = cores
        .iter()
        .map(|&core| IndexedSeed {
            member: 0,
            contig: 0,
            seed: SharedSeed {
                core,
                context: 0,
                flags: 0,
                position: u64::from(core) + 32,
            },
        })
        .collect::<Vec<_>>();
    let source = directory.path().join("singletons.shared");
    let packed = directory.path().join("singletons-packed.shared");
    let compact = directory.path().join("singletons-compact.shared");
    write_shared_index(&reference, &source, 64, &mut seeds).unwrap();
    crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
    crate::shared_pack::repack_shared_cores(&packed, &compact).unwrap();
    (directory, compact, cores)
}

#[allow(clippy::type_complexity)]
fn serial_posting_evidence(
    reader: &SharedReader,
    entries: &[(u64, Option<crate::trace_index::TraceSeed>)],
) -> Vec<(u64, Vec<(u32, Vec<(u32, u64, bool)>)>)> {
    entries
        .iter()
        .filter_map(|&(key, seed)| {
            let crate::trace_index::TraceSeed::Shared(group) = seed? else {
                unreachable!()
            };
            let members = reader
                .members(group)
                .unwrap()
                .into_iter()
                .map(|member| {
                    let occurrences = reader
                        .member_occurrences(group, member)
                        .unwrap()
                        .into_iter()
                        .map(|value| (value.contig_id, value.position, value.canonical_orientation))
                        .collect();
                    (member.metagenome_id, occurrences)
                })
                .collect();
            Some((key, members))
        })
        .collect()
}

#[allow(clippy::type_complexity)]
fn retained_posting_evidence(
    entries: &[(u64, Option<crate::trace_index::TraceSeed>)],
    postings: &[Option<crate::trace_batch::BatchPosting>],
) -> Vec<(u64, Vec<(u32, Vec<(u32, u64, bool)>)>)> {
    entries
        .iter()
        .zip(postings)
        .filter_map(|(&(key, _), posting)| {
            let posting = posting.as_ref()?;
            let occurrences = posting.occurrences.as_ref().unwrap();
            Some((
                key,
                posting
                    .documents
                    .iter()
                    .zip(occurrences)
                    .map(|(document, values)| {
                        let values = values
                            .iter()
                            .map(|value| {
                                (value.contig_id, value.position, value.canonical_orientation)
                            })
                            .collect();
                        (document.metagenome_id(), values)
                    })
                    .collect(),
            ))
        })
        .collect()
}

#[test]
fn posting_operation_fills_admitted_member_and_occurrence_storage() {
    let (_directory, reader, core, last_position) = compact_member_prefix_fixture(16);
    let group = reader.find(SharedKey::core(core)).unwrap().unwrap();
    let expected_members = reader.members(group).unwrap();
    assert_eq!(expected_members.len(), 17);
    let operation = reader.posting_operation().unwrap();

    let mut empty = Vec::new();
    operation
        .append_member_range(group, 0, 0, &mut empty)
        .unwrap();
    assert!(empty.is_empty());
    let mut partial = Vec::with_capacity(4);
    operation
        .append_member_range(group, 3, 4, &mut partial)
        .unwrap();
    assert_eq!(partial, expected_members[3..7]);
    let mut exact = Vec::with_capacity(expected_members.len());
    operation
        .append_member_range(group, 0, expected_members.len(), &mut exact)
        .unwrap();
    assert_eq!(exact, expected_members);

    let mut unadmitted = Vec::new();
    assert!(matches!(
        operation.append_member_range(group, 0, 1, &mut unadmitted),
        Err(SharedError::ResourceLimit)
    ));
    assert!(unadmitted.is_empty());
    assert!(matches!(
        operation.append_member_range(group, group.member_count(), 1, &mut exact),
        Err(SharedError::Invalid("member range"))
    ));

    let member = *expected_members.last().unwrap();
    let expected = reader.member_occurrences(group, member).unwrap();
    assert_eq!(expected.len(), 5000);
    let blank = crate::jidx_reader::SeedOccurrence {
        contig_id: 0,
        position: 0,
        canonical_orientation: false,
    };
    let mut first = vec![blank; 4096];
    operation
        .fill_occurrence_block(group, member, 0, &mut first)
        .unwrap();
    assert_eq!(first, expected[..4096]);
    let mut tail = vec![blank; 904];
    operation
        .fill_occurrence_block(group, member, 4096, &mut tail)
        .unwrap();
    assert_eq!(tail, expected[4096..]);
    assert_eq!(tail.last().unwrap().position, last_position);
    let mut too_long = vec![blank; 905];
    assert!(matches!(
        operation.fill_occurrence_block(group, member, 4096, &mut too_long),
        Err(SharedError::Invalid("occurrence block"))
    ));
    operation.finish().unwrap();
}

#[test]
fn posting_operation_rejects_foreign_and_reopened_generation_handles() {
    let (directory, reader, core, _) = compact_member_prefix_fixture(2);
    let path = directory.path().join("members-compact.shared");
    let group = reader.find(SharedKey::core(core)).unwrap().unwrap();
    let member = reader.members(group).unwrap()[0];
    let other = SharedReader::open(&path).unwrap();
    let other_operation = other.posting_operation().unwrap();
    let mut members = Vec::with_capacity(1);
    assert!(matches!(
        other_operation.append_member_range(group, 0, 1, &mut members),
        Err(SharedError::Invalid("group handle"))
    ));
    let mut occurrence = [crate::jidx_reader::SeedOccurrence {
        contig_id: 0,
        position: 0,
        canonical_orientation: false,
    }];
    assert!(matches!(
        other_operation.fill_occurrence_block(group, member, 0, &mut occurrence),
        Err(SharedError::Invalid("group handle"))
    ));
    other_operation.finish().unwrap();

    let stale = reader.posting_operation().unwrap();
    mutate_and_resign(&path, |bytes, header| {
        let occurrence = header.section(Section::Occurrences).offset as usize;
        let flags = crate::jidx::read_u32(bytes, occurrence + 8) ^ 1;
        bytes[occurrence + 8..occurrence + 12].copy_from_slice(&flags.to_le_bytes());
    });
    assert!(matches!(stale.finish(), Err(SharedError::SourceChanged)));
    let reopened = SharedReader::open(&path).unwrap();
    let reopened_operation = reopened.posting_operation().unwrap();
    assert!(matches!(
        reopened_operation.append_member_range(group, 0, 1, &mut members),
        Err(SharedError::Invalid("group handle"))
    ));
    reopened_operation.finish().unwrap();
}

#[test]
fn posting_operation_keeps_only_valid_member_prefix_on_corruption() {
    for (name, corrupt_ordinal, expected_prefix) in [("early", 0usize, 0), ("late", 8, 8)] {
        let (directory, reader, core, _) = compact_member_prefix_fixture(16);
        drop(reader);
        let path = directory.path().join("members-compact.shared");
        mutate_and_resign(&path, |bytes, header| {
            let width = header.row_bytes(Section::Members) as usize;
            let count = header.id_bytes() + 4;
            let at = header.section(Section::Members).offset as usize + corrupt_ordinal * width;
            bytes[at + count..at + count + 4].fill(0);
        });
        let reader = SharedReader::open(&path).unwrap();
        let group = reader.find(SharedKey::core(core)).unwrap().unwrap();
        let operation = reader.posting_operation().unwrap();
        let mut members = Vec::with_capacity(group.member_count() as usize);
        assert!(matches!(
            operation.append_member_range(group, 0, group.member_count() as usize, &mut members),
            Err(SharedError::Invalid("member row"))
        ));
        assert_eq!(members.len(), expected_prefix, "{name}");
        assert_eq!(
            members
                .iter()
                .map(|member| member.metagenome_id)
                .collect::<Vec<_>>(),
            (0..expected_prefix as u32).collect::<Vec<_>>(),
            "{name}"
        );
        operation.finish().unwrap();
    }
}

#[test]
fn bounded_posting_executor_is_stable_and_parallel_for_many_singletons() {
    let (_directory, path, cores) = singleton_posting_fixture(2048);
    let mut expected_plan = None;
    for threads in [1, 4, 8, 16] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let (actual, plan) = pool.install(|| {
            let reader = SharedReader::open_observed(&path).unwrap();
            let keys = cores
                .iter()
                .copied()
                .map(SharedKey::core)
                .collect::<Vec<_>>();
            let entries = keys
                .iter()
                .zip(reader.find_many(&keys).unwrap())
                .map(|(key, group)| {
                    (
                        key.packed().unwrap(),
                        group.map(crate::trace_index::TraceSeed::Shared),
                    )
                })
                .collect::<Vec<_>>();
            let expected = serial_posting_evidence(&reader, &entries);
            let mut postings = (0..entries.len()).map(|_| None).collect::<Vec<_>>();
            let execution = crate::trace_postings::prepare_shared_postings(
                &reader,
                &entries,
                &mut postings,
                128 * 1024 * 1024,
                true,
            )
            .unwrap()
            .unwrap();
            assert!(execution.complete);
            assert_eq!(execution.admitted_member_rows, 2048);
            assert_eq!(execution.admitted_position_rows, 2048);
            assert_eq!(execution.member_copies, 2048);
            assert_eq!(
                execution.peak_bytes,
                execution.scratch_bytes + execution.retained_bytes
            );
            assert!(execution.member_tasks > 1 && execution.position_tasks > 1);
            if threads >= 4 {
                assert!(execution.peak_parallel_tasks > 1, "{threads} threads");
            }
            let actual = retained_posting_evidence(&entries, &postings);
            assert_eq!(actual, expected);
            let plan = (
                execution.member_tasks,
                execution.position_tasks,
                execution.task_hash,
                execution.retained_bytes,
                execution.peak_bytes,
                execution.scratch_bytes,
                execution.admitted_member_rows,
                execution.admitted_position_rows,
            );
            (actual, plan)
        });
        if let Some(expected) = expected_plan {
            assert_eq!(plan, expected, "{threads} threads");
        } else {
            expected_plan = Some(plan);
        }
        assert_eq!(actual.len(), 2048);
    }

    let reader = SharedReader::open(&path).unwrap();
    let keys = cores
        .iter()
        .copied()
        .map(SharedKey::core)
        .collect::<Vec<_>>();
    let entries = keys
        .iter()
        .zip(reader.find_many(&keys).unwrap())
        .map(|(key, group)| {
            (
                key.packed().unwrap(),
                group.map(crate::trace_index::TraceSeed::Shared),
            )
        })
        .collect::<Vec<_>>();
    let mut refused = (0..entries.len()).map(|_| None).collect::<Vec<_>>();
    assert!(
        crate::trace_postings::prepare_shared_postings(&reader, &entries, &mut refused, 0, false)
            .unwrap()
            .is_none()
    );
    assert!(refused.iter().all(Option::is_none));
    let mut empty = [];
    assert!(
        crate::trace_postings::prepare_shared_postings(&reader, &[], &mut empty, 0, false)
            .unwrap()
            .unwrap()
            .complete
    );
}

#[test]
fn bounded_posting_executor_matches_context_and_long_member_streams() {
    let (directory, source) = grouped_lookup_fixture(0);
    let packed = directory.path().join("postings-packed.shared");
    let compact = directory.path().join("postings-compact.shared");
    crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
    crate::shared_pack::repack_shared_cores(&packed, &compact).unwrap();
    let reader = SharedReader::open(&compact).unwrap();
    let mut keys = [
        SharedKey::core(TARGET_CORE),
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT >> 20,
            length: 21,
        },
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT,
            length: 31,
        },
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT ^ 1,
            length: 31,
        },
    ];
    keys.sort_unstable_by_key(|key| key.packed().unwrap());
    let entries = keys
        .iter()
        .zip(reader.find_many(&keys).unwrap())
        .map(|(key, group)| {
            (
                key.packed().unwrap(),
                group.map(crate::trace_index::TraceSeed::Shared),
            )
        })
        .collect::<Vec<_>>();
    let member_counts = entries
        .iter()
        .filter_map(|entry| match entry.1 {
            Some(crate::trace_index::TraceSeed::Shared(group)) => Some(group.member_count()),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert!(member_counts.contains(&1) && member_counts.iter().any(|&count| count > 1));
    let expected = serial_posting_evidence(&reader, &entries);
    let mut postings = (0..entries.len()).map(|_| None).collect::<Vec<_>>();
    let execution = crate::trace_postings::prepare_shared_postings(
        &reader,
        &entries,
        &mut postings,
        128 * 1024 * 1024,
        true,
    )
    .unwrap()
    .unwrap();
    assert!(execution.complete);
    assert_eq!(retained_posting_evidence(&entries, &postings), expected);

    let (_directory, reader, core, last_position) = compact_member_prefix_fixture(0);
    let group = reader.find(SharedKey::core(core)).unwrap().unwrap();
    assert_eq!((group.member_count(), group.occurrence_count()), (1, 5000));
    let entries = [(
        SharedKey::core(core).packed().unwrap(),
        Some(crate::trace_index::TraceSeed::Shared(group)),
    )];
    let expected = serial_posting_evidence(&reader, &entries);
    let mut postings = [None];
    let execution = crate::trace_postings::prepare_shared_postings(
        &reader,
        &entries,
        &mut postings,
        128 * 1024 * 1024,
        true,
    )
    .unwrap()
    .unwrap();
    assert!(execution.complete);
    assert_eq!(execution.admitted_position_rows, 5000);
    let actual = retained_posting_evidence(&entries, &postings);
    assert_eq!(actual, expected);
    assert_eq!(actual[0].1[0].1.last().unwrap().1, last_position);
}

#[test]
fn bounded_posting_executor_selects_earliest_error_and_publishes_nothing() {
    let (directory, reader, _) = fixture(16);
    drop(reader);
    let source = directory.path().join("fixture.shared");
    let packed = directory.path().join("postings-packed.shared");
    let compact = directory.path().join("postings-compact.shared");
    crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
    crate::shared_pack::repack_shared_cores(&packed, &compact).unwrap();
    mutate_and_resign(&compact, |bytes, header| {
        let occurrences = header.section(Section::Occurrences);
        let occurrence_width = header.row_bytes(Section::Occurrences) as usize;
        let mut changed_position = false;
        for row in bytes
            [occurrences.offset as usize..(occurrences.offset + occurrences.length) as usize]
            .chunks_exact_mut(occurrence_width)
        {
            if crate::jidx::read_u64(row, 16) == 0 {
                let flags = crate::jidx::read_u32(row, 8) | 8;
                row[8..12].copy_from_slice(&flags.to_le_bytes());
                changed_position = true;
                break;
            }
        }
        assert!(changed_position);

        let members = header.section(Section::Members);
        let member_width = header.row_bytes(Section::Members) as usize;
        let count = header.id_bytes() + 4;
        let mut changed_member = false;
        for row in bytes[members.offset as usize..(members.offset + members.length) as usize]
            .chunks_exact_mut(member_width)
        {
            if row[0] == 1 {
                row[count..count + 4].fill(0);
                changed_member = true;
                break;
            }
        }
        assert!(changed_member);
    });
    let reader = SharedReader::open(&compact).unwrap();
    let keys = [SharedKey::core(0), SharedKey::core(TARGET_CORE)];
    let entries = keys
        .iter()
        .zip(reader.find_many(&keys).unwrap())
        .map(|(key, group)| {
            (
                key.packed().unwrap(),
                group.map(crate::trace_index::TraceSeed::Shared),
            )
        })
        .collect::<Vec<_>>();
    let mut postings = [None, None];
    assert!(matches!(
        crate::trace_postings::prepare_shared_postings(
            &reader,
            &entries,
            &mut postings,
            128 * 1024 * 1024,
            true,
        ),
        Err(crate::trace::TraceError::Shared(SharedError::Invalid(
            "occurrence row"
        )))
    ));
    assert!(postings.iter().all(Option::is_none));
}

#[test]
fn bounded_posting_executor_rejects_foreign_and_stale_entries() {
    let (directory, reader, _) = fixture(0);
    drop(reader);
    let source = directory.path().join("fixture.shared");
    let packed = directory.path().join("postings-packed.shared");
    let compact = directory.path().join("postings-compact.shared");
    crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
    crate::shared_pack::repack_shared_cores(&packed, &compact).unwrap();
    let first = SharedReader::open(&compact).unwrap();
    let group = first.find(SharedKey::core(TARGET_CORE)).unwrap().unwrap();
    let entries = [(
        SharedKey::core(TARGET_CORE).packed().unwrap(),
        Some(crate::trace_index::TraceSeed::Shared(group)),
    )];
    let second = SharedReader::open(&compact).unwrap();
    let mut postings = [None];
    assert!(matches!(
        crate::trace_postings::prepare_shared_postings(
            &second,
            &entries,
            &mut postings,
            128 * 1024 * 1024,
            false,
        ),
        Err(crate::trace::TraceError::Shared(SharedError::Invalid(
            "group handle"
        )))
    ));
    assert!(postings.iter().all(Option::is_none));

    mutate_and_resign(&compact, |bytes, header| {
        let occurrence = header.section(Section::Occurrences).offset as usize;
        let flags = crate::jidx::read_u32(bytes, occurrence + 8) ^ 1;
        bytes[occurrence + 8..occurrence + 12].copy_from_slice(&flags.to_le_bytes());
    });
    assert!(matches!(
        crate::trace_postings::prepare_shared_postings(
            &first,
            &entries,
            &mut postings,
            128 * 1024 * 1024,
            false,
        ),
        Err(crate::trace::TraceError::Shared(SharedError::SourceChanged))
    ));
    assert!(postings.iter().all(Option::is_none));
}

#[test]
fn selected_member_and_late_occurrence_do_not_replay_prefixes() {
    for preceding in [16u32, 4096, 16_384] {
        let (_directory, reader, core, last_position) = member_prefix_fixture(preceding);
        let group = reader.find(SharedKey::core(core)).unwrap().unwrap();
        let before = reader.stats();
        let member = reader.member(group, preceding).unwrap().unwrap();
        let work = delta(reader.stats(), before);
        let logarithmic_bound = u64::from(u32::BITS - (preceding + 1).leading_zeros()) + 1;
        assert!(work[2] <= logarithmic_bound);
        assert_eq!(&work[3..], &[0, 0]);
        assert_eq!(member.occurrence_count(), 5000);
        assert!(reader.occurrence_block(group, member, 0, 4097).is_err());
        assert!(reader.occurrence_block(group, member, 5001, 1).is_err());
        assert_eq!(
            reader.metagenome(preceding).unwrap().unwrap().name,
            "target"
        );
        assert_eq!(
            reader.contig(preceding).unwrap().unwrap().name,
            "target-contig"
        );

        let before = reader.stats();
        let occurrence = reader.occurrence_block(group, member, 4999, 1).unwrap();
        assert_eq!(
            (occurrence[0].contig_id, occurrence[0].position),
            (preceding, last_position)
        );
        assert_eq!(delta(reader.stats(), before), [0, 0, 0, 1, 1]);
    }
}

fn mutate_and_resign(path: &std::path::Path, mutate: impl FnOnce(&mut [u8], &SharedHeader)) {
    let mut bytes = std::fs::read(path).unwrap();
    let mut header = SharedHeader::decode(&bytes[..HEADER_BYTES], bytes.len() as u64).unwrap();
    mutate(&mut bytes, &header);
    let checksum_start = header.section(Section::Checksums).offset;
    let data_pages = checksum_start / PAGE_BYTES - 1;
    let mut hashes = (0..data_pages)
        .map(|page| {
            let start = ((page + 1) * PAGE_BYTES) as usize;
            sha256(&bytes[start..start + PAGE_BYTES as usize])
        })
        .collect::<Vec<_>>();
    for level in checksum_layout(data_pages).unwrap() {
        let mut parents = Vec::new();
        for (page, chunk) in hashes.chunks(128).enumerate() {
            let mut contents = [0; PAGE_BYTES as usize];
            for (index, hash) in chunk.iter().enumerate() {
                contents[index * 32..index * 32 + 32].copy_from_slice(hash);
            }
            let start = (checksum_start + level.offset + page as u64 * PAGE_BYTES) as usize;
            bytes[start..start + PAGE_BYTES as usize].copy_from_slice(&contents);
            parents.push(sha256(&contents));
        }
        hashes = parents;
    }
    header.checksum_root_sha256 = hashes[0];
    header.body_sha256 = sha256(&bytes[HEADER_BYTES..]);
    bytes[..HEADER_BYTES].copy_from_slice(&header.encode().unwrap());
    std::fs::write(path, bytes).unwrap();
}

#[test]
fn authenticated_context_and_reference_corruption_is_rejected() {
    let (directory, reader, _) = fixture(0);
    drop(reader);
    let path = directory.path().join("fixture.shared");
    let native = SharedReader::open(&path).unwrap().stats();
    assert!(!native.observed);
    assert!(!native.file.observed);
    mutate_and_resign(&path, |bytes, header| {
        let occurrence = header.section(Section::Occurrences).offset as usize;
        bytes[occurrence..occurrence + 4]
            .copy_from_slice(&(TARGET_CONTEXT ^ (1 << 20)).to_le_bytes());
    });
    let reader = SharedReader::open_observed(&path).unwrap();
    let key = SharedKey {
        core: TARGET_CORE,
        context: TARGET_CONTEXT >> 20,
        length: 21,
    };
    let groups = reader.find_many(&[key, key]).unwrap();
    assert_eq!(groups[0], groups[1]);
    let group = groups[0].unwrap();
    let member = reader.member(group, 1).unwrap().unwrap();
    assert!(matches!(
        reader.member_occurrences(group, member),
        Err(SharedError::Invalid("occurrence context"))
    ));

    let (directory, reader, _) = fixture(0);
    drop(reader);
    let path = directory.path().join("fixture.shared");
    mutate_and_resign(&path, |bytes, header| {
        let reference = header.section(Section::References).offset as usize;
        bytes[reference..reference + 8].copy_from_slice(&u64::MAX.to_le_bytes());
    });
    let reader = SharedReader::open_observed(path).unwrap();
    let key = SharedKey::core(TARGET_CORE);
    let groups = reader.find_many(&[key, key]).unwrap();
    assert_eq!(groups[0], groups[1]);
    let group = groups[0].unwrap();
    let member = reader.member(group, 1).unwrap().unwrap();
    assert!(matches!(
        reader.member_occurrences(group, member),
        Err(SharedError::Invalid("occurrence reference"))
    ));
}

#[test]
fn resealed_malformed_rows_fail_scalar_and_grouped_lookup() {
    let (directory, reader, _) = fixture(0);
    let key = SharedKey::core(TARGET_CORE);
    let core_ordinal = reader.find(key).unwrap().unwrap().core_ordinal();
    drop(reader);
    let path = directory.path().join("fixture.shared");
    mutate_and_resign(&path, |bytes, header| {
        let core = header.section(Section::Cores).offset as usize
            + core_ordinal as usize * header.row_bytes(Section::Cores) as usize;
        let tagged = crate::jidx::read_u32(bytes, core) | (1 << 30);
        bytes[core..core + 4].copy_from_slice(&tagged.to_le_bytes());
    });
    let reader = SharedReader::open(&path).unwrap();
    assert!(matches!(
        reader.find(key),
        Err(SharedError::Invalid("core row"))
    ));
    assert!(matches!(
        reader.find_many(&[key, key]),
        Err(SharedError::Invalid("core row"))
    ));

    let (directory, reader, _) = fixture(0);
    drop(reader);
    let path = directory.path().join("fixture.shared");
    mutate_and_resign(&path, |bytes, header| {
        let core = header.section(Section::Cores).offset as usize;
        bytes[core + 4..core + 8].fill(0);
    });
    let key = SharedKey::core(TARGET_CORE);
    let reader = SharedReader::open(&path).unwrap();
    assert!(matches!(
        reader.find(key),
        Err(SharedError::Invalid("repeated core"))
    ));
    assert!(matches!(
        reader.find_many(&[key, key]),
        Err(SharedError::Invalid("repeated core"))
    ));

    let (directory, reader, _) = fixture(0);
    drop(reader);
    let path = directory.path().join("fixture.shared");
    mutate_and_resign(&path, |bytes, header| {
        let groups = header.section(Section::Groups);
        let start = groups.offset as usize;
        let end = (groups.offset + groups.length) as usize;
        for group in (start..end).step_by(32) {
            bytes[group + 28..group + 32].copy_from_slice(&1u32.to_le_bytes());
        }
    });
    let reader = SharedReader::open(&path).unwrap();
    assert!(matches!(
        reader.find(key),
        Err(SharedError::Invalid("group row"))
    ));
    assert!(matches!(
        reader.find_many(&[key, key]),
        Err(SharedError::Invalid("group row"))
    ));
}

#[test]
fn packed_groups_preserve_contexts_counts_and_late_direct_access() {
    for preceding in [0, 16_384] {
        let (directory, reference, _) = fixture(preceding);
        let source = directory.path().join("fixture.shared");
        let packed = directory.path().join("packed.shared");
        let stats = crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
        assert!(stats.build.index_bytes <= std::fs::metadata(&source).unwrap().len());
        assert_eq!(stats.logical_references, u64::from(preceding) * 2 + 9);
        let reader = SharedReader::open_observed(&packed).unwrap();
        reader.verify_checksum().unwrap();
        for key in [
            SharedKey::core(TARGET_CORE),
            SharedKey {
                core: TARGET_CORE,
                context: TARGET_CONTEXT >> 20,
                length: 21,
            },
            SharedKey {
                core: TARGET_CORE,
                context: TARGET_CONTEXT,
                length: 31,
            },
            SharedKey::core(TARGET_CORE + 1),
            SharedKey::core(CORE_LIMIT - 1),
        ] {
            assert_eq!(
                group_evidence(&reader, reader.find(key).unwrap()),
                group_evidence(&reference, reference.find(key).unwrap())
            );
        }
        let before = reader.stats();
        let group = reader.find(SharedKey::core(TARGET_CORE)).unwrap().unwrap();
        assert_eq!(
            reader.stats().physical_positions_decoded,
            before.physical_positions_decoded
        );
        let member = reader.member(group, 2).unwrap().unwrap();
        let before = reader.stats();
        let last = reader.occurrence_block(group, member, 0, 1).unwrap();
        assert_eq!((last[0].contig_id, last[0].position), (2, 300));
        assert_eq!(reader.stats().references_decoded, before.references_decoded);
        assert_eq!(
            reader.stats().physical_positions_decoded - before.physical_positions_decoded,
            1
        );
        assert!(crate::shared_pack::repack_shared_index(&source, &packed).is_err());
        assert!(
            crate::shared_pack::repack_shared_index(&packed, directory.path().join("again"))
                .is_err()
        );
    }
    let (directory, source) = grouped_lookup_fixture(4096);
    let packed = directory.path().join("packed.shared");
    crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
    let reference = SharedReader::open(&source).unwrap();
    let reader = SharedReader::open(&packed).unwrap();
    let keys = [
        SharedKey::core(TARGET_CORE),
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT >> 20,
            length: 21,
        },
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT,
            length: 31,
        },
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT ^ 1,
            length: 31,
        },
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT ^ 2,
            length: 31,
        },
        SharedKey::core(CORE_LIMIT - 1),
    ];
    for chunk in [1, 3, 6] {
        for keys in keys.chunks(chunk) {
            let actual = reader.find_many(keys).unwrap();
            for (&key, group) in keys.iter().zip(actual) {
                assert_eq!(
                    group_evidence(&reader, group),
                    group_evidence(&reference, reference.find(key).unwrap())
                );
            }
        }
    }
    let retained = reader.find(keys[0]).unwrap().unwrap();
    mutate_and_resign(&packed, |bytes, header| {
        let at = header.section(Section::Occurrences).offset as usize;
        bytes[at + 16] ^= 1;
    });
    assert!(matches!(
        reader.members(retained),
        Err(SharedError::SourceChanged)
    ));
}

#[test]
fn compact_cores_preserve_evidence_and_bound_prefix_search() {
    let (directory, baseline, _) = fixture(16);
    let source = directory.path().join("fixture.shared");
    let packed = directory.path().join("packed.shared");
    let compact = directory.path().join("compact.shared");
    crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
    let stats = crate::shared_pack::repack_shared_cores(&packed, &compact).unwrap();
    assert_eq!(stats.core_payload_bytes, 13);
    assert_eq!(stats.hot_core_bytes, stats.build.core_count * 4);
    assert_eq!(stats.cold_core_bytes, stats.build.core_count * 13);
    assert_eq!(
        stats.core_prefix_bytes,
        crate::shared_format::CORE_PREFIX_BOUNDARIES as u64 * 4
    );

    let reader = SharedReader::open_observed(&compact).unwrap();
    let keys = [
        SharedKey::core(0),
        SharedKey::core(TARGET_CORE),
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT >> 20,
            length: 21,
        },
        SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT,
            length: 31,
        },
        SharedKey::core(TARGET_CORE + 1),
        SharedKey::core(CORE_LIMIT - 1),
    ];
    let expected = keys
        .iter()
        .map(|&key| group_evidence(&baseline, baseline.find(key).unwrap()))
        .collect::<Vec<_>>();
    let before = reader.stats();
    let actual = reader
        .find_many(&keys)
        .unwrap()
        .into_iter()
        .map(|group| group_evidence(&reader, group))
        .collect::<Vec<_>>();
    assert_eq!(actual, expected);
    let after = reader.stats();
    assert_eq!(
        after.core_resolutions_absent - before.core_resolutions_absent,
        1
    );
    assert!(
        after.core_descriptor_inspections - before.core_descriptor_inspections < keys.len() as u64
    );
    reader.verify_checksum().unwrap();
}

#[test]
fn compact_prefix_endpoints_and_cold_payloads_fail_closed() {
    let (directory, reader, _) = fixture(16);
    drop(reader);
    let source = directory.path().join("fixture.shared");
    let packed = directory.path().join("packed.shared");
    let compact = directory.path().join("compact.shared");
    crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
    let malformed_directory = directory.path().join("malformed-directory.shared");
    crate::shared_pack::repack_shared_cores(&packed, &malformed_directory).unwrap();
    mutate_and_resign(&malformed_directory, |bytes, header| {
        let prefix = usize::try_from(TARGET_CORE >> 14).unwrap();
        let directory = header.section(Section::CorePrefixes).offset as usize;
        let high = crate::jidx::read_u32(bytes, directory + (prefix + 1) * 4);
        bytes[directory + prefix * 4..directory + prefix * 4 + 4]
            .copy_from_slice(&high.to_le_bytes());
    });
    let reader = SharedReader::open(&malformed_directory).unwrap();
    assert!(matches!(
        reader.find(SharedKey::core(TARGET_CORE)),
        Err(SharedError::Invalid("core prefix membership"))
    ));

    crate::shared_pack::repack_shared_cores(&packed, &compact).unwrap();
    mutate_and_resign(&compact, |bytes, header| {
        let prefix = usize::try_from(TARGET_CORE >> 14).unwrap();
        let directory = header.section(Section::CorePrefixes).offset as usize;
        bytes[directory + prefix * 4..directory + prefix * 4 + 4]
            .copy_from_slice(&0u32.to_le_bytes());
    });
    assert!(matches!(
        SharedReader::open(&compact),
        Err(SharedError::Invalid("core prefix directory"))
    ));

    let malformed_cold = directory.path().join("malformed-cold.shared");
    crate::shared_pack::repack_shared_cores(&packed, &malformed_cold).unwrap();
    let core_ordinal = SharedReader::open(&malformed_cold)
        .unwrap()
        .find(SharedKey::core(TARGET_CORE))
        .unwrap()
        .unwrap()
        .core_ordinal();
    mutate_and_resign(&malformed_cold, |bytes, header| {
        let cold = header.section(Section::CorePayloads).offset as usize
            + core_ordinal as usize * usize::from(header.core_payload_bytes);
        bytes[cold + 12] = 1;
    });
    let reader = SharedReader::open(&malformed_cold).unwrap();
    assert!(
        reader
            .find(SharedKey::core(CORE_LIMIT - 1))
            .unwrap()
            .is_none()
    );
    assert!(matches!(
        reader.find(SharedKey::core(TARGET_CORE)),
        Err(SharedError::Invalid("repeated core"))
    ));
    assert!(matches!(
        reader.verify_checksum(),
        Err(SharedError::Invalid("repeated core"))
    ));
}

#[test]
fn compact_sorted_core_errors_clear_sparse_output() {
    #[derive(Clone, Copy)]
    enum Corruption {
        Prefix,
        Tag,
        ForeignPrefix,
        Cold(u32),
    }

    let (directory, reader, _) = fixture(16);
    drop(reader);
    let source = directory.path().join("fixture.shared");
    let packed = directory.path().join("packed.shared");
    crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
    for (name, corruption, prior_successes) in [
        ("prefix", Corruption::Prefix, 1),
        ("tag", Corruption::Tag, 0),
        ("foreign-prefix", Corruption::ForeignPrefix, 1),
        ("cold-early", Corruption::Cold(0), 0),
        ("cold-late", Corruption::Cold(TARGET_CORE), 1),
    ] {
        let compact = directory.path().join(format!("corrupt-{name}.shared"));
        crate::shared_pack::repack_shared_cores(&packed, &compact).unwrap();
        let selected = match corruption {
            Corruption::Cold(core) => core,
            _ => TARGET_CORE,
        };
        let ordinal = SharedReader::open(&compact)
            .unwrap()
            .find(SharedKey::core(selected))
            .unwrap()
            .unwrap()
            .core_ordinal();
        mutate_and_resign(&compact, |bytes, header| match corruption {
            Corruption::Prefix => {
                let prefix = usize::try_from(TARGET_CORE >> 14).unwrap();
                let directory = header.section(Section::CorePrefixes).offset as usize;
                let high = crate::jidx::read_u32(bytes, directory + (prefix + 1) * 4);
                bytes[directory + prefix * 4..directory + prefix * 4 + 4]
                    .copy_from_slice(&high.to_le_bytes());
            }
            Corruption::Tag => {
                let hot = header.section(Section::Cores).offset as usize + ordinal as usize * 4;
                let tagged = crate::jidx::read_u32(bytes, hot) | (1 << 30);
                bytes[hot..hot + 4].copy_from_slice(&tagged.to_le_bytes());
            }
            Corruption::ForeignPrefix => {
                let hot = header.section(Section::Cores).offset as usize + ordinal as usize * 4;
                let tagged = crate::jidx::read_u32(bytes, hot);
                let foreign = (tagged & (1 << 31)) | (TARGET_CORE - (1 << 14));
                bytes[hot..hot + 4].copy_from_slice(&foreign.to_le_bytes());
            }
            Corruption::Cold(_) => {
                let cold = header.section(Section::CorePayloads).offset as usize
                    + ordinal as usize * usize::from(header.core_payload_bytes);
                bytes[cold + 12] = 1;
            }
        });
        let reader = SharedReader::open_observed(&compact).unwrap();
        let before = reader.stats();
        let mut output = Vec::new();
        assert!(
            reader
                .resolve_sorted_cores_into(&[0, TARGET_CORE], &mut output)
                .is_err(),
            "{name}"
        );
        let after = reader.stats();
        assert!(output.is_empty(), "{name}");
        assert_eq!(
            after.core_resolutions_present - before.core_resolutions_present,
            prior_successes,
            "{name}"
        );
    }
}

#[test]
fn compact_sorted_core_task_boundaries_repeat_boundary_evidence() {
    let (directory, reader, _) = fixture(16);
    drop(reader);
    let source = directory.path().join("fixture.shared");
    let packed = directory.path().join("packed.shared");
    let compact = directory.path().join("compact.shared");
    crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
    crate::shared_pack::repack_shared_cores(&packed, &compact).unwrap();
    let reader = SharedReader::open_observed(&compact).unwrap();
    let reference = SharedReader::open(&compact).unwrap();
    let mut boundary_hits = 0;
    for task in [&[0, 1, 2][..], &[2, 3, TARGET_CORE, TARGET_CORE + 1][..]] {
        let mut output = Vec::new();
        reader.resolve_sorted_cores_into(task, &mut output).unwrap();
        let expected = reference
            .find_many(
                &task
                    .iter()
                    .copied()
                    .map(SharedKey::core)
                    .collect::<Vec<_>>(),
            )
            .unwrap()
            .into_iter()
            .flatten()
            .map(|group| (group.key(), group.member_count(), group.occurrence_count()))
            .collect::<Vec<_>>();
        let actual = output
            .iter()
            .map(|group| (group.key(), group.member_count(), group.occurrence_count()))
            .collect::<Vec<_>>();
        assert_eq!(actual, expected);
        boundary_hits += output.iter().filter(|group| group.key().core == 2).count();
    }
    assert_eq!(boundary_hits, 2);
    assert!(reader.stats().core_view_creations >= 2);
}

#[test]
fn compact_cold_row_crossing_a_page_preserves_the_selected_core() {
    let (directory, baseline, _) = fixture(316);
    let source = directory.path().join("fixture.shared");
    let packed = directory.path().join("packed.shared");
    let compact = directory.path().join("compact.shared");
    crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
    crate::shared_pack::repack_shared_cores(&packed, &compact).unwrap();
    let reader = SharedReader::open(&compact).unwrap();
    let key = SharedKey::core(315);
    let group = reader.find(key).unwrap().unwrap();
    assert_eq!(group.core_ordinal(), 315);
    assert_eq!(315 * 13 % PAGE_BYTES as usize, PAGE_BYTES as usize - 1);
    assert_eq!(
        group_evidence(&reader, Some(group)),
        group_evidence(&baseline, baseline.find(key).unwrap())
    );
}

#[test]
fn packed_noninitial_groups_and_payloads_reject_corruption() {
    for group_corruption in [false, true] {
        let (directory, source, _) = fixture(16_384);
        drop(source);
        let original = directory.path().join("fixture.shared");
        let packed = directory.path().join("packed.shared");
        crate::shared_pack::repack_shared_index(&original, &packed).unwrap();
        mutate_and_resign(&packed, |bytes, header| {
            if group_corruption {
                let groups = header.section(Section::Groups);
                let at =
                    (groups.offset + groups.length - header.row_bytes(Section::Groups)) as usize;
                assert!(at > groups.offset as usize + 4096);
                bytes[at + 4] = 0xff;
            } else {
                let occurrences = header.section(Section::Occurrences);
                let at = (occurrences.offset + occurrences.length - 24) as usize;
                assert!(at > occurrences.offset as usize + 4096);
                bytes[at + 12] = 1;
            }
        });
        let reader = SharedReader::open(&packed).unwrap();
        let key = SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT,
            length: 31,
        };
        if group_corruption {
            assert!(reader.find_many(&[key]).is_err());
        } else {
            let group = reader.find_many(&[key]).unwrap()[0].unwrap();
            let member = reader.member(group, 2).unwrap().unwrap();
            assert!(reader.member_occurrences(group, member).is_err());
        }
    }
}

#[test]
fn packed_pilot_shaped_lists_preserve_every_synthetic_context() {
    let directory = tempfile::tempdir().unwrap();
    let metadata = directory.path().join("metadata.jidx");
    let mut writer = JidxWriter::new(
        &metadata,
        &JidxInput {
            k: 15,
            rescue_k15: false,
            minimizer_window: 64,
            jam_sha256: [41; 32],
            manifest_sha256: [42; 32],
        },
    )
    .unwrap();
    for id in 0..12 {
        writer
            .begin_metagenome(MetagenomeInput {
                name: format!("sample-{id:02}"),
                bgzf_uri: format!("sample-{id:02}.bgz"),
                bgzf_bytes: 100,
                bgzf_sha256: [43; 32],
                gzi: vec![0; 8],
            })
            .unwrap();
        writer
            .begin_contig(ContigInput {
                name: format!("contig-{id}"),
                length: 100_000,
                fasta_offset: 4,
                line_bases: 80,
                line_width: 81,
            })
            .unwrap();
    }
    writer.finish().unwrap();
    let reference = JidxReader::open(&metadata).unwrap();
    let mut seeds = Vec::new();
    let mut keys = Vec::new();
    // 21.7% repeated cores; pilot common bins, with its rare 3282-placement tail oversampled.
    for core in 0..5000u32 {
        let count = if core >= 1085 {
            1
        } else if core == 0 {
            3282
        } else {
            match core % 1000 {
                0..=699 => 2,
                700..=894 => 3,
                895..=991 => 4,
                _ => 8,
            }
        };
        for ordinal in 0..count {
            let member = (core + ordinal % 3) % 12;
            let mut seed = SharedSeed {
                core,
                context: (((core * 17 + ordinal * 101) & 4095) << 20)
                    | ((core * 31 + ordinal * 7) & 0xfffff),
                flags: 6 | (ordinal & 1) as u8,
                position: 100 + u64::from(ordinal),
            };
            if (core + ordinal) % 250 == 0 {
                seed.context = 0;
                seed.flags &= 1;
            } else if (core + ordinal) % 100 == 0 {
                seed.context &= !0xfffff;
                seed.flags &= 3;
            }
            for length in [15, 21, 31] {
                if let Some(key) = seed.key(length) {
                    keys.push(key);
                }
            }
            seeds.push(IndexedSeed {
                member,
                contig: member,
                seed,
            });
        }
    }
    let source = directory.path().join("source.shared");
    let packed = directory.path().join("packed.shared");
    let original = write_shared_index(&reference, &source, 64, &mut seeds).unwrap();
    let stats = crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
    assert_eq!(original.singleton_cores, 3915);
    assert_eq!(stats.build.core_count, 5000);
    assert_eq!(original.occurrence_references, stats.logical_references);
    assert_eq!(original.member_descriptors, stats.logical_memberships);
    assert!(stats.inline_member_groups > 0 && stats.direct_placement_members > 0);
    assert!(stats.build.index_bytes < original.index_bytes);
    keys.sort_unstable_by_key(|key| (key.core, key.length, key.context));
    keys.dedup();
    let baseline = SharedReader::open(&source).unwrap();
    let reader = SharedReader::open(&packed).unwrap();
    for chunk in keys.chunks(127) {
        let groups = reader.find_many(chunk).unwrap();
        for (&key, actual) in chunk.iter().zip(groups) {
            assert_eq!(
                group_evidence(&reader, actual),
                group_evidence(&baseline, baseline.find(key).unwrap())
            );
        }
    }
}

#[test]
fn packed_direct_locators_are_bounded_by_physical_placements() {
    for inline in [false, true] {
        let (directory, source) = grouped_lookup_fixture(0);
        let packed = directory.path().join("packed.shared");
        crate::shared_pack::repack_shared_index(&source, &packed).unwrap();
        mutate_and_resign(&packed, |bytes, header| {
            let width = header.id_bytes();
            let physical = header.section(Section::Occurrences).length / 24;
            assert!(header.section(Section::References).length / 4 > physical);
            let section = if inline {
                Section::Groups
            } else {
                Section::Members
            };
            let range = header.section(section);
            let mut changed = 0;
            for row in bytes[range.offset as usize..(range.offset + range.length) as usize]
                .chunks_exact_mut(header.row_bytes(section) as usize)
            {
                let (first, count) = if inline {
                    (5 + width, 9 + width)
                } else {
                    (width, width + 4)
                };
                if (!inline || row[4] & 4 != 0) && crate::jidx::read_u32(row, count) == 1 {
                    row[first..first + 4].copy_from_slice(&(physical as u32).to_le_bytes());
                    changed += 1;
                }
            }
            assert!(changed > 0);
        });
        let reader = SharedReader::open(&packed).unwrap();
        if inline {
            let key = SharedKey {
                core: TARGET_CORE,
                context: TARGET_CONTEXT ^ 1,
                length: 31,
            };
            assert!(reader.find(key).is_err());
        } else {
            let group = reader.find(SharedKey::core(TARGET_CORE)).unwrap().unwrap();
            assert!(reader.members(group).is_err());
        }
    }
}
