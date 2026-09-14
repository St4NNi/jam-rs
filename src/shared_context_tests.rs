use super::*;

#[test]
fn retained_context_members_match_preserved_counts_and_positions() {
    for wide in [false, true] {
        let (_directory, preserved, placed, keys, _) = context_fixture(wide);
        let reference = SharedReader::open(preserved).unwrap();
        let reader = SharedReader::open_observed(&placed).unwrap();
        for workers in [1, 4, 8, 16] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            pool.install(|| {
                // Also request children without their 21 ancestor, with repeated associations.
                for children_only in [false, true] {
                    for same_core in keys.chunk_by(|left, right| left.core == right.core) {
                        let Some(core) = reader.find(SharedKey::core(same_core[0].core)).unwrap()
                        else {
                            continue;
                        };
                        let contexts = same_core
                            .iter()
                            .copied()
                            .filter(|key| !children_only || key.length == 31)
                            .flat_map(|key| [key, key])
                            .collect::<Vec<_>>();
                        let expected = reader.find_in_core(core, &contexts).unwrap();
                        let mut output = vec![None; contexts.len()];
                        let mut retained = Vec::with_capacity(2048);
                        let operation = reader.posting_operation().unwrap();
                        let complete = operation
                            .find_in_core_with_members_into(
                                core,
                                &contexts,
                                &mut output,
                                &mut retained,
                            )
                            .unwrap();
                        assert_eq!(output.as_slice(), expected.as_slice());
                        operation.finish().unwrap();
                        for (ordinal, (&key, group)) in contexts.iter().zip(output).enumerate() {
                            assert_eq!(
                                placed_evidence(&reader, group),
                                placed_evidence(&reference, reference.find(key).unwrap())
                            );
                            let Some(group) = group else {
                                continue;
                            };
                            if !complete {
                                continue;
                            }
                            let first = if ordinal % 2 == 0 {
                                ordinal
                            } else {
                                ordinal - 1
                            };
                            let actual = retained
                                .iter()
                                .filter(|(index, _)| *index == first)
                                .map(|(_, member)| *member)
                                .collect::<Vec<_>>();
                            assert_eq!(actual, reader.members(group).unwrap());
                            let before = reader.stats().placement_bound_searches;
                            let operation = reader.posting_operation().unwrap();
                            for member in actual {
                                let mut positions = vec![
                                    crate::jidx_reader::SeedOccurrence {
                                        contig_id: 0,
                                        position: 0,
                                        canonical_orientation: false,
                                    };
                                    member.occurrence_count() as usize
                                ];
                                operation
                                    .fill_occurrence_block(group, member, 0, &mut positions)
                                    .unwrap();
                                assert_eq!(
                                    positions,
                                    reader.member_occurrences(group, member).unwrap()
                                );
                            }
                            operation.finish().unwrap();
                            assert_eq!(
                                reader.stats().placement_bound_searches,
                                before,
                                "retained ranges must not be searched during position fill"
                            );
                        }
                    }
                }
            });
        }
    }
}

#[test]
fn retained_context_admission_is_complete_and_directory_is_shared() {
    let (_directory, _, placed, _, _) = context_fixture(false);
    let reader = SharedReader::open_observed(&placed).unwrap();
    let core = reader.find(SharedKey::core(TARGET_CORE)).unwrap().unwrap();
    let operation = reader.posting_operation().unwrap();
    let before_empty = reader.stats();
    assert!(
        !operation
            .find_in_core_with_members_into(core, &[], &mut [], &mut Vec::new())
            .unwrap()
    );
    assert_eq!(
        reader.stats().core_descriptor_inspections,
        before_empty.core_descriptor_inspections
    );
    assert_eq!(
        reader.stats().member_descriptor_inspections,
        before_empty.member_descriptor_inspections
    );
    operation.finish().unwrap();
    let contexts = [
        SharedKey {
            core: TARGET_CORE,
            context: 0x5a5,
            length: 21,
        },
        SharedKey {
            core: TARGET_CORE,
            context: (0x5a5 << 20) | 1,
            length: 31,
        },
        SharedKey {
            core: TARGET_CORE,
            context: (0x5a5 << 20) | 2,
            length: 31,
        },
    ];
    let before = reader.stats();
    let expected = reader.find_in_core(core, &contexts).unwrap();
    let baseline =
        reader.stats().member_descriptor_inspections - before.member_descriptor_inspections;
    assert_eq!(expected[0].unwrap().member_count(), 2);
    assert_eq!(expected[0].unwrap().occurrence_count(), 4);
    assert_eq!(
        expected[1].unwrap().member_count() + expected[2].unwrap().member_count(),
        3
    );
    for capacity in [0, 1, 3, 5] {
        let operation = reader.posting_operation().unwrap();
        let mut members = Vec::with_capacity(capacity);
        let mut output = [None; 3];
        let before = reader.stats();
        let complete = operation
            .find_in_core_with_members_into(core, &contexts, &mut output, &mut members)
            .unwrap();
        operation.finish().unwrap();
        assert_eq!(output.as_slice(), expected.as_slice());
        assert_eq!(complete, capacity >= 5);
        assert_eq!(members.len(), if complete { 5 } else { 0 });
        let after = reader.stats();
        assert_eq!(after.context_member_runs - before.context_member_runs, 4);
        assert_eq!(
            after.member_descriptor_inspections - before.member_descriptor_inspections,
            4
        );
        assert_eq!(baseline, 12);
        assert_eq!(
            after.context_member_fallbacks - before.context_member_fallbacks,
            u64::from(!complete)
        );
    }
}

#[test]
fn retained_context_failures_clear_only_unpublished_members() {
    let (directory, _, placed, _, _) = context_fixture(false);
    let reader = SharedReader::open(&placed).unwrap();
    let other = SharedReader::open(&placed).unwrap();
    let core = reader.find(SharedKey::core(30)).unwrap().unwrap();
    let foreign = other.find(core.key()).unwrap().unwrap();
    let key = SharedKey {
        core: 30,
        context: 0x5a5,
        length: 21,
    };
    let operation = reader.posting_operation().unwrap();
    let mut members = Vec::with_capacity(16);
    let mut output = [None];
    assert!(
        operation
            .find_in_core_with_members_into(core, &[key], &mut output, &mut members)
            .unwrap()
    );
    let prefix = members.clone();
    assert!(matches!(
        operation.find_in_core_with_members_into(foreign, &[key], &mut output, &mut members),
        Err(SharedError::Invalid("group handle"))
    ));
    assert_eq!(output, [None]);
    assert_eq!(members, prefix);
    assert!(matches!(
        operation.find_in_core_with_members_into(core, &[], &mut output, &mut members),
        Err(SharedError::Invalid("context result storage"))
    ));
    assert_eq!(members, prefix);
    operation.finish().unwrap();
    let path = directory.path().join("retained-corrupt.shared");
    std::fs::copy(&placed, &path).unwrap();
    mutate_and_resign(&path, |bytes, header| {
        let at = header.section(Section::Members).offset as usize;
        let width = header.row_bytes(Section::Members) as usize;
        bytes.swap(at, at + width);
    });
    let corrupt = SharedReader::open(path).unwrap();
    let core = corrupt.find(core.key()).unwrap().unwrap();
    let operation = corrupt.posting_operation().unwrap();
    members.clear();
    assert!(matches!(
        operation.find_in_core_with_members_into(core, &[key], &mut output, &mut members),
        Err(SharedError::Invalid("member directory"))
    ));
    assert_eq!(output, [None]);
    assert!(members.is_empty());
}

#[test]
fn retained_context_batch_postings_do_not_search_members_again() {
    use crate::trace_batch::{prepare_cores, prepare_lookup_with_cores};
    use crate::trace_index::TraceIndex;
    let (_directory, _, placed, _, _) = context_fixture(false);
    for workers in [1, 4, 8, 16] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        pool.install(|| {
            let index = TraceIndex::Shared(Box::new(SharedReader::open_observed(&placed).unwrap()));
            let cores = prepare_cores(&index, [TARGET_CORE], 1, true)
                .unwrap()
                .unwrap();
            let contexts = [
                (21, 0x5a5),
                (31, (0x5a5 << 20) | 1),
                (31, (0x5a5 << 20) | 2),
            ]
            .map(|(length, context)| {
                SharedKey {
                    core: TARGET_CORE,
                    context,
                    length,
                }
                .packed()
                .unwrap()
            });
            let lookup = prepare_lookup_with_cores(
                &index,
                contexts
                    .into_iter()
                    .flat_map(|key| [(key, 0), (key, 1)])
                    .collect(),
                2,
                true,
                Some(&cores),
            )
            .unwrap()
            .unwrap();
            assert!(lookup.postings_complete);
            assert_eq!(lookup.entries.len(), 3);
            assert_eq!(lookup.query_entries.len(), 6);
            let TraceIndex::Shared(reader) = &index else {
                unreachable!()
            };
            let stats = reader.stats();
            assert_eq!(stats.retained_context_members, 5);
            assert_eq!(stats.context_posting_members_reused, 5);
            assert_eq!(stats.context_posting_members_fallback, 0);
            assert_eq!(
                stats.placement_bound_searches, 24,
                "three contexts times four members times two bounds, with no posting re-search"
            );
            assert_eq!(
                lookup
                    .postings
                    .iter()
                    .flatten()
                    .map(|posting| posting.documents.len())
                    .sum::<usize>(),
                5
            );
        });
    }
}

#[test]
fn checked_contexts_match_public_results_with_constant_identity_checks() {
    let (_directory, reader, _) = fixture(2);
    let cores = [0, TARGET_CORE, TARGET_CORE + 1]
        .map(|core| reader.find(SharedKey::core(core)).unwrap().unwrap());
    let keys = cores.map(|core| {
        [
            core.key(),
            SharedKey {
                core: core.key().core,
                context: TARGET_CONTEXT,
                length: 31,
            },
            SharedKey {
                core: core.key().core,
                context: TARGET_CONTEXT + 1,
                length: 31,
            },
        ]
    });
    let before = reader.stats().file.identity_checks;
    let expected = cores
        .iter()
        .zip(&keys)
        .map(|(&core, keys)| reader.find_in_core(core, keys).unwrap())
        .collect::<Vec<_>>();
    assert_eq!(reader.stats().file.identity_checks - before, 6);
    let before = reader.stats().file.identity_checks;
    let operation = reader.posting_operation().unwrap();
    let mut scratch = [None; 3];
    let mut actual = Vec::new();
    for (&core, keys) in cores.iter().zip(&keys) {
        operation
            .find_in_core_into(core, keys, &mut scratch)
            .unwrap();
        actual.push(scratch.to_vec());
    }
    operation.finish().unwrap();
    assert_eq!(reader.stats().file.identity_checks - before, 2);
    assert_eq!(actual, expected);
    assert!(actual[1][1].is_some());
    assert!(actual[2][1].is_none());
    for (actual, expected) in actual
        .into_iter()
        .flatten()
        .zip(expected.into_iter().flatten())
    {
        assert_eq!(
            group_evidence(&reader, actual),
            group_evidence(&reader, expected)
        );
    }
}

#[test]
fn checked_context_errors_clear_output_and_follow_input_order() {
    let (directory, reader, _) = fixture(0);
    let core = reader.find(SharedKey::core(TARGET_CORE)).unwrap().unwrap();
    let other = SharedReader::open(directory.path().join("fixture.shared")).unwrap();
    let foreign = other.find(core.key()).unwrap().unwrap();
    let operation = reader.posting_operation().unwrap();
    let malformed = SharedKey {
        length: 14,
        ..core.key()
    };
    let wrong_core = SharedKey::core(TARGET_CORE + 1);
    let long = SharedKey {
        context: TARGET_CONTEXT,
        length: 31,
        ..core.key()
    };
    for (keys, message) in [
        ([malformed, wrong_core], "shared key"),
        ([wrong_core, malformed], "core group key"),
        ([long, core.key()], "sorted core contexts"),
    ] {
        let mut output = [Some(core); 2];
        assert!(
            matches!(operation.find_in_core_into(core, &keys, &mut output),
            Err(SharedError::Invalid(actual)) if actual == message)
        );
        assert_eq!(output, [None; 2]);
    }
    let mut output = [Some(core)];
    assert!(matches!(
        operation.find_in_core_into(foreign, &[core.key()], &mut output),
        Err(SharedError::Invalid("group handle"))
    ));
    assert_eq!(output, [None]);
    assert!(matches!(
        operation.find_in_core_into(foreign, &[], &mut []),
        Err(SharedError::Invalid("group handle"))
    ));
    assert!(matches!(
        operation.find_in_core_into(core, &[], &mut output),
        Err(SharedError::Invalid("context result storage"))
    ));
    assert_eq!(output, [None]);
    operation.finish().unwrap();
}

#[test]
fn checked_context_finish_rejects_source_mutation() {
    let (directory, reader, _) = fixture(0);
    let core = reader.find(SharedKey::core(TARGET_CORE)).unwrap().unwrap();
    let operation = reader.posting_operation().unwrap();
    let mut output = [None];
    operation
        .find_in_core_into(core, &[core.key()], &mut output)
        .unwrap();
    assert_eq!(output, [Some(core)]);
    mutate_and_resign(&directory.path().join("fixture.shared"), |bytes, header| {
        let at = header.section(Section::Occurrences).offset as usize + 8;
        let flags = crate::jidx::read_u32(bytes, at) ^ 1;
        bytes[at..at + 4].copy_from_slice(&flags.to_le_bytes());
    });
    assert!(matches!(
        operation.finish(),
        Err(SharedError::SourceChanged)
    ));
}

#[test]
fn checked_context_adapter_preserves_ordinals_and_publishes_only_complete_results() {
    use crate::trace::TraceError;
    use crate::trace_batch::prepare_cores;
    use crate::trace_index::TraceIndex;

    let (directory, reader, _) = fixture(2);
    let other = SharedReader::open(directory.path().join("fixture.shared")).unwrap();
    let foreign = other
        .find(SharedKey::core(TARGET_CORE + 1))
        .unwrap()
        .unwrap();
    let index = TraceIndex::Shared(Box::new(reader));
    let mut cores = prepare_cores(&index, [0, TARGET_CORE, TARGET_CORE + 1], 3, false)
        .unwrap()
        .unwrap();
    let long = SharedKey {
        core: TARGET_CORE,
        context: TARGET_CONTEXT,
        length: 31,
    }
    .packed()
    .unwrap();
    let keys = [
        long,
        u64::from(TARGET_CORE + 1),
        0,
        u64::from(TARGET_CORE - 1),
        u64::from(TARGET_CORE),
        long,
    ];
    let expected = index.find_seeds_batch(&keys).unwrap();
    let TraceIndex::Shared(reader) = &index else {
        unreachable!()
    };
    let before = reader.stats().file.identity_checks;
    assert_eq!(index.find_seeds_in_cores(&keys, &cores).unwrap(), expected);
    assert_eq!(reader.stats().file.identity_checks - before, 3);
    std::sync::Arc::get_mut(&mut cores).unwrap().groups[2] = foreign;
    for malformed in [u64::MAX, 1 << 30, (1 << 62) | (1 << 42)] {
        let before = reader.stats().core_descriptor_inspections;
        assert!(matches!(
            index.find_seeds_in_cores(&[long, malformed], &cores),
            Err(TraceError::Invalid("shared context key"))
        ));
        assert_eq!(reader.stats().core_descriptor_inspections, before);
    }
    assert!(matches!(
        index.find_seeds_in_cores(&[0, u64::from(TARGET_CORE + 1)], &cores),
        Err(TraceError::Shared(SharedError::Invalid("group handle")))
    ));
}

#[test]
fn checked_context_result_admission_clears_unpublished_records() {
    let (_directory, reader, _) = fixture(0);
    let core = reader.find(SharedKey::core(TARGET_CORE)).unwrap().unwrap();
    let count = 64 * 1024 * 1024 / std::mem::size_of::<Option<SharedGroup>>() + 1;
    let keys = vec![core.key(); count];
    let mut output = vec![Some(core); count];
    let operation = reader.posting_operation().unwrap();
    let before = reader.stats().core_descriptor_inspections;
    assert!(matches!(
        operation.find_in_core_into(core, &keys, &mut output),
        Err(SharedError::ResourceLimit)
    ));
    assert!(output.iter().all(Option::is_none));
    assert_eq!(reader.stats().core_descriptor_inspections, before);
    operation.finish().unwrap();
}
