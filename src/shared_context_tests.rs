use super::*;

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
