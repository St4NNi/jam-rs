use crate::jidx::sha256;
use crate::jidx_reader::JidxReader;
use crate::jidx_writer::{ContigInput, JidxInput, JidxWriter, MetagenomeInput};
use crate::owner_format::checksum_layout;
use crate::shared_format::{HEADER_BYTES, PAGE_BYTES, Section, SharedError, SharedHeader};
use crate::shared_reader::{SharedGroup, SharedReadStats, SharedReader};
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
fn exact_counts_and_absence_do_not_decode_payloads_as_prefix_grows() {
    for preceding in [16u32, 4096, 16_384] {
        let (_directory, reader, build) = fixture(preceding);
        assert!(reader.stats().observed);
        assert_eq!(build.repeated_core_positions, u64::from(preceding) * 2 + 3);
        assert_eq!(build.occurrence_references, u64::from(preceding) * 2 + 9);
        let before = reader.stats();
        let group = reader.find(SharedKey::core(TARGET_CORE)).unwrap().unwrap();
        assert_eq!((group.member_count(), group.occurrence_count()), (2, 3));
        let work = delta(reader.stats(), before);
        let logarithmic_bound = u64::from(u32::BITS - (preceding + 1).leading_zeros()) + 1;
        assert!(work[0] <= logarithmic_bound);
        assert!(work[1] <= 3);
        assert_eq!(&work[2..], &[0, 0, 0]);

        let before = reader.stats();
        assert!(
            reader
                .find(SharedKey::core(CORE_LIMIT - 1))
                .unwrap()
                .is_none()
        );
        let work = delta(reader.stats(), before);
        assert!(work[0] <= logarithmic_bound);
        assert_eq!(&work[1..], &[0, 0, 0, 0]);

        let before = reader.stats();
        let context = SharedKey {
            core: TARGET_CORE,
            context: TARGET_CONTEXT >> 20,
            length: 21,
        };
        let group = reader.find(context).unwrap().unwrap();
        assert_eq!((group.member_count(), group.occurrence_count()), (2, 3));
        let work = delta(reader.stats(), before);
        assert!(work[0] <= logarithmic_bound);
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
