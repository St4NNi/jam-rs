use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use jam_rs::alignment::{
    Alignment, AlignmentConfig, AlignmentWorkspace, EditOperation, EditRun, Interval, Strand,
};
use jam_rs::jidx::sha256;
use jam_rs::jidx_writer::{ContigInput, JidxInput, JidxWriter, MetagenomeInput};
use jam_rs::owner_postings;
use jam_rs::shared_pack::{add_shared_core_filter, repack_shared_cores, repack_shared_index};
use jam_rs::shared_reader::SharedReader;
use jam_rs::shared_seed::{HAS_CONTEXT_21, HAS_CONTEXT_31, SharedKey, select_shared_seeds};
use jam_rs::shared_writer::build_shared_index;
use jam_rs::trace::{SearchCompletion, TraceConfig, TraceEngine};
use noodles_bgzf::{self as bgzf, gzi};
use std::fs::File;
use std::hint::black_box;
use std::io::Write;
use std::path::PathBuf;
use std::sync::OnceLock;
use std::time::Duration;
#[cfg(feature = "bench-internals")]
use xorf::{BinaryFuse8Ref, Filter, FilterRef};

fn sequence(length: usize) -> Vec<u8> {
    sequence_from_state(length, 7)
}

fn sequence_from_state(length: usize, mut state: u64) -> Vec<u8> {
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(state >> 62) as usize]
        })
        .collect()
}

fn affine_alignment(criterion: &mut Criterion) {
    let query = sequence(5_000);
    let mut target = query.clone();
    for position in (50..target.len()).step_by(100) {
        target[position] = match target[position] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
    }
    let config = AlignmentConfig::default();
    let mut workspace = AlignmentWorkspace::default();
    assert!(workspace.align(&query, &target, config).unwrap().identity() > 0.98);

    let mut group = criterion.benchmark_group("trace_alignment");
    group.throughput(Throughput::Bytes(query.len() as u64));
    group.bench_function("5kb_99pct", |bencher| {
        bencher.iter(|| {
            workspace
                .align(black_box(&query), black_box(&target), config)
                .unwrap()
        })
    });
    group.finish();
}

fn endpoint_fixture(flank: usize) -> (Vec<u8>, Vec<u8>, Alignment) {
    let core_length = 32;
    let mut query = vec![b'A'; flank];
    query.extend(std::iter::repeat_n(b'C', core_length));
    query.extend(std::iter::repeat_n(b'G', flank));
    let mut target = vec![b'T'; flank];
    target.extend(std::iter::repeat_n(b'C', core_length));
    target.extend(std::iter::repeat_n(b'A', flank));
    let core = Alignment {
        score: (core_length * 2) as i32,
        strand: Strand::Forward,
        query_interval: Interval::new(flank as u64, (flank + core_length) as u64).unwrap(),
        target_interval: Interval::new(flank as u64, (flank + core_length) as u64).unwrap(),
        matches: core_length as u64,
        substitutions: 0,
        insertions: 0,
        deletions: 0,
        cigar: format!("{core_length}="),
        edit_script: vec![EditRun {
            operation: EditOperation::Equal,
            length: core_length as u32,
        }],
    };
    (query, target, core)
}

fn trace_endpoint(criterion: &mut Criterion) {
    let config = AlignmentConfig::default();
    let mut group = criterion.benchmark_group("trace_endpoint");
    group.sample_size(10);
    group.warm_up_time(Duration::from_millis(100));
    group.measurement_time(Duration::from_millis(300));
    for flank in [32, 256] {
        let (query, target, core) = endpoint_fixture(flank);
        group.throughput(Throughput::Bytes((query.len() + target.len()) as u64));
        group.bench_function(format!("{flank}/fresh"), |bencher| {
            bencher.iter(|| {
                AlignmentWorkspace::default()
                    .complete_endpoints(
                        black_box(core.clone()),
                        black_box(&query),
                        black_box(&target),
                        0,
                        flank,
                        config,
                    )
                    .unwrap()
            })
        });
        let mut workspace = AlignmentWorkspace::default();
        group.bench_function(format!("{flank}/reused"), |bencher| {
            bencher.iter(|| {
                workspace
                    .complete_endpoints(
                        black_box(core.clone()),
                        black_box(&query),
                        black_box(&target),
                        0,
                        flank,
                        config,
                    )
                    .unwrap()
            })
        });
    }
    group.finish();
}

fn owner_postings(criterion: &mut Criterion) {
    use owner_postings::{OwnerKey, OwnerMember, OwnerOccurrence};

    let mut key = 0u64;
    let keys = (0..256)
        .map(|ordinal| {
            key += 1 + (ordinal * 65_537) % 1_000_003;
            OwnerKey {
                key,
                members: (0..if ordinal % 64 == 0 { 128 } else { 1 })
                    .map(|document_id| OwnerMember {
                        document_id,
                        occurrences: (0..if ordinal % 64 == 0 { 16 } else { 1 })
                            .map(|position| OwnerOccurrence {
                                local_contig: position / 4,
                                position: u64::from(position % 4) * 17 + 31,
                                canonical_orientation: position % 2 == 0,
                            })
                            .collect(),
                    })
                    .collect(),
            }
        })
        .collect::<Vec<_>>();
    let widths = [32; 128];
    let encoded = owner_postings::encode_block(&keys, &widths, |_, occurrence| {
        Ok(u64::from(occurrence.local_contig) * 1024 + occurrence.position)
    })
    .unwrap();
    assert_eq!(
        owner_postings::decode_block(
            &encoded.hot,
            &encoded.cold,
            &widths,
            64 * 1024 * 1024,
            |_, locus, canonical_orientation| {
                Ok(OwnerOccurrence {
                    local_contig: (locus / 1024) as u32,
                    position: locus % 1024,
                    canonical_orientation,
                })
            },
        )
        .unwrap(),
        keys
    );
    let hot = owner_postings::parse_hot(&encoded.hot).unwrap();
    let selected = owner_postings::find_key(&hot, keys[0].key).unwrap();
    let member = owner_postings::key_members(&hot, selected).unwrap()[0];
    let window = owner_postings::locate_member(&hot, member.member_ordinal).unwrap();
    let start = window.start_bit / 8;
    let end = window.end_bit.div_ceil(8);
    let tiny_key = owner_postings::find_key(&hot, keys[1].key).unwrap();
    let tiny_member = owner_postings::key_members(&hot, tiny_key).unwrap()[0];
    let tiny_window = owner_postings::locate_member(&hot, tiny_member.member_ordinal).unwrap();
    let tiny_start = tiny_window.start_bit / 8;
    let tiny_end = tiny_window.end_bit.div_ceil(8);
    let mut group = criterion.benchmark_group("owner_postings");
    group.sample_size(20);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    group.throughput(Throughput::Bytes(
        (encoded.hot.len() + encoded.cold.len()) as u64,
    ));
    group.bench_function("encode_mixed", |bencher| {
        bencher.iter(|| {
            owner_postings::encode_block(black_box(&keys), black_box(&widths), |_, occurrence| {
                Ok(u64::from(occurrence.local_contig) * 1024 + occurrence.position)
            })
            .unwrap()
        })
    });
    group.bench_function("decode_mixed", |bencher| {
        bencher.iter(|| {
            owner_postings::decode_block(
                black_box(&encoded.hot),
                black_box(&encoded.cold),
                black_box(&widths),
                64 * 1024 * 1024,
                |_, locus, canonical_orientation| {
                    Ok(OwnerOccurrence {
                        local_contig: (locus / 1024) as u32,
                        position: locus % 1024,
                        canonical_orientation,
                    })
                },
            )
            .unwrap()
        })
    });
    group.throughput(Throughput::Bytes(encoded.hot.len() as u64));
    group.bench_function("hot_directory", |bencher| {
        bencher.iter(|| owner_postings::parse_hot(black_box(&encoded.hot)).unwrap())
    });
    group.bench_function("absent_parse_and_find", |bencher| {
        bencher.iter(|| {
            let parsed = owner_postings::parse_hot(black_box(&encoded.hot)).unwrap();
            owner_postings::find_key(&parsed, black_box(key + 1))
                .map(|entry| entry.document_frequency)
        })
    });
    group.bench_function("absent_count", |bencher| {
        bencher.iter(|| owner_postings::find_key(black_box(&hot), black_box(key + 1)))
    });
    group.throughput(Throughput::Elements(member.occurrence_count));
    group.bench_function("selected_member", |bencher| {
        bencher.iter(|| {
            owner_postings::decode_member_window(
                black_box(&encoded.cold[start as usize..end as usize]),
                start,
                window,
                black_box(&hot),
                black_box(&widths),
                64 * 1024 * 1024,
            )
            .unwrap()
        })
    });
    group.throughput(Throughput::Elements(tiny_member.occurrence_count));
    group.bench_function("selected_tiny_member", |bencher| {
        bencher.iter(|| {
            owner_postings::decode_member_window(
                black_box(&encoded.cold[tiny_start as usize..tiny_end as usize]),
                tiny_start,
                tiny_window,
                black_box(&hot),
                black_box(&widths),
                64 * 1024 * 1024,
            )
            .unwrap()
        })
    });
    group.finish();
}

fn synthetic_gzi(bytes: &[u8]) -> gzi::Index {
    let mut entries = Vec::new();
    let mut compressed = 0usize;
    let mut uncompressed = 0u64;
    while compressed < bytes.len() {
        let header = &bytes[compressed..compressed + 18];
        let block_bytes = usize::from(u16::from_le_bytes([header[16], header[17]])) + 1;
        let end = compressed + block_bytes;
        let uncompressed_bytes = u32::from_le_bytes(bytes[end - 4..end].try_into().unwrap());
        if compressed != 0 && uncompressed_bytes != 0 {
            entries.push((compressed as u64, uncompressed));
        }
        uncompressed += u64::from(uncompressed_bytes);
        compressed = end;
    }
    gzi::Index::from(entries)
}

fn shared_lookup_fixture() -> (
    tempfile::TempDir,
    SharedReader,
    Vec<(&'static str, Vec<SharedKey>)>,
) {
    let directory = tempfile::tempdir().unwrap();
    let mut dna = sequence(80_000);
    let repeated = sequence(127);
    while dna.len() < 100_000 {
        dna.extend_from_slice(&repeated);
    }
    dna.truncate(100_000);
    let bgzf_path = directory.path().join("target.bgz");
    let mut fasta = b">target\n".to_vec();
    for line in dna.chunks(80) {
        fasta.extend_from_slice(line);
        fasta.push(b'\n');
    }
    let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
    writer.write_all(&fasta).unwrap();
    writer.finish().unwrap();
    let bgzf_path = std::fs::canonicalize(bgzf_path).unwrap();
    let bgzf_bytes = std::fs::read(&bgzf_path).unwrap();
    let gzi_path = directory.path().join("target.gzi");
    gzi::fs::write(&gzi_path, &synthetic_gzi(&bgzf_bytes)).unwrap();

    let metadata_path = directory.path().join("target.jidx");
    let mut metadata = JidxWriter::new(
        &metadata_path,
        &JidxInput {
            k: 15,
            rescue_k15: false,
            minimizer_window: 64,
            jam_sha256: [1; 32],
            manifest_sha256: [2; 32],
        },
    )
    .unwrap();
    metadata
        .begin_metagenome(MetagenomeInput {
            name: "target".to_owned(),
            bgzf_uri: bgzf_path.to_str().unwrap().to_owned(),
            bgzf_bytes: bgzf_bytes.len() as u64,
            bgzf_sha256: sha256(&bgzf_bytes),
            gzi: std::fs::read(gzi_path).unwrap(),
        })
        .unwrap();
    metadata
        .begin_contig(ContigInput {
            name: "target".to_owned(),
            length: dna.len() as u64,
            fasta_offset: 8,
            line_bases: 80,
            line_width: 81,
        })
        .unwrap();
    metadata.finish().unwrap();
    let shared_path = directory.path().join("target.shared");
    build_shared_index(&metadata_path, &shared_path, 64).unwrap();
    let reader = SharedReader::open(&shared_path).unwrap();

    let seeds = select_shared_seeds(&dna, 64).unwrap();
    let mut nested = seeds
        .iter()
        .filter(|seed| seed.flags & (HAS_CONTEXT_21 | HAS_CONTEXT_31) == 6)
        .take(128)
        .flat_map(|seed| {
            [
                seed.key(15).unwrap(),
                seed.key(21).unwrap(),
                seed.key(31).unwrap(),
            ]
        })
        .collect::<Vec<_>>();
    nested.sort_unstable();
    nested.dedup();
    let mut independent = seeds
        .iter()
        .map(|seed| SharedKey::core(seed.core))
        .collect::<Vec<_>>();
    independent.sort_unstable();
    independent.dedup();
    independent.truncate(512);
    let common = seeds
        .iter()
        .filter(|seed| seed.position >= 80_000 && seed.flags & 6 == 6)
        .max_by_key(|seed| {
            reader
                .find(SharedKey::core(seed.core))
                .unwrap()
                .unwrap()
                .occurrence_count()
        })
        .unwrap();
    let common = [
        common.key(15).unwrap(),
        common.key(21).unwrap(),
        common.key(31).unwrap(),
    ];
    let high_multiplicity = common.into_iter().cycle().take(1_536).collect();
    let mut absent = Vec::new();
    for ordinal in 0u32.. {
        let core = ordinal.wrapping_mul(506_952_113) & ((1 << 30) - 1);
        if reader.find(SharedKey::core(core)).unwrap().is_none() {
            absent.push(SharedKey::core(core));
            if absent.len() == 512 {
                break;
            }
        }
    }
    (
        directory,
        reader,
        vec![
            ("positive_nested", nested),
            ("absent_heavy_spread", absent),
            ("low_reuse", independent),
            ("high_multiplicity", high_multiplicity),
        ],
    )
}

fn shared_lookup(criterion: &mut Criterion) {
    let (_directory, reader, workloads) = shared_lookup_fixture();
    let mut group = criterion.benchmark_group("shared_lookup");
    group.sample_size(20);
    group.nresamples(1_000);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    for (name, keys) in workloads {
        let scalar = keys
            .iter()
            .map(|&key| reader.find(key).unwrap())
            .collect::<Vec<_>>();
        assert_eq!(reader.find_many(&keys).unwrap(), scalar);
        group.throughput(Throughput::Elements(keys.len() as u64));
        group.bench_function(format!("{name}/scalar_calls"), |bencher| {
            bencher.iter(|| {
                black_box(&keys)
                    .iter()
                    .map(|&key| reader.find(key).unwrap())
                    .collect::<Vec<_>>()
            })
        });
        group.bench_function(format!("{name}/grouped_call"), |bencher| {
            bencher.iter(|| reader.find_many(black_box(&keys)).unwrap())
        });
    }
    group.finish();
}

fn absent_cores(present: &[u32], count: usize, prefix: Option<u32>) -> Vec<u32> {
    let mut absent = Vec::with_capacity(count);
    let mut ordinal = 0u32;
    while absent.len() < count {
        let core = match prefix {
            Some(prefix) => (prefix << 14) | (ordinal & ((1 << 14) - 1)),
            None => ordinal.wrapping_mul(506_952_113) & ((1 << 30) - 1),
        };
        if present.binary_search(&core).is_err() {
            absent.push(core);
        }
        ordinal = ordinal.wrapping_add(1);
    }
    absent
}

fn sampled_present_cores(present: &[u32], count: usize) -> Vec<u32> {
    if count == 0 {
        return Vec::new();
    }
    (0..count)
        .map(|ordinal| present[ordinal * (present.len() - 1) / count.saturating_sub(1).max(1)])
        .collect()
}

fn mixed_core_workload(
    present: &[u32],
    present_count: usize,
    prefix: Option<u32>,
) -> Vec<SharedKey> {
    const REQUESTS: usize = 4096;
    let mut cores = sampled_present_cores(present, present_count);
    cores.extend(absent_cores(present, REQUESTS - present_count, prefix));
    cores.sort_unstable();
    cores.dedup();
    assert_eq!(cores.len(), REQUESTS);
    cores.into_iter().map(SharedKey::core).collect()
}

fn prefix_core_workload(present: &[u32], prefix: u32) -> Vec<SharedKey> {
    const REQUESTS: usize = 4096;
    let mut cores = present
        .iter()
        .copied()
        .filter(|core| core >> 14 == prefix)
        .collect::<Vec<_>>();
    assert!(!cores.is_empty() && cores.len() < REQUESTS);
    cores.extend(absent_cores(present, REQUESTS - cores.len(), Some(prefix)));
    cores.sort_unstable();
    assert_eq!(cores.len(), REQUESTS);
    cores.into_iter().map(SharedKey::core).collect()
}

type SharedCoreReader = (&'static str, PathBuf, SharedReader);

struct SharedCoreFixture {
    directory: tempfile::TempDir,
    readers: Vec<SharedCoreReader>,
    workloads: Vec<(&'static str, Vec<SharedKey>)>,
}

fn shared_core_absent_fixture() -> &'static SharedCoreFixture {
    static FIXTURE: OnceLock<SharedCoreFixture> = OnceLock::new();
    FIXTURE.get_or_init(|| {
    let directory = tempfile::tempdir().unwrap();
    let dna = sequence(1_050_000);
    let bgzf_path = directory.path().join("large-target.bgz");
    let mut fasta = b">target\n".to_vec();
    for line in dna.chunks(80) {
        fasta.extend_from_slice(line);
        fasta.push(b'\n');
    }
    let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
    writer.write_all(&fasta).unwrap();
    writer.finish().unwrap();
    let bgzf_path = std::fs::canonicalize(bgzf_path).unwrap();
    let bgzf_bytes = std::fs::read(&bgzf_path).unwrap();
    let gzi_path = directory.path().join("large-target.gzi");
    gzi::fs::write(&gzi_path, &synthetic_gzi(&bgzf_bytes)).unwrap();

    let metadata_path = directory.path().join("large-target.jidx");
    let mut metadata = JidxWriter::new(
        &metadata_path,
        &JidxInput {
            k: 15,
            rescue_k15: false,
            minimizer_window: 1,
            jam_sha256: [3; 32],
            manifest_sha256: [4; 32],
        },
    )
    .unwrap();
    metadata
        .begin_metagenome(MetagenomeInput {
            name: "large-target".to_owned(),
            bgzf_uri: bgzf_path.to_str().unwrap().to_owned(),
            bgzf_bytes: bgzf_bytes.len() as u64,
            bgzf_sha256: sha256(&bgzf_bytes),
            gzi: std::fs::read(gzi_path).unwrap(),
        })
        .unwrap();
    metadata
        .begin_contig(ContigInput {
            name: "large-target".to_owned(),
            length: dna.len() as u64,
            fasta_offset: 8,
            line_bases: 80,
            line_width: 81,
        })
        .unwrap();
    metadata.finish().unwrap();
    let v1_path = directory.path().join("large-target-v1.shared");
    let build = build_shared_index(&metadata_path, &v1_path, 1).unwrap();
    assert!(build.core_count >= 1_000_000, "{} cores", build.core_count);
    let v2_path = directory.path().join("large-target-v2.shared");
    let v3_path = directory.path().join("large-target-v3.shared");
    let v4_path = directory.path().join("large-target-v4.shared");
    let converted = std::time::Instant::now();
    let v2 = repack_shared_index(&v1_path, &v2_path).unwrap();
    let v3 = repack_shared_cores(&v2_path, &v3_path).unwrap();
    let v4 = add_shared_core_filter(&v3_path, &v4_path, usize::MAX).unwrap();
    eprintln!(
        "shared_core_absent kernel setup: conversion {:.3}s, v2 {} bytes, v3 {} bytes, v4 {} bytes, hot {}, cold {}, prefixes {} bytes",
        converted.elapsed().as_secs_f64(),
        v2.build.index_bytes,
        v3.build.index_bytes,
        v4.index_bytes,
        v3.hot_core_bytes,
        v3.cold_core_bytes,
        v3.core_prefix_bytes,
    );

    let mut present = select_shared_seeds(&dna, 1)
        .unwrap()
        .into_iter()
        .map(|seed| seed.core)
        .collect::<Vec<_>>();
    present.sort_unstable();
    present.dedup();
    assert_eq!(present.len() as u64, build.core_count);
    let skew_prefix = present[present.len() / 2] >> 14;
    let workloads = vec![
        ("spread/absent", mixed_core_workload(&present, 0, None)),
        (
            "spread/10pct_present",
            mixed_core_workload(&present, 410, None),
        ),
        (
            "spread/50pct_present",
            mixed_core_workload(&present, 2048, None),
        ),
        ("spread/present", mixed_core_workload(&present, 4096, None)),
        (
            "one_prefix/available_present",
            prefix_core_workload(&present, skew_prefix),
        ),
    ];
    SharedCoreFixture {
        directory,
        readers: vec![
            ("v2", v2_path.clone(), SharedReader::open(v2_path).unwrap()),
            ("v3", v3_path.clone(), SharedReader::open(v3_path).unwrap()),
            ("v4", v4_path.clone(), SharedReader::open(v4_path).unwrap()),
        ],
        workloads,
    }
    })
}

fn shared_core_absent(criterion: &mut Criterion) {
    let fixture_started = std::time::Instant::now();
    let SharedCoreFixture {
        directory: _directory,
        readers,
        workloads,
    } = shared_core_absent_fixture();
    eprintln!(
        "shared_core_absent total kernel-only fixture setup: {:.3}s (excluded from lookup timing)",
        fixture_started.elapsed().as_secs_f64()
    );
    for (name, keys) in workloads {
        let mut expected = None;
        for (version, path, _) in readers {
            let observed = SharedReader::open_observed(path).unwrap();
            let actual = observed
                .find_many(keys)
                .unwrap()
                .into_iter()
                .map(|group| {
                    group.map(|group| (group.key(), group.member_count(), group.occurrence_count()))
                })
                .collect::<Vec<_>>();
            assert_eq!(
                *expected.get_or_insert(actual.clone()),
                actual,
                "{version}/{name}"
            );
            eprintln!(
                "shared_core_absent {version}/{name} observed: {:?}",
                observed.stats()
            );
        }
    }
    let mut group = criterion.benchmark_group("shared_core_absent");
    group.sample_size(20);
    group.nresamples(1_000);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    for (version, _, reader) in readers {
        for (name, keys) in workloads {
            group.throughput(Throughput::Elements(keys.len() as u64));
            group.bench_function(format!("{version}/{name}"), |bencher| {
                bencher.iter(|| reader.find_many(black_box(keys)).unwrap())
            });
        }
    }
    group.finish();
}

#[cfg(feature = "bench-internals")]
fn shared_core_search(criterion: &mut Criterion) {
    let fixture = shared_core_absent_fixture();
    let (_, path, reader) = fixture
        .readers
        .iter()
        .find(|(version, _, _)| *version == "v3")
        .unwrap();
    let mut group = criterion.benchmark_group("shared_core_search");
    group.sample_size(20);
    group.nresamples(1_000);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    for (name, keys) in &fixture.workloads {
        let cores = keys.iter().map(|key| key.core).collect::<Vec<_>>();
        let dense = reader
            .find_many(keys)
            .unwrap()
            .into_iter()
            .flatten()
            .map(|group| (group.key(), group.member_count(), group.occurrence_count()))
            .collect::<Vec<_>>();
        let sparse = reader
            .benchmark_resolve_sorted_cores(&cores)
            .unwrap()
            .into_iter()
            .map(|group| (group.key(), group.member_count(), group.occurrence_count()))
            .collect::<Vec<_>>();
        assert_eq!(sparse, dense, "{name}");

        let fresh = SharedReader::open_observed(path).unwrap();
        let fresh_results = fresh.benchmark_resolve_sorted_cores(&cores).unwrap();
        black_box(fresh_results);
        let stats = fresh.stats();
        eprintln!(
            "shared_core_search fresh-reader {name}: views={} comparisons={} payloads={} requested_bytes={} requested_pages={} authenticated_pages={}",
            stats.core_view_creations,
            stats.core_view_comparisons,
            stats.core_descriptor_inspections,
            stats.file.requested_bytes,
            stats.file.requested_pages,
            stats.file.authenticated_pages,
        );

        group.throughput(Throughput::Elements(keys.len() as u64));
        group.bench_function(format!("{name}/dense_find_many"), |bencher| {
            bencher.iter(|| reader.find_many(black_box(keys)).unwrap())
        });
        group.bench_function(format!("{name}/sparse_core_output"), |bencher| {
            bencher.iter(|| {
                reader
                    .benchmark_resolve_sorted_cores(black_box(&cores))
                    .unwrap()
            })
        });
    }
    group.finish();
}

#[cfg(not(feature = "bench-internals"))]
fn shared_core_search(_: &mut Criterion) {}

#[cfg(feature = "bench-internals")]
fn shared_core_filter(criterion: &mut Criterion) {
    let fixture = shared_core_absent_fixture();
    let (_, _, reader) = fixture
        .readers
        .iter()
        .find(|(version, _, _)| *version == "v3")
        .unwrap();
    let (filter_bytes, filter_keys, end_prefix) = reader.benchmark_core_filter(usize::MAX).unwrap();
    assert_eq!(end_prefix, 65_536);
    let filter = BinaryFuse8Ref::from_dma(&filter_bytes[..20], &filter_bytes[20..]);
    eprintln!(
        "shared_core_filter: {filter_keys} keys, {} serialized bytes",
        filter_bytes.len()
    );

    let mut group = criterion.benchmark_group("shared_core_filter");
    group.sample_size(20);
    group.nresamples(1_000);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    for (name, keys) in &fixture.workloads {
        let cores = keys.iter().map(|key| key.core).collect::<Vec<_>>();
        let mut filtered = cores.clone();
        filtered.sort_unstable();
        filtered.dedup();
        filtered.retain(|core| filter.contains(&u64::from(*core)));
        let control = reader.benchmark_resolve_sorted_cores(&cores).unwrap();
        let filtered_results = reader.benchmark_resolve_sorted_cores(&filtered).unwrap();
        assert_eq!(
            filtered_results
                .iter()
                .map(|group| (group.key(), group.member_count(), group.occurrence_count()))
                .collect::<Vec<_>>(),
            control
                .iter()
                .map(|group| (group.key(), group.member_count(), group.occurrence_count()))
                .collect::<Vec<_>>(),
            "{name}"
        );
        eprintln!(
            "shared_core_filter {name}: {} requested, {} survived, {} exact hits",
            cores.len(),
            filtered.len(),
            filtered_results.len()
        );

        group.throughput(Throughput::Elements(cores.len() as u64));
        group.bench_function(format!("{name}/exact"), |bencher| {
            bencher.iter(|| {
                let mut sorted = black_box(&cores).clone();
                sorted.sort_unstable();
                sorted.dedup();
                reader.benchmark_resolve_sorted_cores(&sorted).unwrap()
            })
        });
        group.bench_function(format!("{name}/filter_exact"), |bencher| {
            bencher.iter(|| {
                let mut filtered = black_box(&cores).clone();
                filtered.sort_unstable();
                filtered.dedup();
                filtered.retain(|core| filter.contains(&u64::from(*core)));
                reader.benchmark_resolve_sorted_cores(&filtered).unwrap()
            })
        });
    }
    group.finish();
}

#[cfg(not(feature = "bench-internals"))]
fn shared_core_filter(_: &mut Criterion) {}

#[cfg(feature = "bench-internals")]
fn shared_core_planning(criterion: &mut Criterion) {
    let fixture = shared_core_absent_fixture();
    let (_, path, _) = fixture
        .readers
        .iter()
        .find(|(version, _, _)| *version == "v4")
        .unwrap();
    let engine = TraceEngine::open_shared(path, None).unwrap();
    let mut group = criterion.benchmark_group("shared_core_planning");
    group.sample_size(20);
    group.nresamples(1_000);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    for (name, source) in fixture.workloads.iter().filter(|(name, _)| {
        matches!(
            *name,
            "spread/absent" | "spread/50pct_present" | "spread/present"
        )
    }) {
        let directory = source.iter().map(|key| key.core).collect::<Vec<_>>();
        for (reuse, repeats) in [("low_reuse_1", 1), ("high_reuse_16", 16)] {
            let keys = (0..repeats)
                .flat_map(|_| directory.iter().copied())
                .collect::<Vec<_>>();
            let (_, early_counts) = engine.benchmark_core_planning(&keys, true).unwrap();
            let (_, late_counts) = engine.benchmark_core_planning(&keys, false).unwrap();
            eprintln!(
                "shared_core_planning {name}/{reuse} [attempted, covered, rejected, uncovered, retained, planned, tasks]: early={early_counts:?}, late={late_counts:?}"
            );
            if *name == "spread/present" {
                assert_eq!(
                    [
                        early_counts[0],
                        early_counts[4],
                        early_counts[5],
                        early_counts[6],
                    ],
                    [
                        late_counts[0],
                        late_counts[4],
                        late_counts[5],
                        late_counts[6],
                    ],
                );
            }
            group.throughput(Throughput::Elements(keys.len() as u64));
            for (mode, early) in [("early", true), ("late", false)] {
                group.bench_function(format!("{name}/{reuse}/{mode}"), |bencher| {
                    bencher.iter(|| {
                        engine
                            .benchmark_core_planning(black_box(&keys), black_box(early))
                            .unwrap()
                    })
                });
            }
        }
    }
    group.finish();
}

#[cfg(not(feature = "bench-internals"))]
fn shared_core_planning(_: &mut Criterion) {}

#[cfg(feature = "bench-internals")]
fn shared_context_task(criterion: &mut Criterion) {
    let fixture = shared_core_absent_fixture();
    let reader = &fixture
        .readers
        .iter()
        .find(|(version, _, _)| *version == "v4")
        .unwrap()
        .2;
    let keys = &fixture
        .workloads
        .iter()
        .find(|(name, _)| *name == "spread/present")
        .unwrap()
        .1;
    let cores = reader.find_many(&keys[..256]).unwrap();
    let requests = cores
        .into_iter()
        .flatten()
        .map(|core| {
            (
                core,
                [
                    core.key(),
                    SharedKey {
                        length: 21,
                        context: 0,
                        ..core.key()
                    },
                    SharedKey {
                        length: 31,
                        context: 0,
                        ..core.key()
                    },
                ],
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        reader.benchmark_context_task(&requests, false).unwrap(),
        reader.benchmark_context_task(&requests, true).unwrap()
    );
    let mut group = criterion.benchmark_group("shared_context_task");
    group.sample_size(20);
    group.nresamples(1_000);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    for (name, batched) in [("public_per_core", false), ("checked_task", true)] {
        group.bench_function(name, |b| {
            b.iter(|| {
                reader
                    .benchmark_context_task(black_box(&requests), batched)
                    .unwrap()
            })
        });
    }
    group.finish();
}

#[cfg(not(feature = "bench-internals"))]
fn shared_context_task(_: &mut Criterion) {}

fn shared_core_filter_open(criterion: &mut Criterion) {
    let fixture = shared_core_absent_fixture();
    let (_, path, resident) = fixture
        .readers
        .iter()
        .find(|(version, _, _)| *version == "v4")
        .unwrap();
    let keys = &fixture
        .workloads
        .iter()
        .find(|(name, _)| *name == "spread/50pct_present")
        .unwrap()
        .1;
    let observed = SharedReader::open_observed(path).unwrap();
    black_box(observed.find_many(keys).unwrap());
    let stats = observed.stats();
    eprintln!(
        "shared_core_filter_open v4: filter_setup_ns={:?}, filter_owned_copies={}, identity_checks={}, reads={:?}",
        stats.filter_setup_ns, stats.filter_owned_copies, stats.file.identity_checks, stats.file,
    );
    let mut group = criterion.benchmark_group("shared_core_filter_open");
    group.sample_size(20);
    group.nresamples(1_000);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    group.bench_function("v4/fresh_open_setup", |bencher| {
        bencher.iter(|| SharedReader::open_observed(black_box(path)).unwrap())
    });
    group.throughput(Throughput::Elements(keys.len() as u64));
    group.bench_function("v4/resident_reads", |bencher| {
        bencher.iter(|| resident.find_many(black_box(keys)).unwrap())
    });
    group.finish();
}

#[cfg(feature = "bench-internals")]
fn shared_posting_preparation(criterion: &mut Criterion) {
    let fixture = shared_core_absent_fixture();
    let (_, path, _) = fixture
        .readers
        .iter()
        .find(|(version, _, _)| *version == "v3")
        .unwrap();
    let present = fixture
        .workloads
        .iter()
        .find(|(name, _)| *name == "spread/present")
        .unwrap()
        .1
        .iter()
        .map(|key| key.core)
        .collect::<Vec<_>>();
    let engine = TraceEngine::open_shared(path, None).unwrap();
    let workloads = [16usize, 256, 1024, 2048, 4096]
        .map(|count| (count, sampled_present_cores(&present, count)));
    eprintln!(
        "shared_posting_preparation boundary: checked key resolution plus posting planning, member fill, position fill, and drop"
    );
    let mut group = criterion.benchmark_group("shared_posting_preparation");
    group.sample_size(20);
    group.nresamples(1_000);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    for (count, cores) in workloads {
        let contexts = cores.into_iter().map(SharedKey::core).collect::<Vec<_>>();
        let keys = contexts
            .iter()
            .map(|key| key.packed().unwrap())
            .collect::<Vec<_>>();
        let fresh_engine = TraceEngine::open_shared(path, None).unwrap();
        let (retained, stats) = fresh_engine
            .benchmark_posting_preparation(&keys, true, true)
            .unwrap();
        black_box(&retained);
        drop(retained);
        let fresh_reader = SharedReader::open_observed(path).unwrap();
        black_box(fresh_reader.find_many(&contexts).unwrap());
        eprintln!(
            "shared_posting_preparation fresh-reader {count}: [member_tasks, position_tasks, member_wall_ns, position_wall_ns, member_cpu_ns, position_cpu_ns, member_rows, position_rows, scratch_bytes, retained_bytes, peak_parallel, plan_hash]={stats:?}; reader={:?}",
            fresh_reader.stats(),
        );

        group.throughput(Throughput::Elements(keys.len() as u64));
        for (mode, parallel) in [("serial", false), ("parallel", true)] {
            group.bench_function(format!("{count}/{mode}"), |bencher| {
                bencher.iter(|| {
                    engine
                        .benchmark_posting_preparation(black_box(&keys), black_box(parallel), false)
                        .unwrap()
                })
            });
        }
    }
    group.finish();
}

#[cfg(not(feature = "bench-internals"))]
fn shared_posting_preparation(_: &mut Criterion) {}

#[cfg(feature = "bench-internals")]
fn shared_prepare(criterion: &mut Criterion) {
    let (directory, _reader, _) = shared_lookup_fixture();
    let v1_path = directory.path().join("target.shared");
    let v2_path = directory.path().join("prepare-packed.shared");
    let v3_path = directory.path().join("prepare-compact.shared");
    let v4_path = directory.path().join("prepare-filtered.shared");
    jam_rs::shared_pack::repack_shared_index(&v1_path, &v2_path).unwrap();
    jam_rs::shared_pack::repack_shared_cores(&v2_path, &v3_path).unwrap();
    jam_rs::shared_pack::add_shared_core_filter(&v3_path, &v4_path, usize::MAX).unwrap();
    let all_v1 = directory.path().join("prepare-all.shared");
    let all_v2 = directory.path().join("prepare-all-packed.shared");
    let all_v3 = directory.path().join("prepare-all-compact.shared");
    let all_v4 = directory.path().join("prepare-all-filtered.shared");
    build_shared_index(directory.path().join("target.jidx"), &all_v1, 1).unwrap();
    jam_rs::shared_pack::repack_shared_index(&all_v1, &all_v2).unwrap();
    jam_rs::shared_pack::repack_shared_cores(&all_v2, &all_v3).unwrap();
    jam_rs::shared_pack::add_shared_core_filter(&all_v3, &all_v4, usize::MAX).unwrap();
    let engines = [
        ("v1", TraceEngine::open_shared(v1_path, None).unwrap()),
        ("v3", TraceEngine::open_shared(v3_path, None).unwrap()),
        ("v4", TraceEngine::open_shared(v4_path, None).unwrap()),
        ("v4all", TraceEngine::open_shared(all_v4, None).unwrap()),
    ];
    let present_2k = sequence(2_000);
    let absent_2k = sequence_from_state(2_000, 101);
    let present_64k = sequence(64_000);
    let mut ambiguous_lowercase_2k = sequence(2_000);
    ambiguous_lowercase_2k.make_ascii_lowercase();
    ambiguous_lowercase_2k[960..1_040].fill(b'N');
    let mut mixed_64k = sequence(32_000);
    mixed_64k.extend(sequence_from_state(32_000, 103));
    let absent_250k = sequence_from_state(250_000, 107);
    let mut mixed_250k = sequence(64_000);
    mixed_250k.extend(sequence_from_state(186_000, 109));
    let boundary = 128 * 1024 * 1024
        / (std::mem::size_of::<u32>()
            + 2 * std::mem::size_of::<jam_rs::shared_reader::SharedGroup>());
    let below_boundary = sequence_from_state(boundary - boundary / 100 + 14, 113);
    let above_boundary = sequence_from_state(boundary + boundary / 100 + 14, 127);
    let repeated_present = sequence(127)
        .into_iter()
        .cycle()
        .take(64_000)
        .collect::<Vec<_>>();
    let workloads = [
        ("64kb/repeated/linear", repeated_present, false),
        ("boundary/below/linear", below_boundary, false),
        ("boundary/above/linear", above_boundary, false),
        ("2kb/present/linear", present_2k, false),
        ("2kb/absent/linear", absent_2k, false),
        (
            "2kb/ambiguous_lowercase/circular",
            ambiguous_lowercase_2k,
            true,
        ),
        ("64kb/mixed/linear", mixed_64k, false),
        ("64kb/present/circular", present_64k, true),
        ("250kb/absent/linear", absent_250k, false),
        ("250kb/mixed/circular", mixed_250k, true),
    ];
    let mut group = criterion.benchmark_group("shared_prepare");
    group.sample_size(20);
    group.nresamples(1_000);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    for (format, engine) in &engines {
        for (name, query, circular) in &workloads {
            if *format == "v1" && name.starts_with("boundary/") {
                continue;
            }
            let config = TraceConfig {
                circular: *circular,
                use_sketch: false,
                ..TraceConfig::default()
            };
            engine.benchmark_prepare(query, config).unwrap();
            group.throughput(Throughput::Bytes(query.len() as u64));
            let benchmark = if *format == "v1" {
                (*name).to_owned()
            } else {
                format!("{format}/{name}")
            };
            if format.starts_with("v4") {
                group.bench_function(format!("{benchmark}/directory"), |bencher| {
                    bencher.iter(|| {
                        engine
                            .benchmark_prepare_directory(black_box(query), black_box(config))
                            .unwrap()
                    })
                });
            }
            group.bench_function(benchmark, |bencher| {
                bencher.iter(|| {
                    engine
                        .benchmark_prepare(black_box(query), black_box(config))
                        .unwrap()
                })
            });
        }
    }
    group.finish();
}

#[cfg(not(feature = "bench-internals"))]
fn shared_prepare(_: &mut Criterion) {}

fn shared_resolved_handle(criterion: &mut Criterion) {
    let (directory, reference, workloads) = shared_lookup_fixture();
    let packed_path = directory.path().join("resolved-packed.shared");
    jam_rs::shared_pack::repack_shared_index(directory.path().join("target.shared"), &packed_path)
        .unwrap();
    let packed = SharedReader::open(packed_path).unwrap();
    let compact_path = directory.path().join("resolved-compact.shared");
    jam_rs::shared_pack::repack_shared_cores(
        directory.path().join("resolved-packed.shared"),
        &compact_path,
    )
    .unwrap();
    let compact = SharedReader::open(compact_path).unwrap();
    let repeated_key = workloads
        .iter()
        .find(|(name, _)| *name == "high_multiplicity")
        .unwrap()
        .1[0];
    let singleton_key = workloads
        .iter()
        .find(|(name, _)| *name == "low_reuse")
        .unwrap()
        .1
        .iter()
        .copied()
        .find(|&key| {
            reference
                .find(key)
                .unwrap()
                .is_some_and(|group| group.occurrence_count() == 1)
        })
        .unwrap();
    let expected = reference
        .find(repeated_key)
        .unwrap()
        .and_then(|group| reference.members(group).unwrap().into_iter().next())
        .map(|member| {
            let group = reference.find(repeated_key).unwrap().unwrap();
            reference.member_occurrences(group, member).unwrap()
        })
        .unwrap();
    let mut group = criterion.benchmark_group("shared_resolved_handle");
    group.sample_size(20);
    group.nresamples(1_000);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    for (format, reader) in [("v1", &reference), ("packed", &packed), ("v3", &compact)] {
        let repeated = reader.find(repeated_key).unwrap().unwrap();
        let members = reader.members(repeated).unwrap();
        assert_eq!(members.len(), 1);
        let member = members[0];
        assert_eq!(
            reader.member_occurrences(repeated, member).unwrap(),
            expected
        );
        assert!(member.occurrence_count() >= 32);
        group.throughput(Throughput::Elements(repeated.member_count() as u64));
        group.bench_function(format!("members/{format}"), |bencher| {
            bencher.iter(|| reader.members(black_box(repeated)).unwrap())
        });
        for (location, start) in [
            ("first", 0),
            ("middle", member.occurrence_count() / 2),
            ("last", member.occurrence_count() - 16),
        ] {
            group.throughput(Throughput::Elements(16));
            group.bench_function(format!("positions/{format}/{location}"), |bencher| {
                bencher.iter(|| {
                    reader
                        .occurrence_block(
                            black_box(repeated),
                            black_box(member),
                            black_box(start),
                            16,
                        )
                        .unwrap()
                })
            });
        }
        let singleton = reader.find(singleton_key).unwrap().unwrap();
        let singleton_member = reader.members(singleton).unwrap()[0];
        group.throughput(Throughput::Elements(1));
        group.bench_function(format!("positions/{format}/singleton"), |bencher| {
            bencher.iter(|| {
                reader
                    .occurrence_block(black_box(singleton), black_box(singleton_member), 0, 1)
                    .unwrap()
            })
        });
    }
    group.finish();
}

fn shared_geometry(criterion: &mut Criterion) {
    let (directory, _reader, _) = shared_lookup_fixture();
    let v1_path = directory.path().join("target.shared");
    let v2_path = directory.path().join("geometry-packed.shared");
    let v3_path = directory.path().join("geometry-compact.shared");
    jam_rs::shared_pack::repack_shared_index(&v1_path, &v2_path).unwrap();
    jam_rs::shared_pack::repack_shared_cores(&v2_path, &v3_path).unwrap();
    let engines = [
        ("v1", TraceEngine::open_shared(v1_path, None).unwrap()),
        ("v2", TraceEngine::open_shared(v2_path, None).unwrap()),
        ("v3", TraceEngine::open_shared(v3_path, None).unwrap()),
    ];
    let repeated = sequence(127);
    let exact = repeated[..64].to_vec();
    let mut mixed = exact.clone();
    mixed.extend(sequence_from_state(32, 113));
    let workloads = [
        ("64bp/repeated/linear", exact.clone(), false),
        ("64bp/repeated/circular", exact, true),
        ("96bp/mixed_repeated/linear", mixed, false),
    ];
    let mut group = criterion.benchmark_group("shared_geometry");
    group.sample_size(20);
    group.nresamples(1_000);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    for (format, engine) in &engines {
        for (name, query, circular) in &workloads {
            let config = TraceConfig {
                circular: *circular,
                use_sketch: false,
                ..TraceConfig::default()
            };
            let expected = engine.search("geometry", query, config).unwrap();
            assert_eq!(expected.completion, SearchCompletion::Complete);
            assert_eq!(expected.query_length, query.len() as u64);
            assert!(
                expected
                    .metagenomes
                    .iter()
                    .any(|metagenome| !metagenome.mosaic.primary.is_empty())
            );
            group.throughput(Throughput::Bytes(query.len() as u64));
            group.bench_function(format!("{format}/{name}"), |bencher| {
                bencher.iter(|| {
                    engine
                        .search("geometry", black_box(query), black_box(config))
                        .unwrap()
                })
            });
        }
    }
    group.finish();
}

#[cfg(feature = "bench-internals")]
fn shared_threads(criterion: &mut Criterion) {
    let SharedCoreFixture { directory, .. } = shared_core_absent_fixture();
    let engine =
        TraceEngine::open_shared(directory.path().join("large-target-v3.shared"), None).unwrap();
    let keys = (0..98_304u32)
        .map(|ordinal| {
            SharedKey::core(ordinal.wrapping_mul(506_952_113) & ((1 << 30) - 1))
                .packed()
                .unwrap()
        })
        .collect::<Vec<_>>();
    let pools = [1, 4, 8, 16].map(|threads| {
        (
            threads,
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap(),
        )
    });
    let mut expected = None;
    for (threads, pool) in &pools {
        let stats = pool.install(|| {
            let (retained, stats) = engine.benchmark_lookup(&keys, true).unwrap();
            black_box(&retained);
            drop(retained);
            stats
        });
        assert!(
            stats[0] >= 3,
            "{threads} threads produced {} tasks",
            stats[0]
        );
        assert_eq!(*expected.get_or_insert(stats), stats, "{threads} threads");
    }
    println!("shared_threads observed [tasks, plan_hash, successes]: {expected:?}");

    let mut group = criterion.benchmark_group("shared_threads");
    group.sample_size(20);
    group.nresamples(1_000);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    group.throughput(Throughput::Elements(keys.len() as u64));
    for (threads, pool) in &pools {
        group.bench_function(format!("{threads}_threads"), |bencher| {
            bencher.iter(|| {
                pool.install(|| {
                    let (retained, stats) =
                        engine.benchmark_lookup(black_box(&keys), false).unwrap();
                    black_box(&retained);
                    drop(retained);
                    stats
                })
            })
        });
    }
    group.finish();
}

#[cfg(not(feature = "bench-internals"))]
fn shared_threads(_: &mut Criterion) {}

fn shared_packed(criterion: &mut Criterion) {
    let (directory, reference, workloads) = shared_lookup_fixture();
    let path = directory.path().join("packed.shared");
    jam_rs::shared_pack::repack_shared_index(directory.path().join("target.shared"), &path)
        .unwrap();
    let packed = SharedReader::open(path).unwrap();
    let mut group = criterion.benchmark_group("shared_packed");
    group.sample_size(20);
    group.nresamples(1_000);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    for (name, keys) in &workloads {
        let counts = |reader: &SharedReader| {
            reader
                .find_many(keys)
                .unwrap()
                .into_iter()
                .map(|group| {
                    group.map(|group| (group.key(), group.member_count(), group.occurrence_count()))
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(counts(&reference), counts(&packed));
        for (format, reader) in [("v1", &reference), ("packed", &packed)] {
            group.bench_function(format!("{name}/{format}"), |bencher| {
                bencher.iter(|| reader.find_many(black_box(keys)).unwrap())
            });
        }
    }
    let common = workloads
        .iter()
        .find(|(name, _)| *name == "high_multiplicity")
        .unwrap()
        .1[0];
    let left = reference.find(common).unwrap().unwrap();
    let right = packed.find(common).unwrap().unwrap();
    let member = reference.members(left).unwrap()[0];
    let other = packed.members(right).unwrap()[0];
    let start = member.occurrence_count().saturating_sub(16);
    assert_eq!(
        reference.occurrence_block(left, member, start, 16).unwrap(),
        packed.occurrence_block(right, other, start, 16).unwrap()
    );
    for (format, reader, selected, member) in [
        ("v1", &reference, left, member),
        ("packed", &packed, right, other),
    ] {
        group.bench_function(format!("late_placement_block/{format}"), |bencher| {
            bencher.iter(|| {
                reader
                    .occurrence_block(selected, member, black_box(start), 16)
                    .unwrap()
            })
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    affine_alignment,
    trace_endpoint,
    owner_postings,
    shared_lookup,
    shared_core_absent,
    shared_core_search,
    shared_core_filter,
    shared_core_planning,
    shared_context_task,
    shared_core_filter_open,
    shared_posting_preparation,
    shared_prepare,
    shared_resolved_handle,
    shared_geometry,
    shared_threads,
    shared_packed
);
criterion_main!(benches);
