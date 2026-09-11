use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use jam_rs::alignment::{AlignmentConfig, AlignmentWorkspace};
use jam_rs::jidx::sha256;
use jam_rs::jidx_writer::{ContigInput, JidxInput, JidxWriter, MetagenomeInput};
use jam_rs::owner_postings;
use jam_rs::shared_reader::SharedReader;
use jam_rs::shared_seed::{HAS_CONTEXT_21, HAS_CONTEXT_31, SharedKey, select_shared_seeds};
use jam_rs::shared_writer::build_shared_index;
#[cfg(feature = "bench-internals")]
use jam_rs::trace::{TraceConfig, TraceEngine};
use noodles_bgzf::{self as bgzf, gzi};
use std::fs::File;
use std::hint::black_box;
use std::io::Write;
use std::time::Duration;

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

fn shared_core_absent_fixture() -> (
    tempfile::TempDir,
    SharedReader,
    Vec<(&'static str, Vec<SharedKey>)>,
) {
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
    let shared_path = directory.path().join("large-target.shared");
    let build = build_shared_index(&metadata_path, &shared_path, 1).unwrap();
    assert!(build.core_count >= 1_000_000, "{} cores", build.core_count);

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
    (
        directory,
        SharedReader::open(shared_path).unwrap(),
        workloads,
    )
}

fn shared_core_absent(criterion: &mut Criterion) {
    let fixture_started = std::time::Instant::now();
    let (_directory, reader, workloads) = shared_core_absent_fixture();
    eprintln!(
        "shared_core_absent fixture setup: {:.3}s",
        fixture_started.elapsed().as_secs_f64()
    );
    let mut group = criterion.benchmark_group("shared_core_absent");
    group.sample_size(20);
    group.nresamples(1_000);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    for (name, keys) in workloads {
        let expected = keys
            .iter()
            .filter(|&&key| reader.find(key).unwrap().is_some())
            .count();
        assert_eq!(
            reader
                .find_many(&keys)
                .unwrap()
                .into_iter()
                .flatten()
                .count(),
            expected
        );
        group.throughput(Throughput::Elements(keys.len() as u64));
        group.bench_function(name, |bencher| {
            bencher.iter(|| reader.find_many(black_box(&keys)).unwrap())
        });
    }
    group.finish();
}

#[cfg(feature = "bench-internals")]
fn shared_prepare(criterion: &mut Criterion) {
    let (directory, _reader, _) = shared_lookup_fixture();
    let engine = TraceEngine::open_shared(directory.path().join("target.shared"), None).unwrap();
    let present_2k = sequence(2_000);
    let absent_2k = sequence_from_state(2_000, 101);
    let present_64k = sequence(64_000);
    let mut mixed_64k = sequence(32_000);
    mixed_64k.extend(sequence_from_state(32_000, 103));
    let absent_250k = sequence_from_state(250_000, 107);
    let mut mixed_250k = sequence(64_000);
    mixed_250k.extend(sequence_from_state(186_000, 109));
    let workloads = [
        ("2kb/present/linear", present_2k, false),
        ("2kb/absent/linear", absent_2k, false),
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
    for (name, query, circular) in workloads {
        let config = TraceConfig {
            circular,
            use_sketch: false,
            ..TraceConfig::default()
        };
        engine.benchmark_prepare(&query, config).unwrap();
        group.throughput(Throughput::Bytes(query.len() as u64));
        group.bench_function(name, |bencher| {
            bencher.iter(|| {
                engine
                    .benchmark_prepare(black_box(&query), black_box(config))
                    .unwrap()
            })
        });
    }
    group.finish();
}

#[cfg(not(feature = "bench-internals"))]
fn shared_prepare(_: &mut Criterion) {}

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
    owner_postings,
    shared_lookup,
    shared_core_absent,
    shared_prepare,
    shared_packed
);
criterion_main!(benches);
