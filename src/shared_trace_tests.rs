use crate::alignment::{Alignment, AlignmentConfig, AlignmentWorkspace, EditOperation, Strand};
use crate::cli::handlers::{TraceArgs, TraceInput, handle_trace_command};
use crate::jidx::sha256;
use crate::jidx_writer::{ContigInput, JidxInput, JidxWriter, MetagenomeInput};
use crate::mosaic::Fragment;
use crate::shared_seed::{SharedSeed, context_seed, select_shared_seeds};
use crate::shared_writer::build_shared_index;
use crate::trace::{MetagenomeTrace, TraceConfig, TraceEngine};
use needletail::Sequence;
use noodles_bgzf::{self as bgzf, gzi};
use std::fs::File;
use std::io::Write;

fn dna(mut state: u64, length: usize) -> Vec<u8> {
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            b"ACGT"[(state >> 62) as usize]
        })
        .collect()
}

fn xorshift_dna(mut state: u64, length: usize) -> Vec<u8> {
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state & 3) as usize]
        })
        .collect()
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

fn periodic(sequence: &[u8], phase: usize) -> Vec<u8> {
    let mut output = sequence.to_vec();
    for position in (phase..output.len()).step_by(20) {
        output[position] = match output[position] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
    }
    output
}

fn write_bgzf(directory: &std::path::Path, name: &str, sequence: &[u8]) -> MetagenomeInput {
    let path = directory.join(format!("{name}.bgz"));
    let mut raw = b">contig\n".to_vec();
    for line in sequence.chunks(80) {
        raw.extend_from_slice(line);
        raw.push(b'\n');
    }
    let mut writer = bgzf::io::Writer::new(File::create(&path).unwrap());
    writer.write_all(&raw).unwrap();
    writer.finish().unwrap();
    let gzi_path = directory.join(format!("{name}.gzi"));
    gzi::fs::write(&gzi_path, &gzi::Index::default()).unwrap();
    let bytes = std::fs::read(&path).unwrap();
    MetagenomeInput {
        name: name.to_owned(),
        bgzf_uri: path.to_str().unwrap().to_owned(),
        bgzf_bytes: bytes.len() as u64,
        bgzf_sha256: sha256(&bytes),
        gzi: std::fs::read(gzi_path).unwrap(),
    }
}

fn write_queries(path: &std::path::Path, rows: &[(&str, &[u8])]) {
    let mut writer = File::create(path).unwrap();
    for (header, sequence) in rows {
        writeln!(writer, ">{header}").unwrap();
        writer.write_all(sequence).unwrap();
        writer.write_all(b"\n").unwrap();
    }
}

fn run_topology_cli(
    shared: &std::path::Path,
    query: &std::path::Path,
    output: &std::path::Path,
) -> anyhow::Result<()> {
    handle_trace_command(TraceArgs {
        query: query.to_owned(),
        input: TraceInput::Shared {
            path: shared.to_owned(),
            read_stats: None,
            query_topology_header: true,
        },
        audit_index: false,
        output: output.to_owned(),
        query_id: None,
        config: TraceConfig {
            use_sketch: false,
            circular: true,
            ..TraceConfig::default()
        },
        s3: None,
        force: false,
    })
}

fn without_json_read_accounting(mut result: serde_json::Value) -> serde_json::Value {
    for metagenome in result["metagenomes"].as_array_mut().unwrap() {
        for field in [
            "compressed_bytes_read",
            "range_requests",
            "bgzf_blocks_decoded",
        ] {
            metagenome[field] = serde_json::json!(0);
        }
    }
    result
}

#[test]
fn default_band_retains_known_gap256_bounded_counterexample() {
    let query = xorshift_dna(0x9e37_79b9_7f4a_7c15 ^ 1_000, 1_000);
    let mut homolog = query.clone();
    homolog.splice(500..500, xorshift_dna(0xabcd_dcba ^ 256, 256));
    let mut target = xorshift_dna(0x1234_5678_9abc_def0 ^ 4_096, 4_096);
    target[1_420..2_676].copy_from_slice(&homolog);
    let align = |band_width| {
        AlignmentWorkspace::default()
            .align(
                &query,
                &target,
                AlignmentConfig {
                    diagonal_offset: 1_420,
                    band_width,
                    max_cells: 20_000_000,
                    ..AlignmentConfig::default()
                },
            )
            .unwrap()
    };
    let default_band = align(128);
    assert_eq!(default_band.score, 1_000);
    assert_eq!(default_band.query_interval.end, 500);
    assert_eq!(default_band.target_interval.end, 1_920);
    assert_eq!(default_band.cigar, "500=");

    let wider_bounded = align(512);
    assert_eq!(wider_bounded.score, 1_739);
    assert_eq!(wider_bounded.query_interval.end, 1_000);
    assert_eq!(wider_bounded.target_interval.end, 2_676);
    assert_eq!(wider_bounded.cigar, "500=256I500=");
    assert!(default_band.score < wider_bounded.score);
}

fn matching_anchors(query: &[u8], seeds: &[SharedSeed]) -> (usize, usize) {
    let mut total = 0;
    let mut strong = 0;
    for seed in seeds {
        let position = seed.position as usize;
        let Some((_, core, orientation)) = query
            .get(position..position.saturating_add(15))
            .and_then(|word| word.bit_kmers(15, true).next())
        else {
            continue;
        };
        if core.0 as u32 != seed.core || orientation != seed.canonical_orientation() {
            continue;
        }
        let Some(query_seed) = context_seed(query, position, core.0 as u32, orientation, false)
        else {
            continue;
        };
        if query_seed.key(15) != seed.key(15) {
            continue;
        }
        total += 1;
        strong += usize::from([31, 21].into_iter().any(|length| {
            query_seed
                .key(length)
                .is_some_and(|key| Some(key) == seed.key(length))
        }));
    }
    (total, strong)
}

fn mixed_query(target: &[u8], seeds: &[SharedSeed]) -> (Vec<u8>, usize, usize) {
    for phase in 0..20 {
        let changed = periodic(target, phase);
        for seed in seeds {
            let position = seed.position as usize;
            if position < 3 || position + 18 > target.len() {
                continue;
            }
            let mut query = changed.clone();
            query[position - 3..position + 18]
                .copy_from_slice(&target[position - 3..position + 18]);
            let (total, strong) = matching_anchors(&query, seeds);
            if total >= 3 && strong == 1 {
                return (query, total, strong);
            }
        }
    }
    panic!("no deterministic mixed-strength query");
}

fn fragments(trace: &MetagenomeTrace) -> impl Iterator<Item = &Fragment> {
    trace
        .mosaic
        .primary
        .iter()
        .map(|selected| &selected.fragment)
        .chain(trace.mosaic.alternatives.iter())
}

fn alignment<'a>(result: &'a crate::trace::TraceResult, name: &str) -> &'a Alignment {
    fragments(
        result
            .metagenomes
            .iter()
            .find(|trace| trace.name == name)
            .unwrap(),
    )
    .max_by_key(|fragment| fragment.alignment.score)
    .map(|fragment| &fragment.alignment)
    .unwrap()
}

fn validate_score(alignment: &Alignment) {
    alignment.validate_cigar().unwrap();
    let score = alignment
        .edit_script
        .iter()
        .map(|run| match run.operation {
            EditOperation::Equal => 2 * run.length as i32,
            EditOperation::Substitution => -3 * run.length as i32,
            EditOperation::Insertion | EditOperation::Deletion => -5 - run.length as i32,
        })
        .sum::<i32>();
    assert_eq!(alignment.score, score);
}

fn without_read_accounting(mut result: crate::trace::TraceResult) -> crate::trace::TraceResult {
    for metagenome in &mut result.metagenomes {
        metagenome.compressed_bytes_read = 0;
        metagenome.range_requests = 0;
        metagenome.bgzf_blocks_decoded = 0;
    }
    result
}

#[test]
fn shared_index_traces_strong_weak_mixed_reverse_and_circular_queries() {
    shared_trace_fixture(false);
    shared_trace_fixture(true);
}

fn shared_trace_fixture(packed: bool) {
    let directory = tempfile::tempdir().unwrap();
    let exact_target = dna(11, 800);
    let mixed_target = dna(29, 800);
    let weak_target = periodic(&exact_target, 0);
    let absent_target = exact_target[..720].to_vec();
    let mut embedded_target = dna(41, 1_600);
    embedded_target[400..1_200].copy_from_slice(&exact_target);
    let reverse_target = reverse_complement(&embedded_target);
    let target_rows = [
        ("absent", absent_target.as_slice()),
        ("embedded", embedded_target.as_slice()),
        ("exact", exact_target.as_slice()),
        ("mixed", mixed_target.as_slice()),
        ("reverse", reverse_target.as_slice()),
        ("weak", weak_target.as_slice()),
    ];
    let reference = directory.path().join("reference.jidx");
    let mut writer = JidxWriter::new(
        &reference,
        &JidxInput {
            k: 15,
            rescue_k15: false,
            minimizer_window: 64,
            jam_sha256: [1; 32],
            manifest_sha256: [2; 32],
        },
    )
    .unwrap();
    for (name, sequence) in target_rows {
        writer
            .begin_metagenome(write_bgzf(directory.path(), name, sequence))
            .unwrap();
        writer
            .begin_contig(ContigInput {
                name: "contig".into(),
                length: sequence.len() as u64,
                fasta_offset: 8,
                line_bases: 80,
                line_width: 81,
            })
            .unwrap();
    }
    writer.finish().unwrap();
    let shared = directory.path().join("targets.shared");
    let stats = build_shared_index(&reference, &shared, 64).unwrap();
    let shared = if packed {
        let output = directory.path().join("targets.packed.shared");
        crate::shared_pack::repack_shared_index(&shared, &output).unwrap();
        output
    } else {
        shared
    };
    assert_eq!(stats.source_bases, 6_320);
    std::fs::remove_file(&reference).unwrap();
    assert!(!reference.exists());
    for name in ["absent", "embedded", "exact", "mixed", "reverse", "weak"] {
        let gzi = directory.path().join(format!("{name}.gzi"));
        std::fs::remove_file(&gzi).unwrap();
        assert!(!gzi.exists());
    }

    let exact_query = exact_target.clone();
    let reverse_query = reverse_complement(&exact_target);
    let circular_query = [exact_target[600..].to_vec(), exact_target[..600].to_vec()].concat();
    let mixed_seeds = select_shared_seeds(&mixed_target, 64).unwrap();
    let (mixed_query, mixed_total, mixed_strong) = mixed_query(&mixed_target, &mixed_seeds);
    assert!(mixed_total >= 3);
    assert_eq!(mixed_strong, 1);
    let weak_seeds = select_shared_seeds(&weak_target, 64).unwrap();
    let (weak_total, weak_strong) = matching_anchors(&exact_query, &weak_seeds);
    assert!(weak_total >= 2);
    assert_eq!(weak_strong, 0);

    let engine = TraceEngine::open_shared(&shared, None).unwrap();
    let index = crate::trace_index::TraceIndex::Shared(Box::new(
        crate::shared_reader::SharedReader::open(&shared).unwrap(),
    ));
    assert!(crate::trace_batch::lookup_bytes(&index, usize::MAX, 1).is_none());
    assert!(crate::trace_batch::lookup_bytes(&index, 1, usize::MAX).is_none());
    engine.verify_index().unwrap();
    let linear = TraceConfig {
        use_sketch: false,
        circular: false,
        ..TraceConfig::default()
    };
    let batch = engine
        .search_batch(
            &[
                ("exact".into(), exact_query.clone()),
                ("mixed".into(), mixed_query),
            ],
            linear,
        )
        .unwrap();
    let exact = alignment(&batch[0], "exact");
    assert_eq!(
        exact.target_interval,
        crate::alignment::Interval::new(0, 800).unwrap()
    );
    assert_eq!(exact.strand, Strand::Forward);
    assert_eq!(exact.identity(), 1.0);
    validate_score(exact);
    let embedded = alignment(&batch[0], "embedded");
    assert_eq!(
        embedded.target_interval,
        crate::alignment::Interval::new(400, 1_200).unwrap()
    );
    assert_eq!(embedded.strand, Strand::Forward);
    validate_score(embedded);
    let embedded_reverse = alignment(&batch[0], "reverse");
    assert_eq!(
        embedded_reverse.target_interval,
        crate::alignment::Interval::new(400, 1_200).unwrap()
    );
    assert_eq!(embedded_reverse.strand, Strand::Reverse);
    validate_score(embedded_reverse);
    let weak = alignment(&batch[0], "weak");
    assert_eq!(weak.strand, Strand::Forward);
    assert!(weak.query_interval.len() >= 760);
    assert!(weak.identity() >= 0.9);
    validate_score(weak);
    let mixed = alignment(&batch[1], "mixed");
    assert!(mixed.query_interval.len() >= 760);
    assert!(mixed.identity() >= 0.9);
    validate_score(mixed);

    let reverse = engine.search("reverse", &reverse_query, linear).unwrap();
    let reverse = alignment(&reverse, "exact");
    assert_eq!(
        reverse.target_interval,
        crate::alignment::Interval::new(0, 800).unwrap()
    );
    assert_eq!(reverse.strand, Strand::Reverse);
    validate_score(reverse);
    let circular = engine
        .search(
            "circular",
            &circular_query,
            TraceConfig {
                circular: true,
                ..linear
            },
        )
        .unwrap();
    let trace = circular
        .metagenomes
        .iter()
        .find(|trace| trace.name == "exact")
        .unwrap();
    let fragment = fragments(trace)
        .max_by_key(|fragment| fragment.alignment.score)
        .unwrap();
    assert_eq!(
        fragment.alignment.target_interval,
        crate::alignment::Interval::new(0, 800).unwrap()
    );
    assert_eq!(fragment.query_segments.len(), 2);
    validate_score(&fragment.alignment);
    for (name, strand) in [("embedded", Strand::Forward), ("reverse", Strand::Reverse)] {
        let embedded = circular
            .metagenomes
            .iter()
            .find(|trace| trace.name == name)
            .unwrap();
        let embedded = fragments(embedded)
            .max_by_key(|fragment| fragment.alignment.score)
            .unwrap();
        assert_eq!(
            embedded.alignment.target_interval,
            crate::alignment::Interval::new(400, 1_200).unwrap()
        );
        assert_eq!(embedded.alignment.strand, strand);
        assert_eq!(embedded.query_segments.len(), 2);
        validate_score(&embedded.alignment);
    }
    let absent = circular
        .metagenomes
        .iter()
        .find(|trace| trace.name == "absent")
        .unwrap();
    let absent = fragments(absent)
        .max_by_key(|fragment| fragment.alignment.score)
        .unwrap();
    assert_eq!(
        absent.alignment.target_interval,
        crate::alignment::Interval::new(0, 720).unwrap()
    );
    assert_eq!(absent.alignment.query_interval.len(), 720);
    validate_score(&absent.alignment);

    let topology_queries = [
        ("topology-linear".to_owned(), exact_query.clone()),
        ("topology-circular".to_owned(), circular_query.clone()),
    ];
    let topology_flags = [false, true];
    let expected_topologies = vec![
        engine
            .search("topology-linear", &exact_query, linear)
            .unwrap(),
        engine
            .search(
                "topology-circular",
                &circular_query,
                TraceConfig {
                    circular: true,
                    ..linear
                },
            )
            .unwrap(),
    ]
    .into_iter()
    .map(without_read_accounting)
    .collect::<Vec<_>>();
    for threads in [1, 4, 8] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let actual = pool
            .install(|| {
                assert!(
                    crate::trace_batch::lookup_bytes(
                        &index,
                        3 * crate::cli::handlers::SHARED_BATCH_QUERY_BASES,
                        64,
                    )
                    .unwrap()
                        > crate::trace_batch::lookup_budget(&index)
                );
                let oversized =
                    Vec::with_capacity(3 * crate::cli::handlers::SHARED_BATCH_QUERY_BASES);
                assert!(
                    crate::trace_batch::prepare_lookup(&index, oversized, 64, false)
                        .unwrap()
                        .is_none()
                );
                TraceEngine::open_shared(&shared, None)
                    .unwrap()
                    .search_batch_topologies(&topology_queries, linear, &topology_flags)
                    .unwrap()
            })
            .into_iter()
            .map(without_read_accounting)
            .collect::<Vec<_>>();
        assert_eq!(actual, expected_topologies, "thread count {threads}");
    }
    for (((id, sequence), &circular), expected) in topology_queries
        .iter()
        .zip(&topology_flags)
        .zip(&expected_topologies)
    {
        let actual = engine
            .search_without_batch(id, sequence, TraceConfig { circular, ..linear })
            .unwrap();
        assert_eq!(&without_read_accounting(actual), expected);
    }

    let cli_query = directory.path().join("mixed-topologies.fa");
    let cli_output = directory.path().join("mixed-topologies.jsonl");
    write_queries(
        &cli_query,
        &[
            ("topology-linear topology=linear", &exact_query),
            ("topology-circular topology=circular", &circular_query),
        ],
    );
    run_topology_cli(&shared, &cli_query, &cli_output).unwrap();
    let actual = std::fs::read_to_string(&cli_output)
        .unwrap()
        .lines()
        .map(|line| {
            without_json_read_accounting(serde_json::from_str::<serde_json::Value>(line).unwrap())
        })
        .collect::<Vec<_>>();
    let expected = expected_topologies
        .iter()
        .map(|result| without_json_read_accounting(serde_json::to_value(result).unwrap()))
        .collect::<Vec<_>>();
    assert_eq!(actual, expected);

    for (case, header, message) in [
        ("missing", "missing", "requires topology"),
        ("invalid", "invalid topology=unknown", "requires topology"),
        (
            "duplicate",
            "duplicate topology=linear topology=circular",
            "Duplicate query topology",
        ),
    ] {
        let query = directory.path().join(format!("{case}-topology.fa"));
        let output = directory.path().join(format!("{case}-topology.jsonl"));
        write_queries(&query, &[(header, &exact_query)]);
        let error = run_topology_cli(&shared, &query, &output).unwrap_err();
        assert!(error.to_string().contains(message));
        assert!(!output.exists());
    }

    let queries = [
        ("reuse".to_owned(), exact_query.clone()),
        ("reuse".to_owned(), exact_query),
    ];
    let batch_engine = TraceEngine::open_shared_observed(&shared, None, true).unwrap();
    let batch_results = batch_engine.search_batch(&queries, linear).unwrap();
    let batch_stats = batch_engine.batch_stats();
    let batch_reads = batch_engine.shared_read_stats().unwrap();
    let mut separate_results = Vec::new();
    let mut separate_core_inspections = 0;
    let mut separate_positions = 0;
    let mut separate_bgzf_decodes = 0;
    let mut separate_unique_keys = 0;
    let mut separate_cached_positions = 0;
    for (id, query) in &queries {
        let engine = TraceEngine::open_shared_observed(&shared, None, true).unwrap();
        separate_results.push(engine.search(id, query, linear).unwrap());
        let reads = engine.shared_read_stats().unwrap();
        let stats = engine.batch_stats();
        separate_core_inspections += reads.core_descriptor_inspections;
        separate_positions += reads.physical_positions_decoded;
        separate_bgzf_decodes += stats.bgzf_blocks_decoded;
        separate_unique_keys += stats.unique_keys;
        separate_cached_positions += stats.cached_positions;
    }
    assert_eq!(
        batch_results
            .into_iter()
            .map(without_read_accounting)
            .collect::<Vec<_>>(),
        separate_results
            .into_iter()
            .map(without_read_accounting)
            .collect::<Vec<_>>()
    );
    assert!(batch_reads.core_descriptor_inspections < separate_core_inspections);
    assert!(batch_reads.physical_positions_decoded < separate_positions);
    assert!(batch_stats.bgzf_blocks_decoded < separate_bgzf_decodes);
    assert!(batch_stats.unique_keys < separate_unique_keys);
    assert!(batch_stats.cached_positions < separate_cached_positions);
    assert!(batch_stats.bgzf_cache_hits > 0);
}
