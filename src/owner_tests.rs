use crate::alignment::Strand;
use crate::jidx_builder::{JidxBuildConfig, build_local_jidx};
use crate::jidx_reader::{JidxReader, SeedOccurrence};
use crate::owner_format::{
    HAS_METADATA, OWNER_CONTIG_SIZE, OWNER_DOCUMENT_SIZE, OWNER_HEADER_SIZE, OWNER_PAGE_SIZE,
    OWNER_SECTION_COUNT, OwnerHeader, OwnerSection, OwnerSectionDescriptor,
};
use crate::owner_postings::{OwnerKey, OwnerMember, OwnerOccurrence};
use crate::owner_reader::OwnerReader;
use crate::owner_writer::{
    OwnerContigInput, OwnerKeyRange, OwnerMetadata, OwnerMetagenomeInput, OwnerWriteInput,
    publish_owner_manifest, write_owner,
};
use crate::trace::{MetagenomeTrace, TraceConfig, TraceEngine, TraceResult};
use crate::writer::{BuildConfig, build};
use noodles_bgzf::{self as bgzf, gzi};
use serde_json::json;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

pub(crate) struct Fixture {
    _directory: Option<tempfile::TempDir>,
    pub(crate) root: PathBuf,
    pub(crate) query: Vec<u8>,
    jam: PathBuf,
    pub(crate) jidx: PathBuf,
    source_manifest: PathBuf,
    pub(crate) owner_manifest: PathBuf,
    auxiliary: Vec<PathBuf>,
}
#[test]
fn owner_generation_preserves_evidence_trace_and_standalone_deployment() {
    let fixture = build_fixture();
    let old_reader = JidxReader::open(&fixture.jidx).unwrap();
    let owner_reader = OwnerReader::open(&fixture.owner_manifest).unwrap();
    owner_reader.verify_checksum().unwrap();
    let old_evidence = old_evidence(&old_reader);
    let owner_evidence = owner_evidence(&owner_reader, &old_evidence);
    assert_eq!(owner_evidence, old_evidence);
    assert!(old_evidence.iter().any(|(_, members)| members.len() > 1));
    assert!(
        old_evidence
            .iter()
            .any(|(_, members)| { members.iter().any(|(_, occurrences)| occurrences.len() > 1) })
    );
    let absent = (0..1u64 << 30)
        .find(|key| old_reader.find_seed(*key).unwrap().is_none())
        .unwrap();
    assert_eq!(owner_reader.find_seeds_batch(&[absent]).unwrap(), [None]);
    let config = TraceConfig {
        use_sketch: false,
        min_seed_hits: 2,
        diagonal_bin_bases: 32,
        flank_bases: 48,
        min_identity: 0.95,
        min_aligned_bases: 128,
        endpoint_bases: 48,
        circular: false,
        ..TraceConfig::default()
    };
    let old = TraceEngine::open(&fixture.jam, &fixture.jidx, &fixture.source_manifest, None)
        .unwrap()
        .search("synthetic-plasmid", &fixture.query, config)
        .unwrap();
    assert!(old.metagenomes.iter().any(|entry| {
        entry.name == "doc-b"
            && fragments(entry)
                .iter()
                .any(|fragment| fragment.alignment.strand == Strand::Reverse)
    }));
    assert!(old.metagenomes.iter().all(|entry| entry.name != "doc-c"));
    drop(owner_reader);
    drop(old_reader);
    let unavailable = fixture.root.join("legacy-unavailable");
    std::fs::create_dir(&unavailable).unwrap();
    for path in std::iter::once(&fixture.jam)
        .chain(std::iter::once(&fixture.jidx))
        .chain(std::iter::once(&fixture.source_manifest))
        .chain(fixture.auxiliary.iter())
    {
        std::fs::rename(path, unavailable.join(path.file_name().unwrap())).unwrap();
    }
    let owner = TraceEngine::open_owner(&fixture.owner_manifest, None).unwrap();
    owner.verify_index().unwrap();
    let actual = owner
        .search("synthetic-plasmid", &fixture.query, config)
        .unwrap();
    validate_cigars(&actual, &fixture.query, &fixture.owner_manifest, config);
    assert_biological_equal(&old, &actual);
}
#[test]
fn owner_header_supports_more_than_two_billion_contigs_without_allocation() {
    let contig_count = 3_000_000_000u32;
    let lengths = [
        0,
        u64::from(OWNER_DOCUMENT_SIZE),
        u64::from(contig_count) * u64::from(OWNER_CONTIG_SIZE),
        0,
        0,
        0,
        0,
    ];
    let mut offset = OWNER_HEADER_SIZE as u64;
    let mut sections = Vec::with_capacity(OWNER_SECTION_COUNT);
    for (kind, length) in OwnerSection::ALL.into_iter().take(7).zip(lengths) {
        offset = offset.next_multiple_of(OWNER_PAGE_SIZE);
        sections.push(OwnerSectionDescriptor {
            kind,
            record_size: kind.record_size(),
            offset,
            length,
        });
        offset = offset.checked_add(length).unwrap();
    }
    offset = offset.next_multiple_of(OWNER_PAGE_SIZE);
    let checksum_length = crate::owner_format::checksum_layout(offset / OWNER_PAGE_SIZE - 1)
        .unwrap()
        .iter()
        .map(|level| level.page_count * OWNER_PAGE_SIZE)
        .sum();
    sections.push(OwnerSectionDescriptor {
        kind: OwnerSection::PageChecksums,
        record_size: OwnerSection::PageChecksums.record_size(),
        offset,
        length: checksum_length,
    });
    let header = OwnerHeader {
        flags: HAS_METADATA,
        k: 15,
        rescue_k15: false,
        minimizer_window: 4,
        owner_ordinal: 0,
        owner_count: 1,
        first_key: 0,
        last_key: u64::MAX,
        key_count: 0,
        occurrence_count: 0,
        document_count: 1,
        contig_count,
        generation_id: [1; 32],
        body_sha256: [2; 32],
        checksum_root_sha256: [3; 32],
        metadata_sha256: [4; 32],
        sections: sections.try_into().unwrap(),
    };
    let bytes = header.encode().unwrap();
    let decoded = OwnerHeader::decode(&bytes, offset + checksum_length).unwrap();
    assert_eq!(decoded.contig_count, contig_count);
}

#[test]
fn sparse_metadata_preserves_high_original_contig_ids_and_extent_boundaries() {
    let fixture = build_fixture_at(None);
    let source = JidxReader::open(&fixture.jidx).unwrap();
    let source_metadata = owner_metadata(&source);
    let source_keys = owner_keys(&source);
    let source_member = source_keys
        .iter()
        .flat_map(|key| key.members.iter().map(move |member| (key.key, member)))
        .find(|(_, member)| {
            member.document_id == 1
                && member
                    .occurrences
                    .iter()
                    .any(|entry| entry.local_contig == 0)
        })
        .unwrap();
    let occurrence = *source_member
        .1
        .occurrences
        .iter()
        .find(|entry| entry.local_contig == 0)
        .unwrap();
    let mut empty = source_metadata.metagenomes[0].clone();
    empty.name = "empty-low-extent".to_string();
    empty.original_contig_start = 0;
    empty.original_contig_count = 4_000_000_000;
    empty.locus_bits = 1;
    empty.contigs.clear();
    let mut high = source_metadata.metagenomes[1].clone();
    high.name = "selected-high-extent".to_string();
    high.original_contig_start = 4_000_000_000;
    high.original_contig_count = 1;
    high.contigs.truncate(1);
    high.contigs[0].local_contig = 0;
    let metadata = OwnerMetadata {
        metagenomes: vec![empty, high],
    };
    let keys = vec![OwnerKey {
        key: source_member.0,
        members: vec![OwnerMember {
            document_id: 1,
            occurrences: vec![occurrence],
        }],
    }];
    let owner = fixture.root.join("sparse-high.jowner");
    write_owner(
        &owner,
        OwnerWriteInput {
            owner_ordinal: 0,
            owner_count: 1,
            range: OwnerKeyRange {
                first: 0,
                last: u64::MAX,
                complete: true,
            },
            k: source.header().k,
            rescue_k15: source.header().rescue_k15,
            minimizer_window: source.header().minimizer_window,
            generation_id: [13; 32],
            document_count: 2,
            contig_count: 4_000_000_001,
            keys: &keys,
            loci: &metadata,
            metadata: Some(&metadata),
        },
    )
    .unwrap();
    let manifest = fixture.root.join("sparse-high.json");
    publish_owner_manifest(&manifest, &[owner], 0, true).unwrap();
    let reader = OwnerReader::open(manifest).unwrap();
    reader.verify_checksum().unwrap();
    assert_eq!(
        reader.metagenome(0).unwrap().unwrap().contig_count,
        4_000_000_000
    );
    assert!(reader.contig(0).unwrap().is_none());
    let contig = reader.contig(4_000_000_000).unwrap().unwrap();
    assert_eq!(contig.metagenome_id, 1);
    let seed = reader.find_seeds_batch(&[keys[0].key]).unwrap()[0].unwrap();
    let document = reader.seed_documents(seed).unwrap()[0];
    let decoded = reader.seed_document_occurrences(seed, document).unwrap();
    assert_eq!(decoded[0].contig_id, 4_000_000_000);
    assert_eq!(decoded[0].position, occurrence.position);
}
fn build_fixture() -> Fixture {
    build_fixture_at(None)
}

pub(crate) fn build_fixture_at(output: Option<PathBuf>) -> Fixture {
    let temp_root = std::env::var_os("TMPDIR")
        .map(PathBuf::from)
        .expect("TMPDIR must be explicit for owner fixtures");
    let directory = output.is_none().then(|| {
        tempfile::Builder::new()
            .prefix("jam-owner-trace-")
            .tempdir_in(temp_root)
            .unwrap()
    });
    let root = output.unwrap_or_else(|| directory.as_ref().unwrap().path().to_path_buf());
    if directory.is_none() {
        std::fs::create_dir(&root).unwrap();
    }
    let query = random_sequence(17, 224);
    let reverse = reverse_complement(&query);
    let targets = [
        (
            "doc-a",
            [
                random_sequence(31, 71),
                query.clone(),
                random_sequence(37, 83),
                query.clone(),
                random_sequence(41, 67),
            ]
            .concat(),
        ),
        (
            "doc-b",
            [random_sequence(43, 97), reverse, random_sequence(47, 109)].concat(),
        ),
        ("doc-c", random_sequence(10_007, 640)),
    ];
    let mut fasta_inputs = Vec::new();
    let mut manifest_entries = Vec::new();
    let mut auxiliary = Vec::new();
    for (name, sequence) in &targets {
        let fasta = root.join(format!("{name}.fa"));
        write_fasta(&fasta, name, sequence, 64);
        fasta_inputs.push(fasta);
        let (bgzf, fai, gzi) = write_bgzf(&root, name, sequence);
        manifest_entries.push(json!({
            "name": name,
            "bgzf": bgzf,
            "fai": fai,
            "gzi": gzi,
        }));
        auxiliary.extend([fai, gzi]);
    }
    let jam = root.join("synthetic.jam");
    build(
        &fasta_inputs,
        &jam,
        &BuildConfig {
            kmer_size: 15,
            fscale: 1,
            singleton: true,
            memory: 1,
            ..BuildConfig::default()
        },
    )
    .unwrap();
    auxiliary.extend(fasta_inputs);
    let source_manifest = root.join("source.json");
    std::fs::write(
        &source_manifest,
        serde_json::to_vec(&json!({"metagenomes": manifest_entries})).unwrap(),
    )
    .unwrap();
    let jidx = root.join("synthetic.jidx");
    build_local_jidx(
        &jam,
        &source_manifest,
        &jidx,
        JidxBuildConfig {
            k: 15,
            minimizer_window: 4,
            rescue_k15: false,
        },
    )
    .unwrap();
    let source = JidxReader::open(&jidx).unwrap();
    let keys = owner_keys(&source);
    let split = keys.len() / 2;
    let boundary = keys[split - 1].key;
    let metadata = owner_metadata(&source);
    let generation_id = [7; 32];
    let owner_paths = [root.join("owner-000.jowner"), root.join("owner-001.jowner")];
    let ranges = [
        OwnerKeyRange {
            first: 0,
            last: boundary,
            complete: true,
        },
        OwnerKeyRange {
            first: boundary + 1,
            last: u64::MAX,
            complete: true,
        },
    ];
    let mut owner_entries = Vec::new();
    for ordinal in 0..2 {
        let owned_keys = if ordinal == 0 {
            &keys[..split]
        } else {
            &keys[split..]
        };
        let stats = write_owner(
            &owner_paths[ordinal],
            OwnerWriteInput {
                owner_ordinal: ordinal as u32,
                owner_count: 2,
                range: ranges[ordinal],
                k: 15,
                rescue_k15: false,
                minimizer_window: 4,
                generation_id,
                document_count: source.header().document_count,
                contig_count: source.header().contig_count,
                keys: owned_keys,
                loci: &metadata,
                metadata: (ordinal == 0).then_some(&metadata),
            },
        )
        .unwrap();
        owner_entries.push(json!({
            "path": owner_paths[ordinal].file_name().unwrap().to_str().unwrap(),
            "header_sha256": hex(&stats.header_sha256),
        }));
    }
    let owner_manifest = root.join("owners.json");
    std::fs::write(
        &owner_manifest,
        serde_json::to_vec(&json!({
            "version": 1,
            "complete": true,
            "metadata_owner": 0,
            "owners": owner_entries,
        }))
        .unwrap(),
    )
    .unwrap();
    drop(source);
    Fixture {
        _directory: directory,
        root,
        query,
        jam,
        jidx,
        source_manifest,
        owner_manifest,
        auxiliary,
    }
}

pub(crate) fn owner_keys(source: &JidxReader) -> Vec<OwnerKey> {
    (0..source.header().seed_count)
        .map(|ordinal| {
            let seed = crate::jidx_postings::entry(source, ordinal).unwrap();
            let members = source
                .seed_documents(seed)
                .unwrap()
                .into_iter()
                .map(|document| {
                    let metagenome = source.metagenome(document.metagenome_id).unwrap().unwrap();
                    let occurrences = source
                        .seed_document_occurrences(seed, document)
                        .unwrap()
                        .into_iter()
                        .map(|occurrence| OwnerOccurrence {
                            local_contig: occurrence.contig_id - metagenome.contig_start,
                            position: occurrence.position,
                            canonical_orientation: occurrence.canonical_orientation,
                        })
                        .collect();
                    OwnerMember {
                        document_id: document.metagenome_id,
                        occurrences,
                    }
                })
                .collect();
            OwnerKey {
                key: seed.packed_key,
                members,
            }
        })
        .collect()
}

pub(crate) fn owner_metadata(source: &JidxReader) -> OwnerMetadata {
    OwnerMetadata {
        metagenomes: (0..source.header().document_count)
            .map(|document_id| {
                let metagenome = source.metagenome(document_id).unwrap().unwrap();
                let contigs = (metagenome.contig_start
                    ..metagenome.contig_start + metagenome.contig_count)
                    .map(|contig_id| {
                        let contig = source.contig(contig_id).unwrap().unwrap();
                        OwnerContigInput {
                            local_contig: contig_id - metagenome.contig_start,
                            name: contig.name.to_string(),
                            length: contig.length,
                            fasta_offset: contig.fasta_offset,
                            line_bases: contig.line_bases,
                            line_width: contig.line_width,
                        }
                    })
                    .collect::<Vec<_>>();
                let maximum_byte = contigs
                    .iter()
                    .map(|contig| {
                        let position = contig.length - 1;
                        contig.fasta_offset
                            + position / u64::from(contig.line_bases) * u64::from(contig.line_width)
                            + position % u64::from(contig.line_bases)
                    })
                    .max()
                    .unwrap();
                OwnerMetagenomeInput {
                    name: metagenome.name.to_string(),
                    bgzf_uri: metagenome.bgzf_uri.to_string(),
                    bgzf_bytes: metagenome.bgzf_bytes,
                    bgzf_sha256: metagenome.bgzf_sha256,
                    gzi: metagenome.gzi.to_vec(),
                    original_contig_start: metagenome.contig_start,
                    original_contig_count: metagenome.contig_count,
                    locus_bits: (64 - maximum_byte.leading_zeros()) as u8,
                    contigs,
                }
            })
            .collect(),
    }
}

type Evidence = Vec<(u64, Vec<(u32, Vec<SeedOccurrence>)>)>;

fn old_evidence(reader: &JidxReader) -> Evidence {
    (0..reader.header().seed_count)
        .map(|ordinal| {
            let seed = crate::jidx_postings::entry(reader, ordinal).unwrap();
            let members = reader
                .seed_documents(seed)
                .unwrap()
                .into_iter()
                .map(|document| {
                    (
                        document.metagenome_id,
                        reader.seed_document_occurrences(seed, document).unwrap(),
                    )
                })
                .collect();
            (seed.packed_key, members)
        })
        .collect()
}

fn owner_evidence(reader: &OwnerReader, reference: &Evidence) -> Evidence {
    reference
        .iter()
        .map(|(key, _)| {
            let seed = reader.find_seeds_batch(&[*key]).unwrap()[0].unwrap();
            let members = reader
                .seed_documents(seed)
                .unwrap()
                .into_iter()
                .map(|document| {
                    (
                        document.metagenome_id,
                        reader.seed_document_occurrences(seed, document).unwrap(),
                    )
                })
                .collect();
            (*key, members)
        })
        .collect()
}

fn assert_biological_equal(expected: &TraceResult, actual: &TraceResult) {
    assert_eq!(actual.completion, expected.completion);
    assert_eq!(actual.metagenomes.len(), expected.metagenomes.len());
    for (actual, expected) in actual.metagenomes.iter().zip(&expected.metagenomes) {
        assert_eq!(actual.name, expected.name);
        assert_eq!(actual.shared_hashes, expected.shared_hashes);
        assert_eq!(actual.containment, expected.containment);
        assert_eq!(actual.exact_seed_hits, expected.exact_seed_hits);
        assert_eq!(actual.contigs, expected.contigs);
        assert_eq!(actual.mosaic, expected.mosaic);
    }
}

fn validate_cigars(result: &TraceResult, query: &[u8], manifest: &Path, config: TraceConfig) {
    use crate::alignment::EditOperation;
    let index = OwnerReader::open(manifest).unwrap();
    for metagenome in &result.metagenomes {
        let source = index.metagenome(metagenome.metagenome_id).unwrap().unwrap();
        let mut sequence = crate::bgzf::BgzfReader::open(source, None, true).unwrap();
        for fragment in fragments(metagenome) {
            let alignment = &fragment.alignment;
            alignment.validate_cigar().unwrap();
            let query_bases = fragment
                .query_segments
                .iter()
                .flat_map(|interval| &query[interval.start as usize..interval.end as usize])
                .copied()
                .collect::<Vec<_>>();
            let contig = index.contig(fragment.contig_id).unwrap().unwrap();
            let target = sequence
                .read_contig_range(
                    contig,
                    alignment.target_interval.start,
                    alignment.target_interval.end,
                )
                .unwrap();
            let target = if alignment.strand == Strand::Reverse {
                reverse_complement(&target)
            } else {
                target
            };
            let (mut q, mut t, mut score) = (0usize, 0usize, 0i64);
            for run in &alignment.edit_script {
                let length = run.length as usize;
                match run.operation {
                    EditOperation::Equal | EditOperation::Substitution => {
                        for (left, right) in query_bases[q..q + length]
                            .iter()
                            .zip(&target[t..t + length])
                        {
                            assert_eq!(
                                left.eq_ignore_ascii_case(right),
                                run.operation == EditOperation::Equal
                            );
                        }
                        score += i64::from(if run.operation == EditOperation::Equal {
                            config.alignment.match_score
                        } else {
                            config.alignment.mismatch_score
                        }) * length as i64;
                        q += length;
                        t += length;
                    }
                    EditOperation::Insertion => {
                        t += length;
                        score += i64::from(config.alignment.gap_open_score)
                            + i64::from(config.alignment.gap_extend_score) * length as i64;
                    }
                    EditOperation::Deletion => {
                        q += length;
                        score += i64::from(config.alignment.gap_open_score)
                            + i64::from(config.alignment.gap_extend_score) * length as i64;
                    }
                }
            }
            assert_eq!((q, t), (query_bases.len(), target.len()));
            assert_eq!(score, i64::from(alignment.score));
        }
    }
}

fn fragments(trace: &MetagenomeTrace) -> Vec<&crate::mosaic::Fragment> {
    trace
        .mosaic
        .primary
        .iter()
        .map(|selected| &selected.fragment)
        .chain(trace.mosaic.alternatives.iter())
        .collect()
}

fn write_bgzf(directory: &Path, name: &str, sequence: &[u8]) -> (PathBuf, PathBuf, PathBuf) {
    let contig = format!("{name}-contig");
    let bgzf_path = directory.join(format!("{name}.bgz"));
    let fai_path = directory.join(format!("{name}.bgz.fai"));
    let gzi_path = directory.join(format!("{name}.bgz.gzi"));
    let mut raw = format!(">{contig}\n").into_bytes();
    let offset = raw.len();
    for line in sequence.chunks(64) {
        raw.extend_from_slice(line);
        raw.push(b'\n');
    }
    let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
    writer.write_all(&raw).unwrap();
    writer.finish().unwrap();
    std::fs::write(
        &fai_path,
        format!("{contig}\t{}\t{offset}\t64\t65\n", sequence.len()),
    )
    .unwrap();
    gzi::fs::write(&gzi_path, &gzi::Index::default()).unwrap();
    (bgzf_path, fai_path, gzi_path)
}

fn write_fasta(path: &Path, name: &str, sequence: &[u8], width: usize) {
    let mut file = File::create(path).unwrap();
    writeln!(file, ">{name}").unwrap();
    for line in sequence.chunks(width) {
        file.write_all(line).unwrap();
        file.write_all(b"\n").unwrap();
    }
}

pub(crate) fn random_sequence(mut state: u64, length: usize) -> Vec<u8> {
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            b"ACGT"[(state >> 62) as usize]
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
            _ => unreachable!(),
        })
        .collect()
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}
