use crate::jidx::{sha256, sha256_reader};
use crate::jidx_reader::JidxReader;
use crate::owner_format::{
    HAS_METADATA, OWNER_HEADER_SIZE, OWNER_PAGE_SIZE, OwnerHeader, OwnerReaderError, OwnerSection,
    put_u32,
};
use crate::owner_reader::OwnerReader;
use crate::owner_tests::{build_fixture_at, owner_keys, owner_metadata, random_sequence};
use crate::owner_writer::{
    MAX_PROTOTYPE_ENCODED_BYTES, OwnerKeyRange, OwnerWriteInput, publish_owner_manifest,
    validate_prototype_payload_bytes, write_owner,
};
use crate::trace::{TraceConfig, TraceEngine, TraceResult};
use serde_json::Value;
use std::fs::{File, OpenOptions};
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

fn owner_paths(root: &Path) -> [PathBuf; 2] {
    [root.join("owner-000.jowner"), root.join("owner-001.jowner")]
}

fn read_header(path: &Path) -> OwnerHeader {
    let mut file = File::open(path).unwrap();
    let length = file.metadata().unwrap().len();
    let mut bytes = [0; OWNER_HEADER_SIZE];
    file.read_exact(&mut bytes).unwrap();
    OwnerHeader::decode(&bytes, length).unwrap()
}

fn write_header(path: &Path, header: &OwnerHeader) {
    let mut file = OpenOptions::new().write(true).open(path).unwrap();
    file.write_all(&header.encode().unwrap()).unwrap();
    file.sync_all().unwrap();
}

fn update_manifest_header(manifest: &Path, owner: usize, header: &[u8]) {
    let mut value: Value = serde_json::from_slice(&std::fs::read(manifest).unwrap()).unwrap();
    value["owners"][owner]["header_sha256"] = Value::String(hex(&sha256(header)));
    std::fs::write(manifest, serde_json::to_vec(&value).unwrap()).unwrap();
}

fn first_stored_key(path: &Path) -> u64 {
    let header = read_header(path);
    let directory = header.section(OwnerSection::BlockDirectory);
    let mut file = File::open(path).unwrap();
    file.seek(SeekFrom::Start(directory.offset)).unwrap();
    let mut bytes = [0; 8];
    file.read_exact(&mut bytes).unwrap();
    u64::from_le_bytes(bytes)
}

#[test]
fn coordinated_data_and_leaf_corruption_fails_trusted_root() {
    let fixture = build_fixture_at(None);
    let path = owner_paths(&fixture.root)[0].clone();
    let header = read_header(&path);
    let checksums = header.section(OwnerSection::PageChecksums);
    let page = header.section(OwnerSection::HotPostings).offset / OWNER_PAGE_SIZE;
    assert!(page > 0);
    let mut bytes = std::fs::read(&path).unwrap();
    let start = usize::try_from(page * OWNER_PAGE_SIZE).unwrap();
    bytes[start] ^= 1;
    let digest = sha256(&bytes[start..start + OWNER_PAGE_SIZE as usize]);
    let leaf = usize::try_from(checksums.offset + (page - 1) * 32).unwrap();
    bytes[leaf..leaf + 32].copy_from_slice(&digest);
    std::fs::write(&path, bytes).unwrap();
    assert!(matches!(
        OwnerReader::open(&fixture.owner_manifest),
        Err(OwnerReaderError::ChecksumMismatch)
    ));
}

#[test]
fn lazy_checksums_authenticate_deep_leaf_pages() {
    let fixture = build_fixture_at(None);
    let source = JidxReader::open(&fixture.jidx).unwrap();
    let keys = owner_keys(&source);
    let mut metadata = owner_metadata(&source);
    metadata.metagenomes[0].name = "x".repeat(130 * OWNER_PAGE_SIZE as usize);
    let path = fixture.root.join("deep.jowner");
    write_owner(
        &path,
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
            generation_id: [11; 32],
            document_count: source.header().document_count,
            contig_count: source.header().contig_count,
            keys: &keys,
            loci: &metadata,
            metadata: Some(&metadata),
        },
    )
    .unwrap();
    let manifest = fixture.root.join("deep.json");
    publish_owner_manifest(&manifest, std::slice::from_ref(&path), 0, true).unwrap();
    let header = read_header(&path);
    let checksums = header.section(OwnerSection::PageChecksums);
    assert!(
        crate::owner_format::checksum_layout(checksums.offset / OWNER_PAGE_SIZE - 1)
            .unwrap()
            .len()
            > 1
    );
    let page = header.section(OwnerSection::HotPostings).offset / OWNER_PAGE_SIZE;
    let mut bytes = std::fs::read(&path).unwrap();
    let start = (page * OWNER_PAGE_SIZE) as usize;
    bytes[start] ^= 1;
    let digest = sha256(&bytes[start..start + OWNER_PAGE_SIZE as usize]);
    let leaf = (checksums.offset + (page - 1) * 32) as usize;
    bytes[leaf..leaf + 32].copy_from_slice(&digest);
    std::fs::write(&path, bytes).unwrap();
    assert!(matches!(
        OwnerReader::open(manifest),
        Err(OwnerReaderError::ChecksumMismatch)
    ));
}

#[test]
fn cached_owner_rejects_in_place_file_mutation() {
    let fixture = build_fixture_at(None);
    let path = owner_paths(&fixture.root)[0].clone();
    let key = first_stored_key(&path);
    let reader = OwnerReader::open(&fixture.owner_manifest).unwrap();
    let seed = reader.find_seeds_batch(&[key]).unwrap()[0].unwrap();
    let document = reader.seed_documents(seed).unwrap()[0];
    reader.seed_document_occurrences(seed, document).unwrap();
    assert!(reader.find_seeds_batch(&[0]).unwrap()[0].is_none());
    let offset = read_header(&path).section(OwnerSection::Contigs).offset;
    let mut byte = [0];
    let mut input = File::open(owner_paths(&fixture.root)[0].clone()).unwrap();
    input.seek(SeekFrom::Start(offset)).unwrap();
    input.read_exact(&mut byte).unwrap();
    byte[0] ^= 1;
    let mut file = OpenOptions::new().write(true).open(path).unwrap();
    file.seek(SeekFrom::Start(offset)).unwrap();
    file.write_all(&byte).unwrap();
    file.sync_all().unwrap();
    assert!(reader.find_seeds_batch(&[key]).is_err());
    assert!(reader.find_seeds_batch(&[0]).is_err());
    assert!(reader.seed_document_occurrences(seed, document).is_err());
}

#[test]
fn complete_manifest_rejects_gap_and_overlap() {
    for overlap in [false, true] {
        let fixture = build_fixture_at(None);
        let owners = owner_paths(&fixture.root);
        let first = read_header(&owners[0]);
        let mut second = read_header(&owners[1]);
        second.first_key = if overlap {
            first.last_key
        } else {
            first.last_key.checked_add(2).unwrap()
        };
        write_header(&owners[1], &second);
        let manifest = fixture
            .root
            .join(if overlap { "overlap.json" } else { "gap.json" });
        assert!(publish_owner_manifest(&manifest, &owners, 0, true).is_err());
        assert!(!manifest.exists());
    }
}

#[test]
fn incomplete_range_never_turns_unconverted_keys_into_absence() {
    let fixture = build_fixture_at(None);
    let owners = owner_paths(&fixture.root);
    let mut last = read_header(&owners[1]);
    last.last_key = u64::MAX - 1;
    write_header(&owners[1], &last);
    let manifest = fixture.root.join("partial.json");
    publish_owner_manifest(&manifest, &owners, 0, false).unwrap();
    let reader = OwnerReader::open(manifest).unwrap();
    assert!(matches!(
        reader.find_seeds_batch(&[u64::MAX]),
        Err(OwnerReaderError::KeyNotCovered(u64::MAX))
    ));
}

#[test]
fn missing_metadata_owner_flag_and_truncation_are_rejected() {
    let fixture = build_fixture_at(None);
    let owners = owner_paths(&fixture.root);
    let mut bytes = std::fs::read(&owners[0]).unwrap();
    let flags = u32::from_le_bytes(bytes[12..16].try_into().unwrap()) & !HAS_METADATA;
    put_u32(&mut bytes, 12, flags);
    std::fs::write(&owners[0], &bytes).unwrap();
    update_manifest_header(&fixture.owner_manifest, 0, &bytes[..OWNER_HEADER_SIZE]);
    assert!(OwnerReader::open(&fixture.owner_manifest).is_err());

    let fixture = build_fixture_at(None);
    let path = owner_paths(&fixture.root)[1].clone();
    let length = std::fs::metadata(&path).unwrap().len();
    OpenOptions::new()
        .write(true)
        .open(path)
        .unwrap()
        .set_len(length - 1)
        .unwrap();
    assert!(OwnerReader::open(&fixture.owner_manifest).is_err());
}

#[test]
fn owners_with_different_locus_metadata_cannot_share_a_generation() {
    let fixture = build_fixture_at(None);
    let source = JidxReader::open(&fixture.jidx).unwrap();
    let keys = owner_keys(&source);
    let split = keys.len() / 2;
    let mut loci = owner_metadata(&source);
    let contig = loci
        .metagenomes
        .last_mut()
        .unwrap()
        .contigs
        .last_mut()
        .unwrap();
    contig.fasta_offset += 1;
    let mismatched = fixture.root.join("owner-mismatched-loci.jowner");
    write_owner(
        &mismatched,
        OwnerWriteInput {
            owner_ordinal: 1,
            owner_count: 2,
            range: OwnerKeyRange {
                first: keys[split - 1].key + 1,
                last: u64::MAX,
                complete: true,
            },
            k: source.header().k,
            rescue_k15: source.header().rescue_k15,
            minimizer_window: source.header().minimizer_window,
            generation_id: [7; 32],
            document_count: source.header().document_count,
            contig_count: source.header().contig_count,
            keys: &keys[split..],
            loci: &loci,
            metadata: None,
        },
    )
    .unwrap();
    let manifest = fixture.root.join("mismatched-loci.json");
    assert!(
        publish_owner_manifest(
            &manifest,
            &[owner_paths(&fixture.root)[0].clone(), mismatched],
            0,
            true,
        )
        .is_err()
    );
    assert!(!manifest.exists());
}

#[test]
fn invalid_occurrence_contig_and_position_are_rejected() {
    for invalid_position in [false, true] {
        let fixture = build_fixture_at(None);
        let source = JidxReader::open(&fixture.jidx).unwrap();
        let mut keys = owner_keys(&source);
        let metadata = owner_metadata(&source);
        let occurrence = &mut keys[0].members[0].occurrences[0];
        if invalid_position {
            occurrence.position = u64::MAX;
        } else {
            occurrence.local_contig = u32::MAX;
        }
        let owner = fixture.root.join(if invalid_position {
            "invalid-position.jowner"
        } else {
            "invalid-contig.jowner"
        });
        assert!(
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
                    generation_id: [9; 32],
                    document_count: source.header().document_count,
                    contig_count: source.header().contig_count,
                    keys: &keys,
                    loci: &metadata,
                    metadata: Some(&metadata),
                },
            )
            .is_err()
        );
        assert!(!owner.exists());
    }
}

#[test]
fn manifest_publication_is_atomic_and_refuses_replacement() {
    let fixture = build_fixture_at(None);
    let owners = owner_paths(&fixture.root);
    let manifest = fixture.root.join("published.json");
    let first = publish_owner_manifest(&manifest, &owners, 0, true).unwrap();
    assert_eq!(
        first,
        sha256_reader(File::open(&manifest).unwrap()).unwrap()
    );
    assert!(publish_owner_manifest(&manifest, &owners, 0, true).is_err());
}

#[test]
fn prototype_encoded_byte_budget_checks_boundaries_without_allocation() {
    assert_eq!(
        validate_prototype_payload_bytes(&[MAX_PROTOTYPE_ENCODED_BYTES]).unwrap(),
        MAX_PROTOTYPE_ENCODED_BYTES
    );
    assert!(validate_prototype_payload_bytes(&[MAX_PROTOTYPE_ENCODED_BYTES, 1]).is_err());
    assert!(validate_prototype_payload_bytes(&[u64::MAX, 1]).is_err());
}

#[test]
fn owner_shared_batch_matches_sequential_low_reuse_and_repeated_queries() {
    fn io(results: &[TraceResult]) -> (u64, u64, u64) {
        results
            .iter()
            .flat_map(|result| &result.metagenomes)
            .fold((0, 0, 0), |totals, trace| {
                (
                    totals.0 + trace.compressed_bytes_read,
                    totals.1 + trace.range_requests,
                    totals.2 + trace.bgzf_blocks_decoded,
                )
            })
    }
    fn normalize_io(results: &mut [TraceResult]) {
        for trace in results
            .iter_mut()
            .flat_map(|result| &mut result.metagenomes)
        {
            trace.compressed_bytes_read = 0;
            trace.range_requests = 0;
            trace.bgzf_blocks_decoded = 0;
        }
    }

    let fixture = build_fixture_at(None);
    let engine = TraceEngine::open_owner(&fixture.owner_manifest, None).unwrap();
    let queries = vec![
        ("repeat".to_string(), fixture.query.clone()),
        ("repeat".to_string(), fixture.query.clone()),
        (
            "low-a".to_string(),
            random_sequence(800_021, fixture.query.len()),
        ),
        (
            "low-b".to_string(),
            random_sequence(800_023, fixture.query.len()),
        ),
    ];
    let config = TraceConfig {
        use_sketch: false,
        circular: false,
        ..TraceConfig::default()
    };
    let mut expected = queries
        .iter()
        .map(|(id, sequence)| engine.search(id, sequence, config).unwrap())
        .collect::<Vec<_>>();
    let mut actual = engine.search_batch(&queries, config).unwrap();
    let expected_total_io = io(&expected);
    let actual_total_io = io(&actual);
    assert!(actual_total_io.0 <= expected_total_io.0);
    assert!(actual_total_io.1 <= expected_total_io.1);
    assert!(actual_total_io.2 <= expected_total_io.2);
    let expected_io = io(&expected[..2]);
    let actual_io = io(&actual[..2]);
    assert!(actual_io.0 < expected_io.0);
    assert!(actual_io.1 < expected_io.1);
    assert!(actual_io.2 < expected_io.2);
    normalize_io(&mut expected);
    normalize_io(&mut actual);
    assert_eq!(actual, expected);
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}
