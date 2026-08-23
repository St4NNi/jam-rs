#[path = "support/remote/mod.rs"]
mod remote;

use jam_rs::jma::ArchiveReader;
use jam_rs::jma::SequenceRange;
use jam_rs::jma::builder::{ArchiveBuildConfig, write_archive_from_fasta};
use jam_rs::jma::format::SECTION_SEQUENCES;
use jam_rs::jma::index::{bucket_for, build_index, encode_index};
use jam_rs::jma::reader::JmaReader;
use jam_rs::resource::ResourceOpenOptions;
use jam_rs::resource::object::ObjectResource;
use jam_rs::trace::config::SeedSensitivity;
use jam_rs::trace::raw::open_resource;
use remote::HttpFixture;
use std::fs;

struct ArchiveFixture {
    archive: Vec<u8>,
    index: Vec<u8>,
    sequence: Vec<u8>,
}

fn archive_fixture() -> ArchiveFixture {
    let directory = tempfile::Builder::new()
        .prefix("trace-jma-range-")
        .tempdir_in("target")
        .expect("temporary JMA fixture directory");
    let sequence: Vec<u8> = (0..2_048)
        .map(|index| b"ACGT"[(index * 37 + index / 7) % 4])
        .collect();
    let input = directory.path().join("assembly.fa");
    let output = directory.path().join("assembly.jma");
    fs::write(
        &input,
        format!(
            ">candidate-contig\n{}\n",
            String::from_utf8_lossy(&sequence)
        ),
    )
    .expect("write JMA input");
    write_archive_from_fasta(
        &input,
        &output,
        ArchiveBuildConfig {
            block_bases: 32,
            k31_scale: 1,
            k21_scale: Some(1),
            ..ArchiveBuildConfig::default()
        },
    )
    .expect("build JMA fixture");
    let archive = fs::read(output).expect("read JMA fixture");
    let index = encode_index(&build_index(&archive).expect("build JMA sidecar")).unwrap();
    ArchiveFixture {
        archive,
        index,
        sequence,
    }
}

fn archive_options() -> ResourceOpenOptions {
    ResourceOpenOptions {
        // Two-byte blocks align the fixed JMA offsets while keeping this test
        // small enough to inspect every request without fetching payloads.
        cache_block_bytes: 2,
        max_cache_bytes: 4 * 1024 * 1024,
        max_retries: 0,
        ..ResourceOpenOptions::default()
    }
}

fn index_options() -> ResourceOpenOptions {
    ResourceOpenOptions {
        cache_block_bytes: 4096,
        max_cache_bytes: 4 * 1024 * 1024,
        max_retries: 0,
        ..ResourceOpenOptions::default()
    }
}

fn request_ranges_cover_exact(records: &[remote::RequestRecord], offset: u64, length: u64) -> bool {
    let Some(expected_end) = offset.checked_add(length) else {
        return false;
    };
    let mut ranges = records
        .iter()
        .filter(|record| record.method == "GET")
        .filter_map(|record| {
            let (start, end) = record.range?;
            let end = u64::try_from(end).ok()?.checked_add(1)?;
            Some((u64::try_from(start).ok()?, end))
        })
        .collect::<Vec<_>>();
    if ranges.is_empty() {
        return false;
    }
    ranges.sort_unstable();
    let mut covered_end = offset;
    for (start, end) in ranges {
        if start < offset || end > expected_end || start > covered_end {
            return false;
        }
        covered_end = covered_end.max(end);
    }
    covered_end == expected_end
}

#[test]
fn indexed_http_reader_fetches_only_header_directory_index_buckets_and_blocks() {
    let fixture = archive_fixture();
    let archive_index = build_index(&fixture.archive).unwrap();
    let archive_server = HttpFixture::new(fixture.archive.clone());
    let index_server = HttpFixture::new(fixture.index.clone());
    let archive_url = archive_server.url("/candidate.jma?token=archive-secret");
    let index_url = index_server.url("/candidate.jma.idx.json?token=index-secret");
    let reader = JmaReader::open_indexed(
        ObjectResource::open(archive_url, archive_options()).unwrap(),
        ObjectResource::open(index_url, index_options()).unwrap(),
    )
    .expect("open indexed remote JMA");

    let sequence_section = reader
        .sections()
        .iter()
        .find(|section| section.kind == SECTION_SEQUENCES)
        .expect("sequence section descriptor");
    let seed_sections = reader
        .sections()
        .iter()
        .filter(|section| section.kind >= 3)
        .map(|section| (section.offset, section.offset + section.length))
        .collect::<Vec<_>>();

    // Opening an indexed reader reads only metadata from the archive. Every
    // archive GET is a range request and none intersects a payload section.
    let opening_requests = archive_server.range_requests();
    assert!(!opening_requests.is_empty());
    assert!(opening_requests.iter().all(|request| {
        let Some((start, end)) = request.range else {
            return false;
        };
        let start = start as u64;
        let end = end as u64 + 1;
        end <= sequence_section.offset
            || seed_sections
                .iter()
                .all(|(section_start, section_end)| end <= *section_start || start >= *section_end)
    }));
    assert!(
        opening_requests
            .iter()
            .all(|request| request.path == "/candidate.jma")
    );
    assert!(
        archive_server
            .requests()
            .iter()
            .all(|request| !request.path.contains("secret"))
    );
    assert!(
        archive_server
            .requests()
            .iter()
            .all(|request| request.range.is_some() || request.method == "HEAD")
    );
    assert!(
        index_server
            .requests()
            .iter()
            .all(|request| request.path == "/candidate.jma.idx.json")
    );

    let query_level = jam_rs::trace::seeds::extract_seed_level(
        &fixture.sequence,
        SeedSensitivity {
            k: 31,
            scale: 1,
            max_occurrences: 128,
        },
    )
    .unwrap();
    let query_seed = *query_level.seeds.first().expect("retained query seed");
    let bucket = bucket_for(&archive_index, 31, 1, query_seed.hash)
        .expect("query hash bucket")
        .clone();
    let request_count_before_seed = archive_server.range_requests().len();
    let occurrences = reader
        .seed_occurrences_at_scale(query_seed.query(31), 1)
        .unwrap();
    assert!(!occurrences.is_empty());
    let seed_requests = archive_server.range_requests()[request_count_before_seed..].to_vec();
    assert!(!seed_requests.is_empty());
    assert!(request_ranges_cover_exact(
        &seed_requests,
        bucket.offset,
        bucket.length
    ));
    assert!(
        seed_requests
            .iter()
            .all(|request| request.path == "/candidate.jma")
    );

    let sequence_index = archive_index
        .sequence_blocks
        .iter()
        .find(|block| block.start < 47 && block.start + block.base_length > 35)
        .expect("intersecting sequence block")
        .clone();
    let request_count_before_sequence = archive_server.range_requests().len();
    let sequence = reader
        .read_sequence(0, SequenceRange::new(35, 47).unwrap())
        .unwrap();
    assert_eq!(sequence, fixture.sequence[35..47]);
    let sequence_requests =
        archive_server.range_requests()[request_count_before_sequence..].to_vec();
    assert!(!sequence_requests.is_empty());
    assert!(request_ranges_cover_exact(
        &sequence_requests,
        sequence_index.offset,
        sequence_index.length
    ));
    assert!(
        sequence_requests
            .iter()
            .all(|request| request.path == "/candidate.jma")
    );

    let all_archive_gets = archive_server.range_requests();
    assert!(
        all_archive_gets
            .iter()
            .all(|request| { request.range != Some((0, fixture.archive.len().saturating_sub(1))) })
    );
    assert!(reader.metrics().seed_buckets_read > 0);
    assert!(reader.metrics().sequence_blocks_read > 0);
}

#[test]
fn stale_sidecar_identity_is_rejected_before_payload_reads() {
    let fixture = archive_fixture();
    let archive_server = HttpFixture::new(fixture.archive.clone());
    let index_server = HttpFixture::new(fixture.index.clone());
    let reader = JmaReader::open_indexed(
        ObjectResource::open(archive_server.url("/candidate.jma"), archive_options()).unwrap(),
        ObjectResource::open(index_server.url("/candidate.jma.idx.json"), index_options()).unwrap(),
    );
    assert!(reader.is_ok(), "fresh archive and sidecar must bind");

    let mut changed = fixture.archive;
    changed.push(0);
    archive_server.set_body(changed, "fixture-v2");
    let stale = JmaReader::open_indexed(
        ObjectResource::open(archive_server.url("/candidate.jma"), archive_options()).unwrap(),
        ObjectResource::open(index_server.url("/candidate.jma.idx.json"), index_options()).unwrap(),
    )
    .err()
    .expect("stale sidecar must fail closed");
    assert!(stale.to_string().contains("sidecar"));
}

#[test]
fn local_and_remote_indexed_jma_readers_return_equal_sequence_and_seed_evidence() {
    let fixture = archive_fixture();
    let directory = tempfile::Builder::new()
        .prefix("trace-jma-local-remote-")
        .tempdir_in("target")
        .unwrap();
    let archive_path = directory.path().join("candidate.jma");
    let index_path = directory.path().join("candidate.jma.idx.json");
    fs::write(&archive_path, &fixture.archive).unwrap();
    fs::write(&index_path, &fixture.index).unwrap();

    let local = JmaReader::open_indexed(
        open_resource(archive_path.to_string_lossy(), archive_options()).unwrap(),
        open_resource(index_path.to_string_lossy(), index_options()).unwrap(),
    )
    .unwrap();
    let remote_archive = HttpFixture::new(fixture.archive.clone());
    let remote_index = HttpFixture::new(fixture.index.clone());
    let remote = JmaReader::open_indexed(
        open_resource(remote_archive.url("/candidate.jma"), archive_options()).unwrap(),
        open_resource(remote_index.url("/candidate.jma.idx.json"), index_options()).unwrap(),
    )
    .unwrap();

    let query = jam_rs::trace::seeds::extract_seed_level(
        &fixture.sequence,
        SeedSensitivity {
            k: 31,
            scale: 1,
            max_occurrences: 128,
        },
    )
    .unwrap()
    .seeds[0];
    assert_eq!(local.header(), remote.header());
    assert_eq!(local.contigs(), remote.contigs());
    assert_eq!(
        local.seed_occurrences_at_scale(query.query(31), 1).unwrap(),
        remote
            .seed_occurrences_at_scale(query.query(31), 1)
            .unwrap()
    );
    assert_eq!(
        local
            .read_sequence(0, SequenceRange::new(100, 180).unwrap())
            .unwrap(),
        remote
            .read_sequence(0, SequenceRange::new(100, 180).unwrap())
            .unwrap()
    );
}
