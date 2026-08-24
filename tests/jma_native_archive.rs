use jam_rs::archive::{NativeJmaArchive, SeedKey, SequenceRequest, TraceArchive};
use jam_rs::jma::builder::{ArchiveBuildConfig, write_archive_from_fasta};
use jam_rs::resource::ResourceOpenOptions;
use jam_rs::resource::local::LocalResource;
use jam_rs::trace::config::SeedSensitivity;
use jam_rs::trace::seeds::extract_seed_level;
use std::collections::BTreeMap;
use std::fs;

fn dna(mut state: u64, length: usize) -> Vec<u8> {
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            b"ACGT"[((state >> 62) & 3) as usize]
        })
        .collect()
}

#[test]
fn native_archive_uses_mmap_and_batches_sequence_ranges() {
    let directory = tempfile::Builder::new()
        .prefix("native-jma-archive-")
        .tempdir_in("target")
        .unwrap();
    let sequence = dna(77, 20_000);
    let fasta = directory.path().join("assembly.fa");
    fs::write(
        &fasta,
        format!(">contig\n{}\n", String::from_utf8_lossy(&sequence)),
    )
    .unwrap();
    let archive = directory.path().join("assembly.jma");
    write_archive_from_fasta(
        &fasta,
        &archive,
        ArchiveBuildConfig {
            sequence_policy: jam_rs::sequence::SequenceBlockPolicy::Fixed { block_bases: 4096 },
            k31_scale: 1,
            k21_scale: Some(1),
            ..ArchiveBuildConfig::default()
        },
    )
    .unwrap();
    let archive_size = fs::metadata(&archive).unwrap().len();
    let native = NativeJmaArchive::from_resource(
        LocalResource::from_path(&archive, ResourceOpenOptions::default()).unwrap(),
    )
    .unwrap();

    let metadata = native.metadata().unwrap();
    assert_eq!(metadata.contigs.len(), 1);
    assert_eq!(metadata.total_bases, sequence.len() as u64);
    let opening_metrics = native.metrics();
    assert_eq!(opening_metrics.mapped_bytes, archive_size);
    assert!(opening_metrics.resident_bytes > 0);
    assert!(opening_metrics.resource.returned_bytes < archive_size);

    let scheme = native
        .available_seed_schemes()
        .unwrap()
        .into_iter()
        .find(|scheme| scheme.span == 31 && scheme.density_parameter == 1)
        .unwrap();
    let seeds = extract_seed_level(
        &sequence,
        SeedSensitivity {
            k: 31,
            scale: 1,
            max_occurrences: 256,
        },
    )
    .unwrap()
    .seeds;
    let mut by_prefix = BTreeMap::<u32, Vec<_>>::new();
    for seed in seeds {
        let prefix = jam_rs::jma::format::hash_prefix(seed.hash, scheme.bucket_bits).unwrap();
        by_prefix.entry(prefix).or_default().push(seed);
    }
    let selected = by_prefix.values().find(|seeds| seeds.len() >= 2).unwrap();
    let keys = selected[..2]
        .iter()
        .map(|seed| SeedKey {
            digest: seed.hash,
            verification: seed.canonical_kmer.to_be_bytes().to_vec(),
        })
        .collect::<Vec<_>>();
    let lookup = native
        .lookup_seeds(jam_rs::archive::SeedSchemeId(scheme.scheme_id), &keys)
        .unwrap();
    assert_eq!(lookup.matches.len(), 2);
    assert!(
        lookup
            .matches
            .iter()
            .all(|seed_match| !seed_match.occurrences.is_empty())
    );
    assert_eq!(lookup.metrics.pages_read, 1);

    let requests = [
        SequenceRequest::new(0, 50, 150, false).unwrap(),
        SequenceRequest::new(0, 10, 100, true).unwrap(),
        SequenceRequest::new(0, 10, 100, false).unwrap(),
    ];
    let slices = native.read_sequences(&requests).unwrap();
    assert_eq!(slices[0].bases, sequence[50..150]);
    assert_eq!(slices[2].bases, sequence[10..100]);
    let mut reverse = sequence[10..100].to_vec();
    reverse.reverse();
    for base in &mut reverse {
        *base = jam_rs::sequence::complement_base(*base).unwrap();
    }
    assert_eq!(slices[1].bases, reverse);
    assert_eq!(native.metrics().resource.full_object_fallbacks, 0);
}
