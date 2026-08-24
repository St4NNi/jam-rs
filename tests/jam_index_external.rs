use flate2::Compression;
use flate2::write::GzEncoder;
use jam_rs::jam_index::{
    JamIndexBuildConfig, JamIndexTraceConfig, MetagenomeSource, build_jam_index, trace_index,
};
use jam_rs::trace::config::{SensitivityConfig, SensitivityProfile};
use std::fs;
use std::io::Write;
use tempfile::Builder;

fn gzip_block(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

fn write_fai(path: &std::path::Path) {
    fs::write(
        format!("{}.fai", path.display()),
        "background\t5000\t12\t5000\t5001\ntarget\t160\t5021\t160\t161\n",
    )
    .unwrap();
}

#[test]
fn external_modes_trace() {
    let sources = Builder::new()
        .prefix("index-external-source-")
        .tempdir_in("target")
        .unwrap();
    let outputs = Builder::new()
        .prefix("index-external-output-")
        .tempdir_in("target")
        .unwrap();
    let query = (0..160)
        .map(|index| b"ACGT"[(index * 13 + index / 7) & 3])
        .collect::<Vec<_>>();
    let fasta = format!(
        ">background\n{}\n>target\n{}\n",
        "A".repeat(5000),
        String::from_utf8_lossy(&query)
    );

    let plain = sources.path().join("plain.fa");
    fs::write(&plain, &fasta).unwrap();
    write_fai(&plain);

    let gzip = sources.path().join("stream.fa.gz");
    fs::write(&gzip, gzip_block(fasta.as_bytes())).unwrap();

    let bgzf = sources.path().join("indexed.fa.gz");
    let split = 5013usize;
    let first = gzip_block(&fasta.as_bytes()[..split]);
    let second = gzip_block(&fasta.as_bytes()[split..]);
    let mut compressed = first.clone();
    compressed.extend_from_slice(&second);
    fs::write(&bgzf, compressed).unwrap();
    write_fai(&bgzf);
    let mut gzi = 1u64.to_le_bytes().to_vec();
    gzi.extend_from_slice(&u64::try_from(first.len()).unwrap().to_le_bytes());
    gzi.extend_from_slice(&u64::try_from(split).unwrap().to_le_bytes());
    fs::write(format!("{}.gzi", bgzf.display()), gzi).unwrap();

    let unrelated = sources.path().join("unrelated.fa.gz");
    fs::write(
        &unrelated,
        gzip_block(format!(">noise\n{}\n", "C".repeat(5000)).as_bytes()),
    )
    .unwrap();
    let index = outputs.path().join("index");
    build_jam_index(
        &index,
        &[
            MetagenomeSource {
                metagenome_id: "plain".to_string(),
                sequence_path: plain.clone(),
            },
            MetagenomeSource {
                metagenome_id: "bgzf".to_string(),
                sequence_path: bgzf.clone(),
            },
            MetagenomeSource {
                metagenome_id: "gzip".to_string(),
                sequence_path: gzip.clone(),
            },
            MetagenomeSource {
                metagenome_id: "unrelated".to_string(),
                sequence_path: unrelated.clone(),
            },
        ],
        &JamIndexBuildConfig {
            target_parts: 4,
            parallel_parts: 4,
            source_manifest_sha256: "a".repeat(64),
            ..JamIndexBuildConfig::default()
        },
    )
    .unwrap();
    fs::remove_file(unrelated).unwrap();

    let mut config = JamIndexTraceConfig::default();
    config.runner.sensitivity = SensitivityConfig::for_profile(SensitivityProfile::Sensitive);
    config.runner.threads = 1;
    config.runner.io_concurrency = 1;
    config.parallel_candidates = 1;
    let output = trace_index(&index, "query", &query, &config).unwrap();
    assert_eq!(output.metagenomes.len(), 3);
    for result in &output.metagenomes {
        assert_eq!(
            result.coverage.as_ref().unwrap().supported_bases,
            query.len() as u64
        );
    }
    let bytes = |id: &str| {
        output
            .metagenomes
            .iter()
            .find(|result| result.metagenome_id == id)
            .unwrap()
            .archive_metrics
            .as_ref()
            .unwrap()
            .sequence_bytes_read
    };
    assert_eq!(bytes("plain"), 160);
    assert!(bytes("bgzf") < fs::metadata(bgzf).unwrap().len());
    assert_eq!(bytes("gzip"), fs::metadata(gzip).unwrap().len());
}
