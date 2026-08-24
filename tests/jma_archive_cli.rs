use jam_rs::jma::reader::JmaReader;
use jam_rs::jma::{ArchiveReader, SequenceRange};
use jam_rs::resource::ResourceOpenOptions;
use jam_rs::resource::local::LocalResource;
use jam_rs::sequence::BlockCodec;
use std::fs;
use std::process::Command;

#[test]
fn archive_cli_builds_one_object_for_fixed_and_gear_codecs() {
    let directory = tempfile::Builder::new()
        .prefix("jma-archive-cli-")
        .tempdir_in("target")
        .unwrap();
    let input = directory.path().join("assembly.fa");
    let sequence = b"ACGTRYKMBVDHN-ACGT".repeat(1000);
    fs::write(
        &input,
        format!(">contig\n{}\n", String::from_utf8_lossy(&sequence)),
    )
    .unwrap();

    for (name, policy, codec) in [
        ("fixed-raw.jma", "fixed", "raw2bit"),
        ("fixed-zstd.jma", "fixed", "zstd2bit"),
        ("gear-zstd.jma", "gear", "zstd2bit"),
    ] {
        let output = directory.path().join(name);
        let result = Command::new(env!("CARGO_BIN_EXE_jam"))
            .args(["--silent", "archive", "--input"])
            .arg(&input)
            .args(["--output"])
            .arg(&output)
            .args([
                "--sequence-block-policy",
                policy,
                "--sequence-block-codec",
                codec,
                "--block-bases",
                "4096",
                "--gear-min-bases",
                "1024",
                "--gear-target-bases",
                "2048",
                "--gear-max-bases",
                "4096",
                "--gear-table",
                "dinucleotide",
                "--primary-scale",
                "1",
                "--rescue-scale",
                "1",
            ])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "archive failed: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert!(output.is_file());
        assert!(!output.with_extension("jma.idx.json").exists());

        let reader = JmaReader::from_resource(
            LocalResource::from_path(&output, ResourceOpenOptions::default()).unwrap(),
        )
        .unwrap();
        assert_eq!(
            reader
                .read_sequence(0, SequenceRange::new(0, sequence.len() as u64).unwrap())
                .unwrap(),
            sequence
        );
        let expected_codec = if codec == "raw2bit" {
            BlockCodec::Raw2Bit
        } else {
            BlockCodec::Zstd2Bit
        };
        assert!(
            reader
                .sequence_index()
                .iter()
                .all(|record| record.codec == expected_codec)
        );
    }
}
