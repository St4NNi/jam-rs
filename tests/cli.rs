//! End-to-end tests of the `jam` binary on small generated entries.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn jam(arguments: &[&std::ffi::OsStr]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_jam"))
        .args(arguments)
        .output()
        .unwrap()
}

fn random_bases(seed: u64, length: usize) -> Vec<u8> {
    let mut rng = StdRng::seed_from_u64(seed);
    (0..length)
        .map(|_| b"ACGT"[rng.random_range(0..4)])
        .collect()
}

fn write_fasta(path: &Path, name: &str, sequence: &[u8]) {
    let mut text = format!(">{name} test\n").into_bytes();
    for line in sequence.chunks(70) {
        text.extend_from_slice(line);
        text.push(b'\n');
    }
    std::fs::write(path, text).unwrap();
}

/// Three entries of 20 kb, and a query that copies a stretch of the second one with one
/// substitution.
struct Fixture {
    directory: tempfile::TempDir,
    shard: PathBuf,
    query: PathBuf,
}

fn fixture() -> Fixture {
    let directory = tempfile::tempdir().unwrap();
    let mut inputs = Vec::new();
    let mut second = Vec::new();
    for (number, seed) in [11u64, 22, 33].into_iter().enumerate() {
        let sequence = random_bases(seed, 20_000);
        let path = directory.path().join(format!("plasmid_{number}.fa"));
        write_fasta(&path, &format!("contig_{number}"), &sequence);
        if number == 1 {
            second = sequence;
        }
        inputs.push(path);
    }
    let mut copy = second[5_000..8_000].to_vec();
    copy[1_500] = if copy[1_500] == b'A' { b'C' } else { b'A' };
    let query = directory.path().join("query.fa");
    write_fasta(&query, "planted", &copy);
    let shard = directory.path().join("shard");
    let mut arguments: Vec<&std::ffi::OsStr> = vec![
        "index".as_ref(),
        "build".as_ref(),
        "-o".as_ref(),
        shard.as_os_str(),
        "-t".as_ref(),
        "2".as_ref(),
    ];
    arguments.extend(inputs.iter().map(|path| path.as_os_str()));
    let output = jam(&arguments);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    Fixture {
        directory,
        shard,
        query,
    }
}

fn search(fixture: &Fixture, extra: &[&str]) -> Vec<String> {
    let mut arguments: Vec<&std::ffi::OsStr> = vec![
        "search".as_ref(),
        fixture.shard.as_os_str(),
        fixture.query.as_os_str(),
        "-t".as_ref(),
        "2".as_ref(),
    ];
    arguments.extend(extra.iter().map(std::ffi::OsStr::new));
    let output = jam(&arguments);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect()
}

fn json(line: &str) -> serde_json::Value {
    serde_json::from_str(line).unwrap()
}

#[test]
fn planted_query_confirmed() {
    let fixture = fixture();
    let lines = search(&fixture, &[]);
    let confirmed: Vec<_> = lines
        .iter()
        .map(|line| json(line))
        .filter(|record| record["kind"] == "confirmed")
        .collect();
    assert_eq!(confirmed.len(), 1);
    let record = &confirmed[0];
    assert_eq!(record["entry"], "plasmid_1");
    assert_eq!(record["query"], "planted");
    let decision = &record["decision"];
    assert!(decision["query_coverage"].as_f64().unwrap() > 0.99);
    assert!(decision["identity"].as_f64().unwrap() > 0.99);
    let fragment = &decision["fragments"][0];
    assert_eq!(fragment["contig"], "contig_1");
    assert_eq!(fragment["contig_interval"]["start"], 5_000);
    assert_eq!(fragment["contig_interval"]["end"], 8_000);
    let summary = json(lines.last().unwrap());
    assert_eq!(summary["kind"], "summary");
    assert_eq!(summary["ending"], "complete");
    assert_eq!(summary["query_length"], 3_000);
}

#[test]
fn paf_fragment_lines() {
    let fixture = fixture();
    let lines = search(&fixture, &["--format", "paf"]);
    assert_eq!(lines.len(), 1);
    let columns: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(
        &columns[..6],
        ["planted", "3000", "0", "3000", "+", "contig_1"]
    );
    assert_eq!(&columns[6..9], ["20000", "5000", "8000"]);
    assert_eq!(columns[9], "2999");
    assert_eq!(columns[10], "3000");
    assert_eq!(columns[11], "255");
    assert!(columns.contains(&"tp:A:P"));
    assert!(columns.contains(&"NM:i:1"));
    assert!(columns.contains(&"cg:Z:1500=1X1499="));
}

#[test]
fn sketch_level_possible() {
    let fixture = fixture();
    let lines = search(&fixture, &["--level", "sketch"]);
    let records: Vec<_> = lines.iter().map(|line| json(line)).collect();
    let (summary, possible) = records.split_last().unwrap();
    assert_eq!(summary["kind"], "summary");
    assert!(possible.iter().all(|record| record["kind"] == "possible"));
    // An unrelated entry is possible only by a false seed hit, so its estimate stays low.
    let estimate = |record: &serde_json::Value| record["containment_estimate"].as_f64().unwrap();
    for record in possible {
        let planted = record["entry"] == "plasmid_1";
        assert!(planted == (estimate(record) > 0.9), "{lines:?}");
    }
    assert!(possible.iter().any(|record| record["entry"] == "plasmid_1"));
    let without = summary["sketch"]["unseeded_entries"].as_u64().unwrap();
    assert_eq!(without + possible.len() as u64, 3);
}

#[test]
fn stats_and_errors() {
    let fixture = fixture();
    let output = jam(&[
        "stats".as_ref(),
        fixture.shard.as_os_str(),
        "--json".as_ref(),
    ]);
    let stats: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(stats["entries"], 3);
    assert_eq!(stats["bases"], 60_000);
    let missing = fixture.directory.path().join("missing");
    let output = jam(&["stats".as_ref(), missing.as_os_str()]);
    assert!(!output.status.success());
    let message = String::from_utf8(output.stderr).unwrap();
    assert_eq!(message.lines().count(), 1);
    assert!(message.starts_with("jam: "));
    let output = jam(&["screen".as_ref()]);
    assert!(!output.status.success());
}

#[test]
fn records_name_shards() {
    let fixture = fixture();
    let other_input = fixture.directory.path().join("other.fa");
    write_fasta(&other_input, "other_contig", &random_bases(44, 20_000));
    let other = fixture.directory.path().join("other_shard");
    let output = jam(&[
        "index".as_ref(),
        "build".as_ref(),
        "-o".as_ref(),
        other.as_os_str(),
        other_input.as_os_str(),
    ]);
    assert!(output.status.success());
    let output = jam(&[
        "search".as_ref(),
        fixture.shard.as_os_str(),
        other.as_os_str(),
        fixture.query.as_os_str(),
        "-t".as_ref(),
        "2".as_ref(),
    ]);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).unwrap();
    let records: Vec<_> = text.lines().map(json).collect();
    assert!(records.iter().all(|record| record["shard"].is_string()));
    let confirmed: Vec<_> = records
        .iter()
        .filter(|record| record["kind"] == "confirmed")
        .collect();
    assert_eq!(confirmed.len(), 1);
    assert_eq!(confirmed[0]["shard"], "shard");
    let summaries: Vec<_> = records
        .iter()
        .filter(|record| record["kind"] == "summary")
        .map(|record| {
            (
                record["shard"].as_str().unwrap(),
                record["ending"].as_str().unwrap(),
            )
        })
        .collect();
    assert_eq!(
        summaries,
        [("shard", "complete"), ("other_shard", "complete")]
    );
}

#[test]
fn incomplete_search_fails() {
    let fixture = fixture();
    let output = jam(&[
        "search".as_ref(),
        fixture.shard.as_os_str(),
        fixture.query.as_os_str(),
        "--deadline".as_ref(),
        "0.000000001".as_ref(),
    ]);
    assert!(!output.status.success());
    let text = String::from_utf8(output.stdout).unwrap();
    assert_eq!(json(text.lines().last().unwrap())["ending"], "deadline");
    let message = String::from_utf8(output.stderr).unwrap();
    assert_eq!(message.trim(), "jam: 1 searches did not complete");
}
