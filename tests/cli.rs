//! End-to-end tests of the `jam` binary on small generated entries.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn jam(arguments: &[&std::ffi::OsStr]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_jam"))
        .args(arguments)
        .output()
        .unwrap()
}

fn random_bases(seed: u64, length: usize) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 33) as usize % 4]
        })
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
fn a_planted_query_is_confirmed_with_its_alignment() {
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
fn paf_output_has_one_standard_line_per_fragment() {
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
fn the_sketch_level_names_possible_entries_only() {
    let fixture = fixture();
    let lines = search(&fixture, &["--level", "sketch"]);
    let records: Vec<_> = lines.iter().map(|line| json(line)).collect();
    let (summary, possible) = records.split_last().unwrap();
    assert_eq!(summary["kind"], "summary");
    assert!(possible.iter().all(|record| record["kind"] == "possible"));
    // A shard this small treats every seed as frequent, so unrelated entries stay possible.
    let planted = possible
        .iter()
        .find(|record| record["entry"] == "plasmid_1")
        .unwrap();
    assert!(planted["containment_estimate"].as_f64().unwrap() > 0.9);
}

#[test]
fn stats_and_errors_are_plain() {
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
