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
    let shard = directory.path().join("shard.jam");
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

/// The standard output and error of a successful search.
fn search_output(fixture: &Fixture, extra: &[&str]) -> (String, String) {
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
    let text = |bytes| String::from_utf8(bytes).unwrap();
    (text(output.stdout), text(output.stderr))
}

fn search(fixture: &Fixture, extra: &[&str]) -> Vec<String> {
    let (stdout, _) = search_output(fixture, extra);
    stdout.lines().map(str::to_string).collect()
}

fn json(line: &str) -> serde_json::Value {
    serde_json::from_str(line).unwrap()
}

#[test]
fn planted_query_confirmed() {
    let fixture = fixture();
    let lines = search(&fixture, &[]);
    assert_eq!(lines.len(), 1, "{lines:?}");
    let record = json(&lines[0]);
    assert_eq!(
        (record["q"].as_str(), record["a"].as_str()),
        (Some("planted"), Some("plasmid_1"))
    );
    assert_eq!(record["s"], 3);
    assert!(record["cov"].as_f64().unwrap() > 0.99);
    assert!(record["id"].as_f64().unwrap() > 0.99);
    assert_eq!(record["bases"], 3_000);
    let fragment = serde_json::json!(["contig_1", 0, 3000, 5000, 8000, "+", "1500=1X1499="]);
    assert_eq!(record["f"], serde_json::json!([fragment]));
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

/// A --max-miss below the seed miss bound turns the pruning off, which -v marks; the planted
/// copy is still found.
#[test]
fn max_miss_unpruned() {
    let fixture = fixture();
    let (stdout, stderr) = search_output(&fixture, &["--max-miss", "0.000000000000001", "-v"]);
    assert_eq!(json(stdout.lines().next().unwrap())["a"], "plasmid_1");
    assert!(stderr.contains("; pruning off\n"), "{stderr}");
    assert!(
        stderr.contains("jam: 1 of 1 queries are above the miss target --max-miss"),
        "{stderr}"
    );
    assert!(!stdout.contains("jam:"), "{stdout}");
    let (_, stderr) = search_output(&fixture, &["-v"]);
    assert!(
        !stderr.contains("pruning off") && !stderr.contains("jam:"),
        "{stderr}"
    );
    assert!(stderr.contains("; miss "), "{stderr}");
}

#[test]
fn tsv_lines() {
    let fixture = fixture();
    let lines = search(&fixture, &["--format", "tsv"]);
    assert_eq!(
        lines,
        [
            "query\tassembly\tcoverage\tidentity\tbases\tfragments",
            "planted\tplasmid_1\t1\t0.9997\t3000\t1"
        ]
    );
    let output = jam(&[
        "search".as_ref(),
        fixture.shard.as_os_str(),
        fixture.query.as_os_str(),
        "--format".as_ref(),
        "tsv".as_ref(),
        "--level".as_ref(),
        "sketch".as_ref(),
    ]);
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
}

#[test]
fn sketch_level_possible() {
    let fixture = fixture();
    let lines = search(&fixture, &["--level", "sketch"]);
    let records: Vec<_> = lines.iter().map(|line| json(line)).collect();
    let possible = &records;
    assert!(possible.iter().all(|record| record["s"] == 1));
    // An unrelated entry is possible only by a false seed hit, so its estimate stays low.
    let estimate = |record: &serde_json::Value| record["c"].as_f64().unwrap();
    for record in possible {
        let planted = record["a"] == "plasmid_1";
        assert!(planted == (estimate(record) > 0.9), "{lines:?}");
    }
    assert!(possible.iter().any(|record| record["a"] == "plasmid_1"));
}

#[test]
fn locate_level_windows() {
    let fixture = fixture();
    let lines = search(&fixture, &["--level", "locate"]);
    let records: Vec<_> = lines.iter().map(|line| json(line)).collect();
    let planted = records.iter().find(|record| record["a"] == "plasmid_1");
    let planted = planted.unwrap();
    assert_eq!(planted["s"], 2);
    assert!(planted["c"].as_f64().unwrap() > 0.9);
    let windows = planted["w"].as_array().unwrap();
    let covers = |window: &serde_json::Value| {
        window[0] == "contig_1"
            && window[1].as_u64().unwrap() <= 5_000
            && window[2].as_u64().unwrap() >= 8_000
    };
    assert!(windows.iter().any(covers), "{lines:?}");
}

/// A query without a hit writes nothing to the standard output.
#[test]
fn no_hit_silent() {
    let fixture = fixture();
    let other = fixture.directory.path().join("other.fa");
    write_fasta(&other, "stranger", &random_bases(77, 3_000));
    let output = jam(&[
        "search".as_ref(),
        fixture.shard.as_os_str(),
        other.as_os_str(),
        "-t".as_ref(),
        "2".as_ref(),
    ]);
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert!(output.stderr.is_empty());
}

/// The standard output holds only matches; -v adds a progress line and -vv the debug lines,
/// both on the standard error.
#[test]
fn verbose_on_stderr() {
    let fixture = fixture();
    let (quiet, nothing) = search_output(&fixture, &[]);
    assert!(nothing.is_empty(), "{nothing}");
    let (stdout, stderr) = search_output(&fixture, &["-v"]);
    assert_eq!(stdout.lines().count(), quiet.lines().count());
    let progress: Vec<&str> = stderr.lines().collect();
    assert_eq!(progress.len(), 1, "{stderr}");
    assert!(progress[0].starts_with("planted: complete; "), "{stderr}");
    assert!(progress[0].contains(" 1 confirmed; miss "), "{stderr}");
    assert!(progress[0].contains("; sketch "), "{stderr}");
    let (stdout, stderr) = search_output(&fixture, &["-vv"]);
    assert!(stdout.lines().all(|line| json(line).get("kind").is_none()));
    let (debug, progress): (Vec<&str>, Vec<&str>) =
        stderr.lines().partition(|line| line.starts_with('{'));
    assert_eq!(progress.len(), 1, "{stderr}");
    let kinds: Vec<String> = debug
        .iter()
        .map(|line| json(line)["kind"].to_string())
        .collect();
    assert!(kinds.contains(&"\"possible\"".to_string()), "{stderr}");
    assert!(kinds.contains(&"\"located\"".to_string()), "{stderr}");
    assert!(!kinds.contains(&"\"confirmed\"".to_string()), "{stderr}");
    assert_eq!(kinds.last().map(String::as_str), Some("\"summary\""));
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
fn shards_share_summary() {
    let fixture = fixture();
    let other_input = fixture.directory.path().join("other.fa");
    write_fasta(&other_input, "other_contig", &random_bases(44, 20_000));
    let other = fixture.directory.path().join("other_shard.jam");
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
    assert_eq!(records.len(), 1, "{text}");
    assert_eq!(records[0]["a"], "plasmid_1");
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
    assert!(output.stdout.is_empty());
    // The problem line comes without -v, then the error.
    let message = String::from_utf8(output.stderr).unwrap();
    let mut lines = message.lines();
    assert_eq!(
        lines.next(),
        Some("planted: stopped (deadline)"),
        "{message}"
    );
    assert!(
        lines
            .next()
            .unwrap()
            .starts_with("jam: 1 queries did not complete"),
        "{message}"
    );
}

/// Writes a FASTA file with one record per named sequence.
fn write_records(path: &Path, records: &[(&str, Vec<u8>)]) {
    let mut text = Vec::new();
    for (name, sequence) in records {
        text.extend_from_slice(format!(">{name} contig\n").as_bytes());
        text.extend_from_slice(sequence);
        text.push(b'\n');
    }
    std::fs::write(path, text).unwrap();
}

/// Two assemblies, one holding a copy of the second entry with one substitution, and a list
/// path that does not exist: the copy is found, and the missing assembly fails on its own.
#[test]
fn screen_finds_copy() {
    let fixture = fixture();
    let directory = fixture.directory.path();
    let mut copy = random_bases(22, 20_000);
    copy[10_000] = if copy[10_000] == b'A' { b'C' } else { b'A' };
    let first = directory.join("first.fa");
    write_records(
        &first,
        &[("chromosome", random_bases(55, 30_000)), ("planted", copy)],
    );
    let second = directory.join("second.fa");
    write_records(&second, &[("other", random_bases(66, 30_000))]);
    let list = directory.join("assemblies.txt");
    let missing = directory.join("missing.fa");
    let paths = [&first, &second, &missing].map(|path| path.display().to_string());
    std::fs::write(&list, paths.join("\n")).unwrap();
    let output = jam(&[
        "screen".as_ref(),
        fixture.shard.as_os_str(),
        "--assemblies".as_ref(),
        list.as_os_str(),
        "-t".as_ref(),
        "2".as_ref(),
    ]);
    assert!(!output.status.success(), "the missing assembly must fail");
    let stdout = String::from_utf8(output.stdout).unwrap();
    let records: Vec<_> = stdout.lines().map(json).collect();
    let of = |kind: &str| -> Vec<&serde_json::Value> {
        records
            .iter()
            .filter(|record| record["type"] == kind)
            .collect()
    };
    let entries = of("entry");
    assert_eq!(entries.len(), 1, "{stdout}");
    let entry = entries[0];
    assert_eq!(
        (entry["entry"].as_str(), entry["assembly"].as_str()),
        (Some("plasmid_1"), Some("first.fa"))
    );
    assert_eq!(entry["contigs"][0]["name"], "planted");
    let alignment = &entry["primary"][0]["fragment"]["alignment"];
    assert_eq!(alignment["cigar"], "10000=1X9999=");
    assert_eq!(of("assembly").len(), 2, "{stdout}");
    let failed = of("failed");
    assert_eq!(failed.len(), 1, "{stdout}");
    assert_eq!(failed[0]["assembly"], "missing.fa");
}

/// Every option that `--help` of a command shows is described in the README.
#[test]
fn readme_lists_options() {
    let readme = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"));
    let readme = readme.unwrap();
    let commands: [&[&str]; 8] = [
        &["index", "build"],
        &["index", "reshard"],
        &["index", "warm"],
        &["index", "copy"],
        &["stats"],
        &["search"],
        &["serve"],
        &["screen"],
    ];
    let mut missing = Vec::new();
    for command in commands {
        assert!(readme.contains(&format!("### jam {}", command.join(" "))));
        let mut arguments: Vec<&std::ffi::OsStr> = command.iter().map(|w| w.as_ref()).collect();
        arguments.push("--help".as_ref());
        let help = String::from_utf8(jam(&arguments).stdout).unwrap();
        let options = help
            .lines()
            .map(str::trim_start)
            .filter(|line| line.starts_with('-'))
            .flat_map(|line| line.split([' ', ',']))
            .filter(|word| word.starts_with("--") && *word != "--help");
        for option in options.map(|word| word.trim_end_matches('.')) {
            if !mentions(&readme, option) {
                missing.push(format!("jam {} {option}", command.join(" ")));
            }
        }
    }
    assert!(missing.is_empty(), "not in the README: {missing:?}");
}

/// Whether `text` names `option` as a whole word, so `--k` does not match `--key`.
fn mentions(text: &str, option: &str) -> bool {
    text.match_indices(option).any(|(at, _)| {
        let next = text[at + option.len()..].chars().next();
        !next.is_some_and(|c| c.is_ascii_alphanumeric() || c == '-')
    })
}
