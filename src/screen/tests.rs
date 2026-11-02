//! End-to-end tests of the screen on small generated entries and assemblies.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::align::{Strand, format_cigar};
use crate::decide::{Fragment, QueryCoverage};
use crate::index::{EntryInput, Shard, build_shard};
use crate::seeds::reverse_complement;
use crate::test_support::{build_options, no_stop, random_bases, rng, with_substitutions};
use std::path::Path;

#[test]
fn low_complexity_walks_off_by_default() {
    assert!(!ScreenConfig::default().low_complexity_walks);
}

/// An index with one entry per item, each a list of contigs.
pub(super) fn index_of(directory: &Path, entries: &[(&str, &[&[u8]])]) -> Shard {
    let inputs: Vec<EntryInput> = entries
        .iter()
        .map(|(name, contigs)| {
            let path = directory.join(format!("{name}.fa"));
            let mut text = Vec::new();
            for (index, contig) in contigs.iter().enumerate() {
                text.extend(format!(">c{index}\n").bytes());
                text.extend_from_slice(contig);
                text.push(b'\n');
            }
            std::fs::write(&path, text).unwrap();
            EntryInput {
                name: (*name).into(),
                path,
            }
        })
        .collect();
    let options = build_options(directory, 2);
    let output = directory.join("index.jam");
    build_shard(&inputs, &output, &options).unwrap();
    Shard::open(&output).unwrap()
}

pub(super) fn plasmids(directory: &Path, plasmids: &[(&str, &[u8])]) -> Shard {
    let entries: Vec<(&str, Vec<&[u8]>)> = plasmids
        .iter()
        .map(|&(name, sequence)| (name, vec![sequence]))
        .collect();
    let entries: Vec<(&str, &[&[u8]])> = entries
        .iter()
        .map(|(name, contigs)| (*name, contigs.as_slice()))
        .collect();
    index_of(directory, &entries)
}

pub(super) fn write_assembly(directory: &Path, contigs: &[(&str, Vec<u8>)]) -> std::path::PathBuf {
    let mut fasta = Vec::new();
    for (name, sequence) in contigs {
        fasta.extend(format!(">{name} extra words\n").bytes());
        fasta.extend(sequence);
        fasta.push(b'\n');
    }
    let path = directory.join("assembly.fa");
    std::fs::write(&path, fasta).unwrap();
    path
}

pub(super) fn screen_path(
    shard: &Shard,
    path: &Path,
    config: &ScreenConfig,
) -> (Vec<ScreenEntry>, AssemblySummary) {
    let control = no_stop();
    let screener = Screener::new(shard, config.clone()).unwrap();
    let mut reported = Vec::new();
    let summary = screener
        .screen_assembly(path, "assembly", control, &mut |report| {
            reported.push(report);
            Ok(())
        })
        .unwrap();
    (reported, summary)
}

pub(super) fn screen(
    shard: &Shard,
    directory: &Path,
    contigs: &[(&str, Vec<u8>)],
    config: &ScreenConfig,
) -> (Vec<ScreenEntry>, AssemblySummary) {
    let path = write_assembly(directory, contigs);
    screen_path(shard, &path, config)
}

fn close_lineage() -> ScreenConfig {
    ScreenConfig::new(DecisionRule {
        min_identity: 0.99,
        min_fragment_identity: 0.0,
        min_query_coverage: "0.8".parse().unwrap(),
        min_contig_coverage: 0.9,
        min_contig_identity: 0.99,
        ..DecisionRule::default()
    })
}

fn names(reported: &[ScreenEntry]) -> Vec<String> {
    reported.iter().map(|report| report.entry.clone()).collect()
}

/// The primary fragments of a reported entry: entry interval, assembly contig name, contig
/// interval, strand and CIGAR, sorted.
fn placed(report: &ScreenEntry) -> Vec<(u64, u64, String, u64, u64, Strand, String)> {
    let mut placed: Vec<_> = report
        .primary
        .iter()
        .map(|selected| &selected.fragment)
        .map(|fragment: &Fragment| {
            let alignment = &fragment.alignment;
            alignment.validate().unwrap();
            let contig = report
                .contigs
                .iter()
                .find(|contig| contig.id == fragment.contig_index)
                .unwrap();
            (
                alignment.query_interval.start,
                alignment.query_interval.end,
                contig.name.clone(),
                alignment.contig_interval.start,
                alignment.contig_interval.end,
                alignment.strand,
                format_cigar(&alignment.runs),
            )
        })
        .collect();
    placed.sort();
    placed
}

#[test]
fn planted_plasmid_found() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_bases(&mut rng(1), 6_000);
    let other = random_bases(&mut rng(2), 6_000);
    let shard = plasmids(
        directory.path(),
        &[("plasmid", &plasmid), ("other", &other)],
    );
    // The copy in the assembly lacks plasmid bases 3,000 to 3,002 (bases only in the entry, `I`)
    // and holds two extra bases after plasmid base 4,500 (bases only in the contig, `D`), and
    // sits between 2,000 bases of something else on each side.
    let copy = [
        &plasmid[..3_000],
        &plasmid[3_003..4_500],
        b"GT",
        &plasmid[4_500..],
    ]
    .concat();
    let contig = [
        &random_bases(&mut rng(3), 2_000)[..],
        &copy,
        &random_bases(&mut rng(4), 2_000),
    ]
    .concat();
    let contigs = [
        ("chromosome", contig),
        ("unrelated", random_bases(&mut rng(5), 5_000)),
    ];
    let (reported, summary) = screen(&shard, directory.path(), &contigs, &ScreenConfig::default());
    assert_eq!(names(&reported), ["plasmid"]);
    let report = &reported[0];
    assert_eq!((report.entry_length, report.covered_bases), (6_000, 5_997));
    assert_eq!(
        (report.insertions, report.deletions, report.substitutions),
        (3, 2, 0)
    );
    assert_eq!(report.identity, 5_997.0 / 6_002.0);
    assert_eq!(placed(report).len(), 1);
    let (start, end, name, contig_start, contig_end, strand, cigar) = &placed(report)[0];
    assert_eq!((*start, *end, name.as_str()), (0, 6_000, "chromosome"));
    assert_eq!((*contig_start, *contig_end), (2_000, 2_000 + 5_999));
    assert_eq!(*strand, Strand::Forward);
    // The gap that a base of the run of `T` and `G` may fall into shifts the CIGAR a little.
    assert!(
        cigar.starts_with("3000=3I") && cigar.contains("2D"),
        "{cigar}"
    );
    assert_eq!(report.contigs.len(), 1);
    let found = &report.contigs[0];
    assert_eq!((found.id, found.length, found.aligned), (0, 9_999, 5_999));
    assert_eq!(
        (found.score, found.best_score, found.entries),
        (found.score, found.score, 1)
    );
    assert_eq!(
        (
            summary.contigs,
            summary.reported_entries,
            summary.reported_contigs
        ),
        (2, 1, 1)
    );
}

#[test]
fn both_strands_placed() {
    let directory = tempfile::tempdir().unwrap();
    let first = random_bases(&mut rng(5), 3_000);
    let second = random_bases(&mut rng(6), 3_000);
    let shard = plasmids(directory.path(), &[("first", &first), ("second", &second)]);
    let contigs = [
        ("forward", first[200..1_800].to_vec()),
        ("unrelated", random_bases(&mut rng(7), 2_000)),
        ("reverse", reverse_complement(&first[1_900..2_900])),
        // A third of the second entry: below the query coverage asked for.
        ("partial", second[..1_000].to_vec()),
    ];
    for batch_bases in [500_000, 2_500] {
        let config = ScreenConfig {
            batch_bases,
            ..ScreenConfig::default()
        };
        let (reported, _) = screen(&shard, directory.path(), &contigs, &config);
        assert_eq!(names(&reported), ["first"], "batch {batch_bases}");
        let report = &reported[0];
        assert_eq!(report.entry_length, 3_000);
        let placed = placed(report);
        assert_eq!(
            placed,
            [
                (
                    200,
                    1_800,
                    "forward".into(),
                    0,
                    1_600,
                    Strand::Forward,
                    "1600=".into()
                ),
                (
                    1_900,
                    2_900,
                    "reverse".into(),
                    0,
                    1_000,
                    Strand::Reverse,
                    "1000=".into()
                ),
            ],
            "batch {batch_bases}"
        );
        assert_eq!(report.covered_bases, 2_600);
    }
}

#[test]
fn split_plasmid_counts() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_bases(&mut rng(31), 6_000);
    let shard = plasmids(directory.path(), &[("plasmid", &plasmid)]);
    let contigs = [
        ("one", plasmid[..2_100].to_vec()),
        ("two", plasmid[2_200..4_200].to_vec()),
        ("three", plasmid[4_300..].to_vec()),
    ];
    let (reported, _) = screen(&shard, directory.path(), &contigs, &close_lineage());
    assert_eq!(reported.len(), 1);
    assert_eq!(reported[0].contigs.len(), 3);
    assert_eq!(reported[0].covered_bases, 5_800);
    // One of the three and a contig from another place cover less than asked for.
    let few = [
        contigs[0].clone(),
        ("odd", random_bases(&mut rng(32), 2_000)),
    ];
    assert!(
        screen(&shard, directory.path(), &few, &close_lineage())
            .0
            .is_empty()
    );
}

#[test]
fn spacer_blocks_joins() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_bases(&mut rng(15), 3_000);
    let shard = plasmids(directory.path(), &[("plasmid", &plasmid)]);
    // The two contigs continue each other on the plasmid, so an alignment of the batch could run
    // through the spacer.
    let contigs = [
        ("left", plasmid[..1_500].to_vec()),
        ("right", plasmid[1_500..].to_vec()),
    ];
    let (reported, _) = screen(&shard, directory.path(), &contigs, &ScreenConfig::default());
    assert_eq!(reported.len(), 1);
    assert_eq!(
        placed(&reported[0]),
        [
            (
                0,
                1_500,
                "left".into(),
                0,
                1_500,
                Strand::Forward,
                "1500=".into()
            ),
            (
                1_500,
                3_000,
                "right".into(),
                0,
                1_500,
                Strand::Forward,
                "1500=".into()
            ),
        ]
    );
    assert_eq!(reported[0].identity, 1.0);
}

#[test]
fn batch_border_whole() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_bases(&mut rng(81), 4_000);
    let shard = plasmids(directory.path(), &[("plasmid", &plasmid)]);
    // A contig longer than a batch goes in pieces; the plasmid sits over the cut at 5,000.
    let long = [
        &random_bases(&mut rng(82), 3_000)[..],
        &plasmid,
        &random_bases(&mut rng(83), 5_000),
    ]
    .concat();
    // A contig of its own right after it starts the next batch.
    let contigs = [("long", long), ("next", plasmid[..2_000].to_vec())];
    let config = ScreenConfig {
        batch_bases: 5_000,
        ..ScreenConfig::default()
    };
    let (reported, summary) = screen(&shard, directory.path(), &contigs, &config);
    assert!(summary.batches >= 3, "{summary:?}");
    assert_eq!(reported.len(), 1);
    // The pieces overlap by 1,250 bases; each keeps its part up to the middle of the overlap.
    let report = &reported[0];
    let long = report.contigs.iter().find(|c| c.name == "long").unwrap().id;
    let fragments =
        (report.primary.iter().map(|selected| &selected.fragment)).chain(&report.alternatives);
    let mut on_long: Vec<(u64, u64, u64)> = fragments
        .filter(|fragment| fragment.contig_index == long)
        .map(|fragment| {
            let alignment = &fragment.alignment;
            let entry = alignment.query_interval;
            (entry.start, entry.end, alignment.contig_interval.start)
        })
        .collect();
    on_long.sort_unstable();
    assert_eq!(on_long, [(0, 1_375, 3_000), (1_375, 4_000, 4_375)]);
    assert_eq!(report.covered_bases, 4_000);
}

#[test]
fn entry_contigs_separate() {
    let directory = tempfile::tempdir().unwrap();
    let (first, second) = (
        random_bases(&mut rng(91), 3_000),
        random_bases(&mut rng(92), 2_000),
    );
    let shard = index_of(directory.path(), &[("pair", &[&first, &second])]);
    let contigs = [("copy", second.clone())];
    let (reported, _) = screen(&shard, directory.path(), &contigs, &ScreenConfig::default());
    assert_eq!(names(&reported), ["pair/c1"]);
    assert_eq!(
        (reported[0].entry_length, reported[0].covered_bases),
        (2_000, 2_000)
    );
}

#[test]
fn sliver_not_aligned() {
    let directory = tempfile::tempdir().unwrap();
    let first = random_bases(&mut rng(5), 3_000);
    let second = random_bases(&mut rng(6), 3_000);
    let shard = plasmids(directory.path(), &[("first", &first), ("second", &second)]);
    let contigs = [
        ("most", first[..2_500].to_vec()),
        ("sliver", second[..300].to_vec()),
    ];
    let (reported, summary) = screen(&shard, directory.path(), &contigs, &close_lineage());
    assert_eq!(names(&reported), ["first"]);
    assert_eq!(summary.surveyed_entries, 1);
}

#[test]
fn distant_relative_skipped() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_bases(&mut rng(11), 6_000);
    let relative = with_substitutions(&mut rng(12), &plasmid, 0.96);
    let shard = plasmids(
        directory.path(),
        &[("plasmid", &plasmid), ("relative", &relative)],
    );
    let contigs = [("contig", plasmid.clone())];
    let (reported, summary) = screen(&shard, directory.path(), &contigs, &ScreenConfig::default());
    assert_eq!(reported.len(), 2);
    assert_eq!(summary.surveyed_entries, 2);
    // At 0.99 the seeds that it shares are too few for it to be aligned at all.
    let (reported, summary) = screen(&shard, directory.path(), &contigs, &close_lineage());
    assert_eq!(names(&reported), ["plasmid"]);
    assert_eq!(summary.surveyed_entries, 1);
}

#[test]
fn margin_picks_closest() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_bases(&mut rng(21), 6_000);
    let cousin = with_substitutions(&mut rng(22), &plasmid, 0.996);
    let shard = plasmids(
        directory.path(),
        &[("plasmid", &plasmid), ("cousin", &cousin)],
    );
    let contigs = [("contig", plasmid.clone())];
    let both = close_lineage();
    let reported = screen(&shard, directory.path(), &contigs, &both).0;
    assert_eq!(names(&reported), ["plasmid", "cousin"]);
    let (plasmid, cousin) = (&reported[0].contigs[0], &reported[1].contigs[0]);
    assert_eq!(plasmid.score, plasmid.best_score);
    assert_eq!(cousin.best_score, plasmid.score);
    assert!(
        cousin.score < plasmid.score - 0.001,
        "{cousin:?} {plasmid:?}"
    );
    assert!(plasmid.identity > cousin.identity);
    assert_eq!((plasmid.entries, cousin.entries), (2, 2));
    let margin = ScreenConfig {
        closest_margin: Some(0.001),
        ..close_lineage()
    };
    let reported = screen(&shard, directory.path(), &contigs, &margin).0;
    assert_eq!(names(&reported), ["plasmid"]);
    assert_eq!(reported[0].contigs[0].entries, 1);
    let wide = ScreenConfig {
        closest_margin: Some(0.05),
        ..close_lineage()
    };
    let reported = screen(&shard, directory.path(), &contigs, &wide).0;
    assert_eq!(names(&reported), ["plasmid", "cousin"]);
}

#[test]
fn multimer_keeps_plasmid() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_bases(&mut rng(51), 3_000);
    let double = [&plasmid[..], &plasmid[..]].concat();
    let shard = plasmids(
        directory.path(),
        &[("single", &plasmid), ("double", &double)],
    );
    let contigs = [("contig", plasmid.clone())];
    let config = ScreenConfig {
        closest_margin: Some(0.001),
        ..close_lineage()
    };
    let reported = screen(&shard, directory.path(), &contigs, &config).0;
    assert!(names(&reported).contains(&"single".to_string()));
}

#[test]
fn shared_genes_ignored() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_bases(&mut rng(41), 6_000);
    let shard = plasmids(directory.path(), &[("plasmid", &plasmid)]);
    // Four chromosome contigs of 20 kb, each with a different quarter of the plasmid in it.
    let contigs: Vec<(String, Vec<u8>)> = (0..4)
        .map(|quarter| {
            let gene = &plasmid[quarter * 1_500..(quarter + 1) * 1_500];
            let chromosome = random_bases(&mut rng(42 + quarter as u64), 20_000);
            (
                format!("chromosome{quarter}"),
                [&chromosome[..8_000], gene, &chromosome[8_000..]].concat(),
            )
        })
        .collect();
    let contigs: Vec<(&str, Vec<u8>)> = contigs
        .iter()
        .map(|(name, sequence)| (name.as_str(), sequence.clone()))
        .collect();
    let (reported, summary) = screen(&shard, directory.path(), &contigs, &close_lineage());
    assert!(reported.is_empty());
    assert_eq!(summary.surveyed_entries, 0);
    // The plasmid does come out when the contigs only have to cover it.
    let reported = screen(&shard, directory.path(), &contigs, &ScreenConfig::default()).0;
    assert_eq!(reported.len(), 1);
    assert_eq!(reported[0].contigs.len(), 4);
}

#[test]
fn contig_rule_limits() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_bases(&mut rng(5), 3_000);
    let shard = plasmids(directory.path(), &[("plasmid", &plasmid)]);
    let carrier = [&plasmid[2_000..], &random_bases(&mut rng(8), 9_000)[..]].concat();
    let contigs = [("pure", plasmid[..2_000].to_vec()), ("carrier", carrier)];
    let spans = |report: &ScreenEntry| -> Vec<(String, u64, u64)> {
        report
            .contigs
            .iter()
            .map(|c| (c.name.clone(), c.length, c.aligned))
            .collect()
    };
    let reported = screen(&shard, directory.path(), &contigs, &ScreenConfig::default()).0;
    assert_eq!(reported.len(), 1);
    assert_eq!(
        spans(&reported[0]),
        [
            ("pure".into(), 2_000, 2_000),
            ("carrier".into(), 10_000, 1_000)
        ]
    );
    let strict = ScreenConfig::new(DecisionRule {
        min_fragment_identity: 0.0,
        min_contig_coverage: 0.5,
        ..DecisionRule::default()
    });
    let reported = screen(&shard, directory.path(), &contigs, &strict).0;
    assert_eq!(reported.len(), 1);
    assert_eq!(spans(&reported[0]), [("pure".into(), 2_000, 2_000)]);
    let mut more = strict.clone();
    more.rule.min_query_coverage = QueryCoverage::Bases(2_700);
    assert!(
        screen(&shard, directory.path(), &contigs, &more)
            .0
            .is_empty()
    );
}

#[test]
fn matched_contigs_written() {
    use std::io::Read;
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_bases(&mut rng(61), 3_000);
    let shard = plasmids(directory.path(), &[("plasmid", &plasmid)]);
    let read = |path: &Path| {
        let mut text = String::new();
        flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        text
    };
    let atlas = directory.path().join("atlas");
    let config = ScreenConfig {
        contig_dir: Some(atlas.clone()),
        survey_only: true,
        ..close_lineage()
    };
    let contigs = [
        ("copy", plasmid[200..2_800].to_vec()),
        ("other", random_bases(&mut rng(62), 2_000)),
    ];
    let (reported, summary) = screen(&shard, directory.path(), &contigs, &config);
    assert!(reported.is_empty());
    assert_eq!(summary.alignments, 0);
    assert_eq!((summary.matched_contigs, summary.matched_bases), (1, 2_600));
    let text = read(&atlas.join("assembly.gz"));
    assert!(text.starts_with(">copy extra words\n"), "{text}");
    assert!(!text.contains("other"));
    assert_eq!(
        text.lines().nth(1).unwrap().as_bytes(),
        &plasmid[200..2_800]
    );
    // Without survey only, the file holds the contigs of the reported.
    let hits_only = directory.path().join("reported");
    let config = ScreenConfig {
        contig_dir: Some(hits_only.clone()),
        survey_only: false,
        ..config
    };
    let (reported, summary) = screen(&shard, directory.path(), &contigs, &config);
    assert_eq!(reported.len(), 1);
    assert_eq!(
        (summary.reported_contigs, summary.reported_bases),
        (1, 2_600)
    );
    let text = read(&hits_only.join("assembly.gz"));
    assert!(text.starts_with(">copy extra words\n") && !text.contains("other"));
    // An assembly without a matching contig leaves no file.
    let empty = directory.path().join("empty");
    let config = ScreenConfig {
        contig_dir: Some(empty.clone()),
        ..config
    };
    let unrelated = [("other", random_bases(&mut rng(63), 2_000))];
    let (_, summary) = screen(&shard, directory.path(), &unrelated, &config);
    assert_eq!(summary.matched_contigs, 0);
    assert!(!empty.exists());
}

#[test]
fn assembly_files_read() {
    use std::io::Write;
    let directory = tempfile::tempdir().unwrap();
    let first = random_bases(&mut rng(5), 3_000);
    let shard = plasmids(directory.path(), &[("first", &first)]);
    let fasta = [
        b">gap\n".as_slice(),
        &[b'N'; 5_000],
        b"\n>c1 words\n",
        &first[200..1_800],
        b"\n>short\nACGT\n",
    ]
    .concat();
    let compressed = directory.path().join("assembly.fa.gz");
    let mut encoder = flate2::write::GzEncoder::new(
        std::fs::File::create(&compressed).unwrap(),
        flate2::Compression::fast(),
    );
    encoder.write_all(&fasta).unwrap();
    encoder.finish().unwrap();
    let (reported, summary) = screen_path(&shard, &compressed, &ScreenConfig::default());
    assert_eq!((summary.contigs, summary.skipped_contigs), (1, 2));
    assert_eq!(reported.len(), 1);
    assert_eq!(reported[0].contigs[0].name, "c1");
    let empty = directory.path().join("empty.fa");
    std::fs::write(&empty, b"").unwrap();
    let (reported, summary) = screen_path(&shard, &empty, &ScreenConfig::default());
    assert!(reported.is_empty());
    assert_eq!(summary.contigs, 0);
}

#[test]
fn records_plain_json() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_bases(&mut rng(71), 3_000);
    let shard = plasmids(directory.path(), &[("plasmid", &plasmid)]);
    let contigs = [("copy", reverse_complement(&plasmid))];
    let (reported, summary) = screen(&shard, directory.path(), &contigs, &ScreenConfig::default());
    let report = serde_json::to_value(&reported[0]).unwrap();
    assert_eq!(report["type"], "entry");
    assert_eq!(report["entry"], "plasmid");
    assert_eq!(report["contigs"][0]["name"], "copy");
    let fragment = &report["primary"][0]["fragment"];
    assert_eq!(fragment["contig_index"], 0);
    assert_eq!(fragment["alignment"]["cigar"], "3000=");
    assert_eq!(fragment["alignment"]["strand"], "reverse");
    let summary = serde_json::to_value(&summary).unwrap();
    assert_eq!(
        (summary["type"].clone(), summary["reported_entries"].clone()),
        ("assembly".into(), 1.into())
    );
}

/// Half of an entry lies across the cut of the first two pieces of a long contig (20 bases
/// before it, 480 after), and another entry lies whole in the overlap of the pieces. The
/// overlap keeps the first match whole, and the second is kept once.
#[test]
fn cut_matches_kept() {
    let directory = tempfile::tempdir().unwrap();
    let half = random_bases(&mut rng(0x19), 1_000);
    let whole = random_bases(&mut rng(0x20), 3_000);
    let shard = plasmids(directory.path(), &[("half", &half), ("whole", &whole)]);
    let mut contig = vec![b'N'; 1_000_000];
    contig[499_980..500_480].copy_from_slice(&half[..500]);
    contig[490_000..493_000].copy_from_slice(&whole);
    let config = ScreenConfig::default();
    let (reported, _) = screen(&shard, directory.path(), &[("long", contig)], &config);
    assert_eq!(names(&reported), ["half", "whole"]);
    let (half, whole) = (&reported[0], &reported[1]);
    assert_eq!(half.covered_bases, 500);
    assert_eq!(placed(half)[0].3..placed(half)[0].4, 499_980..500_480);
    assert_eq!((whole.covered_bases, whole.matches), (3_000, 3_000));
    assert_eq!(whole.primary.len() + whole.alternatives.len(), 1);
    assert_eq!(whole.contigs[0].aligned, 3_000);
}
