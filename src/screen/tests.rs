//! End-to-end tests of the screen on small generated entries and assemblies.

use super::*;
use crate::align::{Strand, complement, format_cigar};
use crate::decide::{Fragment, QueryCoverage};
use crate::index::{BuildOptions, EntryInput, Parameters, Shard, build_shard};
use crate::search::SearchControl;
use std::path::Path;
use std::sync::atomic::AtomicBool;

pub(super) fn random_sequence(length: usize, seed: u64) -> Vec<u8> {
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state % 4) as usize]
        })
        .collect()
}

pub(super) fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|&base| complement(base))
        .collect()
}

/// A copy with about one base in `one_in` changed.
pub(super) fn diverged(sequence: &[u8], one_in: u64, mut state: u64) -> Vec<u8> {
    sequence
        .iter()
        .map(|&base| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            if state.is_multiple_of(one_in) {
                match base {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                }
            } else {
                base
            }
        })
        .collect()
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
    let options = BuildOptions {
        parameters: Parameters {
            k: 21,
            window: 20,
            bin_bits: 8,
            ..Parameters::default()
        },
        bucket_bits: None,
        threads: 2,
        temp_dir: directory.to_path_buf(),
        entry_per_record: false,
    };
    let output = directory.join("index");
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
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
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
fn a_planted_plasmid_is_found_with_its_contig_coverage_identity_and_cigar() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_sequence(6_000, 1);
    let other = random_sequence(6_000, 2);
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
        &random_sequence(2_000, 3)[..],
        &copy,
        &random_sequence(2_000, 4),
    ]
    .concat();
    let contigs = [
        ("chromosome", contig),
        ("unrelated", random_sequence(5_000, 5)),
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
fn contigs_in_both_strands_are_placed_on_the_entry() {
    let directory = tempfile::tempdir().unwrap();
    let first = random_sequence(3_000, 5);
    let second = random_sequence(3_000, 6);
    let shard = plasmids(directory.path(), &[("first", &first), ("second", &second)]);
    let contigs = [
        ("forward", first[200..1_800].to_vec()),
        ("unrelated", random_sequence(2_000, 7)),
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
fn fragments_of_a_plasmid_in_several_whole_contigs_count_together() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_sequence(6_000, 31);
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
    let few = [contigs[0].clone(), ("odd", random_sequence(2_000, 32))];
    assert!(
        screen(&shard, directory.path(), &few, &close_lineage())
            .0
            .is_empty()
    );
}

#[test]
fn adjacent_contigs_in_one_batch_are_not_joined_across_the_spacer() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_sequence(3_000, 15);
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
fn a_plasmid_across_the_border_of_two_batches_is_found_whole() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_sequence(4_000, 81);
    let shard = plasmids(directory.path(), &[("plasmid", &plasmid)]);
    // A contig longer than a batch goes in pieces; the plasmid sits over the cut at 5,000.
    let long = [
        &random_sequence(3_000, 82)[..],
        &plasmid,
        &random_sequence(5_000, 83),
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
    let placed = placed(&reported[0]);
    let on_long: Vec<_> = placed.iter().filter(|p| p.2 == "long").collect();
    assert_eq!(on_long.len(), 2, "{placed:?}");
    assert_eq!(
        (on_long[0].0, on_long[0].1, on_long[0].3),
        (0, 2_000, 3_000)
    );
    assert_eq!(
        (on_long[1].0, on_long[1].1, on_long[1].3),
        (2_000, 4_000, 5_000)
    );
    assert_eq!(reported[0].covered_bases, 4_000);
}

#[test]
fn an_entry_with_several_contigs_is_screened_contig_by_contig() {
    let directory = tempfile::tempdir().unwrap();
    let (first, second) = (random_sequence(3_000, 91), random_sequence(2_000, 92));
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
fn an_entry_that_shares_only_a_sliver_is_not_aligned() {
    let directory = tempfile::tempdir().unwrap();
    let first = random_sequence(3_000, 5);
    let second = random_sequence(3_000, 6);
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
fn a_relative_far_below_the_identity_is_not_aligned() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_sequence(6_000, 11);
    let relative = diverged(&plasmid, 16, 12);
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
fn a_contig_belongs_to_the_closest_plasmids_when_a_margin_is_set() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_sequence(6_000, 21);
    let cousin = diverged(&plasmid, 250, 22);
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
fn a_plasmid_is_found_when_the_index_also_holds_its_multimer() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_sequence(3_000, 51);
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
fn contigs_that_share_single_genes_do_not_add_up_to_a_plasmid() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_sequence(6_000, 41);
    let shard = plasmids(directory.path(), &[("plasmid", &plasmid)]);
    // Four chromosome contigs of 20 kb, each with a different quarter of the plasmid in it.
    let contigs: Vec<(String, Vec<u8>)> = (0..4)
        .map(|quarter| {
            let gene = &plasmid[quarter * 1_500..(quarter + 1) * 1_500];
            let chromosome = random_sequence(20_000, 42 + quarter as u64);
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
fn a_contig_with_a_long_rest_counts_only_when_the_contig_rule_allows() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_sequence(3_000, 5);
    let shard = plasmids(directory.path(), &[("plasmid", &plasmid)]);
    let carrier = [&plasmid[2_000..], &random_sequence(9_000, 8)[..]].concat();
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
fn matched_contigs_are_written_for_a_smaller_set_of_assemblies() {
    use std::io::Read;
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_sequence(3_000, 61);
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
        ("other", random_sequence(2_000, 62)),
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
    let unrelated = [("other", random_sequence(2_000, 63))];
    let (_, summary) = screen(&shard, directory.path(), &unrelated, &config);
    assert_eq!(summary.matched_contigs, 0);
    assert!(!empty.exists());
}

#[test]
fn contigs_of_unknown_bases_are_skipped_and_gzip_and_empty_files_are_read() {
    use std::io::Write;
    let directory = tempfile::tempdir().unwrap();
    let first = random_sequence(3_000, 5);
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
fn records_are_plain_json() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_sequence(3_000, 71);
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
