//! End-to-end tests of the screen on small generated entries and assemblies.

use super::*;
use crate::align::{Strand, complement, format_cigar};
use crate::decide::Fragment;
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
