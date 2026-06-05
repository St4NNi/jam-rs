//! End-to-end tests of stage 3 on small generated shards. Alignment windows come from exact
//! k-mer matches between the query and each contig, as stage 2 would find them from seeds.

use super::*;
use crate::align::{Interval, Strand, complement, format_cigar};
use crate::decide::{Fragment, QueryCoverage};
use crate::index::{BuildOptions, EntryInput, Parameters, build_shard};
use crate::search::ContigWindows;
use std::collections::HashSet;
use std::ops::Range;
use std::path::Path;
use std::sync::atomic::AtomicBool;

pub(super) fn random_sequence(length: usize, mut state: u64) -> Vec<u8> {
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state % 4) as usize]
        })
        .collect()
}

/// Substitutes a base whenever a roll of 1,000 falls below `per_thousand`.
pub(super) fn substitute(sequence: &[u8], per_thousand: u64, mut state: u64) -> Vec<u8> {
    let mut copy = sequence.to_vec();
    for base in &mut copy {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        if state % 1000 < per_thousand {
            *base = match *base {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
    }
    copy
}

pub(super) fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| complement(*base))
        .collect()
}

/// Conserved stretches of `kept` bases between diverged stretches of `diverged` bases.
fn clustered(sequence: &[u8], kept: usize, diverged: usize, per_thousand: u64) -> Vec<u8> {
    let mut copy = Vec::with_capacity(sequence.len());
    for (index, block) in sequence.chunks(kept + diverged).enumerate() {
        copy.extend_from_slice(&block[..kept.min(block.len())]);
        if block.len() > kept {
            copy.extend(substitute(&block[kept..], per_thousand, 77 + index as u64));
        }
    }
    copy
}

/// An entry: its name and its contigs.
pub(super) type EntryContigs = (String, Vec<(String, Vec<u8>)>);

/// Builds a shard with one FASTA file per entry in `directory`.
pub(super) fn build(directory: &Path, entries: &[EntryContigs]) -> Shard {
    let mut inputs = Vec::new();
    for (name, contigs) in entries {
        let mut text = Vec::new();
        for (contig, sequence) in contigs {
            text.extend_from_slice(format!(">{contig}\n").as_bytes());
            text.extend_from_slice(sequence);
            text.push(b'\n');
        }
        let path = directory.join(format!("{name}.fa"));
        std::fs::write(&path, text).unwrap();
        inputs.push(EntryInput {
            name: name.clone(),
            path,
        });
    }
    let options = BuildOptions {
        parameters: Parameters::default(),
        bucket_bits: Some(8),
        threads: 2,
        temp_dir: directory.to_path_buf(),
        entry_per_record: false,
    };
    let root = directory.join("shard");
    build_shard(&inputs, &root, &options).unwrap();
    Shard::open(&root).unwrap()
}

/// Entries of one contig each.
pub(super) fn single_contig_entries(entries: &[(&str, &[u8])]) -> Vec<EntryContigs> {
    entries
        .iter()
        .map(|(name, sequence)| (name.to_string(), vec![("c".to_string(), sequence.to_vec())]))
        .collect()
}

/// The windows stage 2 would give: the 1,024-base bins that hold a 21-mer of the query on
/// either strand, widened by `flank` bases and merged per contig.
pub(super) fn located(shard: &Shard, entry: u32, query: &[u8], flank: u64) -> LocatedEntry {
    const K: usize = 21;
    let query = query.to_ascii_uppercase();
    let mut kmers = HashSet::new();
    for strand in [query.clone(), reverse_complement(&query)] {
        let scanned = [&strand[..], &strand[..(K - 1).min(strand.len())]].concat();
        kmers.extend(scanned.windows(K).map(<[u8]>::to_vec));
    }
    let packed = shard.packed(entry).unwrap();
    let mut windows = Vec::new();
    for contig in 0..packed.contig_count() {
        let range = packed.contig_range(contig).unwrap();
        let length = range.end - range.start;
        let mut bases = Vec::new();
        packed.read(contig, 0..length, &mut bases).unwrap();
        let mut ranges: Vec<Range<u64>> = bases
            .windows(K)
            .enumerate()
            .filter(|(_, kmer)| kmers.contains(*kmer))
            .map(|(at, _)| {
                let bin = at as u64 >> 10 << 10;
                bin.saturating_sub(flank)..(bin + 1024 + flank).min(length)
            })
            .collect();
        window::merge_spans(&mut ranges);
        if !ranges.is_empty() {
            windows.push(ContigWindows {
                contig: contig as u32,
                ranges,
            });
        }
    }
    LocatedEntry {
        entry,
        coverage_bound: query.len() as u64,
        containment_estimate: 1.0,
        windows,
    }
}

pub(super) fn rule(coverage: &str) -> DecisionRule {
    DecisionRule {
        min_query_coverage: coverage.parse::<QueryCoverage>().unwrap(),
        ..DecisionRule::default()
    }
}

/// Runs stage 3 for one located entry.
pub(super) fn align_located(
    shard: &Shard,
    query: &[u8],
    circular: bool,
    located: &LocatedEntry,
    rule: &DecisionRule,
) -> AlignedEntry {
    let config = AlignStageConfig::default();
    let prepared = PreparedQuery::new(query, circular, &config).unwrap();
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let mut workspace = AlignStageWorkspace::default();
    align_entry(
        shard,
        &prepared,
        located,
        &config,
        rule,
        &mut workspace,
        control,
    )
    .unwrap()
}

/// Locates and aligns entry `entry`, with windows of 300 flanking bases.
pub(super) fn align(
    shard: &Shard,
    entry: u32,
    query: &[u8],
    circular: bool,
    rule: &DecisionRule,
) -> AlignedEntry {
    let located = located(shard, entry, query, 300);
    align_located(shard, query, circular, &located, rule)
}

/// Every fragment of the mosaic, primary and alternative.
pub(super) fn fragments(aligned: &AlignedEntry) -> Vec<&Fragment> {
    let mosaic = &aligned.decision.mosaic;
    let primary = mosaic.primary.iter().map(|selected| &selected.fragment);
    primary.chain(mosaic.alternatives.iter()).collect()
}

/// The entries of jam 0.10's search tests, in this order: s1 to s9.
fn fixture_entries(query: &[u8]) -> Vec<EntryContigs> {
    let flank = random_sequence;
    let named = |name: &str, contigs: Vec<(&str, Vec<u8>)>| {
        let contigs = contigs
            .into_iter()
            .map(|(contig, bases)| (contig.to_string(), bases))
            .collect();
        (name.to_string(), contigs)
    };
    vec![
        named(
            "s1",
            vec![(
                "k141_1",
                [flank(400, 1), query.to_vec(), flank(300, 2)].concat(),
            )],
        ),
        named(
            "s2",
            vec![(
                "k141_2",
                reverse_complement(&clustered(query, 500, 500, 250)),
            )],
        ),
        named(
            "s3",
            vec![
                ("a", [query[..1400].to_vec(), flank(200, 3)].concat()),
                ("b", reverse_complement(&query[1500..])),
            ],
        ),
        named(
            "s4",
            vec![(
                "x",
                [flank(2000, 4), query[1000..1400].to_vec(), flank(1000, 5)].concat(),
            )],
        ),
        named(
            "s6",
            vec![
                ("high", substitute(&query[..1600], 50, 6)),
                ("low", clustered(&query[1600..], 100, 300, 370)),
            ],
        ),
        named(
            "s7",
            vec![
                ("p", query[..2000].to_vec()),
                ("q", substitute(&query[..2000], 10, 7)),
                ("r", query[2000..].to_vec()),
            ],
        ),
        named(
            "s8",
            vec![("ring", [&query[2200..], &query[..2200]].concat())],
        ),
        (
            "s9".to_string(),
            query
                .chunks(100)
                .enumerate()
                .map(|(index, piece)| (format!("short{index}"), piece.to_vec()))
                .collect(),
        ),
    ]
}

struct Fixture {
    _directory: tempfile::TempDir,
    shard: Shard,
    query: Vec<u8>,
}

fn fixture() -> Fixture {
    let directory = tempfile::tempdir().unwrap();
    let query = random_sequence(3000, 11);
    let shard = build(directory.path(), &fixture_entries(&query));
    Fixture {
        _directory: directory,
        shard,
        query,
    }
}

/// The entry index of a fixture entry by name.
fn entry(name: &str) -> u32 {
    ["s1", "s2", "s3", "s4", "s6", "s7", "s8", "s9"]
        .iter()
        .position(|entry| *entry == name)
        .unwrap() as u32
}

#[test]
fn a_planted_query_is_confirmed_with_its_location_and_cigar() {
    let fixture = fixture();
    let s1 = align(
        &fixture.shard,
        entry("s1"),
        &fixture.query,
        false,
        &rule("0.5"),
    );
    assert!(s1.decision.qualifies);
    assert_eq!(s1.entry, entry("s1"));
    let only = fragments(&s1);
    assert_eq!(only.len(), 1);
    assert_eq!(only[0].contig_index, 0);
    assert_eq!(only[0].alignment.strand, Strand::Forward);
    assert_eq!(
        only[0].alignment.contig_interval,
        Interval::new(400, 3400).unwrap()
    );
    assert_eq!(only[0].query_segments, [Interval::new(0, 3000).unwrap()]);
    assert_eq!(format_cigar(&only[0].alignment.runs), "3000=");
    assert_eq!(s1.decision.identity(), 1.0);
    assert_eq!(s1.decision.mosaic.covered_bases, 3000);
    assert!(s1.work.windows >= 1 && s1.work.bases_read >= 3000 && s1.work.alignments >= 1);
}

#[test]
fn a_diverged_reverse_copy_is_confirmed_on_the_reverse_strand() {
    let fixture = fixture();
    let s2 = align(
        &fixture.shard,
        entry("s2"),
        &fixture.query,
        false,
        &rule("0.5"),
    );
    assert!(s2.decision.qualifies);
    assert!(
        fragments(&s2)
            .iter()
            .all(|fragment| fragment.alignment.strand == Strand::Reverse)
    );
    let identity = s2.decision.identity();
    assert!((0.84..0.91).contains(&identity), "{identity}");
    assert!(s2.decision.mosaic.covered_bases >= 2950);
}

#[test]
fn split_entries_and_contigs_count_together() {
    let fixture = fixture();
    let s3 = align(
        &fixture.shard,
        entry("s3"),
        &fixture.query,
        false,
        &rule("0.5"),
    );
    let contigs: HashSet<u32> = fragments(&s3).iter().map(|f| f.contig_index).collect();
    assert_eq!(contigs.len(), 2);
    assert!((2890..=2910).contains(&s3.decision.mosaic.covered_bases));
    let s7 = align(
        &fixture.shard,
        entry("s7"),
        &fixture.query,
        false,
        &rule("0.5"),
    );
    assert_eq!(s7.decision.mosaic.covered_bases, 3000);
    assert_eq!(s7.decision.mosaic.primary.len(), 2);
    assert_eq!(s7.decision.mosaic.alternatives.len(), 1);
    let s9 = align(
        &fixture.shard,
        entry("s9"),
        &fixture.query,
        false,
        &rule("1"),
    );
    assert!(s9.decision.qualifies);
    assert_eq!(s9.decision.mosaic.covered_bases, 3000);
    let contigs: HashSet<u32> = fragments(&s9).iter().map(|f| f.contig_index).collect();
    assert_eq!(contigs.len(), 30);
    // At full coverage the gap between the two contigs of s3 leaves it out.
    let s3 = align(
        &fixture.shard,
        entry("s3"),
        &fixture.query,
        false,
        &rule("1"),
    );
    assert!(!s3.decision.qualifies);
}

#[test]
fn an_entry_below_the_limits_is_not_reported() {
    let fixture = fixture();
    // s4 holds 400 query bases.
    let s4 = align(
        &fixture.shard,
        entry("s4"),
        &fixture.query,
        false,
        &rule("0.5"),
    );
    assert!(!s4.decision.qualifies);
    assert_eq!(s4.decision.mosaic.covered_bases, 400);
}

#[test]
fn mixed_identities_qualify_combined() {
    let fixture = fixture();
    let s6 = align(
        &fixture.shard,
        entry("s6"),
        &fixture.query,
        false,
        &rule("0.75"),
    );
    let identities: Vec<f64> = fragments(&s6)
        .iter()
        .map(|fragment| fragment.alignment.identity())
        .collect();
    assert!(
        identities
            .iter()
            .any(|identity| (0.7..0.8).contains(identity)),
        "{identities:?}"
    );
    assert!(s6.decision.qualifies && s6.decision.identity() >= 0.8);
    let per_fragment = DecisionRule {
        min_fragment_identity: 0.8,
        ..rule("0.75")
    };
    let s6 = align(
        &fixture.shard,
        entry("s6"),
        &fixture.query,
        false,
        &per_fragment,
    );
    assert!(!s6.decision.qualifies);
}

#[test]
fn a_window_widens_while_the_alignment_reaches_its_edge() {
    let directory = tempfile::tempdir().unwrap();
    let query = random_sequence(20_000, 42);
    let contig = [
        random_sequence(5000, 43),
        query.clone(),
        random_sequence(5000, 44),
    ]
    .concat();
    let shard = build(
        directory.path(),
        &single_contig_entries(&[("long", &contig)]),
    );
    // Stage 2 found only the first bin of the copy.
    let first_bin = 5120..6144;
    let located = LocatedEntry {
        entry: 0,
        coverage_bound: 20_000,
        containment_estimate: 1.0,
        windows: vec![ContigWindows {
            contig: 0,
            ranges: vec![first_bin],
        }],
    };
    let aligned = align_located(&shard, &query, false, &located, &rule("1"));
    assert!(aligned.decision.qualifies);
    assert!(aligned.work.widenings > 0);
    let only = fragments(&aligned);
    assert_eq!(only.len(), 1);
    assert_eq!(
        only[0].alignment.contig_interval,
        Interval::new(5000, 25_000).unwrap()
    );
    assert_eq!(format_cigar(&only[0].alignment.runs), "20000=");
}
