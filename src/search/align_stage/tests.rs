//! End-to-end tests of stage 3 on small generated shards. Alignment windows come from exact
//! k-mer matches between the query and each contig, as stage 2 would find them from seeds.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::align::{Interval, Strand, format_cigar};
use crate::decide::{Fragment, QueryCoverage};
use crate::index::{BuildOptions, EntryInput, Parameters, Shard, build_shard};
use crate::search::ContigWindows;
use crate::seeds::complement;
use crate::storage::MappedFile;
use std::collections::HashSet;
use std::ops::Range;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

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
        temporary_directory: directory.to_path_buf(),
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
fn a_circular_query_aligns_over_its_origin_in_one_fragment() {
    let fixture = fixture();
    let ring = align(
        &fixture.shard,
        entry("s8"),
        &fixture.query,
        true,
        &rule("0.5"),
    );
    let only = fragments(&ring);
    assert_eq!(only.len(), 1);
    assert_eq!(
        only[0].query_segments,
        [
            Interval::new(2200, 3000).unwrap(),
            Interval::new(0, 2200).unwrap()
        ]
    );
    assert_eq!(format_cigar(&only[0].alignment.runs), "3000=");
    assert_eq!(ring.decision.mosaic.covered_bases, 3000);
    let linear = align(
        &fixture.shard,
        entry("s8"),
        &fixture.query,
        false,
        &rule("0.5"),
    );
    assert_eq!(fragments(&linear).len(), 2);
    assert_eq!(linear.decision.mosaic.covered_bases, 3000);
}

#[test]
fn contigs_that_are_mostly_not_the_query_do_not_count() {
    let fixture = fixture();
    let with_contig_rule = |coverage: f64| DecisionRule {
        min_contig_coverage: coverage,
        ..rule("0.5")
    };
    // The contig of s1 is 3,700 bases long and pairs 3,000 of them (81 percent).
    let s1 = |coverage| {
        align(
            &fixture.shard,
            entry("s1"),
            &fixture.query,
            false,
            &with_contig_rule(coverage),
        )
    };
    assert!(s1(0.5).decision.qualifies);
    assert!(!s1(0.9).decision.qualifies);
}

#[test]
fn a_contig_that_only_shares_a_repeat_is_left_out_by_the_contig_rule() {
    let directory = tempfile::tempdir().unwrap();
    let query = random_sequence(3000, 21);
    // Contig 0 holds the query with a few differences; contig 1 holds 600 query bases in
    // 10,000 unrelated ones, at a lower identity.
    let host = [
        random_sequence(4000, 22),
        substitute(&query[1000..1600], 100, 23),
        random_sequence(5400, 24),
    ]
    .concat();
    let entries = vec![(
        "assembly".to_string(),
        vec![
            ("plasmid".to_string(), substitute(&query, 5, 25)),
            ("chromosome".to_string(), host),
        ],
    )];
    let shard = build(directory.path(), &entries);
    let all = align(&shard, 0, &query, false, &rule("0.5"));
    let contigs: HashSet<u32> = fragments(&all).iter().map(|f| f.contig_index).collect();
    assert_eq!(contigs, HashSet::from([0, 1]));
    let strict = DecisionRule {
        min_contig_coverage: 0.5,
        ..rule("0.5")
    };
    let counted = align(&shard, 0, &query, false, &strict);
    assert!(counted.decision.qualifies);
    assert!(fragments(&counted).iter().all(|f| f.contig_index == 0));
    assert!(counted.decision.identity() > all.decision.identity());
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

#[test]
fn insertions_are_query_bases_and_deletions_are_contig_bases() {
    let directory = tempfile::tempdir().unwrap();
    let query = random_sequence(2000, 31);
    let extra = b"GGGGGGGGGG";
    // The contig lacks ten query bases, and then has ten bases the query lacks.
    let contig = [&query[..800], &query[810..1400], extra, &query[1400..]].concat();
    let reverse = reverse_complement(&contig);
    let entries = single_contig_entries(&[("forward", &contig), ("reverse", &reverse)]);
    let shard = build(directory.path(), &entries);
    for entry in [0, 1] {
        // Query-only bases are not covered: 1,990 of 2,000.
        let aligned = align(&shard, entry, &query, false, &rule("0.99"));
        assert!(aligned.decision.qualifies, "entry {entry}");
        let only = fragments(&aligned);
        assert_eq!(only.len(), 1, "entry {entry}");
        let counts = only[0].alignment.counts();
        assert_eq!((counts.insertions, counts.deletions), (10, 10));
        assert_eq!(aligned.decision.mosaic.covered_bases, 1990);
        assert_eq!(only[0].alignment.contig_interval.len(), 2000);
        let strand = [Strand::Forward, Strand::Reverse][entry as usize];
        assert_eq!(only[0].alignment.strand, strand);
        let cigar = format_cigar(&only[0].alignment.runs);
        assert!(cigar.contains("10I") && cigar.contains("10D"), "{cigar}");
    }
}

/// A packed sequence source that sets the stop flag when it is first read after `armed`.
struct StopOnRead {
    file: MappedFile,
    armed: AtomicBool,
    stop: Arc<AtomicBool>,
    reads: AtomicUsize,
}

impl crate::storage::ReadAt for StopOnRead {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> std::io::Result<()> {
        if self.armed.load(Ordering::Relaxed) {
            self.reads.fetch_add(1, Ordering::Relaxed);
            self.stop.store(true, Ordering::Relaxed);
        }
        self.file.read_exact_at(offset, buffer)
    }

    fn size(&self) -> u64 {
        self.file.size()
    }
}

#[test]
fn a_stop_flag_stops_between_windows() {
    use crate::storage::ReadAt;
    let fixture = fixture();
    let s3 = entry("s3");
    let located = located(&fixture.shard, s3, &fixture.query, 300);
    assert_eq!(located.windows.len(), 2);
    let stop = Arc::new(AtomicBool::new(false));
    let source = Arc::new(StopOnRead {
        file: MappedFile::open(&fixture.shard.packed_path(s3).unwrap(), false).unwrap(),
        armed: AtomicBool::new(false),
        stop: stop.clone(),
        reads: AtomicUsize::new(0),
    });
    let packed = PackedSequence::open(source.clone()).unwrap();
    source.armed.store(true, Ordering::Relaxed);
    let config = AlignStageConfig::default();
    let query = PreparedQuery::new(&fixture.query, false, &config).unwrap();
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let result = align_packed_entry(
        &packed,
        &query,
        &located,
        &config,
        &rule("0.5"),
        &mut AlignStageWorkspace::default(),
        control,
    );
    assert!(matches!(
        result,
        Err(AlignStageError::Stopped(Stopped::Asked))
    ));
    // The first window was read; the second was not.
    assert_eq!(source.reads.load(Ordering::Relaxed), 1);
    assert!(source.size() > 0);
}

#[test]
fn windows_outside_their_contig_and_other_anchor_settings_are_rejected() {
    let fixture = fixture();
    let mut located = located(&fixture.shard, entry("s1"), &fixture.query, 300);
    let beyond = 0..1_000_000;
    located.windows[0].ranges = vec![beyond];
    let config = AlignStageConfig::default();
    let query = PreparedQuery::new(&fixture.query, false, &config).unwrap();
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let mut workspace = AlignStageWorkspace::default();
    let run = |config: &AlignStageConfig, workspace: &mut AlignStageWorkspace| {
        align_entry(
            &fixture.shard,
            &query,
            &located,
            config,
            &rule("0.5"),
            workspace,
            control,
        )
    };
    assert!(matches!(
        run(&config, &mut workspace),
        Err(AlignStageError::WindowOutsideContig { .. })
    ));
    let mut other = config;
    other.anchors.a = 17;
    assert!(matches!(
        run(&other, &mut workspace),
        Err(AlignStageError::InvalidConfig(_))
    ));
}
