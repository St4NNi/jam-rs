//! End-to-end tests of stage 3 on small generated shards. Alignment windows come from exact
//! k-mer matches between the query and each contig, as stage 2 would find them from seeds.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::align::{Interval, Strand, format_cigar};
use crate::decide::{Fragment, QueryCoverage};
use crate::index::{BuildOptions, EntryInput, Parameters, Shard, build_shard};
use crate::search::ContigWindows;
use crate::seeds::reverse_complement;
use crate::storage::{MappedFile, SliceSource};
use crate::test_support::{no_stop, random_bases, rng, with_substitutions};
use std::collections::HashSet;
use std::ops::Range;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

/// Conserved stretches of `kept` bases between diverged stretches of `diverged` bases.
fn clustered(sequence: &[u8], kept: usize, diverged: usize, per_thousand: u64) -> Vec<u8> {
    let mut copy = Vec::with_capacity(sequence.len());
    for (index, block) in sequence.chunks(kept + diverged).enumerate() {
        copy.extend_from_slice(&block[..kept.min(block.len())]);
        if block.len() > kept {
            let identity = 1.0 - per_thousand as f64 / 1000.0;
            copy.extend(with_substitutions(
                &mut rng(77 + index as u64),
                &block[kept..],
                identity,
            ));
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
        singleton: false,
        sketch_only: false,
        chunks: 1,
    };
    let root = directory.join("shard.jam");
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
        alignment_windows: windows,
    }
}

pub(super) fn rule(coverage: &str) -> DecisionRule {
    DecisionRule {
        min_query_coverage: coverage.parse::<QueryCoverage>().unwrap(),
        min_identity: 0.8,
        min_fragment_identity: 0.0,
        min_fragment_columns: 100,
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
    let control = no_stop();
    let workspace = AlignStageWorkspace::default();
    align_entry(
        shard, &prepared, located, &config, rule, &workspace, control,
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

/// The nine entries of the stage 3 tests, named s1 to s9.
fn fixture_entries(query: &[u8]) -> Vec<EntryContigs> {
    let flank = |length: usize, seed: u64| random_bases(&mut rng(seed), length);
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
                (
                    "high",
                    with_substitutions(&mut rng(6), &query[..1600], 0.95),
                ),
                ("low", clustered(&query[1600..], 100, 300, 370)),
            ],
        ),
        named(
            "s7",
            vec![
                ("p", query[..2000].to_vec()),
                ("q", with_substitutions(&mut rng(7), &query[..2000], 0.99)),
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
    let query = random_bases(&mut rng(11), 3000);
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
fn planted_query_placed() {
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
fn reverse_copy_confirmed() {
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
    // The copy ends in a diverged block, so an alignment end may stop a few dozen bases early.
    let covered = s2.decision.mosaic.covered_bases;
    assert!(covered >= 2900, "{covered}");
}

#[test]
fn split_parts_count() {
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
fn low_entry_unreported() {
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
    // A few flank bases may match the query by chance and extend the alignment.
    let covered = s4.decision.mosaic.covered_bases;
    assert!((400..410).contains(&covered), "{covered}");
}

#[test]
fn mixed_identities_combine() {
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
fn circular_origin_fragment() {
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
fn foreign_contigs_ignored() {
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
fn repeat_contig_excluded() {
    let directory = tempfile::tempdir().unwrap();
    let query = random_bases(&mut rng(21), 3000);
    // Contig 0 holds the query with a few differences; contig 1 holds 600 query bases in
    // 10,000 unrelated ones, at a lower identity.
    let host = [
        random_bases(&mut rng(22), 4000),
        with_substitutions(&mut rng(23), &query[1000..1600], 0.9),
        random_bases(&mut rng(24), 5400),
    ]
    .concat();
    let entries = vec![(
        "assembly".to_string(),
        vec![
            (
                "plasmid".to_string(),
                with_substitutions(&mut rng(25), &query, 0.995),
            ),
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
    // The repeat copy is only an alternative, so it does not change the identity.
    let alternatives = &all.decision.mosaic.alternatives;
    assert!(alternatives.iter().any(|f| f.contig_index == 1));
    assert_eq!(counted.decision.identity(), all.decision.identity());
}

#[test]
fn exact_pieces_chained() {
    // Thirteen contigs of 40 query bases each: the anchors of each start at window positions
    // 0, 4, ..., 24 and score 39, below the chain minimum of 40 bases.
    let directory = tempfile::tempdir().unwrap();
    let query = random_bases(&mut rng(20261023), 1000);
    let contigs = (0..13)
        .map(|i| (format!("c{i}"), query[80 * i..80 * i + 40].to_vec()))
        .collect();
    let shard = build(directory.path(), &[("pieces".to_string(), contigs)]);
    let rule = DecisionRule {
        min_query_coverage: "0.5".parse::<QueryCoverage>().unwrap(),
        min_identity: 0.95,
        ..DecisionRule::default()
    };
    let aligned = align(&shard, 0, &query, false, &rule);
    assert!(aligned.decision.qualifies);
    assert_eq!(aligned.decision.mosaic.covered_bases, 520);
    assert_eq!(fragments(&aligned).len(), 13);
}

#[test]
fn window_widens_edge() {
    let directory = tempfile::tempdir().unwrap();
    let query = random_bases(&mut rng(42), 20_000);
    let contig = [
        random_bases(&mut rng(43), 5000),
        query.clone(),
        random_bases(&mut rng(44), 5000),
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
        alignment_windows: vec![ContigWindows {
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
fn insertion_deletion_sides() {
    let directory = tempfile::tempdir().unwrap();
    let query = random_bases(&mut rng(31), 2000);
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
    file: SliceSource,
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
fn stop_between_windows() {
    use crate::storage::ReadAt;
    let fixture = fixture();
    let s3 = entry("s3");
    let located = located(&fixture.shard, s3, &fixture.query, 300);
    assert_eq!(located.alignment_windows.len(), 2);
    let stop = Arc::new(AtomicBool::new(false));
    let file = MappedFile::open(&fixture.shard.sequence_path(s3).unwrap(), false).unwrap();
    let entries = &fixture.shard.manifest().entries;
    let range = entries[s3 as usize].sequence.as_ref().unwrap().range();
    let source = Arc::new(StopOnRead {
        file: SliceSource::new(Arc::new(file), range).unwrap(),
        armed: AtomicBool::new(false),
        stop: stop.clone(),
        reads: AtomicUsize::new(0),
    });
    let packed = PackedSequence::open(source.clone(), fixture.shard.contig_table(s3).unwrap());
    // The first read checks the header; it happens before the source is armed.
    packed.read_bases(0..1, &mut Vec::new()).unwrap();
    source.armed.store(true, Ordering::Relaxed);
    let config = AlignStageConfig::default();
    let query = PreparedQuery::new(&fixture.query, false, &config).unwrap();
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    // The two windows lie on two contigs, aligned on several threads; one thread keeps their order.
    let one_thread = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let result = one_thread.install(|| {
        align_packed_entry(
            &packed,
            &query,
            &located,
            &config,
            &rule("0.5"),
            &AlignStageWorkspace::default(),
            control,
        )
    });
    assert!(matches!(
        result,
        Err(AlignStageError::Stopped(Stopped::Asked))
    ));
    // The first window was read; the second was not.
    assert_eq!(source.reads.load(Ordering::Relaxed), 1);
    assert!(source.size() > 0);
}

#[test]
fn bad_windows_rejected() {
    let fixture = fixture();
    let mut located = located(&fixture.shard, entry("s1"), &fixture.query, 300);
    let beyond = 0..1_000_000;
    located.alignment_windows[0].ranges = vec![beyond];
    let config = AlignStageConfig::default();
    let query = PreparedQuery::new(&fixture.query, false, &config).unwrap();
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let workspace = AlignStageWorkspace::default();
    let run = |config: &AlignStageConfig, workspace: &AlignStageWorkspace| {
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
        run(&config, &workspace),
        Err(AlignStageError::WindowOutsideContig { .. })
    ));
    let mut other = config;
    other.anchors.a = 17;
    assert!(matches!(
        run(&other, &workspace),
        Err(AlignStageError::InvalidConfig(_))
    ));
}
