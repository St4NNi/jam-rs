//! Frozen synthetic suite for the fragment island study. The cases and seeds were fixed before
//! any island result was inspected. The suite reports every lost, gained, shortened, extended or
//! reassigned interval; it asserts only the rules that must hold for any island outcome. The
//! 40-base 80% identity boundary is tested on island windows in `trace.rs`, because a whole-query
//! case always holds an anchor block that outscores it.

use crate::alignment::{Interval, Strand};
use crate::jidx_writer::{ContigInput, JidxInput, JidxWriter, MetagenomeInput};
use crate::trace::{TraceConfig, TraceEngine, TraceResult};
use crate::trace_islands::RegionPolicy;
use noodles_bgzf::{self as bgzf, gzi};
use std::fs::File;
use std::io::Write;

/// Seed used while the policy is developed.
const DEVELOPMENT_SEEDS: [u64; 1] = [0x0d15_ea5e];
/// Independent seeds, run only after the policy is selected.
const CHECK_SEEDS: [u64; 3] = [0x5eed_c0de_0001, 0x5eed_c0de_0002, 0x5eed_c0de_0003];

const LONG_CONTIG: usize = 200_000;
const SHORT_CONTIG: usize = 50_000;
/// Offset of the duplicated mobile element inside the long contig.
const ELEMENT_COPIES: [usize; 2] = [150_000, 180_000];
const ELEMENT_BASES: usize = 1_200;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn dna(&mut self, length: usize) -> Vec<u8> {
        (0..length)
            .map(|_| b"ACGT"[(self.next() & 3) as usize])
            .collect()
    }
}

fn substitute(base: u8) -> u8 {
    match base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    }
}

/// Substitutes every `period`-th base starting at `phase`; no exact run reaches `period` bases.
fn periodic_mutation(sequence: &[u8], period: usize, phase: usize) -> Vec<u8> {
    let mut output = sequence.to_vec();
    for position in (phase..output.len()).step_by(period) {
        output[position] = substitute(output[position]);
    }
    output
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|&base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

struct Case {
    name: &'static str,
    query: Vec<u8>,
    circular: bool,
}

/// Exact blocks of 30 bases hold one full 16-base minimizer window, so each yields an anchor but
/// cannot pass the 40-base acceptance boundary alone.
const ANCHOR_BLOCK: usize = 30;

fn plant(query: &mut [u8], at: usize, piece: &[u8]) {
    query[at..at + piece.len()].copy_from_slice(piece);
}

fn cases(seed: u64, long: &[u8], short: &[u8]) -> Vec<Case> {
    let mut rng = Rng(seed | 1);
    let mut cases = Vec::new();
    let mut background = |length| rng.dna(length);

    // 1. Two accidental collinear anchors separated by an unrelated 30 kb interval.
    let mut query = background(45_000);
    plant(&mut query, 5_000, &long[65_000..65_000 + ANCHOR_BLOCK]);
    plant(&mut query, 35_000, &long[95_000..95_000 + ANCHOR_BLOCK]);
    cases.push(Case {
        name: "accidental_collinear_pair",
        query,
        circular: false,
    });

    // 2. One real 400-base fragment at the first anchor; the second anchor is unrelated.
    let mut query = background(45_000);
    plant(&mut query, 5_000, &long[65_000..65_400]);
    plant(&mut query, 35_000, &long[95_000..95_000 + ANCHOR_BLOCK]);
    cases.push(Case {
        name: "short_fragment_near_one_anchor",
        query,
        circular: false,
    });

    // 3. Two real separate 300-base fragments with unrelated sequence between them.
    let mut query = background(45_000);
    plant(&mut query, 5_000, &long[65_000..65_300]);
    plant(&mut query, 35_000, &long[95_000..95_300]);
    cases.push(Case {
        name: "two_separate_fragments",
        query,
        circular: false,
    });

    // 4. A true 12 kb fragment at 6/7 identity with exact anchor blocks only at its ends.
    let mut query = background(30_000);
    let mut fragment = periodic_mutation(&long[70_000..82_000], 7, 3);
    fragment[..ANCHOR_BLOCK].copy_from_slice(&long[70_000..70_000 + ANCHOR_BLOCK]);
    let tail = fragment.len() - ANCHOR_BLOCK;
    fragment[tail..].copy_from_slice(&long[82_000 - ANCHOR_BLOCK..82_000]);
    plant(&mut query, 9_000, &fragment);
    cases.push(Case {
        name: "long_fragment_few_anchors",
        query: query.clone(),
        circular: false,
    });

    // 5. Same fragment on the reverse strand.
    cases.push(Case {
        name: "long_fragment_few_anchors_reverse",
        query: reverse_complement(&query),
        circular: false,
    });

    // 6. A qualifying 600-base internal fragment without anchors, away from both anchors.
    let mut query = background(45_000);
    plant(&mut query, 5_000, &long[65_000..65_000 + ANCHOR_BLOCK]);
    plant(
        &mut query,
        20_000,
        &periodic_mutation(&long[80_000..80_600], 7, 3),
    );
    plant(&mut query, 35_000, &long[95_000..95_000 + ANCHOR_BLOCK]);
    cases.push(Case {
        name: "internal_fragment_without_anchors",
        query,
        circular: false,
    });

    // 7. A 6 kb fragment at 6/7 identity with a 40-base deletion in the query at its middle.
    let mut query = background(30_000);
    let mut fragment = periodic_mutation(&long[100_000..106_000], 7, 2);
    fragment[..ANCHOR_BLOCK].copy_from_slice(&long[100_000..100_000 + ANCHOR_BLOCK]);
    let tail = fragment.len() - ANCHOR_BLOCK;
    fragment[tail..].copy_from_slice(&long[106_000 - ANCHOR_BLOCK..106_000]);
    fragment.drain(3_000..3_040);
    plant(&mut query, 8_000, &fragment);
    cases.push(Case {
        name: "long_deletion_inside_fragment",
        query,
        circular: false,
    });

    // 8. Repeated mobile element: the query holds the element twice, 30 kb apart, collinear with
    // the two target copies, so each copy also pairs with the other on distant diagonals.
    let element = &long[ELEMENT_COPIES[0]..ELEMENT_COPIES[0] + ELEMENT_BASES];
    let mut query = background(45_000);
    plant(&mut query, 5_000, element);
    plant(&mut query, 35_000, element);
    cases.push(Case {
        name: "repeated_element_collinear_copies",
        query,
        circular: false,
    });

    // 9. Circular query whose 8 kb weak fragment crosses the origin, anchors 3 kb apart per side.
    let mut fragment = periodic_mutation(&long[40_000..48_000], 7, 1);
    for at in [0, 3_000, 5_000, 8_000 - ANCHOR_BLOCK] {
        fragment[at..at + ANCHOR_BLOCK]
            .copy_from_slice(&long[40_000 + at..40_000 + at + ANCHOR_BLOCK]);
    }
    let mut query = background(20_000);
    plant(&mut query, 16_000, &fragment[..4_000]);
    plant(&mut query, 0, &fragment[4_000..]);
    cases.push(Case {
        name: "circular_origin_crossing_fragment",
        query,
        circular: true,
    });

    // 10. Accidental collinear pair on a contig short enough for the whole-contig window.
    let mut query = background(45_000);
    plant(&mut query, 5_000, &short[5_000..5_000 + ANCHOR_BLOCK]);
    plant(&mut query, 35_000, &short[35_000..35_000 + ANCHOR_BLOCK]);
    cases.push(Case {
        name: "short_contig_pair",
        query,
        circular: false,
    });

    // 11. Weak fragment running off the query end, with an ambiguous run and a sparse anchor.
    let mut query = background(20_000);
    let mut fragment = periodic_mutation(&long[120_000..128_000], 7, 4);
    for at in [0, 5_000] {
        fragment[at..at + ANCHOR_BLOCK]
            .copy_from_slice(&long[120_000 + at..120_000 + at + ANCHOR_BLOCK]);
    }
    fragment[2_000..2_050].fill(b'N');
    let kept = query.len() - 14_000;
    plant(&mut query, 14_000, &fragment[..kept]);
    cases.push(Case {
        name: "missing_end_with_ambiguity",
        query,
        circular: false,
    });

    cases
}

fn write_target(directory: &std::path::Path, contigs: &[(&str, &[u8])]) -> std::path::PathBuf {
    let bgzf_path = directory.join("target.bgz");
    let mut raw = Vec::new();
    let mut offsets = Vec::new();
    for (name, sequence) in contigs {
        raw.extend_from_slice(format!(">{name}\n").as_bytes());
        offsets.push(raw.len() as u64);
        for line in sequence.chunks(80) {
            raw.extend_from_slice(line);
            raw.push(b'\n');
        }
    }
    let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
    let mut blocks = Vec::new();
    for (ordinal, chunk) in raw.chunks(32_000).enumerate() {
        if ordinal > 0 {
            blocks.push((writer.position(), (ordinal * 32_000) as u64));
        }
        writer.write_all(chunk).unwrap();
        writer.flush().unwrap();
    }
    writer.finish().unwrap();
    let gzi_path = directory.join("target.gzi");
    gzi::fs::write(&gzi_path, &gzi::Index::from(blocks)).unwrap();
    let bytes = std::fs::read(&bgzf_path).unwrap();
    let reference = directory.join("reference.jidx");
    let mut jidx = JidxWriter::new(
        &reference,
        &JidxInput {
            k: 15,
            rescue_k15: false,
            minimizer_window: 16,
            jam_sha256: [1; 32],
            manifest_sha256: [2; 32],
        },
    )
    .unwrap();
    jidx.begin_metagenome(MetagenomeInput {
        name: "target".into(),
        bgzf_uri: bgzf_path.to_str().unwrap().to_owned(),
        bgzf_bytes: bytes.len() as u64,
        bgzf_sha256: crate::jidx::sha256(&bytes),
        gzi: std::fs::read(gzi_path).unwrap(),
    })
    .unwrap();
    for ((name, sequence), offset) in contigs.iter().zip(offsets) {
        jidx.begin_contig(ContigInput {
            name: (*name).into(),
            length: sequence.len() as u64,
            fasta_offset: offset,
            line_bases: 80,
            line_width: 81,
        })
        .unwrap();
    }
    jidx.finish().unwrap();
    let shared = directory.join("target.shared");
    crate::shared_writer::build_shared_index(&reference, &shared, 16).unwrap();
    shared
}

#[derive(Clone, Debug, PartialEq)]
struct Placed {
    contig_id: u32,
    strand: Strand,
    query: Vec<Interval>,
    target: Interval,
}

fn placed(result: &TraceResult) -> Vec<Placed> {
    let mut output = result
        .metagenomes
        .iter()
        .flat_map(|metagenome| {
            metagenome
                .mosaic
                .primary
                .iter()
                .map(|selected| &selected.fragment)
                .chain(&metagenome.mosaic.alternatives)
        })
        .map(|fragment| Placed {
            contig_id: fragment.contig_id,
            strand: fragment.alignment.strand,
            query: fragment.query_segments.clone(),
            target: fragment.alignment.target_interval,
        })
        .collect::<Vec<_>>();
    output.sort_by_key(|item| (item.query[0].start, item.target.start));
    output
}

fn overlaps(first: &[Interval], second: &[Interval]) -> bool {
    first
        .iter()
        .any(|a| second.iter().any(|b| a.start < b.end && b.start < a.end))
}

/// Classifies each reference interval against the candidate: equal, shortened, extended,
/// changed (both), reassigned (same query, other target) or lost; unmatched candidate
/// intervals are gained.
fn classify(reference: &[Placed], candidate: &[Placed]) -> serde_json::Value {
    let mut used = vec![false; candidate.len()];
    let mut rows = Vec::new();
    for item in reference {
        let bases = |segments: &[Interval]| segments.iter().map(|s| s.end - s.start).sum::<u64>();
        let matched = candidate
            .iter()
            .enumerate()
            .find(|(index, other)| !used[*index] && overlaps(&item.query, &other.query));
        let class = match matched {
            Some((index, other)) => {
                used[index] = true;
                if other == item {
                    "equal"
                } else if other.contig_id != item.contig_id
                    || other.strand != item.strand
                    || !(other.target.start < item.target.end
                        && item.target.start < other.target.end)
                {
                    "reassigned"
                } else if bases(&other.query) < bases(&item.query) {
                    "shortened"
                } else if bases(&other.query) > bases(&item.query) {
                    "extended"
                } else {
                    "changed"
                }
            }
            None => "lost",
        };
        rows.push(serde_json::json!({"class": class, "query": item.query, "target": item.target}));
    }
    for (index, item) in candidate.iter().enumerate() {
        if !used[index] {
            rows.push(
                serde_json::json!({"class": "gained", "query": item.query, "target": item.target}),
            );
        }
    }
    serde_json::Value::Array(rows)
}

fn run_suite(seed: u64) -> Vec<serde_json::Value> {
    let mut rng = Rng(seed.rotate_left(17) | 1);
    let mut long = rng.dna(LONG_CONTIG);
    let element = long[ELEMENT_COPIES[0]..ELEMENT_COPIES[0] + ELEMENT_BASES].to_vec();
    long[ELEMENT_COPIES[1]..ELEMENT_COPIES[1] + ELEMENT_BASES].copy_from_slice(&element);
    let short = rng.dna(SHORT_CONTIG);
    let directory = tempfile::tempdir().unwrap();
    let shared = write_target(directory.path(), &[("long", &long), ("short", &short)]);
    let mut rows = Vec::new();
    for case in cases(seed, &long, &short) {
        let config = TraceConfig {
            use_sketch: false,
            circular: case.circular,
            ..TraceConfig::default()
        };
        let run = |policy| {
            let mut engine = TraceEngine::open_shared(&shared, None).unwrap();
            engine.observed = true;
            engine.region_policy = policy;
            let result = engine.search(case.name, &case.query, config).unwrap();
            (result, engine.batch_stats())
        };
        let (reference, reference_stats) = run(RegionPolicy::Parent);
        let (candidate, candidate_stats) = run(RegionPolicy::Islands);
        let (reference_placed, candidate_placed) = (placed(&reference), placed(&candidate));
        let intervals = classify(&reference_placed, &candidate_placed);
        let classes = intervals
            .as_array()
            .unwrap()
            .iter()
            .map(|row| row["class"].as_str().unwrap())
            .collect::<Vec<_>>();
        let tasks = (
            reference_stats.alignment_tasks,
            candidate_stats.alignment_tasks,
        );
        // A policy that planned exactly the parent tasks must return exactly the parent result.
        if tasks.0 == tasks.1
            && reference_stats.alignment_work.local_cells
                == candidate_stats.alignment_work.local_cells
        {
            assert_eq!(reference_placed, candidate_placed, "{}", case.name);
        }
        rows.push(serde_json::json!({
            "seed": seed,
            "case": case.name,
            "reference_tasks": tasks.0,
            "candidate_tasks": tasks.1,
            "reference_local_cells": reference_stats.alignment_work.local_cells,
            "candidate_local_cells": candidate_stats.alignment_work.local_cells,
            "reference_recomputed_cells": reference_stats.alignment_work.local_recomputed_cells,
            "candidate_recomputed_cells": candidate_stats.alignment_work.local_recomputed_cells,
            "reference_target_bytes": reference_stats.alignment_target_bytes,
            "candidate_target_bytes": candidate_stats.alignment_target_bytes,
            "all_equal": classes.iter().all(|class| *class == "equal"),
            "intervals": intervals,
        }));
    }
    rows
}

#[test]
fn island_suite_development_seed_reports_interval_changes() {
    for seed in DEVELOPMENT_SEEDS {
        for row in run_suite(seed) {
            println!("{row}");
        }
    }
}

#[test]
#[ignore = "independent seeds run once after the island policy is selected"]
fn island_suite_check_seeds_report_interval_changes() {
    for seed in CHECK_SEEDS {
        for row in run_suite(seed) {
            println!("{row}");
        }
    }
}
