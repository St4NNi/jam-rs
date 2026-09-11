use anyhow::{Context, Result, ensure};
use jam_rs::shared_reader::{SharedReadStats, SharedReader};
use jam_rs::shared_seed::{SharedKey, context_seed};
use needletail::{Sequence, parse_fastx_file};
use serde::Serialize;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::time::Instant;

const LOOKUP_CHUNK_REQUESTS: usize = 8192;

#[derive(Serialize)]
struct Phase {
    elapsed_ns: u64,
    requests: usize,
    found: usize,
    intentional_chunk_core_repeats: u64,
    stats: SharedReadStats,
}

#[derive(Serialize)]
struct Report {
    query_records: usize,
    linear_records: usize,
    circular_records: usize,
    core_occurrences: u64,
    extraction_ns: u64,
    grouping_ns: u64,
    classification_ns_excluded: u64,
    distinct_cores: usize,
    distinct_context_requests: usize,
    absent_cores: usize,
    present_cores: usize,
    absent: Vec<Phase>,
    present_contexts: Vec<Phase>,
}

fn topology(header: &[u8]) -> Result<bool> {
    let mut tokens = header
        .split(|byte| byte.is_ascii_whitespace())
        .skip(1)
        .filter(|token| token.starts_with(b"topology="));
    let circular = match tokens.next() {
        Some(b"topology=linear") => false,
        Some(b"topology=circular") => true,
        _ => anyhow::bail!("each query header requires topology=linear or topology=circular"),
    };
    ensure!(tokens.next().is_none(), "duplicate query topology token");
    Ok(circular)
}

fn query_keys(path: &Path) -> Result<(Vec<SharedKey>, usize, usize, u64, u64, u64)> {
    let extraction = Instant::now();
    let mut input = parse_fastx_file(path)?;
    let mut keys = Vec::new();
    let mut linear = 0;
    let mut circular_count = 0;
    let mut core_occurrences = 0u64;
    while let Some(record) = input.next() {
        let record = record?;
        let circular = topology(record.id())?;
        linear += usize::from(!circular);
        circular_count += usize::from(circular);
        let query = record.seq().normalize(false).into_owned();
        if query.len() < 15 {
            continue;
        }
        let mut sequence = query.clone();
        if circular {
            sequence.extend_from_slice(&query[..14]);
        }
        for (position, word, orientation) in sequence.bit_kmers(15, true) {
            if position >= query.len() {
                break;
            }
            core_occurrences += 1;
            let context = context_seed(&query, position, word.0 as u32, orientation, circular)
                .context("query core context")?;
            for length in [15, 21, 31] {
                if let Some(key) = context.key(length) {
                    keys.push(key);
                }
            }
        }
    }
    let extraction_ns = extraction.elapsed().as_nanos() as u64;
    let grouping = Instant::now();
    keys.sort_unstable_by_key(|key| (key.core, key.length, key.context));
    keys.dedup();
    let grouping_ns = grouping.elapsed().as_nanos() as u64;
    Ok((
        keys,
        linear,
        circular_count,
        core_occurrences,
        extraction_ns,
        grouping_ns,
    ))
}

fn core_aligned_chunks(keys: &[SharedKey]) -> (Vec<Range<usize>>, u64) {
    let mut chunks = Vec::new();
    let mut repeats = 0;
    let mut start = 0;
    while start < keys.len() {
        let mut end = (start + LOOKUP_CHUNK_REQUESTS).min(keys.len());
        if end < keys.len() && keys[end - 1].core == keys[end].core {
            let boundary = end;
            while end > start && keys[end - 1].core == keys[boundary].core {
                end -= 1;
            }
            if end == start {
                end = boundary;
                repeats += 1;
            }
        }
        chunks.push(start..end);
        start = end;
    }
    (chunks, repeats)
}

fn replay(index: &Path, keys: &[SharedKey]) -> Result<Phase> {
    let reader = SharedReader::open_observed(index)?;
    let (chunks, intentional_chunk_core_repeats) = core_aligned_chunks(keys);
    let start = Instant::now();
    let mut found = 0;
    for chunk in chunks {
        found += reader
            .find_many(&keys[chunk])?
            .into_iter()
            .flatten()
            .count();
    }
    let elapsed_ns = start.elapsed().as_nanos() as u64;
    let stats = reader.stats();
    ensure!(
        stats.member_descriptor_inspections == 0,
        "membership read during lookup replay"
    );
    ensure!(
        stats.references_decoded == 0,
        "reference read during lookup replay"
    );
    ensure!(
        stats.physical_positions_decoded == 0,
        "position read during lookup replay"
    );
    Ok(Phase {
        elapsed_ns,
        requests: keys.len(),
        found,
        intentional_chunk_core_repeats,
        stats,
    })
}

fn main() -> Result<()> {
    let mut args = std::env::args_os().skip(1);
    let index = PathBuf::from(
        args.next()
            .context("usage: shared_lookup INDEX QUERY [REPETITIONS]")?,
    );
    let query = PathBuf::from(
        args.next()
            .context("usage: shared_lookup INDEX QUERY [REPETITIONS]")?,
    );
    let repetitions = args
        .next()
        .map(|value| value.to_string_lossy().parse::<usize>())
        .transpose()?
        .unwrap_or(3);
    ensure!(
        repetitions > 0 && args.next().is_none(),
        "invalid arguments"
    );

    let (keys, linear, circular, core_occurrences, extraction_ns, grouping_ns) =
        query_keys(&query)?;
    let mut cores = keys.iter().map(|key| key.core).collect::<Vec<_>>();
    cores.dedup();
    let classifier = SharedReader::open(&index)?;
    let core_keys = cores
        .iter()
        .copied()
        .map(SharedKey::core)
        .collect::<Vec<_>>();
    let classification = Instant::now();
    let mut absent = Vec::new();
    let mut present = Vec::new();
    for (core_chunk, key_chunk) in cores
        .chunks(LOOKUP_CHUNK_REQUESTS)
        .zip(core_keys.chunks(LOOKUP_CHUNK_REQUESTS))
    {
        for (&core, group) in core_chunk.iter().zip(classifier.find_many(key_chunk)?) {
            if group.is_some() {
                present.push(core);
            } else {
                absent.push(SharedKey::core(core));
            }
        }
    }
    let classification_ns_excluded = classification.elapsed().as_nanos() as u64;
    let present_contexts = keys
        .iter()
        .copied()
        .filter(|key| present.binary_search(&key.core).is_ok())
        .collect::<Vec<_>>();

    let mut absent_phases = Vec::with_capacity(repetitions);
    let mut present_phases = Vec::with_capacity(repetitions);
    for _ in 0..repetitions {
        let phase = replay(&index, &absent)?;
        ensure!(
            phase.found == 0,
            "classified absent core resolved as present"
        );
        absent_phases.push(phase);
        present_phases.push(replay(&index, &present_contexts)?);
    }
    serde_json::to_writer_pretty(
        std::io::stdout(),
        &Report {
            query_records: linear + circular,
            linear_records: linear,
            circular_records: circular,
            core_occurrences,
            extraction_ns,
            grouping_ns,
            classification_ns_excluded,
            distinct_cores: cores.len(),
            distinct_context_requests: keys.len(),
            absent_cores: absent.len(),
            present_cores: present.len(),
            absent: absent_phases,
            present_contexts: present_phases,
        },
    )?;
    println!();
    Ok(())
}
