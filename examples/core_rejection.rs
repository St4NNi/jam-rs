#[cfg(not(feature = "bench-internals"))]
fn main() {
    panic!("requires bench-internals");
}

#[cfg(feature = "bench-internals")]
fn main() -> anyhow::Result<()> {
    use anyhow::{Context, ensure};
    use jam_rs::shared_reader::{SharedGroup, SharedReader};
    use needletail::{Sequence, parse_fastx_file};
    use rayon::prelude::*;
    use serde_json::json;
    use std::time::Instant;
    use xorf::{BinaryFuse8Ref, Filter, FilterRef};
    let args = std::env::args().collect::<Vec<_>>();
    ensure!(args.len() >= 4, "INDEX OUTPUT MAX_KEYS [QUERY ...]");
    rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build_global()?;
    let started = Instant::now();
    let builder = SharedReader::open(&args[1])?;
    let (bytes, count, end_prefix) = builder.benchmark_core_filter(args[3].parse()?)?;
    let build_ns = started.elapsed().as_nanos();
    ensure!(!bytes.is_empty(), "empty test target");
    let filter = BinaryFuse8Ref::from_dma(&bytes[..20], &bytes[20..]);
    std::fs::write(&args[2], &bytes)?;
    println!(
        "{}",
        json!({"build_ns":build_ns,"covered_cores":count,
        "end_prefix_exclusive":end_prefix,"descriptor_bytes":20,
        "fingerprint_bytes":bytes.len()-20,"serialized_bytes":bytes.len(),
        "descriptor":bytes[..20],"sha256":jam_rs::jidx::sha256(&bytes)})
    );
    drop(builder);
    for query in &args[4..] {
        let mut input = parse_fastx_file(query)?;
        let mut requests = Vec::new();
        let mut records = 0;
        while let Some(record) = input.next() {
            let record = record?;
            let circular = record
                .id()
                .split(|b| b.is_ascii_whitespace())
                .find(|s| s.starts_with(b"topology="))
                .context("topology missing")?
                == b"topology=circular";
            let mut sequence = record.seq().normalize(false).into_owned();
            let length = sequence.len();
            if circular && length >= 15 {
                sequence.extend_from_within(..14);
            }
            requests.extend(
                sequence
                    .bit_kmers(15, true)
                    .take_while(|(position, _, _)| *position < length)
                    .map(|(_, word, _)| word.0 as u32),
            );
            records += 1;
        }
        let run = |enabled: bool, observed: bool| -> anyhow::Result<_> {
            let started = Instant::now();
            let reader = if observed {
                SharedReader::open_observed(&args[1])?
            } else {
                SharedReader::open(&args[1])?
            };
            let mut cores = requests.clone();
            cores.par_sort_unstable();
            cores.dedup();
            let distinct = cores.len();
            let mut ranges = Vec::new();
            let mut start = 0;
            while start < cores.len() {
                let mut end = (start + 32_768).min(cores.len());
                if end < cores.len() {
                    let prefix = cores[end] >> 14;
                    while end > start && cores[end - 1] >> 14 == prefix {
                        end -= 1;
                    }
                    if end == start {
                        end = (start + 32_768).min(cores.len());
                    }
                }
                ranges.push(start..end);
                start = end;
            }
            let mut found = Vec::<SharedGroup>::new();
            let mut survivors = 0;
            for wave in ranges.chunks(8) {
                let chunks = wave
                    .par_iter()
                    .map(|range| {
                        let keys = &cores[range.clone()];
                        let filtered;
                        let keys = if enabled {
                            filtered = keys
                                .iter()
                                .copied()
                                .filter(|core| {
                                    *core >> 14 >= end_prefix || filter.contains(&u64::from(*core))
                                })
                                .collect::<Vec<_>>();
                            &filtered
                        } else {
                            keys
                        };
                        Ok::<_, anyhow::Error>((
                            keys.len(),
                            reader.benchmark_resolve_sorted_cores(keys)?,
                        ))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                for (n, groups) in chunks {
                    survivors += n;
                    found.extend(groups);
                }
            }
            let elapsed_ns = started.elapsed().as_nanos();
            let hits = found
                .iter()
                .map(|g| {
                    (
                        g.key(),
                        g.core_ordinal(),
                        g.member_count(),
                        g.occurrence_count(),
                    )
                })
                .collect::<Vec<_>>();
            Ok((
                json!({"enabled":enabled,"observed":observed,"elapsed_ns":elapsed_ns,
                "distinct":distinct,"survivors":survivors,"hits":hits.len(),"stats":reader.stats()}),
                hits,
            ))
        };
        let (_, expected) = run(false, false)?;
        let mut samples = Vec::new();
        for repetition in 0..5 {
            for enabled in if repetition % 2 == 0 {
                [false, true]
            } else {
                [true, false]
            } {
                let (mut sample, hits) = run(enabled, false)?;
                ensure!(hits == expected, "filter changed exact results");
                sample["repetition"] = json!(repetition);
                samples.push(sample);
            }
        }
        for enabled in [false, true] {
            let (sample, hits) = run(enabled, true)?;
            ensure!(hits == expected, "observer changed results");
            samples.push(sample);
        }
        let covered = requests
            .iter()
            .filter(|core| **core >> 14 < end_prefix)
            .count();
        println!(
            "{}",
            json!({"query":query,"records":records,"occurrences":requests.len(),
            "covered_occurrences":covered,"samples":samples})
        );
    }
    Ok(())
}
