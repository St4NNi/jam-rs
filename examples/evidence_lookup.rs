use anyhow::{Context, Result, ensure};
use jam_rs::jidx::RESCUE_K15_TAG;
use jam_rs::jidx_builder::{JidxBuildConfig, build_local_jidx};
use jam_rs::jidx_reader::JidxReader;
use jam_rs::jidx_writer::{ContigInput, JidxInput, JidxWriter, MetagenomeInput, SelectedSeed};
use needletail::{Sequence, parse_fastx_file};
use serde_json::json;
use std::collections::BTreeMap;
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

fn probe(index: &Path, query: &Path, output: &Path, circular_query: bool) -> Result<()> {
    let reader = JidxReader::open(index)?;
    let mut queries = parse_fastx_file(query)?;
    let mut output = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)?,
    );
    while let Some(record) = queries.next() {
        let record = record?;
        let sequence = record.seq().normalize(false).into_owned();
        let mut keys = BTreeMap::<u64, Vec<(usize, bool)>>::new();
        let mut words = [0u64; 2];
        for (scheme, k) in [(0, reader.header().k), (1, 15)] {
            if (scheme == 1 && !reader.header().rescue_k15) || sequence.len() < k as usize {
                continue;
            }
            let mut circular = sequence.clone();
            if circular_query {
                circular.extend_from_slice(&sequence[..k as usize - 1]);
            }
            for (position, word, strand) in circular.bit_kmers(k, true) {
                if position >= sequence.len() {
                    break;
                }
                let key = word.0 | if scheme == 1 { RESCUE_K15_TAG } else { 0 };
                keys.entry(key).or_default().push((position, strand));
                words[scheme] += 1;
            }
        }
        let began = Instant::now();
        let mut present = [0u64; 2];
        let mut memberships = [0u64; 2];
        let mut positions = [0u64; 2];
        let mut placements = Vec::new();
        for (key, query_positions) in &keys {
            let scheme = usize::from(key & RESCUE_K15_TAG != 0);
            let Some(seed) = reader.find_seed(*key)? else {
                continue;
            };
            present[scheme] += 1;
            let documents = reader.seed_documents(seed)?;
            ensure!(documents.len() == seed.document_frequency as usize);
            memberships[scheme] += documents.len() as u64;
            for document in documents {
                let sample = reader
                    .metagenome(document.metagenome_id)?
                    .context("sample")?;
                let occurrences = reader.seed_document_occurrences(seed, document)?;
                ensure!(occurrences.len() as u64 == document.occurrence_count);
                positions[scheme] += occurrences.len() as u64;
                for occurrence in occurrences {
                    let contig = reader.contig(occurrence.contig_id)?.context("contig")?;
                    placements.push(json!({
                        "k": if scheme == 1 { 15 } else { reader.header().k },
                        "sample": sample.name, "contig": contig.name,
                        "target_position": occurrence.position,
                        "target_canonical": occurrence.canonical_orientation,
                        "query_positions": query_positions,
                    }));
                }
            }
        }
        let result = json!({
            "query_id": std::str::from_utf8(record.id())?,
            "query_words_k21_k15": words, "dictionary_probes": keys.len(),
            "present_keys_k21_k15": present,
            "decoded_memberships_k21_k15": memberships,
            "decoded_positions_k21_k15": positions,
            "probe_seconds": began.elapsed().as_secs_f64(), "placements": placements,
            "circular": circular_query,
            "semantics": "exhaustive dense lookup audit; separate from native search timing",
        });
        serde_json::to_writer(&mut output, &result)?;
        writeln!(output)?;
    }
    output.flush()?;
    Ok(())
}

fn prefix(root: &Path) -> Result<()> {
    let mut rows = Vec::new();
    for unrelated in [0u32, 127, 1023] {
        let path = root.join(format!("prefix-{unrelated}.jidx"));
        let mut writer = JidxWriter::new(
            &path,
            &JidxInput {
                k: 21,
                rescue_k15: false,
                minimizer_window: 16,
                jam_sha256: [1; 32],
                manifest_sha256: [2; 32],
            },
        )?;
        for id in 0..=unrelated {
            writer.begin_metagenome(MetagenomeInput {
                name: if id == unrelated {
                    "target".into()
                } else {
                    format!("prefix-{id:04}")
                },
                bgzf_uri: "unused.bgz".into(),
                bgzf_bytes: 100,
                bgzf_sha256: [3; 32],
                gzi: vec![0; 8],
            })?;
            let contig = writer.begin_contig(ContigInput {
                name: "contig".into(),
                length: 100,
                fasta_offset: 0,
                line_bases: 100,
                line_width: 101,
            })?;
            writer.add_seeds(
                contig,
                &[SelectedSeed {
                    packed_key: 7,
                    position: 17,
                    canonical_orientation: true,
                }],
            )?;
        }
        let stats = writer.finish()?;
        let reader = JidxReader::open(path)?;
        reader.verify_checksum()?;
        ensure!(reader.metagenome(unrelated)?.context("target")?.name == "target");
        let mut absent_ns = Vec::new();
        let mut count_ns = Vec::new();
        let mut late_ns = Vec::new();
        for _ in 0..25 {
            let began = Instant::now();
            ensure!(reader.find_seed(8)?.is_none());
            absent_ns.push(began.elapsed().as_nanos());
            let began = Instant::now();
            let seed = reader.find_seed(7)?.context("present seed")?;
            ensure!(seed.document_frequency == unrelated + 1);
            count_ns.push(began.elapsed().as_nanos());
            let began = Instant::now();
            let documents = reader.seed_documents(seed)?;
            let document = *documents
                .iter()
                .find(|d| d.metagenome_id == unrelated)
                .context("late member")?;
            let found = reader.seed_document_occurrences(seed, document)?;
            ensure!(found.len() == 1 && found[0].position == 17 && found[0].canonical_orientation);
            late_ns.push(began.elapsed().as_nanos());
        }
        rows.push(json!({
            "unrelated_prefix_members": unrelated, "file_bytes": stats.file_bytes,
            "requested_member": "target", "requested_contig": "contig",
            "returned_positions": [17], "members_decoded_for_late_answer": unrelated + 1,
            "absent_nanoseconds": absent_ns, "document_count_nanoseconds": count_ns,
            "late_member_nanoseconds": late_ns,
            "bounded_late_access": unrelated == 0,
            "total_occurrence_count_in_key_record": false,
        }));
    }
    println!("{}", serde_json::to_string_pretty(&rows)?);
    Ok(())
}

fn fixture(root: &Path) -> Result<()> {
    let root = root.canonicalize()?;
    let mut state = 982_451_653u64;
    let sequence = (0..5000)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(state >> 62) as usize]
        })
        .collect::<Vec<_>>();
    let mut sources = Vec::new();
    let mut targets = Vec::new();
    for i in 0..2 {
        let name = format!("sample{i}");
        let mut main = if i == 0 {
            sequence.clone()
        } else {
            sequence.reverse_complement()
        };
        main[2500..2530].fill(b'N');
        let contigs = [
            ("main", main),
            ("short15", sequence[..15].to_vec()),
            ("short20", sequence[80..100].to_vec()),
            ("repeat", b"ACGT".repeat(100)),
        ];
        let mut raw = Vec::new();
        let mut fai = String::new();
        for (contig, bases) in &contigs {
            writeln!(raw, ">{contig}")?;
            let offset = raw.len();
            raw.extend_from_slice(bases);
            raw.push(b'\n');
            fai.push_str(&format!(
                "{contig}\t{}\t{offset}\t{}\t{}\n",
                bases.len(),
                bases.len(),
                bases.len() + 1
            ));
        }
        let source = root.join(format!("{name}.fa"));
        std::fs::write(&source, &raw)?;
        let bgzf = root.join(format!("{name}.bgz"));
        let fai_path = root.join(format!("{name}.bgz.fai"));
        let gzi = root.join(format!("{name}.bgz.gzi"));
        let mut writer = noodles_bgzf::io::Writer::new(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&bgzf)?,
        );
        writer.write_all(&raw)?;
        writer.finish()?;
        std::fs::write(&fai_path, fai)?;
        noodles_bgzf::gzi::fs::write(&gzi, &noodles_bgzf::gzi::Index::default())?;
        sources.push(source);
        targets
            .push(json!({"name": format!("{name}.fa"), "bgzf": bgzf, "fai": fai_path, "gzi": gzi}));
    }
    let database = root.join("database.jam");
    jam_rs::writer::build(
        &sources,
        &database,
        &jam_rs::writer::BuildConfig {
            kmer_size: 21,
            fscale: 100,
            singleton: false,
            memory: 1,
            ..jam_rs::writer::BuildConfig::default()
        },
    )?;
    let manifest = root.join("manifest.json");
    std::fs::write(
        &manifest,
        serde_json::to_vec(&json!({"metagenomes": targets}))?,
    )?;
    for minimizer_window in [16, 32, 64] {
        build_local_jidx(
            &database,
            &manifest,
            root.join(format!("w{minimizer_window}.jidx")),
            JidxBuildConfig {
                k: 21,
                rescue_k15: true,
                minimizer_window,
            },
        )?;
    }
    let query = &sequence[321..721];
    std::fs::write(
        root.join("query.fa"),
        format!(
            ">forward\n{}\n>reverse\n{}\n",
            std::str::from_utf8(query)?,
            std::str::from_utf8(&query.reverse_complement())?
        ),
    )?;
    Ok(())
}

fn main() -> Result<()> {
    let args = std::env::args().collect::<Vec<_>>();
    match args.get(1).map(String::as_str) {
        Some("probe" | "probe-linear") if args.len() == 5 => probe(
            Path::new(&args[2]),
            Path::new(&args[3]),
            Path::new(&args[4]),
            args[1] == "probe",
        ),
        Some("prefix") if args.len() == 3 => prefix(Path::new(&args[2])),
        Some("fixture") if args.len() == 3 => fixture(Path::new(&args[2])),
        _ => anyhow::bail!(
            "usage: evidence_lookup probe[-linear] INDEX QUERY OUTPUT | prefix NEW_OUTPUT_DIRECTORY"
        ),
    }
}
