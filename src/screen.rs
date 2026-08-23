use crate::format::VERSION;
use crate::provenance;
use crate::query::{QueryEngine, QuerySketch, SampleMatch};
use crate::reader::JamReader;
use anyhow::{Context, Result};
use serde_json::json;
use std::cmp::Ordering;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

pub const CONTIG_HEADER: &str = "schema_version\tassembly\tquery_contig\treference\tshared_hashes\tquery_hashes\treference_hashes\tquery_containment\treference_containment\tretained_query_containment\tretained_reference_containment\tbias_weighted_query_containment\tuniform_hash_e_value\trank\tscore_mode\tbias_table_id";
pub const SUMMARY_HEADER: &str = "schema_version\tassembly\treference\tsupporting_contigs\tshared_hashes_union\treference_hashes\taggregate_reference_containment\tmax_contig_query_containment\tretained_aggregate_reference_containment\tmax_contig_retained_query_containment\tmax_contig_bias_weighted_query_containment\trank\tscore_mode\tbias_table_id";

#[derive(Debug)]
pub struct ScreenConfig {
    pub input: PathBuf,
    pub database: PathBuf,
    pub output: PathBuf,
    pub summary: PathBuf,
    pub metadata: Option<PathBuf>,
    pub assembly_name: Option<String>,
    pub min_shared: u32,
    pub min_query_containment: f64,
    pub min_reference_containment: f64,
    pub top_per_contig: usize,
    pub top_references: usize,
    pub threads: usize,
    pub memory_gb: usize,
    pub force: bool,
}

#[derive(Debug, Clone, Copy)]
pub struct ScreenStats {
    pub contig_count: usize,
    pub contig_rows: usize,
    pub summary_rows: usize,
}

#[derive(Debug)]
struct ContigHit {
    sample_id: u32,
    reference: String,
    shared_hashes: u32,
    query_hashes: usize,
    reference_hashes: u64,
    query_containment: Option<f64>,
    reference_containment: Option<f64>,
    retained_query_containment: Option<f64>,
    retained_reference_containment: Option<f64>,
    bias_weighted_query_containment: Option<f64>,
    uniform_hash_e_value: Option<f64>,
}

#[derive(Debug, Default)]
struct AssemblyEvidence {
    supporting_contigs: usize,
    max_contig_query_containment: f64,
    max_contig_retained_query_containment: f64,
    max_contig_bias_weighted_query_containment: f64,
}

#[derive(Debug)]
struct SummaryHit {
    sample_id: u32,
    reference: String,
    supporting_contigs: usize,
    shared_hashes_union: u64,
    reference_hashes: u64,
    aggregate_reference_containment: Option<f64>,
    max_contig_query_containment: Option<f64>,
    retained_aggregate_reference_containment: Option<f64>,
    max_contig_retained_query_containment: Option<f64>,
    max_contig_bias_weighted_query_containment: Option<f64>,
}

pub fn run(config: &ScreenConfig) -> Result<ScreenStats> {
    validate_config(config)?;
    let started = provenance::unix_time_seconds();
    let engine = QueryEngine::open(&config.database)?;
    let manifest = provenance::load_database_manifest(&config.database)?;
    validate_manifest(manifest.as_ref(), engine.reader())?;

    let has_bias = engine.has_bias_table();
    let score_mode = if has_bias { "bias" } else { "uniform" };
    let bias_table_id = bias_table_id(manifest.as_ref(), &engine);
    let assembly = config
        .assembly_name
        .clone()
        .unwrap_or_else(|| default_assembly_name(&config.input));
    let sketch = QuerySketch::from_fasta(&config.input, engine.reader(), true)
        .map_err(|error| anyhow::anyhow!(error))?;

    let budget_bytes = config.memory_gb.saturating_mul(1024 * 1024 * 1024);
    let sketch_bytes =
        sketch.total_entries().saturating_mul(16) + sketch.sample_count().saturating_mul(128);
    if sketch_bytes > budget_bytes.saturating_mul(7) / 10 {
        anyhow::bail!(
            "query sketch requires approximately {} bytes, exceeding 70% of the {} GiB memory limit; increase --memory or use a larger --fscale database",
            sketch_bytes,
            config.memory_gb
        );
    }

    let mut contig_writer = output_writer(&config.output)?;
    writeln!(contig_writer, "{CONTIG_HEADER}")?;

    let db_sample_count = engine.reader().stats().sample_count as usize;
    let mut assembly_evidence: Vec<AssemblyEvidence> = (0..db_sample_count)
        .map(|_| AssemblyEvidence::default())
        .collect();
    let mut contig_rows = 0usize;

    if sketch.sample_count() > 0 {
        let chunk_size = screen_chunk_size(
            sketch.sample_count(),
            budget_bytes.saturating_mul(7) / 10,
            db_sample_count,
        );
        for chunk_start in (0..sketch.sample_count()).step_by(chunk_size) {
            let chunk_end = (chunk_start + chunk_size).min(sketch.sample_count());
            let results = engine.query_sketch_chunked(&sketch, chunk_start..chunk_end, 0.0);
            for (offset, result) in results.iter().enumerate() {
                let query_idx = chunk_start + offset;
                let mut hits: Vec<ContigHit> = result
                    .matches
                    .iter()
                    .filter_map(|matched| {
                        make_contig_hit(matched, result.query_size, engine.reader(), has_bias)
                    })
                    .filter(|hit| passes_filters(hit, config, has_bias))
                    .collect();

                hits.sort_by(|left, right| compare_contig_hits(left, right, has_bias));
                for hit in &hits {
                    let evidence = &mut assembly_evidence[hit.sample_id as usize];
                    evidence.supporting_contigs += 1;
                    if let Some(value) = hit.query_containment {
                        evidence.max_contig_query_containment =
                            evidence.max_contig_query_containment.max(value);
                    }
                    if let Some(value) = hit.retained_query_containment {
                        evidence.max_contig_retained_query_containment =
                            evidence.max_contig_retained_query_containment.max(value);
                    }
                    if let Some(value) = hit.bias_weighted_query_containment {
                        evidence.max_contig_bias_weighted_query_containment = evidence
                            .max_contig_bias_weighted_query_containment
                            .max(value);
                    }
                }

                for (index, hit) in hits.iter().take(config.top_per_contig).enumerate() {
                    write_contig_row(
                        &mut contig_writer,
                        &assembly,
                        &sketch.sample_names[query_idx],
                        hit,
                        index + 1,
                        score_mode,
                        &bias_table_id,
                    )?;
                    contig_rows += 1;
                }
            }
        }
    }
    contig_writer.flush()?;

    let union_counts = shared_hash_union_counts(&sketch, engine.reader());
    let mut summaries =
        make_summaries(&assembly_evidence, &union_counts, engine.reader(), has_bias);
    summaries.sort_by(|left, right| compare_summary_hits(left, right, has_bias));
    summaries.truncate(config.top_references);

    let mut summary_writer = output_writer(&config.summary)?;
    writeln!(summary_writer, "{SUMMARY_HEADER}")?;
    for (index, hit) in summaries.iter().enumerate() {
        write_summary_row(
            &mut summary_writer,
            &assembly,
            hit,
            index + 1,
            score_mode,
            &bias_table_id,
        )?;
    }
    summary_writer.flush()?;

    let stats = ScreenStats {
        contig_count: sketch.sample_count(),
        contig_rows,
        summary_rows: summaries.len(),
    };
    if let Some(path) = config.metadata.as_deref() {
        write_run_metadata(
            path,
            config,
            &assembly,
            score_mode,
            &bias_table_id,
            manifest.as_ref(),
            stats,
            started,
        )?;
    }
    Ok(stats)
}

fn validate_config(config: &ScreenConfig) -> Result<()> {
    if !config.input.is_file() {
        anyhow::bail!("assembly input is not a file: {}", config.input.display());
    }
    if !config.database.is_file() {
        anyhow::bail!("database is not a file: {}", config.database.display());
    }
    if config.threads == 0 {
        anyhow::bail!("--threads must be > 0");
    }
    if config.memory_gb == 0 {
        anyhow::bail!("--memory must be > 0");
    }
    if config.min_shared == 0 {
        anyhow::bail!("--min-shared must be > 0");
    }
    validate_containment("--min-query-containment", config.min_query_containment)?;
    validate_containment(
        "--min-reference-containment",
        config.min_reference_containment,
    )?;
    if config.top_per_contig == 0 {
        anyhow::bail!("--top-per-contig must be > 0");
    }
    if config.top_references == 0 {
        anyhow::bail!("--top-references must be > 0");
    }

    let mut outputs = vec![&config.output, &config.summary];
    if let Some(metadata) = &config.metadata {
        outputs.push(metadata);
    }
    for (index, path) in outputs.iter().enumerate() {
        if **path == config.input || **path == config.database {
            anyhow::bail!("output path would overwrite an input: {}", path.display());
        }
        if outputs[..index].contains(path) {
            anyhow::bail!("screen output paths must be distinct: {}", path.display());
        }
        if path.exists() {
            if !path.is_file() {
                anyhow::bail!("output path is not a file: {}", path.display());
            }
            if !config.force {
                anyhow::bail!(
                    "output file {} already exists; use --force to overwrite",
                    path.display()
                );
            }
        }
    }
    Ok(())
}

fn validate_containment(name: &str, value: f64) -> Result<()> {
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        anyhow::bail!("{name} must be finite and between 0.0 and 1.0, got {value}");
    }
    Ok(())
}

fn validate_manifest(
    manifest: Option<&provenance::DatabaseManifest>,
    reader: &JamReader,
) -> Result<()> {
    let Some(manifest) = manifest else {
        return Ok(());
    };
    if manifest.database_format_version != VERSION {
        anyhow::bail!(
            "database manifest format version {} does not match readable .jam version {}",
            manifest.database_format_version,
            VERSION
        );
    }
    if manifest.hash_id != provenance::HASH_ID {
        anyhow::bail!(
            "database hash_id '{}' does not match required '{}'",
            manifest.hash_id,
            provenance::HASH_ID
        );
    }
    if manifest.kmer_size != reader.kmer_size()
        || manifest.hash_threshold != reader.threshold()
        || manifest.entropy_threshold.to_bits() != reader.min_entropy().to_bits()
    {
        anyhow::bail!("database manifest parameters do not match the .jam header");
    }
    if manifest.bias.is_some() != reader.has_bias_table() {
        anyhow::bail!("database manifest bias mode does not match the .jam database");
    }
    let stats = reader.stats();
    if manifest.sample_count != stats.sample_count
        || manifest.entry_count != stats.entry_count
        || manifest.unique_hash_count != stats.unique_hash_count
        || manifest.database_size_bytes != stats.file_size
    {
        anyhow::bail!("database manifest counts or file size do not match the .jam database");
    }
    Ok(())
}

fn bias_table_id(manifest: Option<&provenance::DatabaseManifest>, engine: &QueryEngine) -> String {
    if let Some(id) = manifest
        .and_then(|manifest| manifest.bias.as_ref())
        .map(|bias| bias.table_id.clone())
    {
        return id;
    }
    engine
        .bias_table()
        .map(|table| format!("sha256:{}", provenance::sha256_bytes(&table.to_bytes())))
        .unwrap_or_else(|| "NA".to_string())
}

fn default_assembly_name(path: &Path) -> String {
    path.file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("assembly")
        .to_string()
}

fn output_writer(path: &Path) -> Result<BufWriter<File>> {
    let file = File::create(path)
        .with_context(|| format!("failed to create screen output {}", path.display()))?;
    Ok(BufWriter::with_capacity(1024 * 1024, file))
}

fn screen_chunk_size(total: usize, budget_bytes: usize, db_sample_count: usize) -> usize {
    let per_query_bytes = db_sample_count.max(1).saturating_mul(40);
    (budget_bytes / per_query_bytes.max(1)).clamp(1, total)
}

fn make_contig_hit(
    matched: &SampleMatch,
    query_hashes: usize,
    reader: &JamReader,
    has_bias: bool,
) -> Option<ContigHit> {
    let reference_hashes = reader.sample_size(matched.sample_id)?;
    let retained_query = if query_hashes > 0 {
        matched.hit_count as f64 / query_hashes as f64
    } else {
        0.0
    };
    let retained_reference = if reference_hashes > 0 {
        matched.hit_count as f64 / reference_hashes as f64
    } else {
        0.0
    };
    Some(ContigHit {
        sample_id: matched.sample_id,
        reference: reader
            .sample_name(matched.sample_id)
            .unwrap_or("unknown")
            .to_string(),
        shared_hashes: matched.hit_count,
        query_hashes,
        reference_hashes,
        query_containment: (!has_bias).then_some(retained_query),
        reference_containment: (!has_bias).then_some(retained_reference),
        retained_query_containment: has_bias.then_some(retained_query),
        retained_reference_containment: has_bias.then_some(retained_reference),
        bias_weighted_query_containment: has_bias.then_some(matched.containment),
        uniform_hash_e_value: (!has_bias).then_some(matched.e_value),
    })
}

fn passes_filters(hit: &ContigHit, config: &ScreenConfig, has_bias: bool) -> bool {
    let query = if has_bias {
        hit.retained_query_containment.unwrap_or(0.0)
    } else {
        hit.query_containment.unwrap_or(0.0)
    };
    let reference = if has_bias {
        hit.retained_reference_containment.unwrap_or(0.0)
    } else {
        hit.reference_containment.unwrap_or(0.0)
    };
    hit.shared_hashes >= config.min_shared
        && query >= config.min_query_containment
        && reference >= config.min_reference_containment
}

fn compare_contig_hits(left: &ContigHit, right: &ContigHit, has_bias: bool) -> Ordering {
    let left_primary = if has_bias {
        left.bias_weighted_query_containment.unwrap_or(0.0)
    } else {
        left.query_containment.unwrap_or(0.0)
    };
    let right_primary = if has_bias {
        right.bias_weighted_query_containment.unwrap_or(0.0)
    } else {
        right.query_containment.unwrap_or(0.0)
    };
    let left_reference = if has_bias {
        left.retained_reference_containment.unwrap_or(0.0)
    } else {
        left.reference_containment.unwrap_or(0.0)
    };
    let right_reference = if has_bias {
        right.retained_reference_containment.unwrap_or(0.0)
    } else {
        right.reference_containment.unwrap_or(0.0)
    };
    right_primary
        .total_cmp(&left_primary)
        .then_with(|| right.shared_hashes.cmp(&left.shared_hashes))
        .then_with(|| right_reference.total_cmp(&left_reference))
        .then_with(|| left.reference.cmp(&right.reference))
        .then_with(|| left.sample_id.cmp(&right.sample_id))
}

fn write_contig_row(
    writer: &mut impl Write,
    assembly: &str,
    query_contig: &str,
    hit: &ContigHit,
    rank: usize,
    score_mode: &str,
    bias_table_id: &str,
) -> Result<()> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        provenance::OUTPUT_SCHEMA_VERSION,
        tsv_field(assembly),
        tsv_field(query_contig),
        tsv_field(&hit.reference),
        hit.shared_hashes,
        hit.query_hashes,
        hit.reference_hashes,
        format_optional(hit.query_containment),
        format_optional(hit.reference_containment),
        format_optional(hit.retained_query_containment),
        format_optional(hit.retained_reference_containment),
        format_optional(hit.bias_weighted_query_containment),
        format_e_value(hit.uniform_hash_e_value),
        rank,
        score_mode,
        bias_table_id,
    )?;
    Ok(())
}

fn shared_hash_union_counts(sketch: &QuerySketch, reader: &JamReader) -> Vec<u64> {
    let mut counts = vec![0u64; reader.stats().sample_count as usize];
    for bucket_idx in 0..crate::format::BUCKET_COUNT {
        let query = sketch.bucket(bucket_idx);
        let database = reader.bucket_entries(bucket_idx);
        let mut query_index = 0usize;
        let mut database_index = 0usize;
        while query_index < query.len() && database_index < database.len() {
            let query_hash = query[query_index].0;
            while query_index < query.len() && query[query_index].0 == query_hash {
                query_index += 1;
            }
            while database_index < database.len() && database[database_index].hash < query_hash {
                database_index += 1;
            }
            let mut last_sample = None;
            while database_index < database.len() && database[database_index].hash == query_hash {
                let sample_id = database[database_index].sample_id;
                if last_sample != Some(sample_id) {
                    counts[sample_id as usize] += 1;
                    last_sample = Some(sample_id);
                }
                database_index += 1;
            }
        }
        reader.release_bucket(bucket_idx);
    }
    counts
}

fn make_summaries(
    evidence: &[AssemblyEvidence],
    union_counts: &[u64],
    reader: &JamReader,
    has_bias: bool,
) -> Vec<SummaryHit> {
    evidence
        .iter()
        .enumerate()
        .filter(|(_, evidence)| evidence.supporting_contigs > 0)
        .map(|(sample_id, evidence)| {
            let reference_hashes = reader.sample_size(sample_id as u32).unwrap_or(0);
            let union = union_counts.get(sample_id).copied().unwrap_or(0);
            let containment = if reference_hashes > 0 {
                union as f64 / reference_hashes as f64
            } else {
                0.0
            };
            SummaryHit {
                sample_id: sample_id as u32,
                reference: reader
                    .sample_name(sample_id as u32)
                    .unwrap_or("unknown")
                    .to_string(),
                supporting_contigs: evidence.supporting_contigs,
                shared_hashes_union: union,
                reference_hashes,
                aggregate_reference_containment: (!has_bias).then_some(containment),
                max_contig_query_containment: (!has_bias)
                    .then_some(evidence.max_contig_query_containment),
                retained_aggregate_reference_containment: has_bias.then_some(containment),
                max_contig_retained_query_containment: has_bias
                    .then_some(evidence.max_contig_retained_query_containment),
                max_contig_bias_weighted_query_containment: has_bias
                    .then_some(evidence.max_contig_bias_weighted_query_containment),
            }
        })
        .collect()
}

fn compare_summary_hits(left: &SummaryHit, right: &SummaryHit, has_bias: bool) -> Ordering {
    let left_primary = if has_bias {
        left.retained_aggregate_reference_containment.unwrap_or(0.0)
    } else {
        left.aggregate_reference_containment.unwrap_or(0.0)
    };
    let right_primary = if has_bias {
        right
            .retained_aggregate_reference_containment
            .unwrap_or(0.0)
    } else {
        right.aggregate_reference_containment.unwrap_or(0.0)
    };
    let left_query = if has_bias {
        left.max_contig_bias_weighted_query_containment
            .unwrap_or(0.0)
    } else {
        left.max_contig_query_containment.unwrap_or(0.0)
    };
    let right_query = if has_bias {
        right
            .max_contig_bias_weighted_query_containment
            .unwrap_or(0.0)
    } else {
        right.max_contig_query_containment.unwrap_or(0.0)
    };
    right_primary
        .total_cmp(&left_primary)
        .then_with(|| right.shared_hashes_union.cmp(&left.shared_hashes_union))
        .then_with(|| right.supporting_contigs.cmp(&left.supporting_contigs))
        .then_with(|| right_query.total_cmp(&left_query))
        .then_with(|| left.reference.cmp(&right.reference))
        .then_with(|| left.sample_id.cmp(&right.sample_id))
}

fn write_summary_row(
    writer: &mut impl Write,
    assembly: &str,
    hit: &SummaryHit,
    rank: usize,
    score_mode: &str,
    bias_table_id: &str,
) -> Result<()> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        provenance::OUTPUT_SCHEMA_VERSION,
        tsv_field(assembly),
        tsv_field(&hit.reference),
        hit.supporting_contigs,
        hit.shared_hashes_union,
        hit.reference_hashes,
        format_optional(hit.aggregate_reference_containment),
        format_optional(hit.max_contig_query_containment),
        format_optional(hit.retained_aggregate_reference_containment),
        format_optional(hit.max_contig_retained_query_containment),
        format_optional(hit.max_contig_bias_weighted_query_containment),
        rank,
        score_mode,
        bias_table_id,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_run_metadata(
    path: &Path,
    config: &ScreenConfig,
    assembly: &str,
    score_mode: &str,
    bias_table_id: &str,
    manifest: Option<&provenance::DatabaseManifest>,
    stats: ScreenStats,
    started: u64,
) -> Result<()> {
    let database_manifest_path = provenance::sidecar_path(&config.database);
    let database_sha256 = match manifest {
        Some(manifest) => manifest.database_sha256.clone(),
        None => provenance::sha256_file(&config.database)?,
    };
    let database_manifest_sha256 = database_manifest_path
        .exists()
        .then(|| provenance::sha256_file(&database_manifest_path))
        .transpose()?;
    let query = provenance::file_identity(&config.input)?;
    let metadata = json!({
        "schema_version": provenance::OUTPUT_SCHEMA_VERSION,
        "jam_rs_version": env!("CARGO_PKG_VERSION"),
        "source_commit": provenance::source_commit(),
        "source_dirty": provenance::source_dirty(),
        "hash_id": provenance::HASH_ID,
        "hash_zero_policy": provenance::HASH_ZERO_POLICY,
        "score_mode": score_mode,
        "bias_table_id": bias_table_id,
        "assembly": assembly,
        "database": {
            "path": config.database.display().to_string(),
            "sha256": database_sha256,
            "format_version": VERSION,
            "manifest_path": database_manifest_path.exists().then(|| database_manifest_path.display().to_string()),
            "manifest_sha256": database_manifest_sha256,
        },
        "query": query,
        "parameters": {
            "min_shared": config.min_shared,
            "min_query_containment": config.min_query_containment,
            "min_reference_containment": config.min_reference_containment,
            "top_per_contig": config.top_per_contig,
            "top_references": config.top_references,
            "threads": config.threads,
            "memory_gb": config.memory_gb,
        },
        "results": {
            "contig_count": stats.contig_count,
            "contig_candidate_rows": stats.contig_rows,
            "assembly_candidate_rows": stats.summary_rows,
            "contig_output": config.output.display().to_string(),
            "assembly_summary": config.summary.display().to_string(),
        },
        "command": provenance::command_line(),
        "started_unix_seconds": started,
        "completed_unix_seconds": provenance::unix_time_seconds(),
    });
    provenance::write_json(path, &metadata)
}

fn format_optional(value: Option<f64>) -> String {
    value
        .map(|value| format!("{value:.6}"))
        .unwrap_or_else(|| "NA".to_string())
}

fn format_e_value(value: Option<f64>) -> String {
    value
        .map(|value| format!("{value:.6e}"))
        .unwrap_or_else(|| "NA".to_string())
}

fn tsv_field(value: &str) -> String {
    value
        .chars()
        .map(|character| match character {
            '\t' | '\r' | '\n' => ' ',
            other => other,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chunk_size_respects_small_memory_budget() {
        assert_eq!(screen_chunk_size(1000, 4000, 100), 1);
    }

    #[test]
    fn tsv_fields_are_single_line() {
        assert_eq!(tsv_field("a\tb\nc"), "a b c");
    }
}
