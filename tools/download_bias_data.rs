use anyhow::{bail, Context, Result};
use clap::Parser;
use flate2::read::GzDecoder;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Duration;

const READ_BUFFER_SIZE: usize = 50 * 1024 * 1024;
const WRITE_BUFFER_SIZE: usize = 50 * 1024 * 1024;
const ZSTD_COMPRESSION_LEVEL: i32 = 3;

const PLASMIDSCOPE_URLS: &[(&str, &str)] = &[
    ("PLSDB", "https://plasmidapi.deepomics.org/api/database/files/PLSDB/PLSDB.fasta.tar.gz"),
    ("IMG-PR", "https://plasmidapi.deepomics.org/api/database/files/IMG-PR/IMG-PR.fasta.tar.gz"),
    ("COMPASS", "https://plasmidapi.deepomics.org/api/database/files/COMPASS/COMPASS.fasta.tar.gz"),
    ("mMGE", "https://plasmidapi.deepomics.org/api/database/files/mMGE/mMGE.fasta.tar.gz"),
    ("GenBank", "https://plasmidapi.deepomics.org/api/database/files/GenBank/GenBank.fasta.tar.gz"),
    ("RefSeq", "https://plasmidapi.deepomics.org/api/database/files/RefSeq/RefSeq.fasta.tar.gz"),
    ("ENA", "https://plasmidapi.deepomics.org/api/database/files/ENA/ENA.fasta.tar.gz"),
    ("Kraken2", "https://plasmidapi.deepomics.org/api/database/files/Kraken2/Kraken2.fasta.tar.gz"),
    ("DDBJ", "https://plasmidapi.deepomics.org/api/database/files/DDBJ/DDBJ.fasta.tar.gz"),
    ("TPA", "https://plasmidapi.deepomics.org/api/database/files/TPA/TPA.fasta.tar.gz"),
];

const NCBI_ASSEMBLY_SUMMARY: &str =
    "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt";

#[derive(Parser)]
#[command(name = "download_bias_data")]
#[command(about = "Download plasmid/chromosome data for bias table generation")]
struct Args {
    #[arg(short, long)]
    output: PathBuf,

    #[arg(short, long)]
    temp: Option<PathBuf>,

    #[arg(short = 'j', long, default_value = "16")]
    threads: usize,

    #[arg(short = 'm', long)]
    max_chromosomes: Option<usize>,

    #[arg(long)]
    skip_plasmids: bool,

    #[arg(long)]
    skip_chromosomes: bool,

    #[arg(long)]
    keep_temp: bool,

    #[arg(long, default_value = "3")]
    retries: u32,
}

#[derive(Debug)]
struct DownloadResult {
    name: String,
    success: bool,
    plasmids_filtered: usize,
    error: Option<String>,
}

impl DownloadResult {
    fn ok(name: impl Into<String>, plasmids_filtered: usize) -> Self {
        Self {
            name: name.into(),
            success: true,
            plasmids_filtered,
            error: None,
        }
    }

    fn err(name: impl Into<String>, error: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            success: false,
            plasmids_filtered: 0,
            error: Some(error.into()),
        }
    }
}

#[derive(Debug, Default)]
struct FirstFailure {
    name: String,
    error: String,
}

impl FirstFailure {
    fn set(&mut self, name: impl Into<String>, error: impl Into<String>) {
        if self.name.is_empty() {
            self.name = name.into();
            self.error = error.into();
        }
    }

    fn format_abort_reason(&self) -> String {
        format!("aborted due to failure in '{}': {}", self.name, self.error)
    }
}

#[derive(Debug, Default)]
struct AssemblyReportInfo {
    chromosome_accessions: HashSet<String>,
    plasmid_accessions: HashSet<String>,
    other_accessions: HashSet<String>,
}

#[derive(Debug, Default)]
struct FilterStats {
    chromosomes_kept: usize,
    plasmids_filtered: usize,
}

fn ftp_to_https(url: &str) -> String {
    if url.starts_with("ftp://") {
        url.replacen("ftp://", "https://", 1)
    } else {
        url.to_string()
    }
}

fn build_client() -> Result<reqwest::blocking::Client> {
    reqwest::blocking::Client::builder()
        .timeout(Duration::from_secs(600))
        .connect_timeout(Duration::from_secs(30))
        .gzip(true)
        .build()
        .context("Failed to build HTTP client")
}

fn download_file_with_retry(
    client: &reqwest::blocking::Client,
    url: &str,
    output: &Path,
    pb: Option<&ProgressBar>,
    retries: u32,
) -> Result<()> {
    let https_url = ftp_to_https(url);
    let mut last_error = None;

    for attempt in 0..=retries {
        if attempt > 0 {
            let delay = Duration::from_secs(1 << (attempt - 1));
            std::thread::sleep(delay);
            if let Some(pb) = pb {
                pb.set_message(format!("retry {}/{}...", attempt, retries));
            }
        }

        match download_file_inner(client, &https_url, output, pb) {
            Ok(()) => return Ok(()),
            Err(e) => {
                last_error = Some(e);
                if let Some(pb) = pb {
                    pb.set_position(0);
                }
            }
        }
    }

    Err(last_error.unwrap())
}

fn download_file_inner(
    client: &reqwest::blocking::Client,
    url: &str,
    output: &Path,
    pb: Option<&ProgressBar>,
) -> Result<()> {
    let resp = client.get(url).send().context("Request failed")?;

    if !resp.status().is_success() {
        bail!("HTTP {}: {}", resp.status(), url);
    }

    let total = resp.content_length().unwrap_or(0);
    if let Some(pb) = pb {
        pb.set_length(total);
    }

    let mut reader = BufReader::with_capacity(READ_BUFFER_SIZE, resp);
    let file = File::create(output).context("Failed to create output file")?;
    let mut writer = BufWriter::with_capacity(WRITE_BUFFER_SIZE, file);

    let mut buf = vec![0u8; 64 * 1024];
    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }
        writer.write_all(&buf[..n])?;
        if let Some(pb) = pb {
            pb.inc(n as u64);
        }
    }
    writer.flush()?;

    Ok(())
}

fn download_text_with_retry(
    client: &reqwest::blocking::Client,
    url: &str,
    retries: u32,
) -> Result<String> {
    let https_url = ftp_to_https(url);
    let mut last_error = None;

    for attempt in 0..=retries {
        if attempt > 0 {
            let delay = Duration::from_secs(1 << (attempt - 1));
            std::thread::sleep(delay);
        }

        match download_text_inner(client, &https_url) {
            Ok(text) => return Ok(text),
            Err(e) => last_error = Some(e),
        }
    }

    Err(last_error.unwrap())
}

fn download_text_inner(client: &reqwest::blocking::Client, url: &str) -> Result<String> {
    let resp = client.get(url).send().context("Request failed")?;
    if !resp.status().is_success() {
        bail!("HTTP {}: {}", resp.status(), url);
    }
    resp.text().context("Failed to read response text")
}

fn parse_assembly_report(content: &str) -> AssemblyReportInfo {
    let mut info = AssemblyReportInfo::default();
    let mut col_indices: HashMap<&str, usize> = HashMap::new();

    for line in content.lines() {
        if line.starts_with("# Sequence-Name") {
            let headers: Vec<&str> = line[2..].split('\t').collect();
            for (i, h) in headers.iter().enumerate() {
                col_indices.insert(h.trim(), i);
            }
            continue;
        }

        if line.starts_with('#') || line.trim().is_empty() || col_indices.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();

        let mol_type_idx = col_indices.get("Assigned-Molecule-Location/Type");
        let refseq_idx = col_indices.get("RefSeq-Accn");
        let genbank_idx = col_indices.get("GenBank-Accn");

        let Some(&mol_idx) = mol_type_idx else {
            continue;
        };

        if fields.len() <= mol_idx {
            continue;
        }

        let mol_type = fields[mol_idx].trim().to_lowercase();

        let accession = refseq_idx
            .and_then(|&i| fields.get(i))
            .map(|s| s.trim())
            .filter(|s| !s.is_empty() && *s != "na")
            .or_else(|| {
                genbank_idx
                    .and_then(|&i| fields.get(i))
                    .map(|s| s.trim())
                    .filter(|s| !s.is_empty() && *s != "na")
            });

        let Some(acc) = accession else {
            continue;
        };

        let add_accession = |set: &mut HashSet<String>, acc: &str| {
            set.insert(acc.to_string());
            if let Some(pos) = acc.find('.') {
                set.insert(acc[..pos].to_string());
            }
        };

        match mol_type.as_str() {
            "chromosome" => add_accession(&mut info.chromosome_accessions, acc),
            "plasmid" => add_accession(&mut info.plasmid_accessions, acc),
            _ => add_accession(&mut info.other_accessions, acc),
        }
    }

    info
}

fn filter_fasta_uncompressed(
    input_gz: &Path,
    output: &Path,
    report_info: &AssemblyReportInfo,
) -> Result<FilterStats> {
    let file = File::open(input_gz)?;
    let reader = BufReader::with_capacity(READ_BUFFER_SIZE, GzDecoder::new(file));

    let out_file = File::create(output)?;
    let mut writer = BufWriter::with_capacity(WRITE_BUFFER_SIZE, out_file);

    let mut stats = FilterStats::default();
    let mut is_chromosome = false;

    for line in reader.lines() {
        let line = line?;

        if let Some(header) = line.strip_prefix('>') {
            let seq_id = header.split_whitespace().next().unwrap_or("");
            let base_id = seq_id.split('.').next().unwrap_or("");

            is_chromosome = report_info.chromosome_accessions.contains(seq_id)
                || report_info.chromosome_accessions.contains(base_id);

            if is_chromosome {
                stats.chromosomes_kept += 1;
                writeln!(writer, "{}", line)?;
            } else if report_info.plasmid_accessions.contains(seq_id)
                || report_info.plasmid_accessions.contains(base_id)
            {
                stats.plasmids_filtered += 1;
            }
        } else if is_chromosome {
            writeln!(writer, "{}", line)?;
        }
    }

    writer.flush()?;
    Ok(stats)
}

fn download_plasmidscope_db(
    client: &reqwest::blocking::Client,
    name: &str,
    url: &str,
    temp_dir: &Path,
    output_dir: &Path,
    pb: &ProgressBar,
    retries: u32,
) -> DownloadResult {
    let output_file = output_dir.join(format!("{}.fasta", name));

    if output_file.exists() {
        pb.finish_with_message(format!("{} (cached)", name));
        return DownloadResult::ok(name, 0);
    }

    let work_dir = temp_dir.join(name);
    let _ = fs::create_dir_all(&work_dir);
    let archive = work_dir.join(format!("{}.tar.gz", name));

    let result = (|| -> Result<DownloadResult> {
        pb.set_message(format!("{} downloading...", name));
        download_file_with_retry(client, url, &archive, Some(pb), retries)?;

        pb.set_message(format!("{} extracting...", name));
        pb.set_length(0);
        pb.set_position(0);

        let file = File::open(&archive)?;
        let reader = BufReader::with_capacity(READ_BUFFER_SIZE, file);
        let gz = GzDecoder::new(reader);
        let mut tar = tar::Archive::new(gz);
        tar.unpack(&work_dir)?;

        let _ = fs::remove_file(&archive);

        let fastas: Vec<PathBuf> = walkdir(&work_dir)?
            .into_iter()
            .filter(|p| {
                let ext = p.extension().and_then(|e| e.to_str()).unwrap_or("");
                ext == "fasta" || ext == "fa" || ext == "fna"
            })
            .collect();

        if fastas.is_empty() {
            bail!("no FASTA files in archive");
        }

        pb.set_message(format!("{} concatenating {} files...", name, fastas.len()));

        let out_file = File::create(&output_file)?;
        let mut writer = BufWriter::with_capacity(WRITE_BUFFER_SIZE, out_file);

        for fasta in &fastas {
            let mut file = File::open(fasta)?;
            std::io::copy(&mut file, &mut writer)?;
        }

        writer.flush()?;
        pb.finish_with_message(format!("{}: done", name));
        Ok(DownloadResult::ok(name, 0))
    })();

    let _ = fs::remove_dir_all(&work_dir);

    match result {
        Ok(r) => r,
        Err(e) => {
            pb.finish_with_message(format!("{}: FAILED", name));
            DownloadResult::err(name, e.to_string())
        }
    }
}

fn download_chromosome(
    client: &reqwest::blocking::Client,
    ftp_path: &str,
    temp_dir: &Path,
    output_dir: &Path,
    retries: u32,
) -> DownloadResult {
    let name = ftp_path.rsplit('/').next().unwrap_or("unknown");
    let output_file = output_dir.join(format!("{}.chromosome.fna", name));

    if output_file.exists() {
        return DownloadResult::ok(name, 0);
    }

    let work_dir = temp_dir.join(name);
    let _ = fs::create_dir_all(&work_dir);
    let fasta_file = work_dir.join("genome.fna.gz");

    let result = (|| -> Result<DownloadResult> {
        let report_url = format!("{}/{}_assembly_report.txt", ftp_path, name);
        let report_content = download_text_with_retry(client, &report_url, retries)?;
        let report_info = parse_assembly_report(&report_content);

        if report_info.chromosome_accessions.is_empty() {
            bail!(
                "no chromosome sequences found (found {} plasmids, {} other)",
                report_info.plasmid_accessions.len() / 2,
                report_info.other_accessions.len() / 2
            );
        }

        let fasta_url = format!("{}/{}_genomic.fna.gz", ftp_path, name);
        download_file_with_retry(client, &fasta_url, &fasta_file, None, retries)?;

        let stats = filter_fasta_uncompressed(&fasta_file, &output_file, &report_info)?;

        if stats.chromosomes_kept == 0 {
            let _ = fs::remove_file(&output_file);
            bail!(
                "no chromosome sequences matched (filtered {} plasmids)",
                stats.plasmids_filtered
            );
        }

        Ok(DownloadResult::ok(name, stats.plasmids_filtered))
    })();

    let _ = fs::remove_dir_all(&work_dir);

    match result {
        Ok(r) => r,
        Err(e) => DownloadResult::err(name, e.to_string()),
    }
}

fn walkdir(dir: &Path) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    walkdir_inner(dir, &mut files)?;
    Ok(files)
}

fn walkdir_inner(dir: &Path, files: &mut Vec<PathBuf>) -> Result<()> {
    if dir.is_dir() {
        for entry in fs::read_dir(dir)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_dir() {
                walkdir_inner(&path, files)?;
            } else {
                files.push(path);
            }
        }
    }
    Ok(())
}

fn download_all_plasmidscope(
    temp_dir: &Path,
    output_dir: &Path,
    threads: usize,
    abort_flag: &AtomicBool,
    first_failure: &Mutex<FirstFailure>,
    retries: u32,
    multi: &MultiProgress,
) -> Result<Vec<DownloadResult>> {
    fs::create_dir_all(output_dir)?;
    fs::create_dir_all(temp_dir)?;

    let style = ProgressStyle::default_bar()
        .template(
            "{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {bytes}/{total_bytes} {msg}",
        )
        .unwrap()
        .progress_chars("##-");

    let client = build_client()?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.min(PLASMIDSCOPE_URLS.len()))
        .build()?;

    let results: Vec<DownloadResult> = pool.install(|| {
        PLASMIDSCOPE_URLS
            .par_iter()
            .map(|(name, url)| {
                if abort_flag.load(Ordering::Relaxed) {
                    let reason = first_failure.lock().unwrap().format_abort_reason();
                    return DownloadResult::err(*name, reason);
                }

                let pb = multi.add(ProgressBar::new(0));
                pb.set_style(style.clone());
                pb.set_message(format!("[plasmid] {}", name));

                let result =
                    download_plasmidscope_db(&client, name, url, temp_dir, output_dir, &pb, retries);

                if !result.success {
                    if let Some(ref err) = result.error {
                        first_failure.lock().unwrap().set(*name, err);
                    }
                    abort_flag.store(true, Ordering::Relaxed);
                }

                result
            })
            .collect()
    });

    for r in &results {
        if !r.success {
            bail!(
                "{}: {}",
                r.name,
                r.error.as_deref().unwrap_or("unknown error")
            );
        }
    }

    Ok(results)
}

fn download_all_chromosomes(
    temp_dir: &Path,
    output_dir: &Path,
    threads: usize,
    max_genomes: Option<usize>,
    abort_flag: &AtomicBool,
    first_failure: &Mutex<FirstFailure>,
    retries: u32,
    multi: &MultiProgress,
) -> Result<(Vec<DownloadResult>, usize)> {
    fs::create_dir_all(output_dir)?;
    fs::create_dir_all(temp_dir)?;

    let client = build_client()?;

    let summary_pb = multi.add(ProgressBar::new_spinner());
    summary_pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} {msg}")
            .unwrap(),
    );
    summary_pb.set_message("[chrom] Downloading NCBI assembly summary...");

    let summary_path = temp_dir.join("assembly_summary.txt");
    download_file_with_retry(&client, NCBI_ASSEMBLY_SUMMARY, &summary_path, None, retries)?;
    summary_pb.finish_with_message("[chrom] Assembly summary downloaded");

    let content = fs::read_to_string(&summary_path)?;
    let mut ftp_paths: Vec<String> = Vec::new();

    for line in content.lines() {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() > 19
            && fields[11] == "Complete Genome"
            && fields[10] == "latest"
            && fields[19] != "na"
        {
            ftp_paths.push(fields[19].to_string());
            if max_genomes.is_some_and(|max| ftp_paths.len() >= max) {
                break;
            }
        }
    }

    let pb = multi.add(ProgressBar::new(ftp_paths.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta}) {msg}",
            )
            .unwrap()
            .progress_chars("##-"),
    );
    pb.set_message(format!("[chrom] Downloading {} genomes...", ftp_paths.len()));

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?;

    let total_plasmids_filtered = AtomicUsize::new(0);

    let results: Vec<DownloadResult> = pool.install(|| {
        ftp_paths
            .par_iter()
            .map(|ftp_path| {
                if abort_flag.load(Ordering::Relaxed) {
                    let name = ftp_path.rsplit('/').next().unwrap_or("unknown");
                    let reason = first_failure.lock().unwrap().format_abort_reason();
                    return DownloadResult::err(name, reason);
                }

                let result = download_chromosome(&client, ftp_path, temp_dir, output_dir, retries);

                if result.success {
                    total_plasmids_filtered.fetch_add(result.plasmids_filtered, Ordering::Relaxed);
                } else {
                    if let Some(ref err) = result.error {
                        let name = ftp_path.rsplit('/').next().unwrap_or("unknown");
                        first_failure.lock().unwrap().set(name, err);
                    }
                    abort_flag.store(true, Ordering::Relaxed);
                }

                pb.inc(1);
                result
            })
            .collect()
    });

    pb.finish_with_message("[chrom] Done");

    let plasmids_filtered = total_plasmids_filtered.load(Ordering::Relaxed);

    for r in &results {
        if !r.success {
            bail!(
                "{}: {}",
                r.name,
                r.error.as_deref().unwrap_or("unknown error")
            );
        }
    }

    Ok((results, plasmids_filtered))
}

fn combine_fastas_zstd(
    input_dir: &Path,
    output_file: &Path,
    pattern: &str,
    desc: &str,
) -> Result<(usize, usize)> {
    let files: Vec<PathBuf> = fs::read_dir(input_dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .map(|n| glob_match(pattern, n))
                .unwrap_or(false)
        })
        .collect();

    if files.is_empty() {
        return Ok((0, 0));
    }

    let pb = ProgressBar::new(files.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} {msg}",
            )
            .unwrap()
            .progress_chars("##-"),
    );
    pb.set_message(desc.to_string());

    let out_file = File::create(output_file)?;
    let buf_writer = BufWriter::with_capacity(WRITE_BUFFER_SIZE, out_file);
    let mut writer = zstd::stream::Encoder::new(buf_writer, ZSTD_COMPRESSION_LEVEL)?;

    let mut seqs = 0usize;
    let mut bases = 0usize;

    for file in &files {
        let f = File::open(file)?;
        let reader = BufReader::with_capacity(READ_BUFFER_SIZE, f);
        for line in reader.lines() {
            let line = line?;
            writeln!(writer, "{}", line)?;
            if line.starts_with('>') {
                seqs += 1;
            } else {
                bases += line.trim().len();
            }
        }
        pb.inc(1);
    }

    writer.finish()?;
    pb.finish_with_message(format!("{}: {} seqs", desc, seqs));

    Ok((seqs, bases))
}

fn glob_match(pattern: &str, name: &str) -> bool {
    if let Some(ext) = pattern.strip_prefix("*.") {
        name.ends_with(&format!(".{}", ext))
    } else if pattern.contains('*') {
        let parts: Vec<&str> = pattern.split('*').collect();
        if parts.len() == 2 {
            let prefix = parts[0];
            let suffix = parts[1];
            (prefix.is_empty() || name.starts_with(prefix))
                && (suffix.is_empty() || name.ends_with(suffix))
        } else {
            name == pattern
        }
    } else {
        name == pattern
    }
}

fn fmt_size(bytes: u64) -> String {
    const UNITS: &[&str] = &["B", "KB", "MB", "GB", "TB"];
    let mut size = bytes as f64;
    for unit in UNITS {
        if size < 1024.0 {
            return format!("{:.2} {}", size, unit);
        }
        size /= 1024.0;
    }
    format!("{:.2} PB", size)
}

fn fmt_bases(bases: usize) -> String {
    let b = bases as f64;
    if b >= 1e12 {
        format!("{:.2} Tbp", b / 1e12)
    } else if b >= 1e9 {
        format!("{:.2} Gbp", b / 1e9)
    } else if b >= 1e6 {
        format!("{:.2} Mbp", b / 1e6)
    } else if b >= 1e3 {
        format!("{:.2} Kbp", b / 1e3)
    } else {
        format!("{} bp", bases)
    }
}

fn main() -> Result<()> {
    let args = Args::parse();

    fs::create_dir_all(&args.output)?;

    let temp_base = args
        .temp
        .unwrap_or_else(|| std::env::temp_dir().join(format!("bias_data_{}", std::process::id())));
    fs::create_dir_all(&temp_base)?;

    println!("Output:  {}", args.output.display());
    println!("Temp:    {}", temp_base.display());
    println!("Threads: {}", args.threads);
    println!("Retries: {}\n", args.retries);

    let plasmid_temp = temp_base.join("plasmids");
    let chrom_temp = temp_base.join("chromosomes");

    let abort_flag = Arc::new(AtomicBool::new(false));
    let first_failure = Arc::new(Mutex::new(FirstFailure::default()));
    let multi = MultiProgress::new();

    println!("=== Downloading Data ===");

    let total_plasmids_filtered = std::thread::scope(|s| -> Result<usize> {
        let plasmid_handle = if !args.skip_plasmids {
            let abort = &abort_flag;
            let failure = &first_failure;
            let temp = &plasmid_temp;
            let threads = args.threads;
            let retries = args.retries;
            let mp = &multi;
            Some(s.spawn(move || {
                download_all_plasmidscope(temp, temp, threads, abort, failure, retries, mp)
            }))
        } else {
            None
        };

        let chrom_handle = if !args.skip_chromosomes {
            let abort = &abort_flag;
            let failure = &first_failure;
            let temp = &chrom_temp;
            let threads = args.threads;
            let max = args.max_chromosomes;
            let retries = args.retries;
            let mp = &multi;
            Some(s.spawn(move || {
                download_all_chromosomes(temp, temp, threads, max, abort, failure, retries, mp)
            }))
        } else {
            None
        };

        if let Some(h) = plasmid_handle {
            h.join().expect("plasmid thread panicked")?;
        }

        let filtered = if let Some(h) = chrom_handle {
            let (_, filtered) = h.join().expect("chromosome thread panicked")?;
            filtered
        } else {
            0
        };

        Ok(filtered)
    })?;

    println!();

    let plasmid_out = args.output.join("all_plasmids.fasta.zst");
    let chrom_out = args.output.join("all_chromosomes.fasta.zst");

    let (p_seqs, p_bases) = if !args.skip_plasmids {
        println!("=== Combining Plasmid Files ===");
        combine_fastas_zstd(
            &plasmid_temp,
            &plasmid_out,
            "*.fasta",
            "Combining plasmids",
        )?
    } else {
        (0, 0)
    };

    let (c_seqs, c_bases) = if !args.skip_chromosomes {
        println!("=== Combining Chromosome Files ===");
        combine_fastas_zstd(
            &chrom_temp,
            &chrom_out,
            "*.chromosome.fna",
            "Combining chromosomes",
        )?
    } else {
        (0, 0)
    };

    if !args.keep_temp {
        let _ = fs::remove_dir_all(&temp_base);
    }

    println!("\n=== Summary ===");
    if !args.skip_plasmids && plasmid_out.exists() {
        let size = fs::metadata(&plasmid_out).map(|m| m.len()).unwrap_or(0);
        println!(
            "Plasmids (positive): {} seqs, {}, {}",
            p_seqs,
            fmt_bases(p_bases),
            fmt_size(size)
        );
    }
    if !args.skip_chromosomes && chrom_out.exists() {
        let size = fs::metadata(&chrom_out).map(|m| m.len()).unwrap_or(0);
        println!(
            "Chromosomes (negative): {} seqs, {}, {}",
            c_seqs,
            fmt_bases(c_bases),
            fmt_size(size)
        );
        println!(
            "Plasmids filtered from chromosome set: {} sequences",
            total_plasmids_filtered
        );
    }

    println!(
        "\njam bias -p {} -n {} -o plasmid.bias",
        plasmid_out.display(),
        chrom_out.display()
    );

    Ok(())
}
