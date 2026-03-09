use anyhow::{Context, Result};
use jam_rs::bias::{BiasCreateConfig, CMSConfig, HashBiasTable};
use jam_rs::jamhash_u64;
use needletail::{Sequence, parse_fastx_file};
use std::io::Write;
use std::path::Path;
use std::time::Instant;
use tempfile::NamedTempFile;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

const PLASMID_PATH: &str = "tests/testfiles/plasmidscope_clean_50MB.fasta";
const CHROMO_PATH: &str = "tests/testfiles/chromosomes_small.fasta";
const TRAIN_FRACTION: f64 = 0.7;

// ---------------------------------------------------------------------------
// Parameter sweep types
// ---------------------------------------------------------------------------

struct TrialResult {
    alpha: f32,
    target_fscale: u64,
    negative_fscale: u64,
    pos_retention: f32,
    neg_retention: f32,
    fold_enrichment: f32,
    max_fold_enrichment: f32,
    threshold: i8,
    weight_min: f32,
    weight_max: f32,
    weight_mean: f32,
    weight_std: f32,
    pct_positive: f64,
    pct_above_threshold: f64,
    elapsed_secs: f64,
}

fn run_trial(
    alpha: f32,
    target_fscale: Option<u64>,
    negative_fscale: Option<u64>,
    fscale: u64,
    k: u8,
) -> Result<TrialResult> {
    let config = BiasCreateConfig {
        cms: CMSConfig {
            k,
            fscale,
            ..Default::default()
        },
        alpha,
        target_fscale,
        negative_fscale,
        unseen_fscale: None,
    };

    let plasmid = Path::new(PLASMID_PATH);
    let chromo = Path::new(CHROMO_PATH);

    let start = Instant::now();
    let table = HashBiasTable::create(&[plasmid], &[chromo], &config, None)?;
    let elapsed = start.elapsed().as_secs_f64();

    let (min, max, mean, std, positive, above_threshold) = table.weight_stats();
    let total_cells = config.cms.width * config.cms.depth;

    Ok(TrialResult {
        alpha,
        target_fscale: target_fscale.unwrap_or(0),
        negative_fscale: negative_fscale.unwrap_or(0),
        pos_retention: table.positive_retention,
        neg_retention: table.negative_retention,
        fold_enrichment: table.fold_enrichment(),
        max_fold_enrichment: table.max_fold_enrichment,
        threshold: table.threshold,
        weight_min: min,
        weight_max: max,
        weight_mean: mean,
        weight_std: std,
        pct_positive: positive as f64 / total_cells as f64 * 100.0,
        pct_above_threshold: above_threshold as f64 / total_cells as f64 * 100.0,
        elapsed_secs: elapsed,
    })
}

fn print_header() {
    println!(
        "{:<8} {:<12} {:<12} {:>8} {:>8} {:>8} {:>10} {:>6} {:>8} {:>8} {:>8} {:>8} {:>8} {:>10} {:>8}",
        "alpha", "tgt_fscale", "neg_fscale",
        "pos_ret", "neg_ret", "fold_e", "max_fold",
        "thr", "w_min", "w_max", "w_mean", "w_std",
        "%pos", "%>=thr", "time_s"
    );
    println!("{}", "-".repeat(160));
}

fn print_row(r: &TrialResult) {
    println!(
        "{:<8.2} {:<12} {:<12} {:>7.2}% {:>7.2}% {:>8.2} {:>10.2} {:>6} {:>8.2} {:>8.2} {:>8.4} {:>8.4} {:>7.1}% {:>9.1}% {:>8.1}",
        r.alpha, r.target_fscale, r.negative_fscale,
        r.pos_retention * 100.0, r.neg_retention * 100.0,
        r.fold_enrichment, r.max_fold_enrichment,
        r.threshold,
        r.weight_min, r.weight_max, r.weight_mean, r.weight_std,
        r.pct_positive, r.pct_above_threshold,
        r.elapsed_secs
    );
}

// ---------------------------------------------------------------------------
// FASTA splitting & k-mer retention measurement
// ---------------------------------------------------------------------------

/// One FASTA record stored in memory (header + sequence bytes).
struct FastaRecord {
    header: Vec<u8>,
    seq: Vec<u8>,
}

/// Read all records from a FASTA file into memory.
fn read_fasta_records(path: &Path) -> Result<Vec<FastaRecord>> {
    let mut records = Vec::new();
    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("Failed to open {}", path.display()))?;
    while let Some(rec) = reader.next() {
        let rec = rec.context("Failed to parse record")?;
        records.push(FastaRecord {
            header: rec.id().to_vec(),
            seq: rec.normalize(false).to_vec(),
        });
    }
    Ok(records)
}

/// Write a slice of records to a temporary FASTA file; returns the temp file handle.
fn write_temp_fasta(records: &[FastaRecord]) -> Result<NamedTempFile> {
    let mut tmp = NamedTempFile::new()?;
    for rec in records {
        write!(tmp, ">")?;
        tmp.write_all(&rec.header)?;
        writeln!(tmp)?;
        tmp.write_all(&rec.seq)?;
        writeln!(tmp)?;
    }
    tmp.flush()?;
    Ok(tmp)
}

/// Count k-mers and how many pass the bias filter (hard threshold).
struct RetentionStats {
    total_kmers: u64,
    retained_kmers: u64,
    /// k-mers that even pass the fscale subsampling (i.e., are in sketch space).
    sketch_kmers: u64,
    /// Of the sketch-space k-mers, how many pass the bias filter.
    sketch_retained: u64,
}

impl RetentionStats {
    fn retention_pct(&self) -> f64 {
        if self.sketch_kmers == 0 {
            return 0.0;
        }
        self.sketch_retained as f64 / self.sketch_kmers as f64 * 100.0
    }

    fn overall_pct(&self) -> f64 {
        if self.total_kmers == 0 {
            return 0.0;
        }
        self.retained_kmers as f64 / self.total_kmers as f64 * 100.0
    }
}

/// Measure k-mer retention for a set of FASTA records against a bias table.
fn measure_retention(records: &[FastaRecord], table: &HashBiasTable, k: u8) -> RetentionStats {
    let frac_max = u64::MAX / table.fscale();
    let mut stats = RetentionStats {
        total_kmers: 0,
        retained_kmers: 0,
        sketch_kmers: 0,
        sketch_retained: 0,
    };

    for rec in records {
        if rec.seq.len() < k as usize {
            continue;
        }
        for (_, kmer, _) in rec.seq.bit_kmers(k, true) {
            let hash = jamhash_u64(kmer.0);
            stats.total_kmers += 1;
            if hash < frac_max {
                stats.sketch_kmers += 1;
                if table.passes_filter(hash) {
                    stats.sketch_retained += 1;
                    stats.retained_kmers += 1;
                }
            }
        }
    }
    stats
}

/// Measure per-record retention and return (mean, std, min, max) of per-record retention %.
fn measure_per_record_retention(
    records: &[FastaRecord],
    table: &HashBiasTable,
    k: u8,
) -> (f64, f64, f64, f64) {
    let frac_max = u64::MAX / table.fscale();
    let mut retentions = Vec::with_capacity(records.len());

    for rec in records {
        if rec.seq.len() < k as usize {
            continue;
        }
        let mut sketch = 0u64;
        let mut passed = 0u64;
        for (_, kmer, _) in rec.seq.bit_kmers(k, true) {
            let hash = jamhash_u64(kmer.0);
            if hash < frac_max {
                sketch += 1;
                if table.passes_filter(hash) {
                    passed += 1;
                }
            }
        }
        if sketch > 0 {
            retentions.push(passed as f64 / sketch as f64 * 100.0);
        }
    }

    if retentions.is_empty() {
        return (0.0, 0.0, 0.0, 0.0);
    }

    let n = retentions.len() as f64;
    let mean = retentions.iter().sum::<f64>() / n;
    let var = retentions.iter().map(|r| (r - mean).powi(2)).sum::<f64>() / n;
    let std = var.sqrt();
    let min = retentions.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = retentions.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    (mean, std, min, max)
}

/// Run the generalization test: train/test split on plasmids, measure retention
/// on train plasmids, test (held-out) plasmids, and chromosomes.
fn run_generalization_test(alpha: f32, fscale: u64, k: u8) -> Result<()> {
    let plasmid_path = Path::new(PLASMID_PATH);
    let chromo_path = Path::new(CHROMO_PATH);

    // Load all records
    println!("  Loading FASTA records...");
    let all_plasmids = read_fasta_records(plasmid_path)?;
    let chromosomes = read_fasta_records(chromo_path)?;

    let n_total = all_plasmids.len();
    let n_train = (n_total as f64 * TRAIN_FRACTION).round() as usize;
    let n_test = n_total - n_train;

    println!(
        "  Plasmids: {} total, {} train ({:.0}%), {} test ({:.0}%)",
        n_total,
        n_train,
        TRAIN_FRACTION * 100.0,
        n_test,
        (1.0 - TRAIN_FRACTION) * 100.0
    );
    println!("  Chromosomes: {} records", chromosomes.len());

    let train_plasmids = &all_plasmids[..n_train];
    let test_plasmids = &all_plasmids[n_train..];

    // Write train split to temp file for bias table creation
    println!("  Writing train split to temp file...");
    let train_tmp = write_temp_fasta(train_plasmids)?;

    // Build bias table on train plasmids vs chromosomes
    let config = BiasCreateConfig {
        cms: CMSConfig {
            k,
            fscale,
            ..Default::default()
        },
        alpha,
        target_fscale: None,
        negative_fscale: None,
        unseen_fscale: None,
    };

    println!("  Building bias table (train plasmids vs chromosomes)...");
    let start = Instant::now();
    let table = HashBiasTable::create(&[train_tmp.path()], &[chromo_path], &config, None)?;
    let build_time = start.elapsed().as_secs_f64();
    println!("  Built in {build_time:.1}s");
    println!();

    // Also build a "full" table using all plasmids for comparison
    println!("  Building full bias table (all plasmids vs chromosomes)...");
    let full_table = HashBiasTable::create(&[plasmid_path], &[chromo_path], &config, None)?;

    // --- Aggregate retention ---
    println!();
    println!("  === Aggregate k-mer Retention (train-based bias table) ===");
    println!(
        "  {:.<30} {:>12} {:>12} {:>10} {:>10}",
        "Population", "sketch_kmers", "retained", "sketch_%", "overall_%"
    );
    println!("  {}", "-".repeat(78));

    for (label, records) in [
        ("Train plasmids (seen)", train_plasmids),
        ("Test plasmids (unseen)", test_plasmids),
        ("Chromosomes (negative)", chromosomes.as_slice()),
    ] {
        let s = measure_retention(records, &table, k);
        println!(
            "  {:.<30} {:>12} {:>12} {:>9.2}% {:>9.4}%",
            label,
            s.sketch_kmers,
            s.sketch_retained,
            s.retention_pct(),
            s.overall_pct()
        );
    }

    // --- Per-record retention distribution ---
    println!();
    println!("  === Per-Record Retention Distribution (train-based bias table) ===");
    println!(
        "  {:.<30} {:>10} {:>10} {:>10} {:>10}",
        "Population", "mean_%", "std_%", "min_%", "max_%"
    );
    println!("  {}", "-".repeat(54));

    for (label, records) in [
        ("Train plasmids (seen)", train_plasmids),
        ("Test plasmids (unseen)", test_plasmids),
        ("Chromosomes (negative)", chromosomes.as_slice()),
    ] {
        let (mean, std, min, max) = measure_per_record_retention(records, &table, k);
        println!(
            "  {:.<30} {:>9.2}% {:>9.2}% {:>9.2}% {:>9.2}%",
            label, mean, std, min, max
        );
    }

    // --- Comparison: train-only vs full table ---
    println!();
    println!("  === Comparison: Train-only vs Full Table (on test plasmids) ===");
    let test_train_only = measure_retention(test_plasmids, &table, k);
    let test_full = measure_retention(test_plasmids, &full_table, k);
    let chromo_train_only = measure_retention(&chromosomes, &table, k);
    let chromo_full = measure_retention(&chromosomes, &full_table, k);

    println!(
        "  {:.<40} {:>10} {:>10}",
        "", "train-only", "full"
    );
    println!("  {}", "-".repeat(64));
    println!(
        "  {:.<40} {:>9.2}% {:>9.2}%",
        "Test plasmid retention",
        test_train_only.retention_pct(),
        test_full.retention_pct()
    );
    println!(
        "  {:.<40} {:>9.2}% {:>9.2}%",
        "Chromosome retention",
        chromo_train_only.retention_pct(),
        chromo_full.retention_pct()
    );
    let fold_train = if chromo_train_only.retention_pct() > 0.0 {
        test_train_only.retention_pct() / chromo_train_only.retention_pct()
    } else {
        f64::INFINITY
    };
    let fold_full = if chromo_full.retention_pct() > 0.0 {
        test_full.retention_pct() / chromo_full.retention_pct()
    } else {
        f64::INFINITY
    };
    println!(
        "  {:.<40} {:>9.2}x {:>9.2}x",
        "Fold enrichment (test/chromo)", fold_train, fold_full
    );

    // --- Generalization gap ---
    let train_ret = measure_retention(train_plasmids, &table, k);
    let gap = train_ret.retention_pct() - test_train_only.retention_pct();
    println!();
    println!("  === Generalization Gap ===");
    println!(
        "  Train plasmid retention:  {:.2}%",
        train_ret.retention_pct()
    );
    println!(
        "  Test plasmid retention:   {:.2}%",
        test_train_only.retention_pct()
    );
    println!("  Gap (train - test):       {:.2} pp", gap);
    println!(
        "  Interpretation: {}",
        if gap.abs() < 2.0 {
            "Excellent generalization (gap < 2pp)"
        } else if gap.abs() < 5.0 {
            "Good generalization (gap < 5pp)"
        } else if gap.abs() < 10.0 {
            "Moderate generalization (gap < 10pp)"
        } else {
            "Poor generalization (gap >= 10pp) - possible overfitting"
        }
    );

    Ok(())
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() -> Result<()> {
    let k: u8 = 21;
    let fscale: u64 = 1000;

    // Check files exist
    for p in [PLASMID_PATH, CHROMO_PATH] {
        if !Path::new(p).exists() {
            anyhow::bail!("Missing test file: {p}");
        }
    }

    println!("Bias Table Effectiveness Benchmark");
    println!("===================================");
    println!("  Positive (plasmid):    {PLASMID_PATH}");
    println!("  Negative (chromosome): {CHROMO_PATH}");
    println!("  k={k}, fscale={fscale}");
    println!();

    // --- Hard filter mode: sweep alpha ---
    println!("=== Hard Filter Mode (alpha sweep) ===");
    print_header();

    let alphas = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0];
    for &alpha in &alphas {
        match run_trial(alpha, None, None, fscale, k) {
            Ok(r) => print_row(&r),
            Err(e) => println!("alpha={alpha:.2}: ERROR: {e}"),
        }
    }

    println!();

    // --- Soft filter mode: sweep alpha x target_fscale x negative_fscale ---
    println!("=== Soft Filter Mode (alpha x fscale sweep) ===");
    print_header();

    let soft_alphas = [0.5, 1.0, 2.0, 5.0];
    let target_fscales: Vec<u64> = vec![2000, 5000, 10000];
    let neg_fscale_multipliers: Vec<u64> = vec![2, 5, 10];

    for &alpha in &soft_alphas {
        for &tf in &target_fscales {
            for &mult in &neg_fscale_multipliers {
                let nf = tf * mult;
                match run_trial(alpha, Some(tf), Some(nf), fscale, k) {
                    Ok(r) => print_row(&r),
                    Err(e) => println!("alpha={alpha:.2} tf={tf} nf={nf}: ERROR: {e}"),
                }
            }
        }
    }

    println!();

    // --- Generalization test ---
    println!("=== Generalization Test (train/test split, {:.0}%/{:.0}%) ===",
        TRAIN_FRACTION * 100.0, (1.0 - TRAIN_FRACTION) * 100.0);
    println!();

    for &alpha in &[0.5, 1.0, 2.0, 5.0] {
        println!("--- alpha={alpha:.1} ---");
        run_generalization_test(alpha, fscale, k)?;
        println!();
    }

    // --- Detailed stats for reference ---
    println!("=== Detailed Stats (alpha=1.0, hard filter, all plasmids) ===");
    let config = BiasCreateConfig {
        cms: CMSConfig {
            k,
            fscale,
            ..Default::default()
        },
        alpha: 1.0,
        target_fscale: None,
        negative_fscale: None,
        unseen_fscale: None,
    };

    let table = HashBiasTable::create(
        &[Path::new(PLASMID_PATH)],
        &[Path::new(CHROMO_PATH)],
        &config,
        None,
    )?;
    table.print_stats();

    Ok(())
}
