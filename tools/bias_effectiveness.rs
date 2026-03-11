use anyhow::{Context, Result};
use jam_rs::bias::{BiasCreateConfig, CMSConfig, HashBiasTable};
use jam_rs::jamhash_u64;
use needletail::{Sequence, parse_fastx_file};
use rayon::prelude::*;
use std::path::Path;
use std::time::Instant;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

const PLASMID_TRAIN: &str = "tests/testfiles/plasmidscope_clean_50MB_train.fasta";
const PLASMID_TEST: &str = "tests/testfiles/plasmidscope_clean_50MB_test.fasta";
const CHROMO_TRAIN: &str = "tests/testfiles/all_chromosomes_1GB_train.fasta";
const CHROMO_TEST: &str = "tests/testfiles/all_chromosomes_1GB_test.fasta";

const ALPHA: f32 = 1.0;
const K: u8 = 21;
const FSCALE: u64 = 100;
const WIDTH: usize = 1 << 23; // 8M
const DEPTH: usize = 7;

struct FastaRecord {
    seq: Vec<u8>,
}

fn read_fasta_records(path: &Path) -> Result<Vec<FastaRecord>> {
    let mut records = Vec::new();
    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("Failed to open {}", path.display()))?;
    while let Some(rec) = reader.next() {
        let rec = rec.context("Failed to parse record")?;
        records.push(FastaRecord {
            seq: rec.normalize(false).to_vec(),
        });
    }
    Ok(records)
}

struct RetentionStats {
    sketch_kmers: u64,
    sketch_retained: u64,
}

impl RetentionStats {
    fn retention_pct(&self) -> f64 {
        if self.sketch_kmers == 0 {
            return 0.0;
        }
        self.sketch_retained as f64 / self.sketch_kmers as f64 * 100.0
    }

    fn effective_fscale(&self, base_fscale: u64) -> f64 {
        if self.sketch_retained == 0 {
            return f64::INFINITY;
        }
        base_fscale as f64 * self.sketch_kmers as f64 / self.sketch_retained as f64
    }
}

fn measure_retention(records: &[FastaRecord], table: &HashBiasTable, k: u8) -> RetentionStats {
    let frac_max = u64::MAX / table.fscale();
    let (sketch_kmers, sketch_retained) = records
        .par_iter()
        .filter(|rec| rec.seq.len() >= k as usize)
        .map(|rec| {
            let mut sk = 0u64;
            let mut sr = 0u64;
            for (_, kmer, _) in rec.seq.bit_kmers(k, true) {
                let hash = jamhash_u64(kmer.0);
                if hash < frac_max {
                    sk += 1;
                    if table.passes_filter(hash) {
                        sr += 1;
                    }
                }
            }
            (sk, sr)
        })
        .reduce(|| (0, 0), |(a1, b1), (a2, b2)| (a1 + a2, b1 + b2));

    RetentionStats { sketch_kmers, sketch_retained }
}

fn measure_weight_distribution(
    records: &[FastaRecord],
    table: &HashBiasTable,
    k: u8,
) -> [u64; 255] {
    let frac_max = u64::MAX / table.fscale();
    let histograms: Vec<[u64; 255]> = records
        .par_iter()
        .filter(|rec| rec.seq.len() >= k as usize)
        .map(|rec| {
            let mut hist = [0u64; 255];
            for (_, kmer, _) in rec.seq.bit_kmers(k, true) {
                let hash = jamhash_u64(kmer.0);
                if hash < frac_max {
                    let w = table.weight(hash);
                    let idx = (w as i16 + 127) as usize;
                    hist[idx] += 1;
                }
            }
            hist
        })
        .collect();

    let mut total = [0u64; 255];
    for h in &histograms {
        for i in 0..255 {
            total[i] += h[i];
        }
    }
    total
}

fn print_weight_distribution(label: &str, hist: &[u64; 255], table: &HashBiasTable) {
    let total: u64 = hist.iter().sum();
    if total == 0 {
        println!("  {label}: no k-mers");
        return;
    }

    println!("  Weight distribution for: {label}");
    println!(
        "    {:>7} {:>7} {:>12} {:>8} {:>12}",
        "weight", "index", "count", "pct", "eff_fscale"
    );

    let ranges = [
        (-127i8, -101i8, "[-127,-101]"),
        (-100, -51, "[-100, -51]"),
        (-50, -11, "[ -50, -11]"),
        (-10, -1, "[ -10,  -1]"),
        (0, 0, "[   0,   0]"),
        (1, 10, "[   1,  10]"),
        (11, 50, "[  11,  50]"),
        (51, 100, "[  51, 100]"),
        (101, 127, "[ 101, 127]"),
    ];

    for (lo, hi, label_range) in &ranges {
        let lo_idx = (*lo as i16 + 127) as usize;
        let hi_idx = (*hi as i16 + 127) as usize;
        let count: u64 = hist[lo_idx..=hi_idx].iter().sum();
        let pct = count as f64 / total as f64 * 100.0;
        let mid = ((*lo as i16 + *hi as i16) / 2) as i8;
        let fs = table.effective_fscale_at(mid);
        println!(
            "    {:<11} {:>5}-{:<5} {:>12} {:>7.2}% {:>12.0}",
            label_range, lo_idx, hi_idx, count, pct, fs
        );
    }
    println!(
        "    {:.<11} {:>13} {:>12}",
        "TOTAL", "", total
    );
    println!();
}

fn print_oracle_analysis(
    test_pos_hist: &[u64; 255],
    test_neg_hist: &[u64; 255],
    table: &HashBiasTable,
    base_fscale: u64,
    target_fscale: u64,
    negative_fscale: u64,
) {
    let pos_total: u64 = test_pos_hist.iter().sum();
    let neg_total: u64 = test_neg_hist.iter().sum();
    if pos_total == 0 || neg_total == 0 {
        return;
    }

    let ranges: &[(i8, i8, &str)] = &[
        (-127, -101, "[-127,-101]"),
        (-100, -51, "[-100, -51]"),
        (-50, -11, "[ -50, -11]"),
        (-10, -1, "[ -10,  -1]"),
        (0, 0, "[   0,   0]"),
        (1, 10, "[   1,  10]"),
        (11, 50, "[  11,  50]"),
        (51, 100, "[  51, 100]"),
        (101, 127, "[ 101, 127]"),
    ];

    println!("  === Per-Bin Analysis (test/unseen data) ===");
    println!(
        "    {:>11} {:>10} {:>10} {:>8} {:>8} {:>8} {:>10} {:>10}",
        "weight", "pos_cnt", "neg_cnt", "pos_%", "neg_%", "ratio", "cur_fscale", "ideal_fs"
    );

    // Compute per-bin probabilities on test data
    let pt = pos_total as f64;
    let nt = neg_total as f64;
    let base = base_fscale as f64;
    let target_ret = base / target_fscale as f64;

    // For each range, compute enrichment ratio and "ideal" fscale
    struct BinInfo {
        label: String,
        pos_count: u64,
        neg_count: u64,
        pos_frac: f64,
        neg_frac: f64,
        ratio: f64,
        cur_fscale: f64,
    }

    let mut bins: Vec<BinInfo> = Vec::new();

    for &(lo, hi, label) in ranges {
        let lo_idx = (lo as i16 + 127) as usize;
        let hi_idx = (hi as i16 + 127) as usize;
        let pos_count: u64 = test_pos_hist[lo_idx..=hi_idx].iter().sum();
        let neg_count: u64 = test_neg_hist[lo_idx..=hi_idx].iter().sum();
        let pos_frac = pos_count as f64 / pt;
        let neg_frac = neg_count as f64 / nt;
        let ratio = if neg_frac > 1e-12 { pos_frac / neg_frac } else { f64::INFINITY };
        let mid = (lo as i16 + hi as i16) / 2;
        let cur_fscale = table.effective_fscale_at(mid as i8);

        bins.push(BinInfo {
            label: label.to_string(), pos_count, neg_count, pos_frac, neg_frac, ratio, cur_fscale,
        });
    }

    // Compute oracle LUT via greedy: sort bins by pos/neg ratio descending,
    // assign max retention (base_fscale) to highest-ratio bins until we exhaust
    // the target retention budget on a uniform prior, then assign min retention.
    // Use fine-grained per-index bins for the oracle.
    let mut idx_bins: Vec<(usize, f64, f64)> = (0..255)
        .map(|i| {
            let pf = test_pos_hist[i] as f64 / pt;
            let nf = test_neg_hist[i] as f64 / nt;
            (i, pf, nf)
        })
        .collect();

    // Sort by pos/neg ratio descending (highest enrichment first)
    idx_bins.sort_by(|a, b| {
        let ra = if a.2 > 1e-15 { a.1 / a.2 } else { f64::INFINITY };
        let rb = if b.2 > 1e-15 { b.1 / b.2 } else { f64::INFINITY };
        rb.partial_cmp(&ra).unwrap()
    });

    // Use p_target as average of pos and neg priors (unseen mix)
    let mut oracle_fscale = [negative_fscale as f64; 255];
    let min_ret = base / negative_fscale as f64;
    let max_ret = 1.0f64; // retention at base_fscale

    // Start with everything at max fscale (min retention)
    // Budget: how much more retention can we add?
    let baseline_target_ret: f64 = (0..255)
        .map(|i| {
            let tf = (test_pos_hist[i] as f64 / pt + test_neg_hist[i] as f64 / nt) / 2.0;
            tf * min_ret
        })
        .sum();

    let mut remaining_budget = target_ret - baseline_target_ret;

    for &(i, _pf, _nf) in &idx_bins {
        if remaining_budget <= 0.0 {
            break;
        }
        let tf = (test_pos_hist[i] as f64 / pt + test_neg_hist[i] as f64 / nt) / 2.0;
        let extra_ret = tf * (max_ret - min_ret);
        if extra_ret <= 0.0 {
            oracle_fscale[i] = base_fscale as f64;
            continue;
        }
        if extra_ret <= remaining_budget {
            oracle_fscale[i] = base_fscale as f64;
            remaining_budget -= extra_ret;
        } else {
            // Partial: find fscale that uses exactly the remaining budget
            let needed_ret = remaining_budget / tf + min_ret;
            let fs = base / needed_ret;
            oracle_fscale[i] = fs.clamp(base, negative_fscale as f64);
            remaining_budget = 0.0;
        }
    }

    // Compute oracle fold enrichment on test data
    let oracle_pos_ret: f64 = (0..255)
        .map(|i| {
            let pf = test_pos_hist[i] as f64 / pt;
            pf * (base / oracle_fscale[i]).min(1.0)
        })
        .sum();
    let oracle_neg_ret: f64 = (0..255)
        .map(|i| {
            let nf = test_neg_hist[i] as f64 / nt;
            nf * (base / oracle_fscale[i]).min(1.0)
        })
        .sum();
    let oracle_fold = if oracle_neg_ret > 1e-30 { oracle_pos_ret / oracle_neg_ret } else { f64::INFINITY };
    let oracle_pos_eff = if oracle_pos_ret > 0.0 { base / oracle_pos_ret } else { f64::INFINITY };
    let oracle_neg_eff = if oracle_neg_ret > 0.0 { base / oracle_neg_ret } else { f64::INFINITY };

    // Compute oracle fscale per display range (weighted average)
    let range_oracle: Vec<f64> = ranges.iter().map(|&(lo, hi, _)| {
        let lo_idx = (lo as i16 + 127) as usize;
        let hi_idx = (hi as i16 + 127) as usize;
        let total_weight: f64 = (lo_idx..=hi_idx)
            .map(|i| test_pos_hist[i] as f64 + test_neg_hist[i] as f64)
            .sum();
        if total_weight > 0.0 {
            let weighted_fs: f64 = (lo_idx..=hi_idx)
                .map(|i| (test_pos_hist[i] as f64 + test_neg_hist[i] as f64) * oracle_fscale[i])
                .sum();
            weighted_fs / total_weight
        } else {
            negative_fscale as f64
        }
    }).collect();

    // Compute current actual fold using per-index fscale (not range midpoints)
    let cur_pos_ret: f64 = (0..255).map(|i| {
        let pf = test_pos_hist[i] as f64 / pt;
        let w = (i as i16 - 127) as i8;
        let fs = table.effective_fscale_at(w);
        pf * (base / fs).min(1.0)
    }).sum();
    let cur_neg_ret: f64 = (0..255).map(|i| {
        let nf = test_neg_hist[i] as f64 / nt;
        let w = (i as i16 - 127) as i8;
        let fs = table.effective_fscale_at(w);
        nf * (base / fs).min(1.0)
    }).sum();
    let cur_fold = if cur_neg_ret > 1e-30 { cur_pos_ret / cur_neg_ret } else { f64::INFINITY };

    // Compute weighted-average current fscale per range (using per-index values)
    let range_cur: Vec<f64> = ranges.iter().map(|&(lo, hi, _)| {
        let lo_idx = (lo as i16 + 127) as usize;
        let hi_idx = (hi as i16 + 127) as usize;
        let total_weight: f64 = (lo_idx..=hi_idx)
            .map(|i| test_pos_hist[i] as f64 + test_neg_hist[i] as f64)
            .sum();
        if total_weight > 0.0 {
            let weighted_fs: f64 = (lo_idx..=hi_idx)
                .map(|i| {
                    let w = (i as i16 - 127) as i8;
                    let fs = table.effective_fscale_at(w);
                    (test_pos_hist[i] as f64 + test_neg_hist[i] as f64) * fs
                })
                .sum();
            weighted_fs / total_weight
        } else {
            negative_fscale as f64
        }
    }).collect();

    for (j, b) in bins.iter().enumerate() {
        println!(
            "    {:<11} {:>10} {:>10} {:>7.2}% {:>7.2}% {:>7.1}x {:>10.0} {:>10.0}",
            b.label, b.pos_count, b.neg_count,
            b.pos_frac * 100.0, b.neg_frac * 100.0,
            b.ratio, range_cur[j], range_oracle[j]
        );
    }

    let cur_pos_eff = if cur_pos_ret > 0.0 { base / cur_pos_ret } else { f64::INFINITY };
    let cur_neg_eff = if cur_neg_ret > 0.0 { base / cur_neg_ret } else { f64::INFINITY };

    println!();
    println!("    Current:  pos_eff={:.0}, neg_eff={:.0}, fold={:.2}x",
        cur_pos_eff, cur_neg_eff, cur_fold);
    println!("    Oracle:   pos_eff={:.0}, neg_eff={:.0}, fold={:.2}x  (theoretical max at target={})",
        oracle_pos_eff, oracle_neg_eff, oracle_fold, target_fscale);
    println!("    Gap:      {:.1}x fold left on the table", oracle_fold / cur_fold);
    println!();
}

fn measure_per_record_retention(
    records: &[FastaRecord],
    table: &HashBiasTable,
    k: u8,
) -> (f64, f64, f64, f64) {
    let frac_max = u64::MAX / table.fscale();
    let retentions: Vec<f64> = records
        .par_iter()
        .filter_map(|rec| {
            if rec.seq.len() < k as usize {
                return None;
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
                Some(passed as f64 / sketch as f64 * 100.0)
            } else {
                None
            }
        })
        .collect();

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

struct SoftTrialResult {
    target_fscale: u64,
    negative_fscale: u64,
    unseen_fscale: u64,
    train_pos_ret: f64,
    test_pos_ret: f64,
    train_neg_ret: f64,
    test_neg_ret: f64,
    fold_train: f64,
    fold_test: f64,
    eff_fs_test_pos: f64,
    eff_fs_test_neg: f64,
    pos_gap: f64,
    neg_gap: f64,
    build_secs: f64,
}

fn run_soft_trial(
    target_fscale: u64,
    negative_fscale: u64,
    unseen_fscale: Option<u64>,
    train_plasmids: &[FastaRecord],
    test_plasmids: &[FastaRecord],
    train_chromosomes: &[FastaRecord],
    test_chromosomes: &[FastaRecord],
) -> Result<SoftTrialResult> {
    let config = BiasCreateConfig {
        cms: CMSConfig {
            k: K,
            fscale: FSCALE,
            width: WIDTH,
            depth: DEPTH,
        },
        alpha: ALPHA,
        target_fscale: Some(target_fscale),
        negative_fscale: Some(negative_fscale),
        unseen_fscale,
    };

    let start = Instant::now();

    let table = HashBiasTable::create(
        &[Path::new(PLASMID_TRAIN)],
        &[Path::new(CHROMO_TRAIN)],
        &config,
        None,
    )?;
    let build_secs = start.elapsed().as_secs_f64();

    let tp = measure_retention(train_plasmids, &table, K);
    let tsp = measure_retention(test_plasmids, &table, K);
    let tn = measure_retention(train_chromosomes, &table, K);
    let tsn = measure_retention(test_chromosomes, &table, K);

    let fold_train = if tn.retention_pct() > 0.0 {
        tp.retention_pct() / tn.retention_pct()
    } else {
        f64::INFINITY
    };
    let fold_test = if tsn.retention_pct() > 0.0 {
        tsp.retention_pct() / tsn.retention_pct()
    } else {
        f64::INFINITY
    };

    Ok(SoftTrialResult {
        target_fscale,
        negative_fscale,
        unseen_fscale: unseen_fscale.unwrap_or(0),
        train_pos_ret: tp.retention_pct(),
        test_pos_ret: tsp.retention_pct(),
        train_neg_ret: tn.retention_pct(),
        test_neg_ret: tsn.retention_pct(),
        fold_train,
        fold_test,
        eff_fs_test_pos: tsp.effective_fscale(FSCALE),
        eff_fs_test_neg: tsn.effective_fscale(FSCALE),
        pos_gap: tp.retention_pct() - tsp.retention_pct(),
        neg_gap: tn.retention_pct() - tsn.retention_pct(),
        build_secs,
    })
}

fn print_summary_header() {
    println!(
        "{:<10} {:<10} {:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>8} {:>8} {:>6}",
        "tgt_fs", "neg_fs", "unseen_fs",
        "trn_pos", "tst_pos", "trn_neg", "tst_neg",
        "fold_trn", "fold_tst",
        "eff_pos", "eff_neg",
        "pos_gap", "neg_gap", "time"
    );
    println!("{}", "-".repeat(148));
}

fn print_summary_row(r: &SoftTrialResult) {
    println!(
        "{:<10} {:<10} {:<10} {:>9.2}% {:>9.2}% {:>9.2}% {:>9.2}% {:>9.2}x {:>9.2}x {:>10.0} {:>10.0} {:>7.1}pp {:>7.1}pp {:>5.1}s",
        r.target_fscale, r.negative_fscale,
        if r.unseen_fscale > 0 { format!("{}", r.unseen_fscale) } else { "-".to_string() },
        r.train_pos_ret, r.test_pos_ret,
        r.train_neg_ret, r.test_neg_ret,
        r.fold_train, r.fold_test,
        r.eff_fs_test_pos, r.eff_fs_test_neg,
        r.pos_gap, r.neg_gap,
        r.build_secs
    );
}

fn run_detailed(
    label: &str,
    target_fscale: u64,
    negative_fscale: u64,
    unseen_fscale: Option<u64>,
    train_plasmids: &[FastaRecord],
    test_plasmids: &[FastaRecord],
    train_chromosomes: &[FastaRecord],
    test_chromosomes: &[FastaRecord],
) -> Result<()> {
    let config = BiasCreateConfig {
        cms: CMSConfig {
            k: K,
            fscale: FSCALE,
            width: WIDTH,
            depth: DEPTH,
        },
        alpha: ALPHA,
        target_fscale: Some(target_fscale),
        negative_fscale: Some(negative_fscale),
        unseen_fscale,
    };

    println!("--- {} (tgt={}, neg={}, unseen={}) ---",
        label, target_fscale, negative_fscale,
        unseen_fscale.map_or("-".to_string(), |v| v.to_string()));

    let start = Instant::now();
    let table = HashBiasTable::create(
        &[Path::new(PLASMID_TRAIN)],
        &[Path::new(CHROMO_TRAIN)],
        &config,
        None,
    )?;
    let build_time = start.elapsed().as_secs_f64();
    println!("  Built in {build_time:.1}s");

    println!(
        "  {:.<35} {:>12} {:>12} {:>10} {:>12}",
        "Population", "sketch_kmers", "retained", "sketch_%", "eff_fscale"
    );
    println!("  {}", "-".repeat(87));

    for (lbl, records) in [
        ("Train plasmids (seen)", train_plasmids),
        ("Test plasmids (unseen)", test_plasmids),
        ("Train chromosomes (seen neg)", train_chromosomes),
        ("Test chromosomes (unseen neg)", test_chromosomes),
    ] {
        let s = measure_retention(records, &table, K);
        println!(
            "  {:.<35} {:>12} {:>12} {:>9.2}% {:>12.0}",
            lbl, s.sketch_kmers, s.sketch_retained, s.retention_pct(),
            s.effective_fscale(FSCALE)
        );
    }

    println!();
    println!(
        "  {:.<35} {:>10} {:>10} {:>10} {:>10}",
        "Population", "mean_%", "std_%", "min_%", "max_%"
    );
    println!("  {}", "-".repeat(59));

    for (lbl, records) in [
        ("Train plasmids (seen)", train_plasmids),
        ("Test plasmids (unseen)", test_plasmids),
        ("Train chromosomes (seen neg)", train_chromosomes),
        ("Test chromosomes (unseen neg)", test_chromosomes),
    ] {
        let (mean, std, min, max) = measure_per_record_retention(records, &table, K);
        println!(
            "  {:.<35} {:>9.2}% {:>9.2}% {:>9.2}% {:>9.2}%",
            lbl, mean, std, min, max
        );
    }
    println!();

    println!("  === Weight Distributions ===");
    let mut test_pos_hist = [0u64; 255];
    let mut test_neg_hist = [0u64; 255];
    for (lbl, records, capture) in [
        ("Train plasmids (seen)", train_plasmids, None),
        ("Test plasmids (unseen)", test_plasmids, Some(&mut test_pos_hist as &mut [u64; 255])),
        ("Train chromosomes (seen neg)", train_chromosomes, None),
        ("Test chromosomes (unseen neg)", test_chromosomes, Some(&mut test_neg_hist as &mut [u64; 255])),
    ] {
        let hist = measure_weight_distribution(records, &table, K);
        print_weight_distribution(lbl, &hist, &table);
        if let Some(dst) = capture {
            *dst = hist;
        }
    }

    print_oracle_analysis(
        &test_pos_hist, &test_neg_hist, &table,
        FSCALE, target_fscale, negative_fscale,
    );

    Ok(())
}

fn main() -> Result<()> {
    for p in [PLASMID_TRAIN, PLASMID_TEST, CHROMO_TRAIN, CHROMO_TEST] {
        if !Path::new(p).exists() {
            anyhow::bail!("Missing test file: {p}");
        }
    }

    println!("Bias Table Effectiveness — Soft Filter Sweep");
    println!("==============================================");
    println!("  Train positive: {PLASMID_TRAIN}");
    println!("  Test positive:  {PLASMID_TEST}");
    println!("  Train negative: {CHROMO_TRAIN}");
    println!("  Test negative:  {CHROMO_TEST}");
    println!("  k={K}, fscale={FSCALE}, alpha={ALPHA}, width={}M, depth={DEPTH}", WIDTH >> 20);
    println!();

    println!("Loading FASTA records...");
    let train_plasmids = read_fasta_records(Path::new(PLASMID_TRAIN))?;
    let test_plasmids = read_fasta_records(Path::new(PLASMID_TEST))?;
    let train_chromosomes = read_fasta_records(Path::new(CHROMO_TRAIN))?;
    let test_chromosomes = read_fasta_records(Path::new(CHROMO_TEST))?;
    println!(
        "  Train: {} plasmids, {} chromosomes",
        train_plasmids.len(), train_chromosomes.len()
    );
    println!(
        "  Test:  {} plasmids, {} chromosomes",
        test_plasmids.len(), test_chromosomes.len()
    );
    println!();

    println!("=== Sweep: target_fscale x negative_fscale (unseen = target) ===");
    print_summary_header();

    let target_fscales: Vec<u64> = vec![500, 1000, 2000, 5000];
    let neg_fscales: Vec<u64> = vec![10000, 50000, 100000, 500000];

    let mut configs: Vec<(u64, u64)> = Vec::new();
    for &tf in &target_fscales {
        for &nf in &neg_fscales {
            if nf > tf {
                configs.push((tf, nf));
            }
        }
    }

    let results: Vec<_> = configs
        .par_iter()
        .map(|&(tf, nf)| {
            (tf, nf, run_soft_trial(tf, nf, Some(tf),
                &train_plasmids, &test_plasmids, &train_chromosomes, &test_chromosomes))
        })
        .collect();

    for (tf, nf, res) in &results {
        match res {
            Ok(r) => print_summary_row(r),
            Err(e) => println!("tgt={tf} neg={nf}: ERROR: {e}"),
        }
    }

    println!();
    println!("=== Detailed Results with Weight Distributions ===");
    println!();

    run_detailed("User config (100/1000/1000/100000)", 1000, 100000, Some(1000),
        &train_plasmids, &test_plasmids, &train_chromosomes, &test_chromosomes)?;

    Ok(())
}
