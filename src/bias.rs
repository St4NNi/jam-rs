use anyhow::{Context, Result};
use indicatif::ProgressBar;
use jamhash::jamhash_u64;
use needletail::{Sequence, parse_fastx_file};
use std::io::{Read, Write};
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

const BIAS_MAGIC: &[u8; 4] = b"BIA1";
const BIAS_VERSION: u32 = 1;

const DEFAULT_CMS_WIDTH: usize = 1 << 20;
const DEFAULT_CMS_DEPTH: usize = 5;
const QUANTIZATION_SCALE: f32 = 10.0;
const MAX_SAMPLE_HASHES: usize = 100_000;
const HISTOGRAM_EPS: f64 = 1e-6;

#[derive(Debug, Clone)]
pub struct CMSConfig {
    pub width: usize,
    pub depth: usize,
    pub k: u8,
    pub fscale: u64,
}

impl Default for CMSConfig {
    fn default() -> Self {
        Self {
            width: DEFAULT_CMS_WIDTH,
            depth: DEFAULT_CMS_DEPTH,
            k: 21,
            fscale: 1000,
        }
    }
}

#[derive(Debug, Clone)]
pub struct CountMinSketch {
    width: usize,
    depth: usize,
    seeds: Vec<u64>,
    counts: Vec<u64>,
}

impl CountMinSketch {
    pub fn new(width: usize, depth: usize) -> Self {
        let seeds: Vec<u64> = (0..depth).map(|i| jamhash_u64(i as u64)).collect();
        let counts = vec![0u64; width * depth];
        Self {
            width,
            depth,
            seeds,
            counts,
        }
    }

    pub fn with_seeds(width: usize, depth: usize, seeds: Vec<u64>) -> Self {
        assert_eq!(seeds.len(), depth);
        let counts = vec![0u64; width * depth];
        Self {
            width,
            depth,
            seeds,
            counts,
        }
    }

    #[inline]
    fn index(&self, row: usize, hash: u64) -> usize {
        let mixed = hash.wrapping_mul(self.seeds[row]);
        row * self.width + (mixed as usize % self.width)
    }

    #[inline]
    pub fn increment(&mut self, hash: u64) {
        for row in 0..self.depth {
            let idx = self.index(row, hash);
            self.counts[idx] = self.counts[idx].saturating_add(1);
        }
    }

    #[inline]
    pub fn estimate(&self, hash: u64) -> u64 {
        (0..self.depth)
            .map(|row| self.counts[self.index(row, hash)])
            .min()
            .unwrap_or(0)
    }

    pub fn width(&self) -> usize {
        self.width
    }
    pub fn depth(&self) -> usize {
        self.depth
    }
    pub fn seeds(&self) -> &[u64] {
        &self.seeds
    }
    pub fn counts(&self) -> &[u64] {
        &self.counts
    }

    pub fn cell_stats(&self) -> (u64, u64, f64, f64, usize) {
        let min = *self.counts.iter().min().unwrap_or(&0);
        let max = *self.counts.iter().max().unwrap_or(&0);
        let sum: u64 = self.counts.iter().sum();
        let mean = sum as f64 / self.counts.len() as f64;
        let variance: f64 = self
            .counts
            .iter()
            .map(|&c| {
                let d = c as f64 - mean;
                d * d
            })
            .sum::<f64>()
            / self.counts.len() as f64;
        let non_zero = self.counts.iter().filter(|&&c| c > 0).count();
        (min, max, mean, variance.sqrt(), non_zero)
    }
}

#[derive(Debug, Clone)]
pub(crate) struct RawHashCounts {
    pub(crate) config: CMSConfig,
    pub(crate) cms: CountMinSketch,
    pub(crate) total: u64,
    pub(crate) samples: Vec<u64>,
}

impl RawHashCounts {
    pub(crate) fn new(config: CMSConfig) -> Self {
        let cms = CountMinSketch::new(config.width, config.depth);
        Self {
            config,
            cms,
            total: 0,
            samples: Vec::with_capacity(MAX_SAMPLE_HASHES),
        }
    }

    pub(crate) fn build(
        paths: &[&Path],
        config: CMSConfig,
        record_counter: &AtomicU64,
        hash_counter: &AtomicU64,
    ) -> Result<Self> {
        let frac_max = u64::MAX / config.fscale;
        let k = config.k;

        let mut raw = RawHashCounts::new(config);
        for path in paths {
            process_path(&mut raw, path, k, frac_max, record_counter, hash_counter)?;
        }

        if raw.samples.len() > MAX_SAMPLE_HASHES {
            downsample_samples(&mut raw.samples);
        }

        Ok(raw)
    }
}

fn process_path(
    raw: &mut RawHashCounts,
    path: &Path,
    k: u8,
    frac_max: u64,
    record_counter: &AtomicU64,
    hash_counter: &AtomicU64,
) -> Result<()> {
    let mut reader = match parse_fastx_file(path) {
        Ok(reader) => reader,
        Err(e) if e.kind == needletail::errors::ParseErrorKind::EmptyFile => {
            return Ok(());
        }
        Err(e) => {
            return Err(e).with_context(|| format!("Failed to parse: {}", path.display()));
        }
    };

    while let Some(record) = reader.next() {
        let record = record.context("Failed to parse sequence record")?;
        let seq = record.normalize(false);
        record_counter.fetch_add(1, Ordering::Relaxed);

        if seq.len() < k as usize {
            continue;
        }

        for (_, kmer, _) in seq.bit_kmers(k, true) {
            let hash = jamhash_u64(kmer.0);
            if hash < frac_max {
                raw.cms.increment(hash);
                raw.total += 1;
                if raw.samples.len() < MAX_SAMPLE_HASHES {
                    raw.samples.push(hash);
                } else {
                    let seen = raw.total;
                    let pick = (jamhash_u64(hash ^ seen) % seen) as usize;
                    if pick < MAX_SAMPLE_HASHES {
                        raw.samples[pick] = hash;
                    }
                }
                hash_counter.fetch_add(1, Ordering::Relaxed);
            }
        }
    }

    Ok(())
}

fn downsample_samples(samples: &mut Vec<u64>) {
    if samples.len() <= MAX_SAMPLE_HASHES {
        return;
    }
    samples.sort_unstable_by_key(|&hash| jamhash_u64(hash));
    samples.truncate(MAX_SAMPLE_HASHES);
}

#[derive(Debug, Clone)]
pub struct BiasCreateConfig {
    pub cms: CMSConfig,
    pub alpha: f32,
    pub target_fscale: Option<u64>,
    pub negative_fscale: Option<u64>,
    pub unseen_fscale: Option<u64>,
}

#[derive(Debug, Clone, Copy)]
pub struct CalibrationResult {
    pub threshold: i8,
    pub positive_retention: f32,
    pub negative_retention: f32,
    pub fold_enrichment: f32,
    pub max_fold_enrichment: f32,
}

#[derive(Debug, Clone)]
pub struct HashBiasTable {
    pub config: CMSConfig,
    seeds: Vec<u64>,
    weights: Vec<i8>,
    pub alpha: f32,
    pub threshold: i8,
    pub positive_retention: f32,
    pub negative_retention: f32,
    pub max_fold_enrichment: f32,
    pub target_fscale: u64,
    pub negative_fscale: u64,
    pub unseen_fscale: u64,
    pub fscale_lut: [u64; 255],
    frac_max_lut: [u64; 255],
    cms_pos_counts: Option<Vec<u32>>,
    cms_neg_counts: Option<Vec<u32>>,
    cms_pos_total: f64,
    cms_neg_total: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct SoftFilterReferencePoint {
    pub weight: i8,
    pub weight_f32: f32,
    pub effective_fscale: f64,
    pub retention_pct: f64,
    pub vs_base: f64,
}

fn validate_cms_compatibility(positive: &RawHashCounts, negative: &RawHashCounts) -> Result<()> {
    if positive.config.k != negative.config.k {
        anyhow::bail!(
            "k-mer size mismatch: positive={}, negative={}",
            positive.config.k,
            negative.config.k
        );
    }
    if positive.config.fscale != negative.config.fscale {
        anyhow::bail!(
            "fscale mismatch: positive={}, negative={}",
            positive.config.fscale,
            negative.config.fscale
        );
    }
    if positive.config.width != negative.config.width
        || positive.config.depth != negative.config.depth
    {
        anyhow::bail!(
            "CMS dimensions mismatch: positive={}x{}, negative={}x{}",
            positive.config.width,
            positive.config.depth,
            negative.config.width,
            negative.config.depth
        );
    }
    Ok(())
}

fn histogram_weights(hashes: &[u64], bias_table: &HashBiasTable) -> ([f64; 255], [u64; 255]) {
    let mut counts = [0u64; 255];
    for &h in hashes {
        let w = bias_table.weight(h);
        let idx = HashBiasTable::weight_to_index(w);
        counts[idx] += 1;
    }
    let total: f64 = counts.iter().sum::<u64>() as f64 + 255.0 * HISTOGRAM_EPS;
    let mut hist = [0.0f64; 255];
    for (h, &c) in hist.iter_mut().zip(counts.iter()) {
        *h = (c as f64 + HISTOGRAM_EPS) / total;
    }
    (hist, counts)
}

/// Fill unreliable (sparse) bins via log-space linear interpolation between reliable neighbors.
#[allow(clippy::needless_range_loop)]
fn empirical_prior_from_weights(weights: &[i8]) -> [f64; 255] {
    let mut counts = [0u64; 255];
    for &w in weights {
        counts[HashBiasTable::weight_to_index(w)] += 1;
    }

    let total = weights.len() as f64 + 255.0 * HISTOGRAM_EPS;
    let mut prior = [0.0f64; 255];
    for i in 0..255 {
        prior[i] = (counts[i] as f64 + HISTOGRAM_EPS) / total;
    }
    prior
}

const DROP_MODE_LOG_CAP: f64 = 1e18;


#[inline]
fn retention_from_log_fscale(log_fs: f64, base_fs: f64) -> f64 {
    (base_fs / log_fs.exp()).min(1.0)
}

fn compute_retentions(
    x: &[f64; 255],
    p_pos: &[f64; 255],
    p_neg: &[f64; 255],
    p_target: &[f64; 255],
    base_fs: f64,
) -> (f64, f64, f64) {
    let mut pos_ret = 0.0f64;
    let mut neg_ret = 0.0f64;
    let mut target_ret = 0.0f64;
    for i in 0..255 {
        let r = retention_from_log_fscale(x[i], base_fs);
        pos_ret += p_pos[i] * r;
        neg_ret += p_neg[i] * r;
        target_ret += p_target[i] * r;
    }
    (pos_ret, neg_ret, target_ret)
}

#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
fn target_fscale_lut(
    p_pos: &[f64; 255],
    p_neg: &[f64; 255],
    _p_target: &[f64; 255],
    _pos_counts: &[u64; 255],
    _neg_counts: &[u64; 255],
    base_fscale: u64,
    target_fscale: u64,
    negative_fscale: u64,
    unseen_fscale: u64,
) -> Result<[u64; 255]> {
    if base_fscale == 0 || target_fscale == 0 || negative_fscale == 0 {
        anyhow::bail!("fscale values must be > 0");
    }
    if target_fscale < base_fscale {
        anyhow::bail!(
            "target_fscale ({}) must be >= base fscale ({})",
            target_fscale,
            base_fscale
        );
    }
    if negative_fscale != u64::MAX && negative_fscale < target_fscale {
        anyhow::bail!(
            "negative_fscale ({}) must be >= target_fscale ({})",
            negative_fscale,
            target_fscale
        );
    }

    let base_fs = base_fscale as f64;
    let neg_fs = if negative_fscale == u64::MAX {
        DROP_MODE_LOG_CAP
    } else {
        negative_fscale as f64
    };
    let f0 = if unseen_fscale > 0 {
        unseen_fscale as f64
    } else {
        target_fscale as f64
    };
    let w0_idx = HashBiasTable::weight_to_index(0);

    let r_max = 1.0f64;
    let r_min = base_fs / neg_fs;

    let budget_fscale = (target_fscale as f64 * 1.5).min(neg_fs);
    let target_ret = base_fs / budget_fscale;

    let fixed_neg_ret: f64 = [0usize, w0_idx, 254]
        .iter()
        .map(|&i| {
            let fs = match i {
                0 => neg_fs,
                127 => f0,
                254 => base_fs,
                _ => unreachable!(),
            };
            p_neg[i] * (base_fs / fs).min(1.0)
        })
        .sum();

    let neg_budget = (target_ret - fixed_neg_ret).max(0.0);

    struct BinEntry {
        idx: usize,
        ratio: f64,
        neg_cost: f64,
    }

    let mut bins: Vec<BinEntry> = Vec::new();
    for i in 0..255 {
        if i == 0 || i == w0_idx || i == 254 {
            continue;
        }
        let neg_cost = p_neg[i] * (r_max - r_min);
        let ratio = if p_neg[i] > 1e-30 {
            p_pos[i] / p_neg[i]
        } else if p_pos[i] > 1e-30 {
            f64::INFINITY
        } else {
            0.0
        };
        bins.push(BinEntry { idx: i, ratio, neg_cost });
    }

    bins.sort_by(|a, b| b.ratio.partial_cmp(&a.ratio).unwrap_or(std::cmp::Ordering::Equal));

    let all_min_neg_ret: f64 = bins.iter().map(|b| p_neg[b.idx] * r_min).sum();
    let mut remaining_budget = neg_budget - all_min_neg_ret;

    let mut retention = [0.0f64; 255];
    retention[0] = (base_fs / neg_fs).min(1.0);
    retention[w0_idx] = (base_fs / f0).min(1.0);
    retention[254] = 1.0;

    for b in &bins {
        retention[b.idx] = r_min;
    }

    for b in &bins {
        if remaining_budget <= 1e-30 {
            break;
        }
        if b.neg_cost <= remaining_budget {
            retention[b.idx] = r_max;
            remaining_budget -= b.neg_cost;
        } else {
            let frac = remaining_budget / b.neg_cost;
            retention[b.idx] = r_min + frac * (r_max - r_min);
            remaining_budget = 0.0;
        }
    }

    let mut lut = [0u64; 255];
    for i in 0..255 {
        let r = retention[i].clamp(r_min, r_max);
        let fs = base_fs / r;
        lut[i] = (fs.round() as u64).clamp(base_fscale, if negative_fscale == u64::MAX { u64::MAX } else { negative_fscale });
    }

    lut[0] = if negative_fscale == u64::MAX {
        u64::MAX
    } else {
        negative_fscale
    };
    lut[w0_idx] = if unseen_fscale > 0 {
        unseen_fscale
    } else {
        target_fscale
    };
    lut[254] = base_fscale;

    let (final_pos_ret, final_neg_ret, _) = {
        let mut x = [0.0f64; 255];
        for i in 0..255 {
            let fs = if lut[i] == u64::MAX { neg_fs } else { lut[i] as f64 };
            x[i] = fs.ln();
        }
        compute_retentions(&x, p_pos, p_neg, _p_target, base_fs)
    };
    let final_fold = if final_neg_ret > 1e-30 {
        final_pos_ret / final_neg_ret
    } else {
        f64::INFINITY
    };
    let pos_eff = base_fs / final_pos_ret;
    let neg_eff = base_fs / final_neg_ret;
    eprintln!(
        "  Knapsack LUT: pos_eff={:.0}, neg_eff={:.0}, fold={:.2}x, neg_ret={:.4} (budget={:.4})",
        pos_eff, neg_eff, final_fold, final_neg_ret, target_ret
    );

    Ok(lut)
}

impl HashBiasTable {
    pub fn create(
        positive_paths: &[&Path],
        negative_paths: &[&Path],
        config: &BiasCreateConfig,
        progress: Option<ProgressBar>,
    ) -> Result<Self> {
        // Validate fscale parameters early, before the expensive sketching step.
        let target_fscale = config.target_fscale.unwrap_or(0);
        let negative_fscale = config.negative_fscale.unwrap_or(0);
        let unseen_fscale = config.unseen_fscale.unwrap_or(0);
        let global_fscale = config.cms.fscale;

        if target_fscale > 0 || negative_fscale > 0 {
            if target_fscale == 0 || negative_fscale == 0 {
                anyhow::bail!("target_fscale and negative_fscale must be set together");
            }
            if target_fscale < global_fscale {
                anyhow::bail!(
                    "target_fscale ({}) must be >= global fscale ({})",
                    target_fscale,
                    global_fscale
                );
            }
            if negative_fscale <= target_fscale {
                anyhow::bail!(
                    "negative_fscale ({}) must be > target_fscale ({})",
                    negative_fscale,
                    target_fscale
                );
            }
            if unseen_fscale > 0 {
                if unseen_fscale < target_fscale {
                    anyhow::bail!(
                        "unseen_fscale ({}) must be >= target_fscale ({})",
                        unseen_fscale,
                        target_fscale
                    );
                }
                if negative_fscale != u64::MAX && unseen_fscale > negative_fscale {
                    anyhow::bail!(
                        "unseen_fscale ({}) must be <= negative_fscale ({})",
                        unseen_fscale,
                        negative_fscale
                    );
                }
            }
        } else if unseen_fscale > 0 {
            anyhow::bail!(
                "--unseen-fscale requires soft filter mode (set --target-fscale and --negative-fscale)"
            );
        }

        let record_counter = Arc::new(AtomicU64::new(0));
        let hash_counter = Arc::new(AtomicU64::new(0));
        let stop_flag = Arc::new(AtomicBool::new(false));

        let update_handle = progress.as_ref().map(|pb| {
            let pb = pb.clone();
            let record_counter = Arc::clone(&record_counter);
            let hash_counter = Arc::clone(&hash_counter);
            let stop_flag = Arc::clone(&stop_flag);

            std::thread::spawn(move || {
                loop {
                    if stop_flag.load(Ordering::Relaxed) || pb.is_finished() {
                        break;
                    }
                    let records = record_counter.load(Ordering::Relaxed);
                    let hashes = hash_counter.load(Ordering::Relaxed);
                    pb.set_message(format!(
                        "{} records, {} hashes",
                        format_number(records),
                        format_number(hashes)
                    ));
                    std::thread::sleep(std::time::Duration::from_millis(100));
                }
            })
        });

        let (pos_raw, neg_raw) = rayon::join(
            || {
                RawHashCounts::build(
                    positive_paths,
                    config.cms.clone(),
                    &record_counter,
                    &hash_counter,
                )
            },
            || {
                RawHashCounts::build(
                    negative_paths,
                    config.cms.clone(),
                    &record_counter,
                    &hash_counter,
                )
            },
        );

        stop_flag.store(true, Ordering::Relaxed);
        if let Some(handle) = update_handle {
            let _ = handle.join();
        }

        let pos_raw = pos_raw?;
        let neg_raw = neg_raw?;

        if let Some(ref pb) = progress {
            pb.set_message("Computing bias weights...");
        }

        let table = Self::build(
            &pos_raw,
            &neg_raw,
            config.alpha,
            config.target_fscale.unwrap_or(0),
            config.negative_fscale.unwrap_or(0),
            config.unseen_fscale.unwrap_or(0),
        )?;

        if let Some(ref pb) = progress {
            pb.finish();
        }

        Ok(table)
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn build(
        positive: &RawHashCounts,
        negative: &RawHashCounts,
        alpha: f32,
        target_fscale: u64,
        negative_fscale: u64,
        unseen_fscale: u64,
    ) -> Result<Self> {
        validate_cms_compatibility(positive, negative)?;

        if !alpha.is_finite() || alpha <= 0.0 {
            anyhow::bail!("alpha must be finite and > 0, got {}", alpha);
        }

        if target_fscale > 0 || negative_fscale > 0 {
            if target_fscale == 0 || negative_fscale == 0 {
                anyhow::bail!("target_fscale and negative_fscale must be set together");
            }
            if target_fscale < positive.config.fscale {
                anyhow::bail!(
                    "target_fscale ({}) must be >= global fscale ({})",
                    target_fscale,
                    positive.config.fscale
                );
            }
            if negative_fscale <= target_fscale {
                anyhow::bail!(
                    "negative_fscale ({}) must be > target_fscale ({})",
                    negative_fscale,
                    target_fscale
                );
            }
            if unseen_fscale > 0 {
                if unseen_fscale < target_fscale {
                    anyhow::bail!(
                        "unseen_fscale ({}) must be >= target_fscale ({})",
                        unseen_fscale,
                        target_fscale
                    );
                }
                if negative_fscale != u64::MAX && unseen_fscale > negative_fscale {
                    anyhow::bail!(
                        "unseen_fscale ({}) must be <= negative_fscale ({})",
                        unseen_fscale,
                        negative_fscale
                    );
                }
            }
        } else if unseen_fscale > 0 {
            anyhow::bail!(
                "--unseen-fscale requires soft filter mode (set --target-fscale and --negative-fscale)"
            );
        }

        let width = positive.config.width;
        let depth = positive.config.depth;
        let seeds = positive.cms.seeds().to_vec();

        let pos_counts_raw = positive.cms.counts();
        let neg_counts_raw = negative.cms.counts();

        let pos_total = positive.total as f64;
        let neg_total = negative.total as f64;

        let cms_pos_counts: Vec<u32> = pos_counts_raw
            .iter()
            .map(|&c| c.min(u32::MAX as u64) as u32)
            .collect();
        let cms_neg_counts: Vec<u32> = neg_counts_raw
            .iter()
            .map(|&c| c.min(u32::MAX as u64) as u32)
            .collect();

        let mut weights = vec![0i8; width * depth];

        if pos_total > 0.0 && neg_total > 0.0 {
            let scale = pos_total.max(neg_total);

            for i in 0..(width * depth) {
                let norm_pos = (pos_counts_raw[i] as f64 / pos_total) * scale;
                let norm_neg = (neg_counts_raw[i] as f64 / neg_total) * scale;
                let adj_pos = (norm_pos - norm_neg).max(0.0) as f32;
                let adj_neg = (norm_neg - norm_pos).max(0.0) as f32;

                let log_ratio = ((adj_pos + alpha) / (adj_neg + alpha)).ln();
                let quantized = (log_ratio * QUANTIZATION_SCALE).clamp(-127.0, 127.0) as i8;
                weights[i] = quantized;
            }
        }

        let calibration = calibrate_threshold(positive, negative, &weights, &seeds, width, None)?;

        let base_fscale = positive.config.fscale;
        let fscale_lut = [base_fscale; 255];
        let frac_max_lut = [0u64; 255];

        let mut table = Self {
            config: positive.config.clone(),
            seeds,
            weights,
            alpha,
            threshold: calibration.threshold,
            positive_retention: calibration.positive_retention,
            negative_retention: calibration.negative_retention,
            max_fold_enrichment: calibration.max_fold_enrichment,
            target_fscale,
            negative_fscale,
            unseen_fscale,
            fscale_lut,
            frac_max_lut,
            cms_pos_counts: Some(cms_pos_counts),
            cms_neg_counts: Some(cms_neg_counts),
            cms_pos_total: pos_total,
            cms_neg_total: neg_total,
        };

        if table.is_soft_filter() {
            let fscale = table.config.fscale;
            let frac_max = u64::MAX / fscale;

            let sample_hashes = |raw: &RawHashCounts, max_samples: usize| -> Vec<u64> {
                if raw.samples.len() <= max_samples {
                    return raw.samples.clone();
                }
                let step = raw.samples.len() / max_samples;
                raw.samples
                    .iter()
                    .step_by(step)
                    .take(max_samples)
                    .copied()
                    .collect()
            };

            let pos_hashes = sample_hashes(positive, MAX_SAMPLE_HASHES);
            let neg_hashes = sample_hashes(negative, MAX_SAMPLE_HASHES);

            let pos_in_range: Vec<u64> = pos_hashes
                .iter()
                .copied()
                .filter(|&h| h < frac_max)
                .collect();
            let neg_in_range: Vec<u64> = neg_hashes
                .iter()
                .copied()
                .filter(|&h| h < frac_max)
                .collect();

            let (p_pos, pos_bin_counts) = histogram_weights(&pos_in_range, &table);
            let (p_neg, neg_bin_counts) = histogram_weights(&neg_in_range, &table);
            let p_target = if table.cms_pos_counts.is_some() {
                let mut combined: Vec<u64> = Vec::with_capacity(pos_in_range.len() + neg_in_range.len());
                combined.extend_from_slice(&pos_in_range);
                combined.extend_from_slice(&neg_in_range);
                histogram_weights(&combined, &table).0
            } else {
                empirical_prior_from_weights(&table.weights)
            };

            let lut = target_fscale_lut(
                &p_pos,
                &p_neg,
                &p_target,
                &pos_bin_counts,
                &neg_bin_counts,
                table.config.fscale,
                target_fscale,
                negative_fscale,
                unseen_fscale,
            )?;

            table.fscale_lut = lut;

            table.recompute_frac_max_lut();

            let (pos_ret, neg_ret) = Self::compute_retention(&table, &pos_in_range, &neg_in_range);
            table.positive_retention = pos_ret;
            table.negative_retention = neg_ret;
            if table.negative_retention > 0.0 {
                table.max_fold_enrichment = table.positive_retention / table.negative_retention;
            } else {
                table.max_fold_enrichment = f32::INFINITY;
            }
        }

        table.recompute_frac_max_lut();
        Ok(table)
    }

    fn compute_retention(table: &Self, pos_in_range: &[u64], neg_in_range: &[u64]) -> (f32, f32) {
        let pos_passing = pos_in_range
            .iter()
            .filter(|&&h| table.passes_soft_filter(h))
            .count();
        let neg_passing = neg_in_range
            .iter()
            .filter(|&&h| table.passes_soft_filter(h))
            .count();

        let pos_ret = if pos_in_range.is_empty() {
            0.0
        } else {
            pos_passing as f32 / pos_in_range.len() as f32
        };
        let neg_ret = if neg_in_range.is_empty() {
            0.0
        } else {
            neg_passing as f32 / neg_in_range.len() as f32
        };
        (pos_ret, neg_ret)
    }

    #[inline]
    fn index(&self, row: usize, hash: u64) -> usize {
        let mixed = hash.wrapping_mul(self.seeds[row]);
        row * self.config.width + (mixed as usize % self.config.width)
    }

    #[inline]
    pub fn weight(&self, hash: u64) -> i8 {
        if let (Some(pos), Some(neg)) = (&self.cms_pos_counts, &self.cms_neg_counts) {
            self.cms_query_weight(hash, pos, neg)
        } else {
            self.cell_min_weight(hash)
        }
    }

    #[inline]
    fn cell_min_weight(&self, hash: u64) -> i8 {
        (0..self.config.depth)
            .map(|row| self.weights[self.index(row, hash)])
            .min()
            .unwrap_or(0)
    }

    #[inline]
    fn cms_query_weight(&self, hash: u64, pos: &[u32], neg: &[u32]) -> i8 {
        let mut min_pos = u32::MAX;
        let mut min_neg = u32::MAX;
        for row in 0..self.config.depth {
            let idx = self.index(row, hash);
            min_pos = min_pos.min(pos[idx]);
            min_neg = min_neg.min(neg[idx]);
        }

        let scale = self.cms_pos_total.max(self.cms_neg_total);
        let norm_pos = (min_pos as f64 / self.cms_pos_total) * scale;
        let norm_neg = (min_neg as f64 / self.cms_neg_total) * scale;
        let adj_pos = (norm_pos - norm_neg).max(0.0) as f32;
        let adj_neg = (norm_neg - norm_pos).max(0.0) as f32;

        let log_ratio = ((adj_pos + self.alpha) / (adj_neg + self.alpha)).ln();
        (log_ratio * QUANTIZATION_SCALE).clamp(-127.0, 127.0) as i8
    }

    #[inline]
    pub fn passes_filter(&self, hash: u64) -> bool {
        if !self.is_soft_filter() {
            return self.weight(hash) >= self.threshold;
        }
        self.passes_soft_filter(hash)
    }

    #[inline]
    pub fn is_soft_filter(&self) -> bool {
        self.target_fscale > 0 && self.negative_fscale > 0
    }

    #[inline]
    fn weight_to_index(w: i8) -> usize {
        (w as i16 + 127) as usize
    }

    fn recompute_frac_max_lut(&mut self) {
        for (frac, &fs) in self.frac_max_lut.iter_mut().zip(self.fscale_lut.iter()) {
            *frac = u64::MAX.checked_div(fs).unwrap_or(u64::MAX);
        }
    }

    #[inline]
    fn effective_fscale(&self, w: i8) -> f64 {
        self.fscale_lut[Self::weight_to_index(w)] as f64
    }

    #[inline]
    fn passes_soft_filter(&self, hash: u64) -> bool {
        let w = self.weight(hash);
        hash < self.frac_max_lut[Self::weight_to_index(w)]
    }

    pub fn k(&self) -> u8 {
        self.config.k
    }
    pub fn fscale(&self) -> u64 {
        self.config.fscale
    }

    pub fn fold_enrichment(&self) -> f32 {
        if self.negative_retention > 0.0 {
            self.positive_retention / self.negative_retention
        } else {
            f32::INFINITY
        }
    }

    pub fn effective_fscale_at(&self, w: i8) -> f64 {
        self.effective_fscale(w)
    }

    pub fn effective_fscale_on_population(&self, retention: f32) -> f64 {
        if retention > 0.0 {
            self.config.fscale as f64 / retention as f64
        } else {
            f64::INFINITY
        }
    }

    pub fn effective_fscale_combined(&self) -> f64 {
        let pos = self.positive_retention as f64;
        let neg = self.negative_retention as f64;
        if pos > 0.0 && neg > 0.0 {
            self.config.fscale as f64 / (pos * neg).sqrt()
        } else {
            f64::INFINITY
        }
    }

    fn effective_fscale_on_bucket_distribution(&self, dist: &[f64; 255]) -> f64 {
        let base = self.config.fscale as f64;
        let mut retention = 0.0f64;
        for (i, &p) in dist.iter().enumerate() {
            let eff = self.fscale_lut[i] as f64;
            let r = (base / eff).min(1.0);
            retention += p * r;
        }

        if retention > 0.0 {
            base / retention
        } else {
            f64::INFINITY
        }
    }

    pub fn effective_fscale_target_prior(&self) -> f64 {
        let prior = empirical_prior_from_weights(&self.weights);
        self.effective_fscale_on_bucket_distribution(&prior)
    }

    pub fn negative_fscale_label(&self) -> String {
        if self.negative_fscale == u64::MAX {
            "drop".to_string()
        } else {
            self.negative_fscale.to_string()
        }
    }

    pub fn save(&self, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)
            .with_context(|| format!("Failed to create bias table file: {}", path.display()))?;

        file.write_all(BIAS_MAGIC)?;
        file.write_all(&BIAS_VERSION.to_le_bytes())?;
        file.write_all(&[self.config.k])?;
        file.write_all(&self.config.fscale.to_le_bytes())?;
        file.write_all(&(self.config.width as u32).to_le_bytes())?;
        file.write_all(&[self.config.depth as u8])?;
        file.write_all(&self.alpha.to_le_bytes())?;
        file.write_all(&[self.threshold as u8])?;
        file.write_all(&self.positive_retention.to_le_bytes())?;
        file.write_all(&self.negative_retention.to_le_bytes())?;
        file.write_all(&self.target_fscale.to_le_bytes())?;
        file.write_all(&self.negative_fscale.to_le_bytes())?;
        file.write_all(&0.0f32.to_le_bytes())?; // legacy min_positive_retention slot
        file.write_all(&[0u8; 8])?; // reserved
        file.write_all(&self.unseen_fscale.to_le_bytes())?;

        for &f in &self.fscale_lut {
            file.write_all(&f.to_le_bytes())?;
        }

        for &seed in &self.seeds {
            file.write_all(&seed.to_le_bytes())?;
        }
        for &w in &self.weights {
            file.write_all(&[w as u8])?;
        }

        Ok(())
    }

    pub fn load(path: &Path) -> Result<Self> {
        let mut file = std::fs::File::open(path)
            .with_context(|| format!("Failed to open bias table file: {}", path.display()))?;

        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)?;

        if &magic != BIAS_MAGIC {
            anyhow::bail!(
                "Invalid bias table file (bad magic, expected BIA1): {}",
                path.display()
            );
        }

        let mut buf4 = [0u8; 4];
        file.read_exact(&mut buf4)?;
        let version = u32::from_le_bytes(buf4);
        if version != BIAS_VERSION {
            anyhow::bail!(
                "Unsupported bias table version {} (expected {})",
                version,
                BIAS_VERSION
            );
        }

        let mut k_buf = [0u8; 1];
        file.read_exact(&mut k_buf)?;
        let k = k_buf[0];

        let mut buf8 = [0u8; 8];
        file.read_exact(&mut buf8)?;
        let fscale = u64::from_le_bytes(buf8);

        file.read_exact(&mut buf4)?;
        let width = u32::from_le_bytes(buf4) as usize;

        let mut depth_buf = [0u8; 1];
        file.read_exact(&mut depth_buf)?;
        let depth = depth_buf[0] as usize;

        file.read_exact(&mut buf4)?;
        let alpha = f32::from_le_bytes(buf4);

        let mut threshold_buf = [0u8; 1];
        file.read_exact(&mut threshold_buf)?;
        let threshold = threshold_buf[0] as i8;

        file.read_exact(&mut buf4)?;
        let positive_retention = f32::from_le_bytes(buf4);

        file.read_exact(&mut buf4)?;
        let negative_retention = f32::from_le_bytes(buf4);

        file.read_exact(&mut buf8)?;
        let target_fscale = u64::from_le_bytes(buf8);

        file.read_exact(&mut buf8)?;
        let negative_fscale = u64::from_le_bytes(buf8);

        file.read_exact(&mut buf4)?;
        let _min_positive_retention = f32::from_le_bytes(buf4); // legacy, discarded

        let mut reserved = [0u8; 8];
        file.read_exact(&mut reserved)?;

        file.read_exact(&mut buf8)?;
        let unseen_fscale = u64::from_le_bytes(buf8);

        let mut fscale_lut = [0u64; 255];
        for entry in &mut fscale_lut {
            file.read_exact(&mut buf8)?;
            *entry = u64::from_le_bytes(buf8);
        }

        let mut seeds = Vec::with_capacity(depth);
        for _ in 0..depth {
            file.read_exact(&mut buf8)?;
            seeds.push(u64::from_le_bytes(buf8));
        }

        let mut weights = vec![0i8; width * depth];
        let mut weight_buf = vec![0u8; width * depth];
        file.read_exact(&mut weight_buf)?;
        for (i, &b) in weight_buf.iter().enumerate() {
            weights[i] = b as i8;
        }

        let config = CMSConfig {
            width,
            depth,
            k,
            fscale,
        };

        let max_fold_enrichment = if negative_retention > 0.0 {
            positive_retention / negative_retention
        } else {
            f32::INFINITY
        };

        let mut table = Self {
            config,
            seeds,
            weights,
            alpha,
            threshold,
            positive_retention,
            negative_retention,
            max_fold_enrichment,
            target_fscale,
            negative_fscale,
            unseen_fscale,
            fscale_lut,
            frac_max_lut: [0u64; 255],
            cms_pos_counts: None,
            cms_neg_counts: None,
            cms_pos_total: 0.0,
            cms_neg_total: 0.0,
        };
        table.recompute_frac_max_lut();
        Ok(table)
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let header_size = 4 + 4 + 1 + 8 + 4 + 1 + 4 + 1 + 4 + 4 + 8 + 8 + 4 + 8 + 8 + 255 * 8;
        let seeds_size = self.config.depth * 8;
        let weights_size = self.config.width * self.config.depth;
        let total_size = header_size + seeds_size + weights_size;

        let mut out = Vec::with_capacity(total_size);
        out.extend_from_slice(BIAS_MAGIC);
        out.extend_from_slice(&BIAS_VERSION.to_le_bytes());
        out.push(self.config.k);
        out.extend_from_slice(&self.config.fscale.to_le_bytes());
        out.extend_from_slice(&(self.config.width as u32).to_le_bytes());
        out.push(self.config.depth as u8);
        out.extend_from_slice(&self.alpha.to_le_bytes());
        out.push(self.threshold as u8);
        out.extend_from_slice(&self.positive_retention.to_le_bytes());
        out.extend_from_slice(&self.negative_retention.to_le_bytes());
        out.extend_from_slice(&self.target_fscale.to_le_bytes());
        out.extend_from_slice(&self.negative_fscale.to_le_bytes());
        out.extend_from_slice(&0.0f32.to_le_bytes()); // legacy min_positive_retention slot
        out.extend_from_slice(&[0u8; 8]); // reserved
        out.extend_from_slice(&self.unseen_fscale.to_le_bytes());

        for &f in &self.fscale_lut {
            out.extend_from_slice(&f.to_le_bytes());
        }

        for &seed in &self.seeds {
            out.extend_from_slice(&seed.to_le_bytes());
        }
        for &w in &self.weights {
            out.push(w as u8);
        }

        out
    }

    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        const HEADER_SIZE: usize = 2111;
        if data.len() < HEADER_SIZE {
            anyhow::bail!(
                "Bias table data too small: {} bytes (expected >= {})",
                data.len(),
                HEADER_SIZE
            );
        }

        let magic: [u8; 4] = data[0..4].try_into().unwrap();
        if &magic != BIAS_MAGIC {
            anyhow::bail!("Invalid bias table magic bytes (expected BIA1)");
        }

        let version = u32::from_le_bytes(data[4..8].try_into().unwrap());
        if version != BIAS_VERSION {
            anyhow::bail!(
                "Unsupported bias table version {} (expected {})",
                version,
                BIAS_VERSION
            );
        }

        let k = data[8];
        let fscale = u64::from_le_bytes(data[9..17].try_into().unwrap());
        let width = u32::from_le_bytes(data[17..21].try_into().unwrap()) as usize;
        let depth = data[21] as usize;
        let alpha = f32::from_le_bytes(data[22..26].try_into().unwrap());
        let threshold = data[26] as i8;
        let positive_retention = f32::from_le_bytes(data[27..31].try_into().unwrap());
        let negative_retention = f32::from_le_bytes(data[31..35].try_into().unwrap());
        let target_fscale = u64::from_le_bytes(data[35..43].try_into().unwrap());
        let negative_fscale = u64::from_le_bytes(data[43..51].try_into().unwrap());
        let _min_positive_retention = f32::from_le_bytes(data[51..55].try_into().unwrap()); // legacy, discarded
        // bytes 55..63 reserved
        let unseen_fscale = u64::from_le_bytes(data[63..71].try_into().unwrap());

        let mut fscale_lut = [0u64; 255];
        for (i, slot) in fscale_lut.iter_mut().enumerate() {
            let offset = 71 + i * 8;
            *slot = u64::from_le_bytes(data[offset..offset + 8].try_into().unwrap());
        }

        let seeds_start = HEADER_SIZE;
        let seeds_end = seeds_start + depth * 8;
        let weights_start = seeds_end;
        let weights_end = weights_start + width * depth;

        if data.len() < weights_end {
            anyhow::bail!(
                "Bias table data truncated: expected {} bytes, got {}",
                weights_end,
                data.len()
            );
        }

        let mut seeds = Vec::with_capacity(depth);
        for i in 0..depth {
            let offset = seeds_start + i * 8;
            seeds.push(u64::from_le_bytes(
                data[offset..offset + 8].try_into().unwrap(),
            ));
        }

        let mut weights = vec![0i8; width * depth];
        for (i, &b) in data[weights_start..weights_end].iter().enumerate() {
            weights[i] = b as i8;
        }

        let config = CMSConfig {
            width,
            depth,
            k,
            fscale,
        };

        let max_fold_enrichment = if negative_retention > 0.0 {
            positive_retention / negative_retention
        } else {
            f32::INFINITY
        };

        let mut table = Self {
            config,
            seeds,
            weights,
            alpha,
            threshold,
            positive_retention,
            negative_retention,
            max_fold_enrichment,
            target_fscale,
            negative_fscale,
            unseen_fscale,
            fscale_lut,
            frac_max_lut: [0u64; 255],
            cms_pos_counts: None,
            cms_neg_counts: None,
            cms_pos_total: 0.0,
            cms_neg_total: 0.0,
        };
        table.recompute_frac_max_lut();
        Ok(table)
    }

    pub fn weight_stats(&self) -> (f32, f32, f32, f32, usize, usize) {
        let min = *self.weights.iter().min().unwrap_or(&0) as f32 / QUANTIZATION_SCALE;
        let max = *self.weights.iter().max().unwrap_or(&0) as f32 / QUANTIZATION_SCALE;
        let sum: i64 = self.weights.iter().map(|&w| w as i64).sum();
        let mean = sum as f32 / self.weights.len() as f32 / QUANTIZATION_SCALE;
        let variance: f32 = self
            .weights
            .iter()
            .map(|&w| {
                let d = w as f32 / QUANTIZATION_SCALE - mean;
                d * d
            })
            .sum::<f32>()
            / self.weights.len() as f32;
        let positive = self.weights.iter().filter(|&&w| w > 0).count();
        let above_threshold = self
            .weights
            .iter()
            .filter(|&&w| w >= self.threshold)
            .count();
        (min, max, mean, variance.sqrt(), positive, above_threshold)
    }

    pub fn memory_usage(&self) -> usize {
        self.weights.len() + self.seeds.len() * 8
    }

    pub fn threshold_f32(&self) -> f32 {
        self.threshold as f32 / QUANTIZATION_SCALE
    }

    pub fn soft_filter_reference_points(&self) -> Vec<SoftFilterReferencePoint> {
        if !self.is_soft_filter() {
            return Vec::new();
        }

        const WEIGHTS: [i8; 9] = [-127, -64, -32, -16, 0, 16, 32, 64, 127];
        let base = self.config.fscale as f64;

        WEIGHTS
            .iter()
            .map(|&w| {
                let effective_fscale = self.effective_fscale_at(w);
                SoftFilterReferencePoint {
                    weight: w,
                    weight_f32: w as f32 / QUANTIZATION_SCALE,
                    effective_fscale,
                    retention_pct: 100.0 / effective_fscale,
                    vs_base: base / effective_fscale,
                }
            })
            .collect()
    }

    pub fn print_stats(&self) {
        let (min, max, mean, std, positive, above_threshold) = self.weight_stats();
        let total_cells = self.config.width * self.config.depth;
        eprintln!("Hash Bias Table (v5)");
        eprintln!("====================");
        eprintln!("  k-mer size:     {}", self.config.k);
        eprintln!("  fscale:         {}", self.config.fscale);
        eprintln!(
            "  CMS dimensions: {} x {}",
            self.config.width, self.config.depth
        );
        eprintln!("  Smoothing (alpha): {:.1}", self.alpha);
        if self.is_soft_filter() {
            eprintln!("  Filter mode: enrichment LUT");
        } else {
            eprintln!("  Filter mode: hard cutoff");
        }
        eprintln!();

        eprintln!("Calibration:");
        if self.is_soft_filter() {
            eprintln!("  requested target:    {}", self.target_fscale);
            eprintln!("  negative_fscale:     {}", self.negative_fscale_label());
        }
        eprintln!(
            "  positive retention:  {:.2}%",
            self.positive_retention * 100.0
        );
        eprintln!(
            "  negative retention:  {:.2}%",
            self.negative_retention * 100.0
        );
        eprintln!("  fold enrichment:     {:.2}x", self.fold_enrichment());
        if self.is_soft_filter() {
            let target_fs = self.target_fscale as f64;
            let eff_target = self.effective_fscale_target_prior();
            let eff_pos = self.effective_fscale_on_population(self.positive_retention);
            let eff_neg = self.effective_fscale_on_population(self.negative_retention);
            let delta_pct = |achieved: f64| {
                if target_fs > 0.0 {
                    (achieved / target_fs - 1.0) * 100.0
                } else {
                    0.0
                }
            };
            eprintln!(
                "  eff. fscale (target prior): {:.0} ({:+.1}%)",
                eff_target,
                delta_pct(eff_target)
            );
            eprintln!(
                "  eff. fscale (pos):      {:.0} ({:+.1}%)",
                eff_pos,
                delta_pct(eff_pos)
            );
            eprintln!(
                "  eff. fscale (neg):      {:.0} ({:+.1}%)",
                eff_neg,
                delta_pct(eff_neg)
            );
            eprintln!(
                "  eff. fscale (combined): {:.0}",
                self.effective_fscale_combined()
            );
            if self.unseen_fscale > 0 {
                eprintln!("  unseen_fscale:       {}", self.unseen_fscale);
            }

            eprintln!("  reference points:");
            for p in self.soft_filter_reference_points() {
                eprintln!(
                    "    w={:>5.2}: eff={:>10.1}, ret={:>9.4}%, vs_base={:>7.2}x",
                    p.weight_f32, p.effective_fscale, p.retention_pct, p.vs_base
                );
            }
        }
        eprintln!(
            "  threshold:           {:.2} (quantized: {})",
            self.threshold_f32(),
            self.threshold
        );
        eprintln!();

        eprintln!("Weight distribution (clamped to +/-12.70):");
        eprintln!("  min:    {:.2}", min);
        eprintln!("  max:    {:.2}", max);
        eprintln!("  mean:   {:.2}", mean);
        eprintln!("  std:    {:.2}", std);
        eprintln!(
            "  >0:             {} cells ({:.1}%)",
            positive,
            positive as f64 / total_cells as f64 * 100.0
        );
        eprintln!(
            "  >= threshold:   {} cells ({:.1}%)",
            above_threshold,
            above_threshold as f64 / total_cells as f64 * 100.0
        );

        if self.is_soft_filter() {
            self.print_lut_curve(min, max);
        }
    }

    pub fn print_lut_curve(&self, min: f32, max: f32) {
        let base = self.config.fscale as f64;
        let threshold_q = self.threshold as i32;

        let distinct: std::collections::BTreeSet<u64> = self.fscale_lut.iter().copied().collect();
        let lut_min = self.fscale_lut.iter().copied().min().unwrap_or(0);
        let lut_max = self.fscale_lut.iter().copied().max().unwrap_or(0);

        eprintln!();
        eprintln!(
            "LUT response curve (enrichment, {} distinct fscale values):",
            distinct.len()
        );
        eprintln!("  fscale range: {} \u{2014} {}", lut_min, lut_max);

        let start_q = (min * QUANTIZATION_SCALE).floor() as i32;
        let end_q = (max * QUANTIZATION_SCALE).ceil() as i32;

        let mut points: Vec<i32> = (start_q..=end_q).map(|q| q.clamp(-127, 127)).collect();
        if !points.contains(&0) {
            points.push(0);
        }
        if threshold_q != 0 && !points.contains(&threshold_q) {
            points.push(threshold_q);
        }
        points.sort_unstable();
        points.dedup();

        struct Group {
            start_q: i32,
            end_q: i32,
            eff: f64,
        }
        let mut groups: Vec<Group> = Vec::new();
        for &q in &points {
            let w = q.clamp(-128, 127) as i8;
            let eff = self.effective_fscale_at(w);
            match groups.last_mut() {
                Some(g) if (g.eff - eff).abs() < 0.01 => g.end_q = q,
                _ => groups.push(Group {
                    start_q: q,
                    end_q: q,
                    eff,
                }),
            }
        }

        eprintln!();
        eprintln!(
            "  {:>15}  {:>11}  {:>10}  {:>8}",
            "Weight", "Eff.fscale", "Retention", "vs Base"
        );
        eprintln!(
            "  {:>15}  {:>11}  {:>10}  {:>8}",
            "---------------", "----------", "---------", "-------"
        );

        for g in &groups {
            let retention_pct = 100.0 / g.eff;
            let vs_base = base / g.eff;
            let w_start = g.start_q as f64 / 10.0;
            let w_end = g.end_q as f64 / 10.0;

            let weight_col = if g.start_q == g.end_q {
                format!("{:>7.2}       ", w_start)
            } else {
                format!("{:>7.2}→{:<5.2}", w_start, w_end)
            };

            let mut markers: Vec<String> = Vec::new();
            if g.start_q <= 0 && g.end_q >= 0 {
                markers.push("unseen".to_string());
            }
            if threshold_q != 0 && g.start_q <= threshold_q && g.end_q >= threshold_q {
                markers.push("threshold".to_string());
            }
            let count = (g.end_q - g.start_q + 1) as usize;
            if g.start_q != g.end_q {
                markers.push(format!("{} steps", count));
            }

            let marker = if markers.is_empty() {
                String::new()
            } else {
                format!("  <- {}", markers.join(", "))
            };

            eprintln!(
                "  {}  {:>11.1}  {:>9.4}%  {:>7.2}x{}",
                weight_col, g.eff, retention_pct, vs_base, marker
            );
        }
    }
}

fn calibrate_threshold(
    positive: &RawHashCounts,
    negative: &RawHashCounts,
    weights: &[i8],
    seeds: &[u64],
    width: usize,
    target_fold_enrichment: Option<f32>,
) -> Result<CalibrationResult> {
    let sample_hashes = |raw: &RawHashCounts, max_samples: usize| -> Vec<u64> {
        if raw.samples.len() <= max_samples {
            return raw.samples.clone();
        }
        let step = raw.samples.len() / max_samples;
        raw.samples
            .iter()
            .step_by(step)
            .take(max_samples)
            .copied()
            .collect()
    };

    let estimate_weight = |hash: u64| -> i8 {
        let depth = seeds.len();
        (0..depth)
            .map(|row| {
                let mixed = hash.wrapping_mul(seeds[row]);
                let idx = row * width + (mixed as usize % width);
                weights[idx]
            })
            .min()
            .unwrap_or(0)
    };

    let pos_sample_weights: Vec<i8> = sample_hashes(positive, 100_000)
        .iter()
        .map(|&h| estimate_weight(h))
        .collect();
    let neg_sample_weights: Vec<i8> = sample_hashes(negative, 100_000)
        .iter()
        .map(|&h| estimate_weight(h))
        .collect();

    if pos_sample_weights.is_empty() || neg_sample_weights.is_empty() {
        return Ok(CalibrationResult {
            threshold: 0,
            positive_retention: 1.0,
            negative_retention: 1.0,
            fold_enrichment: 1.0,
            max_fold_enrichment: 1.0,
        });
    }

    let mut max_enrichment = 0.0f32;
    let mut max_threshold = 0i8;
    let mut max_pos_ret = 1.0f32;
    let mut max_neg_ret = 1.0f32;

    for t in -127i8..=127i8 {
        let pos_passing = pos_sample_weights.iter().filter(|&&w| w >= t).count();
        let neg_passing = neg_sample_weights.iter().filter(|&&w| w >= t).count();

        let pos_ret = pos_passing as f32 / pos_sample_weights.len() as f32;
        let neg_ret = neg_passing as f32 / neg_sample_weights.len().max(1) as f32;

        if neg_ret < 1e-6 {
            continue;
        }

        let enrichment = pos_ret / neg_ret;
        if enrichment > max_enrichment {
            max_enrichment = enrichment;
            max_threshold = t;
            max_pos_ret = pos_ret;
            max_neg_ret = neg_ret;
        }
    }

    match target_fold_enrichment {
        None => Ok(CalibrationResult {
            threshold: max_threshold,
            positive_retention: max_pos_ret,
            negative_retention: max_neg_ret,
            fold_enrichment: max_enrichment,
            max_fold_enrichment: max_enrichment,
        }),
        Some(target) => {
            if target > max_enrichment {
                return Ok(CalibrationResult {
                    threshold: max_threshold,
                    positive_retention: max_pos_ret,
                    negative_retention: max_neg_ret,
                    fold_enrichment: max_enrichment,
                    max_fold_enrichment: max_enrichment,
                });
            }

            let mut best_threshold = 0i8;
            let mut best_diff = f32::MAX;
            let mut best_pos_ret = 1.0f32;
            let mut best_neg_ret = 1.0f32;

            for t in -127i8..=127i8 {
                let pos_passing = pos_sample_weights.iter().filter(|&&w| w >= t).count();
                let neg_passing = neg_sample_weights.iter().filter(|&&w| w >= t).count();

                let pos_ret = pos_passing as f32 / pos_sample_weights.len() as f32;
                let neg_ret = neg_passing as f32 / neg_sample_weights.len().max(1) as f32;

                if neg_ret < 1e-6 {
                    continue;
                }

                let enrichment = pos_ret / neg_ret;
                let diff = (enrichment - target).abs();

                if diff < best_diff {
                    best_diff = diff;
                    best_threshold = t;
                    best_pos_ret = pos_ret;
                    best_neg_ret = neg_ret;
                }
            }

            Ok(CalibrationResult {
                threshold: best_threshold,
                positive_retention: best_pos_ret,
                negative_retention: best_neg_ret,
                fold_enrichment: if best_neg_ret > 0.0 {
                    best_pos_ret / best_neg_ret
                } else {
                    f32::INFINITY
                },
                max_fold_enrichment: max_enrichment,
            })
        }
    }
}

fn format_number(n: u64) -> String {
    if n >= 1_000_000_000 {
        format!("{:.2}G", n as f64 / 1_000_000_000.0)
    } else if n >= 1_000_000 {
        format!("{:.2}M", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.2}K", n as f64 / 1_000.0)
    } else {
        format!("{}", n)
    }
}

pub fn format_bp(bp: u64) -> String {
    if bp >= 1_000_000_000 {
        format!("{:.2} Gbp", bp as f64 / 1_000_000_000.0)
    } else if bp >= 1_000_000 {
        format!("{:.2} Mbp", bp as f64 / 1_000_000.0)
    } else if bp >= 1_000 {
        format!("{:.2} Kbp", bp as f64 / 1_000.0)
    } else {
        format!("{} bp", bp)
    }
}

pub const BIAS_TABLE_SERIALIZED_SIZE: usize =
    2111 + DEFAULT_CMS_DEPTH * 8 + DEFAULT_CMS_WIDTH * DEFAULT_CMS_DEPTH;

impl PartialEq for HashBiasTable {
    fn eq(&self, other: &Self) -> bool {
        self.config.k == other.config.k
            && self.config.fscale == other.config.fscale
            && self.config.width == other.config.width
            && self.config.depth == other.config.depth
            && self.alpha == other.alpha
            && self.threshold == other.threshold
            && self.positive_retention == other.positive_retention
            && self.negative_retention == other.negative_retention
            && self.target_fscale == other.target_fscale
            && self.negative_fscale == other.negative_fscale
            && self.unseen_fscale == other.unseen_fscale
            && self.fscale_lut == other.fscale_lut
            && self.seeds == other.seeds
            && self.weights == other.weights
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn assert_approx_eq(actual: f64, expected: f64, tol: f64) {
        let diff = (actual - expected).abs();
        assert!(
            diff <= tol,
            "expected {expected}, got {actual} (diff {diff} > tol {tol})"
        );
    }

    fn create_fasta(sequences: &[&str]) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        for (i, seq) in sequences.iter().enumerate() {
            writeln!(file, ">seq_{}", i).unwrap();
            writeln!(file, "{}", seq).unwrap();
        }
        file
    }

    #[test]
    fn test_cms_basic() {
        let mut cms = CountMinSketch::new(1024, 5);
        let hash = 0x12345678u64;

        assert_eq!(cms.estimate(hash), 0);

        cms.increment(hash);
        assert_eq!(cms.estimate(hash), 1);

        for _ in 0..9 {
            cms.increment(hash);
        }
        assert_eq!(cms.estimate(hash), 10);
    }

    #[test]
    fn test_cms_collision_handling() {
        let mut cms = CountMinSketch::new(16, 5);

        for i in 0..100u64 {
            cms.increment(i);
        }

        for i in 0..100u64 {
            assert!(cms.estimate(i) >= 1);
        }
    }

    #[test]
    fn test_raw_hash_counts_build() {
        let fasta = create_fasta(&["ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"]);
        let config = CMSConfig {
            width: 1024,
            depth: 3,
            k: 11,
            fscale: 1,
        };

        let raw = RawHashCounts::build(
            &[fasta.path()],
            config,
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();
        assert!(raw.total > 0);
    }

    #[test]
    fn test_hash_bias_table_build() {
        let pos = create_fasta(&[
            "ATATATATATATATATATATATATATATATATATATATAT",
            "TATATATATATATATATATATATATATATATATATATAT",
        ]);
        let neg = create_fasta(&[
            "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC",
            "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG",
        ]);

        let config = CMSConfig {
            width: 1024,
            depth: 3,
            k: 11,
            fscale: 1,
        };

        let pos_raw = RawHashCounts::build(
            &[pos.path()],
            config.clone(),
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();
        let neg_raw = RawHashCounts::build(
            &[neg.path()],
            config,
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 0, 0).unwrap();
        assert!(table.threshold >= -127);
    }

    #[test]
    fn test_hash_bias_table_save_load() {
        let pos = create_fasta(&["ATATATATATATATATATATATATATATATATATATATAT"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"]);

        let config = CMSConfig {
            width: 1024,
            depth: 3,
            k: 11,
            fscale: 10,
        };

        let pos_raw = RawHashCounts::build(
            &[pos.path()],
            config.clone(),
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();
        let neg_raw = RawHashCounts::build(
            &[neg.path()],
            config,
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 0, 0).unwrap();

        let output = NamedTempFile::new().unwrap();
        table.save(output.path()).unwrap();

        let loaded = HashBiasTable::load(output.path()).unwrap();
        assert_eq!(table.config.k, loaded.config.k);
        assert_eq!(table.threshold, loaded.threshold);
        assert_eq!(table.weights, loaded.weights);
    }

    #[test]
    fn test_hash_bias_table_bytes_roundtrip() {
        let pos = create_fasta(&["ATATATATATATATATATATATATATATATATATATATAT"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"]);

        let config = CMSConfig {
            width: 512,
            depth: 3,
            k: 11,
            fscale: 10,
        };

        let pos_raw = RawHashCounts::build(
            &[pos.path()],
            config.clone(),
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();
        let neg_raw = RawHashCounts::build(
            &[neg.path()],
            config,
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 0, 0).unwrap();

        let bytes = table.to_bytes();
        let loaded = HashBiasTable::from_bytes(&bytes).unwrap();

        assert_eq!(table, loaded);
    }

    #[test]
    fn test_passes_filter() {
        let pos = create_fasta(&["ATATATATATATATATATATATATATATATATATATATAT"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"]);

        let config = CMSConfig {
            width: 1024,
            depth: 3,
            k: 11,
            fscale: 1,
        };

        let pos_raw = RawHashCounts::build(
            &[pos.path()],
            config.clone(),
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();
        let neg_raw = RawHashCounts::build(
            &[neg.path()],
            config,
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 0, 0).unwrap();

        let mut passed = 0;
        let mut failed = 0;
        for h in 0..1000u64 {
            if table.passes_filter(h) {
                passed += 1;
            } else {
                failed += 1;
            }
        }

        assert!(passed > 0 || failed > 0);
    }

    #[test]
    fn test_maximize_fold_enrichment() {
        let pos = create_fasta(&[
            "ATATATATATATATATATATATATATATATATATATATAT",
            "TATATATATATATATATATATATATATATATATATATAT",
        ]);
        let neg = create_fasta(&[
            "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC",
            "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG",
        ]);

        let config = CMSConfig {
            width: 1024,
            depth: 3,
            k: 11,
            fscale: 1,
        };

        let pos_raw = RawHashCounts::build(
            &[pos.path()],
            config.clone(),
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();
        let neg_raw = RawHashCounts::build(
            &[neg.path()],
            config,
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 0, 0).unwrap();
        assert!(table.threshold >= -127);
        assert!(table.fold_enrichment() >= 1.0);
    }

    #[test]
    fn test_create_unified() {
        let pos = create_fasta(&["ATATATATATATATATATATATATATATATATATATATAT"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"]);

        let config = BiasCreateConfig {
            cms: CMSConfig {
                width: 1024,
                depth: 3,
                k: 11,
                fscale: 1,
            },
            alpha: 1.0,
            target_fscale: None,
            negative_fscale: None,
            unseen_fscale: None,
        };

        let table = HashBiasTable::create(&[pos.path()], &[neg.path()], &config, None).unwrap();

        assert!(table.threshold >= -127);
        assert!(table.fold_enrichment() >= 1.0);
    }

    #[test]
    fn test_hard_cutoff_equivalence() {
        let pos = create_fasta(&[
            "ATATATATATATATATATATATATATATATATATATATAT",
            "TATATATATATATATATATATATATATATATATATATAT",
        ]);
        let neg = create_fasta(&[
            "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC",
            "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG",
        ]);

        let config = CMSConfig {
            width: 1024,
            depth: 3,
            k: 11,
            fscale: 1,
        };

        let pos_raw = RawHashCounts::build(
            &[pos.path()],
            config.clone(),
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();
        let neg_raw = RawHashCounts::build(
            &[neg.path()],
            config,
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 0, 0).unwrap();
        assert!(!table.is_soft_filter());

        for h in 0..10_000u64 {
            let expected = table.weight(h) >= table.threshold;
            assert_eq!(
                table.passes_filter(h),
                expected,
                "Mismatch for hash {} (weight={}, threshold={})",
                h,
                table.weight(h),
                table.threshold
            );
        }
    }

    #[test]
    fn test_soft_filter_deterministic() {
        let pos = create_fasta(&["ATATATATATATATATATATATATATATATATATATATAT"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"]);

        let config = CMSConfig {
            width: 1024,
            depth: 3,
            k: 11,
            fscale: 1,
        };

        let pos_raw = RawHashCounts::build(
            &[pos.path()],
            config.clone(),
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();
        let neg_raw = RawHashCounts::build(
            &[neg.path()],
            config,
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 1, 1000, 0).unwrap();
        assert!(table.is_soft_filter());

        for h in 0..10_000u64 {
            let first = table.passes_filter(h);
            let second = table.passes_filter(h);
            assert_eq!(first, second, "Non-deterministic for hash {}", h);
        }
    }

    #[test]
    fn test_soft_filter_bytes_roundtrip() {
        let pos = create_fasta(&["ATATATATATATATATATATATATATATATATATATATAT"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"]);

        let config = CMSConfig {
            width: 512,
            depth: 3,
            k: 11,
            fscale: 10,
        };

        let pos_raw = RawHashCounts::build(
            &[pos.path()],
            config.clone(),
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();
        let neg_raw = RawHashCounts::build(
            &[neg.path()],
            config,
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 10, 5000, 0).unwrap();
        assert!(table.is_soft_filter());

        let bytes = table.to_bytes();
        let loaded = HashBiasTable::from_bytes(&bytes).unwrap();
        assert_eq!(table, loaded);
    }

    #[test]
    fn test_w0_bucket_data_driven_without_unbiased() {
        let pos = create_fasta(&["ATATATATATATATATATATATATATATATATATATATAT"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"]);

        let config = CMSConfig {
            width: 1024,
            depth: 3,
            k: 11,
            fscale: 10,
        };

        let pos_raw = RawHashCounts::build(
            &[pos.path()],
            config.clone(),
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();
        let neg_raw = RawHashCounts::build(
            &[neg.path()],
            config,
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 100, 10_000, 0).unwrap();
        assert!(table.is_soft_filter());

        // When unseen_fscale is not set (0), w=0 should be anchored at target_fscale
        let at_0 = table.effective_fscale_at(0);
        assert_approx_eq(at_0, 100.0, 1e-6);
    }

    #[test]
    fn test_build_rejects_invalid_invariants() {
        let pos = create_fasta(&["ATATATATATATATATATATATATATATATATATATATAT"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"]);

        let config = CMSConfig {
            width: 512,
            depth: 3,
            k: 11,
            fscale: 100,
        };

        let pos_raw = RawHashCounts::build(
            &[pos.path()],
            config.clone(),
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();
        let neg_raw = RawHashCounts::build(
            &[neg.path()],
            config,
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();

        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 50, 1000, 0);
        assert!(result.is_err());

        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 100, 0, 0);
        assert!(result.is_err());

        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 1000, 0);
        assert!(result.is_err());

        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 100, 100, 0);
        assert!(result.is_err());

        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 100, 5000, 0);
        assert!(result.is_ok());

        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 100, 5000, 50);
        assert!(result.is_err());

        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 100, 5000, 6000);
        assert!(result.is_err());

        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 0, 500);
        assert!(result.is_err());
    }

    #[test]
    fn test_custom_unbiased_fscale() {
        let pos = create_fasta(&["ATATATATATATATATATATATATATATATATATATATAT"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"]);

        let config = CMSConfig {
            width: 1024,
            depth: 3,
            k: 11,
            fscale: 1,
        };

        let pos_raw = RawHashCounts::build(
            &[pos.path()],
            config.clone(),
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();
        let neg_raw = RawHashCounts::build(
            &[neg.path()],
            config,
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 10, 100, 10).unwrap();
        assert!(table.is_soft_filter());
        assert_eq!(table.unseen_fscale, 10);

        assert_approx_eq(table.effective_fscale_at(0), 10.0, 1e-6);

        let geom_mean = (10.0_f64 * 100.0_f64).sqrt();
        assert!((table.effective_fscale_at(0) - geom_mean).abs() > 1.0);
    }

    #[test]
    fn test_lut_direct_lookup() {
        let mut fscale_lut = [0u64; 255];
        for (i, slot) in fscale_lut.iter_mut().enumerate() {
            *slot = 10_000u64 - ((10_000u64 - 100) * i as u64 / 254);
        }
        fscale_lut[127] = 1_000;

        let mut table = HashBiasTable {
            config: CMSConfig {
                width: 1,
                depth: 1,
                k: 11,
                fscale: 100,
            },
            seeds: vec![1],
            weights: vec![0],
            alpha: 1.0,
            threshold: 10,
            positive_retention: 0.0,
            negative_retention: 0.0,
            max_fold_enrichment: 0.0,
            target_fscale: 100,
            negative_fscale: 10_000,
            unseen_fscale: 1_000,
            fscale_lut,
            frac_max_lut: [0u64; 255],
            cms_pos_counts: None,
            cms_neg_counts: None,
            cms_pos_total: 0.0,
            cms_neg_total: 0.0,
        };
        table.recompute_frac_max_lut();

        assert_approx_eq(table.effective_fscale_at(0), 1_000.0, 1e-6);

        assert_approx_eq(table.effective_fscale_at(-127), 10_000.0, 1e-6);

        assert_approx_eq(table.effective_fscale_at(127), 100.0, 1.0);

        let below = table.effective_fscale_at(-50);
        let above = table.effective_fscale_at(50);
        assert!(
            below > above,
            "expected monotonic decrease: below={}, above={}",
            below,
            above
        );
    }

    #[test]
    fn test_lut_bounds_respected() {
        let pos = create_fasta(&["ATATATATATATATATATATATAT"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGC"]);
        let config = CMSConfig {
            width: 1024,
            depth: 3,
            k: 11,
            fscale: 100,
        };
        let pos_raw = RawHashCounts::build(
            &[pos.path()],
            config.clone(),
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();
        let neg_raw = RawHashCounts::build(
            &[neg.path()],
            config,
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 100, 10_000, 0).unwrap();
        assert!(table.is_soft_filter());

        for (i, &fs) in table.fscale_lut.iter().enumerate() {
            assert!(
                fs >= table.config.fscale && fs <= table.negative_fscale,
                "fscale_lut[{}] = {} out of bounds [{}, {}]",
                i,
                fs,
                table.config.fscale,
                table.negative_fscale
            );
        }
    }

    #[test]
    fn test_tracks_target_effective_fscale() {
        let mut p_pos = [HISTOGRAM_EPS; 255];
        let mut p_neg = [HISTOGRAM_EPS; 255];
        let pos_counts = [100u64; 255];
        let neg_counts = [100u64; 255];

        for i in 200..255 {
            p_pos[i] = 1.0 / 55.0;
        }
        for i in 0..55 {
            p_neg[i] = 1.0 / 55.0;
        }

        let pos_sum: f64 = p_pos.iter().sum();
        let neg_sum: f64 = p_neg.iter().sum();
        let mut p_pos_norm = p_pos;
        let mut p_neg_norm = p_neg;
        for v in p_pos_norm.iter_mut() {
            *v /= pos_sum;
        }
        for v in p_neg_norm.iter_mut() {
            *v /= neg_sum;
        }

        let base_fscale = 100u64;
        let target_fscale = 500u64;
        let negative_fscale = 10000u64;
        let unseen_fscale = 500u64;
        let p_target = [1.0 / 255.0; 255];

        let lut = target_fscale_lut(
            &p_pos_norm,
            &p_neg_norm,
            &p_target,
            &pos_counts,
            &neg_counts,
            base_fscale,
            target_fscale,
            negative_fscale,
            unseen_fscale,
        )
        .unwrap();

        let actual_target_ret: f64 = (0..255)
            .map(|i| {
                let r = (base_fscale as f64 / lut[i] as f64).min(1.0);
                p_target[i] * r
            })
            .sum();
        let actual_eff_target = base_fscale as f64 / actual_target_ret;
        let rel_err = ((actual_eff_target - target_fscale as f64) / target_fscale as f64).abs();
        assert!(
            rel_err <= 0.03,
            "target-prior effective fscale too far from target: target={}, got={:.2}, rel_err={:.4}",
            target_fscale,
            actual_eff_target,
            rel_err
        );

        // Positive buckets should have lower fscale than negative buckets
        let avg_pos_fscale: f64 = (200..255).map(|i| lut[i] as f64).sum::<f64>() / 55.0;
        let avg_neg_fscale: f64 = (0..55).map(|i| lut[i] as f64).sum::<f64>() / 55.0;
        assert!(
            avg_pos_fscale < avg_neg_fscale,
            "positive buckets should have lower fscale: pos={avg_pos_fscale}, neg={avg_neg_fscale}"
        );

        // Endpoints should be anchored
        assert_eq!(lut[254], base_fscale);
        assert_eq!(lut[0], negative_fscale);
        assert_eq!(lut[127], unseen_fscale);

        // Should have many distinct fscale values (anti-regression for flat LUT)
        let distinct: std::collections::HashSet<u64> = lut.iter().copied().collect();
        assert!(
            distinct.len() >= 10,
            "LUT should have many distinct values, got {}",
            distinct.len()
        );
    }

    #[test]
    fn test_target_fscale_lut_integration() {
        let pos = create_fasta(&[
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
            "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
            "AACGTAACGTAACGTAACGTAACGTAACGTAACGTAACGT",
            "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT",
        ]);
        let neg = create_fasta(&[
            "TTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAAT",
            "CCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGC",
            "GGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCG",
            "AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTT",
        ]);

        let config = CMSConfig {
            width: 1024,
            depth: 3,
            k: 11,
            fscale: 1,
        };

        let pos_raw = RawHashCounts::build(
            &[pos.path()],
            config.clone(),
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();
        let neg_raw = RawHashCounts::build(
            &[neg.path()],
            config,
            &AtomicU64::new(0),
            &AtomicU64::new(0),
        )
        .unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 2, 100, 0).unwrap();
        assert!(table.is_soft_filter());
        assert!(table.positive_retention > 0.0);
        assert!(
            table.fold_enrichment() >= 1.05,
            "fold enrichment should be > 1.05, got {}",
            table.fold_enrichment()
        );
    }

    #[test]
    fn test_lut_smooth_gradient_max_step() {
        let mut p_pos = [HISTOGRAM_EPS; 255];
        let mut p_neg = [HISTOGRAM_EPS; 255];
        let pos_counts = [100u64; 255];
        let neg_counts = [100u64; 255];

        for i in 180..255 {
            p_pos[i] = 1.0 / 75.0;
        }
        for i in 0..75 {
            p_neg[i] = 1.0 / 75.0;
        }

        let pos_sum: f64 = p_pos.iter().sum();
        let neg_sum: f64 = p_neg.iter().sum();
        for v in p_pos.iter_mut() {
            *v /= pos_sum;
        }
        for v in p_neg.iter_mut() {
            *v /= neg_sum;
        }

        let base_fscale = 100u64;
        let negative_fscale = 10000u64;
        let target_fscale = 1000u64;
        let p_target = [1.0 / 255.0; 255];

        let lut = target_fscale_lut(
            &p_pos,
            &p_neg,
            &p_target,
            &pos_counts,
            &neg_counts,
            base_fscale,
            target_fscale,
            negative_fscale,
            target_fscale,
        )
        .unwrap();

        let ln_base = (base_fscale as f64).ln();
        let ln_neg = (negative_fscale as f64).ln();
        let s_max = (ln_neg - ln_base) / LUT_MAX_LOG_STEP_DIVISOR;
        let w0_idx = HashBiasTable::weight_to_index(0);

        for i in 1..255 {
            if i == w0_idx || i == w0_idx + 1 {
                continue;
            }
            let log_prev = (lut[i - 1] as f64).ln();
            let log_cur = (lut[i] as f64).ln();
            let diff = (log_prev - log_cur).abs();
            assert!(
                diff <= s_max + 0.05,
                "Adjacent log-step too large at index {}: |ln({}) - ln({})| = {:.4}, s_max = {:.4}",
                i,
                lut[i - 1],
                lut[i],
                diff,
                s_max
            );
        }

        assert_eq!(lut[254], base_fscale);
        assert_eq!(lut[0], negative_fscale);

        let avg_pos_fscale: f64 = (180..255).map(|i| lut[i] as f64).sum::<f64>() / 75.0;
        let avg_neg_fscale: f64 = (0..75).map(|i| lut[i] as f64).sum::<f64>() / 75.0;
        assert!(
            avg_pos_fscale < avg_neg_fscale,
            "positive bins should have lower fscale: pos={avg_pos_fscale}, neg={avg_neg_fscale}"
        );

        let distinct: std::collections::HashSet<u64> = lut.iter().copied().collect();
        assert!(
            distinct.len() >= 10,
            "LUT should have many distinct values, got {}",
            distinct.len()
        );
    }
}
