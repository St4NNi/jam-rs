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
    pub target_positive_retention: Option<f32>,
    pub target_negative_retention: Option<f32>,
    pub positive_fscale: Option<u64>,
    pub negative_fscale: Option<u64>,
    pub unbiased_fscale: Option<u64>,
    pub lambda: Option<f32>,
    pub gamma: f32,
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
    pub positive_fscale: u64,
    pub negative_fscale: u64,
    pub lambda: f32,
    pub gamma: f32,
    pub unbiased_fscale: u64,
    pub fscale_lut: [u64; 255],
    frac_max_lut: [u64; 255],
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

// ── LUT Optimizer ──────────────────────────────────────────────────────────

/// Build a histogram of hashes per weight bucket, normalized with Laplace smoothing.
/// Returns `[f64; 255]` where index `i` corresponds to quantized weight `i - 127`.
fn histogram_weights(hashes: &[u64], bias_table: &HashBiasTable) -> [f64; 255] {
    let mut counts = [0u64; 255];
    for &h in hashes {
        let w = bias_table.weight(h);
        let idx = HashBiasTable::weight_to_index(w);
        counts[idx] += 1;
    }
    let total: f64 = counts.iter().sum::<u64>() as f64 + 255.0; // Laplace smoothing
    let mut hist = [0.0f64; 255];
    for i in 0..255 {
        hist[i] = (counts[i] as f64 + 1.0) / total;
    }
    hist
}

/// Pool Adjacent Violators Algorithm — enforce monotonically non-decreasing.
/// Operates in-place on a slice. Clamps values to `[lo, hi]`.
fn pava_increasing(r: &mut [f64], lo: f64, hi: f64) {
    let n = r.len();
    if n == 0 {
        return;
    }
    // Clamp first
    for v in r.iter_mut() {
        *v = v.clamp(lo, hi);
    }
    // Weighted PAVA using block representation
    // Each block stores (sum_of_values, count)
    let mut blocks: Vec<(f64, usize)> = r.iter().map(|&v| (v, 1)).collect();
    let mut i = 0;
    while i < blocks.len() - 1 {
        let mean_i = blocks[i].0 / blocks[i].1 as f64;
        let mean_next = blocks[i + 1].0 / blocks[i + 1].1 as f64;
        if mean_i > mean_next {
            // Pool blocks i and i+1
            blocks[i].0 += blocks[i + 1].0;
            blocks[i].1 += blocks[i + 1].1;
            blocks.remove(i + 1);
            // Step back
            if i > 0 {
                i -= 1;
            }
        } else {
            i += 1;
        }
    }
    // Expand blocks back into r
    let mut idx = 0;
    for (sum, count) in &blocks {
        let val = (sum / *count as f64).clamp(lo, hi);
        for _ in 0..*count {
            r[idx] = val;
            idx += 1;
        }
    }
}

/// Optimize the LUT via projected gradient ascent.
///
/// - `p_pos[i]`, `p_neg[i]`: histogram weights for positive/negative samples
/// - `lambda`: penalty on negative retention
/// - `gamma`: smoothness penalty coefficient
/// - `r_min`, `r_max`: retention bounds (1/negative_fscale, 1/positive_fscale)
/// - `anchor_idx`: if Some, fix that index to `anchor_val` (w=0 pinned to unbiased_fscale)
/// - Returns: `[u64; 255]` fscale values (fscale = round(1/r))
fn optimize_lut(
    p_pos: &[f64; 255],
    p_neg: &[f64; 255],
    lambda: f64,
    gamma: f64,
    r_min: f64,
    r_max: f64,
    anchor_idx: Option<usize>,
    anchor_val: f64,
) -> [u64; 255] {
    let mut r = [0.0f64; 255];

    // Initialize linearly from r_min to r_max
    for i in 0..255 {
        r[i] = r_min + (r_max - r_min) * (i as f64 / 254.0);
    }

    // Pin anchor
    if let Some(idx) = anchor_idx {
        r[idx] = anchor_val;
    }

    let lr = 0.01;
    let max_iters = 2000;

    for _iter in 0..max_iters {
        let mut grad = [0.0f64; 255];

        // Objective gradient: P+(i) - lambda * P-(i)
        for i in 0..255 {
            grad[i] = p_pos[i] - lambda * p_neg[i];
        }

        // Smoothness penalty gradient: -2*gamma*(2*r[i] - r[i-1] - r[i+1])
        if gamma > 0.0 {
            for i in 0..255 {
                let left = if i > 0 { r[i - 1] } else { r[i] };
                let right = if i < 254 { r[i + 1] } else { r[i] };
                grad[i] -= 2.0 * gamma * (2.0 * r[i] - left - right);
            }
        }

        // Zero anchor gradient
        if let Some(idx) = anchor_idx {
            grad[idx] = 0.0;
        }

        // Gradient step
        let mut delta = 0.0f64;
        for i in 0..255 {
            let old = r[i];
            r[i] += lr * grad[i];
            r[i] = r[i].clamp(r_min, r_max);
            delta += (r[i] - old).abs();
        }

        // Restore anchor
        if let Some(idx) = anchor_idx {
            r[idx] = anchor_val;
        }

        // PAVA projection: run on left/right halves independently with anchor-aware bounds
        if let Some(idx) = anchor_idx {
            // Left half: indices 0..=idx, must be in [r_min, anchor_val]
            pava_increasing(&mut r[..=idx], r_min, anchor_val);
            // Right half: indices idx..255, must be in [anchor_val, r_max]
            pava_increasing(&mut r[idx..], anchor_val, r_max);
            // Restore anchor (PAVA may have nudged it)
            r[idx] = anchor_val;
        } else {
            pava_increasing(&mut r[..], r_min, r_max);
        }

        // Convergence check
        if delta < 1e-10 {
            break;
        }
    }

    // Convert retention → fscale: fscale = round(1/r), clamped
    let fscale_min = (1.0 / r_max).round() as u64;
    let fscale_max = (1.0 / r_min).round() as u64;
    let mut lut = [0u64; 255];
    for i in 0..255 {
        let fs = if r[i] <= 0.0 {
            fscale_max
        } else {
            (1.0 / r[i]).round() as u64
        };
        lut[i] = fs.clamp(fscale_min, fscale_max);
    }
    lut
}

/// Two-phase lambda search: coarse grid scan + binary search.
/// Returns lambda that achieves the target retention for the given sample type.
///
/// `target_fn` computes retention given a lambda value. It should return the
/// retention metric to match against `target`.
fn search_lambda<F>(target: f64, max_lambda: f64, mut target_fn: F) -> f64
where
    F: FnMut(f64) -> f64,
{
    // Phase A: coarse grid scan to bracket the target
    let grid_steps = 50;
    let mut best_lo = 0.0f64;
    let mut best_hi = max_lambda;
    let mut val_lo = target_fn(best_lo);
    let mut found_bracket = false;

    for i in 1..=grid_steps {
        let lam = max_lambda * i as f64 / grid_steps as f64;
        let val = target_fn(lam);
        if (val_lo >= target && val <= target) || (val_lo <= target && val >= target) {
            best_lo = max_lambda * (i - 1) as f64 / grid_steps as f64;
            best_hi = lam;
            found_bracket = true;
            break;
        }
        val_lo = val;
    }

    if !found_bracket {
        // If no bracket found, return the endpoint closer to target
        let val_0 = target_fn(0.0);
        let val_max = target_fn(max_lambda);
        return if (val_0 - target).abs() < (val_max - target).abs() {
            0.0
        } else {
            max_lambda
        };
    }

    // Phase B: binary search within bracket
    for _ in 0..50 {
        let mid = (best_lo + best_hi) / 2.0;
        let val_mid = target_fn(mid);
        let val_at_lo = target_fn(best_lo);

        if (val_at_lo >= target) == (val_mid >= target) {
            best_lo = mid;
        } else {
            best_hi = mid;
        }

        if (best_hi - best_lo) < 1e-6 {
            break;
        }
    }

    (best_lo + best_hi) / 2.0
}

impl HashBiasTable {
    pub fn create(
        positive_paths: &[&Path],
        negative_paths: &[&Path],
        config: &BiasCreateConfig,
        progress: Option<ProgressBar>,
    ) -> Result<Self> {
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
            config.target_positive_retention,
            config.target_negative_retention,
            config.positive_fscale.unwrap_or(0),
            config.negative_fscale.unwrap_or(0),
            config.unbiased_fscale.unwrap_or(0),
            config.lambda,
            config.gamma,
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
        target_positive_retention: Option<f32>,
        target_negative_retention: Option<f32>,
        positive_fscale: u64,
        negative_fscale: u64,
        unbiased_fscale: u64,
        lambda_arg: Option<f32>,
        gamma: f32,
    ) -> Result<Self> {
        validate_cms_compatibility(positive, negative)?;

        if !alpha.is_finite() || alpha <= 0.0 {
            anyhow::bail!("alpha must be finite and > 0, got {}", alpha);
        }

        if let Some(target) = target_positive_retention
            && (!target.is_finite() || !(0.0..=1.0).contains(&target))
        {
            anyhow::bail!("target_positive_retention must be in [0,1], got {}", target);
        }

        if let Some(target) = target_negative_retention
            && (!target.is_finite() || !(0.0..=1.0).contains(&target))
        {
            anyhow::bail!("target_negative_retention must be in [0,1], got {}", target);
        }

        if let Some(lam) = lambda_arg {
            if !lam.is_finite() || lam < 0.0 {
                anyhow::bail!("lambda must be finite and >= 0, got {}", lam);
            }
        }

        if !gamma.is_finite() || gamma < 0.0 {
            anyhow::bail!("gamma must be finite and >= 0, got {}", gamma);
        }

        // Enforce soft-filter invariants
        if positive_fscale > 0 || negative_fscale > 0 {
            if positive_fscale == 0 || negative_fscale == 0 {
                anyhow::bail!("positive_fscale and negative_fscale must be set together");
            }
            if positive_fscale < positive.config.fscale {
                anyhow::bail!(
                    "positive_fscale ({}) must be >= global fscale ({})",
                    positive_fscale,
                    positive.config.fscale
                );
            }
            if negative_fscale <= positive_fscale {
                anyhow::bail!(
                    "negative_fscale ({}) must be > positive_fscale ({})",
                    negative_fscale,
                    positive_fscale
                );
            }
            if unbiased_fscale > 0 {
                if unbiased_fscale < positive_fscale {
                    anyhow::bail!(
                        "unbiased_fscale ({}) must be >= positive_fscale ({})",
                        unbiased_fscale,
                        positive_fscale
                    );
                }
                if negative_fscale != u64::MAX && unbiased_fscale > negative_fscale {
                    anyhow::bail!(
                        "unbiased_fscale ({}) must be <= negative_fscale ({})",
                        unbiased_fscale,
                        negative_fscale
                    );
                }
            }
        } else if unbiased_fscale > 0 {
            anyhow::bail!(
                "--unbiased-fscale requires soft filter mode (set --positive-fscale and --negative-fscale)"
            );
        }

        let width = positive.config.width;
        let depth = positive.config.depth;
        let seeds = positive.cms.seeds().to_vec();

        let pos_counts = positive.cms.counts();
        let neg_counts = negative.cms.counts();

        let pos_total = positive.total as f64;
        let neg_total = negative.total as f64;

        let mut weights = vec![0i8; width * depth];

        if pos_total > 0.0 && neg_total > 0.0 {
            let scale = pos_total.max(neg_total);

            for i in 0..(width * depth) {
                let norm_pos = (pos_counts[i] as f64 / pos_total) * scale;
                let norm_neg = (neg_counts[i] as f64 / neg_total) * scale;
                let adj_neg = (norm_neg - norm_pos).max(0.0) as f32;
                let norm_pos_f32 = norm_pos as f32;

                let log_ratio = ((norm_pos_f32 + alpha) / (adj_neg + alpha)).ln();
                let quantized = (log_ratio * QUANTIZATION_SCALE).clamp(-127.0, 127.0) as i8;
                weights[i] = quantized;
            }
        }

        // Calibrate threshold (informational for soft filter, functional for hard cutoff)
        let calibration = calibrate_threshold(
            positive, negative, &weights, &seeds, width,
            None, // always maximize — threshold is informational for soft mode
        )?;

        // Initialize with default LUT (all base fscale)
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
            positive_fscale,
            negative_fscale,
            lambda: lambda_arg.unwrap_or(0.0),
            gamma,
            unbiased_fscale,
            fscale_lut,
            frac_max_lut,
        };

        // LUT optimization when soft filter is active
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

            // Histogram positive/negative hashes by weight bucket
            let p_pos = histogram_weights(&pos_in_range, &table);
            let p_neg = histogram_weights(&neg_in_range, &table);

            // Retention bounds from fscale endpoints
            let r_min = 1.0 / negative_fscale as f64; // lowest retention (most filtered)
            let r_max = 1.0 / positive_fscale as f64; // highest retention (least filtered)

            // Anchor w=0 to unbiased_fscale if set, else geometric mean
            let anchor_idx = Some(127usize); // w=0
            let anchor_fscale = if unbiased_fscale > 0 {
                unbiased_fscale as f64
            } else {
                (positive_fscale as f64 * negative_fscale as f64).sqrt()
            };
            let anchor_val = 1.0 / anchor_fscale;

            // Determine lambda
            let gamma_f64 = gamma as f64;
            let resolved_lambda = if let Some(lam) = lambda_arg {
                lam as f64
            } else if let Some(target) = target_negative_retention {
                // Auto-derive lambda to hit negative retention target
                let target_f64 = target as f64;
                search_lambda(target_f64, 1000.0, |lam| {
                    let lut = optimize_lut(
                        &p_pos, &p_neg, lam, gamma_f64, r_min, r_max, anchor_idx, anchor_val,
                    );
                    let mut trial = table.clone();
                    trial.fscale_lut = lut;
                    trial.recompute_frac_max_lut();
                    let (_, neg_ret) =
                        Self::compute_retention(&trial, &pos_in_range, &neg_in_range);
                    neg_ret as f64
                })
            } else if let Some(target) = target_positive_retention {
                // Auto-derive lambda to hit positive retention target
                let target_f64 = target as f64;
                search_lambda(target_f64, 1000.0, |lam| {
                    let lut = optimize_lut(
                        &p_pos, &p_neg, lam, gamma_f64, r_min, r_max, anchor_idx, anchor_val,
                    );
                    let mut trial = table.clone();
                    trial.fscale_lut = lut;
                    trial.recompute_frac_max_lut();
                    let (pos_ret, _) =
                        Self::compute_retention(&trial, &pos_in_range, &neg_in_range);
                    pos_ret as f64
                })
            } else {
                // Default lambda
                10.0
            };

            table.lambda = resolved_lambda as f32;

            // Run optimizer with resolved lambda
            table.fscale_lut = optimize_lut(
                &p_pos,
                &p_neg,
                resolved_lambda,
                gamma_f64,
                r_min,
                r_max,
                anchor_idx,
                anchor_val,
            );

            table.recompute_frac_max_lut();

            // Compute actual retention
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

    /// Compute (positive_retention, negative_retention) using the soft filter.
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
        (0..self.config.depth)
            .map(|row| self.weights[self.index(row, hash)])
            .min()
            .unwrap_or(0)
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
        self.positive_fscale > 0 && self.negative_fscale > 0
    }

    /// Map a quantized weight to a LUT index.
    #[inline]
    fn weight_to_index(w: i8) -> usize {
        (w as i16 + 127) as usize
    }

    /// Recompute the precomputed frac_max LUT from fscale_lut.
    fn recompute_frac_max_lut(&mut self) {
        for i in 0..255 {
            self.frac_max_lut[i] = if self.fscale_lut[i] == 0 {
                u64::MAX
            } else {
                u64::MAX / self.fscale_lut[i]
            };
        }
    }

    /// Compute the effective fscale for a given weight via LUT lookup.
    #[inline]
    fn effective_fscale(&self, w: i8) -> f64 {
        self.fscale_lut[Self::weight_to_index(w)] as f64
    }

    /// Soft LUT-based filter. Uses precomputed frac_max for zero floating-point overhead.
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

    /// Public wrapper for effective_fscale. Takes a quantized i8 weight.
    pub fn effective_fscale_at(&self, w: i8) -> f64 {
        self.effective_fscale(w)
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
        file.write_all(&self.positive_fscale.to_le_bytes())?;
        file.write_all(&self.negative_fscale.to_le_bytes())?;
        file.write_all(&self.lambda.to_le_bytes())?;
        file.write_all(&self.gamma.to_le_bytes())?;
        file.write_all(&[0u8; 4])?; // reserved
        file.write_all(&self.unbiased_fscale.to_le_bytes())?;

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
        let positive_fscale = u64::from_le_bytes(buf8);

        file.read_exact(&mut buf8)?;
        let negative_fscale = u64::from_le_bytes(buf8);

        file.read_exact(&mut buf4)?;
        let lambda = f32::from_le_bytes(buf4);

        file.read_exact(&mut buf4)?;
        let gamma = f32::from_le_bytes(buf4);

        let mut reserved = [0u8; 4];
        file.read_exact(&mut reserved)?;
        if reserved != [0u8; 4] {
            anyhow::bail!(
                "Reserved bytes are non-zero — file may be corrupt or unsupported format"
            );
        }

        file.read_exact(&mut buf8)?;
        let unbiased_fscale = u64::from_le_bytes(buf8);

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
            positive_fscale,
            negative_fscale,
            lambda,
            gamma,
            unbiased_fscale,
            fscale_lut,
            frac_max_lut: [0u64; 255],
        };
        table.recompute_frac_max_lut();
        Ok(table)
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let header_size = 4 + 4 + 1 + 8 + 4 + 1 + 4 + 1 + 4 + 4 + 8 + 8 + 4 + 4 + 4 + 8 + 255 * 8; // 2111 bytes
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
        out.extend_from_slice(&self.positive_fscale.to_le_bytes());
        out.extend_from_slice(&self.negative_fscale.to_le_bytes());
        out.extend_from_slice(&self.lambda.to_le_bytes());
        out.extend_from_slice(&self.gamma.to_le_bytes());
        out.extend_from_slice(&[0u8; 4]); // reserved
        out.extend_from_slice(&self.unbiased_fscale.to_le_bytes());

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
        let positive_fscale = u64::from_le_bytes(data[35..43].try_into().unwrap());
        let negative_fscale = u64::from_le_bytes(data[43..51].try_into().unwrap());
        let lambda = f32::from_le_bytes(data[51..55].try_into().unwrap());
        let gamma = f32::from_le_bytes(data[55..59].try_into().unwrap());
        let reserved: [u8; 4] = data[59..63].try_into().unwrap();
        if reserved != [0u8; 4] {
            anyhow::bail!(
                "Reserved bytes are non-zero — file may be corrupt or unsupported format"
            );
        }
        let unbiased_fscale = u64::from_le_bytes(data[63..71].try_into().unwrap());

        let mut fscale_lut = [0u64; 255];
        for i in 0..255 {
            let offset = 71 + i * 8;
            fscale_lut[i] = u64::from_le_bytes(data[offset..offset + 8].try_into().unwrap());
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
            positive_fscale,
            negative_fscale,
            lambda,
            gamma,
            unbiased_fscale,
            fscale_lut,
            frac_max_lut: [0u64; 255],
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
            eprintln!(
                "  Filter mode: optimized LUT (lambda={:.2}, gamma={:.2})",
                self.lambda, self.gamma
            );
        } else {
            eprintln!("  Filter mode: hard cutoff");
        }
        eprintln!();

        // Primary outcomes
        eprintln!("Calibration:");
        if self.is_soft_filter() {
            eprintln!("  positive_fscale:     {}", self.positive_fscale);
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
            eprintln!("  lambda:              {:.4}", self.lambda);
            eprintln!("  gamma:               {:.4}", self.gamma);
            if self.unbiased_fscale > 0 {
                eprintln!("  unbiased_fscale:     {}", self.unbiased_fscale);
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

        // Weight distribution
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

        // Collect distinct fscale values for summary
        let distinct: std::collections::BTreeSet<u64> = self.fscale_lut.iter().copied().collect();
        let lut_min = self.fscale_lut.iter().copied().min().unwrap_or(0);
        let lut_max = self.fscale_lut.iter().copied().max().unwrap_or(0);

        eprintln!();
        eprintln!(
            "LUT response curve (optimized, {} distinct fscale values):",
            distinct.len()
        );
        eprintln!("  fscale range: {} \u{2014} {}", lut_min, lut_max);

        let mut points = std::collections::BTreeSet::new();
        let start_q = (min * QUANTIZATION_SCALE).floor() as i32;
        let end_q = (max * QUANTIZATION_SCALE).ceil() as i32;
        for q in start_q..=end_q {
            points.insert(q.clamp(-127, 127));
        }
        points.insert(0);
        points.insert(threshold_q);

        eprintln!();
        eprintln!(
            "  {:>7}  {:>11}  {:>10}  {:>8}",
            "Weight", "Eff.fscale", "Retention", "vs Base"
        );
        eprintln!(
            "  {:>7}  {:>11}  {:>10}  {:>8}",
            "------", "----------", "---------", "-------"
        );

        for q in &points {
            let w = (*q).clamp(-128, 127) as i8;
            let eff = self.effective_fscale_at(w);
            let retention_pct = 100.0 / eff;
            let vs_base = base / eff;
            let weight_f = *q as f64 / 10.0;

            let mut markers: Vec<&str> = Vec::new();
            if *q == 0 {
                markers.push("unseen");
            }
            if *q == threshold_q && *q != 0 {
                markers.push("threshold (informational)");
            }

            let marker = if markers.is_empty() {
                String::new()
            } else {
                format!("  <- {}", markers.join(", "))
            };

            eprintln!(
                "  {:>7.2}  {:>11.1}  {:>9.4}%  {:>7.2}x{}",
                weight_f, eff, retention_pct, vs_base, marker
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
            && self.positive_fscale == other.positive_fscale
            && self.negative_fscale == other.negative_fscale
            && self.lambda == other.lambda
            && self.gamma == other.gamma
            && self.unbiased_fscale == other.unbiased_fscale
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

        let table =
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 0, 0, 0, None, 1.0).unwrap();
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

        let table =
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 0, 0, 0, None, 1.0).unwrap();

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

        let table =
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 0, 0, 0, None, 1.0).unwrap();

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

        let table =
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 0, 0, 0, None, 1.0).unwrap();

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

        let table =
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 0, 0, 0, None, 1.0).unwrap();
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
            target_positive_retention: None,
            target_negative_retention: None,
            positive_fscale: None,
            negative_fscale: None,
            unbiased_fscale: None,
            lambda: None,
            gamma: 1.0,
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

        // Hard cutoff mode (positive_fscale=0, negative_fscale=0)
        let table =
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 0, 0, 0, None, 1.0).unwrap();
        assert!(!table.is_soft_filter());

        // Verify passes_filter matches weight >= threshold for many hashes
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

        let table =
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 1, 1000, 0, None, 1.0)
                .unwrap();
        assert!(table.is_soft_filter());

        // Same hash always gives the same result
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

        let table =
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 10, 5000, 0, None, 1.0)
                .unwrap();
        assert!(table.is_soft_filter());

        let bytes = table.to_bytes();
        let loaded = HashBiasTable::from_bytes(&bytes).unwrap();
        assert_eq!(table, loaded);
    }

    #[test]
    fn test_unseen_midpoint_is_geometric_mean() {
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

        let table = HashBiasTable::build(
            &pos_raw, &neg_raw, 1.0, None, None, 100, 10_000, 0, None, 1.0,
        )
        .unwrap();
        assert!(table.is_soft_filter());

        let expected = (100.0_f64 * 10_000.0_f64).sqrt();
        let actual = table.effective_fscale_at(0);
        assert_approx_eq(actual, expected, 1e-6);
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

        // positive_fscale < global fscale → error
        let result =
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 50, 1000, 0, None, 1.0);
        assert!(result.is_err());

        // positive_fscale without negative_fscale → error
        let result =
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 100, 0, 0, None, 1.0);
        assert!(result.is_err());

        // negative_fscale without positive_fscale → error
        let result =
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 0, 1000, 0, None, 1.0);
        assert!(result.is_err());

        // negative_fscale <= positive_fscale → error
        let result =
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 100, 100, 0, None, 1.0);
        assert!(result.is_err());

        // Valid soft filter
        let result =
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 100, 5000, 0, None, 1.0);
        assert!(result.is_ok());

        // unbiased_fscale < positive_fscale → error
        let result = HashBiasTable::build(
            &pos_raw, &neg_raw, 1.0, None, None, 100, 5000, 50, None, 1.0,
        );
        assert!(result.is_err());

        // unbiased_fscale > negative_fscale → error
        let result = HashBiasTable::build(
            &pos_raw, &neg_raw, 1.0, None, None, 100, 5000, 6000, None, 1.0,
        );
        assert!(result.is_err());

        // unbiased_fscale without soft filter → error
        let result =
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 0, 0, 500, None, 1.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_soft_filter_monotonicity() {
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

        let table = HashBiasTable::build(
            &pos_raw, &neg_raw, 1.0, None, None, 100, 10_000, 0, None, 1.0,
        )
        .unwrap();
        assert!(table.is_soft_filter());

        let pos_fscale = table.positive_fscale as f64;
        let neg_fscale = table.negative_fscale as f64;

        // Positive side: effective_fscale should monotonically decrease toward positive_fscale
        let mut prev = table.effective_fscale_at(1);
        for w in 2..=127i8 {
            let curr = table.effective_fscale_at(w);
            assert!(
                curr <= prev + 1e-6,
                "Positive side not monotonic at w={}: prev={}, curr={}",
                w,
                prev,
                curr
            );
            prev = curr;
        }
        // At w=127, should be at or below midpoint (anchor is at geometric mean)
        let midpoint = (pos_fscale * neg_fscale).sqrt();
        let at_127 = table.effective_fscale_at(127);
        assert!(
            at_127 <= midpoint + 1e-6,
            "w=127 should be at or below midpoint={}. got={}",
            midpoint,
            at_127
        );

        // Negative side: effective_fscale should monotonically increase toward negative_fscale
        let mut prev = table.effective_fscale_at(-1);
        for w in (-127..=-2).rev() {
            let w = w as i8;
            let curr = table.effective_fscale_at(w);
            assert!(
                curr >= prev - 1e-6,
                "Negative side not monotonic at w={}: prev={}, curr={}",
                w,
                prev,
                curr
            );
            prev = curr;
        }
        // At w=-127, should be closer to negative_fscale than to midpoint
        let at_neg127 = table.effective_fscale_at(-127);
        assert!(
            at_neg127 > midpoint,
            "w=-127 should be above midpoint={}. got={}",
            midpoint,
            at_neg127
        );
    }

    #[test]
    fn test_custom_unbiased_fscale() {
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

        // Custom unbiased_fscale = 500 (between 100 and 10_000)
        let table = HashBiasTable::build(
            &pos_raw, &neg_raw, 1.0, None, None, 100, 10_000, 500, None, 1.0,
        )
        .unwrap();
        assert!(table.is_soft_filter());
        assert_eq!(table.unbiased_fscale, 500);

        // w=0 should return unbiased_fscale, not geometric mean
        assert_approx_eq(table.effective_fscale_at(0), 500.0, 1e-6);

        // Geometric mean would be sqrt(100 * 10_000) = 1000
        // Our custom value 500 is different
        let geom_mean = (100.0_f64 * 10_000.0_f64).sqrt();
        assert!((table.effective_fscale_at(0) - geom_mean).abs() > 1.0);
    }

    #[test]
    fn test_lut_direct_lookup() {
        // Build a LUT that linearly decreases fscale from 10000 to 100
        let mut fscale_lut = [0u64; 255];
        for i in 0..255 {
            // fscale decreases as weight increases (retention increases)
            fscale_lut[i] = 10_000u64 - ((10_000u64 - 100) * i as u64 / 254);
        }
        // Pin w=0 (index 127) to unbiased_fscale
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
            positive_fscale: 100,
            negative_fscale: 10_000,
            lambda: 10.0,
            gamma: 1.0,
            unbiased_fscale: 1_000,
            fscale_lut,
            frac_max_lut: [0u64; 255],
        };
        table.recompute_frac_max_lut();

        // w=0 (index 127) returns the pinned value
        assert_approx_eq(table.effective_fscale_at(0), 1_000.0, 1e-6);

        // w=-127 (index 0) should be at max fscale
        assert_approx_eq(table.effective_fscale_at(-127), 10_000.0, 1e-6);

        // w=127 (index 254) should be at min fscale
        assert_approx_eq(table.effective_fscale_at(127), 100.0, 1.0);

        // Monotonicity: higher weight → lower fscale
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
        // All fscale_lut entries should be in [positive_fscale, negative_fscale]
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

        let table = HashBiasTable::build(
            &pos_raw, &neg_raw, 1.0, None, None, 100, 10_000, 0, None, 1.0,
        )
        .unwrap();
        assert!(table.is_soft_filter());

        for (i, &fs) in table.fscale_lut.iter().enumerate() {
            assert!(
                fs >= table.positive_fscale && fs <= table.negative_fscale,
                "fscale_lut[{}] = {} out of bounds [{}, {}]",
                i,
                fs,
                table.positive_fscale,
                table.negative_fscale
            );
        }
    }

    #[test]
    fn test_gamma_zero_produces_step() {
        // gamma=0 should produce at most 3 distinct fscale values
        // (left extreme, anchor, right extreme)
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

        let table = HashBiasTable::build(
            &pos_raw, &neg_raw, 1.0, None, None, 100, 10_000, 0, None, 0.0,
        )
        .unwrap();
        assert!(table.is_soft_filter());

        let distinct: std::collections::BTreeSet<u64> = table.fscale_lut.iter().copied().collect();
        assert!(
            distinct.len() <= 3,
            "gamma=0 should produce at most 3 distinct fscale values, got {}",
            distinct.len()
        );
    }

    #[test]
    fn test_gamma_positive_produces_graded() {
        // With synthetic histogram data, gamma>0 should produce a graded curve
        let mut p_pos = [0.0f64; 255];
        let mut p_neg = [0.0f64; 255];
        // Positive hashes cluster at high weights, negative at low
        for i in 128..255 {
            p_pos[i] = 1.0 / 127.0;
        }
        for i in 0..127 {
            p_neg[i] = 1.0 / 127.0;
        }

        let r_min = 1.0 / 10_000.0;
        let r_max = 1.0 / 100.0;
        let anchor_val = 1.0 / 1_000.0;

        let lut = optimize_lut(
            &p_pos,
            &p_neg,
            10.0,
            5.0,
            r_min,
            r_max,
            Some(127),
            anchor_val,
        );

        let distinct: std::collections::BTreeSet<u64> = lut.iter().copied().collect();
        assert!(
            distinct.len() > 3,
            "gamma>0 should produce more than 3 distinct fscale values, got {}",
            distinct.len()
        );
    }

    #[test]
    fn test_pava_basic() {
        // Simple test: enforce monotonically non-decreasing
        let mut r = [3.0, 1.0, 2.0, 5.0, 4.0];
        pava_increasing(&mut r, 0.0, 10.0);
        for i in 0..r.len() - 1 {
            assert!(
                r[i] <= r[i + 1] + 1e-10,
                "PAVA failed monotonicity at {}: {} > {}",
                i,
                r[i],
                r[i + 1]
            );
        }

        // All values should be clamped to bounds
        let mut r2 = [5.0, 3.0, 1.0];
        pava_increasing(&mut r2, 2.0, 4.0);
        for &v in &r2 {
            assert!(v >= 2.0 - 1e-10 && v <= 4.0 + 1e-10);
        }
        for i in 0..r2.len() - 1 {
            assert!(r2[i] <= r2[i + 1] + 1e-10);
        }
    }

    #[test]
    fn test_optimizer_convergence() {
        // Create a simple histogram and verify optimizer produces valid output
        let mut p_pos = [0.0f64; 255];
        let mut p_neg = [0.0f64; 255];

        // Positive hashes cluster at high weights, negative at low weights
        for i in 128..255 {
            p_pos[i] = 1.0 / 127.0;
        }
        for i in 0..127 {
            p_neg[i] = 1.0 / 127.0;
        }

        let r_min = 1.0 / 10_000.0;
        let r_max = 1.0 / 100.0;
        let anchor_val = 1.0 / 1_000.0;

        let lut = optimize_lut(
            &p_pos,
            &p_neg,
            10.0,
            1.0,
            r_min,
            r_max,
            Some(127),
            anchor_val,
        );

        // Verify bounds
        for (i, &fs) in lut.iter().enumerate() {
            assert!(
                fs >= 100 && fs <= 10_000,
                "lut[{}] = {} out of bounds [100, 10000]",
                i,
                fs
            );
        }

        // Verify monotonicity (fscale should decrease as index increases)
        for i in 0..254 {
            assert!(
                lut[i] >= lut[i + 1],
                "LUT not monotonic at {}: {} < {}",
                i,
                lut[i],
                lut[i + 1]
            );
        }

        // Verify anchor
        assert_eq!(
            lut[127], 1000,
            "Anchor at index 127 should be 1000, got {}",
            lut[127]
        );
    }
}
