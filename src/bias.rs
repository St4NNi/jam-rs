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
const MIN_BIN_FRACTION: f64 = 0.01;

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
    pub positive_fscale: Option<u64>,
    pub negative_fscale: Option<u64>,
    pub unbiased_fscale: Option<u64>,
    pub min_positive_retention: f32,
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
    pub min_positive_retention: f32,
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

fn fill_unreliable_bins(
    lut: &mut [u64; 255],
    pos_counts: &[u64; 255],
    neg_counts: &[u64; 255],
    w0_idx: usize,
) {
    let total: u64 = pos_counts.iter().sum::<u64>() + neg_counts.iter().sum::<u64>();
    let expected_per_bin = total as f64 / 255.0;
    let threshold = (expected_per_bin * MIN_BIN_FRACTION).max(1.0) as u64;

    let is_reliable =
        |i: usize| -> bool { i == w0_idx || pos_counts[i] + neg_counts[i] >= threshold };

    let mut last_reliable_fscale = lut[w0_idx];
    for i in (0..w0_idx).rev() {
        if is_reliable(i) {
            last_reliable_fscale = lut[i];
        } else {
            lut[i] = last_reliable_fscale;
        }
    }

    last_reliable_fscale = lut[w0_idx];
    for i in (w0_idx + 1)..255 {
        if is_reliable(i) {
            last_reliable_fscale = lut[i];
        } else {
            lut[i] = last_reliable_fscale;
        }
    }
}

fn enrichment_lut(
    p_pos: &[f64; 255],
    p_neg: &[f64; 255],
    pos_counts: &[u64; 255],
    neg_counts: &[u64; 255],
    base_fscale: u64,
    positive_fscale: u64,
    negative_fscale: u64,
    unbiased_fscale: u64,
    min_positive_retention: f64,
) -> Result<[u64; 255]> {
    if base_fscale == 0 || positive_fscale == 0 || negative_fscale == 0 {
        anyhow::bail!("fscale values must be > 0");
    }

    if !min_positive_retention.is_finite() || !(0.0..=1.0).contains(&min_positive_retention) {
        anyhow::bail!(
            "min_positive_retention must be in [0,1], got {}",
            min_positive_retention
        );
    }

    let base_fs = base_fscale as f64;
    let pos_fs = positive_fscale as f64;
    let neg_fs = negative_fscale as f64;
    let unbiased_fs = unbiased_fscale as f64;
    let w0_idx = HashBiasTable::weight_to_index(0);
    let retention_tol = 1e-6;

    let mut e = [0.0f64; 255];
    for i in 0..255 {
        e[i] = p_pos[i] / p_neg[i].max(1e-30);
    }

    let e_min = e.iter().copied().fold(f64::INFINITY, f64::min);
    let e_max = e.iter().copied().fold(0.0f64, f64::max);
    let c_lo = (pos_fs * e_min).max(1e-30).ln();
    let c_hi = (neg_fs * e_max).max(1e-30).ln();

    let retention_for_fscale = |fs: f64| -> f64 { (base_fs / fs).min(1.0) };

    let max_pos_ret = (0..255)
        .map(|i| {
            let fs = if i == w0_idx && unbiased_fscale > 0 {
                unbiased_fs
            } else {
                pos_fs
            };
            p_pos[i] * retention_for_fscale(fs)
        })
        .sum::<f64>();

    if max_pos_ret + retention_tol < min_positive_retention {
        anyhow::bail!(
            "min_positive_retention={} is infeasible (max achievable {:.6}) for fscale range [{}, {}]",
            min_positive_retention,
            max_pos_ret,
            positive_fscale,
            negative_fscale
        );
    }

    let compute_score = |c: f64| -> (f64, f64, f64) {
        let mut pos_ret = 0.0f64;
        let mut neg_ret = 0.0f64;
        for i in 0..255 {
            let fs = if i == w0_idx && unbiased_fscale > 0 {
                unbiased_fs
            } else {
                (c / e[i]).clamp(pos_fs, neg_fs)
            };
            let r = retention_for_fscale(fs);
            pos_ret += p_pos[i] * r;
            neg_ret += p_neg[i] * r;
        }
        let score = if neg_ret > 0.0 {
            pos_ret * pos_ret / neg_ret
        } else {
            f64::INFINITY
        };
        (score, pos_ret, neg_ret)
    };

    // Coarse grid search
    let grid_steps = 500;
    let mut best_c = c_lo.exp();
    let mut best_score = f64::NEG_INFINITY;
    let mut found_feasible = false;

    for step in 0..=grid_steps {
        let c = (c_lo + (c_hi - c_lo) * step as f64 / grid_steps as f64).exp();
        let (score, pos_ret, _) = compute_score(c);
        if pos_ret + retention_tol >= min_positive_retention && score > best_score {
            found_feasible = true;
            best_score = score;
            best_c = c;
        }
    }

    if !found_feasible {
        anyhow::bail!(
            "no feasible LUT found for min_positive_retention={} in fscale range [{}, {}]",
            min_positive_retention,
            positive_fscale,
            negative_fscale
        );
    }

    // Ternary search refinement
    let grid_width = (c_hi - c_lo) / grid_steps as f64;
    let mut lo = (best_c.ln() - grid_width).max(c_lo);
    let mut hi = (best_c.ln() + grid_width).min(c_hi);
    for _ in 0..50 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        let (s1, pr1, _) = compute_score(m1.exp());
        let (s2, pr2, _) = compute_score(m2.exp());
        let s1_valid = if pr1 + retention_tol >= min_positive_retention {
            s1
        } else {
            f64::NEG_INFINITY
        };
        let s2_valid = if pr2 + retention_tol >= min_positive_retention {
            s2
        } else {
            f64::NEG_INFINITY
        };
        if s1_valid < s2_valid {
            lo = m1;
            if s2_valid > best_score {
                best_score = s2_valid;
                best_c = m2.exp();
            }
        } else {
            hi = m2;
            if s1_valid > best_score {
                best_score = s1_valid;
                best_c = m1.exp();
            }
        }
    }

    let mut lut = [0u64; 255];
    for (i, slot) in lut.iter_mut().enumerate() {
        *slot = if i == w0_idx && unbiased_fscale > 0 {
            unbiased_fscale
        } else {
            (best_c / e[i]).round().clamp(pos_fs, neg_fs) as u64
        };
    }

    fill_unreliable_bins(&mut lut, pos_counts, neg_counts, w0_idx);

    Ok(lut)
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
            config.positive_fscale.unwrap_or(0),
            config.negative_fscale.unwrap_or(0),
            config.unbiased_fscale.unwrap_or(0),
            config.min_positive_retention,
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
        positive_fscale: u64,
        negative_fscale: u64,
        unbiased_fscale: u64,
        min_positive_retention: f32,
    ) -> Result<Self> {
        validate_cms_compatibility(positive, negative)?;

        if !alpha.is_finite() || alpha <= 0.0 {
            anyhow::bail!("alpha must be finite and > 0, got {}", alpha);
        }

        if !min_positive_retention.is_finite() || !(0.0..=1.0).contains(&min_positive_retention) {
            anyhow::bail!(
                "min_positive_retention must be in [0,1], got {}",
                min_positive_retention
            );
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

        let calibration = calibrate_threshold(
            positive, negative, &weights, &seeds, width,
            None,
        )?;

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
            min_positive_retention,
            unbiased_fscale,
            fscale_lut,
            frac_max_lut,
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

            table.fscale_lut = enrichment_lut(
                &p_pos,
                &p_neg,
                &pos_bin_counts,
                &neg_bin_counts,
                table.config.fscale,
                positive_fscale,
                negative_fscale,
                unbiased_fscale,
                min_positive_retention as f64,
            )?;

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
        file.write_all(&self.min_positive_retention.to_le_bytes())?;
        file.write_all(&[0u8; 8])?; // reserved
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
        let min_positive_retention = f32::from_le_bytes(buf4);

        let mut reserved = [0u8; 8];
        file.read_exact(&mut reserved)?;

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
            min_positive_retention,
            unbiased_fscale,
            fscale_lut,
            frac_max_lut: [0u64; 255],
        };
        table.recompute_frac_max_lut();
        Ok(table)
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let header_size = 4 + 4 + 1 + 8 + 4 + 1 + 4 + 1 + 4 + 4 + 8 + 8 + 4 + 8 + 8 + 255 * 8; // 2111 bytes
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
        out.extend_from_slice(&self.min_positive_retention.to_le_bytes());
        out.extend_from_slice(&[0u8; 8]); // reserved
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
        let min_positive_retention = f32::from_le_bytes(data[51..55].try_into().unwrap());
        // bytes 55..63 reserved
        let unbiased_fscale = u64::from_le_bytes(data[63..71].try_into().unwrap());

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
            positive_fscale,
            negative_fscale,
            min_positive_retention,
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
            eprintln!("  Filter mode: enrichment LUT");
        } else {
            eprintln!("  Filter mode: hard cutoff");
        }
        eprintln!();

        eprintln!("Calibration:");
        if self.is_soft_filter() {
            eprintln!("  positive_fscale:     {}", self.positive_fscale);
            eprintln!("  negative_fscale:     {}", self.negative_fscale_label());
            eprintln!("  min_pos_retention:   {:.2}", self.min_positive_retention);
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
            eprintln!(
                "  eff. fscale (pos):   {:.0}",
                self.effective_fscale_on_population(self.positive_retention)
            );
            eprintln!(
                "  eff. fscale (neg):   {:.0}",
                self.effective_fscale_on_population(self.negative_retention)
            );
            eprintln!(
                "  eff. fscale (combined): {:.0}",
                self.effective_fscale_combined()
            );
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
            && self.positive_fscale == other.positive_fscale
            && self.negative_fscale == other.negative_fscale
            && self.min_positive_retention == other.min_positive_retention
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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 0, 0, 0.1).unwrap();
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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 0, 0, 0.1).unwrap();

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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 0, 0, 0.1).unwrap();

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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 0, 0, 0.1).unwrap();

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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 0, 0, 0.1).unwrap();
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
            positive_fscale: None,
            negative_fscale: None,
            unbiased_fscale: None,
            min_positive_retention: 0.1,
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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 0, 0, 0.1).unwrap();
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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 1, 1000, 0, 0.1).unwrap();
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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 10, 5000, 0, 0.1).unwrap();
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

        // Without unbiased_fscale, w=0 is determined by enrichment data
        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 100, 10_000, 0, 0.1).unwrap();
        assert!(table.is_soft_filter());

        let at_0 = table.effective_fscale_at(0);
        assert!((100.0..=10_000.0).contains(&at_0));
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
        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 50, 1000, 0, 0.1);
        assert!(result.is_err());

        // positive_fscale without negative_fscale → error
        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 100, 0, 0, 0.1);
        assert!(result.is_err());

        // negative_fscale without positive_fscale → error
        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 1000, 0, 0.1);
        assert!(result.is_err());

        // negative_fscale <= positive_fscale → error
        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 100, 100, 0, 0.1);
        assert!(result.is_err());

        // Valid soft filter
        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 100, 5000, 0, 0.1);
        assert!(result.is_ok());

        // unbiased_fscale < positive_fscale → error
        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 100, 5000, 50, 0.1);
        assert!(result.is_err());

        // unbiased_fscale > negative_fscale → error
        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 100, 5000, 6000, 0.1);
        assert!(result.is_err());

        // unbiased_fscale without soft filter → error
        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 0, 0, 500, 0.1);
        assert!(result.is_err());
    }

    #[test]
    fn test_enrichment_lut_bounds() {
        let mut p_pos = [HISTOGRAM_EPS; 255];
        let mut p_neg = [HISTOGRAM_EPS; 255];
        // Clear separation: positive at high indices, negative at low
        for v in p_pos[180..255].iter_mut() {
            *v = 1.0 / 75.0;
        }
        for v in p_neg[0..75].iter_mut() {
            *v = 1.0 / 75.0;
        }
        // All bins reliable (high counts)
        let counts_high = [1000u64; 255];

        let lut = enrichment_lut(&p_pos, &p_neg, &counts_high, &counts_high, 100, 100, 10_000, 0, 0.0).unwrap();

        for (i, &fs) in lut.iter().enumerate() {
            assert!(
                (100..=10_000).contains(&fs),
                "lut[{}] = {} out of bounds [100, 10000]",
                i,
                fs
            );
        }

        // High-enrichment and low-enrichment buckets should differ
        let has_low_fscale = lut.iter().any(|&fs| fs <= 500);
        let has_high_fscale = lut.iter().any(|&fs| fs >= 5000);
        assert!(
            has_low_fscale && has_high_fscale,
            "LUT should have both low and high fscale values"
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
        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 100, 10_000, 500, 0.05).unwrap();
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
        for (i, slot) in fscale_lut.iter_mut().enumerate() {
            // fscale decreases as weight increases (retention increases)
            *slot = 10_000u64 - ((10_000u64 - 100) * i as u64 / 254);
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
            min_positive_retention: 0.1,
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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 100, 10_000, 0, 0.1).unwrap();
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
    fn test_enrichment_lut_separates_signal() {
        let mut p_pos = [HISTOGRAM_EPS; 255];
        let mut p_neg = [HISTOGRAM_EPS; 255];
        for v in p_pos[200..255].iter_mut() {
            *v = 1.0 / 55.0;
        }
        for v in p_neg[0..55].iter_mut() {
            *v = 1.0 / 55.0;
        }
        let counts_high = [1000u64; 255];

        let lut = enrichment_lut(&p_pos, &p_neg, &counts_high, &counts_high, 100, 100, 10_000, 0, 0.0).unwrap();

        // High-enrichment buckets (200+) should have low fscale
        let high_avg: f64 = lut[200..255].iter().map(|&f| f as f64).sum::<f64>() / 55.0;
        // Low-enrichment buckets (0..55) should have high fscale
        let low_avg: f64 = lut[0..55].iter().map(|&f| f as f64).sum::<f64>() / 55.0;

        assert!(
            low_avg > high_avg,
            "negative-enriched buckets should have higher fscale: low_avg={}, high_avg={}",
            low_avg,
            high_avg
        );
    }

    #[test]
    fn test_enrichment_lut_min_retention_floor() {
        let mut p_pos = [HISTOGRAM_EPS; 255];
        let mut p_neg = [HISTOGRAM_EPS; 255];
        // Very few positive hashes
        p_pos[200] = 1.0;
        // Lots of negative everywhere
        for v in p_neg.iter_mut() {
            *v = 1.0 / 255.0;
        }

        let counts_high = [1000u64; 255];
        let lut_low = enrichment_lut(&p_pos, &p_neg, &counts_high, &counts_high, 100, 100, 10_000, 0, 0.0).unwrap();
        let lut_high = enrichment_lut(&p_pos, &p_neg, &counts_high, &counts_high, 100, 100, 10_000, 0, 0.5).unwrap();

        // With high retention floor, more buckets should have lower fscale (higher retention)
        let open_low = lut_low.iter().filter(|&&f| f == 100).count();
        let open_high = lut_high.iter().filter(|&&f| f == 100).count();
        assert!(
            open_high >= open_low,
            "higher min retention should open more buckets: floor=0.0 has {} open, floor=0.5 has {} open",
            open_low,
            open_high
        );
    }

    #[test]
    fn test_enrichment_lut_uses_full_range() {
        let mut p_pos = [HISTOGRAM_EPS; 255];
        let mut p_neg = [HISTOGRAM_EPS; 255];
        // Positive hashes at high indices, negative at low, with overlap
        p_pos[80..255].fill(1.0 / 175.0);
        p_neg[0..180].fill(1.0 / 180.0);
        let counts_high = [1000u64; 255];

        let lut = enrichment_lut(&p_pos, &p_neg, &counts_high, &counts_high, 100, 100, 10_000, 0, 0.0).unwrap();

        let has_low = lut.iter().any(|&fs| fs <= 500);
        let has_high = lut.iter().any(|&fs| fs >= 3000);
        assert!(
            has_low && has_high,
            "LUT should use the fscale range: has_low={}, has_high={}",
            has_low,
            has_high
        );
    }
}
