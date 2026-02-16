use anyhow::{Context, Result};
use indicatif::ProgressBar;
use jamhash::jamhash_u64;
use needletail::{Sequence, parse_fastx_file};
use std::io::{Read, Write};
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

const BIAS_MAGIC: &[u8; 4] = b"BIA4";
const BIAS_VERSION: u32 = 4;

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
    pub curve_power: f32,
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
    pub k_pos: f32,
    pub k_neg: f32,
    pub curve_power: f32,
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
            config.curve_power,
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
        curve_power: f32,
    ) -> Result<Self> {
        validate_cms_compatibility(positive, negative)?;

        if !alpha.is_finite() || alpha <= 0.0 {
            anyhow::bail!("alpha must be finite and > 0, got {}", alpha);
        }

        if !curve_power.is_finite() || curve_power <= 0.0 {
            anyhow::bail!("curve_power must be finite and > 0, got {}", curve_power);
        }

        if let Some(target) = target_positive_retention {
            if !target.is_finite() || !(0.0..=1.0).contains(&target) {
                anyhow::bail!("target_positive_retention must be in [0,1], got {}", target);
            }
        }

        if let Some(target) = target_negative_retention {
            if !target.is_finite() || !(0.0..=1.0).contains(&target) {
                anyhow::bail!("target_negative_retention must be in [0,1], got {}", target);
            }
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
            k_pos: 0.0,
            k_neg: 0.0,
            curve_power,
        };

        // Auto-derive k_pos and k_neg when soft filter is active
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

            // Derive k_pos
            table.k_pos = if let Some(target) = target_positive_retention {
                Self::search_k_for_target_retention(&table, &pos_in_range, target, true)
            } else {
                Self::sweep_k_for_max_enrichment(&table, &pos_in_range, &neg_in_range, true)
            };

            // Derive k_neg (with k_pos already set)
            table.k_neg = if let Some(target) = target_negative_retention {
                Self::search_k_for_target_retention(&table, &neg_in_range, target, false)
            } else {
                Self::sweep_k_for_max_enrichment(&table, &pos_in_range, &neg_in_range, false)
            };

            // Final retention computation with both k values set
            let (pos_ret, neg_ret) = Self::compute_retention(&table, &pos_in_range, &neg_in_range);
            table.positive_retention = pos_ret;
            table.negative_retention = neg_ret;
            if table.negative_retention > 0.0 {
                table.max_fold_enrichment = table.positive_retention / table.negative_retention;
            } else {
                table.max_fold_enrichment = f32::INFINITY;
            }
        }

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

    /// Binary search for k that achieves the target retention.
    /// `for_positive_side`: true = searching k_pos (affects w>=0), false = searching k_neg (affects w<0).
    fn search_k_for_target_retention(
        table: &Self,
        hashes: &[u64],
        target: f32,
        for_positive_side: bool,
    ) -> f32 {
        if hashes.is_empty() {
            return 0.0;
        }

        let mut lo: f32 = 0.001;
        let mut hi: f32 = 5.0;

        for _ in 0..50 {
            let mid = (lo + hi) / 2.0;
            let mut probe = table.clone();
            if for_positive_side {
                probe.k_pos = mid;
            } else {
                probe.k_neg = mid;
            }

            let retention = hashes
                .iter()
                .filter(|&&h| probe.passes_soft_filter(h))
                .count() as f32
                / hashes.len() as f32;

            if for_positive_side {
                // Higher k_pos → higher positive retention
                if retention < target {
                    lo = mid;
                } else {
                    hi = mid;
                }
            } else {
                // Higher k_neg → lower negative retention
                if retention > target {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
        }

        (lo + hi) / 2.0
    }

    /// Sweep k over a log-scale grid to maximize fold enrichment.
    /// `for_positive_side`: true = sweeping k_pos, false = sweeping k_neg.
    fn sweep_k_for_max_enrichment(
        table: &Self,
        pos_in_range: &[u64],
        neg_in_range: &[u64],
        for_positive_side: bool,
    ) -> f32 {
        const STEPS: usize = 200;
        const MIN_RETENTION: f64 = 1e-6;
        let log_lo = 0.001_f64.ln();
        let log_hi = 5.0_f64.ln();

        let mut best_k: f32 = 0.001;
        let mut best_enrichment: f64 = 0.0;

        for i in 0..=STEPS {
            let k = (log_lo + (log_hi - log_lo) * i as f64 / STEPS as f64).exp() as f32;
            let mut probe = table.clone();
            if for_positive_side {
                probe.k_pos = k;
            } else {
                probe.k_neg = k;
            }

            let (pos_ret, neg_ret) = Self::compute_retention(&probe, pos_in_range, neg_in_range);

            let neg_ret = (neg_ret as f64).max(MIN_RETENTION);
            let enrichment = pos_ret as f64 / neg_ret;

            if enrichment > best_enrichment && enrichment.is_finite() {
                best_enrichment = enrichment;
                best_k = k;
            }
        }

        best_k
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

    /// Compute the effective fscale using piecewise sigmoid centered at 0.
    /// w >= 0 uses k_pos, w < 0 uses k_neg. At w=0: s=0.5 → geometric mean (midpoint anchor).
    /// When curve_power > 1.0, a power transform flattens the sigmoid's middle region,
    /// so only strongly biased k-mers get preferential/penalized treatment.
    #[inline]
    fn effective_fscale(&self, w: i8) -> f64 {
        let k = if w >= 0 { self.k_pos } else { self.k_neg };
        let w_norm = w as f64 / 127.0;
        let w_curved = w_norm.signum() * w_norm.abs().powf(self.curve_power as f64);
        let x = k as f64 * w_curved * 127.0;
        let s = 1.0 / (1.0 + (-x).exp());
        let log_pos = (self.positive_fscale as f64).ln();
        let log_neg = (self.negative_fscale as f64).ln();
        (log_pos * s + log_neg * (1.0 - s)).exp()
    }

    /// Soft sigmoid-based filter. Uses the hash value directly against an adjusted frac_max.
    #[inline]
    fn passes_soft_filter(&self, hash: u64) -> bool {
        let w = self.weight(hash);
        let eff = self.effective_fscale(w);
        hash < (u64::MAX as f64 / eff) as u64
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
        file.write_all(&self.k_pos.to_le_bytes())?;
        file.write_all(&self.k_neg.to_le_bytes())?;
        file.write_all(&self.curve_power.to_le_bytes())?;

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
                "Invalid bias table file (bad magic, expected BIA4): {}",
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
        let k_pos = f32::from_le_bytes(buf4);

        file.read_exact(&mut buf4)?;
        let k_neg = f32::from_le_bytes(buf4);

        file.read_exact(&mut buf4)?;
        let curve_power = f32::from_le_bytes(buf4);

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

        Ok(Self {
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
            k_pos,
            k_neg,
            curve_power,
        })
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let header_size = 4 + 4 + 1 + 8 + 4 + 1 + 4 + 1 + 4 + 4 + 8 + 8 + 4 + 4 + 4; // 63 bytes
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
        out.extend_from_slice(&self.k_pos.to_le_bytes());
        out.extend_from_slice(&self.k_neg.to_le_bytes());
        out.extend_from_slice(&self.curve_power.to_le_bytes());

        for &seed in &self.seeds {
            out.extend_from_slice(&seed.to_le_bytes());
        }
        for &w in &self.weights {
            out.push(w as u8);
        }

        out
    }

    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        if data.len() < 63 {
            anyhow::bail!(
                "Bias table data too small: {} bytes (expected >= 63)",
                data.len()
            );
        }

        let magic: [u8; 4] = data[0..4].try_into().unwrap();
        if &magic != BIAS_MAGIC {
            anyhow::bail!("Invalid bias table magic bytes (expected BIA4)");
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
        let k_pos = f32::from_le_bytes(data[51..55].try_into().unwrap());
        let k_neg = f32::from_le_bytes(data[55..59].try_into().unwrap());
        let curve_power = f32::from_le_bytes(data[59..63].try_into().unwrap());

        let seeds_start = 63;
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

        Ok(Self {
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
            k_pos,
            k_neg,
            curve_power,
        })
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
        eprintln!("Hash Bias Table (v4)");
        eprintln!("====================");
        eprintln!("  k-mer size:     {}", self.config.k);
        eprintln!("  fscale:         {}", self.config.fscale);
        eprintln!(
            "  CMS dimensions: {} x {}",
            self.config.width, self.config.depth
        );
        eprintln!("  Smoothing (alpha): {:.1}", self.alpha);
        if self.is_soft_filter() {
            eprintln!("  Filter mode: soft sigmoid (piecewise, centered at 0)");
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
            eprintln!("  k_pos:               {:.4}", self.k_pos);
            eprintln!("  k_neg:               {:.4}", self.k_neg);
            eprintln!("  curve_power:         {:.2}", self.curve_power);

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
            self.print_sigmoid_curve(min, max);
        }
    }

    pub fn print_sigmoid_curve(&self, min: f32, max: f32) {
        let base = self.config.fscale as f64;
        let k_pos = self.k_pos as f64;
        let k_neg = self.k_neg as f64;
        let threshold_q = self.threshold as i32;
        let ln19 = 19.0_f64.ln();

        // Asymmetric 5%/95% bounds centered at 0
        let neg_5pct_qi = if k_neg > 0.0 {
            (-ln19 / k_neg).round() as i32
        } else {
            -127
        };
        let pos_95pct_qi = if k_pos > 0.0 {
            (ln19 / k_pos).round() as i32
        } else {
            127
        };

        eprintln!();
        eprintln!("Sigmoid response curve (piecewise, centered at 0):");
        eprintln!(
            "  Bounds: {:.2} (neg 5%) \u{2014} 0.00 (unseen/midpoint) \u{2014} {:.2} (pos 95%)",
            neg_5pct_qi as f64 / 10.0,
            pos_95pct_qi as f64 / 10.0
        );

        let mut points = std::collections::BTreeSet::new();
        let start_q = (min.floor() as i32) * 10;
        let end_q = (max.ceil() as i32) * 10;
        for q in (start_q..=end_q).step_by(10) {
            points.insert(q.clamp(-127, 127));
        }
        points.insert(0);
        points.insert(threshold_q);
        points.insert(neg_5pct_qi.clamp(-127, 127));
        points.insert(pos_95pct_qi.clamp(-127, 127));

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
                markers.push("threshold");
            }
            if *q == neg_5pct_qi && *q != 0 {
                markers.push("neg 5%");
            }
            if *q == pos_95pct_qi && *q != 0 {
                markers.push("pos 95%");
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
    63 + DEFAULT_CMS_DEPTH * 8 + DEFAULT_CMS_WIDTH * DEFAULT_CMS_DEPTH;

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
            && self.k_pos == other.k_pos
            && self.k_neg == other.k_neg
            && self.curve_power == other.curve_power
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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 0, 0, 1.0).unwrap();
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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 0, 0, 1.0).unwrap();

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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 0, 0, 1.0).unwrap();

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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 0, 0, 1.0).unwrap();

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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 0, 0, 1.0).unwrap();
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
            curve_power: 1.0,
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
        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 0, 0, 1.0).unwrap();
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
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 1, 1000, 1.0).unwrap();
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
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 10, 5000, 1.0).unwrap();
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

        let table =
            HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 100, 10_000, 1.0).unwrap();
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
        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 50, 1000, 1.0);
        assert!(result.is_err());

        // positive_fscale without negative_fscale → error
        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 100, 0, 1.0);
        assert!(result.is_err());

        // negative_fscale without positive_fscale → error
        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 0, 1000, 1.0);
        assert!(result.is_err());

        // negative_fscale <= positive_fscale → error
        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 100, 100, 1.0);
        assert!(result.is_err());

        // Valid soft filter
        let result = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None, None, 100, 5000, 1.0);
        assert!(result.is_ok());
    }
}
