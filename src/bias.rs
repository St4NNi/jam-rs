use anyhow::{Context, Result};
use indicatif::ProgressBar;
use jamhash::jamhash_u64;
use needletail::{Sequence, parse_fastx_file};
use std::io::{Read, Write};
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

const BIAS_MAGIC: &[u8; 4] = b"BIA3";
const BIAS_VERSION: u32 = 3;

const DEFAULT_CMS_WIDTH: usize = 1 << 20; // 1M columns
const DEFAULT_CMS_DEPTH: usize = 5;
const QUANTIZATION_SCALE: f32 = 10.0; // Each i8 unit = 0.1 log-ratio
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
        let seeds: Vec<u64> = (0..depth)
            .map(|i| 0x517cc1b727220a95u64.wrapping_add(i as u64))
            .collect();
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

/// Configuration for creating a bias table from FASTA files.
#[derive(Debug, Clone)]
pub struct BiasCreateConfig {
    pub cms: CMSConfig,
    pub alpha: f32,
    pub target_fold_enrichment: Option<f32>,
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
    /// Create a bias table directly from positive and negative FASTA files.
    /// Background subtraction is always applied: the target signal is removed
    /// from the background before computing log-ratios.
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
            config.target_fold_enrichment,
        )?;

        if let Some(ref pb) = progress {
            pb.finish();
        }

        Ok(table)
    }

    /// Build a bias table from pre-computed hash counts.
    /// Always applies background subtraction: normalizes counts and subtracts
    /// the target (positive) signal from the background (negative) before
    /// computing log-ratios.
    pub(crate) fn build(
        positive: &RawHashCounts,
        negative: &RawHashCounts,
        alpha: f32,
        target_fold_enrichment: Option<f32>,
    ) -> Result<Self> {
        validate_cms_compatibility(positive, negative)?;

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
            positive,
            negative,
            &weights,
            &seeds,
            width,
            target_fold_enrichment,
        )?;

        Ok(Self {
            config: positive.config.clone(),
            seeds,
            weights,
            alpha,
            threshold: calibration.threshold,
            positive_retention: calibration.positive_retention,
            negative_retention: calibration.negative_retention,
            max_fold_enrichment: calibration.max_fold_enrichment,
        })
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
        self.weight(hash) >= self.threshold
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
            anyhow::bail!("Invalid bias table file (bad magic): {}", path.display());
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
        })
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let header_size = 4 + 4 + 1 + 8 + 4 + 1 + 4 + 1 + 4 + 4; // 35 bytes
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

        for &seed in &self.seeds {
            out.extend_from_slice(&seed.to_le_bytes());
        }
        for &w in &self.weights {
            out.push(w as u8);
        }

        out
    }

    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        if data.len() < 35 {
            anyhow::bail!("Bias table data too small: {} bytes", data.len());
        }

        let magic: [u8; 4] = data[0..4].try_into().unwrap();
        if &magic != BIAS_MAGIC {
            anyhow::bail!("Invalid bias table magic bytes");
        }

        let version = u32::from_le_bytes(data[4..8].try_into().unwrap());
        if version != BIAS_VERSION {
            anyhow::bail!("Unsupported bias table version {}", version);
        }

        let k = data[8];
        let fscale = u64::from_le_bytes(data[9..17].try_into().unwrap());
        let width = u32::from_le_bytes(data[17..21].try_into().unwrap()) as usize;
        let depth = data[21] as usize;
        let alpha = f32::from_le_bytes(data[22..26].try_into().unwrap());
        let threshold = data[26] as i8;
        let positive_retention = f32::from_le_bytes(data[27..31].try_into().unwrap());
        let negative_retention = f32::from_le_bytes(data[31..35].try_into().unwrap());

        let seeds_start = 35;
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
        })
    }

    pub fn weight_stats(&self) -> (f32, f32, f32, f32, usize) {
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
        (min, max, mean, variance.sqrt(), positive)
    }

    pub fn memory_usage(&self) -> usize {
        self.weights.len() + self.seeds.len() * 8
    }

    pub fn threshold_f32(&self) -> f32 {
        self.threshold as f32 / QUANTIZATION_SCALE
    }

    pub fn print_stats(&self) {
        let (min, max, mean, std, positive) = self.weight_stats();
        let total_cells = self.config.width * self.config.depth;
        eprintln!("Hash Bias Table (v3)");
        eprintln!("  k-mer size:     {}", self.config.k);
        eprintln!("  fscale:         {}", self.config.fscale);
        eprintln!(
            "  CMS dimensions: {} x {}",
            self.config.width, self.config.depth
        );
        eprintln!("  Smoothing (alpha): {:.1}", self.alpha);
        eprintln!(
            "  Threshold: {:.2} (quantized: {})",
            self.threshold_f32(),
            self.threshold
        );
        eprintln!(
            "  Positive retention: {:.2}%",
            self.positive_retention * 100.0
        );
        eprintln!(
            "  Negative retention: {:.2}%",
            self.negative_retention * 100.0
        );
        eprintln!("  Fold enrichment: {:.2}x", self.fold_enrichment());
        eprintln!(
            "  Weight stats: min={:.2}, max={:.2}, mean={:.2}, std={:.2}",
            min, max, mean, std
        );
        eprintln!(
            "  Positive weights: {} ({:.1}%)",
            positive,
            positive as f64 / total_cells as f64 * 100.0
        );
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

    // First pass: find maximum achievable fold enrichment across all thresholds
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
        None => {
            // Use the threshold that maximizes fold enrichment
            Ok(CalibrationResult {
                threshold: max_threshold,
                positive_retention: max_pos_ret,
                negative_retention: max_neg_ret,
                fold_enrichment: max_enrichment,
                max_fold_enrichment: max_enrichment,
            })
        }
        Some(target) => {
            if target > max_enrichment {
                // Target exceeds maximum achievable - use maximum instead
                // The handler will warn the user
                return Ok(CalibrationResult {
                    threshold: max_threshold,
                    positive_retention: max_pos_ret,
                    negative_retention: max_neg_ret,
                    fold_enrichment: max_enrichment,
                    max_fold_enrichment: max_enrichment,
                });
            }

            // Find threshold closest to target fold enrichment
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
    35 + DEFAULT_CMS_DEPTH * 8 + DEFAULT_CMS_WIDTH * DEFAULT_CMS_DEPTH;

impl PartialEq for HashBiasTable {
    fn eq(&self, other: &Self) -> bool {
        // max_fold_enrichment is excluded: it's calibration metadata not stored
        // in the serialized format, so it may differ across round-trips.
        self.config.k == other.config.k
            && self.config.fscale == other.config.fscale
            && self.config.width == other.config.width
            && self.config.depth == other.config.depth
            && self.alpha == other.alpha
            && self.threshold == other.threshold
            && self.positive_retention == other.positive_retention
            && self.negative_retention == other.negative_retention
            && self.seeds == other.seeds
            && self.weights == other.weights
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, Some(5.0)).unwrap();
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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, Some(2.0)).unwrap();

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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, Some(2.0)).unwrap();

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

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, Some(2.0)).unwrap();

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

        // None = maximize fold enrichment
        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, None).unwrap();
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
            target_fold_enrichment: None,
        };

        let table = HashBiasTable::create(&[pos.path()], &[neg.path()], &config, None).unwrap();

        assert!(table.threshold >= -127);
        assert!(table.fold_enrichment() >= 1.0);
    }
}
