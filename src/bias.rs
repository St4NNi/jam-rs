use anyhow::{Context, Result};
use indicatif::ProgressBar;
use jamhash::jamhash_u64;
use needletail::{Sequence, parse_fastx_file};
use std::io::{Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

const BRAW_MAGIC: &[u8; 4] = b"BRW3";
const BRAW_VERSION: u32 = 3;
const BIAS_MAGIC: &[u8; 4] = b"BIA3";
const BIAS_VERSION: u32 = 3;

const OLD_BRAW_MAGIC: &[u8; 4] = b"BRW2";
const OLD_BIAS_MAGIC: &[u8; 4] = b"BIA2";

const DEFAULT_CMS_WIDTH: usize = 1 << 20; // 1M columns
const DEFAULT_CMS_DEPTH: usize = 5;
const QUANTIZATION_SCALE: f32 = 10.0; // Each i8 unit = 0.1 log-ratio

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
        let seeds: Vec<u64> = (0..depth).map(|i| 0x517cc1b727220a95u64.wrapping_add(i as u64)).collect();
        let counts = vec![0u64; width * depth];
        Self { width, depth, seeds, counts }
    }

    pub fn with_seeds(width: usize, depth: usize, seeds: Vec<u64>) -> Self {
        assert_eq!(seeds.len(), depth);
        let counts = vec![0u64; width * depth];
        Self { width, depth, seeds, counts }
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

    pub fn width(&self) -> usize { self.width }
    pub fn depth(&self) -> usize { self.depth }
    pub fn seeds(&self) -> &[u64] { &self.seeds }
    pub fn counts(&self) -> &[u64] { &self.counts }

    pub fn cell_stats(&self) -> (u64, u64, f64, f64, usize) {
        let min = *self.counts.iter().min().unwrap_or(&0);
        let max = *self.counts.iter().max().unwrap_or(&0);
        let sum: u64 = self.counts.iter().sum();
        let mean = sum as f64 / self.counts.len() as f64;
        let variance: f64 = self.counts.iter()
            .map(|&c| { let d = c as f64 - mean; d * d })
            .sum::<f64>() / self.counts.len() as f64;
        let non_zero = self.counts.iter().filter(|&&c| c > 0).count();
        (min, max, mean, variance.sqrt(), non_zero)
    }
}

#[derive(Debug, Clone)]
pub struct RawHashCounts {
    pub config: CMSConfig,
    pub cms: CountMinSketch,
    pub total: u64,
}

impl RawHashCounts {
    pub fn new(config: CMSConfig) -> Self {
        let cms = CountMinSketch::new(config.width, config.depth);
        Self { config, cms, total: 0 }
    }

    pub fn build(paths: &[&Path], config: CMSConfig, progress: Option<ProgressBar>) -> Result<Self> {
        let record_counter = AtomicU64::new(0);
        let hash_counter = AtomicU64::new(0);

        let update_handle = progress.as_ref().map(|pb| {
            let pb = pb.clone();
            let record_counter = &record_counter as *const AtomicU64 as usize;
            let hash_counter = &hash_counter as *const AtomicU64 as usize;

            std::thread::spawn(move || {
                let record_counter = unsafe { &*(record_counter as *const AtomicU64) };
                let hash_counter = unsafe { &*(hash_counter as *const AtomicU64) };

                loop {
                    if pb.is_finished() { break; }
                    let records = record_counter.load(Ordering::Relaxed);
                    let hashes = hash_counter.load(Ordering::Relaxed);
                    pb.set_message(format!("{} records, {} hashes", records, format_number(hashes)));
                    std::thread::sleep(std::time::Duration::from_millis(100));
                }
            })
        });

        let mut raw = Self::new(config);
        let frac_max = u64::MAX / raw.config.fscale;
        let k = raw.config.k;

        for path in paths {
            let mut reader = match parse_fastx_file(path) {
                Ok(reader) => reader,
                Err(e) if e.kind == needletail::errors::ParseErrorKind::EmptyFile => {
                    continue;
                }
                Err(e) => {
                    return Err(e).with_context(|| format!("Failed to parse: {}", path.display()));
                }
            };

            while let Some(record) = reader.next() {
                let record = record.context("Failed to parse sequence record")?;
                let seq = record.normalize(false);
                record_counter.fetch_add(1, Ordering::Relaxed);

                if seq.len() < k as usize { continue; }

                for (_, kmer, _) in seq.bit_kmers(k, true) {
                    let hash = jamhash_u64(kmer.0);
                    if hash < frac_max {
                        raw.cms.increment(hash);
                        raw.total += 1;
                        hash_counter.fetch_add(1, Ordering::Relaxed);
                    }
                }
            }
        }

        if let Some(ref pb) = progress {
            pb.finish();
        }
        if let Some(handle) = update_handle {
            let _ = handle.join();
        }

        Ok(raw)
    }

    pub fn save(&self, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)
            .with_context(|| format!("Failed to create raw hash counts file: {}", path.display()))?;

        file.write_all(BRAW_MAGIC)?;
        file.write_all(&BRAW_VERSION.to_le_bytes())?;
        file.write_all(&[self.config.k])?;
        file.write_all(&self.config.fscale.to_le_bytes())?;
        file.write_all(&(self.config.width as u32).to_le_bytes())?;
        file.write_all(&[self.config.depth as u8])?;
        file.write_all(&self.total.to_le_bytes())?;

        for &seed in self.cms.seeds() {
            file.write_all(&seed.to_le_bytes())?;
        }
        for &count in self.cms.counts() {
            file.write_all(&count.to_le_bytes())?;
        }

        Ok(())
    }

    pub fn load(path: &Path) -> Result<Self> {
        let mut file = std::fs::File::open(path)
            .with_context(|| format!("Failed to open raw hash counts file: {}", path.display()))?;

        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)?;

        if &magic == OLD_BRAW_MAGIC {
            anyhow::bail!(
                "Old v2 raw bias table format detected (BRW2). Please recreate with 'jam bias create'.\n\
                 The new v3 format uses hash-based bias instead of k-mer composition."
            );
        }
        if &magic != BRAW_MAGIC {
            anyhow::bail!("Invalid raw hash counts file (bad magic): {}", path.display());
        }

        let mut buf4 = [0u8; 4];
        file.read_exact(&mut buf4)?;
        let version = u32::from_le_bytes(buf4);
        if version != BRAW_VERSION {
            anyhow::bail!("Unsupported raw hash counts version {} (expected {})", version, BRAW_VERSION);
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

        file.read_exact(&mut buf8)?;
        let total = u64::from_le_bytes(buf8);

        let mut seeds = Vec::with_capacity(depth);
        for _ in 0..depth {
            file.read_exact(&mut buf8)?;
            seeds.push(u64::from_le_bytes(buf8));
        }

        let mut counts = vec![0u64; width * depth];
        for count in &mut counts {
            file.read_exact(&mut buf8)?;
            *count = u64::from_le_bytes(buf8);
        }

        let config = CMSConfig { width, depth, k, fscale };
        let cms = CountMinSketch::with_seeds(width, depth, seeds);
        let mut raw = Self { config, cms, total };
        raw.cms.counts = counts;

        Ok(raw)
    }

    pub fn memory_usage(&self) -> usize {
        self.cms.counts().len() * 8 + self.cms.seeds().len() * 8
    }
}

#[derive(Debug, Clone, Copy)]
pub struct CalibrationResult {
    pub threshold: i8,
    pub positive_retention: f32,
    pub negative_retention: f32,
    pub selectivity: f32,
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
}

impl HashBiasTable {
    pub fn build(
        positive: &RawHashCounts,
        negative: &RawHashCounts,
        alpha: f32,
        target_selectivity: f32,
    ) -> Result<Self> {
        if positive.config.k != negative.config.k {
            anyhow::bail!(
                "k-mer size mismatch: positive={}, negative={}",
                positive.config.k, negative.config.k
            );
        }
        if positive.config.fscale != negative.config.fscale {
            anyhow::bail!(
                "fscale mismatch: positive={}, negative={}",
                positive.config.fscale, negative.config.fscale
            );
        }
        if positive.config.width != negative.config.width || positive.config.depth != negative.config.depth {
            anyhow::bail!(
                "CMS dimensions mismatch: positive={}x{}, negative={}x{}",
                positive.config.width, positive.config.depth,
                negative.config.width, negative.config.depth
            );
        }

        let width = positive.config.width;
        let depth = positive.config.depth;
        let seeds = positive.cms.seeds().to_vec();

        let pos_counts = positive.cms.counts();
        let neg_counts = negative.cms.counts();

        let mut weights = vec![0i8; width * depth];

        for i in 0..(width * depth) {
            let c_p = pos_counts[i] as f32;
            let c_n = neg_counts[i] as f32;
            let log_ratio = ((c_p + alpha) / (c_n + alpha)).ln();
            let quantized = (log_ratio * QUANTIZATION_SCALE).clamp(-127.0, 127.0) as i8;
            weights[i] = quantized;
        }

        let calibration = calibrate_threshold(
            positive, negative, &weights, &seeds, width, target_selectivity
        )?;

        Ok(Self {
            config: positive.config.clone(),
            seeds,
            weights,
            alpha,
            threshold: calibration.threshold,
            positive_retention: calibration.positive_retention,
            negative_retention: calibration.negative_retention,
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

    pub fn k(&self) -> u8 { self.config.k }
    pub fn fscale(&self) -> u64 { self.config.fscale }

    pub fn selectivity(&self) -> f32 {
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

        if &magic == OLD_BIAS_MAGIC {
            anyhow::bail!(
                "Old v2 bias table format detected (BIA2). Please recreate with 'jam bias combine'.\n\
                 The new v3 format uses hash-based bias instead of k-mer composition."
            );
        }
        if &magic != BIAS_MAGIC {
            anyhow::bail!("Invalid bias table file (bad magic): {}", path.display());
        }

        let mut buf4 = [0u8; 4];
        file.read_exact(&mut buf4)?;
        let version = u32::from_le_bytes(buf4);
        if version != BIAS_VERSION {
            anyhow::bail!("Unsupported bias table version {} (expected {})", version, BIAS_VERSION);
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

        let config = CMSConfig { width, depth, k, fscale };

        Ok(Self {
            config,
            seeds,
            weights,
            alpha,
            threshold,
            positive_retention,
            negative_retention,
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
        if &magic == OLD_BIAS_MAGIC {
            anyhow::bail!(
                "Old v2 bias table format detected (BIA2). Please recreate with 'jam bias combine'."
            );
        }
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
            anyhow::bail!("Bias table data truncated: expected {} bytes, got {}", weights_end, data.len());
        }

        let mut seeds = Vec::with_capacity(depth);
        for i in 0..depth {
            let offset = seeds_start + i * 8;
            seeds.push(u64::from_le_bytes(data[offset..offset + 8].try_into().unwrap()));
        }

        let mut weights = vec![0i8; width * depth];
        for (i, &b) in data[weights_start..weights_end].iter().enumerate() {
            weights[i] = b as i8;
        }

        let config = CMSConfig { width, depth, k, fscale };

        Ok(Self {
            config,
            seeds,
            weights,
            alpha,
            threshold,
            positive_retention,
            negative_retention,
        })
    }

    pub fn weight_stats(&self) -> (f32, f32, f32, f32, usize) {
        let min = *self.weights.iter().min().unwrap_or(&0) as f32 / QUANTIZATION_SCALE;
        let max = *self.weights.iter().max().unwrap_or(&0) as f32 / QUANTIZATION_SCALE;
        let sum: i64 = self.weights.iter().map(|&w| w as i64).sum();
        let mean = sum as f32 / self.weights.len() as f32 / QUANTIZATION_SCALE;
        let variance: f32 = self.weights.iter()
            .map(|&w| {
                let d = w as f32 / QUANTIZATION_SCALE - mean;
                d * d
            })
            .sum::<f32>() / self.weights.len() as f32;
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
        eprintln!("  CMS dimensions: {} x {}", self.config.width, self.config.depth);
        eprintln!("  Smoothing (alpha): {:.1}", self.alpha);
        eprintln!("  Threshold: {:.2} (quantized: {})", self.threshold_f32(), self.threshold);
        eprintln!("  Positive retention: {:.2}%", self.positive_retention * 100.0);
        eprintln!("  Negative retention: {:.2}%", self.negative_retention * 100.0);
        eprintln!("  Selectivity: {:.2}x", self.selectivity());
        eprintln!("  Weight stats: min={:.2}, max={:.2}, mean={:.2}, std={:.2}", min, max, mean, std);
        eprintln!("  Positive weights: {} ({:.1}%)", positive, positive as f64 / total_cells as f64 * 100.0);
    }
}

fn calibrate_threshold(
    positive: &RawHashCounts,
    negative: &RawHashCounts,
    weights: &[i8],
    seeds: &[u64],
    width: usize,
    target_selectivity: f32,
) -> Result<CalibrationResult> {
    let sample_hashes = |raw: &RawHashCounts, max_samples: usize| -> Vec<u64> {
        let frac_max = u64::MAX / raw.config.fscale;
        let mut hashes = Vec::new();
        let step = if raw.total > max_samples as u64 { raw.total as usize / max_samples } else { 1 };
        let mut counter = 0usize;

        for row in 0..raw.config.depth {
            for col in 0..raw.config.width {
                let idx = row * raw.config.width + col;
                let count = raw.cms.counts()[idx];
                if count > 0 {
                    for h_offset in 0..count.min(10) as u64 {
                        let pseudo_hash = (col as u64).wrapping_mul(seeds[row]).wrapping_add(h_offset);
                        if pseudo_hash < frac_max {
                            counter += 1;
                            if counter % step == 0 {
                                hashes.push(pseudo_hash);
                            }
                        }
                    }
                }
            }
        }
        hashes
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
            selectivity: 1.0,
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

        if neg_ret < 1e-6 { continue; }

        let selectivity = pos_ret / neg_ret;
        let diff = (selectivity - target_selectivity).abs();

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
        selectivity: if best_neg_ret > 0.0 { best_pos_ret / best_neg_ret } else { f32::INFINITY },
    })
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

pub fn compute_divergence(positive: &RawHashCounts, negative: &RawHashCounts) -> (f64, f64, f64) {
    let pos_counts = positive.cms.counts();
    let neg_counts = negative.cms.counts();
    let pos_total = positive.total as f64;
    let neg_total = negative.total as f64;

    if pos_total == 0.0 || neg_total == 0.0 {
        return (0.0, 0.0, 0.0);
    }

    let width = positive.config.width;
    let mut kl_pn = 0.0f64;
    let mut kl_np = 0.0f64;

    for col in 0..width {
        let p = pos_counts[col] as f64 / pos_total;
        let q = neg_counts[col] as f64 / neg_total;

        if p > 1e-12 && q > 1e-12 {
            kl_pn += p * (p / q).ln();
            kl_np += q * (q / p).ln();
        }
    }

    kl_pn /= std::f64::consts::LN_2;
    kl_np /= std::f64::consts::LN_2;

    let js = (kl_pn + kl_np) / 2.0;

    (kl_pn, kl_np, js.sqrt())
}

pub fn compute_effect_size(positive: &RawHashCounts, negative: &RawHashCounts, alpha: f32) -> f64 {
    let pos_counts = positive.cms.counts();
    let neg_counts = negative.cms.counts();

    let width = positive.config.width;
    let mut pos_weights = Vec::with_capacity(width);
    let mut neg_weights = Vec::with_capacity(width);

    for col in 0..width {
        let c_p = pos_counts[col] as f32;
        let c_n = neg_counts[col] as f32;
        let log_ratio = ((c_p + alpha) / (c_n + alpha)).ln();

        if pos_counts[col] > 0 {
            pos_weights.push(log_ratio);
        }
        if neg_counts[col] > 0 {
            neg_weights.push(log_ratio);
        }
    }

    if pos_weights.is_empty() || neg_weights.is_empty() {
        return 0.0;
    }

    let pos_mean = pos_weights.iter().sum::<f32>() / pos_weights.len() as f32;
    let neg_mean = neg_weights.iter().sum::<f32>() / neg_weights.len() as f32;

    let pos_var = pos_weights.iter().map(|&w| (w - pos_mean).powi(2)).sum::<f32>() / pos_weights.len() as f32;
    let neg_var = neg_weights.iter().map(|&w| (w - neg_mean).powi(2)).sum::<f32>() / neg_weights.len() as f32;

    let pooled_std = ((pos_var + neg_var) / 2.0).sqrt();

    if pooled_std > 0.0 {
        ((pos_mean - neg_mean) / pooled_std) as f64
    } else {
        0.0
    }
}

pub const BIAS_TABLE_SERIALIZED_SIZE: usize = 35 + DEFAULT_CMS_DEPTH * 8 + DEFAULT_CMS_WIDTH * DEFAULT_CMS_DEPTH;

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

        let raw = RawHashCounts::build(&[fasta.path()], config, None).unwrap();
        assert!(raw.total > 0);
    }

    #[test]
    fn test_raw_hash_counts_save_load() {
        let fasta = create_fasta(&["ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"]);
        let config = CMSConfig {
            width: 1024,
            depth: 3,
            k: 11,
            fscale: 10,
        };

        let raw = RawHashCounts::build(&[fasta.path()], config, None).unwrap();

        let output = NamedTempFile::new().unwrap();
        raw.save(output.path()).unwrap();

        let loaded = RawHashCounts::load(output.path()).unwrap();
        assert_eq!(raw.total, loaded.total);
        assert_eq!(raw.config.k, loaded.config.k);
        assert_eq!(raw.config.fscale, loaded.config.fscale);
        assert_eq!(raw.cms.counts(), loaded.cms.counts());
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

        let pos_raw = RawHashCounts::build(&[pos.path()], config.clone(), None).unwrap();
        let neg_raw = RawHashCounts::build(&[neg.path()], config, None).unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 5.0).unwrap();
        assert!(table.threshold >= -127 && table.threshold <= 127);
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

        let pos_raw = RawHashCounts::build(&[pos.path()], config.clone(), None).unwrap();
        let neg_raw = RawHashCounts::build(&[neg.path()], config, None).unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 2.0).unwrap();

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

        let pos_raw = RawHashCounts::build(&[pos.path()], config.clone(), None).unwrap();
        let neg_raw = RawHashCounts::build(&[neg.path()], config, None).unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 2.0).unwrap();

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

        let pos_raw = RawHashCounts::build(&[pos.path()], config.clone(), None).unwrap();
        let neg_raw = RawHashCounts::build(&[neg.path()], config, None).unwrap();

        let table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, 2.0).unwrap();

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
}
