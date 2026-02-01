use anyhow::{Context, Result};
use indicatif::ProgressBar;
use needletail::{Sequence, parse_fastx_file};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

const HEXAMER_COUNT: usize = 4096; // 4^6
const BIAS_MAGIC: &[u8; 4] = b"BIAS";
const BIAS_VERSION: u32 = 1;
const BRAW_MAGIC: &[u8; 4] = b"BRAW";
const BRAW_VERSION: u32 = 1;

/// Serialized size: magic(4) + version(4) + threshold(4) + scores(4096 * 4) = 16,396 bytes
pub const BIAS_TABLE_SERIALIZED_SIZE: usize = 4 + 4 + 4 + (HEXAMER_COUNT * 4);

/// Serialized size: magic(4) + version(4) + metadata(24) + counts(4096 * 8) = 32,800 bytes
pub const RAW_BIAS_TABLE_SERIALIZED_SIZE: usize = 4 + 4 + 24 + (HEXAMER_COUNT * 8);

#[derive(Debug, Clone, Default)]
pub struct FileStats {
    pub num_records: u64,
    pub total_bases: u64,
    pub gc_count: u64,
}

impl FileStats {
    pub fn gc_content(&self) -> f64 {
        if self.total_bases == 0 {
            0.0
        } else {
            self.gc_count as f64 / self.total_bases as f64 * 100.0
        }
    }
}

pub fn idx_to_hexamer(idx: usize) -> String {
    const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
    (0..6)
        .rev()
        .map(|i| BASES[(idx >> (2 * i)) & 0x3])
        .collect()
}

pub fn reverse_complement_idx(idx: usize) -> usize {
    let mut rc = 0usize;
    for i in 0..6 {
        let base = (idx >> (2 * i)) & 0x3;
        let comp = 3 - base; // A<->T, C<->G in 2-bit: 0<->3, 1<->2
        rc = (rc << 2) | comp;
    }
    rc
}

#[derive(Debug, Clone)]
pub struct RawBiasTable {
    counts: [u64; HEXAMER_COUNT],
    pub stats: FileStats,
}

impl RawBiasTable {
    pub fn build(path: &Path, progress: Option<ProgressBar>) -> Result<Self> {
        let record_counter = AtomicU64::new(0);
        let bases_counter = AtomicU64::new(0);

        let update_handle = progress.as_ref().map(|pb| {
            let pb = pb.clone();
            let record_counter = &record_counter as *const AtomicU64 as usize;
            let bases_counter = &bases_counter as *const AtomicU64 as usize;

            std::thread::spawn(move || {
                let record_counter = unsafe { &*(record_counter as *const AtomicU64) };
                let bases_counter = unsafe { &*(bases_counter as *const AtomicU64) };

                loop {
                    if pb.is_finished() {
                        break;
                    }

                    let records = record_counter.load(Ordering::Relaxed);
                    let bases = bases_counter.load(Ordering::Relaxed);

                    pb.set_message(format!(
                        "{} records, {}",
                        records,
                        format_bp(bases)
                    ));

                    std::thread::sleep(std::time::Duration::from_millis(100));
                }
            })
        });

        let (counts, stats) = count_hexamers(path, &record_counter, &bases_counter)
            .with_context(|| format!("Failed to count hexamers in: {}", path.display()))?;

        if let Some(ref pb) = progress {
            pb.finish();
        }

        if let Some(handle) = update_handle {
            let _ = handle.join();
        }

        Ok(Self { counts, stats })
    }

    pub fn save(&self, path: &Path) -> Result<()> {
        use std::io::Write;

        let mut file = std::fs::File::create(path)
            .with_context(|| format!("Failed to create raw bias table file: {}", path.display()))?;

        file.write_all(BRAW_MAGIC)?;
        file.write_all(&BRAW_VERSION.to_le_bytes())?;
        file.write_all(&self.stats.num_records.to_le_bytes())?;
        file.write_all(&self.stats.total_bases.to_le_bytes())?;
        file.write_all(&self.stats.gc_count.to_le_bytes())?;

        for count in &self.counts {
            file.write_all(&count.to_le_bytes())?;
        }

        Ok(())
    }

    pub fn load(path: &Path) -> Result<Self> {
        use std::io::Read;

        let mut file = std::fs::File::open(path)
            .with_context(|| format!("Failed to open raw bias table file: {}", path.display()))?;

        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)?;
        if &magic != BRAW_MAGIC {
            return Err(anyhow::anyhow!(
                "Invalid raw bias table file (bad magic bytes): {}",
                path.display()
            ));
        }

        let mut version_bytes = [0u8; 4];
        file.read_exact(&mut version_bytes)?;
        let version = u32::from_le_bytes(version_bytes);
        if version != BRAW_VERSION {
            return Err(anyhow::anyhow!(
                "Unsupported raw bias table version {} (expected {})",
                version,
                BRAW_VERSION
            ));
        }

        let mut num_records_bytes = [0u8; 8];
        file.read_exact(&mut num_records_bytes)?;
        let num_records = u64::from_le_bytes(num_records_bytes);

        let mut total_bases_bytes = [0u8; 8];
        file.read_exact(&mut total_bases_bytes)?;
        let total_bases = u64::from_le_bytes(total_bases_bytes);

        let mut gc_count_bytes = [0u8; 8];
        file.read_exact(&mut gc_count_bytes)?;
        let gc_count = u64::from_le_bytes(gc_count_bytes);

        let mut counts = [0u64; HEXAMER_COUNT];
        for count in &mut counts {
            let mut buf = [0u8; 8];
            file.read_exact(&mut buf)?;
            *count = u64::from_le_bytes(buf);
        }

        Ok(Self {
            counts,
            stats: FileStats {
                num_records,
                total_bases,
                gc_count,
            },
        })
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(RAW_BIAS_TABLE_SERIALIZED_SIZE);
        out.extend_from_slice(BRAW_MAGIC);
        out.extend_from_slice(&BRAW_VERSION.to_le_bytes());
        out.extend_from_slice(&self.stats.num_records.to_le_bytes());
        out.extend_from_slice(&self.stats.total_bases.to_le_bytes());
        out.extend_from_slice(&self.stats.gc_count.to_le_bytes());
        for count in &self.counts {
            out.extend_from_slice(&count.to_le_bytes());
        }
        out
    }

    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        if data.len() < RAW_BIAS_TABLE_SERIALIZED_SIZE {
            return Err(anyhow::anyhow!(
                "Raw bias table data too small: {} bytes (expected {})",
                data.len(),
                RAW_BIAS_TABLE_SERIALIZED_SIZE
            ));
        }

        let magic: [u8; 4] = data[0..4].try_into().unwrap();
        if &magic != BRAW_MAGIC {
            return Err(anyhow::anyhow!(
                "Invalid raw bias table magic bytes: {:?}",
                magic
            ));
        }

        let version = u32::from_le_bytes(data[4..8].try_into().unwrap());
        if version != BRAW_VERSION {
            return Err(anyhow::anyhow!(
                "Unsupported raw bias table version {} (expected {})",
                version,
                BRAW_VERSION
            ));
        }

        let num_records = u64::from_le_bytes(data[8..16].try_into().unwrap());
        let total_bases = u64::from_le_bytes(data[16..24].try_into().unwrap());
        let gc_count = u64::from_le_bytes(data[24..32].try_into().unwrap());

        let mut counts = [0u64; HEXAMER_COUNT];
        for (i, count) in counts.iter_mut().enumerate() {
            let offset = 32 + i * 8;
            *count = u64::from_le_bytes(data[offset..offset + 8].try_into().unwrap());
        }

        Ok(Self {
            counts,
            stats: FileStats {
                num_records,
                total_bases,
                gc_count,
            },
        })
    }

    pub fn total_hexamers(&self) -> u64 {
        self.counts.iter().sum()
    }

    pub fn frequency(&self, idx: usize) -> f64 {
        let total = self.total_hexamers();
        if total == 0 {
            0.0
        } else {
            self.counts[idx] as f64 / total as f64
        }
    }

    pub fn count(&self, idx: usize) -> u64 {
        self.counts[idx]
    }

    pub fn top_hexamers(&self, n: usize) -> Vec<(String, u64, f64)> {
        let total = self.total_hexamers();
        let mut indexed: Vec<(usize, u64)> = self.counts.iter().copied().enumerate().collect();
        indexed.sort_by(|a, b| b.1.cmp(&a.1));

        indexed
            .into_iter()
            .take(n)
            .map(|(idx, count)| {
                let freq = if total > 0 {
                    count as f64 / total as f64
                } else {
                    0.0
                };
                (idx_to_hexamer(idx), count, freq)
            })
            .collect()
    }

    pub fn shannon_entropy(&self) -> f64 {
        let total = self.total_hexamers();
        if total == 0 {
            return 0.0;
        }

        let mut entropy = 0.0;
        for &count in &self.counts {
            if count > 0 {
                let p = count as f64 / total as f64;
                entropy -= p * p.log2();
            }
        }
        entropy
    }

    pub fn non_zero_count(&self) -> usize {
        self.counts.iter().filter(|&&c| c > 0).count()
    }

    pub fn low_count_hexamers(&self, threshold: u64) -> usize {
        self.counts
            .iter()
            .filter(|&&c| c > 0 && c < threshold)
            .count()
    }

    pub fn combine(&self, other: &RawBiasTable, threshold: f32) -> Result<BiasTable> {
        if !(0.0..=1.0).contains(&threshold) {
            return Err(anyhow::anyhow!(
                "Threshold must be between 0.0 and 1.0, got {}",
                threshold
            ));
        }

        let pos_total = self.total_hexamers();
        let neg_total = other.total_hexamers();

        if pos_total == 0 {
            return Err(anyhow::anyhow!("No hexamers in positive table"));
        }
        if neg_total == 0 {
            return Err(anyhow::anyhow!("No hexamers in negative table"));
        }

        let mut scores = [0.0f32; HEXAMER_COUNT];
        for (i, score) in scores.iter_mut().enumerate() {
            let f_pos = self.counts[i] as f64 / pos_total as f64;
            let f_neg = other.counts[i] as f64 / neg_total as f64;
            let raw_score = (f_pos / (f_pos + f_neg + 1e-12)) as f32;
            *score = raw_score - threshold;
        }

        Ok(BiasTable { scores, threshold })
    }
}

#[derive(Debug, Clone)]
pub struct CompareStats {
    pub kl_divergence: f64,
    pub gc_diff: f64,
    pub threshold_sweep: Vec<ThresholdPoint>,
    pub top_positive_enriched: Vec<EnrichedHexamer>,
    pub top_negative_enriched: Vec<EnrichedHexamer>,
    pub hexamers_favoring_positive: usize,
    pub hexamers_favoring_negative: usize,
}

#[derive(Debug, Clone)]
pub struct ThresholdPoint {
    pub threshold: f32,
    pub pos_retain: f64,
    pub neg_retain: f64,
    pub selectivity: f64,
}

#[derive(Debug, Clone)]
pub struct EnrichedHexamer {
    pub hexamer: String,
    pub rev_comp: String,
    pub p_pos_given_hex: f64,
    pub pos_freq: f64,
    pub neg_freq: f64,
    pub fold_change: f64,
}

impl RawBiasTable {
    pub fn compare(&self, other: &RawBiasTable, threshold: f32) -> CompareStats {
        let pos_total = self.total_hexamers() as f64;
        let neg_total = other.total_hexamers() as f64;

        let kl_divergence = if pos_total > 0.0 && neg_total > 0.0 {
            let mut kl = 0.0;
            for i in 0..HEXAMER_COUNT {
                let p = self.counts[i] as f64 / pos_total;
                let q = other.counts[i] as f64 / neg_total;
                if p > 0.0 && q > 0.0 {
                    kl += p * (p / q).ln();
                }
            }
            kl / std::f64::consts::LN_2
        } else {
            0.0
        };

        let gc_diff = self.stats.gc_content() - other.stats.gc_content();

        let thresholds = [0.3f32, 0.4, 0.5, 0.6, 0.7];
        let threshold_sweep: Vec<ThresholdPoint> = thresholds
            .iter()
            .map(|&t| {
                let (pos_retain, neg_retain) =
                    self.estimate_retention(other, t, pos_total, neg_total);
                let selectivity = if neg_retain > 0.0 {
                    pos_retain / neg_retain
                } else {
                    f64::INFINITY
                };
                ThresholdPoint {
                    threshold: t,
                    pos_retain: pos_retain * 100.0,
                    neg_retain: neg_retain * 100.0,
                    selectivity,
                }
            })
            .collect();

        let mut enrichment: Vec<(usize, f64, f64, f64, f64)> = (0..HEXAMER_COUNT)
            .map(|i| {
                let pos_freq = if pos_total > 0.0 {
                    self.counts[i] as f64 / pos_total
                } else {
                    0.0
                };
                let neg_freq = if neg_total > 0.0 {
                    other.counts[i] as f64 / neg_total
                } else {
                    0.0
                };
                let p_pos = pos_freq / (pos_freq + neg_freq + 1e-12);
                let fold = if neg_freq > 0.0 {
                    pos_freq / neg_freq
                } else if pos_freq > 0.0 {
                    f64::INFINITY
                } else {
                    1.0
                };
                (i, p_pos, pos_freq, neg_freq, fold)
            })
            .collect();

        enrichment.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        let top_positive_enriched: Vec<EnrichedHexamer> = enrichment
            .iter()
            .take(10)
            .map(|&(idx, p_pos, pos_freq, neg_freq, fold)| {
                let rc_idx = reverse_complement_idx(idx);
                EnrichedHexamer {
                    hexamer: idx_to_hexamer(idx),
                    rev_comp: idx_to_hexamer(rc_idx),
                    p_pos_given_hex: p_pos,
                    pos_freq: pos_freq * 100.0,
                    neg_freq: neg_freq * 100.0,
                    fold_change: fold,
                }
            })
            .collect();

        let top_negative_enriched: Vec<EnrichedHexamer> = enrichment
            .iter()
            .rev()
            .take(10)
            .map(|&(idx, p_pos, pos_freq, neg_freq, fold)| {
                let rc_idx = reverse_complement_idx(idx);
                EnrichedHexamer {
                    hexamer: idx_to_hexamer(idx),
                    rev_comp: idx_to_hexamer(rc_idx),
                    p_pos_given_hex: p_pos,
                    pos_freq: pos_freq * 100.0,
                    neg_freq: neg_freq * 100.0,
                    fold_change: fold,
                }
            })
            .collect();

        let hexamers_favoring_positive = enrichment
            .iter()
            .filter(|(_, p_pos, _, _, _)| *p_pos > threshold as f64)
            .count();
        let hexamers_favoring_negative = HEXAMER_COUNT - hexamers_favoring_positive;

        CompareStats {
            kl_divergence,
            gc_diff,
            threshold_sweep,
            top_positive_enriched,
            top_negative_enriched,
            hexamers_favoring_positive,
            hexamers_favoring_negative,
        }
    }

    fn estimate_retention(
        &self,
        other: &RawBiasTable,
        threshold: f32,
        pos_total: f64,
        neg_total: f64,
    ) -> (f64, f64) {
        let mut pos_retain_sum = 0.0;
        let mut neg_retain_sum = 0.0;

        for i in 0..HEXAMER_COUNT {
            let pos_freq = self.counts[i] as f64 / pos_total.max(1.0);
            let neg_freq = other.counts[i] as f64 / neg_total.max(1.0);
            let score = pos_freq / (pos_freq + neg_freq + 1e-12);

            if score >= threshold as f64 {
                pos_retain_sum += pos_freq;
                neg_retain_sum += neg_freq;
            }
        }

        (pos_retain_sum, neg_retain_sum)
    }
}

#[derive(Debug, Clone)]
pub struct BiasResult {
    pub table: BiasTable,
    pub positive_stats: FileStats,
    pub negative_stats: FileStats,
}

impl BiasResult {
    pub fn print_stats(&self) {
        self.print_file_stats();
        self.table.print_stats();
    }

    pub fn print_file_stats(&self) {
        eprintln!("Input file statistics:");
        eprintln!(
            "  Positive: {} records, {} bp, {:.1}% GC",
            self.positive_stats.num_records,
            format_bp(self.positive_stats.total_bases),
            self.positive_stats.gc_content()
        );
        eprintln!(
            "  Negative: {} records, {} bp, {:.1}% GC",
            self.negative_stats.num_records,
            format_bp(self.negative_stats.total_bases),
            self.negative_stats.gc_content()
        );
    }
}

fn format_bp(bp: u64) -> String {
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

#[derive(Debug, Clone)]
pub struct BiasTable {
    /// Pre-adjusted hexamer scores (raw_score - threshold).
    /// Indexed directly by 2-bit encoded hexamer value.
    /// At runtime: sum(scores) >= 0 means k-mer passes.
    scores: [f32; HEXAMER_COUNT],
    pub threshold: f32,
}

impl BiasTable {
    pub fn build(
        pos_path: &Path,
        neg_path: &Path,
        threshold: f32,
        progress: Option<ProgressBar>,
    ) -> Result<BiasResult> {
        if !(0.0..=1.0).contains(&threshold) {
            return Err(anyhow::anyhow!(
                "Threshold must be between 0.0 and 1.0, got {}",
                threshold
            ));
        }

        let pos_records = AtomicU64::new(0);
        let neg_records = AtomicU64::new(0);
        let pos_bases = AtomicU64::new(0);
        let neg_bases = AtomicU64::new(0);

        let update_handle = progress.as_ref().map(|pb| {
            let pb = pb.clone();
            let pos_records = &pos_records as *const AtomicU64 as usize;
            let neg_records = &neg_records as *const AtomicU64 as usize;
            let pos_bases = &pos_bases as *const AtomicU64 as usize;
            let neg_bases = &neg_bases as *const AtomicU64 as usize;

            std::thread::spawn(move || {
                // SAFETY: pointers are valid for the duration of build()
                let pos_records = unsafe { &*(pos_records as *const AtomicU64) };
                let neg_records = unsafe { &*(neg_records as *const AtomicU64) };
                let pos_bases = unsafe { &*(pos_bases as *const AtomicU64) };
                let neg_bases = unsafe { &*(neg_bases as *const AtomicU64) };

                loop {
                    let pr = pos_records.load(Ordering::Relaxed);
                    let nr = neg_records.load(Ordering::Relaxed);
                    let pb_val = pos_bases.load(Ordering::Relaxed);
                    let nb = neg_bases.load(Ordering::Relaxed);

                    pb.set_message(format!(
                        "pos: {} records, {} | neg: {} records, {}",
                        pr,
                        format_bp(pb_val),
                        nr,
                        format_bp(nb)
                    ));

                    std::thread::sleep(std::time::Duration::from_millis(100));

                    if pr > 0 && nr > 0 && pb.is_finished() {
                        break;
                    }
                }
            })
        });

        let (pos_result, neg_result) = rayon::join(
            || {
                count_hexamers(pos_path, &pos_records, &pos_bases).with_context(|| {
                    format!(
                        "Failed to count hexamers in positive file: {}",
                        pos_path.display()
                    )
                })
            },
            || {
                count_hexamers(neg_path, &neg_records, &neg_bases).with_context(|| {
                    format!(
                        "Failed to count hexamers in negative file: {}",
                        neg_path.display()
                    )
                })
            },
        );

        if let Some(ref pb) = progress {
            pb.set_message("Computing bias scores...");
        }

        if let Some(handle) = update_handle {
            let _ = handle.join();
        }

        let (pos_counts, positive_stats) = pos_result?;
        let (neg_counts, negative_stats) = neg_result?;

        let pos_total: u64 = pos_counts.iter().sum();
        let neg_total: u64 = neg_counts.iter().sum();

        if pos_total == 0 {
            return Err(anyhow::anyhow!("No hexamers found in positive file"));
        }
        if neg_total == 0 {
            return Err(anyhow::anyhow!("No hexamers found in negative file"));
        }

        let mut scores = [0.0f32; HEXAMER_COUNT];
        for i in 0..HEXAMER_COUNT {
            let f_pos = pos_counts[i] as f64 / pos_total as f64;
            let f_neg = neg_counts[i] as f64 / neg_total as f64;

            // Normalized score: P(positive | hexamer) in [0, 1]
            let raw_score = (f_pos / (f_pos + f_neg + 1e-12)) as f32;
            scores[i] = raw_score - threshold;
        }

        let table = Self { scores, threshold };
        Ok(BiasResult {
            table,
            positive_stats,
            negative_stats,
        })
    }

    #[inline(always)]
    pub fn passes_filter(&self, kmer: u64, kmer_len: u8) -> bool {
        let n_hexamers = kmer_len.saturating_sub(5) as usize;
        if n_hexamers == 0 {
            return true;
        }

        let mut sum: f32 = 0.0;
        for i in 0..n_hexamers {
            let hexamer = ((kmer >> (2 * i)) & 0xFFF) as usize;
            sum += self.scores[hexamer];
        }
        sum >= 0.0
    }

    pub fn score(&self, idx: usize) -> f32 {
        self.scores[idx]
    }

    pub fn save(&self, path: &Path) -> Result<()> {
        use std::io::Write;

        let mut file = std::fs::File::create(path)
            .with_context(|| format!("Failed to create bias table file: {}", path.display()))?;

        file.write_all(BIAS_MAGIC)?;
        file.write_all(&BIAS_VERSION.to_le_bytes())?;
        file.write_all(&self.threshold.to_le_bytes())?;

        for score in &self.scores {
            file.write_all(&score.to_le_bytes())?;
        }

        Ok(())
    }

    pub fn load(path: &Path) -> Result<Self> {
        use std::io::Read;

        let mut file = std::fs::File::open(path)
            .with_context(|| format!("Failed to open bias table file: {}", path.display()))?;

        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)?;
        if &magic != BIAS_MAGIC {
            return Err(anyhow::anyhow!(
                "Invalid bias table file (bad magic bytes): {}",
                path.display()
            ));
        }

        let mut version_bytes = [0u8; 4];
        file.read_exact(&mut version_bytes)?;
        let version = u32::from_le_bytes(version_bytes);
        if version != BIAS_VERSION {
            return Err(anyhow::anyhow!(
                "Unsupported bias table version {} (expected {})",
                version,
                BIAS_VERSION
            ));
        }

        let mut threshold_bytes = [0u8; 4];
        file.read_exact(&mut threshold_bytes)?;
        let threshold = f32::from_le_bytes(threshold_bytes);

        let mut scores = [0.0f32; HEXAMER_COUNT];
        for score in &mut scores {
            let mut buf = [0u8; 4];
            file.read_exact(&mut buf)?;
            *score = f32::from_le_bytes(buf);
        }

        Ok(Self { scores, threshold })
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(BIAS_TABLE_SERIALIZED_SIZE);
        out.extend_from_slice(BIAS_MAGIC);
        out.extend_from_slice(&BIAS_VERSION.to_le_bytes());
        out.extend_from_slice(&self.threshold.to_le_bytes());
        for score in &self.scores {
            out.extend_from_slice(&score.to_le_bytes());
        }
        out
    }

    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        if data.len() < BIAS_TABLE_SERIALIZED_SIZE {
            return Err(anyhow::anyhow!(
                "Bias table data too small: {} bytes (expected {})",
                data.len(),
                BIAS_TABLE_SERIALIZED_SIZE
            ));
        }

        let magic: [u8; 4] = data[0..4].try_into().unwrap();
        if &magic != BIAS_MAGIC {
            return Err(anyhow::anyhow!(
                "Invalid bias table magic bytes: {:?}",
                magic
            ));
        }

        let version = u32::from_le_bytes(data[4..8].try_into().unwrap());
        if version != BIAS_VERSION {
            return Err(anyhow::anyhow!(
                "Unsupported bias table version {} (expected {})",
                version,
                BIAS_VERSION
            ));
        }

        let threshold = f32::from_le_bytes(data[8..12].try_into().unwrap());

        let mut scores = [0.0f32; HEXAMER_COUNT];
        for (i, score) in scores.iter_mut().enumerate() {
            let offset = 12 + i * 4;
            *score = f32::from_le_bytes(data[offset..offset + 4].try_into().unwrap());
        }

        Ok(Self { scores, threshold })
    }

    pub fn print_stats(&self) {
        let mut positive_count = 0usize;
        let mut negative_count = 0usize;

        for &score in &self.scores {
            if score > 0.0 {
                positive_count += 1;
            } else if score < 0.0 {
                negative_count += 1;
            }
        }

        eprintln!("Bias table statistics:");
        eprintln!("  Threshold: {:.2}", self.threshold);
        eprintln!(
            "  Hexamers above threshold: {}/{} ({:.1}%)",
            positive_count,
            HEXAMER_COUNT,
            positive_count as f64 / HEXAMER_COUNT as f64 * 100.0
        );
        eprintln!(
            "  Hexamers below threshold: {}/{} ({:.1}%)",
            negative_count,
            HEXAMER_COUNT,
            negative_count as f64 / HEXAMER_COUNT as f64 * 100.0
        );
    }
}

impl PartialEq for BiasTable {
    fn eq(&self, other: &Self) -> bool {
        self.threshold == other.threshold && self.scores == other.scores
    }
}

fn count_hexamers(
    path: &Path,
    record_counter: &AtomicU64,
    bases_counter: &AtomicU64,
) -> Result<([u64; HEXAMER_COUNT], FileStats)> {
    let mut reader = match parse_fastx_file(path) {
        Ok(reader) => reader,
        Err(e) if e.kind == needletail::errors::ParseErrorKind::EmptyFile => {
            anyhow::bail!("Empty file: {}", path.display());
        }
        Err(e) => {
            return Err(e)
                .with_context(|| format!("Failed to parse FASTA file: {}", path.display()));
        }
    };

    let mut counts = [0u64; HEXAMER_COUNT];
    let mut stats = FileStats::default();

    while let Some(record) = reader.next() {
        let record = record.context("Failed to parse sequence record")?;
        let seq = record.normalize(false);

        bases_counter.fetch_add(seq.len() as u64, Ordering::Relaxed);
        record_counter.fetch_add(1, Ordering::Relaxed);

        for (_, kmer, _) in seq.bit_kmers(6, true) {
            let idx = (kmer.0 & 0xFFF) as usize;
            counts[idx] += 1;
        }

        stats.num_records += 1;
        stats.total_bases += seq.len() as u64;
        for &base in seq.iter() {
            if base == b'G' || base == b'C' || base == b'g' || base == b'c' {
                stats.gc_count += 1;
            }
        }
    }

    Ok((counts, stats))
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
    fn test_bias_table_build() {
        let pos = create_fasta(&["AAATATATATATATATATATATAT", "TATATATATATATATATATATATA"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCG", "CGCGCGCGCGCGCGCGCGCGCGC"]);

        let result = BiasTable::build(pos.path(), neg.path(), 0.5, None).unwrap();
        assert_eq!(result.table.threshold, 0.5);
        assert_eq!(result.positive_stats.num_records, 2);
        assert_eq!(result.negative_stats.num_records, 2);
    }

    #[test]
    fn test_save_load_roundtrip() {
        let pos = create_fasta(&["ATCGATCGATCGATCGATCGATCG"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGC"]);

        let result = BiasTable::build(pos.path(), neg.path(), 0.6, None).unwrap();

        let output = NamedTempFile::new().unwrap();
        result.table.save(output.path()).unwrap();

        let loaded = BiasTable::load(output.path()).unwrap();
        assert_eq!(result.table.threshold, loaded.threshold);
        assert_eq!(result.table.scores, loaded.scores);
    }

    #[test]
    fn test_to_bytes_from_bytes_roundtrip() {
        let pos = create_fasta(&["ATCGATCGATCGATCGATCGATCG"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGC"]);

        let result = BiasTable::build(pos.path(), neg.path(), 0.6, None).unwrap();
        let bytes = result.table.to_bytes();

        assert_eq!(bytes.len(), BIAS_TABLE_SERIALIZED_SIZE);

        let loaded = BiasTable::from_bytes(&bytes).unwrap();
        assert_eq!(result.table, loaded);
    }

    #[test]
    fn test_partial_eq() {
        let pos = create_fasta(&["ATCGATCGATCGATCGATCGATCG"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGC"]);

        let result1 = BiasTable::build(pos.path(), neg.path(), 0.6, None).unwrap();
        let result2 = BiasTable::build(pos.path(), neg.path(), 0.6, None).unwrap();
        let result3 = BiasTable::build(pos.path(), neg.path(), 0.7, None).unwrap();

        assert_eq!(result1.table, result2.table);
        assert_ne!(result1.table, result3.table);
    }

    #[test]
    fn test_passes_filter() {
        let pos = create_fasta(&["AAATATATATATATATATATATAT", "AATATAATATAATATAATATAATAT"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCG", "CGCGCGCGCGCGCGCGCGCGCGC"]);

        let result = BiasTable::build(pos.path(), neg.path(), 0.5, None).unwrap();

        // GC-rich k-mer should NOT pass (enriched in negative set)
        // G=0b10, C=0b01 -> GCGCGC... = repeating 01_10 pattern
        let gc_kmer: u64 = 0x6666666666; // GCGCGC pattern in 2-bit
        assert!(!result.table.passes_filter(gc_kmer, 21));

        // Build an ATATAT k-mer (A=00, T=11) which matches positive set
        // ATATAT = 00_11_00_11_00_11 per hexamer = 0xCCC
        let mut at_kmer: u64 = 0;
        for i in 0..21u32 {
            if i % 2 == 1 {
                at_kmer |= 0b11 << (2 * i);
            }
        }
        assert!(result.table.passes_filter(at_kmer, 21));
    }

    #[test]
    fn test_threshold_validation() {
        let pos = create_fasta(&["ATCGATCGATCG"]);
        let neg = create_fasta(&["GCGCGCGCGCGC"]);

        assert!(BiasTable::build(pos.path(), neg.path(), -0.1, None).is_err());
        assert!(BiasTable::build(pos.path(), neg.path(), 1.1, None).is_err());
    }

    #[test]
    fn test_short_kmer_passes() {
        let pos = create_fasta(&["ATCGATCGATCG"]);
        let neg = create_fasta(&["GCGCGCGCGCGC"]);

        let result = BiasTable::build(pos.path(), neg.path(), 0.5, None).unwrap();
        assert!(result.table.passes_filter(0, 5));
        assert!(result.table.passes_filter(0, 4));
    }

    #[test]
    fn test_file_stats() {
        let pos = create_fasta(&["AAATTTCCC", "GGGAAATTT"]); // 18 bp, 6 GC
        let neg = create_fasta(&["GCGCGCGCGC"]); // 10 bp, 10 GC (100%)

        let result = BiasTable::build(pos.path(), neg.path(), 0.5, None).unwrap();

        assert_eq!(result.positive_stats.num_records, 2);
        assert_eq!(result.positive_stats.total_bases, 18);
        assert_eq!(result.positive_stats.gc_count, 6);
        assert!((result.positive_stats.gc_content() - 33.33).abs() < 0.1);

        assert_eq!(result.negative_stats.num_records, 1);
        assert_eq!(result.negative_stats.total_bases, 10);
        assert_eq!(result.negative_stats.gc_count, 10);
        assert!((result.negative_stats.gc_content() - 100.0).abs() < 0.01);
    }

    #[test]
    fn test_idx_to_hexamer() {
        assert_eq!(idx_to_hexamer(0), "AAAAAA");
        assert_eq!(idx_to_hexamer(1), "AAAAAC");
        assert_eq!(idx_to_hexamer(2), "AAAAAG");
        assert_eq!(idx_to_hexamer(3), "AAAAAT");
        assert_eq!(idx_to_hexamer(4095), "TTTTTT");
    }

    #[test]
    fn test_reverse_complement_idx() {
        let aaaaaa_idx = 0;
        let tttttt_idx = 4095;
        assert_eq!(reverse_complement_idx(aaaaaa_idx), tttttt_idx);
        assert_eq!(reverse_complement_idx(tttttt_idx), aaaaaa_idx);

        let acgtac_idx = 0b00_01_10_11_00_01;
        let gtacgt_idx = reverse_complement_idx(acgtac_idx);
        assert_eq!(idx_to_hexamer(gtacgt_idx), "GTACGT");
    }

    #[test]
    fn test_raw_bias_table_build() {
        let fasta = create_fasta(&["ATCGATCGATCGATCGATCGATCG"]);
        let raw = RawBiasTable::build(fasta.path(), None).unwrap();

        assert_eq!(raw.stats.num_records, 1);
        assert_eq!(raw.stats.total_bases, 24);
        assert!(raw.total_hexamers() > 0);
    }

    #[test]
    fn test_raw_bias_table_save_load_roundtrip() {
        let fasta = create_fasta(&["ATCGATCGATCGATCGATCGATCG"]);
        let raw = RawBiasTable::build(fasta.path(), None).unwrap();

        let output = NamedTempFile::new().unwrap();
        raw.save(output.path()).unwrap();

        let loaded = RawBiasTable::load(output.path()).unwrap();
        assert_eq!(raw.stats.num_records, loaded.stats.num_records);
        assert_eq!(raw.stats.total_bases, loaded.stats.total_bases);
        assert_eq!(raw.stats.gc_count, loaded.stats.gc_count);
        assert_eq!(raw.counts, loaded.counts);
    }

    #[test]
    fn test_raw_bias_table_to_bytes_from_bytes_roundtrip() {
        let fasta = create_fasta(&["ATCGATCGATCGATCGATCGATCG"]);
        let raw = RawBiasTable::build(fasta.path(), None).unwrap();
        let bytes = raw.to_bytes();

        assert_eq!(bytes.len(), RAW_BIAS_TABLE_SERIALIZED_SIZE);

        let loaded = RawBiasTable::from_bytes(&bytes).unwrap();
        assert_eq!(raw.stats.num_records, loaded.stats.num_records);
        assert_eq!(raw.stats.total_bases, loaded.stats.total_bases);
        assert_eq!(raw.counts, loaded.counts);
    }

    #[test]
    fn test_raw_bias_table_combine() {
        let pos = create_fasta(&["AAATATATATATATATATATATAT"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCG"]);

        let pos_raw = RawBiasTable::build(pos.path(), None).unwrap();
        let neg_raw = RawBiasTable::build(neg.path(), None).unwrap();

        let combined = pos_raw.combine(&neg_raw, 0.5).unwrap();
        assert_eq!(combined.threshold, 0.5);
    }

    #[test]
    fn test_raw_bias_table_shannon_entropy() {
        let fasta = create_fasta(&["ATCGATCGATCGATCGATCGATCG"]);
        let raw = RawBiasTable::build(fasta.path(), None).unwrap();

        let entropy = raw.shannon_entropy();
        assert!(entropy > 0.0);
        assert!(entropy <= 12.0);
    }

    #[test]
    fn test_raw_bias_table_top_hexamers() {
        let fasta = create_fasta(&["AAAAAAAAAAAAAAAAAAAAAAAAA"]);
        let raw = RawBiasTable::build(fasta.path(), None).unwrap();

        let top = raw.top_hexamers(3);
        assert!(!top.is_empty());
        assert_eq!(top[0].0, "AAAAAA");
    }
}
