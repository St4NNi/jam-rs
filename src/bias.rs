use anyhow::{Context, Result};
use needletail::{Sequence, parse_fastx_file};
use std::path::Path;

const HEXAMER_COUNT: usize = 4096; // 4^6
const BIAS_MAGIC: &[u8; 4] = b"BIAS";
const BIAS_VERSION: u32 = 1;

/// Serialized size: magic(4) + version(4) + threshold(4) + scores(4096 * 4) = 16,396 bytes
pub const BIAS_TABLE_SERIALIZED_SIZE: usize = 4 + 4 + 4 + (HEXAMER_COUNT * 4);

#[derive(Debug, Clone)]
pub struct BiasTable {
    /// Pre-adjusted hexamer scores (raw_score - threshold).
    /// Indexed directly by 2-bit encoded hexamer value.
    /// At runtime: sum(scores) >= 0 means k-mer passes.
    scores: [f32; HEXAMER_COUNT],
    pub threshold: f32,
}

impl BiasTable {
    pub fn build(pos_path: &Path, neg_path: &Path, threshold: f32) -> Result<Self> {
        if !(0.0..=1.0).contains(&threshold) {
            return Err(anyhow::anyhow!(
                "Threshold must be between 0.0 and 1.0, got {}",
                threshold
            ));
        }

        let pos_counts = count_hexamers(pos_path)
            .with_context(|| format!("Failed to count hexamers in positive file: {}", pos_path.display()))?;
        let neg_counts = count_hexamers(neg_path)
            .with_context(|| format!("Failed to count hexamers in negative file: {}", neg_path.display()))?;

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

        Ok(Self { scores, threshold })
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

    /// Serialize the bias table to bytes for embedding in JAM files.
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

    /// Deserialize a bias table from a byte slice (e.g., from mmap).
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

fn count_hexamers(path: &Path) -> Result<[u64; HEXAMER_COUNT]> {
    let mut counts = [0u64; HEXAMER_COUNT];

    let mut reader = parse_fastx_file(path)
        .with_context(|| format!("Failed to parse FASTA file: {}", path.display()))?;

    while let Some(record) = reader.next() {
        let record = record.context("Failed to parse sequence record")?;
        let sequence = record.normalize(false);

        for (_, kmer, _) in sequence.bit_kmers(6, true) {
            let idx = (kmer.0 & 0xFFF) as usize;
            counts[idx] += 1;
        }
    }

    Ok(counts)
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
        let pos = create_fasta(&[
            "AAATATATATATATATATATATAT",
            "TATATATATATATATATATATATA",
        ]);
        let neg = create_fasta(&[
            "GCGCGCGCGCGCGCGCGCGCGCG",
            "CGCGCGCGCGCGCGCGCGCGCGC",
        ]);

        let table = BiasTable::build(pos.path(), neg.path(), 0.5).unwrap();
        assert_eq!(table.threshold, 0.5);
    }

    #[test]
    fn test_save_load_roundtrip() {
        let pos = create_fasta(&["ATCGATCGATCGATCGATCGATCG"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGC"]);

        let table = BiasTable::build(pos.path(), neg.path(), 0.6).unwrap();

        let output = NamedTempFile::new().unwrap();
        table.save(output.path()).unwrap();

        let loaded = BiasTable::load(output.path()).unwrap();
        assert_eq!(table.threshold, loaded.threshold);
        assert_eq!(table.scores, loaded.scores);
    }

    #[test]
    fn test_to_bytes_from_bytes_roundtrip() {
        let pos = create_fasta(&["ATCGATCGATCGATCGATCGATCG"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGC"]);

        let table = BiasTable::build(pos.path(), neg.path(), 0.6).unwrap();
        let bytes = table.to_bytes();

        assert_eq!(bytes.len(), BIAS_TABLE_SERIALIZED_SIZE);

        let loaded = BiasTable::from_bytes(&bytes).unwrap();
        assert_eq!(table, loaded);
    }

    #[test]
    fn test_partial_eq() {
        let pos = create_fasta(&["ATCGATCGATCGATCGATCGATCG"]);
        let neg = create_fasta(&["GCGCGCGCGCGCGCGCGCGCGCGC"]);

        let table1 = BiasTable::build(pos.path(), neg.path(), 0.6).unwrap();
        let table2 = BiasTable::build(pos.path(), neg.path(), 0.6).unwrap();
        let table3 = BiasTable::build(pos.path(), neg.path(), 0.7).unwrap();

        assert_eq!(table1, table2);
        assert_ne!(table1, table3);
    }

    #[test]
    fn test_passes_filter() {
        let pos = create_fasta(&[
            "AAATATATATATATATATATATAT",
            "AATATAATATAATATAATATAATAT",
        ]);
        let neg = create_fasta(&[
            "GCGCGCGCGCGCGCGCGCGCGCG",
            "CGCGCGCGCGCGCGCGCGCGCGC",
        ]);

        let table = BiasTable::build(pos.path(), neg.path(), 0.5).unwrap();

        // GC-rich k-mer should NOT pass (enriched in negative set)
        // G=0b10, C=0b01 -> GCGCGC... = repeating 01_10 pattern
        let gc_kmer: u64 = 0x6666666666; // GCGCGC pattern in 2-bit
        assert!(!table.passes_filter(gc_kmer, 21));

        // Build an ATATAT k-mer (A=00, T=11) which matches positive set
        // ATATAT = 00_11_00_11_00_11 per hexamer = 0xCCC
        let mut at_kmer: u64 = 0;
        for i in 0..21u32 {
            if i % 2 == 1 {
                at_kmer |= 0b11 << (2 * i);
            }
        }
        assert!(table.passes_filter(at_kmer, 21));
    }

    #[test]
    fn test_threshold_validation() {
        let pos = create_fasta(&["ATCGATCGATCG"]);
        let neg = create_fasta(&["GCGCGCGCGCGC"]);

        assert!(BiasTable::build(pos.path(), neg.path(), -0.1).is_err());
        assert!(BiasTable::build(pos.path(), neg.path(), 1.1).is_err());
    }

    #[test]
    fn test_short_kmer_passes() {
        let pos = create_fasta(&["ATCGATCGATCG"]);
        let neg = create_fasta(&["GCGCGCGCGCGC"]);

        let table = BiasTable::build(pos.path(), neg.path(), 0.5).unwrap();
        assert!(table.passes_filter(0, 5));
        assert!(table.passes_filter(0, 4));
    }
}
