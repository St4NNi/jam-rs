use anyhow::Result;
use glob::glob;
use rayon::prelude::*;
use serde::Serialize;
use std::cmp::{max, min};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

const S_LIST: [u64; 4] = [10, 100, 1000, 10000];

// Tune these for your storage system
const READER_CAP_BYTES: usize = 1024 * 1024 * 256; // 256 MiB BufReader internal buffer
const IO_CHUNK_BYTES: usize = 1024 * 1024 * 256; // 256 MiB read chunks

#[derive(Clone)]
struct Stats {
    n: u128,
    minv: u64,
    maxv: u64,
    // byte_hist[pos][byte], pos=0 is least significant byte
    byte_hist: [[u64; 256]; 8],
    // pass counts for thresholds at S_LIST
    pass_counts: [u128; S_LIST.len()],
}

impl Stats {
    fn new() -> Self {
        Self {
            n: 0,
            minv: u64::MAX,
            maxv: 0,
            byte_hist: [[0u64; 256]; 8],
            pass_counts: [0u128; S_LIST.len()],
        }
    }

    #[inline(always)]
    fn update(&mut self, h: u64) {
        self.n += 1;
        self.minv = min(self.minv, h);
        self.maxv = max(self.maxv, h);

        // per-byte histograms (pos 0 = least significant byte)
        for pos in 0..8 {
            let b = ((h >> (8 * pos)) & 0xFF) as usize;
            self.byte_hist[pos][b] += 1;
        }

        // FracMinHash threshold tests
        let hmax = u64::MAX;
        for (i, &s) in S_LIST.iter().enumerate() {
            let t = hmax / s; // floor((2^64-1)/s)
            if h < t {
                self.pass_counts[i] += 1;
            }
        }
    }

    fn merge(mut a: Stats, b: Stats) -> Stats {
        a.n += b.n;
        a.minv = min(a.minv, b.minv);
        a.maxv = max(a.maxv, b.maxv);

        for pos in 0..8 {
            for byte in 0..256 {
                a.byte_hist[pos][byte] += b.byte_hist[pos][byte];
            }
        }
        for i in 0..S_LIST.len() {
            a.pass_counts[i] += b.pass_counts[i];
        }
        a
    }
}

fn list_files(dir: &Path, pattern: &str) -> Result<Vec<PathBuf>> {
    let mut out = Vec::new();
    let g = dir.join(pattern);
    let gstr = g
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Bad path: {:?}", g))?
        .to_string();

    for entry in glob(&gstr)? {
        out.push(entry?);
    }
    out.sort();

    if out.is_empty() {
        return Err(anyhow::anyhow!("No files matched: {}", gstr));
    }
    Ok(out)
}

#[inline(always)]
fn parse_u64_le(bytes: &[u8]) -> u64 {
    // bytes must be length 8
    u64::from_le_bytes([
        bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
    ])
}

fn process_file(path: &Path) -> Result<Stats> {
    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(READER_CAP_BYTES, file);

    let mut st = Stats::new();

    let mut buf = vec![0u8; IO_CHUNK_BYTES];
    let mut tail: Vec<u8> = Vec::with_capacity(8);

    loop {
        let nread = reader.read(&mut buf)?;
        if nread == 0 {
            break;
        }

        // prepend leftover tail (0..7 bytes) from previous read
        if !tail.is_empty() {
            let mut combined = Vec::with_capacity(tail.len() + nread);
            combined.extend_from_slice(&tail);
            combined.extend_from_slice(&buf[..nread]);
            tail.clear();

            let usable = combined.len() / 8 * 8;
            for chunk in combined[..usable].chunks_exact(8) {
                let h = parse_u64_le(chunk);
                st.update(h);
            }

            let rem = combined.len() - usable;
            if rem > 0 {
                tail.extend_from_slice(&combined[usable..]);
            }
            continue;
        }

        // normal path: parse buf[..nread] in 8-byte chunks
        let usable = nread / 8 * 8;
        for chunk in buf[..usable].chunks_exact(8) {
            let h = parse_u64_le(chunk);
            st.update(h);
        }

        let rem = nread - usable;
        if rem > 0 {
            tail.extend_from_slice(&buf[usable..nread]);
        }
    }

    if !tail.is_empty() {
        return Err(anyhow::anyhow!(
            "Trailing {} bytes (file not multiple of 8?): {:?}",
            tail.len(),
            path
        ));
    }

    Ok(st)
}

fn compute(dir: &Path, pattern: &str, label: &str) -> Result<Stats> {
    let files = list_files(dir, pattern)?;
    eprintln!("{}: found {} files", label, files.len());

    let st = files
        .par_iter()
        .map(|p| process_file(p).unwrap_or_else(|e| panic!("Failed {:?}: {}", p, e)))
        .reduce(Stats::new, Stats::merge);

    eprintln!(
        "{}: N={} min=0x{:016x} max=0x{:016x}",
        label, st.n, st.minv, st.maxv
    );
    Ok(st)
}

#[derive(Serialize)]
struct SummaryRow {
    hash: String,
    n: String,
    min_hash: String,
    max_hash: String,
}

#[derive(Serialize)]
struct ThresholdRow {
    hash: String,
    scale: u64,
    threshold_hex: String,
    pass_count: String,
    pass_fraction: f64,
    expected_fraction: f64,
    rel_error: f64,
}

#[derive(Serialize)]
struct ByteHistRow {
    hash: String,
    pos: usize,
    byte: usize,
    count: u64,
}

fn write_csvs(out_dir: &Path, label: &str, st: &Stats) -> Result<()> {
    // summary
    let mut w = csv::Writer::from_path(out_dir.join(format!("{}_summary.csv", label)))?;
    w.serialize(SummaryRow {
        hash: label.to_string(),
        n: st.n.to_string(),
        min_hash: format!("0x{:016x}", st.minv),
        max_hash: format!("0x{:016x}", st.maxv),
    })?;
    w.flush()?;

    // thresholds
    let mut w = csv::Writer::from_path(out_dir.join(format!("{}_thresholds.csv", label)))?;
    let hmax = u64::MAX;
    for (i, &s) in S_LIST.iter().enumerate() {
        let t = hmax / s;
        let pass = st.pass_counts[i];
        let frac = (pass as f64) / (st.n as f64);
        let exp = 1.0 / (s as f64);
        let rel = (frac - exp) / exp;

        w.serialize(ThresholdRow {
            hash: label.to_string(),
            scale: s,
            threshold_hex: format!("0x{:016x}", t),
            pass_count: pass.to_string(),
            pass_fraction: frac,
            expected_fraction: exp,
            rel_error: rel,
        })?;
    }
    w.flush()?;

    // byte hist
    let mut w = csv::Writer::from_path(out_dir.join(format!("{}_bytehist.csv", label)))?;
    for pos in 0..8 {
        for byte in 0..256 {
            w.serialize(ByteHistRow {
                hash: label.to_string(),
                pos,
                byte,
                count: st.byte_hist[pos][byte],
            })?;
        }
    }
    w.flush()?;

    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 || args.len() > 3 {
        return Err(anyhow::anyhow!(
            "Usage: hash_dump_stats <base_dir> [output_dir]"
        ));
    }
    let base = Path::new(&args[1]);
    let out_dir = if args.len() == 3 {
        Path::new(&args[2])
    } else {
        base
    };

    // Input directories
    let jam_dir = base.join("sorted_jam");
    let mur_dir = base.join("sorted_murmur");

    // Input patterns
    let jam_pattern = "unsorted_bin_jam_hashes_*_*.bin";
    let mur_pattern = "unsorted_bin_murmur_hashes_*_*.bin";

    // Compute + write JamHash
    let jam = compute(&jam_dir, jam_pattern, "jamhash")?;
    write_csvs(out_dir, "jamhash", &jam)?;

    // Compute + write Murmur
    let mur = compute(&mur_dir, mur_pattern, "murmur")?;
    write_csvs(out_dir, "murmur", &mur)?;

    Ok(())
}
