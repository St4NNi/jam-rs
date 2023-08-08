use crate::sketcher;
use anyhow::Result;
use needletail::parse_fastx_file;
use std::{
    fs::{self, File},
    io::{Read, Write},
};

pub struct FileHandler {}

impl FileHandler {
    pub fn sketch_file(
        &self,
        input: &str,
        output: &str,
        kmer_length: u8,
        scale: f32,
    ) -> Result<()> {
        let mut x = fs::metadata(input)?.len();
        if input.ends_with(".gz") {
            // Approximate the size of the uncompressed file
            x = x * 3;
        }
        let start = std::time::Instant::now();
        let kmer_num = (x as f64 * scale as f64) as u64;
        let mut sketcher = sketcher::Sketcher::new(kmer_length, kmer_num);
        let mut reader = parse_fastx_file(input)?;
        let mut counter = 0;
        while let Some(record) = reader.next() {
            let seqrec = record?;
            sketcher.process(&seqrec);
            counter += 1;
        }
        let mut o_file = File::create(output)?;
        o_file.write_all(&bincode::serialize(&sketcher.finalize())?)?;
        let elapsed = start.elapsed().as_millis();
        println!(
            "Processed {} with {} records, in {:?} seconds",
            input,
            counter,
            elapsed as f64 / 1000.0,
        );
        Ok(())
    }

    pub fn read_sketch(&self, input: &str) -> Result<sketcher::Sketch> {
        let mut i_file = File::open(input)?;
        let mut buffer = Vec::new();
        i_file.read_to_end(&mut buffer)?;
        Ok(bincode::deserialize(&buffer)?)
    }
}
