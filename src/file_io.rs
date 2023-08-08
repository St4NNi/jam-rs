use crate::sketcher;
use anyhow::Result;
use needletail::parse_fastx_file;
use std::{fs, io::BufReader};

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
            x *= 3;
        }
        let start = std::time::Instant::now();
        let kmer_num = (x as f64 * scale as f64) as u64;
        let mut sketcher = sketcher::Sketcher::new(kmer_length, kmer_num, input.to_string());
        let mut reader = parse_fastx_file(input)?;
        let mut counter = 0;
        while let Some(record) = reader.next() {
            let seqrec = record?;
            sketcher.process(&seqrec);
            counter += 1;
        }
        let o_file = std::fs::File::create(output)?;
        let bufwriter = std::io::BufWriter::new(o_file);

        bincode::serialize_into(bufwriter, &sketcher.finalize())?;
        let elapsed = start.elapsed().as_millis();
        println!(
            "Processed {} with {} records, in {:?} seconds",
            input,
            counter,
            elapsed as f64 / 1000.0,
        );
        Ok(())
    }

    pub fn read_sketches(&self, input: &str) -> Result<Vec<sketcher::Sketch>> {
        let mut vec: Vec<sketcher::Sketch> = vec![];

        let mut reader = BufReader::new(std::fs::File::open(input)?);

        while let Ok(result) = bincode::deserialize_from(&mut reader) {
            vec.push(result);
        }

        Ok(vec)
    }
}
