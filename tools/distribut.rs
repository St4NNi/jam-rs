use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use clap::Parser;
use crossbeam_channel::Sender;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use std::thread::JoinHandle;

#[derive(Parser)]
#[command(name = "hash_dedup")]
#[command(about = "Genomic hash deduplication tool for distributed processing")]
struct Args {
    /// Directory containing jam_hashes_*.bin files
    #[arg(long)]
    input_dir: String,

    /// Output directory for sorted bin files
    #[arg(long)]
    output_dir: String,

    /// Hash name
    #[arg(long, default_value = "jam_hashes")]
    hash_name: String,
}

fn process_hashes(
    hash_name: String,
    output_dir: String,
    input_channel: crossbeam_channel::Receiver<u64>,
) -> JoinHandle<Result<(), Box<dyn std::error::Error + Send + Sync>>> {
    std::thread::spawn(move || {
        let mut bins = [0u64; 65536];

        while let Ok(hash) = input_channel.recv() {
            bins[hash as usize % bins.len()] += 1;
        }

        // Write sorted output and detect duplicates directly from heap
        let output_filename = format!("distrib_bin_{}.bin", hash_name);
        let output_path = Path::new(&output_dir).join(&output_filename);
        let output_file = File::create(&output_path)?;
        let mut writer = BufWriter::new(output_file);

        for x in bins {
            writer.write_u64::<LittleEndian>(x)?;
            writer.write_all(b"\n")?;
        }

        writer.flush()?;

        println!("Finished writing to: {}", output_filename);
        println!("Max bin count: {}", bins.iter().max().unwrap_or(&0));
        println!("Min bin count: {}", bins.iter().min().unwrap_or(&0));
        println!(
            "Average bin count: {}",
            bins.iter().sum::<u64>() as f64 / bins.len() as f64
        );

        Ok(())
    })
}

const THREADS: usize = 16; // Number of threads to use for processing

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let mut process_handles: Vec<_> = vec![];
    let (tx, rx) = crossbeam_channel::bounded(10_000_000);

    process_handles.push(process_hashes(args.hash_name.clone(), args.output_dir, rx));

    let vec = (0..2000).collect::<Vec<_>>();
    let iter = vec.chunks(vec.len() / THREADS);

    println!("Processing files in {} threads...", iter.len());

    for chunk in iter.into_iter() {
        println!("Processing chunk: {:?}", chunk);
        let chunk = chunk.to_vec();
        let tx = tx.clone();
        let hash_name = args.hash_name.clone();
        let input_dir = args.input_dir.clone();
        process_handles.push(std::thread::spawn(move || {
            for file_ids in chunk {
                let filename = format!("{}/{}_{}.bin", input_dir, hash_name, file_ids);
                let filepath = Path::new(&filename);
                if !filepath.exists() {
                    eprintln!("File not found: {}", filepath.display());
                }
                process_file(&filepath, tx.clone())?;
            }
            Ok(())
        }));
    }

    drop(tx); // Close the sender to signal no more messages

    for handle in process_handles {
        handle
            .join()
            .expect("Thread panicked")
            .expect("Failed to process file");
    }

    Ok(())
}

fn process_file(
    filepath: &Path,
    sender: Sender<u64>,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let file = File::open(filepath)?;
    // Read 100 MB at a time for efficiency
    let mut reader = BufReader::with_capacity(1024 * 1024 * 500, file);

    loop {
        match reader.read_u64::<LittleEndian>() {
            Ok(hash) => {
                sender
                    .send(hash)
                    .expect("Failed to send hash to processing thread");
            }
            Err(ref e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e.into()),
        }
    }

    println!("Finished processing file: {}", filepath.display());

    Ok(())
}
