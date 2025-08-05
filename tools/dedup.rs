use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use clap::Parser;
use kanal::Sender;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::Arc;
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

    /// Bindex
    #[arg(long)]
    bindex: usize,
}

// 16 Bindex â 128 Threads == 2048 files
fn process_output(
    bin_id: usize,
    bindex: usize,
    hash_name: String,
    output_dir: String,
    input_channel: kanal::Receiver<u64>,
) -> JoinHandle<Result<[u64; 512], Box<dyn std::error::Error + Send + Sync>>> {
    std::thread::spawn(move || {
        let output_file = format!("unsorted_bin_{}_{}_{}.bin", hash_name, bindex, bin_id);
        let output_path = Path::new(&output_dir).join(&output_file);
        let output_file = File::create(&output_path)?;
        let mut writer = BufWriter::with_capacity(1024 * 1024 * 256, output_file); // 1 GB Buffersize
        let mut category = [0u64; 512];

        while let Ok(hash) = input_channel.recv() {
            writer.write_u64::<LittleEndian>(hash)?;
            category[hash as usize % category.len()] += 1;
        }
        writer.flush()?;

        println!("Completed bin {}:", bin_id);

        Ok(category)
    })
}

const THREADS: usize = 64; // Number of threads to use for processing
const OUTPUT_THREADS: usize = 128; // Number of threads for output processing

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let mut senders = vec![];
    let mut process_handles: Vec<
        JoinHandle<Result<[u64; 512], Box<dyn std::error::Error + Send + Sync>>>,
    > = vec![];
    for bin_id in 0..OUTPUT_THREADS {
        let (sender, receiver) = kanal::bounded(10_000_000);
        senders.push(sender);
        process_handles.push(process_output(
            bin_id,
            args.bindex,
            args.hash_name.clone(),
            args.output_dir.clone(),
            receiver,
        ));
    }

    let arc_senders = std::sync::Arc::new(senders);

    let vec = (0..2000).collect::<Vec<_>>();
    let iter = vec.chunks(vec.len() / THREADS);

    let mut read_threads: Vec<
        std::thread::JoinHandle<Result<(), Box<dyn std::error::Error + Send + Sync>>>,
    > = vec![];

    for chunk in iter.into_iter() {
        let chunk = chunk.to_vec();
        let hash_name = args.hash_name.clone();
        let input_dir = args.input_dir.clone();
        let bindex = args.bindex;
        let senders = arc_senders.clone();
        read_threads.push(std::thread::spawn(move || {
            for file_ids in chunk {
                let senders = senders.clone();
                let filename = format!("{}/{}_{}.bin", input_dir, hash_name, file_ids);
                let filepath = Path::new(&filename);
                if !filepath.exists() {
                    eprintln!("File not found: {}", filepath.display());
                }
                process_file(&filepath, bindex, senders)?;
            }
            Ok(())
        }));
    }

    for handle in read_threads {
        handle
            .join()
            .expect("Thread panicked")
            .expect("Failed to process file");
    }

    arc_senders.iter().for_each(|sender| {
        while !sender.is_empty() {}
        sender.close().unwrap();
    });

    let mut bins = vec![0u64; 512 * OUTPUT_THREADS];

    for (bin_id, handle) in process_handles.into_iter().enumerate() {
        let category = handle
            .join()
            .expect("Thread panicked")
            .expect("Failed to process output");
        for (i, count) in category.iter().enumerate() {
            bins[bin_id * 512 + i] += count;
        }
    }

    // Write the bin counts to a file
    let output_file = format!(
        "{}/{}_bindex_{}_counts.txt",
        args.output_dir, args.hash_name, args.bindex
    );
    let output_path = Path::new(&output_file);
    let output_file = File::create(&output_path)?;
    let mut writer = BufWriter::new(output_file);

    for bin in bins.iter() {
        writer.write_u64::<LittleEndian>(*bin)?;
    }

    writer.flush()?;

    println!("Bin counts written to: {}", output_path.display());
    println!("Total bins: {}", bins.len());
    println!("Total hashes processed: {}", bins.iter().sum::<u64>());
    println!("Max bin count: {}", bins.iter().max().unwrap_or(&0));
    println!("Min bin count: {}", bins.iter().min().unwrap_or(&0));
    println!(
        "Average bin count: {}",
        bins.iter().sum::<u64>() as f64 / bins.len() as f64
    );

    Ok(())
}

fn process_file(
    filepath: &Path,
    our_bindex: usize,
    senders: Arc<Vec<Sender<u64>>>,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let file = File::open(filepath)?;
    let mut reader = BufReader::with_capacity(1024 * 1024 * 256, file);

    loop {
        match reader.read_u64::<LittleEndian>() {
            Ok(hash) => {
                // First nibble is the "bindex"
                let bindex = (hash >> 60) as usize;

                if bindex != our_bindex {
                    continue; // Skip hashes that do not match the bindex
                }
                let index = ((hash >> 53) & 0x7F) as usize;
                senders[index]
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
