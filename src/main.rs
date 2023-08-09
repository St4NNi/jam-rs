use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "jam")]
#[command(bin_name = "jam")]
#[command(version = "0.1.0")]
#[command(
    about = "Just another minhasher, obviously blazingly fast",
    long_about = "This is an heavily optimized minhash implementation that focuses less on accuracy and more on quick scans of large datasets."
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
    /// Number of threads to use
    #[arg(short, long, global = true)]
    threads: Option<usize>,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Sketches one or more files and writes the result to a file
    #[command(arg_required_else_help = true)]
    Sketch {
        /// Input file, directory or file with list of files to be hashed
        #[arg(short, long, required = true)]
        input: PathBuf,
        /// Output file
        #[arg(short, long, required = true)]
        output: PathBuf,
        /// The kmer size to use all sketches must have the same kmer size
        #[arg(short, long, default_value = "21")]
        kmer_size: u8,
        /// The estimated scaling factor to apply
        #[arg(short, long, default_value = "0.1")]
        scale: f32,
    },
    /// Merge multiple input sketches into a single sketch
    #[command(arg_required_else_help = true)]
    Merge {
        /// One or more input sketches
        inputs: Vec<PathBuf>,
        /// Output file
        #[arg(short, long, required = true)]
        output: PathBuf,
    },
    /// Compare a raw file or sketch against one or more sketches as database
    /// Requires all sketches to have the same kmer size
    #[command(arg_required_else_help = true)]
    Compare {
        /// Input sketch or raw file
        #[arg(short, long)]
        input: PathBuf,
        /// Database sketch(es)
        #[arg(short, long)]
        database: PathBuf,
        /// Output to file instead of stdout
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}

fn main() {
    let args = Cli::parse();

    match args.command {
        Commands::Sketch {
            input,
            output,
            kmer_size,
            scale,
        } => todo!(),
        Commands::Merge { inputs, output } => {
            dbg!(inputs, output);
        }
        Commands::Compare {
            input,
            database,
            output,
        } => todo!(),
    }
}
