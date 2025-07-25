// Update to cli.rs - enhanced version with input path expansion
use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "jam")]
#[command(bin_name = "jam")]
#[command(version = "1.0.0")]
#[command(
    about = "Just another (genomic) minhasher (jam), obviously blazingly fast",
    long_about = "An optimized minhash implementation that focuses on quick scans for small sequences in large datasets."
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
    /// Number of threads to use
    #[arg(short, long, global = true, default_value = "1")]
    pub threads: Option<usize>,
    /// Overwrite output files
    #[arg(short, long, global = true, default_value = "false")]
    pub force: bool,
    /// Silent mode, no (additional) output to stdout
    /// Only errors and output files will be printed
    #[arg(short, long, global = true, default_value = "false")]
    pub silent: bool,
}

#[derive(Debug, Subcommand, Clone)]
pub enum Commands {
    /// Sketch one or more files and write the result to an output file
    #[command(arg_required_else_help = true)]
    Sketch {
        /// Input file(s), directories, or file with list of files to be hashed
        #[arg(value_parser = clap::value_parser!(std::path::PathBuf))]
        input: Vec<PathBuf>,
        /// Output file (.lmdb will be appended if not present)
        #[arg(short, long)]
        #[arg(value_parser = clap::value_parser!(std::path::PathBuf))]
        output: PathBuf,
        /// K-mer size, all sketches must have the same size to be compared and below 32
        #[arg(short = 'k', long = "kmer-size", default_value = "21")]
        kmer_size: u8,
        /// Scale the hash space to a minimum fraction of the maximum hash value (FracMinHash)
        #[arg(long)]
        fscale: Option<u64>,
        /// Maximum number of k-mers (per record) to be hashed, top cut-off
        #[arg(long)]
        nmax: Option<u64>,
        /// Complexity cut-off, only hash sequences with complexity above this value
        /// This is created via shannon entropy
        #[arg(long, default_value = "0.0")]
        complexity: f64,
        /// Create a separate sketch for each sequence record
        /// Will increase the size of the output file
        #[arg(long)]
        singleton: bool,
        /// Custom temporary directory for intermediate files during sorting
        #[arg(long)]
        temp_dir: Option<PathBuf>,
    },

    /// Estimate containment of a (small) sketch against a subset of one or more sketches as database.
    /// Requires all sketches to have the same kmer size
    #[command(arg_required_else_help = true)]
    Dist {
        /// Input sketch or raw sequence file
        #[arg(short, long)]
        input: PathBuf,
        /// Database sketch (.lmdb file)
        #[arg(short, long)]
        database: PathBuf,
        /// Output to file instead of stdout
        #[arg(short, long)]
        #[arg(value_parser = clap::value_parser!(std::path::PathBuf))]
        output: Option<PathBuf>,
        /// Cut-off value for similarity/containment
        #[arg(short, long, default_value = "0.0")]
        cutoff: f64,
        /// Singleton mode, process each query sequence separately
        #[arg(short, long, default_value = "false")]
        singleton: bool,
    },

    /// Display statistics about an LMDB database
    #[command(arg_required_else_help = true)]
    Stats {
        /// Input LMDB database
        #[arg(short, long)]
        input: PathBuf,
        /// Short summary only
        #[arg(short, long)]
        short: bool,
    },
}
