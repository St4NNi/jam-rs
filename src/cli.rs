use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "jam")]
#[command(bin_name = "jam")]
#[command(version = "0.1.0")]
#[command(
    about = "Just another minhasher, obviously blazingly fast",
    long_about = "A heavily optimized minhash implementation that focuses less on accuracy and more on quick scans of large datasets."
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
}

#[derive(ValueEnum, Debug, Clone)]
pub enum OutputFormats {
    Bin,
    // Sourmash compatible json
    Sourmash,
}

#[derive(ValueEnum, Debug, Clone)]
pub enum HashAlgorithms {
    Default, // AHash < 32 | Xxhash >= 32
    AHash,
    Xxhash,
    Murmur3,
}

#[derive(Debug, Subcommand, Clone)]
pub enum Commands {
    /// Sketches one or more files and writes the result to an output file
    #[command(arg_required_else_help = true)]
    Sketch {
        /// Input file(s), one directory or one file with list of files to be hashed
        #[arg(value_parser = clap::value_parser!(std::path::PathBuf))]
        input: Vec<PathBuf>,
        /// Output file
        #[arg(short, long)]
        #[arg(value_parser = clap::value_parser!(std::path::PathBuf))]
        output: Option<PathBuf>,
        /// kmer size all sketches to be compared must have the same size
        #[arg(short = 'k', long = "kmer-size", default_value = "21")]
        kmer_size: u8,
        /// Scale the hash spaces to a minimum fraction of the maximum hash value
        #[arg(long, default_value = "0")]
        fscale: u64,
        /// Scale the hash spaces to a minimum fraction of all k-mers
        #[arg(long, default_value = "1000")]
        kscale: u64,
        /// Minimum number of k-mers (per record) to be hashed
        #[arg(long)]
        nmin: Option<u64>,
        /// Maximum number of k-mers (per record) to be hashed
        #[arg(long)]
        nmax: Option<u64>,
        /// Change to other output formats
        #[arg(long, default_value = "bin")]
        format: OutputFormats,
        /// Change the hashing algorithm
        #[arg(long, default_value = "default")]
        algorithm: HashAlgorithms,
        /// Create a separate sketch for each sequence record
        #[arg(long)]
        singleton: bool,
    },
    /// Merge multiple input sketches into a single sketch
    #[command(arg_required_else_help = true)]
    Merge {
        /// One or more input sketches
        #[arg(value_parser = clap::value_parser!(std::path::PathBuf))]
        inputs: Vec<PathBuf>,
        /// Output file
        #[arg(short, long, required = true)]
        #[arg(value_parser = clap::value_parser!(std::path::PathBuf))]
        output: PathBuf,
    },
    /// Calculate distance of a (small) sketch against one or more sketches as database.
    /// Requires all sketches to have the same kmer size
    #[command(arg_required_else_help = true)]
    Dist {
        /// Input sketch or raw file
        #[arg(short, long)]
        input: PathBuf,
        /// Database sketch(es)
        #[arg(short, long)]
        database: PathBuf,
        /// Output to file instead of stdout
        #[arg(short, long)]
        #[arg(value_parser = clap::value_parser!(std::path::PathBuf))]
        output: Option<PathBuf>,
        /// Cut-off value for similarity
        #[arg(short, long, default_value = "0.0")]
        cutoff: f64,
    },
}