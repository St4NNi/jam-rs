use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "jam")]
#[command(bin_name = "jam")]
#[command(version = "0.9.10")]
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
    /// Maximum memory usage in bytes in GB
    #[arg(short, long, global = true, default_value = "2")]
    pub memory: Option<usize>,
}

#[derive(Debug, Subcommand, Clone)]
pub enum Commands {
    /// Sketch one or more files and write the result to an output file
    #[command(arg_required_else_help = true)]
    Sketch {
        /// Input file(s), directories, or file with list of files to be hashed
        #[arg(value_parser = clap::value_parser!(std::path::PathBuf))]
        input: Vec<PathBuf>,
        /// Output file (.jam format)
        #[arg(short, long)]
        #[arg(value_parser = clap::value_parser!(std::path::PathBuf))]
        output: PathBuf,
        /// K-mer size, all sketches must have the same size to be compared and below 32
        #[arg(short = 'k', long = "kmer-size", default_value = "21")]
        kmer_size: u8,
        /// Scale the hash space to a minimum fraction of the maximum hash value (FracMinHash)
        #[arg(long)]
        fscale: Option<u64>,
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
        /// Path to a bias table file (.bias) for hash-based filtering
        #[arg(long)]
        bias_table: Option<PathBuf>,
    },

    /// Estimate containment of a query sequence against a sketch database.
    /// Requires all sketches to have the same kmer size
    #[command(arg_required_else_help = true)]
    Dist {
        /// Input FASTA/FASTQ file to query
        #[arg(short, long)]
        input: PathBuf,
        /// Database sketch (.jam file)
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
        #[arg(long, default_value = "false")]
        singleton: bool,
    },

    /// Build and analyze hash bias tables for filtering
    #[command(arg_required_else_help = true)]
    Bias {
        #[command(subcommand)]
        command: BiasCommands,
    },

    /// Display statistics about a JAM database
    #[command(arg_required_else_help = true)]
    Stats {
        /// Input JAM database (.jam file)
        #[arg(short, long)]
        input: PathBuf,
        /// Short summary only
        #[arg(long)]
        short: bool,
        /// Include the full entry statistics
        #[arg(long)]
        full: bool,
    },
}

#[derive(Debug, Subcommand, Clone)]
pub enum BiasCommands {
    /// Count hashes from reference FASTA files into a raw hash counts file
    #[command(arg_required_else_help = true)]
    Create {
        /// Input FASTA file(s) or glob pattern
        #[arg(required = true)]
        input: Vec<PathBuf>,
        /// Output raw hash counts file (.braw)
        #[arg(short, long)]
        output: PathBuf,
        /// K-mer size (must match sketch k-mer size)
        #[arg(short = 'k', long = "kmer-size", default_value = "21")]
        kmer_size: u8,
        /// FracMinHash scale (must match sketch fscale)
        #[arg(long, default_value = "1000")]
        fscale: u64,
        /// Count-Min Sketch width (columns, power of 2 recommended)
        #[arg(long, default_value = "1048576")]
        cms_width: usize,
        /// Count-Min Sketch depth (number of hash functions)
        #[arg(long, default_value = "5")]
        cms_depth: usize,
    },

    /// Combine positive and negative raw counts into a trained bias table
    #[command(arg_required_else_help = true)]
    Combine {
        /// Positive raw hash counts (.braw) - sequences to enrich
        positive: PathBuf,
        /// Negative raw hash counts (.braw) - sequences to deplete
        negative: PathBuf,
        /// Output combined bias table file (.bias)
        #[arg(short, long)]
        output: PathBuf,
        /// Smoothing parameter for log-ratio computation
        #[arg(long, default_value = "1.0")]
        alpha: f32,
        /// Target selectivity ratio (positive retention / negative retention).
        /// E.g., 10.0 means positive hashes are retained 10x more than negative.
        #[arg(long, default_value = "10.0")]
        selectivity: f32,
    },

    /// Display statistics for a bias table or raw hash counts
    #[command(arg_required_else_help = true)]
    Stats {
        /// Input file (.braw or .bias)
        input: PathBuf,
        /// Optional second raw hash counts (.braw) for comparison
        #[arg(long)]
        compare: Option<PathBuf>,
        /// Output JSON report to file instead of stderr
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}
