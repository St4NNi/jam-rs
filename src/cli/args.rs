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
        /// Path to a bias table file (.bias) for compositional filtering
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

    /// Build and analyze bias tables for compositional filtering
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
    /// Create a raw hexamer frequency table from a FASTA file
    #[command(arg_required_else_help = true)]
    Create {
        /// Input FASTA file
        input: PathBuf,
        /// Output raw bias table file (.braw)
        #[arg(short, long)]
        output: PathBuf,
    },

    /// Combine two raw tables into a final bias table
    #[command(arg_required_else_help = true)]
    Combine {
        /// Positive raw bias table (.braw)
        positive: PathBuf,
        /// Negative raw bias table (.braw)
        negative: PathBuf,
        /// Output combined bias table file (.bias)
        #[arg(short, long)]
        output: PathBuf,
        /// Score threshold (0.0-1.0). Higher = stricter filtering.
        /// 0.5 = neutral, 0.7 = keep hexamers 70%+ likely from positive set
        #[arg(long, default_value = "0.5")]
        threshold: f32,
    },

    /// Display statistics for a bias table
    #[command(arg_required_else_help = true)]
    Stats {
        /// Input bias table file (.braw or .bias)
        input: PathBuf,
        /// Output JSON report to file instead of stderr
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Compare two raw bias tables
    #[command(arg_required_else_help = true)]
    Compare {
        /// First raw bias table (.braw) - treated as positive
        positive: PathBuf,
        /// Second raw bias table (.braw) - treated as negative
        negative: PathBuf,
        /// Output JSON report to file instead of stderr
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Threshold for summary statistics
        #[arg(long, default_value = "0.5")]
        threshold: f32,
    },
}
