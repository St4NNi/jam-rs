use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "jam")]
#[command(bin_name = "jam")]
#[command(version)]
#[command(
    about = "High-speed reference-guided fragment search in metagenomic assemblies",
    long_about = "Find candidate sequence fragments from one query element across metagenomic assemblies. Sequence evidence requires independent biological interpretation."
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
    /// Internal memory target in GiB; process RSS may include additional overhead
    #[arg(
        short = 'm',
        long = "memory-target",
        visible_alias = "memory",
        global = true,
        default_value = "2"
    )]
    pub memory_target: Option<usize>,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum TraceSensitivityArg {
    Fast,
    Balanced,
    Sensitive,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
pub enum QueryKindArg {
    Plasmid,
    Phage,
    Other,
    #[default]
    Unknown,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
pub enum TopologyArg {
    Linear,
    Circular,
    #[default]
    Auto,
    Unknown,
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
        /// Scale the hash space to a minimum fraction of the maximum hash value (FracMinHash, default 100)
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

    /// Screen assembly contigs against a reference catalog and aggregate evidence
    #[command(arg_required_else_help = true)]
    Screen {
        /// Assembly FASTA/FASTQ input (compressed formats supported)
        #[arg(short, long)]
        input: PathBuf,
        /// Reference catalog database (.jam)
        #[arg(short, long)]
        database: PathBuf,
        /// Contig-level candidate TSV
        #[arg(short, long)]
        output: PathBuf,
        /// Assembly-level reference summary TSV
        #[arg(long)]
        summary: PathBuf,
        /// Optional run metadata JSON
        #[arg(long)]
        metadata: Option<PathBuf>,
        /// Assembly identifier (defaults to the input file name)
        #[arg(long)]
        assembly_name: Option<String>,
        /// Minimum distinct shared hashes per contig-reference pair
        #[arg(long, default_value = "3")]
        min_shared: u32,
        /// Minimum query containment (retained-query containment in bias mode)
        #[arg(long, default_value = "0.0")]
        min_query_containment: f64,
        /// Minimum reference containment (retained-reference containment in bias mode)
        #[arg(long, default_value = "0.0")]
        min_reference_containment: f64,
        /// Maximum reported references per contig after deterministic ranking
        #[arg(long, default_value = "10")]
        top_per_contig: usize,
        /// Maximum references in the assembly summary after deterministic ranking
        #[arg(long, default_value = "100")]
        top_references: usize,
    },

    /// Build a JMA v1 positional archive for trace searches
    #[command(arg_required_else_help = true)]
    Archive {
        /// Metagenomic assembly FASTA/FASTQ input
        #[arg(short, long)]
        input: PathBuf,
        /// Output JMA archive
        #[arg(short, long)]
        output: PathBuf,
        /// Maximum decoded bases per packed sequence block
        #[arg(long, default_value = "1048576")]
        block_bases: usize,
        /// Primary k=31 FracMinHash scale
        #[arg(long, default_value = "200")]
        primary_scale: u64,
        /// Rescue k=21 FracMinHash scale
        #[arg(long, default_value = "500")]
        rescue_scale: u64,
        /// Omit the k=21 rescue seed section
        #[arg(long)]
        no_rescue: bool,
        /// Optional Shannon-entropy threshold in bits per base
        #[arg(long)]
        complexity: Option<f64>,
    },

    /// Trace one query sequence element across candidate metagenomic assemblies
    #[command(arg_required_else_help = true)]
    Trace {
        /// FASTA/FASTQ containing exactly one query sequence
        #[arg(
            short = 'q',
            long,
            required_unless_present = "plasmid",
            conflicts_with = "plasmid"
        )]
        query: Option<PathBuf>,
        /// Compatibility alias for --query; implies --query-kind plasmid
        #[arg(short = 'p', long, conflicts_with = "query")]
        plasmid: Option<PathBuf>,
        /// Declared query class; this is metadata, not a classifier result
        #[arg(long, value_enum, default_value = "unknown")]
        query_kind: QueryKindArg,
        /// Requested query-coordinate handling; wrap does not prove circular biology
        #[arg(long, value_enum, default_value = "auto")]
        topology: TopologyArg,
        /// Existing metagenome candidate index (.jam)
        #[arg(short, long)]
        database: String,
        /// TSV or JSON catalog mapping database sample IDs to JMA/raw resources
        #[arg(short = 'c', long = "metagenomes", visible_alias = "catalog")]
        metagenomes: PathBuf,
        /// JSONL output; .zst or .zstd enables Zstandard compression
        #[arg(short, long)]
        output: PathBuf,
        /// Upload the finalized local output once to an HTTP(S) or S3 object
        #[arg(long)]
        upload_to: Option<String>,
        /// Override the query FASTA record identifier
        #[arg(long = "query-id", visible_alias = "plasmid-id")]
        query_id: Option<String>,
        /// Bounded execution profile, not a calibrated biological sensitivity
        #[arg(long, value_enum, default_value = "balanced")]
        sensitivity: TraceSensitivityArg,
        /// Minimum shared hashes for metagenome candidate retrieval
        #[arg(long, default_value = "3")]
        min_shared: u32,
        /// Minimum retained query containment
        #[arg(
            long = "min-query-containment",
            visible_alias = "min-plasmid-containment",
            default_value = "0.0"
        )]
        min_query_containment: f64,
        /// Minimum retained metagenome/reference containment
        #[arg(long, default_value = "0.0")]
        min_metagenome_containment: f64,
        /// Override the profile's deterministic candidate limit
        #[arg(long)]
        top_candidates: Option<usize>,
        /// Maximum retained alignments per candidate metagenome
        #[arg(long, default_value = "256")]
        max_alignments: usize,
        /// Maximum candidate resource tasks in flight
        #[arg(long)]
        io_concurrency: Option<usize>,
        /// Minimum new-base margin required to select one auto topology model
        #[arg(long)]
        topology_margin_bases: Option<u64>,
        /// Directory for identity-checked remote .jam materialization
        #[arg(long)]
        cache_dir: Option<PathBuf>,
        /// In-memory remote range-cache block size in bytes
        #[arg(long, default_value = "1048576")]
        cache_block_bytes: u64,
        /// Per-request timeout in seconds
        #[arg(long, default_value = "60")]
        request_timeout_seconds: u64,
        /// Number of retries after a failed remote request or database download
        #[arg(long, default_value = "3")]
        max_retries: u32,
        /// Reject servers that require a complete-object fallback for range reads
        #[arg(long)]
        no_full_download_fallback: bool,
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
        /// Minimum query containment (bias-weighted query containment in bias mode)
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
        /// Emit stable JSON statistics
        #[arg(long, conflicts_with_all = ["short", "full"])]
        json: bool,
    },
}

#[derive(Debug, Subcommand, Clone)]
pub enum BiasCommands {
    /// Create a bias table from positive (target) and negative (background) FASTA files.
    /// Target signal is always subtracted from background before computing bias weights.
    /// Note: bias table format is experimental and may change between versions.
    #[command(arg_required_else_help = true)]
    Create {
        /// Positive (target) FASTA file(s) - sequences to enrich for
        #[arg(long, required = true, num_args = 1..)]
        positive: Vec<PathBuf>,
        /// Negative (background) FASTA file(s) - sequences to deplete.
        /// Target signal is subtracted from background automatically.
        #[arg(long, required = true, num_args = 1..)]
        negative: Vec<PathBuf>,
        /// Output bias table file (.bias)
        #[arg(short, long)]
        output: PathBuf,
        /// K-mer size (must match sketch k-mer size)
        #[arg(short = 'k', long = "kmer-size", default_value = "21")]
        kmer_size: u8,
        /// FracMinHash scale (must match sketch fscale)
        #[arg(long, default_value = "100")]
        fscale: u64,
        /// Count-Min Sketch width (columns, power of 2 recommended)
        #[arg(long, default_value = "1048576")]
        cms_width: usize,
        /// Count-Min Sketch depth (number of hash functions)
        #[arg(long, default_value = "5")]
        cms_depth: usize,
        /// Smoothing parameter for log-ratio computation
        #[arg(long, default_value = "1.0")]
        alpha: f32,
        /// Target effective fscale under the empirical CMS weight prior.
        /// Enables enrichment LUT filtering. Must be >= base fscale.
        #[arg(long)]
        target_fscale: Option<u64>,
        /// Maximum fscale per bucket (minimum retention). Higher = more suppression.
        /// Use "drop" to allow complete suppression of negatively biased hashes.
        /// Required when --target-fscale is set.
        #[arg(long)]
        max_fscale: Option<String>,
        /// Fixed fscale for weight=0 (unseen/balanced k-mers).
        /// Default: same as --target-fscale.
        #[arg(long)]
        unseen_fscale: Option<u64>,
        /// Number of threads to use for bias sketching
        #[arg(long)]
        threads: Option<usize>,
        /// Reject calibration below this positive-set retention fraction
        #[arg(long, default_value = "0.25")]
        min_positive_retention: f32,
    },

    /// Display statistics for a bias table (.bias file)
    #[command(arg_required_else_help = true)]
    Stats {
        /// Input bias table file (.bias)
        input: PathBuf,
        /// Output JSON report to file instead of stderr
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}
