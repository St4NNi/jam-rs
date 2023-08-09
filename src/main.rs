use clap::{error::ErrorKind, CommandFactory, Parser, Subcommand};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "jam")]
#[command(bin_name = "jam")]
#[command(version = "0.1.0")]
#[command(
    about = "Just another minhasher, obviously blazingly fast",
    long_about = "A heavily optimized minhash implementation that focuses less on accuracy and more on quick scans of large datasets."
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
    /// Number of threads to use
    #[arg(short, long, global = true, default_value = "1")]
    threads: Option<usize>,
    /// Overwrite output files
    #[arg(short, long, global = true, default_value = "false")]
    force: bool,
}

#[derive(Debug, Subcommand, Clone)]
enum Commands {
    /// Sketches one or more files and writes the result to an output file
    #[command(arg_required_else_help = true)]
    Sketch {
        /// Input file, directory or file with list of files to be hashed
        #[arg(short, long, required = true)]
        #[arg(value_parser = clap::value_parser!(std::path::PathBuf))]
        input: PathBuf,
        /// Output file
        #[arg(short, long, required = true)]
        #[arg(value_parser = clap::value_parser!(std::path::PathBuf))]
        output: PathBuf,
        /// kmer size all sketches to be compared must have the same size
        #[arg(short, long, default_value = "21")]
        kmer_size: u8,
        /// The estimated scaling factor to apply
        #[arg(short, long, default_value = "0.01")]
        scale: f32,
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
        #[arg(value_parser = clap::value_parser!(std::path::PathBuf))]
        output: Option<PathBuf>,
    },
}

fn main() {
    let args = Cli::parse();

    match &args.command {
        cmd @ Commands::Sketch {
            input,
            output,
            kmer_size,
            scale,
        } => {
            let mut cmd = Cli::command();
            cmd.error(ErrorKind::ArgumentConflict, "Unknown file type")
                .exit();
        }
        cmd @ Commands::Merge { inputs, output } => {
            dbg!(inputs, output);
        }
        cmd @ Commands::Compare {
            input,
            database,
            output,
        } => todo!(),
    }
}
