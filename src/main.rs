use clap::{error::ErrorKind, CommandFactory, Parser, Subcommand, ValueEnum};
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

#[derive(ValueEnum, Debug, Clone)]
pub enum OutputFormats {
    Bin,
    Sourmash,
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
        #[arg(short = 'k', long = "kmer-size", default_value = "21")]
        kmer_size: u8,
        /// Scale the hash spaces to a minimum fraction of the maximum hash value
        #[arg(long, default_value = "0")]
        fscale: u64,
        /// Scale the hash spaces to a minimum fraction of all k-mers
        #[arg(long, default_value = "1000")]
        kscale: u64,
        /// Minimum number of k-mers (per record) to be hashed
        #[arg(long, default_value = "0")]
        nmin: u64,
        /// Maximum number of k-mers (per record) to be hashed
        #[arg(long, default_value = "0")]
        nmax: u64,
        /// Change to other output formats
        #[arg(long, default_value = "bin")]
        format: OutputFormats,
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

fn main() {
    let args = Cli::parse();

    match args.command {
        Commands::Sketch {
            input,
            output,
            kmer_size,
            fscale: _,
            kscale: _,
            nmin: _,
            nmax: _,
            format: _,
        } => {
            let mut cmd = Cli::command();

            let files = jam_rs::file_io::FileHandler::test_and_collect_files(vec![input], true);
            let fs = match files {
                Ok(f) => f,
                Err(e) => {
                    cmd.error(ErrorKind::ArgumentConflict, e).exit();
                }
            };
            match jam_rs::file_io::FileHandler::sketch_files(
                fs,
                output,
                kmer_size,
                0.0, // FIXME: new_scales
                args.threads.unwrap(),
            ) {
                Ok(_) => {}
                Err(e) => {
                    cmd.error(ErrorKind::ArgumentConflict, e).exit();
                }
            }
        }
        Commands::Merge { inputs, output } => {
            match jam_rs::file_io::FileHandler::concat(inputs, output) {
                Ok(_) => {}
                Err(e) => {
                    Cli::command().error(ErrorKind::ArgumentConflict, e).exit();
                }
            }
        }
        Commands::Dist {
            input,
            database,
            output,
            cutoff,
        } => {
            let mut cmd = Cli::command();
            let database_files =
                jam_rs::file_io::FileHandler::test_and_collect_files(vec![database], false);
            let fs = match database_files {
                Ok(f) => f,
                Err(e) => {
                    cmd.error(ErrorKind::ArgumentConflict, e).exit();
                }
            };

            let mut db_sketches = Vec::new();
            for db_path in fs {
                match jam_rs::file_io::FileHandler::read_sketches(&db_path) {
                    Ok(r) => {
                        db_sketches.extend(r);
                    }
                    Err(e) => {
                        cmd.error(ErrorKind::ArgumentConflict, e).exit();
                    }
                }
            }

            let input_files =
                jam_rs::file_io::FileHandler::test_and_collect_files(vec![input], false);
            let fs_input = match input_files {
                Ok(f) => f,
                Err(e) => {
                    cmd.error(ErrorKind::ArgumentConflict, e).exit();
                }
            };

            let mut input_sketch = Vec::new();
            for db_path in fs_input {
                match jam_rs::file_io::FileHandler::read_sketches(&db_path) {
                    Ok(r) => {
                        input_sketch.extend(r);
                    }
                    Err(e) => {
                        cmd.error(ErrorKind::ArgumentConflict, e).exit();
                    }
                }
            }

            match jam_rs::compare::MultiComp::new(
                input_sketch,
                db_sketches,
                args.threads.unwrap(),
                cutoff,
            ) {
                Ok(mut mc) => {
                    if let Err(e) = mc.compare() {
                        cmd.error(ErrorKind::ArgumentConflict, e).exit();
                    }
                    let result = mc.finalize();
                    match output {
                        Some(o) => {
                            if let Err(e) = jam_rs::file_io::FileHandler::write_result(&result, o) {
                                cmd.error(ErrorKind::ArgumentConflict, e).exit();
                            }
                        }
                        None => {
                            for result in result {
                                println!("{}", result);
                            }
                        }
                    }
                }
                Err(e) => {
                    cmd.error(ErrorKind::ArgumentConflict, e).exit();
                }
            }
        }
    }
}
