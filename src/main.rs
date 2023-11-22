use clap::{error::ErrorKind, CommandFactory, Parser};
use jam_rs::cli::{Cli, Commands};

fn main() {
    let args = jam_rs::cli::Cli::parse();

    match args.command {
        Commands::Sketch { .. } => {
            let mut cmd = Cli::command();
            match jam_rs::file_io::FileHandler::sketch_files(args.command, args.threads) {
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
            stats,
            gc_lower,
            gc_upper,
        } => {
            let mut cmd = Cli::command();
            let database_files =
                jam_rs::file_io::FileHandler::test_and_collect_files(database, false);
            let fs = match database_files {
                Ok(f) => f,
                Err(e) => {
                    cmd.error(ErrorKind::ArgumentConflict, e).exit();
                }
            };

            let mut db_sketches = Vec::new();
            for db_path in fs {
                match jam_rs::file_io::FileHandler::read_signatures(&db_path) {
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

            let gc_bounds = match (gc_lower, gc_upper) {
                (Some(l), Some(u)) => Some((l, u)),
                (None, None) => None,
                _ => {
                    cmd.error(
                        ErrorKind::ArgumentConflict,
                        "Both gc_lower and gc_upper must be set",
                    )
                    .exit();
                }
            };

            let mut input_sketch = Vec::new();
            for db_path in fs_input {
                match jam_rs::file_io::FileHandler::read_signatures(&db_path) {
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
                stats,
                gc_bounds,
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
