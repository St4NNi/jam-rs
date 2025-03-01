use std::path::PathBuf;

use clap::{error::ErrorKind, CommandFactory, Parser};
use indicatif::ProgressIterator;
use jam_rs::{
    cli::{Cli, Commands},
    hash_functions::ahash,
    heed::HeedHandler,
};

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
        // Commands::Merge { inputs, output } => {
        //     match jam_rs::file_io::FileHandler::concat(inputs, output) {
        //         Ok(_) => {}
        //         Err(e) => {
        //             Cli::command().error(ErrorKind::ArgumentConflict, e).exit();
        //         }
        //     }
        // }
        Commands::Dist {
            input,
            database,
            output,
            cutoff,
        } => {
            let mut cmd = Cli::command();

            let input_files =
                jam_rs::file_io::FileHandler::test_and_collect_files(vec![input], false);
            let fs_input = match input_files {
                Ok(f) => f,
                Err(e) => {
                    cmd.error(ErrorKind::ArgumentConflict, e).exit();
                }
            };

            if database.len() == 1 {
                let mut lmdb = false;
                if let Some(first) = database.first() {
                    if first.is_file() {
                        if first.extension() == Some("mdb".as_ref()) {
                            lmdb = true;
                        }
                    }
                    if lmdb {
                        let mut lmdb_comparator = jam_rs::compare::LmdbComparator::new(
                            first.clone(),
                            args.threads.unwrap_or(1),
                            cutoff,
                            args.silent
                        )
                        .unwrap();

                        let mut input_sketch = Vec::new();

                        let iterator:Box<dyn Iterator<Item = PathBuf>> = if args.silent {
                            Box::new(fs_input.into_iter())
                        } else {
                            Box::new(fs_input.into_iter().progress())
                        };

                        for db_path in iterator {
                            // TODO: Remove hardcoded kmer sizes / settings / parse from db
                            match jam_rs::file_io::FileHandler::sketch_file(
                                &db_path,
                                lmdb_comparator.kmer_size,
                                lmdb_comparator.fscale,
                                None,
                                false,
                                jam_rs::hash_functions::Function::Small(&ahash),
                                jam_rs::cli::HashAlgorithms::Ahash,
                                false,
                            ) {
                                Ok(r) => {
                                    input_sketch.push(r);
                                }
                                Err(e) => {
                                    cmd.error(ErrorKind::ArgumentConflict, e).exit();
                                }
                            }
                        }

                        lmdb_comparator.set_signatures(input_sketch);

                        let mut result = match lmdb_comparator.compare() {
                            Ok(r) => r,
                            Err(e) => {
                                cmd.error(ErrorKind::ArgumentConflict, e).exit();
                            }
                        };

                        result.sort_by(|a, b| b.estimated_containment.total_cmp(&a.estimated_containment));

                        match output {
                            Some(o) => {
                                if let Err(e) =
                                    jam_rs::file_io::FileHandler::write_result(&result, o)
                                {
                                    cmd.error(ErrorKind::ArgumentConflict, e).exit();
                                }
                            }
                            None => {
                                for result in result {
                                    println!("{}", result);
                                }
                            }
                        }
                        return;
                    }
                }
            };

            let mut input_sketch = Vec::new();
            eprintln!("Reading input sketches");
            for db_path in fs_input {
                // TODO: Remove hardcoded kmer sizes / settings / parse from db
                match jam_rs::file_io::FileHandler::sketch_file(
                    &db_path,
                    21,
                    None,
                    None,
                    false,
                    jam_rs::hash_functions::Function::Small(&ahash),
                    jam_rs::cli::HashAlgorithms::Ahash,
                    false,
                ) {
                    Ok(r) => {
                        input_sketch.push(r);
                    }
                    Err(e) => {
                        cmd.error(ErrorKind::ArgumentConflict, e).exit();
                    }
                }
            }

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
        Commands::Stats { input, short } => {
            let mut cmd = Cli::command();

            let heed_handler = match HeedHandler::new_ro(input) {
                Ok(heed_handler) => heed_handler,
                Err(e) => {
                    cmd.error(ErrorKind::ArgumentConflict, e).exit();
                }
            };

            if short {
                match heed_handler.summarize_stats() {
                    Ok(_) => {}
                    Err(e) => {
                        cmd.error(ErrorKind::ArgumentConflict, e).exit();
                    }
                }
            } else {
                match heed_handler.detail_sigs() {
                    Ok(_) => {}
                    Err(e) => {
                        cmd.error(ErrorKind::ArgumentConflict, e).exit();
                    }
                }
            }
        }
    }
}
