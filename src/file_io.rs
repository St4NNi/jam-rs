use crate::cli::Commands;
use crate::cli::HashAlgorithms;
use crate::cli::OutputFormats;
use crate::compare::CompareResult;
use crate::hash_functions::Function;
use crate::signature::Signature;
use crate::sketch::Sketch;
use crate::sketcher;
use anyhow::anyhow;
use anyhow::Result;
use byteorder::BigEndian;
use heed::types::SerdeBincode;
use heed::types::U32;
use heed::types::U64;
use heed::DatabaseFlags;
use heed::EnvFlags;
use heed::PutFlags;
use indicatif::MultiProgress;
use indicatif::ParallelProgressIterator;
use indicatif::ProgressBar;
use needletail::parse_fastx_file;
use rayon::prelude::IntoParallelRefIterator;
use rayon::prelude::ParallelIterator;
use serde::Deserialize;
use serde::Serialize;
use sourmash::signature::Signature as SourmashSignature;
use std::collections::BTreeMap;
use std::fs;
use std::fs::remove_file;
use std::io;
use std::io::Write;
use std::path;
use std::sync::mpsc;
use std::sync::mpsc::Receiver;
use std::thread;
use std::{
    ffi::OsStr,
    fs::File,
    io::{BufRead, BufReader},
    path::PathBuf,
};

pub struct FileHandler {}

#[derive(Debug, Serialize, Deserialize)]
pub struct ShortSketchInfo {
    pub file_name: String,
    pub num_hashes: usize,
    pub kmer_size: u8,
    pub fscale: Option<u64>,
}

impl FileHandler {
    pub fn sketch_files(command: Commands, threads: Option<usize>) -> Result<()> {
        match command.to_owned() {
            Commands::Sketch {
                input,
                output,
                kmer_size,
                fscale,
                nmax,
                algorithm,
                format,
                singleton,
            } => {
                let files = FileHandler::test_and_collect_files(input, true)?;
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads.unwrap_or_default())
                    .build()?;

                let function = Function::from_alg(algorithm.clone(), kmer_size);

                let (send, recv) = mpsc::sync_channel(10);

                let multi_bar = MultiProgress::new();
                let multi_bar_clone = multi_bar.clone();

                let is_stdout = output.is_none();
                let handler = thread::spawn(move || {
                    FileHandler::write_output(fscale, output, format, recv, multi_bar_clone)
                });

                let pb = ProgressBar::new(files.len() as u64);
                let pb = multi_bar.add(pb);
                pb.set_style(indicatif::ProgressStyle::default_bar()
                        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta}) {msg}")
                        .unwrap()
                        .progress_chars("#>-"));
                let pb_clone = pb.clone();
                let _ = pool.install(|| {
                    files
                        .par_iter()
                        .progress_with(pb)
                        .try_for_each(|file_path| {
                            pb_clone.set_message(format!("{:?}", file_path.clone()));
                            match FileHandler::sketch_file(
                                file_path,
                                kmer_size,
                                fscale,
                                nmax,
                                singleton,
                                function.clone(),
                                algorithm.clone(),
                                is_stdout,
                            ) {
                                Ok(sig) => {
                                    send.send(sig).map_err(|_| anyhow!("Error while sending"))
                                }
                                Err(_) => {
                                    Err(anyhow!("Error while sketching file {:?}", file_path))
                                }
                            }
                        })
                });

                drop(send);

                Ok(handler
                    .join()
                    .map_err(|_| anyhow!("Unable to join threads"))??)
            }
            _ => Err(anyhow!("Wrong command")),
        }
    }

    pub fn sketch_file(
        input: &PathBuf,
        kmer_length: u8,
        fscale: Option<u64>,
        nmax: Option<u64>,
        singleton: bool,
        function: Function,
        algorithm: HashAlgorithms,
        _stdout: bool,
    ) -> Result<Signature> {
        //let start = std::time::Instant::now();
        let max_hash = if let Some(fscale) = fscale {
            (u64::MAX as f64 / fscale as f64) as u64
        } else {
            u64::MAX
        };
        let mut sketcher = sketcher::Sketcher::new(
            kmer_length,
            input
                .to_str()
                .ok_or_else(|| anyhow!("Unknown path"))?
                .to_string(),
            singleton,
            max_hash,
            nmax,
            function,
            algorithm,
        );
        let mut reader = parse_fastx_file(input)?;
        //let mut counter = 0;
        while let Some(record) = reader.next() {
            sketcher.process(&record?);
        }
        //let elapsed = start.elapsed().as_millis();
        // if !stdout {
        //     println!(
        //         "Processed {:?} with {} records, in {:?} seconds",
        //         input,
        //         counter,
        //         elapsed as f64 / 1000.0,
        //     );
        // }
        Ok(sketcher.finish())
    }

    pub fn write_output(
        fscale: Option<u64>,
        output: Option<PathBuf>,
        output_format: OutputFormats,
        signature_recv: Receiver<Signature>,
        multibar: MultiProgress,
    ) -> Result<()> {
        let stdout = output.is_none();

        match output_format {
            OutputFormats::Sourmash => {
                let mut output: Box<dyn Write> = match output {
                    Some(o) => Box::new(std::io::BufWriter::new(File::create(o)?)),
                    None => Box::new(std::io::BufWriter::new(io::stdout())),
                };
                output.write_all(b"[\n")?;
                let mut first = true;
                while let Ok(sig) = signature_recv.recv() {
                    let sourmash_sig: SourmashSignature = sig.into();
                    if !first {
                        output.write_all(b",\n")?;
                        first = false;
                    }
                    serde_json::to_writer(&mut output, &sourmash_sig)?;
                }
                output.write_all(b"]")?;
            }
            OutputFormats::Lmdb => {
                if stdout {
                    return Err(anyhow!("Output format lmdb is not supported for stdout"));
                }
                let Some(output) = output else {
                    return Err(anyhow!("Output folder is required for lmdb"));
                };
                if !output.is_dir() {
                    return Err(anyhow!(
                        "Output folder {:?} does not exist or is no directory",
                        output
                    ));
                }

                let heed_env = unsafe {
                    heed::EnvOpenOptions::new()
                        .map_size(10 * 1024 * 1024 * 1024 * 1024)
                        .max_dbs(2)
                        .flags(EnvFlags::WRITE_MAP | EnvFlags::MAP_ASYNC)
                        .open(output.clone())?
                };
                {
                    let mut write_txn = heed_env.write_txn()?;

                    let sigs_db = heed_env
                        .create_database::<U32<BigEndian>, SerdeBincode<ShortSketchInfo>>(
                            &mut write_txn,
                            Some("sigs"),
                        )?;
                    let hashes_db = heed_env
                        .database_options()
                        .types::<U64<BigEndian>, U32<BigEndian>>()
                        .name("hashes")
                        .flags(DatabaseFlags::DUP_SORT)
                        .create(&mut write_txn)?;

                    let mut counter: u32 = 0;
                    let mut hashes = BTreeMap::new();
                    while let Ok(sig) = signature_recv.recv() {
                        for sketch in sig.sketches {
                            sigs_db.put(
                                &mut write_txn,
                                &counter,
                                &ShortSketchInfo {
                                    file_name: sketch.name,
                                    num_hashes: sketch.num_kmers,
                                    kmer_size: sig.kmer_size,
                                    fscale,
                                },
                            )?;
                            for hash in sketch.hashes {
                                hashes.entry(hash).or_insert_with(Vec::new).push(counter);
                            }
                            counter += 1;
                        }
                        write_txn.commit()?;
                        write_txn = heed_env.write_txn()?;
                    }
                    let _ = multibar.println("Signatures finished, writing hashes");

                    let bar = multibar.add(ProgressBar::new(hashes.len() as u64));
                    bar.set_style(indicatif::ProgressStyle::default_bar()
                        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta}) {msg}")
                        .unwrap()
                        .progress_chars("#>-"));

                    for (hash, sigs) in hashes {
                        for sig in sigs {
                            hashes_db.put_with_flags(
                                &mut write_txn,
                                PutFlags::APPEND_DUP,
                                &hash,
                                &sig,
                            )?;
                        }
                        bar.inc(1);
                    }
                    write_txn.commit()?;
                }

                heed_env.prepare_for_closing().wait();

                let heed_env = unsafe {
                    heed::EnvOpenOptions::new()
                        .map_size(10 * 1024 * 1024 * 1024 * 1024)
                        .max_dbs(2)
                        .open(output.clone())?
                };

                let canonical_path = fs::canonicalize(format!("{}/", output.to_string_lossy()))?;
                println!(
                    "Compacting database to {:?}/compact.mdb",
                    canonical_path.to_string_lossy()
                );
                heed_env
                    .copy_to_file(
                        format!("{}/compact.mdb", canonical_path.to_string_lossy()),
                        heed::CompactionOption::Enabled,
                    )
                    .map_err(|e| {
                        println!("Error in copy file: {e}");
                        e
                    })?;

                remove_file(format!("{}/data.mdb", output.to_string_lossy())).map_err(|e| {
                    println!("Error deleting data.mdb: {e}");
                    e
                })?;
                remove_file(format!("{}/lock.mdb", output.to_string_lossy())).map_err(|e| {
                    println!("Error deleting lock.mdb: {e}");
                    e
                })?;
            }
        }

        Ok(())
    }

    pub fn read_signatures(input: &PathBuf) -> Result<Vec<Signature>> {
        Ok(
            sourmash::signature::Signature::from_path(path::Path::new(input))?
                .into_iter()
                .map(Signature::from)
                .collect(),
        )
    }

    pub fn concat(inputs: Vec<PathBuf>, output: PathBuf) -> Result<()> {
        let o_file = std::fs::File::create(output)?;
        let mut bufwriter = std::io::BufWriter::new(o_file);

        for input in inputs {
            let mut reader = BufReader::new(std::fs::File::open(input)?);
            while let Ok(result) =
                bincode::deserialize_from::<&mut BufReader<File>, Sketch>(&mut reader)
            {
                bincode::serialize_into(&mut bufwriter, &result)?;
            }
        }
        Ok(())
    }

    pub fn test_and_collect_files(input: Vec<PathBuf>, check_ext: bool) -> Result<Vec<PathBuf>> {
        let mut resulting_paths = Vec::new();
        let mut found_list: Option<PathBuf> = None;
        for path in input {
            if !path.exists() {
                return Err(anyhow::anyhow!("File {:?} does not exist", path));
            }
            if path.is_dir() {
                for p in path.read_dir()? {
                    let p = p?;
                    if p.path().is_file() {
                        if let Some(ext) = p.path().extension() {
                            if test_extension(ext) {
                                resulting_paths.push(p.path());
                            } else {
                                println!("Skipping file with invalid extension: {:?}", p.path());
                            }
                        } else {
                            println!("Skipping file without extension: {:?}", p.path());
                        }
                    } else {
                        println!("Skipping directory: {:?}", p.path());
                    }
                }
            }

            if path.is_file() {
                if let Some(ext) = path.extension() {
                    if test_extension(ext) || !check_ext {
                        resulting_paths.push(path);
                    } else if ext == "list" {
                        if resulting_paths.is_empty() {
                            found_list = Some(path);
                            break;
                        } else {
                            return Err(anyhow::anyhow!("Found multiple list files in {:?}", path));
                        }
                    } else {
                        return Err(anyhow::anyhow!("File with {:?} invalid extension", path));
                    }
                } else {
                    return Err(anyhow::anyhow!(
                        "File {:?} does not have an extension",
                        path
                    ));
                }
            }
        }

        if let Some(list) = found_list {
            let reader = BufReader::new(std::fs::File::open(list)?);
            for line in reader.lines() {
                let as_path_buf = PathBuf::from(line?);
                if as_path_buf.exists()
                    && test_extension(as_path_buf.extension().ok_or_else(|| {
                        anyhow::anyhow!("File {:?} does not have an extension", as_path_buf)
                    })?)
                    || !check_ext
                {
                    resulting_paths.push(as_path_buf);
                }
            }
        }
        Ok(resulting_paths)
    }

    pub fn write_result(result: &Vec<CompareResult>, output: PathBuf) -> Result<()> {
        let o_file = std::fs::File::create(output)?;
        let mut bufwriter = std::io::BufWriter::new(o_file);
        for r in result {
            writeln!(bufwriter, "{}", r)?;
        }
        Ok(())
    }
}

pub fn test_extension(ext: &OsStr) -> bool {
    !(ext != "fasta" && ext != "fa" && ext != "fastq" && ext != "fq" && ext != "gz")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_test_extension() {
        assert!(test_extension(OsStr::new("fasta")));
        assert!(test_extension(OsStr::new("fa")));
        assert!(test_extension(OsStr::new("fastq")));
        assert!(test_extension(OsStr::new("fq")));
        assert!(test_extension(OsStr::new("gz")));
        assert!(!test_extension(OsStr::new("txt")));
        assert!(!test_extension(OsStr::new("list")));
    }
}
