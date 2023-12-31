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
use needletail::parse_fastx_file;
use rayon::prelude::IntoParallelRefIterator;
use rayon::prelude::ParallelIterator;
use sourmash::signature::Signature as SourmashSignature;
use std::io;
use std::io::Write;
use std::sync::mpsc;
use std::sync::mpsc::Receiver;
use std::thread;
use std::{
    ffi::OsStr,
    fs::{self, File},
    io::{BufRead, BufReader},
    path::PathBuf,
};

pub struct FileHandler {}

impl FileHandler {
    pub fn sketch_files(command: Commands, threads: Option<usize>) -> Result<()> {
        match command.to_owned() {
            Commands::Sketch {
                input,
                output,
                kmer_size,
                fscale,
                kscale,
                nmin,
                nmax,
                algorithm,
                format,
                singleton,
                stats,
            } => {
                let files = FileHandler::test_and_collect_files(input, true)?;
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads.unwrap_or_default())
                    .build()?;

                let function = Function::from_alg(algorithm.clone(), kmer_size);

                let (send, recv) = mpsc::channel();

                let is_stdout = output.is_none();
                let handler = thread::spawn(|| {
                    FileHandler::write_output(output, format, recv)
                    // thread code
                });

                let _ = pool.install(|| {
                    files.par_iter().try_for_each(|file_path| {
                        match FileHandler::sketch_file(
                            file_path,
                            kmer_size,
                            fscale,
                            kscale,
                            nmin,
                            nmax,
                            singleton,
                            stats,
                            function.clone(),
                            algorithm.clone(),
                            is_stdout,
                        ) {
                            Ok(sig) => send.send(sig).map_err(|_| anyhow!("Error while sending")),
                            Err(_) => Err(anyhow!("Error while sketching file {:?}", file_path)),
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
        kscale: Option<u64>,
        nmin: Option<u64>,
        nmax: Option<u64>,
        singleton: bool,
        stats: bool,
        function: Function,
        algorithm: HashAlgorithms,
        stdout: bool,
    ) -> Result<Signature> {
        let mut x = fs::metadata(input)?.len();
        if let Some(ext) = input.extension() {
            if let Some(ext_str) = ext.to_str() {
                if ext_str == "gz" {
                    // Approximate the size of uncompressed file
                    x *= 3;
                }
            }
        }
        let start = std::time::Instant::now();
        let kscale = if let Some(kscale) = kscale {
            (x as f64 / kscale as f64) as u64
        } else {
            u64::MAX
        };
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
            stats,
            kscale,
            max_hash,
            nmin,
            nmax,
            function,
            algorithm,
        );
        let mut reader = parse_fastx_file(input)?;
        let mut counter = 0;
        while let Some(record) = reader.next() {
            sketcher.process(&record?);
            counter += 1;
        }
        let elapsed = start.elapsed().as_millis();
        if !stdout {
            println!(
                "Processed {:?} with {} records, in {:?} seconds",
                input,
                counter,
                elapsed as f64 / 1000.0,
            );
        }
        Ok(sketcher.finish())
    }

    pub fn write_output(
        output: Option<PathBuf>,
        output_format: OutputFormats,
        signature_recv: Receiver<Signature>,
    ) -> Result<()> {
        let stdout = output.is_none();
        let mut output: Box<dyn Write> = match output {
            Some(o) => Box::new(std::io::BufWriter::new(File::create(o)?)),
            None => Box::new(std::io::BufWriter::new(io::stdout())),
        };

        match output_format {
            OutputFormats::Bin => {
                while let Ok(sig) = signature_recv.recv() {
                    let name = sig.file_name.clone();
                    let len = sig.sketches.first().unwrap().hashes.len();
                    bincode::serialize_into(&mut output, &vec![sig])?;
                    if !stdout {
                        println!("Wrote signature: {:?} with {:?} hashes.", name, len);
                    }
                }
            }
            OutputFormats::Sourmash => {
                while let Ok(sig) = signature_recv.recv() {
                    let sourmash_sig: SourmashSignature = sig.into();
                    serde_json::to_writer(&mut output, &vec![sourmash_sig])?;
                }
            }
        }

        Ok(())
    }

    pub fn read_signatures(input: &PathBuf) -> Result<Vec<Signature>> {
        let read_to_bytes = std::fs::read(input)?;
        Ok(bincode::deserialize_from(read_to_bytes.as_slice()).unwrap())
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
                            } else if ext == "list" {
                                if resulting_paths.is_empty() {
                                    found_list = Some(p.path());
                                    break;
                                } else {
                                    return Err(anyhow::anyhow!(
                                        "Found multiple list files in {:?}",
                                        path
                                    ));
                                }
                            } else {
                                return Err(anyhow::anyhow!(
                                    "File with {:?} invalid extension",
                                    path
                                ));
                            }
                        } else {
                            return Err(anyhow::anyhow!(
                                "File {:?} does not have an extension",
                                p.path()
                            ));
                        }
                    } else {
                        return Err(anyhow::anyhow!("File {:?} is not a file", p.path()));
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
