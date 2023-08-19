use crate::cli::Commands;
use crate::compare::CompareResult;
use crate::hash_functions::Function;
use crate::sketcher;
use anyhow::anyhow;
use anyhow::Result;
use needletail::parse_fastx_file;
use rayon::prelude::IntoParallelRefIterator;
use rayon::prelude::ParallelIterator;
use std::io::Write;
use std::ops::DerefMut;
use std::sync::Arc;
use std::sync::Mutex;
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
            } => {
                let files = FileHandler::test_and_collect_files(input, true)?;
                let output = Arc::new(Mutex::new(File::create(output.unwrap())?)); //TODO: stdout support
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads.unwrap_or_default())
                    .build()?;

                let function = Function::from_alg(algorithm, kmer_size);

                pool.install(|| {
                    files
                        .par_iter()
                        .try_for_each(|file_path| {
                            let output = output.clone();
                            FileHandler::sketch_file(
                                file_path,
                                output,
                                kmer_size,
                                fscale,
                                kscale,
                                nmin,
                                nmax,
                                singleton,
                                false,
                                function.clone(),
                            )
                        })
                        .unwrap()
                });
                Ok(())
            }
            _ => Err(anyhow!("Wrong command")),
        }
    }

    pub fn sketch_file(
        input: &PathBuf,
        output: Arc<Mutex<File>>,
        kmer_length: u8,
        fscale: u64,
        kscale: u64,
        nmin: Option<u64>,
        nmax: Option<u64>,
        singleton: bool,
        stats: bool,
        function: Function,
    ) -> Result<()> {
        let mut x = fs::metadata(input)?.len();
        if input.ends_with(".gz") {
            // Approximate the size of the uncompressed file
            x *= 3;
        }
        let start = std::time::Instant::now();
        let kscale = (x as f64 / kscale as f64) as u64;
        let max_hash = (u64::MAX as f64 / fscale as f64) as u64;
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
        );
        let mut reader = parse_fastx_file(input)?;
        let mut counter = 0;

        while let Some(record) = reader.next() {
            sketcher.process_small(&record?);
            counter += 1;
        }

        let output = output.clone();
        let mut o_file = output.lock().unwrap();
        let mut bufwriter = std::io::BufWriter::new(o_file.deref_mut());
        bincode::serialize_into(&mut bufwriter, &sketcher.finalize())?;
        let elapsed = start.elapsed().as_millis();
        println!(
            "Processed {:?} with {} records, in {:?} seconds",
            input,
            counter,
            elapsed as f64 / 1000.0,
        );
        bufwriter.flush()?;
        Ok(())
    }

    pub fn read_sketches(input: &PathBuf) -> Result<Vec<sketcher::Sketch>> {
        let mut vec: Vec<sketcher::Sketch> = vec![];

        let mut reader = BufReader::new(std::fs::File::open(input)?);

        while let Ok(result) =
            bincode::deserialize_from::<&mut BufReader<File>, sketcher::Sketch>(&mut reader)
        {
            vec.push(result);
        }

        Ok(vec)
    }

    pub fn concat(inputs: Vec<PathBuf>, output: PathBuf) -> Result<()> {
        let o_file = std::fs::File::create(output)?;
        let mut bufwriter = std::io::BufWriter::new(o_file);

        for input in inputs {
            let mut reader = BufReader::new(std::fs::File::open(input)?);
            while let Ok(result) =
                bincode::deserialize_from::<&mut BufReader<File>, sketcher::Sketch>(&mut reader)
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
