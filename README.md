[![Rust](https://img.shields.io/badge/built_with-Rust-dca282.svg)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/License-MIT-brightgreen.svg)](https://github.com/St4NNi/jam-rs/blob/main/LICENSE)
[![Crates.io](https://img.shields.io/crates/v/jam-rs.svg)](https://crates.io/crates/jam-rs)
[![Codecov](https://codecov.io/github/St4NNi/jam-rs/coverage.svg?branch=main)](https://codecov.io/gh/St4NNi/jam-rs)
[![Dependency status](https://deps.rs/repo/github/St4NNi/jam-rs/status.svg)](https://deps.rs/repo/github/St4NNi/jam-rs)

# jam-rs

Just another minhash (jam) implementation. A high-performance minhashing tool for genomic sequence similarity analysis, specifically optimized for plasmids and other small genomic elements.

Implements FracMinHash algorithm for rapid similarity comparisons with enhanced metadata tracking including GC content and sequence length categories tuned for typical plasmid ranges.

### Installation

Install the latest release from [crates.io](https://crates.io/):

```bash
cargo install jam-rs
```

Or install the development version from git:

```bash
cargo install --git https://github.com/St4NNi/jam-rs
```

### Key Features

- **Plasmid-optimized**: GC content and length categories specifically tuned for plasmid analysis (30-70% GC, 1kB-500kB lengths)
- **Fast sketching**: Entropy-filtered k-mers with optimized hash functions to exclude low-complexity regions
- **Rich metadata**: Enhanced metadata tracking with file index, GC category, and length category for each hash
- **Memory-efficient**: External sorting for processing datasets larger than available RAM
- **LMDB storage**: Fast random access and compact representation with dual database structure
- **Parallel execution**: File-level parallelization with configurable thread count

### Scaling Methods

Multiple scaling methods for different use cases:
- **FracMinHash** (`--fscale`): Restricts hash-space to a fraction of `u64::MAX` / `fscale`
- **Max hashes** (`--nmax`): Limits maximum number of hashes per sequence (memory control)
- **Singleton mode** (`--singleton`): Creates separate sketch per sequence record

### Usage

```console
$ jam --help
Just another (genomic) minhasher (jam), obviously blazingly fast

Usage: jam [OPTIONS] <COMMAND>

Commands:
  sketch  Sketch one or more files and write the result to an output file
  dist    Estimate containment of sequences against a database
  stats   Display statistics about a sketch database
  help    Print this message or the help of the given subcommand(s)

Options:
  -t, --threads <THREADS>  Number of threads to use [default: 1]
  -f, --force              Overwrite output files
  -s, --silent             Silent mode, no output to stdout
  -h, --help               Print help
  -V, --version            Print version
```

#### Sketching

Create sketches from FASTA/FASTQ files. Supports single files, multiple files, or directories.

```console
$ jam sketch --help
Sketch one or more files and write the result to an output file

Usage: jam sketch [OPTIONS] <INPUT>... --output <OUTPUT>

Arguments:
  <INPUT>...  Input file(s), directories, or file with list of files to be hashed

Options:
  -o, --output <OUTPUT>        Output file (.lmdb will be appended if not present)
  -k, --kmer-size <KMER_SIZE>  K-mer size, all sketches must have the same size to be compared and below 32 [default: 21]
      --fscale <FSCALE>        Scale the hash space to a minimum fraction of the maximum hash value (FracMinHash)
      --nmax <NMAX>            Maximum number of k-mers (per record) to be hashed, top cut-off
      --singleton              Create a separate sketch for each sequence record
  -t, --threads <THREADS>      Number of threads to use [default: 1]
  -f, --force                  Overwrite output files
  -h, --help                   Print help
```

Examples:
```bash
# Basic plasmid sketching
jam sketch -i plasmid.fasta -o sketch.lmdb

# Multiple plasmid files with custom k-mer size
jam sketch -i plasmids/ -o plasmid_db.lmdb -k 21 -t 8

# Large collections with memory limits
jam sketch -i large_collection/ -o database.lmdb --nmax 10000 --fscale 1000000

# Separate sketch per plasmid sequence
jam sketch -i multi_plasmids.fasta -o sketches.lmdb --singleton
```

#### Distance Calculation

Compare sequences against a sketch database. Supports both raw sequence files and pre-computed sketches.

```console
$ jam dist --help
Estimate containment of a (small) sketch against a subset of one or more sketches as database

Usage: jam dist [OPTIONS] --input <INPUT> --database <DATABASE>

Options:
  -i, --input <INPUT>        Input sketch or raw sequence file
  -d, --database <DATABASE>  Database sketch (.lmdb file)
  -o, --output <OUTPUT>      Output to file instead of stdout
  -c, --cutoff <CUTOFF>      Cut-off value for similarity/containment [default: 0.0]
  -h, --help                 Print help
```

Examples:
```bash
# Query plasmid against database
jam dist -i query_plasmid.fasta -d plasmid_db.lmdb -c 0.1 -o results.tsv

# Sketch-to-sketch comparison
jam dist -i query.lmdb -d plasmid_db.lmdb -c 0.05
```

#### Statistics

Display database statistics including hash counts and distribution analysis.

```console
$ jam stats --help
Display statistics about an LMDB database

Usage: jam stats [OPTIONS] --input <INPUT>

Options:
  -i, --input <INPUT>  Input LMDB database
  -s, --short          Short summary only
  -h, --help           Print help
```

Examples:
```bash
# Summary statistics
jam stats -i plasmid_db.lmdb --short

# Detailed distributions
jam stats -i plasmid_db.lmdb
```

### Output Format

Distance results are tab-separated with columns:
```
query  target  containment_query_in_target  containment_target_in_query  jaccard  shared_hashes  query_hashes  target_hashes
```

Statistics include hash counts, GC content distribution, and sequence length categories optimized for plasmids and small genomic elements.

### Algorithm

JAM uses entropy-filtered k-mers to exclude low-complexity regions, stores rich metadata (file index, GC category, length category) with each hash, and employs external sorting for memory-efficient processing of large datasets. The categorization system is specifically tuned for plasmid analysis, with fine-grained bins in typical plasmid size and GC content ranges.

### License

This project is licensed under the MIT license. See the [LICENSE](LICENSE) file for more info.

### Feedback & Contributions

If you have any ideas, suggestions, or issues, please don't hesitate to open an issue and/or PR. Contributions to this project are always welcome! We appreciate your help in making this project better.

### Credits

This tool is inspired by [finch-rs](https://github.com/onecodex/finch-rs) and [sourmash](https://github.com/sourmash-bio/sourmash). Check them out if you need a more mature ecosystem with well tested hash functions and more features.