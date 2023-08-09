[![Rust](https://img.shields.io/badge/built_with-Rust-dca282.svg)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/License-MIT-brightgreen.svg)](https://github.com/St4NNi/jam-rs/blob/main/LICENSE)
![CI](https://github.com/St4NNi/jam-rs/actions/workflows/push.yaml/badge.svg)
[![Codecov](https://codecov.io/github/St4NNi/jam-rs/coverage.svg?branch=main)](https://codecov.io/gh/St4NNi/jam-rs)
[![Dependency status](https://deps.rs/repo/github/St4NNi/jam-rs/status.svg)](https://deps.rs/repo/github/St4NNi/jam-rs)
___
# jam-rs

Just another minhash (jam) implementation. A high performance minhash variant to screen extremely large (metagenomic) datasets in a very short timeframe.
Implements parts of the ScaledMinHash / FracMinHash algorithm described in [sourmash](https://joss.theoj.org/papers/10.21105/joss.00027).

Unlike traditional implementations like [sourmash](https://joss.theoj.org/papers/10.21105/joss.00027) or [mash](https://doi.org/10.1186/s13059-016-0997-x) this version focuses on estimating containment of small sequences in a large set. This can be used to screen terabytes of data in just a few seconds / minutes.

### Comparison

- xxhash3 instead of murmurhash
- No jaccard similarity since this is meaningless when comparing small embedd
- Quick estimation of hash fractions from file-size

### Usage

```console
$ jam
Just another minhasher, obviously blazingly fast

Usage: jam [OPTIONS] <COMMAND>

Commands:
  sketch   Sketches one or more files and writes the result to an output file
  merge    Merge multiple input sketches into a single sketch
  compare  Compare a raw file or sketch against one or more sketches as database Requires all sketches to have the same kmer size
  help     Print this message or the help of the given subcommand(s)

Options:
  -t, --threads <THREADS>  Number of threads to use [default: 1]
  -f, --force              Overwrite output files
  -h, --help               Print help (see more with '--help')
  -V, --version            Print version
```

#### Sketching

The easiest way to sketch files is to use the `jam sketch` command. This accepts one or more input files (fastx / fastx.gz) or a `.list` file with a full list of input files. And sketches all inputs to a specific outpuf sketch file.

```console
$ jam sketch
Sketches one or more files and writes the result to an output file

Usage: jam sketch [OPTIONS] --input <INPUT> --output <OUTPUT>

Options:
  -i, --input <INPUT>          Input file, directory or file with list of files to be hashed
  -o, --output <OUTPUT>        Output file
  -k, --kmer-size <KMER_SIZE>  kmer size all sketches to be compared must have the same size [default: 21]
  -s, --scale <SCALE>          The estimated scaling factor to apply [default: 0.001]
  -t, --threads <THREADS>      Number of threads to use [default: 1]
  -f, --force                  Overwrite output files
  -h, --help                   Print help
```

#### Compare

Compare one or more inputs vs. a large set of database sketches. Optionally specify a minimum cutoff in percent of matching kmers. Output is optional if not specified the result will be printed to stdout.

```console
$ jam compare
Compare a raw file or sketch against one or more sketches as database Requires all sketches to have the same kmer size

Usage: jam compare [OPTIONS] --input <INPUT> --database <DATABASE>

Options:
  -i, --input <INPUT>        Input sketch or raw file
  -d, --database <DATABASE>  Database sketch(es)
  -o, --output <OUTPUT>      Output to file instead of stdout
  -c, --cutoff <CUTOFF>      Cut-off value for similarity [default: 0.0]
  -t, --threads <THREADS>    Number of threads to use [default: 1]
  -f, --force                Overwrite output files
  -h, --help                 Print help
```




#### Merge

Merge multiple sketches into one large one.

```console
$ jam merge
Merge multiple input sketches into a single sketch

Usage: jam merge [OPTIONS] --output <OUTPUT> [INPUTS]...

Arguments:
  [INPUTS]...  One or more input sketches

Options:
  -o, --output <OUTPUT>    Output file
  -t, --threads <THREADS>  Number of threads to use [default: 1]
  -f, --force              Overwrite output files
  -h, --help               Print help
```



