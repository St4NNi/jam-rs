mod compare;
mod file_io;
mod hasher;
mod sketcher;

fn main() {
    let filename = "test.fasta.gz";
    let reader = file_io::FileHandler {};
    reader
        .sketch_file(filename, "output.sketch", 21, 0.01)
        .unwrap();
}
