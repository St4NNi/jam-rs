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

    let start = std::time::Instant::now();
    let sketch1 = reader.read_sketch("output.sketch").unwrap();
    let sketch2 = reader.read_sketch("output.sketch").unwrap();

    let mut comparator = compare::Comparator::new(sketch1, sketch2);
    comparator.compare();
    println!("{:?}", comparator.finalize());
    let elapsed = start.elapsed().as_millis();
    println!("Processed in {:?} seconds", elapsed as f64 / 1000.0);
}
