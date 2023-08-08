use jam_rs::{compare, file_io};

#[tokio::main]
async fn main() {
    let filename = "test.fasta.gz";
    let reader = file_io::FileHandler {};
    reader
        .sketch_file(filename, "output.sketch", 21, 0.01)
        .unwrap();

    let start = std::time::Instant::now();
    let sketch1 = reader.read_sketches("output.sketch").unwrap();
    let sketch2 = reader.read_sketches("output.sketch").unwrap();

    let mut comparator = compare::MultiComp::new(sketch1, sketch2);
    comparator.compare();
    println!("{}", comparator.finalize().first().unwrap());
    let elapsed = start.elapsed().as_millis();
    println!("Processed in {:?} seconds", elapsed as f64 / 1000.0);
}
