use jam_rs::{compare, file_io};

fn main() {
    let reader = file_io::FileHandler {};
    let start = std::time::Instant::now();
    let sketch1 = reader.read_sketches("output.sketch").unwrap();
    let sketch2 = reader.read_sketches("combined.sketch").unwrap();

    let mut comparator = compare::MultiComp::new(sketch1, sketch2);
    comparator.compare();
    for x in comparator.finalize().iter() {
        println!("{}", x);
    }
    let elapsed = start.elapsed().as_millis();
    println!("Processed in {:?} seconds", elapsed as f64 / 1000.0);
}
