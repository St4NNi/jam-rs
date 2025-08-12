use jam_rs::hash_functions::jamhash;

#[inline]
pub fn murmur3_u64(kmer: u64) -> u64 {
    fastmurmur3::murmur3_x64_128(&kmer.to_be_bytes(), 42) as u64
}

fn main() {
    let mut vec = Vec::with_capacity(34359738368 / 8);

    for x in 0..(34359738368u64 / 8) {
        let hash = jamhash(x); //jamhash(x);
        if x % 100_000_000 == 0 {
            println!("Processing: {} / {}", x, 34359738368u64 / 8);
        }
        vec.push(hash);
    }

    vec.sort_unstable();
    let len = vec.len();
    println!("Total elements: {}", len);
    vec.dedup();
    let len_after = vec.len();

    println!("Number of duplicates: {}", len - len_after);
}
