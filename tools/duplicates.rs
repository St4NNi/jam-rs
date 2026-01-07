use std::sync::mpsc;
use std::thread;


/// Reverse complement a `BitKmer` (reverses the sequence and swaps A<>T and G<>C)
pub fn reverse_complement(sequence: u64, len: u8) -> u64 {
    // FIXME: this is not going to work with BitKmers of u128 or u32
    // inspired from https://www.biostars.org/p/113640/
    let mut new_kmer = sequence;
    // reverse it
    new_kmer =
        ((new_kmer >> 2) & 0x3333_3333_3333_3333) | ((new_kmer & 0x3333_3333_3333_3333) << 2);
    new_kmer =
        ((new_kmer >> 4) & 0x0F0F_0F0F_0F0F_0F0F) | ((new_kmer & 0x0F0F_0F0F_0F0F_0F0F) << 4);
    new_kmer =
        ((new_kmer >> 8) & 0x00FF_00FF_00FF_00FF) | ((new_kmer & 0x00FF_00FF_00FF_00FF) << 8);
    new_kmer =
        ((new_kmer >> 16) & 0x0000_FFFF_0000_FFFF) | ((new_kmer & 0x0000_FFFF_0000_FFFF) << 16);
    new_kmer =
        ((new_kmer >> 32) & 0x0000_0000_FFFF_FFFF) | ((new_kmer & 0x0000_0000_FFFF_FFFF) << 32);
    // complement it
    new_kmer ^= 0xFFFF_FFFF_FFFF_FFFF;
    // shift it to the right size
    new_kmer >>= 2 * (32 - len);
    new_kmer
}

pub fn bit_kmer_to_string(kmer: u64, len: u8) -> String {
    let mut s = String::with_capacity(len as usize);
    for i in (0..len).rev() {
        let base = (kmer >> (2 * i)) & 0b11;
        s.push(match base {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => unreachable!(),
        });
    }
    s
}

/// Count hash collisions in a sorted vector of hashes
fn count_collisions(iter: impl ExactSizeIterator<Item = u64>) -> usize {
    if iter.len() <= 1 {
        return 0;
    }

    let mut collisions = 0;
    let mut last_elem = 0u64;

    for elem in iter {
        if last_elem == elem {
            collisions += 1;
        }
        last_elem = elem;
    }

    collisions
}

const TOTAL: u64 = 4_398_046_511_103; // 2^42 - 1
const SAMPLE_RATE: u64 = 100;

#[inline]
pub fn murmur3_u64(kmer: u64) -> u64 {
    fastmurmur3::murmur3_x64_128(&kmer.to_be_bytes(), 42) as u64
}

// pub fn foldhash(kmer: u64) -> u64 {
//     let mut hasher = foldhash::quality::FixedState::with_seed(42).build_hasher();
//     hasher.write_u64(kmer);
//     hasher.finish()
// }

// pub fn rapidhash(kmer: u64) -> u64 {
//     let mut hasher = rapidhash::quality::RapidHasher::new(42);
//     hasher.write_u64(kmer);
//     hasher.finish()
// }

pub fn main() {
    let kmer_len = 21;

    println!(
        "Expected RAM usage: ~{:.2} GB",
        ((TOTAL / SAMPLE_RATE) * 8) as f64 / (1024.0 * 1024.0 * 1024.0)
    );

    println!("Starting k-mer hash collision test...");
    println!("Total k-mer space: {}", TOTAL);
    println!("Sample rate: 1/{}", SAMPLE_RATE);

    let start_time = std::time::Instant::now();

    // Calculate work distribution
    let num_threads = thread::available_parallelism().unwrap().get();
    let chunk_size = TOTAL / num_threads as u64;

    println!("Using {} threads", num_threads);

    // Create channel for collecting results
    let (tx, rx) = mpsc::channel::<Vec<u64>>();

    // Spawn worker threads
    let mut handles = Vec::new();

    for thread_id in 0..num_threads {
        let tx_clone = tx.clone();

        let handle = thread::spawn(move || {
            let start = thread_id as u64 * chunk_size;
            let end = if thread_id == num_threads - 1 {
                TOTAL
            } else {
                (thread_id + 1) as u64 * chunk_size
            };

            // Process in smaller batches to avoid memory spikes
            const BATCH_SIZE: usize = 50_000;
            let mut batch = Vec::with_capacity(BATCH_SIZE);

            for kmer in start..end {
                if kmer % SAMPLE_RATE == 0 {
                    let rc_kmer = reverse_complement(kmer, kmer_len);

                    // Skip if RC is smaller and would also be sampled
                    if rc_kmer < kmer && rc_kmer % SAMPLE_RATE == 0 {
                        continue;
                    }

                    let canonical = kmer.min(rc_kmer);

                    //batch.push(murmur3_u64(canonical));
                    batch.push(jamhash::jamhash_u64(canonical));
                    //batch.push(foldhash(canonical));
                    //batch.push(rapidhash(canonical));

                    // Send batch when it's full
                    if batch.len() >= BATCH_SIZE {
                        tx_clone.send(batch).unwrap();
                        batch = Vec::with_capacity(BATCH_SIZE);
                    }
                }
            }

            // Send any remaining items in the final batch
            if !batch.is_empty() {
                tx_clone.send(batch).unwrap();
            }
        });

        handles.push(handle);
    }

    // Drop the original sender so receiver knows when all threads are done
    drop(tx);

    // Estimate total capacity upfront
    let estimated_total = (TOTAL / SAMPLE_RATE) as usize;
    let mut hashes = Vec::with_capacity(estimated_total);

    // Collect batches as they arrive and extend directly
    while let Ok(mut batch) = rx.recv() {
        hashes.extend(batch.drain(..));
    }

    // Wait for all threads to complete
    for handle in handles {
        handle.join().unwrap();
    }

    let collection_time = start_time.elapsed();
    println!("Collected {} hashes in {:?}", hashes.len(), collection_time);

    // Sort hashes (can still use rayon for this part)
    let sort_start = std::time::Instant::now();
    hashes.sort_unstable();
    let sort_time = sort_start.elapsed();
    println!("Sorted hashes in {:?}", sort_time);

    // Count collisions
    let collision_start = std::time::Instant::now();
    let len = hashes.len();
    let collision_count = count_collisions(hashes.into_iter());
    let collision_time = collision_start.elapsed();

    let total_time = start_time.elapsed();

    println!("\n=== Results ===");
    println!("Total canonical k-mers processed: {}", len);
    println!("Hash collisions found: {}", collision_count);
    println!(
        "Collision rate: {:.8}%",
        (collision_count as f64 / len as f64) * 100.0
    );
    println!(
        "Memory usage: ~{:.2} GB",
        (len * 8) as f64 / (1024.0 * 1024.0 * 1024.0)
    );
    println!("\n=== Timing ===");
    println!("Collection time: {:?}", collection_time);
    println!("Sorting time: {:?}", sort_time);
    println!("Collision detection time: {:?}", collision_time);
    println!("Total time: {:?}", total_time);
}
