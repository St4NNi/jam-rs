use std::{
    io::{BufWriter, Write},
    path::PathBuf,
};

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

const TOTAL: u64 = 4_398_046_511_103; // 2^42 - 1
const CHUNKS: u64 = 2000;

fn get_chunk_iter(total: u64, chunks: u64, idx: u64) -> std::ops::Range<u64> {
    if idx >= chunks {
        panic!("Index out of bounds: {idx} >= {chunks}");
    }
    let base = total / chunks;
    let rem = total % chunks;
    let start = if idx < rem {
        idx * (base + 1)
    } else {
        rem * (base + 1) + (idx - rem) * base
    };
    let size = if idx < rem { base + 1 } else { base };
    start..(start + size)
}

fn do_hash(kmer: u64, jam: &mut Vec<u64>) {
    jam.push(jamhash::jamhash_u64(kmer));
}

fn write(name: &str, chunk: u64, path: PathBuf, mut values: Vec<u64>) {
    let file_path = path.join(format!("{name}_{chunk}.bin"));
    let mut file = BufWriter::with_capacity(
        1024 * 1024 * 512, // 512 MB buffer
        std::fs::File::create(file_path).expect("Failed to create file"),
    );

    for value in values.drain(..) {
        file.write(&value.to_le_bytes())
            .expect("Failed to write kmer");
    }

    file.flush().expect("Failed to flush file");
}

pub fn main() {
    let chunk = std::env::args().nth(1).unwrap().parse::<u64>().unwrap();
    let path = PathBuf::from(std::env::args().nth(2).unwrap());

    let kmer_len = 21;

    let mut rc = 0u64;
    let mut kc = 0u64;
    let mut skipped_rc = 0u64;

    let mut jam_hashes: Vec<u64> = Vec::new();

    let range = get_chunk_iter(TOTAL, CHUNKS, chunk);
    let range_start = range.start;
    let range_end = range.end;

    for x in range {
        let kmer = x;
        let rc_kmer = reverse_complement(kmer, kmer_len);

        if kmer < rc_kmer {
            kc += 1;
            do_hash(rc_kmer, &mut jam_hashes);
        } else {
            if rc_kmer >= range_start && rc_kmer < range_end {
                rc += 1;
                do_hash(rc_kmer, &mut jam_hashes);
            } else {
                skipped_rc += 1;
            }
        }
    }

    // Ensure the directory exists and is a dir
    if !path.is_dir() || !path.exists() {
        println!(
            "Path is not a directory or does not exist: {}",
            path.display()
        );
    }

    // Write the results to a file
    write("jam_hashes", chunk, path.clone(), jam_hashes);

    println!("Range: {range_start}..{range_end}");
    println!("Kmers: {kc}");
    println!("Reverse complements: {rc} | Skipped: {skipped_rc}");
    println!("Total: {}", kc + rc);
}
