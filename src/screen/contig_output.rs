//! The contig output: the contigs of an assembly that share enough with an entry, written to one
//! gzip FASTA per assembly.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::batch::read_batches;
use super::{ScreenConfig, ScreenError};
use std::collections::HashSet;
use std::io::Write;
use std::path::Path;

/// Writes the contigs of the assembly at `path` whose ids are in `chosen` to a gzip FASTA in
/// `directory`, named like the assembly with `.gz`. Nothing is written when no contig is chosen.
/// These files form a much smaller set that still holds every close relative of the entries.
pub(super) fn write_contigs(
    path: &Path,
    name: &str,
    config: &ScreenConfig,
    directory: &Path,
    chosen: &HashSet<u32>,
) -> Result<(), ScreenError> {
    if chosen.is_empty() {
        return Ok(());
    }
    std::fs::create_dir_all(directory)?;
    let output = directory.join(if name.ends_with(".gz") {
        name.to_string()
    } else {
        format!("{name}.gz")
    });
    let file = tempfile::NamedTempFile::new_in(directory)?;
    let mut encoder = flate2::write::GzEncoder::new(
        std::io::BufWriter::new(file.as_file()),
        flate2::Compression::fast(),
    );
    let mut sink = |contig: u32, header: &[u8], sequence: &[u8]| -> Result<(), ScreenError> {
        if chosen.contains(&contig) {
            encoder.write_all(b">")?;
            encoder.write_all(header)?;
            encoder.write_all(b"\n")?;
            encoder.write_all(sequence)?;
            encoder.write_all(b"\n")?;
        }
        Ok(())
    };
    read_batches(
        path,
        config.batch_bases,
        config.rule.min_fragment_columns as usize,
        &mut |_, _| Ok(()),
        Some(&mut sink),
    )?;
    encoder.finish()?.flush()?;
    file.persist(&output).map_err(|error| error.error)?;
    Ok(())
}
