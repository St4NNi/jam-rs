//! Joins shards into one. Every entry keeps its positions section and its packed sequence
//! unchanged; only the sketch is built again, with color sets over all entries.

use super::build::{complete_shard, in_partial_directory};
use super::manifest::{Manifest, PACKED_DIRECTORY, Parameters, check_names};
use super::positions::{POSITIONS_FILE, Positions};
use super::{IndexError, MappedFile};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// How shards are joined.
#[derive(Clone, Debug)]
pub struct JoinOptions {
    /// Threads that build the sketch.
    pub threads: usize,
    /// Directory for temporary files.
    pub temp_dir: PathBuf,
    /// Bucket bits to use instead of the ones chosen from the joined shard's size.
    pub bucket_bits: Option<u8>,
    /// Link every packed sequence into the new shard's packed directory. Without it, the
    /// joined shard is opened with the directory that holds the packed sequences of all shards.
    pub link_packed: bool,
}

/// Joins the shards in `inputs`, in order, into the shard directory `output`.
pub fn join_shards(
    inputs: &[PathBuf],
    output: &Path,
    options: &JoinOptions,
) -> Result<Manifest, IndexError> {
    let manifests = inputs
        .iter()
        .map(|input| Manifest::read(input))
        .collect::<Result<Vec<_>, IndexError>>()?;
    let first = manifests
        .first()
        .ok_or_else(|| IndexError::Input("no shard to join".into()))?;
    let parameters = first.parameters;
    // Bucket bits follow each shard's size, and the joined sketch chooses its own.
    let same = |other: &Parameters| {
        Parameters {
            bucket_bits: parameters.bucket_bits,
            ..*other
        } == parameters
    };
    if manifests.iter().any(|manifest| !same(&manifest.parameters)) {
        return Err(IndexError::Input(
            "shards with different parameters cannot be joined".into(),
        ));
    }
    let entries: Vec<_> = manifests
        .iter()
        .flat_map(|manifest| manifest.entries.iter().cloned())
        .collect();
    check_names(entries.iter().map(|entry| entry.name.as_str()))?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.threads.max(1))
        .build()
        .map_err(|error| IndexError::Input(format!("thread pool: {error}")))?;
    in_partial_directory(output, |directory| {
        if options.link_packed {
            link_packed(inputs, &manifests, &directory.join(PACKED_DIRECTORY))?;
        }
        pool.install(|| {
            complete_shard(
                directory,
                parameters,
                entries,
                &options.temp_dir,
                options.bucket_bits,
                |writer| {
                    for input in inputs {
                        let file = MappedFile::open(&input.join(POSITIONS_FILE), false)?;
                        let positions = Positions::open(Arc::new(file))?;
                        for entry in 0..positions.entries() {
                            // The section's CRC-32 is checked as it is read, so damage is not copied.
                            let bytes = positions.section_bytes(entry)?;
                            writer.push_section(positions.section(entry)?, &bytes)?;
                        }
                    }
                    Ok(())
                },
            )
        })
    })
}

/// Links the packed sequences of all inputs into `packed`.
fn link_packed(
    inputs: &[PathBuf],
    manifests: &[Manifest],
    packed: &Path,
) -> Result<(), IndexError> {
    std::fs::create_dir(packed)?;
    for (input, manifest) in inputs.iter().zip(manifests) {
        let source = std::fs::canonicalize(input.join(PACKED_DIRECTORY))?;
        for entry in &manifest.entries {
            std::os::unix::fs::symlink(
                source.join(&entry.packed.path),
                packed.join(&entry.packed.path),
            )?;
        }
    }
    Ok(())
}
