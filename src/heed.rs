use std::path::PathBuf;

use byteorder::BigEndian;
use heed::{
    types::{SerdeBincode, U32, U64},
    DatabaseFlags, EnvFlags,
};

use crate::file_io::ShortSketchInfo;

pub struct HeedHandler {
    heed_env: heed::Env,
    signatures: heed::Database<U32<BigEndian>, SerdeBincode<ShortSketchInfo>>,
    hashes: heed::Database<U64<BigEndian>, U32<BigEndian>>,
}

impl HeedHandler {
    pub fn new_ro(path: PathBuf) -> anyhow::Result<Self> {
        let heed_env = if path.is_dir() {
            unsafe {
                heed::EnvOpenOptions::new()
                    .map_size(10 * 1024 * 1024 * 1024 * 1024)
                    .max_dbs(2)
                    .flags(EnvFlags::READ_ONLY)
                    .open(path.clone())?
            }
        } else {
            unsafe {
                heed::EnvOpenOptions::new()
                    .map_size(10 * 1024 * 1024 * 1024 * 1024)
                    .max_dbs(2)
                    .flags(EnvFlags::READ_ONLY | EnvFlags::NO_SUB_DIR)
                    .open(path.clone())?
            }
        };

        let rtxn = heed_env.read_txn()?;

        let sigs_db = heed_env
            .open_database::<U32<BigEndian>, SerdeBincode<ShortSketchInfo>>(&rtxn, Some("sigs"))?
            .ok_or_else(|| anyhow::anyhow!("Unable to open signatures database"))?;
        let hashes = heed_env
            .database_options()
            .types::<U64<BigEndian>, U32<BigEndian>>()
            .name("hashes")
            .flags(DatabaseFlags::DUP_SORT)
            .open(&rtxn)?
            .ok_or_else(|| anyhow::anyhow!("Unable to open signatures database"))?;
        rtxn.commit()?;
        Ok(HeedHandler {
            heed_env,
            signatures: sigs_db,
            hashes,
        })
    }

    pub fn summarize_stats(&self) -> anyhow::Result<()> {
        let rtxn = self.heed_env.read_txn()?;
        let num_of_sigs = self.signatures.len(&rtxn)?;
        println!("Number of signatures: {}", num_of_sigs);
        let num_of_hashes = self.hashes.len(&rtxn)?;
        println!("Number of hashes: {}", num_of_hashes);
        Ok(())
    }

    pub fn detail_sigs(&self) -> anyhow::Result<()> {
        let rtxn = self.heed_env.read_txn()?;
        for (_, value) in self.signatures.iter(&rtxn)?.enumerate() {
            let (_, value) = value?;
            println!(
                "{},{:?},{},{}",
                value.file_name, value.fscale, value.kmer_size, value.num_hashes
            );
        }
        Ok(())
    }
}
