use anyhow::{Context, Result};
use heed::{Env, EnvFlags, EnvOpenOptions};
use std::path::Path;

pub enum EnvMode {
    ReadOnly,
    ReadWrite,
}

pub fn open_env(path: &Path, mode: EnvMode) -> Result<Env> {
    let map_size = 10 * 1024 * 1024 * 1024 * 1024; // 10TB

    // SAFETY: heed requires unsafe for mmap; we control file access
    let env = unsafe {
        let mut opts = EnvOpenOptions::new();
        opts.max_dbs(3).map_size(map_size);

        match mode {
            EnvMode::ReadOnly => {
                opts.flags(EnvFlags::READ_ONLY | EnvFlags::NO_SUB_DIR);
            }
            EnvMode::ReadWrite => {
                opts.flags(
                    EnvFlags::NO_SUB_DIR | EnvFlags::WRITE_MAP | EnvFlags::MAP_ASYNC,
                );
            }
        }

        opts.open(path)
            .with_context(|| format!("Failed to open database: {}", path.display()))?
    };

    Ok(env)
}
