use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read};
use std::path::{Path, PathBuf};

pub const HASH_ID: &str = "jamhash_u64_v1";
pub const HASH_ZERO_POLICY: &str = "excluded";
pub const OUTPUT_SCHEMA_VERSION: &str = "1.0.0";
pub const MANIFEST_SCHEMA_VERSION: &str = "1.0.0";

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileIdentity {
    pub path: String,
    pub size_bytes: u64,
    pub sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiasDatabaseMetadata {
    pub table_id: String,
    pub source_path: String,
    pub sha256: String,
    pub kmer_size: u8,
    pub base_fscale: u64,
    pub cms_width: usize,
    pub cms_depth: usize,
    pub alpha: f32,
    pub filter_mode: String,
    pub target_fscale: u64,
    pub negative_fscale: String,
    pub unseen_fscale: u64,
    pub positive_retention: f32,
    pub negative_retention: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatabaseManifest {
    pub schema_version: String,
    pub output_schema_version: String,
    pub database_format_version: u32,
    pub jam_rs_version: String,
    pub source_commit: String,
    pub source_dirty: Option<bool>,
    pub hash_id: String,
    pub hash_zero_policy: String,
    pub kmer_size: u8,
    pub fscale: u64,
    pub hash_threshold: u64,
    pub entropy_threshold: f64,
    pub bias: Option<BiasDatabaseMetadata>,
    pub input_catalog_files: Vec<FileIdentity>,
    pub catalog_manifest_sha256: String,
    pub sample_count: u32,
    pub entry_count: u64,
    pub unique_hash_count: u64,
    pub database_file: String,
    pub database_size_bytes: u64,
    pub database_sha256: String,
    pub creation_command: Vec<String>,
    pub creation_time_unix_seconds: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiasTableManifest {
    pub schema_version: String,
    pub jam_rs_version: String,
    pub source_commit: String,
    pub source_dirty: Option<bool>,
    pub hash_id: String,
    pub hash_zero_policy: String,
    pub table_id: String,
    pub table_file: String,
    pub table_size_bytes: u64,
    pub table_sha256: String,
    pub kmer_size: u8,
    pub base_fscale: u64,
    pub cms_width: usize,
    pub cms_depth: usize,
    pub alpha: f32,
    pub filter_mode: String,
    pub target_fscale: u64,
    pub negative_fscale: String,
    pub unseen_fscale: u64,
    pub positive_retention: f32,
    pub negative_retention: f32,
    pub minimum_positive_retention: f32,
    pub positive_files: Vec<FileIdentity>,
    pub chromosome_background_files: Vec<FileIdentity>,
    pub creation_command: Vec<String>,
    pub creation_time_unix_seconds: u64,
}

pub fn sidecar_path(path: &Path) -> PathBuf {
    let mut value = OsString::from(path.as_os_str());
    value.push(".json");
    PathBuf::from(value)
}

pub fn sha256_file(path: &Path) -> Result<String> {
    let file = File::open(path)
        .with_context(|| format!("failed to open {} for checksum", path.display()))?;
    let mut reader = BufReader::with_capacity(1024 * 1024, file);
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 1024 * 1024];
    loop {
        let read = reader.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

pub fn sha256_bytes(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

pub fn file_identity(path: &Path) -> Result<FileIdentity> {
    Ok(FileIdentity {
        path: path.display().to_string(),
        size_bytes: path.metadata()?.len(),
        sha256: sha256_file(path)?,
    })
}

pub fn file_identities(paths: &[PathBuf]) -> Result<Vec<FileIdentity>> {
    paths.iter().map(|path| file_identity(path)).collect()
}

pub fn identities_checksum(identities: &[FileIdentity]) -> Result<String> {
    Ok(sha256_bytes(&serde_json::to_vec(identities)?))
}

pub fn write_json<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    let writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    serde_json::to_writer_pretty(writer, value)?;
    Ok(())
}

pub fn load_database_manifest(database_path: &Path) -> Result<Option<DatabaseManifest>> {
    let path = sidecar_path(database_path);
    if !path.exists() {
        return Ok(None);
    }
    let file = File::open(&path)?;
    let manifest = serde_json::from_reader(BufReader::new(file))
        .with_context(|| format!("failed to parse database manifest {}", path.display()))?;
    Ok(Some(manifest))
}

pub fn source_commit() -> String {
    if let Some(commit) = option_env!("JAM_RS_SOURCE_COMMIT") {
        return commit.to_string();
    }
    std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .ok()
        .filter(|output| output.status.success())
        .and_then(|output| String::from_utf8(output.stdout).ok())
        .map(|commit| commit.trim().to_string())
        .filter(|commit| !commit.is_empty())
        .unwrap_or_else(|| "unknown".to_string())
}

pub fn source_dirty() -> Option<bool> {
    let output = std::process::Command::new("git")
        .args(["status", "--porcelain", "--untracked-files=no"])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .ok()?;
    output.status.success().then_some(!output.stdout.is_empty())
}

pub fn unix_time_seconds() -> u64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs()
}

pub fn command_line() -> Vec<String> {
    std::env::args().collect()
}
