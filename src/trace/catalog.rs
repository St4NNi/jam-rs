//! Catalog records for one-plasmid-to-many-metagenome trace searches.
//!
//! The catalog is deliberately small and boring: it maps the stable sample
//! identifier stored in a `.jam` index to one self-contained JMA resource.
//! It is not a second sequence database. JSON and TSV
//! are accepted so a pipeline can keep the catalog under ordinary version
//! control without introducing another database format.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use thiserror::Error;

/// One row in a trace resource catalog.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CatalogEntry {
    /// Identifier used as the sample name in the existing `.jam` database.
    #[serde(alias = "id", alias = "sample", alias = "sample_id")]
    pub metagenome_id: String,
    /// URI or local path of the single self-contained JMA object.
    pub resource_uri: String,
    /// SHA-256 of the complete JMA object.
    pub sha256: String,
    #[serde(default)]
    pub etag: Option<String>,
    #[serde(default)]
    pub object_version: Option<String>,
    #[serde(default)]
    pub label: Option<String>,
    #[serde(default)]
    pub source_assembly_sha256: Option<String>,
}

impl CatalogEntry {
    #[must_use]
    pub fn resource(&self) -> &str {
        &self.resource_uri
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(untagged)]
enum CatalogDocument {
    Rows(Vec<CatalogEntry>),
    Wrapped { entries: Vec<CatalogEntry> },
}

/// A validated, deterministic catalog indexed by metagenome identifier.
#[derive(Clone, Debug, Default)]
pub struct TraceCatalog {
    source: PathBuf,
    entries: Vec<CatalogEntry>,
    index: HashMap<String, usize>,
}

/// Short compatibility alias used by CLI and pipeline callers.
pub type Catalog = TraceCatalog;

impl TraceCatalog {
    /// Read a JSON or tab-separated catalog from disk.
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, CatalogError> {
        let path = path.as_ref();
        let bytes = fs::read(path).map_err(|source| CatalogError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        let text = std::str::from_utf8(&bytes).map_err(|_| CatalogError::Invalid {
            path: path.to_path_buf(),
            message: "catalog is not valid UTF-8".to_string(),
        })?;
        let mut entries = if matches!(text.trim_start().chars().next(), Some('[' | '{')) {
            parse_json(text, path)?
        } else {
            parse_tsv(text, path)?
        };
        let base = path.parent().unwrap_or_else(|| Path::new("."));
        for entry in &mut entries {
            validate_entry(entry, path)?;
            entry.resource_uri = resolve_resource(base, std::mem::take(&mut entry.resource_uri));
            entry.sha256.make_ascii_lowercase();
            if let Some(source_sha256) = &mut entry.source_assembly_sha256 {
                source_sha256.make_ascii_lowercase();
            }
        }
        let mut index = HashMap::with_capacity(entries.len());
        for (position, entry) in entries.iter().enumerate() {
            if index
                .insert(entry.metagenome_id.clone(), position)
                .is_some()
            {
                return Err(CatalogError::Invalid {
                    path: path.to_path_buf(),
                    message: format!("duplicate metagenome identifier {:?}", entry.metagenome_id),
                });
            }
        }
        Ok(Self {
            source: path.to_path_buf(),
            entries,
            index,
        })
    }

    /// Construct a catalog from already-resolved rows.  This is useful for
    /// callers that receive catalog data from a manifest service.
    pub fn from_entries(mut entries: Vec<CatalogEntry>) -> Result<Self, CatalogError> {
        let mut index = HashMap::with_capacity(entries.len());
        for (position, entry) in entries.iter_mut().enumerate() {
            validate_entry(entry, Path::new("<memory>"))?;
            entry.sha256.make_ascii_lowercase();
            if let Some(source_sha256) = &mut entry.source_assembly_sha256 {
                source_sha256.make_ascii_lowercase();
            }
            if index
                .insert(entry.metagenome_id.clone(), position)
                .is_some()
            {
                return Err(CatalogError::Invalid {
                    path: PathBuf::from("<memory>"),
                    message: format!("duplicate metagenome identifier {:?}", entry.metagenome_id),
                });
            }
        }
        Ok(Self {
            source: PathBuf::from("<memory>"),
            entries,
            index,
        })
    }

    #[must_use]
    pub fn source(&self) -> &Path {
        &self.source
    }

    #[must_use]
    pub fn entries(&self) -> &[CatalogEntry] {
        &self.entries
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    #[must_use]
    pub fn get(&self, metagenome_id: &str) -> Option<&CatalogEntry> {
        self.index
            .get(metagenome_id)
            .and_then(|position| self.entries.get(*position))
    }

    /// Look up every requested candidate in one pass while retaining the
    /// caller's deterministic candidate order.
    pub fn resolve<I>(&self, ids: I) -> Result<Vec<&CatalogEntry>, CatalogError>
    where
        I: IntoIterator,
        I::Item: AsRef<str>,
    {
        ids.into_iter()
            .map(|id| {
                self.get(id.as_ref()).ok_or_else(|| CatalogError::Missing {
                    metagenome_id: id.as_ref().to_string(),
                })
            })
            .collect()
    }
}

fn parse_json(text: &str, path: &Path) -> Result<Vec<CatalogEntry>, CatalogError> {
    let document: CatalogDocument =
        serde_json::from_str(text).map_err(|source| CatalogError::Json {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(match document {
        CatalogDocument::Rows(rows) => rows,
        CatalogDocument::Wrapped { entries } => entries,
    })
}

fn parse_tsv(text: &str, path: &Path) -> Result<Vec<CatalogEntry>, CatalogError> {
    let mut lines = text
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.trim().is_empty());
    let Some((header_line, header)) = lines.next() else {
        return Err(CatalogError::Invalid {
            path: path.to_path_buf(),
            message: "catalog is empty".to_string(),
        });
    };
    let columns = header
        .split('\t')
        .map(|column| column.trim().to_ascii_lowercase())
        .collect::<Vec<_>>();
    let id_column = find_column(&columns, &["metagenome_id", "id", "sample_id", "sample"])
        .ok_or_else(|| CatalogError::Invalid {
            path: path.to_path_buf(),
            message: "TSV catalog requires a metagenome_id column".to_string(),
        })?;
    let resource_column =
        find_column(&columns, &["resource_uri"]).ok_or_else(|| CatalogError::Invalid {
            path: path.to_path_buf(),
            message: "TSV catalog requires a resource_uri column".to_string(),
        })?;
    let sha256_column =
        find_column(&columns, &["sha256"]).ok_or_else(|| CatalogError::Invalid {
            path: path.to_path_buf(),
            message: "TSV catalog requires a sha256 column".to_string(),
        })?;
    let etag_column = find_column(&columns, &["etag"]);
    let object_version_column = find_column(&columns, &["object_version"]);
    let label_column = find_column(&columns, &["label"]);
    let source_sha256_column = find_column(&columns, &["source_assembly_sha256"]);
    let allowed = [
        "metagenome_id",
        "resource_uri",
        "sha256",
        "etag",
        "object_version",
        "label",
        "source_assembly_sha256",
    ];
    if let Some(column) = columns
        .iter()
        .find(|column| !allowed.contains(&column.as_str()))
    {
        return Err(CatalogError::Invalid {
            path: path.to_path_buf(),
            message: format!("unsupported TSV catalog column {column:?}"),
        });
    }
    let mut rows = Vec::new();
    for (line_number, line) in lines {
        let fields = line.split('\t').collect::<Vec<_>>();
        let get = |column: Option<usize>| {
            column
                .and_then(|position| fields.get(position))
                .map(|value| value.trim())
                .filter(|value| !value.is_empty())
                .map(str::to_string)
        };
        let metagenome_id = get(Some(id_column)).ok_or_else(|| CatalogError::Invalid {
            path: path.to_path_buf(),
            message: format!("line {} has an empty metagenome_id", line_number + 1),
        })?;
        rows.push(CatalogEntry {
            metagenome_id,
            resource_uri: get(Some(resource_column)).ok_or_else(|| CatalogError::Invalid {
                path: path.to_path_buf(),
                message: format!("line {} has an empty resource_uri", line_number + 1),
            })?,
            sha256: get(Some(sha256_column)).ok_or_else(|| CatalogError::Invalid {
                path: path.to_path_buf(),
                message: format!("line {} has an empty sha256", line_number + 1),
            })?,
            etag: get(etag_column),
            object_version: get(object_version_column),
            label: get(label_column),
            source_assembly_sha256: get(source_sha256_column),
        });
    }
    if rows.is_empty() && header_line == 0 {
        return Err(CatalogError::Invalid {
            path: path.to_path_buf(),
            message: "TSV catalog has a header but no rows".to_string(),
        });
    }
    Ok(rows)
}

fn find_column(columns: &[String], names: &[&str]) -> Option<usize> {
    columns
        .iter()
        .position(|column| names.iter().any(|name| column == name))
}

fn validate_entry(entry: &CatalogEntry, path: &Path) -> Result<(), CatalogError> {
    if entry.metagenome_id.trim().is_empty() {
        return Err(CatalogError::Invalid {
            path: path.to_path_buf(),
            message: "metagenome_id must not be empty".to_string(),
        });
    }
    if entry.resource_uri.trim().is_empty() {
        return Err(CatalogError::Invalid {
            path: path.to_path_buf(),
            message: format!("catalog row {:?} has no resource_uri", entry.metagenome_id),
        });
    }
    if !valid_sha256(&entry.sha256) {
        return Err(CatalogError::Invalid {
            path: path.to_path_buf(),
            message: format!(
                "catalog row {:?} has an invalid sha256",
                entry.metagenome_id
            ),
        });
    }
    if entry
        .source_assembly_sha256
        .as_deref()
        .is_some_and(|checksum| !valid_sha256(checksum))
    {
        return Err(CatalogError::Invalid {
            path: path.to_path_buf(),
            message: format!(
                "catalog row {:?} has an invalid source_assembly_sha256",
                entry.metagenome_id
            ),
        });
    }
    Ok(())
}

fn valid_sha256(value: &str) -> bool {
    value.len() == 64 && value.bytes().all(|byte| byte.is_ascii_hexdigit())
}

fn resolve_resource(base: &Path, resource: String) -> String {
    if resource.contains("://") || Path::new(&resource).is_absolute() {
        resource
    } else {
        base.join(resource).to_string_lossy().into_owned()
    }
}

#[derive(Debug, Error)]
pub enum CatalogError {
    #[error("catalog I/O failed for {path}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    #[error("catalog JSON parse failed for {path}: {source}")]
    Json {
        path: PathBuf,
        source: serde_json::Error,
    },
    #[error("invalid catalog {path}: {message}")]
    Invalid { path: PathBuf, message: String },
    #[error("candidate metagenome is absent from catalog: {metagenome_id}")]
    Missing { metagenome_id: String },
}
