//! Manifest-backed collection of chromosome-specific OSA databases.

use crate::common::chrom_aliases;
use crate::reader::{AnySaReader, SaCacheFormat};
use anyhow::{Context, Result};
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue, SaMetadata};
use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;

pub const SHARD_MANIFEST_SUFFIX: &str = ".osa-shards.json";

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
struct ShardManifest {
    schema_version: u16,
    shards: Vec<ShardEntry>,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
struct ShardEntry {
    chromosome: String,
    file: String,
}

/// One logical annotation provider backed by one verified OSA per chromosome.
pub struct ShardedSaReader {
    metadata: SaMetadata,
    readers: Vec<AnySaReader>,
    chromosome_to_reader: HashMap<String, usize>,
    fallback_reader: Option<usize>,
}

impl ShardedSaReader {
    pub fn open(manifest_path: &Path) -> Result<Self> {
        let manifest_bytes = std::fs::read(manifest_path)
            .with_context(|| format!("Reading OSA shard manifest {}", manifest_path.display()))?;
        let manifest: ShardManifest = serde_json::from_slice(&manifest_bytes)
            .with_context(|| format!("Parsing OSA shard manifest {}", manifest_path.display()))?;
        if manifest.schema_version != 1 {
            anyhow::bail!(
                "Unsupported OSA shard manifest schemaVersion {} in {}",
                manifest.schema_version,
                manifest_path.display()
            );
        }
        if manifest.shards.is_empty() {
            anyhow::bail!(
                "OSA shard manifest contains no shards: {}",
                manifest_path.display()
            );
        }

        let manifest_dir = manifest_path
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .canonicalize()
            .with_context(|| {
                format!("Resolving shard directory for {}", manifest_path.display())
            })?;
        let mut readers = Vec::with_capacity(manifest.shards.len());
        let mut chromosome_to_reader = HashMap::new();
        let mut fallback_reader = None;
        let mut metadata: Option<SaMetadata> = None;
        let mut format: Option<SaCacheFormat> = None;

        for shard in manifest.shards {
            let relative = Path::new(&shard.file);
            if relative.is_absolute() {
                anyhow::bail!("OSA shard path must be relative: {}", shard.file);
            }
            let shard_path = manifest_dir
                .join(relative)
                .canonicalize()
                .with_context(|| {
                    format!(
                        "Resolving OSA shard '{}' from {}",
                        shard.file,
                        manifest_path.display()
                    )
                })?;
            if !shard_path.starts_with(&manifest_dir) {
                anyhow::bail!("OSA shard escapes its manifest directory: {}", shard.file);
            }

            let reader = AnySaReader::open(&shard_path)
                .with_context(|| format!("Opening OSA shard {}", shard_path.display()))?;
            if let Some(expected) = format {
                if reader.format() != expected {
                    anyhow::bail!("OSA shard manifest mixes OSA v1 and OSA2 files");
                }
            } else {
                format = Some(reader.format());
            }
            if let Some(expected) = &metadata {
                validate_same_source(expected, reader.metadata(), &shard.file)?;
            } else {
                metadata = Some(reader.metadata().clone());
            }

            let reader_index = readers.len();
            if shard.chromosome.eq_ignore_ascii_case("all")
                && fallback_reader.replace(reader_index).is_some()
            {
                anyhow::bail!("Duplicate all-chromosome fallback in OSA shard manifest");
            }
            for alias in chrom_aliases(&shard.chromosome) {
                if chromosome_to_reader
                    .insert(alias.clone(), reader_index)
                    .is_some()
                {
                    anyhow::bail!(
                        "Duplicate chromosome/alias '{}' in OSA shard manifest",
                        alias
                    );
                }
            }
            readers.push(reader);
        }

        Ok(Self {
            metadata: metadata.expect("non-empty manifest establishes metadata"),
            readers,
            chromosome_to_reader,
            fallback_reader,
        })
    }

    fn reader_for(&self, chromosome: &str) -> Option<&AnySaReader> {
        let index = chrom_aliases(chromosome)
            .iter()
            .find_map(|alias| self.chromosome_to_reader.get(alias))
            .copied()
            .or(self.fallback_reader)?;
        self.readers.get(index)
    }
}

fn validate_same_source(expected: &SaMetadata, actual: &SaMetadata, file: &str) -> Result<()> {
    if expected.name != actual.name
        || expected.version != actual.version
        || expected.assembly != actual.assembly
        || expected.json_key != actual.json_key
        || expected.match_by_allele != actual.match_by_allele
        || expected.is_array != actual.is_array
        || expected.record_list != actual.record_list
        || expected.is_positional != actual.is_positional
    {
        anyhow::bail!(
            "OSA shard '{}' has source metadata inconsistent with the first shard",
            file
        );
    }
    Ok(())
}

impl AnnotationProvider for ShardedSaReader {
    fn name(&self) -> &str {
        &self.metadata.name
    }
    fn json_key(&self) -> &str {
        &self.metadata.json_key
    }
    fn metadata(&self) -> &SaMetadata {
        &self.metadata
    }

    fn cache_load_count(&self) -> Option<u64> {
        let mut total = 0_u64;
        for reader in &self.readers {
            total = total.checked_add(reader.cache_load_count()?)?;
        }
        Some(total)
    }

    fn annotate_position(
        &self,
        chrom: &str,
        pos: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Option<AnnotationValue>> {
        match self.reader_for(chrom) {
            Some(reader) => reader.annotate_position(chrom, pos, ref_allele, alt_allele),
            None => Ok(None),
        }
    }

    fn preload(&self, chrom: &str, positions: &[u64]) -> Result<()> {
        match self.reader_for(chrom) {
            Some(reader) => reader.preload(chrom, positions),
            None => Ok(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{AnnotationRecord, SCHEMA_VERSION};
    use crate::index::IndexHeader;
    use crate::writer::SaWriter;
    use crate::writer_v2::{Osa2Metadata, Osa2Record, Osa2Writer};

    fn write_shard(dir: &Path, name: &str, chromosome: &str, chrom_idx: u16, position: u32) {
        let header = IndexHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "dbnsfp".into(),
            name: "dbNSFP".into(),
            version: "4.9a".into(),
            description: "test".into(),
            assembly: "GRCh38".into(),
            match_by_allele: true,
            is_array: false,
            record_list: false,
            is_positional: false,
        };
        let base = dir.join(name);
        let mut writer = SaWriter::new(header);
        writer
            .write_to_files(
                &base,
                std::iter::once(AnnotationRecord {
                    chrom_idx,
                    position,
                    ref_allele: "A".into(),
                    alt_allele: "G".into(),
                    json: format!(r#"{{"chromosome":"{}"}}"#, chromosome),
                }),
                &[chromosome.to_string()],
            )
            .unwrap();
    }

    fn write_osa2_shard(dir: &Path, name: &str, chromosome: &str, position: u32) {
        let metadata = Osa2Metadata {
            format_version: 2,
            name: "dbNSFP".into(),
            version: "4.9a".into(),
            description: "test".into(),
            assembly: "GRCh38".into(),
            json_key: "dbnsfp".into(),
            match_by_allele: true,
            is_array: false,
            record_list: false,
            is_positional: false,
            chunk_bits: 20,
        };
        Osa2Writer::new(metadata, Vec::new())
            .write_all(
                std::fs::File::create(dir.join(name)).unwrap(),
                &[Osa2Record {
                    chrom: chromosome.into(),
                    position,
                    ref_allele: b"A".to_vec(),
                    alt_allele: b"G".to_vec(),
                    values: Vec::new(),
                    json_blob: None,
                }],
            )
            .unwrap();
    }

    #[test]
    fn dispatches_each_chromosome_through_one_logical_provider() {
        let temp = tempfile::tempdir().unwrap();
        let shard_dir = temp.path().join("shards");
        std::fs::create_dir(&shard_dir).unwrap();
        write_shard(&shard_dir, "chr1", "chr1", 0, 100);
        write_shard(&shard_dir, "chr2", "chr2", 0, 200);
        let manifest = temp.path().join("dbnsfp.osa-shards.json");
        std::fs::write(
            &manifest,
            r#"{"schemaVersion":1,"shards":[{"chromosome":"1","file":"shards/chr1.osa"},{"chromosome":"2","file":"shards/chr2.osa"}]}"#,
        )
        .unwrap();

        let reader = ShardedSaReader::open(&manifest).unwrap();
        assert_eq!(reader.name(), "dbNSFP");
        assert_eq!(reader.json_key(), "dbnsfp");
        assert!(reader
            .annotate_position("1", 100, "A", "G")
            .unwrap()
            .is_some());
        assert!(reader
            .annotate_position("chr2", 200, "A", "G")
            .unwrap()
            .is_some());
        assert!(reader
            .annotate_position("chr3", 100, "A", "G")
            .unwrap()
            .is_none());
        assert!(reader
            .annotate_position("chr1", 200, "A", "G")
            .unwrap()
            .is_none());
    }

    #[test]
    fn uses_all_chromosome_shard_as_fallback() {
        let temp = tempfile::tempdir().unwrap();
        let shard_dir = temp.path().join("shards");
        std::fs::create_dir(&shard_dir).unwrap();
        write_shard(&shard_dir, "all", "chr1", 0, 100);
        let manifest = temp.path().join("clinvar.osa-shards.json");
        std::fs::write(
            &manifest,
            r#"{"schemaVersion":1,"shards":[{"chromosome":"all","file":"shards/all.osa"}]}"#,
        )
        .unwrap();

        let reader = ShardedSaReader::open(&manifest).unwrap();
        assert!(reader
            .annotate_position("1", 100, "A", "G")
            .unwrap()
            .is_some());
        assert!(reader
            .annotate_position("chr1", 100, "A", "G")
            .unwrap()
            .is_some());
        assert!(reader
            .annotate_position("2", 100, "A", "G")
            .unwrap()
            .is_none());
    }

    #[test]
    fn rejects_duplicate_chromosomes_and_directory_escape() {
        let temp = tempfile::tempdir().unwrap();
        let shard_dir = temp.path().join("shards");
        std::fs::create_dir(&shard_dir).unwrap();
        write_shard(&shard_dir, "chr1", "chr1", 0, 100);

        let duplicate = temp.path().join("duplicate.osa-shards.json");
        std::fs::write(
            &duplicate,
            r#"{"schemaVersion":1,"shards":[{"chromosome":"1","file":"shards/chr1.osa"},{"chromosome":"chr1","file":"shards/chr1.osa"}]}"#,
        )
        .unwrap();
        assert!(ShardedSaReader::open(&duplicate)
            .err()
            .expect("duplicate chromosomes must fail")
            .to_string()
            .contains("Duplicate chromosome"));

        let outside_dir = tempfile::tempdir().unwrap();
        write_shard(outside_dir.path(), "outside", "chr1", 0, 100);
        let escape = temp.path().join("escape.osa-shards.json");
        let outside = outside_dir.path().join("outside.osa");
        std::fs::write(
            &escape,
            serde_json::json!({
                "schemaVersion": 1,
                "shards": [{"chromosome": "1", "file": outside.to_string_lossy()}]
            })
            .to_string(),
        )
        .unwrap();
        assert!(ShardedSaReader::open(&escape)
            .err()
            .expect("absolute shard path must fail")
            .to_string()
            .contains("relative"));
    }

    #[test]
    fn rejects_mixed_cache_formats() {
        let temp = tempfile::tempdir().unwrap();
        let shard_dir = temp.path().join("shards");
        std::fs::create_dir(&shard_dir).unwrap();
        write_shard(&shard_dir, "chr1", "1", 0, 100);
        write_osa2_shard(&shard_dir, "chr2.osa2", "2", 200);
        let manifest = temp.path().join("mixed.osa-shards.json");
        std::fs::write(
            &manifest,
            r#"{"schemaVersion":1,"shards":[{"chromosome":"1","file":"shards/chr1.osa"},{"chromosome":"2","file":"shards/chr2.osa2"}]}"#,
        )
        .unwrap();

        let error = ShardedSaReader::open(&manifest)
            .err()
            .expect("mixed cache formats must fail");
        assert!(error.to_string().contains("mixes OSA v1 and OSA2"));
    }
}
