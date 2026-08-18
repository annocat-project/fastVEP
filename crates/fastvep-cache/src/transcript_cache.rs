//! Binary transcript cache for fast startup.
//!
//! Serializes fully-built `Vec<Transcript>` (including spliced sequences)
//! to a compact binary format using bincode + zstd compression.
//! Subsequent loads skip GFF3 parsing, FASTA loading, and sequence building.

use anyhow::{Context, Result};
use fastvep_genome::Transcript;
use serde::Serialize;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;
use std::time::SystemTime;

/// Magic header for zstd-compressed caches (current format).
const CACHE_MAGIC_V2: &[u8; 8] = b"FSTVEP02";
/// Magic header for legacy gzip-compressed caches (read-only support).
const CACHE_MAGIC_V1: &[u8; 8] = b"FSTVEP01";

/// Save transcripts to a binary cache file (bincode + zstd).
pub fn save_cache(transcripts: &[Transcript], path: &Path) -> Result<()> {
    let dir = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let tmp = tempfile::Builder::new()
        .prefix(".fastvep-cache-")
        .suffix(".partial")
        .tempfile_in(dir)
        .with_context(|| format!("Creating temporary cache file in {}", dir.display()))?;
    let writer = BufWriter::new(
        tmp.reopen()
            .with_context(|| format!("Reopening temporary cache file {}", tmp.path().display()))?,
    );
    // zstd level 1: fast compression, still much better decompression than gzip
    let mut zst = zstd::Encoder::new(writer, 1)?;
    zst.include_checksum(true)?;

    // Write magic header
    use std::io::Write;
    zst.write_all(CACHE_MAGIC_V2)?;

    // Serialize with bincode
    bincode::serialize_into(&mut zst, transcripts)
        .with_context(|| "Serializing transcripts to cache")?;

    let mut writer = zst.finish()?;
    writer.flush().with_context(|| "Flushing cache writer")?;
    writer
        .get_ref()
        .sync_all()
        .with_context(|| "Syncing cache file to disk")?;
    drop(writer);

    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mode = std::fs::metadata(path)
            .map(|metadata| metadata.permissions().mode() & 0o7777)
            .unwrap_or(0o644);
        std::fs::set_permissions(tmp.path(), std::fs::Permissions::from_mode(mode))
            .with_context(|| format!("Setting cache permissions for {}", path.display()))?;
    }

    tmp.persist(path).map_err(|error| {
        anyhow::anyhow!("Publishing cache file {}: {}", path.display(), error.error)
    })?;
    Ok(())
}

/// Load transcripts from a binary cache file.
/// Supports both zstd (v2) and legacy gzip (v1) formats.
pub fn load_cache(path: &Path) -> Result<Vec<Transcript>> {
    let file =
        File::open(path).with_context(|| format!("Opening cache file: {}", path.display()))?;
    let mut reader = BufReader::new(file);

    // Peek at the first bytes to detect format.
    // zstd frames start with 0x28B52FFD; gzip starts with 0x1F8B.
    use std::io::Read;
    let mut peek = [0u8; 4];
    reader
        .read_exact(&mut peek)
        .with_context(|| "Reading cache header")?;

    // Rewind so the decompressor sees the full stream
    use std::io::Seek;
    reader.seek(std::io::SeekFrom::Start(0))?;

    if peek[0..2] == [0x1F, 0x8B] {
        // Legacy gzip format (v1)
        load_cache_gzip(reader)
    } else {
        // zstd format (v2, or future)
        load_cache_zstd(reader)
    }
}

fn load_cache_zstd<R: std::io::Read>(reader: R) -> Result<Vec<Transcript>> {
    let mut zst = zstd::Decoder::new(reader)?;

    use std::io::Read;
    let mut magic = [0u8; 8];
    zst.read_exact(&mut magic)
        .with_context(|| "Reading cache header")?;
    if &magic != CACHE_MAGIC_V2 {
        anyhow::bail!("Invalid cache file (wrong magic header, expected FSTVEP02)");
    }

    let transcripts: Vec<Transcript> = bincode::deserialize_from(&mut zst)
        .with_context(|| "Deserializing transcripts from cache")?;
    require_decompressed_eof(&mut zst)?;
    Ok(transcripts)
}

fn load_cache_gzip<R: std::io::Read>(reader: R) -> Result<Vec<Transcript>> {
    use flate2::read::GzDecoder;

    let mut gz = GzDecoder::new(reader);

    use std::io::Read;
    let mut magic = [0u8; 8];
    gz.read_exact(&mut magic)
        .with_context(|| "Reading cache header")?;
    if &magic != CACHE_MAGIC_V1 {
        anyhow::bail!("Invalid cache file (wrong magic header, expected FSTVEP01)");
    }

    let transcripts: Vec<Transcript> = bincode::deserialize_from(&mut gz)
        .with_context(|| "Deserializing transcripts from cache")?;
    require_decompressed_eof(&mut gz)?;
    Ok(transcripts)
}

fn require_decompressed_eof(reader: &mut impl std::io::Read) -> Result<()> {
    let mut trailing = [0u8; 1];
    if reader
        .read(&mut trailing)
        .with_context(|| "Finishing transcript cache stream")?
        != 0
    {
        anyhow::bail!("Transcript cache contains trailing decoded data");
    }
    Ok(())
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct TranscriptCacheVerification {
    pub schema_version: u32,
    pub cache_format: &'static str,
    pub cache_bytes: u64,
    pub transcript_count: u64,
    pub coding_transcript_count: u64,
    pub coding_with_sequence_count: u64,
    pub primary_coding_missing_sequence_count: u64,
    pub non_primary_coding_missing_sequence_count: u64,
}

/// Fully decode a transcript cache and validate the structural invariants used
/// by annotation. Installed caches may require sequence-complete coding
/// transcripts on the primary assembly while still permitting annotation
/// records on alt contigs that are absent from a primary-only FASTA.
pub fn verify_cache(
    path: &Path,
    require_primary_coding_sequences: bool,
) -> Result<TranscriptCacheVerification> {
    let cache_bytes = path
        .metadata()
        .with_context(|| format!("Reading cache metadata: {}", path.display()))?
        .len();
    if cache_bytes == 0 {
        anyhow::bail!("Transcript cache is empty");
    }

    let cache_format = cache_format(path)?;
    let transcripts = load_cache(path)?;
    if transcripts.is_empty() {
        anyhow::bail!("Transcript cache contains no transcripts");
    }

    let mut coding_transcript_count = 0u64;
    let mut coding_with_sequence_count = 0u64;
    let mut primary_coding_missing_sequence_count = 0u64;
    let mut non_primary_coding_missing_sequence_count = 0u64;

    for transcript in &transcripts {
        if transcript.stable_id.is_empty()
            || transcript.chromosome.is_empty()
            || transcript.start == 0
            || transcript.end < transcript.start
            || transcript.exons.is_empty()
            || transcript
                .exons
                .iter()
                .any(|exon| exon.start == 0 || exon.end < exon.start)
        {
            anyhow::bail!(
                "Transcript cache contains an invalid transcript record: {}",
                transcript.stable_id
            );
        }

        if transcript.is_coding() {
            coding_transcript_count += 1;
            if transcript
                .spliced_seq
                .as_ref()
                .is_some_and(|value| !value.is_empty())
            {
                coding_with_sequence_count += 1;
            } else if is_primary_chromosome(&transcript.chromosome) {
                primary_coding_missing_sequence_count += 1;
            } else {
                non_primary_coding_missing_sequence_count += 1;
            }
        }
    }

    if require_primary_coding_sequences && primary_coding_missing_sequence_count != 0 {
        anyhow::bail!(
            "Transcript cache has {} primary-assembly coding transcripts without sequence data",
            primary_coding_missing_sequence_count
        );
    }

    Ok(TranscriptCacheVerification {
        schema_version: 1,
        cache_format,
        cache_bytes,
        transcript_count: transcripts.len() as u64,
        coding_transcript_count,
        coding_with_sequence_count,
        primary_coding_missing_sequence_count,
        non_primary_coding_missing_sequence_count,
    })
}

fn cache_format(path: &Path) -> Result<&'static str> {
    use std::io::Read;
    let mut file =
        File::open(path).with_context(|| format!("Opening cache file: {}", path.display()))?;
    let mut prefix = [0u8; 4];
    file.read_exact(&mut prefix)
        .with_context(|| "Reading cache compression header")?;
    Ok(if prefix[0..2] == [0x1F, 0x8B] {
        "FSTVEP01"
    } else {
        "FSTVEP02"
    })
}

fn is_primary_chromosome(chromosome: &str) -> bool {
    let chromosome = chromosome
        .strip_prefix("chr")
        .or_else(|| chromosome.strip_prefix("CHR"))
        .unwrap_or(chromosome);
    matches!(chromosome, "X" | "Y" | "M" | "MT")
        || chromosome
            .parse::<u8>()
            .is_ok_and(|value| (1..=22).contains(&value))
}

/// Check if cache file is newer than source file.
pub fn cache_is_fresh(cache_path: &Path, source_path: &Path) -> bool {
    let cache_mtime = cache_path
        .metadata()
        .and_then(|m| m.modified())
        .unwrap_or(SystemTime::UNIX_EPOCH);
    let source_mtime = source_path
        .metadata()
        .and_then(|m| m.modified())
        .unwrap_or(SystemTime::now());
    cache_mtime > source_mtime
}

/// Get the default cache path for a given GFF3 path.
pub fn default_cache_path(gff3_path: &Path) -> std::path::PathBuf {
    let mut cache_path = gff3_path.to_path_buf();
    let name = cache_path
        .file_name()
        .map(|n| {
            let s = n.to_string_lossy();
            if s.ends_with(".fastvep.cache") {
                s.to_string()
            } else {
                format!("{}.fastvep.cache", s)
            }
        })
        .unwrap_or_else(|| "transcripts.fastvep.cache".to_string());
    cache_path.set_file_name(name);
    cache_path
}

#[cfg(test)]
mod tests {
    use super::*;
    use fastvep_core::Strand;
    use fastvep_genome::{Exon, Gene, Transcript, Translation};
    use std::sync::Arc;
    use tempfile::NamedTempFile;

    fn make_test_transcript() -> Transcript {
        Transcript {
            stable_id: Arc::from("ENST00000001"),
            version: Some(1),
            gene: Gene {
                stable_id: Arc::from("ENSG00000001"),
                symbol: Some(Arc::from("TEST")),
                symbol_source: None,
                hgnc_id: None,
                biotype: Arc::from("protein_coding"),
                chromosome: Arc::from("1"),
                start: 1000,
                end: 5000,
                strand: Strand::Forward,
            },
            biotype: Arc::from("protein_coding"),
            chromosome: Arc::from("1"),
            start: 1000,
            end: 5000,
            strand: Strand::Forward,
            exons: vec![Exon {
                stable_id: "ENSE001".into(),
                start: 1000,
                end: 1200,
                strand: Strand::Forward,
                phase: 0,
                end_phase: -1,
                rank: 1,
            }],
            translation: None,
            cdna_coding_start: Some(1),
            cdna_coding_end: Some(200),
            coding_region_start: Some(1000),
            coding_region_end: Some(1200),
            spliced_seq: Some("ACGTACGT".into()),
            translateable_seq: Some("ACGT".into()),
            peptide: Some("T".into()),
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: Some(1),
            appris: Some("P1".into()),
            ccds: None,
            protein_id: Some("ENSP001".into()),
            protein_version: Some(1),
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
        }
    }

    #[test]
    fn test_cache_roundtrip() {
        let transcripts = vec![make_test_transcript()];
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("transcripts.cache");

        save_cache(&transcripts, &path).unwrap();
        let loaded = load_cache(&path).unwrap();

        assert_eq!(loaded.len(), 1);
        assert_eq!(&*loaded[0].stable_id, "ENST00000001");
        assert_eq!(&**loaded[0].gene.symbol.as_ref().unwrap(), "TEST");
        assert_eq!(loaded[0].spliced_seq.as_deref(), Some("ACGTACGT"));
        assert_eq!(loaded[0].canonical, true);
        assert_eq!(loaded[0].tsl, Some(1));
    }

    #[test]
    fn test_invalid_magic() {
        let tmp = NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), b"NOTVALID").unwrap();
        assert!(load_cache(tmp.path()).is_err());
    }

    #[test]
    fn test_legacy_gzip_cache_loads() {
        // Create a legacy gzip cache and verify it still loads
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;

        let transcripts = vec![make_test_transcript()];
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();

        let file = File::create(path).unwrap();
        let writer = BufWriter::new(file);
        let mut gz = GzEncoder::new(writer, Compression::fast());
        gz.write_all(CACHE_MAGIC_V1).unwrap();
        bincode::serialize_into(&mut gz, &transcripts).unwrap();
        gz.finish().unwrap();

        let loaded = load_cache(path).unwrap();
        assert_eq!(loaded.len(), 1);
        assert_eq!(&*loaded[0].stable_id, "ENST00000001");
    }

    fn make_coding_transcript(chromosome: &str, sequence: bool) -> Transcript {
        let mut transcript = make_test_transcript();
        transcript.chromosome = chromosome.into();
        transcript.gene.chromosome = chromosome.into();
        transcript.translation = Some(Translation {
            stable_id: "ENSP00000001".into(),
            genomic_start: 1000,
            genomic_end: 1200,
            start_exon_rank: 1,
            start_exon_offset: 0,
            end_exon_rank: 1,
            end_exon_offset: 200,
        });
        if !sequence {
            transcript.spliced_seq = None;
            transcript.translateable_seq = None;
            transcript.peptide = None;
        }
        transcript
    }

    #[test]
    fn verifies_a_sequence_complete_primary_cache() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("transcripts.cache");
        save_cache(&[make_coding_transcript("1", true)], &path).unwrap();

        let report = verify_cache(&path, true).unwrap();
        assert_eq!(report.cache_format, "FSTVEP02");
        assert_eq!(report.transcript_count, 1);
        assert_eq!(report.coding_with_sequence_count, 1);
        assert_eq!(report.primary_coding_missing_sequence_count, 0);
    }

    #[test]
    fn rejects_missing_primary_coding_sequences_when_required() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("transcripts.cache");
        save_cache(&[make_coding_transcript("chr1", false)], &path).unwrap();

        assert!(verify_cache(&path, true).is_err());
        let report = verify_cache(&path, false).unwrap();
        assert_eq!(report.primary_coding_missing_sequence_count, 1);
    }

    #[test]
    fn permits_missing_non_primary_sequences() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("transcripts.cache");
        save_cache(&[make_coding_transcript("KI270713.1", false)], &path).unwrap();

        let report = verify_cache(&path, true).unwrap();
        assert_eq!(report.non_primary_coding_missing_sequence_count, 1);
    }

    #[test]
    fn replaces_an_existing_cache_without_leaving_a_partial_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("transcripts.cache");
        save_cache(&[make_test_transcript()], &path).unwrap();

        let mut replacement = make_test_transcript();
        replacement.stable_id = "ENST00000002".into();
        save_cache(&[replacement], &path).unwrap();

        let loaded = load_cache(&path).unwrap();
        assert_eq!(&*loaded[0].stable_id, "ENST00000002");
        assert!(std::fs::read_dir(dir.path()).unwrap().all(|entry| {
            !entry
                .unwrap()
                .file_name()
                .to_string_lossy()
                .ends_with(".partial")
        }));
    }
}
