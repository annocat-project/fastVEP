//! AnnoCAT v1: inspectable provenance followed by a checked zstd/bincode payload.
use anyhow::{ensure, Context, Result};
use fastvep_genome::Transcript;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    fs::File,
    io::{self, Read, Write},
    path::Path,
};

pub const MAGIC: &[u8; 8] = b"ANNOCATC";
pub const FORMAT: &str = "ANNOCATC1";
const VERSION: u32 = 1;
const MAX_HEADER: usize = 64 * 1024;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct SequenceEdit {
    pub kind: String,
    pub start: u64,
    pub end: u64,
    pub replacement: String,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Enrichment {
    pub gene_version: Option<u32>,
    pub mature_mirna_ranges: Vec<(u64, u64)>,
    pub rna_edits: Vec<SequenceEdit>,
    pub translation_edits: Vec<SequenceEdit>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Header {
    pub species: String,
    pub assembly: String,
    pub ensembl_release: u32,
    pub vep_release: String,
    pub capabilities: Vec<String>,
    pub provenance: serde_json::Value,
    pub transcript_count: u64,
    pub coding_transcript_count: u64,
    pub contig_counts: BTreeMap<String, u64>,
    pub payload_bytes: u64,
    pub payload_sha256: String,
    pub semantic_sha256: String,
}

#[derive(Serialize, Deserialize)]
struct PayloadV1 {
    transcripts: Vec<crate::transcript_wire::Transcript>,
    enrichment: BTreeMap<String, Enrichment>,
}

pub struct LoadedCache {
    pub transcripts: Vec<Transcript>,
    pub enrichment: BTreeMap<String, Enrichment>,
    pub header: Option<Header>,
}

pub fn sha256(path: &Path) -> Result<String> {
    let mut file = File::open(path)?;
    let mut digest = Sha256::new();
    let mut buffer = [0u8; 65536];
    loop {
        let n = file.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        digest.update(&buffer[..n]);
    }
    Ok(format!("{:x}", digest.finalize()))
}

struct DigestWriter(Sha256);
impl Write for DigestWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0.update(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

fn semantic_digest(payload: &PayloadV1) -> Result<String> {
    let mut hash = DigestWriter(Sha256::new());
    if std::env::var_os("ANNOCAT_BENCH_CACHE_BUFFERED").is_some() {
        let mut buffered = io::BufWriter::with_capacity(64 * 1024, &mut hash);
        bincode::serialize_into(&mut buffered, payload)?;
        buffered.flush()?;
    } else {
        bincode::serialize_into(&mut hash, payload)?;
    }
    Ok(format!("{:x}", hash.0.finalize()))
}

pub fn read_header(path: &Path) -> Result<Option<Header>> {
    let mut file = File::open(path)?;
    let mut magic = [0; 8];
    file.read_exact(&mut magic)
        .context("Reading transcript cache signature")?;
    if &magic != MAGIC {
        return Ok(None);
    }
    Ok(Some(decode_header(&mut file)?))
}

fn decode_header(file: &mut File) -> Result<Header> {
    let mut word = [0; 4];
    file.read_exact(&mut word)?;
    ensure!(
        u32::from_le_bytes(word) == VERSION,
        "Unsupported AnnoCAT cache version"
    );
    file.read_exact(&mut word)?;
    let length = u32::from_le_bytes(word) as usize;
    ensure!(
        length > 0 && length <= MAX_HEADER,
        "Invalid AnnoCAT cache header length"
    );
    let mut bytes = vec![0; length];
    file.read_exact(&mut bytes)?;
    let header: Header = serde_json::from_slice(&bytes)?;
    ensure!(
        header.transcript_count > 0 && header.payload_bytes > 0,
        "Empty AnnoCAT cache"
    );
    for hash in [&header.payload_sha256, &header.semantic_sha256] {
        ensure!(
            hash.len() == 64 && hash.bytes().all(|b| b.is_ascii_hexdigit()),
            "Invalid cache digest"
        );
    }
    ensure!(
        file.metadata()?.len() == 16 + length as u64 + header.payload_bytes,
        "AnnoCAT cache length does not match its header"
    );
    Ok(header)
}

pub fn load(path: &Path) -> Result<LoadedCache> {
    let mut file = File::open(path)?;
    let mut magic = [0; 8];
    file.read_exact(&mut magic)?;
    if &magic != MAGIC {
        return Ok(LoadedCache {
            transcripts: crate::transcript_cache::load_legacy_cache(path)?,
            enrichment: BTreeMap::new(),
            header: None,
        });
    }
    let header = decode_header(&mut file)?;
    use std::io::{Seek, SeekFrom};
    let start = file.stream_position()?;
    let mut digest = DigestWriter(Sha256::new());
    io::copy(&mut file, &mut digest)?;
    ensure!(
        format!("{:x}", digest.0.finalize()) == header.payload_sha256,
        "Cache payload checksum mismatch"
    );
    file.seek(SeekFrom::Start(start))?;
    let capacity = if std::env::var_os("ANNOCAT_BENCH_CACHE_BUFFERED").is_some() {64 * 1024} else {0};
    let mut decoder = io::BufReader::with_capacity(capacity, zstd::Decoder::new(file)?);
    let decode_started = std::time::Instant::now();
    use bincode::Options;
    let payload: PayloadV1 = bincode::DefaultOptions::new()
        .with_fixint_encoding()
        .with_limit(8 * 1024 * 1024 * 1024)
        .deserialize_from(&mut decoder)?;
    if std::env::var_os("ANNOCAT_BENCH_PROFILE").is_some() {eprintln!("cacheDecodeSeconds={}", decode_started.elapsed().as_secs_f64());}
    let mut trailing = [0];
    ensure!(
        decoder.read(&mut trailing)? == 0,
        "Trailing data in AnnoCAT cache payload"
    );
    ensure!(
        payload.transcripts.len() as u64 == header.transcript_count,
        "Cache transcript count mismatch"
    );
    let semantic_started = std::time::Instant::now();
    ensure!(
        semantic_digest(&payload)? == header.semantic_sha256,
        "Cache semantic digest mismatch"
    );
    if std::env::var_os("ANNOCAT_BENCH_PROFILE").is_some() {eprintln!("cacheSemanticDigestSeconds={}", semantic_started.elapsed().as_secs_f64());}
    let mut transcripts: Vec<Transcript> = payload.transcripts.into_iter().map(Into::into).collect();
    ensure!(
        transcripts.iter().filter(|tr| tr.is_coding()).count() as u64
            == header.coding_transcript_count,
        "Cache coding transcript count mismatch"
    );
    let mut contigs = BTreeMap::new();
    for tr in &transcripts {
        *contigs.entry(tr.chromosome.to_string()).or_insert(0u64) += 1;
    }
    ensure!(
        contigs == header.contig_counts,
        "Cache contig counts mismatch"
    );
    validate_enrichment(&transcripts, &payload.enrichment)?;
    for tr in &mut transcripts {
        if let Some(cds) = &tr.translateable_seq {
            let edits = payload.enrichment.get(&format!("{}:{}", tr.chromosome, tr.stable_id))
                .map(|extra| extra.translation_edits.as_slice()).unwrap_or_default();
            let reference = crate::ensembl_core::translated_peptide(cds, &tr.chromosome, edits, true)?;
            if tr.peptide.as_deref() != Some(reference.as_str()) {
                tr.reference_peptide = Some(reference);
            }
        }
    }
    Ok(LoadedCache {
        transcripts,
        enrichment: payload.enrichment,
        header: Some(header),
    })
}

pub fn validate_enrichment(
    transcripts: &[Transcript],
    enrichment: &BTreeMap<String, Enrichment>,
) -> Result<()> {
    let mut seen = std::collections::BTreeSet::new();
    for tr in transcripts {
        let key = format!("{}:{}", tr.chromosome, tr.stable_id);
        ensure!(seen.insert(key.clone()), "Duplicate cache transcript {key}");
        if let Some(data) = enrichment.get(&key) {
            ensure!(
                data.mature_mirna_ranges
                    .iter()
                    .all(|&(s, e)| tr.biotype.as_ref() == "miRNA"
                        && s > 0
                        && s <= e
                        && e <= tr.cdna_length()),
                "Invalid mature-miRNA coordinates for {key}"
            );
        }
    }
    ensure!(
        enrichment.keys().all(|k| seen.contains(k)),
        "Enrichment refers to an absent transcript"
    );
    Ok(())
}

/// Takes ownership to avoid cloning transcript sequences during serialization.
pub fn save(
    transcripts: Vec<Transcript>,
    enrichment: BTreeMap<String, Enrichment>,
    mut header: Header,
    path: &Path,
) -> Result<Header> {
    validate_enrichment(&transcripts, &enrichment)?;
    ensure!(!transcripts.is_empty(), "Refusing to write an empty cache");
    header.transcript_count = transcripts.len() as u64;
    header.coding_transcript_count = transcripts.iter().filter(|t| t.is_coding()).count() as u64;
    header.contig_counts.clear();
    for tr in &transcripts {
        *header
            .contig_counts
            .entry(tr.chromosome.to_string())
            .or_default() += 1;
    }
    let payload = PayloadV1 {
        transcripts: transcripts.into_iter().map(Into::into).collect(),
        enrichment,
    };
    header.semantic_sha256 = semantic_digest(&payload)?;
    let dir = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let compressed = tempfile::NamedTempFile::new_in(dir)?;
    let mut encoder = zstd::Encoder::new(std::io::BufWriter::new(compressed.reopen()?), 1)?;
    encoder.include_checksum(true)?;
    bincode::serialize_into(&mut encoder, &payload)?;
    encoder.finish()?.flush()?;
    header.payload_bytes = compressed.as_file().metadata()?.len();
    header.payload_sha256 = sha256(compressed.path())?;
    let bytes = serde_json::to_vec(&header)?;
    ensure!(bytes.len() <= MAX_HEADER, "Cache header exceeds size limit");
    let mut pending = tempfile::NamedTempFile::new_in(dir)?;
    pending.write_all(MAGIC)?;
    pending.write_all(&VERSION.to_le_bytes())?;
    pending.write_all(&(bytes.len() as u32).to_le_bytes())?;
    pending.write_all(&bytes)?;
    io::copy(&mut File::open(compressed.path())?, &mut pending)?;
    pending.flush()?;
    pending.as_file().sync_all()?;
    pending
        .persist(path)
        .map_err(|e| anyhow::anyhow!("Publishing cache: {}", e.error))?;
    Ok(header)
}
