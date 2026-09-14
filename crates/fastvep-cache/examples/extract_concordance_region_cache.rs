//! Memory-bounded development fixture extraction; never release qualification.
use anyhow::{ensure, Result};
use fastvep_cache::annocat_cache::{self, Enrichment};
use sha2::{Digest, Sha256};
use std::{collections::{BTreeMap, BTreeSet}, fs::File, io::{self, BufRead, Read, Seek, SeekFrom}, path::Path};

// Reuse the frozen wire layout rather than duplicating it in a fixture tool.
#[path = "../src/transcript_wire.rs"]
mod transcript_wire;

struct Hashed<R>(R, Sha256);

fn same_contig(left: &str, right: &str) -> bool {
    left.trim_start_matches("chr") == right.trim_start_matches("chr")
        || (fastvep_genome::is_mitochondrial(left) && fastvep_genome::is_mitochondrial(right))
}

#[test]
fn region_selection_accepts_mitochondrial_aliases() {
    assert!(same_contig("MT", "chrM"));
    assert!(same_contig("chr1", "1"));
    assert!(!same_contig("MT", "1"));
}
impl<R: Read> Read for Hashed<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.0.read(buf)?;
        self.1.update(&buf[..n]);
        Ok(n)
    }
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    ensure!(args.len() >= 4, "usage: extract_concordance_region_cache CACHE OUTPUT INPUT_VCF...");
    let source = Path::new(&args[1]);
    let output = Path::new(&args[2]);
    ensure!(!output.exists(), "Output already exists");
    let mut regions = BTreeSet::new();
    for input in &args[3..] {
        for line in io::BufReader::new(File::open(input)?).lines() {
            let line = line?;
            if line.starts_with('#') { continue; }
            let cols: Vec<_> = line.split('\t').collect();
            ensure!(cols.len() >= 5, "Invalid VCF row");
            let start: u64 = cols[1].parse()?;
            regions.insert((cols[0].trim_start_matches("chr").to_string(),
                start.saturating_sub(5001), start + cols[3].len() as u64 + 5001));
        }
    }
    let mut header = annocat_cache::read_header(source)?.ok_or_else(|| anyhow::anyhow!("ANNOCATC1 required"))?;
    let mut file = File::open(source)?;
    file.seek(SeekFrom::Start(12))?;
    let mut word = [0; 4]; file.read_exact(&mut word)?;
    let payload_start = 16 + u32::from_le_bytes(word) as u64;
    file.seek(SeekFrom::Start(payload_start))?;
    let mut compressed = Hashed(&mut file, Sha256::new());
    io::copy(&mut compressed, &mut io::sink())?;
    ensure!(format!("{:x}", compressed.1.finalize()) == header.payload_sha256, "Payload checksum mismatch");
    file.seek(SeekFrom::Start(payload_start))?;
    let mut decoder = Hashed(zstd::Decoder::new(file)?, Sha256::new());
    let count: u64 = bincode::deserialize_from(&mut decoder)?;
    ensure!(count == header.transcript_count, "Transcript count mismatch");
    let mut kept = Vec::new();
    let mut coding = 0;
    let mut contigs = BTreeMap::new();
    let mut keys = BTreeSet::new();
    for _ in 0..count {
        let tr: transcript_wire::Transcript = bincode::deserialize_from(&mut decoder)?;
        let tr: fastvep_genome::Transcript = tr.into();
        coding += u64::from(tr.is_coding());
        *contigs.entry(tr.chromosome.to_string()).or_insert(0u64) += 1;
        if regions.iter().any(|(chrom, start, end)|
            same_contig(&tr.chromosome, chrom) && tr.start <= *end && tr.end >= *start) {
            keys.insert(format!("{}:{}", tr.chromosome, tr.stable_id));
            kept.push(tr);
        }
    }
    let mut enrichment: BTreeMap<String, Enrichment> = bincode::deserialize_from(&mut decoder)?;
    let mut trailing = [0];
    ensure!(decoder.read(&mut trailing)? == 0, "Trailing payload data");
    ensure!(format!("{:x}", decoder.1.finalize()) == header.semantic_sha256, "Semantic checksum mismatch");
    ensure!(coding == header.coding_transcript_count && contigs == header.contig_counts, "Cache counts mismatch");
    enrichment.retain(|key, _| keys.contains(key));
    header.provenance = serde_json::json!({"purpose": "Development region subset; not release qualification",
        "parentSha256": annocat_cache::sha256(source)?, "parentProvenance": header.provenance,
        "regions": regions, "inputVcfs": args[3..]});
    let result = annocat_cache::save(kept, enrichment, header, output)?;
    println!("{}", serde_json::to_string_pretty(&result)?);
    Ok(())
}
