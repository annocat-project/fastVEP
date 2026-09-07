use anyhow::{bail, Context, Result};
use fastvep_cache::transcript_cache::{load_cache, verify_cache};
use fastvep_genome::Transcript;
use flate2::{write::GzEncoder, Compression, GzBuilder};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::collections::BTreeSet;
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

fn fingerprint(sequence: Option<&str>, trim_terminal_stop: bool) -> Value {
    match sequence.filter(|sequence| !sequence.is_empty()) {
        Some(sequence) => {
            let sequence = if trim_terminal_stop {
                sequence.strip_suffix('*').unwrap_or(sequence)
            } else {
                sequence
            };
            json!({
                "length": sequence.len(),
                "sha256": format!("{:x}", Sha256::digest(sequence.as_bytes())),
            })
        }
        None => Value::Null,
    }
}

fn strings(values: impl IntoIterator<Item = String>) -> Vec<String> {
    let mut values: Vec<_> = values
        .into_iter()
        .filter(|value| !value.is_empty() && value != "-")
        .collect();
    values.sort();
    values.dedup();
    values
}

fn optional(value: &Option<String>) -> Vec<String> {
    value.iter().cloned().collect()
}

fn normalized_chromosome(chromosome: &str) -> &str {
    match chromosome.strip_prefix("chr").unwrap_or(chromosome) {
        "M" => "MT",
        normalized => normalized,
    }
}

fn record(transcript: &Transcript) -> Value {
    let mut flags = transcript.flags.clone();
    flags.sort();
    flags.dedup();
    json!({
        "recordType": "transcript",
        "transcriptId": transcript.stable_id,
        "version": transcript.version,
        "chromosome": normalized_chromosome(&transcript.chromosome),
        "start": transcript.start,
        "end": transcript.end,
        "strand": transcript.strand.as_int(),
        "biotype": transcript.biotype,
        "source": transcript.source,
        "gene": {
            "stableId": transcript.gene.stable_id,
            "version": Value::Null,
            "symbol": transcript.gene.symbol,
            "symbolSource": transcript.gene.symbol_source,
            "hgncId": transcript.gene.hgnc_id,
            "start": transcript.gene.start,
            "end": transcript.gene.end,
            "strand": transcript.gene.strand.as_int(),
        },
        "exons": transcript.exons.iter().map(|exon| json!({
            "start": exon.start,
            "end": exon.end,
            "strand": exon.strand.as_int(),
            "rank": exon.rank,
        })).collect::<Vec<_>>(),
        "translation": transcript.translation.as_ref().map(|translation| json!({
            "stableId": translation.stable_id,
            "version": transcript.protein_version,
            "genomicStart": translation.genomic_start,
            "genomicEnd": translation.genomic_end,
            "startExonRank": translation.start_exon_rank,
            "startExonOffset": translation.start_exon_offset,
            "endExonRank": translation.end_exon_rank,
            "endExonOffset": translation.end_exon_offset,
        })),
        "coding": {
            "cdnaStart": transcript.cdna_coding_start,
            "cdnaEnd": transcript.cdna_coding_end,
            "genomicStart": transcript.coding_region_start,
            "genomicEnd": transcript.coding_region_end,
            "startPhase": transcript.translation.as_ref().map(|_| transcript.codon_table_start_phase),
            "codonTable": if fastvep_genome::is_mitochondrial(&transcript.chromosome) { 2 } else { 1 },
        },
        "metadata": {
            "canonical": transcript.canonical,
            "gencodeBasic": Value::Null,
            "gencodePrimary": transcript.gencode_primary,
            "flags": flags,
            "maneSelect": optional(&transcript.mane_select),
            "manePlusClinical": optional(&transcript.mane_plus_clinical),
            "tsl": transcript.tsl.map(|value| vec![value.to_string()]).unwrap_or_default(),
            "appris": optional(&transcript.appris),
            "matureMirnaRanges": Value::Null,
            "ncrnaStructures": Value::Null,
            "rnaEdits": Value::Null,
            "refseqComparisonAttributes": Value::Null,
            "otherAttributes": Value::Null,
        },
        "identifiers": {
            "ccds": optional(&transcript.ccds),
            "refseq": optional(&transcript.refseq_id),
            "protein": optional(&transcript.protein_id),
            "swissprot": strings(transcript.swissprot.clone()),
            "trembl": strings(transcript.trembl.clone()),
            "uniparc": strings(transcript.uniparc.clone()),
            "uniprotIsoform": Value::Null,
        },
        "sequences": {
            "translateable": fingerprint(transcript.translateable_seq.as_deref(), false),
            "peptide": fingerprint(transcript.peptide.as_deref(), true),
        },
        "sequenceEdits": Value::Null,
    })
}

fn writer(path: &Path) -> Result<Box<dyn Write>> {
    let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    if path
        .extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("gz"))
    {
        let encoder: GzEncoder<File> = GzBuilder::new().mtime(0).write(file, Compression::new(6));
        Ok(Box::new(BufWriter::new(encoder)))
    } else {
        Ok(Box::new(BufWriter::new(file)))
    }
}

fn parse_args() -> Result<(PathBuf, PathBuf, BTreeSet<String>, u32, String)> {
    let mut args = env::args_os().skip(1);
    let input = args
        .next()
        .map(PathBuf::from)
        .context("missing CACHE argument")?;
    let output = args
        .next()
        .map(PathBuf::from)
        .context("missing OUTPUT argument")?;
    let contigs = args
        .next()
        .context("missing comma-separated CONTIGS argument")?
        .to_string_lossy()
        .split(',')
        .filter(|value| !value.is_empty())
        .map(str::to_owned)
        .collect::<BTreeSet<_>>();
    let release = args
        .next()
        .context("missing RELEASE argument")?
        .to_string_lossy()
        .parse::<u32>()
        .context("RELEASE must be an integer")?;
    let assembly = args
        .next()
        .context("missing ASSEMBLY argument")?
        .to_string_lossy()
        .into_owned();
    if contigs.is_empty() || release == 0 || assembly.is_empty() || args.next().is_some() {
        bail!("usage: export_transcript_inventory CACHE OUTPUT CONTIGS RELEASE ASSEMBLY")
    }
    Ok((input, output, contigs, release, assembly))
}

fn main() -> Result<()> {
    let (input, output, contigs, release, assembly) = parse_args()?;
    let input_bytes = input.metadata()?.len();
    let format = verify_cache(&input, false)?.cache_format;
    let mut transcripts = load_cache(&input)?;
    transcripts.retain(|transcript| contigs.contains(transcript.chromosome.as_ref()));
    transcripts.sort_by(|left, right| {
        left.stable_id
            .cmp(&right.stable_id)
            .then(left.version.cmp(&right.version))
    });
    for pair in transcripts.windows(2) {
        if pair[0].stable_id == pair[1].stable_id {
            bail!("duplicate transcript {} in cache", pair[0].stable_id)
        }
    }

    let mut output = writer(&output)?;
    writeln!(
        output,
        "{}",
        json!({
            "recordType": "manifest",
            "schemaVersion": 1,
            "source": "fastvep-transcript-cache",
            "cacheFormat": format,
            "cacheBytes": input_bytes,
            "release": release,
            "assembly": assembly,
            "contigs": contigs.iter().map(|contig| normalized_chromosome(contig)).collect::<BTreeSet<_>>(),
            "capabilities": [
                "transcript-membership", "transcript-structure", "coding-boundaries",
                "translation-boundaries", "identifiers", "cds-completeness-flags",
                "translateable-sequence", "peptide", "codon-table"
            ],
        })
    )?;
    let mut digest = Sha256::new();
    for transcript in &transcripts {
        let line = serde_json::to_string(&record(transcript))?;
        digest.update(line.as_bytes());
        digest.update(b"\n");
        writeln!(output, "{line}")?;
    }
    writeln!(
        output,
        "{}",
        json!({
            "recordType": "summary",
            "transcriptCount": transcripts.len(),
            "contentSha256": format!("{:x}", digest.finalize()),
        })
    )?;
    output.flush()?;
    eprintln!("exported {} transcripts", transcripts.len());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{fingerprint, normalized_chromosome};

    #[test]
    fn sequence_fingerprint_is_sha256() {
        assert_eq!(
            fingerprint(Some("M"), false)["sha256"],
            "08f271887ce94707da822d5263bae19d5519cb3614e0daedc4c7ce5dab7473f1"
        );
    }

    #[test]
    fn terminal_stop_is_not_part_of_the_peptide_fingerprint() {
        assert_eq!(fingerprint(Some("M*"), true), fingerprint(Some("M"), false));
        assert!(fingerprint(Some(""), false).is_null());
    }

    #[test]
    fn chromosome_aliases_are_normalized() {
        assert_eq!(normalized_chromosome("chr1"), "1");
        assert_eq!(normalized_chromosome("chrM"), "MT");
        assert_eq!(normalized_chromosome("MT"), "MT");
    }
}
