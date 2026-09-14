use anyhow::{bail, Context, Result};
use fastvep_cache::transcript_cache::{load_cache, save_cache, verify_cache};
use fastvep_core::Strand;
use fastvep_genome::{Exon, Gene, Transcript};
use flate2::read::GzDecoder;
use serde::Deserialize;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::sync::Arc;

const PRIMARY_CONTIGS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT",
];

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct InventoryTranscript {
    transcript_id: String,
    version: Option<u32>,
    chromosome: String,
    start: u64,
    end: u64,
    strand: i8,
    biotype: String,
    source: Option<String>,
    gene: InventoryGene,
    exons: Vec<InventoryExon>,
    translation: Option<Value>,
    coding: InventoryCoding,
    metadata: InventoryMetadata,
    identifiers: InventoryIdentifiers,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct InventoryGene {
    stable_id: String,
    symbol: Option<String>,
    symbol_source: Option<String>,
    hgnc_id: Option<String>,
    start: Option<u64>,
    end: Option<u64>,
    strand: Option<i8>,
}

#[derive(Deserialize)]
struct InventoryExon {
    start: u64,
    end: u64,
    strand: i8,
    rank: u32,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct InventoryCoding {
    cdna_start: Option<u64>,
    cdna_end: Option<u64>,
    genomic_start: Option<u64>,
    genomic_end: Option<u64>,
    start_phase: Option<u64>,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct InventoryMetadata {
    canonical: bool,
    gencode_primary: bool,
    flags: Vec<String>,
    mane_select: Vec<String>,
    mane_plus_clinical: Vec<String>,
    tsl: Vec<String>,
    appris: Vec<String>,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct InventoryIdentifiers {
    ccds: Vec<String>,
    protein: Vec<String>,
    refseq: Vec<String>,
    swissprot: Vec<String>,
    trembl: Vec<String>,
    uniparc: Vec<String>,
}

fn sha256(path: &Path) -> Result<String> {
    let mut file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut digest = Sha256::new();
    let mut buffer = vec![0u8; 1024 * 1024];
    loop {
        let read = file.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        digest.update(&buffer[..read]);
    }
    Ok(format!("{:x}", digest.finalize()))
}

fn reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    if path
        .extension()
        .and_then(|value| value.to_str())
        .is_some_and(|value| value.eq_ignore_ascii_case("gz"))
    {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn normalized_chromosome(value: &str) -> &str {
    match value.strip_prefix("chr").unwrap_or(value) {
        "M" => "MT",
        normalized => normalized,
    }
}

fn primary(value: &str) -> bool {
    PRIMARY_CONTIGS.contains(&normalized_chromosome(value))
}

fn cache_chromosome(value: &str, prefixed: bool) -> Arc<str> {
    let normalized = normalized_chromosome(value);
    if !prefixed {
        return Arc::from(normalized);
    }
    Arc::from(if normalized == "MT" {
        "chrM".to_owned()
    } else {
        format!("chr{normalized}")
    })
}

fn optional(values: &[String], label: &str) -> Result<Option<String>> {
    if values.len() > 1 {
        bail!("{label} has {} values but fastVEP stores one", values.len());
    }
    Ok(values.first().cloned())
}

fn tsl(values: &[String]) -> Result<Option<u8>> {
    let value = optional(values, "TSL")?;
    match value {
        None => Ok(None),
        Some(value) if value.starts_with("tslNA") => Ok(None),
        Some(value) => Ok(Some(value.parse().context("parsing TSL")?)),
    }
}

fn added_transcript(value: Value, prefixed: bool) -> Result<Transcript> {
    let record: InventoryTranscript =
        serde_json::from_value(value).context("parsing official-only transcript")?;
    if record.translation.is_some()
        || record.coding.cdna_start.is_some()
        || record.coding.cdna_end.is_some()
        || record.coding.genomic_start.is_some()
        || record.coding.genomic_end.is_some()
    {
        bail!(
            "official-only transcript {} is coding; sequence-complete addition is required",
            record.transcript_id
        );
    }
    let chromosome = cache_chromosome(&record.chromosome, prefixed);
    let strand = Strand::from_int(record.strand);
    let biotype: Arc<str> = Arc::from(record.biotype);
    let exons = record
        .exons
        .into_iter()
        .map(|exon| Exon {
            stable_id: format!("{}:exon:{}", record.transcript_id, exon.rank),
            start: exon.start,
            end: exon.end,
            strand: Strand::from_int(exon.strand),
            phase: -1,
            end_phase: -1,
            rank: exon.rank,
        })
        .collect();
    Ok(Transcript {
        stable_id: Arc::from(record.transcript_id),
        version: record.version,
        gene: Gene {
            stable_id: Arc::from(record.gene.stable_id),
            symbol: record.gene.symbol.map(Arc::from),
            symbol_source: record.gene.symbol_source,
            hgnc_id: record.gene.hgnc_id,
            biotype: biotype.clone(),
            chromosome: chromosome.clone(),
            start: record.gene.start.unwrap_or(record.start),
            end: record.gene.end.unwrap_or(record.end),
            strand: Strand::from_int(record.gene.strand.unwrap_or(record.strand)),
        },
        biotype,
        chromosome,
        start: record.start,
        end: record.end,
        strand,
        exons,
        translation: None,
        cdna_coding_start: None,
        cdna_coding_end: None,
        coding_region_start: None,
        coding_region_end: None,
        spliced_seq: None,
        translateable_seq: None,
        peptide: None,
        canonical: record.metadata.canonical,
        mane_select: optional(&record.metadata.mane_select, "MANE Select")?,
        mane_plus_clinical: optional(&record.metadata.mane_plus_clinical, "MANE Plus Clinical")?,
        tsl: tsl(&record.metadata.tsl)?,
        appris: optional(&record.metadata.appris, "APPRIS")?,
        ccds: optional(&record.identifiers.ccds, "CCDS")?,
        protein_id: optional(&record.identifiers.protein, "protein identifier")?,
        protein_version: None,
        swissprot: record.identifiers.swissprot,
        trembl: record.identifiers.trembl,
        uniparc: record.identifiers.uniparc,
        refseq_id: optional(&record.identifiers.refseq, "RefSeq")?,
        source: record.source,
        gencode_primary: record.metadata.gencode_primary,
        flags: record.metadata.flags,
        codon_table_start_phase: record.coding.start_phase.unwrap_or(0),
        reference_peptide: None,
    })
}

fn strings(value: &Value, label: &str) -> Result<Vec<String>> {
    value
        .as_array()
        .with_context(|| format!("{label} is not an array"))?
        .iter()
        .map(|item| {
            item.as_str()
                .map(str::to_owned)
                .with_context(|| format!("{label} contains a non-string"))
        })
        .collect()
}

fn one_string(value: &Value, label: &str) -> Result<Option<String>> {
    optional(&strings(value, label)?, label)
}

fn current_value(transcript: &Transcript, path: &str) -> Result<Value> {
    let optional_string = |value: Option<&str>| match value {
        Some(value) => Value::String(value.to_owned()),
        None => Value::Null,
    };
    Ok(match path {
        "gene.symbol" => optional_string(transcript.gene.symbol.as_deref()),
        "gene.symbolSource" => optional_string(transcript.gene.symbol_source.as_deref()),
        "gene.hgncId" => optional_string(transcript.gene.hgnc_id.as_deref()),
        "identifiers.refseq" => json!(transcript.refseq_id.iter().collect::<Vec<_>>()),
        "identifiers.swissprot" => json!(transcript.swissprot),
        "identifiers.trembl" => json!(transcript.trembl),
        "identifiers.uniparc" => json!(transcript.uniparc),
        "metadata.appris" => json!(transcript.appris.iter().collect::<Vec<_>>()),
        "metadata.flags" => {
            let mut values = transcript.flags.clone();
            values.sort();
            values.dedup();
            json!(values)
        }
        "metadata.tsl" => json!(transcript
            .tsl
            .map(|value| vec![value.to_string()])
            .unwrap_or_default()),
        "translation.version" => json!(transcript.protein_version),
        _ => bail!("unsupported direct path {path}"),
    })
}

fn set_value(transcript: &mut Transcript, path: &str, value: &Value) -> Result<()> {
    let nullable = |label: &str| -> Result<Option<String>> {
        if value.is_null() {
            Ok(None)
        } else {
            Ok(Some(
                value
                    .as_str()
                    .with_context(|| format!("{label} is not a string"))?
                    .to_owned(),
            ))
        }
    };
    match path {
        "gene.symbol" => transcript.gene.symbol = nullable(path)?.map(Arc::from),
        "gene.symbolSource" => transcript.gene.symbol_source = nullable(path)?,
        "gene.hgncId" => transcript.gene.hgnc_id = nullable(path)?,
        "identifiers.refseq" => transcript.refseq_id = one_string(value, path)?,
        "identifiers.swissprot" => transcript.swissprot = strings(value, path)?,
        "identifiers.trembl" => transcript.trembl = strings(value, path)?,
        "identifiers.uniparc" => transcript.uniparc = strings(value, path)?,
        "metadata.appris" => transcript.appris = one_string(value, path)?,
        "metadata.flags" => transcript.flags = strings(value, path)?,
        "metadata.tsl" => {
            transcript.tsl = one_string(value, path)?
                .map(|value| value.parse().context("parsing TSL"))
                .transpose()?
        }
        "translation.version" => {
            transcript.protein_version = if value.is_null() {
                None
            } else {
                Some(
                    value
                        .as_u64()
                        .context("translation.version is not an unsigned integer")?
                        .try_into()
                        .context("translation.version exceeds u32")?,
                )
            }
        }
        _ => bail!("unsupported direct path {path}"),
    }
    Ok(())
}

fn peptide_fingerprint(transcript: &Transcript) -> (Value, Value) {
    let Some(peptide) = transcript
        .peptide
        .as_deref()
        .filter(|value| !value.is_empty())
    else {
        return (Value::Null, Value::Null);
    };
    let peptide = peptide.strip_suffix('*').unwrap_or(peptide);
    (
        json!(peptide.len()),
        json!(format!("{:x}", Sha256::digest(peptide.as_bytes()))),
    )
}

fn apply_sequence_edits(transcript: &mut Transcript, value: &Value) -> Result<()> {
    let edits = value.as_array().context("sequenceEdits is not an array")?;
    let peptide = transcript
        .peptide
        .as_mut()
        .with_context(|| format!("{} has sequence edits but no peptide", transcript.stable_id))?;
    let mut bytes = peptide.as_bytes().to_vec();
    for edit in edits {
        let start = edit["start"].as_u64().context("sequence edit start")?;
        let end = edit["end"].as_u64().context("sequence edit end")?;
        let alternate = edit["alternate"]
            .as_str()
            .context("sequence edit alternate")?
            .as_bytes();
        if start == 0 || start != end || alternate.len() != 1 {
            bail!("only one-residue sequence edits are supported in the test candidate");
        }
        let target = bytes
            .get_mut(start as usize - 1)
            .with_context(|| format!("{} sequence edit is out of range", transcript.stable_id))?;
        *target = alternate[0];
    }
    *peptide = String::from_utf8(bytes).context("edited peptide is not UTF-8")?;
    Ok(())
}

fn count(map: &mut BTreeMap<String, u64>, key: &str) {
    *map.entry(key.to_owned()).or_default() += 1;
}

fn remember(map: &mut BTreeMap<String, Vec<String>>, key: &str, transcript: &str) {
    let examples = map.entry(key.to_owned()).or_default();
    if examples.len() < 20 {
        examples.push(transcript.to_owned());
    }
}

fn parse_args() -> Result<(PathBuf, PathBuf, PathBuf, String, bool)> {
    let mut args = env::args_os().skip(1);
    let cache = args.next().map(PathBuf::from).context("missing CACHE")?;
    let output = args.next().map(PathBuf::from).context("missing OUTPUT")?;
    let delta = args.next().map(PathBuf::from).context("missing DELTA")?;
    let expected_sha = args
        .next()
        .context("missing EXPECTED_CACHE_SHA256")?
        .to_string_lossy()
        .into_owned();
    let allow_unrepresented = args
        .next()
        .is_some_and(|value| value == "--allow-unrepresented");
    if args.next().is_some() || expected_sha.len() != 64 || cache == output || output.exists() {
        bail!(
            "usage: apply_inventory_candidate CACHE OUTPUT DELTA EXPECTED_CACHE_SHA256 [--allow-unrepresented]"
        );
    }
    Ok((cache, output, delta, expected_sha, allow_unrepresented))
}

fn main() -> Result<()> {
    let (cache, output, delta, expected_sha, allow_unrepresented) = parse_args()?;
    let cache_sha = sha256(&cache)?;
    if !cache_sha.eq_ignore_ascii_case(&expected_sha) {
        bail!("input cache SHA-256 mismatch: expected {expected_sha}, found {cache_sha}");
    }
    let verification = verify_cache(&cache, false)?;
    if verification.cache_format != "FSTVEP02" {
        bail!("inventory candidate requires an FSTVEP02 cache");
    }
    let mut transcripts = load_cache(&cache)?;
    let prefixed = transcripts
        .iter()
        .any(|transcript| transcript.chromosome.as_ref() == "chr1");
    let mut index = HashMap::new();
    for (position, transcript) in transcripts.iter().enumerate() {
        if primary(&transcript.chromosome)
            && index
                .insert(transcript.stable_id.to_string(), position)
                .is_some()
        {
            bail!(
                "duplicate primary-contig transcript {}",
                transcript.stable_id
            );
        }
    }

    let direct_paths = BTreeSet::from([
        "gene.symbol",
        "gene.symbolSource",
        "gene.hgncId",
        "identifiers.refseq",
        "identifiers.swissprot",
        "identifiers.trembl",
        "identifiers.uniparc",
        "metadata.appris",
        "metadata.flags",
        "translation.version",
    ]);
    let mut applied = BTreeMap::new();
    let mut ignored = BTreeMap::new();
    let mut unrepresented = BTreeMap::new();
    let mut unrepresented_examples = BTreeMap::new();
    let mut touched = BTreeSet::new();
    let mut remove = BTreeSet::new();
    let mut additions = Vec::new();

    for (number, line) in reader(&delta)?.lines().enumerate() {
        let line =
            line.with_context(|| format!("reading {} line {}", delta.display(), number + 1))?;
        if line.trim().is_empty() {
            continue;
        }
        let row: Value = serde_json::from_str(&line)
            .with_context(|| format!("parsing {} line {}", delta.display(), number + 1))?;
        let id = row["transcriptId"]
            .as_str()
            .context("delta row has no transcriptId")?;
        let changes = row["changes"]
            .as_array()
            .context("delta row has no changes")?;
        let sequence_edits = changes
            .iter()
            .find(|change| change["path"] == "sequenceEdits");
        let position = index.get(id).copied();
        let peptide_before = position.and_then(|found| transcripts[found].peptide.clone());
        let peptide_fingerprint_before =
            position.map(|found| peptide_fingerprint(&transcripts[found]));
        let mut sequence_edit_failed = false;

        for change in changes {
            let path = change["path"].as_str().context("change has no path")?;
            let impact = change["impact"].as_str().context("change has no impact")?;
            if impact != "output-active" {
                count(&mut ignored, path);
                continue;
            }
            if path == "membership.candidateOnly" {
                let found = position.with_context(|| format!("candidate-only {id} is absent"))?;
                if normalized_chromosome(&transcripts[found].chromosome)
                    != normalized_chromosome(
                        change["candidate"]["chromosome"].as_str().unwrap_or(""),
                    )
                {
                    bail!("candidate-only {id} chromosome mismatch");
                }
                remove.insert(id.to_owned());
                touched.insert(id.to_owned());
                count(&mut applied, path);
                continue;
            }
            if path == "membership.officialOnly" {
                if position.is_some() {
                    bail!("official-only {id} already exists");
                }
                additions.push(added_transcript(change["official"].clone(), prefixed)?);
                touched.insert(id.to_owned());
                count(&mut applied, path);
                continue;
            }
            if path == "identifiers.uniprotIsoform" {
                count(&mut ignored, path);
                continue;
            }
            if path == "metadata.tsl" {
                let found =
                    position.with_context(|| format!("shared transcript {id} is absent"))?;
                let current = current_value(&transcripts[found], path)?;
                if current != change["candidate"] {
                    bail!("{id} {path} cache value does not match the delta candidate");
                }
                let official = strings(&change["official"], path)?;
                let resolved = tsl(&official)?;
                if resolved.is_none()
                    && official
                        .first()
                        .is_some_and(|value| value.starts_with("tslNA"))
                {
                    count(&mut ignored, path);
                } else {
                    transcripts[found].tsl = resolved;
                    touched.insert(id.to_owned());
                    count(&mut applied, path);
                }
                continue;
            }
            if direct_paths.contains(path) {
                let found =
                    position.with_context(|| format!("shared transcript {id} is absent"))?;
                let current = current_value(&transcripts[found], path)?;
                if current != change["candidate"] {
                    bail!(
                        "{id} {path} cache value {} does not match delta candidate {}",
                        current,
                        change["candidate"]
                    );
                }
                set_value(&mut transcripts[found], path, &change["official"])?;
                touched.insert(id.to_owned());
                count(&mut applied, path);
                continue;
            }
            if path == "sequenceEdits" {
                let found = position
                    .with_context(|| format!("sequence-edited transcript {id} is absent"))?;
                if apply_sequence_edits(&mut transcripts[found], &change["official"]).is_err() {
                    sequence_edit_failed = true;
                }
                continue;
            }
            if path == "sequences.peptide.length" || path == "sequences.peptide.sha256" {
                position.with_context(|| format!("peptide transcript {id} is absent"))?;
                let (candidate_length, candidate_sha) = peptide_fingerprint_before
                    .as_ref()
                    .context("peptide fingerprint is absent")?;
                let current = if path.ends_with("length") {
                    candidate_length.clone()
                } else {
                    candidate_sha.clone()
                };
                if current != change["candidate"] {
                    bail!("{id} {path} does not match the delta candidate");
                }
                if sequence_edits.is_none() {
                    count(&mut unrepresented, path);
                    remember(&mut unrepresented_examples, path, id);
                }
                continue;
            }
            count(&mut unrepresented, path);
            remember(&mut unrepresented_examples, path, id);
        }

        if sequence_edits.is_some() {
            let found = position.context("sequence-edited transcript is absent")?;
            let (length, digest) = peptide_fingerprint(&transcripts[found]);
            let has_peptide_check = changes.iter().any(|change| {
                matches!(
                    change["path"].as_str(),
                    Some("sequences.peptide.length" | "sequences.peptide.sha256")
                )
            });
            let mut verified = !sequence_edit_failed && has_peptide_check;
            for change in changes {
                let path = change["path"].as_str().unwrap_or("");
                let actual = if path == "sequences.peptide.length" {
                    Some(&length)
                } else if path == "sequences.peptide.sha256" {
                    Some(&digest)
                } else {
                    None
                };
                if let Some(actual) = actual {
                    if actual != &change["official"] {
                        verified = false;
                    }
                }
            }
            if verified {
                touched.insert(id.to_owned());
                for change in changes {
                    let path = change["path"].as_str().unwrap_or("");
                    if path == "sequenceEdits"
                        || path == "sequences.peptide.length"
                        || path == "sequences.peptide.sha256"
                    {
                        count(&mut applied, path);
                    }
                }
            } else {
                transcripts[found].peptide = peptide_before;
                for change in changes {
                    let path = change["path"].as_str().unwrap_or("");
                    if path == "sequenceEdits"
                        || path == "sequences.peptide.length"
                        || path == "sequences.peptide.sha256"
                    {
                        count(&mut unrepresented, path);
                        remember(&mut unrepresented_examples, path, id);
                    }
                }
            }
        }
    }

    let unrepresented_count: u64 = unrepresented.values().sum();
    if unrepresented_count != 0 && !allow_unrepresented {
        bail!(
            "candidate has {unrepresented_count} unrepresented output-active changes: {}",
            serde_json::to_string(&unrepresented)?
        );
    }

    transcripts.retain(|transcript| {
        !primary(&transcript.chromosome) || !remove.contains(transcript.stable_id.as_ref())
    });
    let additions_count = additions.len();
    transcripts.extend(additions);
    save_cache(&transcripts, &output)?;
    let output_verification = verify_cache(&output, false)?;
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "schemaVersion": 1,
            "testOnly": true,
            "complete": unrepresented_count == 0,
            "inputCache": cache,
            "inputCacheSha256": cache_sha,
            "delta": delta,
            "deltaSha256": sha256(&delta)?,
            "appliedChangesByPath": applied,
            "ignoredChangesByPath": ignored,
            "unrepresentedChangesByPath": unrepresented,
            "unrepresentedTranscriptExamples": unrepresented_examples,
            "transcriptsTouched": touched.len(),
            "transcriptsRemoved": remove.len(),
            "transcriptsAdded": additions_count,
            "outputCache": output,
            "outputCacheSha256": sha256(&output)?,
            "outputCacheBytes": output_verification.cache_bytes,
            "transcriptCount": output_verification.transcript_count,
        }))?
    );
    Ok(())
}
