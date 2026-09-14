//! Public Ensembl Core enrichment for the pinned GRCh38 release-115 GFF3.
use crate::annocat_cache::{self, Enrichment, Header, SequenceEdit};
use crate::providers::SequenceProvider;
use anyhow::{ensure, Context, Result};
use flate2::read::MultiGzDecoder;
use serde::Deserialize;
use std::{
    collections::{BTreeMap, HashMap},
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
    sync::Arc,
};

#[derive(Deserialize)]
struct Source {
    filename: String,
    table: String,
    url: String,
    bytes: u64,
    sha256: String,
    columns: Vec<String>,
}
#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct Manifest {
    schema_version: u32,
    species: String,
    assembly: String,
    ensembl_release: u32,
    vep_release: String,
    sources: Vec<Source>,
    gff3: Input,
    reference: Input,
}
#[derive(Deserialize)]
struct Input {
    bytes: u64,
    sha256: String,
}

fn verify_input(path: &Path, bytes: u64, digest: &str) -> Result<()> {
    ensure!(
        path.metadata()?.len() == bytes && annocat_cache::sha256(path)? == digest,
        "Source size/checksum mismatch: {}",
        path.display()
    );
    Ok(())
}

fn rows(
    manifest: &Manifest,
    directory: &Path,
    table: &str,
    mut visit: impl FnMut(&[&str], &HashMap<&str, usize>) -> Result<()>,
) -> Result<()> {
    let source = manifest
        .sources
        .iter()
        .find(|s| s.table == table)
        .with_context(|| format!("Missing {table} source"))?;
    let columns: HashMap<_, _> = source
        .columns
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    ensure!(!columns.is_empty(), "Missing release schema for {table}");
    let reader = BufReader::new(MultiGzDecoder::new(File::open(
        directory.join(&source.filename),
    )?));
    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let row: Vec<_> = line.split('\t').collect();
        ensure!(
            row.len() == columns.len(),
            "Wrong column count in {table}:{}",
            line_number + 1
        );
        visit(&row, &columns).with_context(|| format!("{table}:{}", line_number + 1))?;
    }
    Ok(())
}
fn field<'a>(row: &[&'a str], columns: &HashMap<&str, usize>, key: &str) -> Result<&'a str> {
    Ok(row[*columns
        .get(key)
        .with_context(|| format!("Missing source column {key}"))?])
}
fn optional(value: &str) -> Option<String> {
    (value != "\\N" && !value.is_empty()).then(|| value.to_owned())
}

#[derive(Default)]
struct CoreTranscript {
    id: String,
    gene: String,
    version: u32,
    translation: Option<String>,
    attributes: BTreeMap<String, Vec<String>>,
    flags: Vec<String>,
}
struct CoreGene {
    stable_id: String,
    version: u32,
    display: String,
    canonical: String,
}
struct CoreTranslation {
    stable_id: String,
    version: u32,
    attributes: Vec<(String, String)>,
}

fn one(attributes: &BTreeMap<String, Vec<String>>, code: &str) -> Result<Option<String>> {
    let Some(values) = attributes.get(code) else {
        return Ok(None);
    };
    ensure!(
        values.len() == 1,
        "Expected one {code} attribute, got {}",
        values.len()
    );
    Ok(values.first().cloned())
}

fn parse_edit(kind: &str, text: &str) -> Result<SequenceEdit> {
    let mut values = text.splitn(3, ' ');
    let start = values.next().context("Missing edit start")?.parse()?;
    let end = values.next().context("Missing edit end")?.parse()?;
    let replacement = values.next().unwrap_or("").to_owned();
    ensure!(
        start > 0 && end + 1 >= start && replacement.is_ascii(),
        "Invalid {kind} edit"
    );
    Ok(SequenceEdit {
        kind: kind.into(),
        start,
        end,
        replacement,
    })
}

fn apply_edits(sequence: &mut String, edits: &[SequenceEdit]) -> Result<()> {
    let mut sorted: Vec<_> = edits.iter().collect();
    sorted.sort_by_key(|e| std::cmp::Reverse(e.start));
    let mut boundary = sequence.len() as u64 + 1;
    for edit in sorted {
        ensure!(
            edit.end < boundary && edit.end <= sequence.len() as u64,
            "Overlapping or out-of-range {} edit",
            edit.kind
        );
        sequence.replace_range(
            (edit.start - 1) as usize..edit.end as usize,
            &edit.replacement,
        );
        boundary = edit.start;
    }
    Ok(())
}

fn consequence_peptide(cds: &str, chromosome: &str, edits: &[SequenceEdit]) -> Result<String> {
    translated_peptide(cds, chromosome, edits, false)
}

pub(crate) fn translated_peptide(cds: &str, chromosome: &str, edits: &[SequenceEdit], reference: bool) -> Result<String> {
    // VEP 115 TranscriptVariationAllele::peptide translates codons without the
    // full Transcript::translate initiator normalization, then applies Core edits.
    let table = if fastvep_genome::is_mitochondrial(chromosome) {
        fastvep_genome::mitochondrial_codon_table()
    } else {
        fastvep_genome::CodonTable::standard()
    };
    let mut peptide = table.translate_seq(cds.as_bytes());
    if reference {
        // Transcript::translate normalizes before modify_translation applies
        // source edits; allele peptide() omits this normalization.
        table.normalize_reference_initiator(&mut peptide, cds.as_bytes());
    }
    let mut peptide = String::from_utf8(peptide)?;
    apply_edits(&mut peptide, edits)?;
    Ok(peptide)
}

pub fn build(
    gff3: &Path,
    fasta: &Path,
    directory: &Path,
    manifest_path: &Path,
    output: &Path,
) -> Result<Header> {
    let provenance: serde_json::Value = serde_json::from_reader(File::open(manifest_path)?)?;
    let manifest: Manifest = serde_json::from_value(provenance.clone())?;
    ensure!(
        manifest.schema_version == 1
            && manifest.species == "homo_sapiens"
            && manifest.assembly == "GRCh38"
            && manifest.ensembl_release == 115
            && manifest.vep_release == "115.2",
        "Unsupported Core source identity"
    );
    verify_input(gff3, manifest.gff3.bytes, &manifest.gff3.sha256)?;
    verify_input(fasta, manifest.reference.bytes, &manifest.reference.sha256)?;
    for source in &manifest.sources {
        ensure!(
            !source.filename.contains(['/', '\\'])
                && source.url.starts_with(
                    "https://ftp.ensembl.org/pub/release-115/mysql/homo_sapiens_core_115_38/"
                ),
            "Invalid Core source path"
        );
        verify_input(
            &directory.join(&source.filename),
            source.bytes,
            &source.sha256,
        )?;
    }
    eprintln!("Verified public Ensembl 115 inputs");
    let mut attribute_types = HashMap::new();
    rows(&manifest, directory, "attrib_type", |r, c| {
        attribute_types.insert(
            field(r, c, "attrib_type_id")?.to_owned(),
            field(r, c, "code")?.to_owned(),
        );
        Ok(())
    })?;
    let mut core = HashMap::new();
    let mut stable_by_id = HashMap::new();
    rows(&manifest, directory, "transcript", |r, c| {
        if field(r, c, "is_current")? != "1" {
            return Ok(());
        }
        let stable = field(r, c, "stable_id")?.to_owned();
        let id = field(r, c, "transcript_id")?.to_owned();
        stable_by_id.insert(id.clone(), stable.clone());
        ensure!(
            core.insert(
                stable,
                CoreTranscript {
                    id,
                    gene: field(r, c, "gene_id")?.into(),
                    version: field(r, c, "version")?.parse()?,
                    translation: optional(field(r, c, "canonical_translation_id")?),
                    attributes: BTreeMap::new(),
                    flags: vec![]
                }
            )
            .is_none(),
            "Duplicate current transcript stable ID"
        );
        Ok(())
    })?;
    rows(&manifest, directory, "transcript_attrib", |r, c| {
        if let Some(stable) = stable_by_id.get(field(r, c, "transcript_id")?) {
            let code = attribute_types
                .get(field(r, c, "attrib_type_id")?)
                .context("Unknown attribute type")?;
            // OutputFactory preserves get_all_Attributes order; do not sort CDS flags.
            if code.starts_with("cds_") {
                core.get_mut(stable).unwrap().flags.push(code.clone());
            }
            if matches!(
                code.as_str(),
                "readthrough_tra"
                    | "cds_start_NF"
                    | "cds_end_NF"
                    | "miRNA"
                    | "appris"
                    | "TSL"
                    | "MANE_Select"
                    | "MANE_Plus_Clinical"
                    | "gencode_primary"
                    | "_rna_edit"
                    | "_transl_start"
                    | "_transl_end"
            ) {
                core.get_mut(stable)
                    .unwrap()
                    .attributes
                    .entry(code.clone())
                    .or_default()
                    .push(field(r, c, "value")?.to_owned());
            }
        }
        Ok(())
    })?;
    let mut genes = HashMap::new();
    rows(&manifest, directory, "gene", |r, c| {
        genes.insert(
            field(r, c, "gene_id")?.to_owned(),
            CoreGene {
                stable_id: field(r, c, "stable_id")?.into(),
                version: field(r, c, "version")?.parse()?,
                display: field(r, c, "display_xref_id")?.into(),
                canonical: field(r, c, "canonical_transcript_id")?.into(),
            },
        );
        Ok(())
    })?;
    let mut databases = HashMap::new();
    rows(&manifest, directory, "external_db", |r, c| {
        databases.insert(
            field(r, c, "external_db_id")?.to_owned(),
            field(r, c, "db_name")?.to_owned(),
        );
        Ok(())
    })?;
    let mut xrefs = HashMap::new();
    let wanted: std::collections::HashSet<_> = genes.values().map(|g| g.display.as_str()).collect();
    rows(&manifest, directory, "xref", |r, c| {
        let db = databases
            .get(field(r, c, "external_db_id")?)
            .context("Unknown xref database")?;
        let id = field(r, c, "xref_id")?;
        if wanted.contains(id) || matches!(db.as_str(), "CCDS" | "RefSeq_gene_name") {
            xrefs.insert(
                id.to_owned(),
                (
                    db.clone(),
                    field(r, c, "dbprimary_acc")?.to_owned(),
                    field(r, c, "display_label")?.to_owned(),
                ),
            );
        }
        Ok(())
    })?;
    let mut ccds = HashMap::new();
    let mut fallback_symbols = HashMap::new();
    rows(&manifest, directory, "object_xref", |r, c| {
        if let Some((db, _, label)) = xrefs.get(field(r, c, "xref_id")?) {
            let id = field(r, c, "ensembl_id")?.to_owned();
            match (field(r, c, "ensembl_object_type")?, db.as_str()) {
                ("Transcript", "CCDS") => {
                    ensure!(
                        ccds.insert(id, label.clone()).is_none(),
                        "Multiple CCDS links require explicit ordering"
                    );
                }
                ("Gene", "RefSeq_gene_name") => {
                    fallback_symbols.entry(id).or_insert_with(|| label.clone());
                }
                _ => {}
            }
        }
        Ok(())
    })?;
    let mut translations = HashMap::new();
    rows(&manifest, directory, "translation", |r, c| {
        translations.insert(
            field(r, c, "translation_id")?.to_owned(),
            CoreTranslation {
                stable_id: field(r, c, "stable_id")?.into(),
                version: field(r, c, "version")?.parse()?,
                attributes: vec![],
            },
        );
        Ok(())
    })?;
    rows(&manifest, directory, "translation_attrib", |r, c| {
        let code = attribute_types
            .get(field(r, c, "attrib_type_id")?)
            .context("Unknown translation attribute")?;
        // VEP Translation::get_all_SeqEdits selects exactly these four codes.
        if matches!(
            code.as_str(),
            "initial_met" | "_selenocysteine" | "amino_acid_sub" | "_stop_codon_rt"
        ) {
            translations
                .get_mut(field(r, c, "translation_id")?)
                .context("Missing translation")?
                .attributes
                .push((code.clone(), field(r, c, "value")?.into()));
        }
        Ok(())
    })?;
    let mut transcripts =
        crate::gff::parse_gff3_with_source(MultiGzDecoder::new(File::open(gff3)?), "Ensembl")?;
    let reference = crate::providers::MmapFastaSequenceProvider::new(
        crate::fasta::MmapFastaReader::open(fasta)?,
    );
    let mut enrichment = BTreeMap::new();
    let mut removed = 0;
    let mut kept = Vec::with_capacity(transcripts.len());
    for mut tr in transcripts.drain(..) {
        let id = tr
            .stable_id
            .strip_suffix("_PAR_Y")
            .unwrap_or(tr.stable_id.as_ref());
        let ct = core
            .get(id)
            .with_context(|| format!("GFF transcript {id} missing from public Core"))?;
        ensure!(
            tr.version == Some(ct.version),
            "Transcript version mismatch for {id}"
        );
        // VEP AnnotationSource/Database/Transcript.pm excludes these at construction.
        if tr.biotype.as_ref() == "artifact" || ct.attributes.contains_key("readthrough_tra") {
            removed += 1;
            continue;
        }
        let gene = genes.get(&ct.gene).context("Missing transcript gene")?;
        ensure!(
            tr.gene
                .stable_id
                .strip_suffix("_PAR_Y")
                .unwrap_or(&tr.gene.stable_id)
                == gene.stable_id,
            "Gene identity mismatch for {id}"
        );
        tr.canonical = ct.id == gene.canonical;
        tr.gene.symbol = None;
        tr.gene.symbol_source = None;
        tr.gene.hgnc_id = None;
        if let Some((db, accession, label)) = xrefs.get(&gene.display) {
            tr.gene.symbol = Some(Arc::from(label.as_str()));
            tr.gene.symbol_source = Some(db.clone());
            if db == "HGNC" {
                tr.gene.hgnc_id = Some(accession.clone());
            }
        } else if let Some(symbol) = fallback_symbols.get(&ct.gene) {
            tr.gene.symbol = Some(Arc::from(symbol.as_str()));
        }
        tr.ccds = ccds.get(&ct.id).cloned();
        tr.appris = one(&ct.attributes, "appris")?
            .map(|v| v.replace("principal", "P").replace("alternative", "A"));
        // VEP OutputFactory extracts tsl(\d+), ignoring previous-version notes.
        tr.tsl = one(&ct.attributes, "TSL")?
            .and_then(|v| {
                v.strip_prefix("tsl").map(|s| {
                    s.chars()
                        .take_while(char::is_ascii_digit)
                        .collect::<String>()
                })
            })
            .filter(|v| !v.is_empty())
            .map(|v| v.parse())
            .transpose()?;
        tr.mane_select = one(&ct.attributes, "MANE_Select")?;
        tr.mane_plus_clinical = one(&ct.attributes, "MANE_Plus_Clinical")?;
        tr.gencode_primary = ct.attributes.contains_key("gencode_primary");
        tr.flags = ct.flags.clone();
        let mut extra = Enrichment {
            gene_version: Some(gene.version),
            ..Default::default()
        };
        if let Some(ranges) = ct.attributes.get("miRNA") {
            for value in ranges {
                let (s, e) = value
                    .split_once('-')
                    .context("Invalid mature-miRNA range")?;
                extra.mature_mirna_ranges.push((s.parse()?, e.parse()?));
            }
            extra.mature_mirna_ranges.sort_unstable();
            extra.mature_mirna_ranges.dedup();
        }
        if let Some(edits) = ct.attributes.get("_rna_edit") {
            for edit in edits {
                extra.rna_edits.push(parse_edit("_rna_edit", edit)?);
            }
        }
        ensure!(
            extra.rna_edits.is_empty(),
            "RNA coordinate edits require an explicit public-source mapper implementation"
        );
        if let Some(translation) = &ct.translation {
            let tl = translations
                .get(translation)
                .context("Missing canonical translation")?;
            ensure!(
                tr.translation.is_some(),
                "Public translation missing from GFF for {id}"
            );
            tr.protein_id = Some(tl.stable_id.clone());
            tr.protein_version = Some(tl.version);
            for (kind, value) in &tl.attributes {
                extra.translation_edits.push(parse_edit(kind, value)?);
            }
        }
        if tr.is_coding() {
            if let Err(error) = tr.build_sequences(|chrom, s, e| {
                reference
                    .fetch_sequence(chrom, s, e)
                    .map_err(|e| e.to_string())
            }) {
                let chrom = tr.chromosome.strip_prefix("chr").unwrap_or(&tr.chromosome);
                let primary = matches!(chrom, "X" | "Y" | "M" | "MT")
                    || chrom.parse::<u8>().is_ok_and(|n| (1..=22).contains(&n));
                ensure!(
                    !primary && error.contains("not found in FASTA"),
                    "Building sequences for {}: {error}",
                    tr.stable_id
                );
            }
            if let Some(cds) = &tr.translateable_seq {
                tr.peptide = Some(consequence_peptide(
                    cds,
                    &tr.chromosome,
                    &extra.translation_edits,
                )?);
            }
        }
        enrichment.insert(format!("{}:{}", tr.chromosome, tr.stable_id), extra);
        kept.push(tr);
    }
    kept.sort_by(|a, b| {
        (&a.chromosome, a.start, &a.stable_id).cmp(&(&b.chromosome, b.start, &b.stable_id))
    });
    eprintln!(
        "Built {} transcripts; excluded {removed} source-declared artifact/readthrough records",
        kept.len()
    );
    let mut provenance = provenance;
    provenance["builderVersion"] = env!("CARGO_PKG_VERSION").into();
    provenance["builderBinarySha256"] = annocat_cache::sha256(&std::env::current_exe()?)?.into();
    let header = Header {
        species: manifest.species,
        assembly: manifest.assembly,
        ensembl_release: 115,
        vep_release: manifest.vep_release,
        capabilities: vec![
            "complete-transcript-membership".into(),
            "vep115-core-metadata".into(),
            "mature-mirna-ranges".into(),
            "translation-sequence-edits".into(),
        ],
        provenance,
        transcript_count: 0,
        coding_transcript_count: 0,
        contig_counts: BTreeMap::new(),
        payload_bytes: 0,
        payload_sha256: String::new(),
        semantic_sha256: String::new(),
    };
    annocat_cache::save(kept, enrichment, header, output)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn codon_peptide_uses_only_source_declared_edits() {
        assert_eq!(consequence_peptide("CTGTGATAA", "1", &[]).unwrap(), "L**");
        assert_eq!(translated_peptide("CTGTGATAA", "1", &[], true).unwrap(), "M**");
        assert_eq!(translated_peptide("CTGTGATAA", "1",
            &[parse_edit("amino_acid_sub", "1 1 A").unwrap()], true).unwrap(), "A**");
        let edits = [
            parse_edit("initial_met", "1 1 M").unwrap(),
            parse_edit("_selenocysteine", "2 2 U").unwrap(),
        ];
        assert_eq!(
            consequence_peptide("CTGTGATAA", "1", &edits).unwrap(),
            "MU*"
        );
        assert_eq!(consequence_peptide("ATTAGATGA", "MT", &[]).unwrap(), "I*W");
        assert!(consequence_peptide(
            "ATG",
            "1",
            &[parse_edit("amino_acid_sub", "2 2 X").unwrap()]
        )
        .is_err());
    }
}
