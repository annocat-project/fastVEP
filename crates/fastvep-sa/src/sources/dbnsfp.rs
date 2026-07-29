//! dbNSFP parser for building .osa annotation files.
//!
//! dbNSFP provides pre-computed functional predictions (SIFT, PolyPhen,
//! REVEL, CADD, etc.) for all possible missense variants.
//!
//! AnnoCat's pinned build retains a versioned, curated dbNSFP 4.9a field set.
//! Values are preserved losslessly as strings (including transcript-aligned
//! semicolon lists) while the iterator keeps memory bounded.

use crate::common::{escape_json, AnnotationRecord};
use crate::writer_v2::Osa2Metadata;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

// Wide transcript-vector records can exceed the reader's 256 MiB
// decompression limit in a 1 MiB genomic window.
const DBNSFP_CHUNK_BITS: u32 = 16;

/// Fields retained by AnnoCat's `dbnsfp-4.9a-annocat-core-v2` contract.
/// Coordinate columns are used as OSA keys and therefore are not duplicated in
/// each record's JSON payload.
pub const CURATED_FIELDS: &[&str] = &[
    "aaref",
    "aaalt",
    "aapos",
    "genename",
    "Ensembl_geneid",
    "Ensembl_transcriptid",
    "Ensembl_proteinid",
    "Uniprot_acc",
    "Uniprot_entry",
    "HGVSc_VEP",
    "HGVSp_VEP",
    "APPRIS",
    "GENCODE_basic",
    "TSL",
    "VEP_canonical",
    "SIFT_score",
    "SIFT_converted_rankscore",
    "SIFT_pred",
    "SIFT4G_score",
    "SIFT4G_converted_rankscore",
    "SIFT4G_pred",
    "Polyphen2_HDIV_score",
    "Polyphen2_HDIV_rankscore",
    "Polyphen2_HDIV_pred",
    "Polyphen2_HVAR_score",
    "Polyphen2_HVAR_rankscore",
    "Polyphen2_HVAR_pred",
    "MutationTaster_score",
    "MutationTaster_converted_rankscore",
    "MutationTaster_pred",
    "MutationAssessor_score",
    "MutationAssessor_rankscore",
    "MutationAssessor_pred",
    "PROVEAN_score",
    "PROVEAN_converted_rankscore",
    "PROVEAN_pred",
    "VEST4_score",
    "VEST4_rankscore",
    "MetaSVM_score",
    "MetaSVM_rankscore",
    "MetaSVM_pred",
    "MetaLR_score",
    "MetaLR_rankscore",
    "MetaLR_pred",
    "MetaRNN_score",
    "MetaRNN_rankscore",
    "MetaRNN_pred",
    "M-CAP_score",
    "M-CAP_rankscore",
    "M-CAP_pred",
    "REVEL_score",
    "REVEL_rankscore",
    "MutPred_score",
    "MutPred_rankscore",
    "MutPred_protID",
    "MutPred_AAchange",
    "MutPred_Top5features",
    "MVP_score",
    "MVP_rankscore",
    "gMVP_score",
    "gMVP_rankscore",
    "MPC_score",
    "MPC_rankscore",
    "PrimateAI_score",
    "PrimateAI_rankscore",
    "PrimateAI_pred",
    "DEOGEN2_score",
    "DEOGEN2_rankscore",
    "DEOGEN2_pred",
    "BayesDel_noAF_score",
    "BayesDel_noAF_rankscore",
    "BayesDel_noAF_pred",
    "ClinPred_score",
    "ClinPred_rankscore",
    "ClinPred_pred",
    "LIST-S2_score",
    "LIST-S2_rankscore",
    "LIST-S2_pred",
    "VARITY_R_score",
    "VARITY_R_rankscore",
    "VARITY_ER_score",
    "VARITY_ER_rankscore",
    "ESM1b_score",
    "ESM1b_rankscore",
    "ESM1b_pred",
    "EVE_score",
    "EVE_rankscore",
    "AlphaMissense_score",
    "AlphaMissense_rankscore",
    "AlphaMissense_pred",
    "PHACTboost_score",
    "PHACTboost_rankscore",
    "MutFormer_score",
    "MutFormer_rankscore",
    "MutScore_score",
    "MutScore_rankscore",
    "Aloft_Fraction_transcripts_affected",
    "Aloft_prob_Tolerant",
    "Aloft_prob_Recessive",
    "Aloft_prob_Dominant",
    "Aloft_pred",
    "Aloft_Confidence",
    "CADD_raw",
    "CADD_raw_rankscore",
    "CADD_phred",
    "DANN_score",
    "DANN_rankscore",
    "fathmm-XF_coding_score",
    "fathmm-XF_coding_rankscore",
    "fathmm-XF_coding_pred",
    "Eigen-raw_coding",
    "Eigen-raw_coding_rankscore",
    "Eigen-phred_coding",
    "Eigen-PC-raw_coding",
    "Eigen-PC-raw_coding_rankscore",
    "Eigen-PC-phred_coding",
    "GERP++_NR",
    "GERP++_RS",
    "GERP++_RS_rankscore",
    "GERP_91_mammals",
    "GERP_91_mammals_rankscore",
    "phyloP100way_vertebrate",
    "phyloP100way_vertebrate_rankscore",
    "phyloP470way_mammalian",
    "phyloP470way_mammalian_rankscore",
    "phastCons100way_vertebrate",
    "phastCons100way_vertebrate_rankscore",
    "phastCons470way_mammalian",
    "phastCons470way_mammalian_rankscore",
    "SiPhy_29way_logOdds",
    "SiPhy_29way_logOdds_rankscore",
    "bStatistic",
    "bStatistic_converted_rankscore",
    "Interpro_domain",
];

fn selected_fields_from(encoded: Option<&str>) -> Result<Vec<&'static str>> {
    let Some(encoded) = encoded else {
        return Ok(CURATED_FIELDS.to_vec());
    };
    let requested: Vec<String> = serde_json::from_str(&encoded)
        .context("ANNOCAT_DBNSFP_FIELDS must be a JSON string array")?;
    if requested.is_empty() || requested.len() > CURATED_FIELDS.len() {
        anyhow::bail!("ANNOCAT_DBNSFP_FIELDS has an invalid field count");
    }
    let requested = requested.into_iter().collect::<std::collections::HashSet<_>>();
    if requested.len() > CURATED_FIELDS.len()
        || requested
            .iter()
            .any(|field| !CURATED_FIELDS.contains(&field.as_str()))
    {
        anyhow::bail!("ANNOCAT_DBNSFP_FIELDS contains an unknown field");
    }
    Ok(CURATED_FIELDS
        .iter()
        .copied()
        .filter(|field| requested.contains(*field))
        .collect())
}

fn selected_fields() -> Result<Vec<&'static str>> {
    let encoded = std::env::var("ANNOCAT_DBNSFP_FIELDS").ok();
    selected_fields_from(encoded.as_deref())
}

/// Standard dbNSFP `.osa2` metadata. dbNSFP's payload is composite prediction
/// strings (`{"sift":"D(0.012)","polyphen":..}`) that don't decompose into
/// numeric u32 columns, so it is stored as a whole-record JSON blob (see
/// [`crate::writer_v2::raw_json_blob_fields`]): byte-identical to v1, with v2's
/// chunk-level zstd shrinking the database.
pub fn dbnsfp_osa2_metadata(assembly: &str) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "dbNSFP".into(),
        version: "latest".into(),
        assembly: assembly.into(),
        json_key: "dbnsfp".into(),
        match_by_allele: true,
        is_array: false,
        record_list: true,
        is_positional: false,
        chunk_bits: DBNSFP_CHUNK_BITS,
        description: format!("dbNSFP SIFT/PolyPhen predictions for {assembly}"),
    }
}

/// Parse a dbNSFP TSV file using the curated schema.
///
/// The header must contain every curated field. This intentionally fails closed
/// if a different dbNSFP layout is supplied.
pub fn parse_dbnsfp<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records: Vec<_> = iter_dbnsfp(reader, chrom_to_idx).collect::<Result<_>>()?;
    records.sort_by(|a, b| {
        a.chrom_idx
            .cmp(&b.chrom_idx)
            .then(a.position.cmp(&b.position))
    });
    Ok(records)
}

/// Stream coordinate-sorted dbNSFP rows without buffering the input.
///
/// Production dbNSFP chromosome files are already coordinate sorted. The OSA
/// writer validates ordering and fails closed if a caller supplies unsorted
/// input. `parse_dbnsfp` remains as a small-fixture compatibility helper.
pub fn iter_dbnsfp<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
) -> DbNsfpRecordIter<'a, R> {
    DbNsfpRecordIter {
        lines: reader.lines(),
        chrom_to_idx,
        col_indices: None,
    }
}

pub struct DbNsfpRecordIter<'a, R: BufRead> {
    lines: std::io::Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    col_indices: Option<DbNsfpColumns>,
}

impl<R: BufRead> Iterator for DbNsfpRecordIter<'_, R> {
    type Item = Result<AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line = match self.lines.next()? {
                Ok(line) => line,
                Err(error) => return Some(Err(error).context("Reading dbNSFP line")),
            };

            if line.starts_with('#') || line.starts_with("chr\t") {
                let header = line.trim_start_matches('#');
                match DbNsfpColumns::from_header(header) {
                    Ok(columns) => self.col_indices = Some(columns),
                    Err(error) => return Some(Err(error)),
                }
                continue;
            }

            if line.is_empty() {
                continue;
            }

            let cols = match &self.col_indices {
                Some(columns) => columns,
                None => continue,
            };

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() <= cols.max_idx() {
                continue;
            }

            let chrom = normalize_chrom(fields[cols.chr]);
            let chrom_idx = match self.chrom_to_idx.get(&chrom) {
                Some(&idx) => idx,
                None => continue,
            };

            let pos: u32 = match fields[cols.pos].parse() {
                Ok(position) => position,
                Err(_) => continue,
            };

            let ref_allele = fields[cols.ref_col].to_string();
            let alt_allele = fields[cols.alt].to_string();

            let mut parts = Vec::with_capacity(cols.curated.len());
            for &(name, index) in &cols.curated {
                let value = fields[index].trim_end_matches('\r');
                if value.is_empty() || value == "." {
                    continue;
                }
                parts.push(format!(
                    "\"{}\":\"{}\"",
                    escape_json(name),
                    escape_json(value)
                ));
            }

            if parts.is_empty() {
                continue;
            }

            return Some(Ok(AnnotationRecord {
                chrom_idx,
                position: pos,
                ref_allele,
                alt_allele,
                json: format!("{{{}}}", parts.join(",")),
            }));
        }
    }
}

#[derive(Debug)]
struct DbNsfpColumns {
    chr: usize,
    pos: usize,
    ref_col: usize,
    alt: usize,
    curated: Vec<(&'static str, usize)>,
}

impl DbNsfpColumns {
    fn from_header(header: &str) -> Result<Self> {
        let fields: Vec<&str> = header.split('\t').collect();
        let find = |names: &[&str]| -> Option<usize> {
            fields.iter().position(|f| {
                let fl = f.to_lowercase();
                names.iter().any(|n| fl == *n)
            })
        };

        let curated = selected_fields()?
            .into_iter()
            .map(|name| {
                fields
                    .iter()
                    .position(|field| *field == name)
                    .map(|index| (name, index))
                    .with_context(|| {
                        format!("dbNSFP 4.9a header is missing curated field '{name}'")
                    })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self {
            chr: find(&["chr", "#chr"]).context("dbNSFP header is missing chr")?,
            pos: find(&["pos(1-based)", "pos", "hg38_pos"])
                .context("dbNSFP header is missing position")?,
            ref_col: find(&["ref", "ref_allele"]).context("dbNSFP header is missing ref")?,
            alt: find(&["alt", "alt_allele"]).context("dbNSFP header is missing alt")?,
            curated,
        })
    }

    fn max_idx(&self) -> usize {
        self.curated
            .iter()
            .fold(self.alt, |maximum, (_, index)| maximum.max(*index))
    }
}

fn normalize_chrom(chrom: &str) -> String {
    if chrom.starts_with("chr") {
        chrom.to_string()
    } else {
        format!("chr{}", chrom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn curated_fields_are_preserved_without_collapsing_transcript_values() {
        let mut header = vec!["#chr", "pos(1-based)", "ref", "alt"];
        header.extend_from_slice(CURATED_FIELDS);
        let mut row = vec!["."; header.len()];
        row[0] = "1";
        row[1] = "10001";
        row[2] = "A";
        row[3] = "G";
        row[header
            .iter()
            .position(|field| *field == "Ensembl_transcriptid")
            .unwrap()] = "ENST1;ENST2";
        row[header
            .iter()
            .position(|field| *field == "SIFT_score")
            .unwrap()] = "0.032;0.450";
        row[header
            .iter()
            .position(|field| *field == "REVEL_score")
            .unwrap()] = "0.91;0.72";
        let data = format!("{}\n{}\n", header.join("\t"), row.join("\t"));
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".into(), 0u16);

        let records = parse_dbnsfp(data.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 1);
        let json: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        assert_eq!(json["Ensembl_transcriptid"], "ENST1;ENST2");
        assert_eq!(json["SIFT_score"], "0.032;0.450");
        assert_eq!(json["REVEL_score"], "0.91;0.72");
        assert!(json.get("CADD_raw").is_none());
    }

    #[test]
    fn curated_schema_fails_closed_when_a_field_is_missing() {
        let error = DbNsfpColumns::from_header("#chr\tpos(1-based)\tref\talt")
            .unwrap_err()
            .to_string();
        assert!(error.contains("missing curated field 'aaref'"));
    }

    #[test]
    fn configured_subset_is_validated_and_keeps_contract_order() {
        let fields = selected_fields_from(Some(r#"["REVEL_score","SIFT_score"]"#)).unwrap();
        assert_eq!(fields, vec!["SIFT_score", "REVEL_score"]);
        assert!(selected_fields_from(Some(r#"["not_a_dbnsfp_field"]"#)).is_err());
    }

    #[test]
    fn osa2_uses_bounded_chunks_for_wide_dbnsfp_records() {
        assert_eq!(dbnsfp_osa2_metadata("GRCh38").chunk_bits, 16);
    }
}
