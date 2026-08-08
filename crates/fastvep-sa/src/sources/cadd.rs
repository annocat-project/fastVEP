//! Streaming parser for CADD score tables.

use crate::common::AnnotationRecord;
use crate::writer_v2::Osa2Metadata;
use anyhow::{anyhow, Context, Result};
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

/// OSA2 metadata that preserves the CADD OSA1 annotation contract.
pub fn cadd_osa2_metadata(assembly: &str) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "CADD".into(),
        version: "1.7".into(),
        assembly: assembly.into(),
        json_key: "cadd".into(),
        match_by_allele: true,
        is_array: false,
        record_list: false,
        is_positional: false,
        chunk_bits: 16,
        description: format!("CADD raw and PHRED scores for {assembly}"),
    }
}

pub fn iter_cadd<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
) -> CaddRecordIter<'a, R> {
    iter_cadd_selected(reader, chrom_to_idx, None)
}

/// Stream CADD records while constructing only the fields selected by the
/// caller. Keeping this at the source adapter avoids serializing a complete
/// JSON object only for the CLI to parse and filter it again.
pub fn iter_cadd_selected<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
    selected_fields: Option<&HashSet<String>>,
) -> CaddRecordIter<'a, R> {
    let include_raw = selected_fields.is_none_or(|fields| fields.contains("raw"));
    let include_phred = selected_fields.is_none_or(|fields| fields.contains("phred"));
    CaddRecordIter {
        reader,
        line: String::new(),
        chrom_to_idx,
        line_number: 0,
        include_raw,
        include_phred,
    }
}

pub struct CaddRecordIter<'a, R: BufRead> {
    reader: R,
    line: String,
    chrom_to_idx: &'a HashMap<String, u16>,
    line_number: u64,
    include_raw: bool,
    include_phred: bool,
}

impl<R: BufRead> Iterator for CaddRecordIter<'_, R> {
    type Item = Result<AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line.clear();
            match self.reader.read_line(&mut self.line) {
                Ok(0) => return None,
                Ok(_) => {}
                Err(error) => return Some(Err(error).context("Reading CADD score line")),
            };
            self.line_number += 1;
            let line = self.line.trim_end_matches(['\r', '\n']);
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let mut fields = line.split('\t');
            let Some(chrom_field) = fields.next() else {
                continue;
            };
            let Some(position_field) = fields.next() else {
                return Some(Err(anyhow!(
                    "CADD line {} has fewer than six columns",
                    self.line_number
                )));
            };
            let Some(ref_allele) = fields.next() else {
                return Some(Err(anyhow!(
                    "CADD line {} has fewer than six columns",
                    self.line_number
                )));
            };
            let Some(alt_allele) = fields.next() else {
                return Some(Err(anyhow!(
                    "CADD line {} has fewer than six columns",
                    self.line_number
                )));
            };
            let Some(raw_field) = fields.next() else {
                return Some(Err(anyhow!(
                    "CADD line {} has fewer than six columns",
                    self.line_number
                )));
            };
            let Some(phred_field) = fields.next() else {
                return Some(Err(anyhow!(
                    "CADD line {} has fewer than six columns",
                    self.line_number
                )));
            };
            let chrom = normalize_chrom(chrom_field);
            let Some(&chrom_idx) = self.chrom_to_idx.get(&chrom) else {
                continue;
            };
            let position = match position_field.parse::<u32>() {
                Ok(position) if position > 0 => position,
                _ => {
                    return Some(Err(anyhow!(
                        "CADD line {} has an invalid position",
                        self.line_number
                    )))
                }
            };
            if ref_allele.is_empty() || alt_allele.is_empty() {
                return Some(Err(anyhow!(
                    "CADD line {} has an empty allele",
                    self.line_number
                )));
            }
            let raw = match raw_field.parse::<f64>() {
                Ok(value) if value.is_finite() => value,
                _ => {
                    return Some(Err(anyhow!(
                        "CADD line {} has an invalid raw score",
                        self.line_number
                    )))
                }
            };
            let phred = match phred_field.parse::<f64>() {
                Ok(value) if value.is_finite() => value,
                _ => {
                    return Some(Err(anyhow!(
                        "CADD line {} has an invalid PHRED score",
                        self.line_number
                    )))
                }
            };
            return Some(Ok(AnnotationRecord {
                chrom_idx,
                position,
                ref_allele: ref_allele.to_string(),
                alt_allele: alt_allele.to_string(),
                json: match (self.include_raw, self.include_phred) {
                    // serde_json's default map order is lexical. Preserve the
                    // existing `phred`, `raw` byte order and serde_json number
                    // formatting while avoiding the old parse/filter/serialize
                    // second pass.
                    (true, true) => serde_json::json!({"phred": phred, "raw": raw}).to_string(),
                    (true, false) => serde_json::json!({"raw": raw}).to_string(),
                    (false, true) => serde_json::json!({"phred": phred}).to_string(),
                    (false, false) => "{}".to_string(),
                },
            }));
        }
    }
}

fn normalize_chrom(chrom: &str) -> String {
    if chrom.starts_with("chr") {
        chrom.to_string()
    } else {
        format!("chr{chrom}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn streams_raw_and_phred_scores_by_allele() {
        let input = b"#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n1\t100\tA\tG\t0.125\t12.4\n";
        let map = HashMap::from([("chr1".to_string(), 0)]);
        let records = iter_cadd(&input[..], &map)
            .collect::<Result<Vec<_>>>()
            .unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].ref_allele, "A");
        assert_eq!(records[0].alt_allele, "G");
        let json: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        assert_eq!(json["raw"], 0.125);
        assert_eq!(json["phred"], 12.4);
    }

    #[test]
    fn malformed_scores_fail_closed() {
        let map = HashMap::from([("chr1".to_string(), 0)]);
        assert!(iter_cadd(&b"1\t100\tA\tG\tbad\t12\n"[..], &map)
            .next()
            .unwrap()
            .is_err());
    }

    #[test]
    fn selected_fields_are_emitted_directly() {
        let input = b"1\t100\tA\tG\t0.125\t12.4\n";
        let map = HashMap::from([("chr1".to_string(), 0)]);
        let selected = HashSet::from(["phred".to_string()]);
        let record = iter_cadd_selected(&input[..], &map, Some(&selected))
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(
            serde_json::from_str::<serde_json::Value>(&record.json).unwrap(),
            serde_json::json!({"phred": 12.4})
        );
    }
}
