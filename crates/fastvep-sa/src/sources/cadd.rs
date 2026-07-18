//! Streaming parser for CADD score tables.

use crate::common::AnnotationRecord;
use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

pub fn iter_cadd<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
) -> CaddRecordIter<'a, R> {
    CaddRecordIter {
        lines: reader.lines(),
        chrom_to_idx,
        line_number: 0,
    }
}

pub struct CaddRecordIter<'a, R: BufRead> {
    lines: std::io::Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    line_number: u64,
}

impl<R: BufRead> Iterator for CaddRecordIter<'_, R> {
    type Item = Result<AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line = match self.lines.next()? {
                Ok(line) => line,
                Err(error) => return Some(Err(error).context("Reading CADD score line")),
            };
            self.line_number += 1;
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            if fields.len() < 6 {
                return Some(Err(anyhow!("CADD line {} has fewer than six columns", self.line_number)));
            }
            let chrom = normalize_chrom(fields[0]);
            let Some(&chrom_idx) = self.chrom_to_idx.get(&chrom) else { continue; };
            let position = match fields[1].parse::<u32>() {
                Ok(position) if position > 0 => position,
                _ => return Some(Err(anyhow!("CADD line {} has an invalid position", self.line_number))),
            };
            if fields[2].is_empty() || fields[3].is_empty() {
                return Some(Err(anyhow!("CADD line {} has an empty allele", self.line_number)));
            }
            let raw = match fields[4].parse::<f64>() {
                Ok(value) if value.is_finite() => value,
                _ => return Some(Err(anyhow!("CADD line {} has an invalid raw score", self.line_number))),
            };
            let phred = match fields[5].parse::<f64>() {
                Ok(value) if value.is_finite() => value,
                _ => return Some(Err(anyhow!("CADD line {} has an invalid PHRED score", self.line_number))),
            };
            return Some(Ok(AnnotationRecord {
                chrom_idx,
                position,
                ref_allele: fields[2].to_string(),
                alt_allele: fields[3].to_string(),
                json: serde_json::json!({"raw": raw, "phred": phred}).to_string(),
            }));
        }
    }
}

fn normalize_chrom(chrom: &str) -> String {
    if chrom.starts_with("chr") { chrom.to_string() } else { format!("chr{chrom}") }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn streams_raw_and_phred_scores_by_allele() {
        let input = b"#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n1\t100\tA\tG\t0.125\t12.4\n";
        let map = HashMap::from([("chr1".to_string(), 0)]);
        let records = iter_cadd(&input[..], &map).collect::<Result<Vec<_>>>().unwrap();
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
        assert!(iter_cadd(&b"1\t100\tA\tG\tbad\t12\n"[..], &map).next().unwrap().is_err());
    }
}
