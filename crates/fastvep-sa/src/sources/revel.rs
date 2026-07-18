//! REVEL score parser for building .osa annotation files.
//!
//! REVEL provides missense pathogenicity predictions as allele-specific scores.
//! Input format: CSV with columns chr, hg19_pos, grch38_pos, ref, alt, REVEL.

use crate::common::AnnotationRecord;
use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Parse a REVEL score file (CSV) into sorted AnnotationRecords.
///
/// REVEL distributes scores as CSV: chr, hg19_pos, grch38_pos, ref, alt, aaref, aaalt, REVEL
/// We use grch38_pos (column index 2) by default.
pub fn parse_revel<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
    pos_column: usize,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = iter_revel(reader, chrom_to_idx, pos_column).collect::<Result<Vec<_>>>()?;

    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

/// Stream coordinate-sorted REVEL CSV rows without retaining a chromosome in memory.
pub fn iter_revel<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
    pos_column: usize,
) -> RevelRecordIter<'a, R> {
    RevelRecordIter { lines: reader.lines(), chrom_to_idx, pos_column, line_number: 0 }
}

pub struct RevelRecordIter<'a, R: BufRead> {
    lines: std::io::Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    pos_column: usize,
    line_number: u64,
}

impl<R: BufRead> Iterator for RevelRecordIter<'_, R> {
    type Item = Result<AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line = match self.lines.next()? {
                Ok(line) => line,
                Err(error) => return Some(Err(error).context("Reading REVEL line")),
            };
            self.line_number += 1;
            if line.starts_with("chr,") || line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let fields = line.split(',').collect::<Vec<_>>();
            if fields.len() < 8 {
                return Some(Err(anyhow!("REVEL line {} has fewer than eight columns", self.line_number)));
            }
            let chrom = normalize_chrom(fields[0]);
            let Some(&chrom_idx) = self.chrom_to_idx.get(&chrom) else { continue; };
            let position = match fields.get(self.pos_column).and_then(|value| value.parse::<u32>().ok()) {
                Some(position) if position > 0 => position,
                _ => return Some(Err(anyhow!("REVEL line {} has an invalid GRCh38 position", self.line_number))),
            };
            if fields[3].is_empty() || fields[4].is_empty() {
                return Some(Err(anyhow!("REVEL line {} has an empty allele", self.line_number)));
            }
            let score = match fields[7].trim().parse::<f64>() {
                Ok(score) if score.is_finite() && (0.0..=1.0).contains(&score) => score,
                _ => return Some(Err(anyhow!("REVEL line {} has an invalid score", self.line_number))),
            };
            let transcript_id = fields.get(8).copied().unwrap_or_default();
            return Some(Ok(AnnotationRecord {
                chrom_idx,
                position,
                ref_allele: fields[3].to_string(),
                alt_allele: fields[4].to_string(),
                json: serde_json::json!({
                    "score": score,
                    "transcriptId": transcript_id,
                    "aaRef": fields.get(5).copied().unwrap_or_default(),
                    "aaAlt": fields.get(6).copied().unwrap_or_default(),
                }).to_string(),
            }));
        }
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
    fn test_parse_revel() {
        let data = "\
chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL
1,35142,35142,G,A,T,M,0.027
1,35142,35142,G,C,T,S,0.035
1,35143,35143,C,A,T,N,0.842
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".into(), 0u16);

        let records = parse_revel(data.as_bytes(), &chrom_map, 2).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].position, 35142);
        assert_eq!(records[0].ref_allele, "G");
        assert_eq!(records[0].alt_allele, "A");
        assert!(records[0].json.contains("0.027"));
        assert_eq!(records[2].position, 35143);
        assert!(records[2].json.contains("0.842"));
    }

    #[test]
    fn malformed_rows_fail_closed() {
        let map = HashMap::from([("chr1".to_string(), 0)]);
        assert!(iter_revel(&b"1,1,bad,G,A,T,M,0.2,ENST1\n"[..], &map, 2)
            .next().unwrap().is_err());
    }
}
