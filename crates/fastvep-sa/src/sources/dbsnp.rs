//! dbSNP VCF parser for building .osa annotation files.
//!
//! Parses dbSNP's VCF release to extract RS IDs and global MAF.
//!
//! The full NCBI dbSNP VCF contains ~800 million records. Use
//! `iter_dbsnp_vcf` for streaming builds; `parse_dbsnp_vcf` is retained for
//! tests and small inputs.

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::borrow::Cow;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fmt::Write;
use std::io::BufRead;

/// Stream a coordinate-sorted dbSNP VCF as `AnnotationRecord`s without
/// buffering the whole file in memory.
///
/// The input must already be sorted by chromosome and position (all standard
/// NCBI dbSNP releases are).
pub fn iter_dbsnp_vcf<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
) -> DbsnpRecordIter<'a, R> {
    iter_dbsnp_vcf_selected(reader, chrom_to_idx, None)
}

/// Stream dbSNP records while extracting and serializing only the requested
/// annotation fields. The record identity is still resolved from ID/RS even
/// when `id` is not retained, preserving the original row-selection behavior.
pub fn iter_dbsnp_vcf_selected<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
    selected_fields: Option<&HashSet<String>>,
) -> DbsnpRecordIter<'a, R> {
    let includes = |field| selected_fields.is_none_or(|fields| fields.contains(field));
    DbsnpRecordIter {
        reader,
        line: String::new(),
        chrom_to_idx,
        pending: VecDeque::new(),
        include_id: includes("id"),
        include_global_maf: includes("globalMaf"),
        include_variant_type: includes("variantType"),
        include_common: includes("common"),
    }
}

pub struct DbsnpRecordIter<'a, R: BufRead> {
    reader: R,
    line: String,
    chrom_to_idx: &'a HashMap<String, u16>,
    pending: VecDeque<AnnotationRecord>,
    include_id: bool,
    include_global_maf: bool,
    include_variant_type: bool,
    include_common: bool,
}

impl<R: BufRead> Iterator for DbsnpRecordIter<'_, R> {
    type Item = Result<AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(record) = self.pending.pop_front() {
                return Some(Ok(record));
            }

            self.line.clear();
            match self.reader.read_line(&mut self.line) {
                Ok(0) => return None,
                Ok(_) => {}
                Err(error) => return Some(Err(error).context("Reading dbSNP VCF line")),
            };
            let line = self.line.trim_end_matches(['\r', '\n']);

            if line.starts_with('#') {
                continue;
            }

            let mut fields = line.splitn(9, '\t');
            let Some(chrom_field) = fields.next() else {
                continue;
            };
            let Some(position_field) = fields.next() else {
                continue;
            };
            let Some(id) = fields.next() else {
                continue;
            };
            let Some(ref_field) = fields.next() else {
                continue;
            };
            let Some(alt_field) = fields.next() else {
                continue;
            };
            // QUAL and FILTER are not needed, but consuming them keeps INFO at
            // the same VCF column as the original parser.
            if fields.next().is_none() || fields.next().is_none() {
                continue;
            }
            let Some(info) = fields.next() else {
                continue;
            };

            let chrom = normalize_chrom(chrom_field);
            let chrom_idx = match self.chrom_to_idx.get(&chrom) {
                Some(&idx) => idx,
                None => continue,
            };

            let pos: u32 = match position_field.parse() {
                Ok(p) => p,
                Err(_) => continue,
            };

            let ref_allele = ref_field.to_string();
            let id_is_rs = id.starts_with("rs");
            let mut rs_info = None;
            let mut caf = None;
            let mut variant_type = None;
            let mut common = false;
            for item in info.split(';') {
                if self.include_common && item == "COMMON" {
                    common = true;
                    continue;
                }
                let Some((key, value)) = item.split_once('=') else {
                    continue;
                };
                match key {
                    "RS" if !id_is_rs => rs_info = Some(value),
                    "CAF" if self.include_global_maf => caf = Some(value),
                    "VC" if self.include_variant_type => variant_type = Some(value),
                    _ => {}
                }
            }

            let rs_id = if id_is_rs {
                Cow::Borrowed(id)
            } else if let Some(rs) = rs_info {
                Cow::Owned(format!("rs{rs}"))
            } else {
                continue;
            };

            // dbSNP's CAF is `ref_freq,alt1_freq,alt2_freq,...` in ALT order
            // (i.e. index i+1 is the frequency for the i-th ALT); index it
            // per-alt below rather than once, or every ALT past the first in
            // a multi-allelic record gets the first ALT's frequency.
            let mut caf_parts = caf.map(|value| value.split(','));
            if let Some(parts) = caf_parts.as_mut() {
                let _ = parts.next(); // CAF starts with the reference frequency.
            }

            for alt in alt_field.split(',') {
                let freq = caf_parts
                    .as_mut()
                    .and_then(|parts| parts.next())
                    .and_then(|value| value.parse::<f64>().ok());
                if alt == "." || alt == "*" {
                    continue;
                }
                let mut json = String::with_capacity(96);
                json.push('{');
                let mut first = true;
                if self.include_id {
                    push_json_string(&mut json, "id", &rs_id, &mut first);
                }
                if self.include_global_maf {
                    if let Some(freq) = freq {
                        push_separator(&mut json, &mut first);
                        let _ = write!(json, "\"globalMaf\":{freq:.6e}");
                    }
                }
                if self.include_variant_type {
                    if let Some(variant_type) = variant_type {
                        push_json_string(&mut json, "variantType", variant_type, &mut first);
                    }
                }
                if self.include_common && common {
                    push_separator(&mut json, &mut first);
                    json.push_str("\"common\":true");
                }
                json.push('}');
                self.pending.push_back(AnnotationRecord {
                    chrom_idx,
                    position: pos,
                    ref_allele: ref_allele.clone(),
                    alt_allele: alt.to_string(),
                    json,
                });
            }
        }
    }
}

/// Parse a dbSNP VCF and produce sorted `AnnotationRecord`s.
///
/// Loads all records into memory — suitable for tests and small inputs.
/// For the full NCBI release use `iter_dbsnp_vcf` via the pipeline instead.
pub fn parse_dbsnp_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records: Vec<_> = iter_dbsnp_vcf(reader, chrom_to_idx).collect::<Result<_>>()?;
    records.sort_by(|a, b| {
        a.chrom_idx
            .cmp(&b.chrom_idx)
            .then(a.position.cmp(&b.position))
    });
    Ok(records)
}

fn push_separator(json: &mut String, first: &mut bool) {
    if *first {
        *first = false;
    } else {
        json.push(',');
    }
}

fn push_json_string(json: &mut String, key: &str, value: &str, first: &mut bool) {
    push_separator(json, first);
    let _ = write!(json, "\"{key}\":{}", serde_json::to_string(value).unwrap());
}

fn normalize_chrom(chrom: &str) -> String {
    // NCBI's dbSNP VCF names contigs by RefSeq accession (`NC_000001.11`).
    // Leave those untouched so the lookup hits the accession keys the builder
    // seeds into the chromosome map; mangling them to `chrNC_000001.11` is what
    // produced "0 records parsed" (issue #51).
    if chrom.starts_with("chr") || fastvep_core::looks_like_refseq_accession(chrom) {
        chrom.to_string()
    } else {
        format!("chr{}", chrom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_dbsnp_vcf() {
        let vcf = "\
##fileformat=VCFv4.0
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t10019\trs775809821\tTA\tT\t.\t.\tRS=775809821;CAF=0.9998,0.0002
1\t10039\trs978760828\tA\tC\t.\t.\tRS=978760828
";

        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_dbsnp_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 2);

        assert_eq!(records[0].position, 10019);
        assert!(records[0].json.contains("rs775809821"));
        assert!(records[0].json.contains("globalMaf"));

        assert_eq!(records[1].position, 10039);
        assert!(records[1].json.contains("rs978760828"));
        assert!(!records[1].json.contains("globalMaf")); // No CAF
    }

    #[test]
    fn test_parse_dbsnp_multiallelic_caf_indexed_per_alt() {
        // CAF=ref_freq,alt1_freq,alt2_freq — each ALT must get its own
        // frequency, not the first ALT's frequency repeated for every ALT.
        let vcf = "\
##fileformat=VCFv4.0
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t10019\trs1\tA\tC,T\t.\t.\tRS=1;CAF=0.90,0.08,0.02
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_dbsnp_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 2);

        let c = records.iter().find(|r| r.alt_allele == "C").unwrap();
        assert!(
            c.json.contains("8.000000e-2"),
            "C should get CAF index 1 (0.08): {}",
            c.json
        );

        let t = records.iter().find(|r| r.alt_allele == "T").unwrap();
        assert!(
            t.json.contains("2.000000e-2"),
            "T should get CAF index 2 (0.02), not C's frequency: {}",
            t.json
        );
    }

    #[test]
    fn test_parse_dbsnp_refseq_accessions() {
        // The real NCBI dbSNP release (GCF_000001405.40) names contigs by
        // RefSeq accession, not `1`/`chr1`. Regression for issue #51: these
        // must resolve when the chromosome map carries the accession key.
        let vcf = "\
##fileformat=VCFv4.0
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
NC_000001.11\t10019\trs775809821\tTA\tT\t.\t.\tRS=775809821;CAF=0.9998,0.0002
NC_000023.11\t100\trs1\tA\tG\t.\t.\tRS=1
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("NC_000001.11".to_string(), 0u16);
        chrom_map.insert("NC_000023.11".to_string(), 22u16);

        let records = parse_dbsnp_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 2, "RefSeq-accession contigs were skipped");
        assert_eq!(records[0].chrom_idx, 0);
        assert!(records[0].json.contains("rs775809821"));
        assert_eq!(records[1].chrom_idx, 22);
        assert!(records[1].json.contains("rs1"));
    }

    #[test]
    fn selected_fields_are_extracted_without_post_filtering() {
        let vcf = "1\t10019\trs1\tA\tC\t.\t.\tRS=1;CAF=0.90,0.10;VC=SNV;COMMON\n";
        let map = HashMap::from([("chr1".to_string(), 0)]);
        let selected = HashSet::from(["id".to_string(), "globalMaf".to_string()]);
        let record = iter_dbsnp_vcf_selected(vcf.as_bytes(), &map, Some(&selected))
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(
            serde_json::from_str::<serde_json::Value>(&record.json).unwrap(),
            serde_json::json!({"id": "rs1", "globalMaf": 0.1})
        );
    }
}
