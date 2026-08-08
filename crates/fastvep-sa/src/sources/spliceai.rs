//! SpliceAI VCF parser for building .osa annotation files.
//!
//! SpliceAI provides splice site effect predictions with delta scores
//! for acceptor/donor gain/loss and their positions.

use crate::common::AnnotationRecord;
use crate::writer_v2::Osa2Metadata;
use anyhow::{Context, Result};
use std::collections::VecDeque;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

/// OSA2 metadata that preserves every gene-specific SpliceAI record per allele.
pub fn spliceai_osa2_metadata(assembly: &str) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "SpliceAI".into(),
        version: "latest".into(),
        assembly: assembly.into(),
        json_key: "spliceAI".into(),
        match_by_allele: true,
        is_array: false,
        record_list: true,
        is_positional: false,
        chunk_bits: 16,
        description: format!("SpliceAI splice site predictions for {assembly}"),
    }
}

/// Parse a SpliceAI VCF and produce sorted AnnotationRecords.
///
/// SpliceAI INFO field format:
/// `SpliceAI=A|GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`
pub fn parse_spliceai_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records: Vec<_> = iter_spliceai_vcf(reader, chrom_to_idx).collect::<Result<_>>()?;
    records.sort_by(|a, b| {
        a.chrom_idx
            .cmp(&b.chrom_idx)
            .then(a.position.cmp(&b.position))
    });
    Ok(records)
}

/// Stream a coordinate-sorted SpliceAI VCF as AnnotationRecords.
///
/// This avoids retaining the whole genome-wide SpliceAI VCF in memory before
/// writing fastSA. The input must already be sorted by chromosome and position.
pub fn iter_spliceai_vcf<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
) -> SpliceAiRecordIter<'a, R> {
    iter_spliceai_vcf_selected(reader, chrom_to_idx, None)
}

/// Stream SpliceAI while emitting only the requested cache fields.
pub fn iter_spliceai_vcf_selected<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
    selected_fields: Option<&HashSet<String>>,
) -> SpliceAiRecordIter<'a, R> {
    let includes = |field| selected_fields.is_none_or(|fields| fields.contains(field));
    SpliceAiRecordIter {
        reader,
        line: String::new(),
        chrom_to_idx,
        pending: VecDeque::new(),
        includes: [
            includes("gene"),
            includes("dsAg"),
            includes("dsAl"),
            includes("dsDg"),
            includes("dsDl"),
            includes("dpAg"),
            includes("dpAl"),
            includes("dpDg"),
            includes("dpDl"),
        ],
    }
}

pub struct SpliceAiRecordIter<'a, R: BufRead> {
    reader: R,
    line: String,
    chrom_to_idx: &'a HashMap<String, u16>,
    pending: VecDeque<AnnotationRecord>,
    includes: [bool; 9],
}

impl<R: BufRead> Iterator for SpliceAiRecordIter<'_, R> {
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
                Err(err) => return Some(Err(err).context("Reading SpliceAI VCF line")),
            }
            let line = self.line.trim_end_matches(['\r', '\n']);
            if line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.splitn(9, '\t').collect();
            if fields.len() < 8 {
                continue;
            }

            let chrom = normalize_chrom(fields[0]);
            let chrom_idx = match self.chrom_to_idx.get(&chrom) {
                Some(&idx) => idx,
                None => continue,
            };

            let pos: u32 = match fields[1].parse() {
                Ok(p) => p,
                Err(_) => continue,
            };

            let ref_allele = fields[3].to_string();
            let info = fields[7];

            for pair in info.split(';') {
                if let Some(val) = pair.strip_prefix("SpliceAI=") {
                    for entry in val.split(',') {
                        if let Some(record) =
                            parse_spliceai_entry(chrom_idx, pos, &ref_allele, entry, &self.includes)
                        {
                            self.pending.push_back(record);
                        }
                    }
                }
            }
        }
    }
}

fn parse_spliceai_entry(
    chrom_idx: u16,
    position: u32,
    ref_allele: &str,
    entry: &str,
    includes: &[bool; 9],
) -> Option<AnnotationRecord> {
    let parts: Vec<&str> = entry.split('|').collect();
    if parts.len() < 10 {
        return None;
    }

    let alt_allele = parts[0].to_string();
    let gene = parts[1];
    let ds_ag: f64 = parts[2].parse().ok()?;
    let ds_al: f64 = parts[3].parse().ok()?;
    let ds_dg: f64 = parts[4].parse().ok()?;
    let ds_dl: f64 = parts[5].parse().ok()?;
    let dp_ag: i32 = parts[6].parse().ok()?;
    let dp_al: i32 = parts[7].parse().ok()?;
    let dp_dg: i32 = parts[8].parse().ok()?;
    let dp_dl: i32 = parts[9].parse().ok()?;

    let mut values = serde_json::Map::new();
    for (include, key, value) in [
        (includes[0], "gene", serde_json::Value::from(gene)),
        (includes[1], "dsAg", serde_json::Value::from(ds_ag)),
        (includes[2], "dsAl", serde_json::Value::from(ds_al)),
        (includes[3], "dsDg", serde_json::Value::from(ds_dg)),
        (includes[4], "dsDl", serde_json::Value::from(ds_dl)),
        (includes[5], "dpAg", serde_json::Value::from(dp_ag)),
        (includes[6], "dpAl", serde_json::Value::from(dp_al)),
        (includes[7], "dpDg", serde_json::Value::from(dp_dg)),
        (includes[8], "dpDl", serde_json::Value::from(dp_dl)),
    ] {
        if include {
            values.insert(key.into(), value);
        }
    }
    let json = serde_json::Value::Object(values).to_string();

    Some(AnnotationRecord {
        chrom_idx,
        position,
        ref_allele: ref_allele.to_string(),
        alt_allele,
        json,
    })
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
    fn test_parse_spliceai_vcf() {
        let vcf = "\
##fileformat=VCFv4.0
##INFO=<ID=SpliceAI,Number=.,Type=String>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t25000\t.\tA\tG\t.\t.\tSpliceAI=G|GENE1|0.01|0.00|0.85|0.00|5|-28|2|-13
1\t30000\t.\tC\tT,A\t.\t.\tSpliceAI=T|GENE2|0.00|0.10|0.00|0.92|3|-5|10|-2,A|GENE2|0.50|0.00|0.00|0.00|7|0|0|0
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".into(), 0u16);

        let records = parse_spliceai_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 3);

        assert_eq!(records[0].position, 25000);
        assert_eq!(records[0].alt_allele, "G");
        let first: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        assert_eq!(first["gene"], "GENE1");
        assert_eq!(first["dsDg"], 0.85);

        assert_eq!(records[1].position, 30000);
        assert_eq!(records[1].alt_allele, "T");
        let second: serde_json::Value = serde_json::from_str(&records[1].json).unwrap();
        assert_eq!(second["dsDl"], 0.92);

        assert_eq!(records[2].position, 30000);
        assert_eq!(records[2].alt_allele, "A");
        let third: serde_json::Value = serde_json::from_str(&records[2].json).unwrap();
        assert_eq!(third["dsAg"], 0.50);
    }

    #[test]
    fn test_iter_spliceai_vcf_streams_records() {
        let vcf = "\
##fileformat=VCFv4.0
##INFO=<ID=SpliceAI,Number=.,Type=String>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t25000\t.\tA\tG\t.\t.\tSpliceAI=G|GENE1|0.01|0.00|0.85|0.00|5|-28|2|-13
1\t30000\t.\tC\tT,A\t.\t.\tSpliceAI=T|GENE2|0.00|0.10|0.00|0.92|3|-5|10|-2,A|GENE2|0.50|0.00|0.00|0.00|7|0|0|0
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".into(), 0u16);

        let records: Vec<_> = iter_spliceai_vcf(vcf.as_bytes(), &chrom_map)
            .collect::<Result<Vec<_>>>()
            .unwrap();

        assert_eq!(records.len(), 3);
        assert_eq!(records[0].position, 25000);
        assert_eq!(records[0].alt_allele, "G");
        assert_eq!(records[1].position, 30000);
        assert_eq!(records[1].alt_allele, "T");
        assert_eq!(records[2].position, 30000);
        assert_eq!(records[2].alt_allele, "A");
    }

    #[test]
    fn test_parse_spliceai_vcf_escapes_gene_for_json() {
        let vcf = "\
##fileformat=VCFv4.0
##INFO=<ID=SpliceAI,Number=.,Type=String>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t25000\t.\tA\tG\t.\t.\tSpliceAI=G|GENE\"1|0.01|0.00|0.85|0.00|5|-28|2|-13
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".into(), 0u16);

        let records = parse_spliceai_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 1);
        let value: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        assert_eq!(value["gene"], "GENE\"1");
    }

    #[test]
    fn selected_fields_are_emitted_directly() {
        let vcf = "1\t25000\t.\tA\tG\t.\t.\tSpliceAI=G|GENE1|0.01|0.00|0.85|0.00|5|-28|2|-13\n";
        let map = HashMap::from([("chr1".to_string(), 0)]);
        let selected = ["gene".to_string(), "dsDg".to_string()]
            .into_iter()
            .collect();
        let record = iter_spliceai_vcf_selected(vcf.as_bytes(), &map, Some(&selected))
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(
            serde_json::from_str::<serde_json::Value>(&record.json).unwrap(),
            serde_json::json!({"gene": "GENE1", "dsDg": 0.85})
        );
    }
}
