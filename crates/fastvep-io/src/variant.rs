use fastvep_core::{
    Allele, Consequence, GeneAnnotation, GenomicPosition, Impact, Strand, SupplementaryAnnotation,
    VariantType,
};
use std::sync::Arc;

/// A variant feature ready for annotation.
#[derive(Debug, Clone)]
pub struct VariationFeature {
    pub position: GenomicPosition,
    /// Allele string in Ensembl format: "REF/ALT1/ALT2"
    pub allele_string: String,
    /// The reference allele after normalization.
    pub ref_allele: Allele,
    /// Alternative alleles after normalization.
    pub alt_alleles: Vec<Allele>,
    /// Variant ID (e.g., rs number) from VCF ID column.
    pub variation_name: Option<String>,
    /// Original VCF fields for reconstruction.
    pub vcf_fields: Option<VcfFields>,
    /// Transcript-level annotations (populated during annotation).
    pub transcript_variations: Vec<TranscriptVariation>,
    /// Co-located known variants (populated during annotation).
    pub existing_variants: Vec<KnownVariant>,
    /// Whether input minimization has completed. Multiallelic input retains
    /// its shared interval; do not repeat the pass after allele case conversion.
    pub minimised: bool,
    /// Most severe consequence across all transcripts/alleles.
    pub most_severe_consequence: Option<Consequence>,
    /// Classified variant type (SNV, insertion, deletion, SV, etc.).
    pub variant_type: VariantType,
    /// For structural variants: the END position from VCF INFO.
    pub sv_end: Option<u64>,
    /// For structural variants: the SVLEN from VCF INFO.
    pub sv_len: Option<i64>,
    /// Supplementary annotations from external data sources (ClinVar, gnomAD, etc.).
    pub supplementary_annotations: Vec<SupplementaryAnnotation>,
    /// Gene-level annotations (OMIM, pLI, etc.).
    pub gene_annotations: Vec<GeneAnnotation>,
}

/// Parsed VCF fields for output reconstruction.
#[derive(Debug, Clone)]
pub struct VcfFields {
    pub chrom: String,
    pub pos: u64,
    pub id: String,
    pub ref_allele: String,
    pub alt: String,
    pub qual: String,
    pub filter: String,
    pub info: String,
    pub rest: Vec<String>,
}

/// Annotation of a variant allele against a specific transcript.
#[derive(Debug, Clone)]
pub struct TranscriptVariation {
    pub transcript_id: Arc<str>,
    pub gene_id: Arc<str>,
    pub gene_symbol: Option<Arc<str>>,
    pub biotype: Arc<str>,
    pub allele_annotations: Vec<AlleleAnnotation>,
    pub canonical: bool,
    pub strand: Strand,
    pub source: Option<String>,
    pub protein_id: Option<String>,
    pub mane_select: Option<String>,
    pub mane_plus_clinical: Option<String>,
    pub tsl: Option<u8>,
    pub appris: Option<String>,
    pub ccds: Option<String>,
    pub gencode_primary: bool,
    pub symbol_source: Option<String>,
    pub hgnc_id: Option<String>,
    /// Flags like "cds_end_NF", "cds_start_NF"
    pub flags: Vec<String>,
}

/// Annotation for a specific allele against a specific transcript.
#[derive(Debug, Clone)]
pub struct AlleleAnnotation {
    pub allele: Allele,
    pub consequences: Vec<Consequence>,
    pub impact: Impact,
    pub cdna_position: PositionRange,
    pub cds_position: PositionRange,
    pub protein_position: PositionRange,
    pub amino_acids: Option<(String, String)>,
    pub codons: Option<(String, String)>,
    /// First exon, last exon, total exons (1-based transcript order).
    pub exon: Option<(u32, u32, u32)>,
    /// First intron, last intron, total introns (1-based transcript order).
    pub intron: Option<(u32, u32, u32)>,
    pub distance: Option<i64>,
    pub hgvsc: Option<String>,
    pub hgvsp: Option<String>,
    pub hgvsg: Option<String>,
    /// HGVS offset: number of bases shifted during 3' normalization.
    pub hgvs_offset: Option<i64>,
    pub existing_variation: Vec<String>,
    pub sift: Option<String>,
    pub polyphen: Option<String>,
    /// Per-allele supplementary annotations as (json_key, json_value) pairs.
    pub supplementary: Vec<(String, String)>,
    /// ACMG-AMP classification result (serialized as serde_json::Value).
    pub acmg_classification: Option<serde_json::Value>,
}

/// A VEP position range whose first or last endpoint may be unknown.
///
/// VEP writes these as `?-N` or `N-?`. Position zero is not valid in these
/// 1-based fields, so zero compactly represents an unknown endpoint.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PositionRange {
    start: u64,
    end: u64,
}

impl PositionRange {
    pub const fn new(start: Option<u64>, end: Option<u64>) -> Self {
        Self {
            start: Self::encode(start),
            end: Self::encode(end),
        }
    }

    pub const fn complete(start: u64, end: u64) -> Self {
        Self::new(Some(start), Some(end))
    }

    const fn encode(position: Option<u64>) -> u64 {
        match position {
            Some(0) => panic!("VEP positions are 1-based"),
            Some(position) => position,
            None => 0,
        }
    }

    pub const fn start(self) -> Option<u64> {
        if self.start == 0 {
            None
        } else {
            Some(self.start)
        }
    }

    pub const fn end(self) -> Option<u64> {
        if self.end == 0 {
            None
        } else {
            Some(self.end)
        }
    }

    pub const fn first_known(self) -> Option<u64> {
        match self.start() {
            Some(start) => Some(start),
            None => self.end(),
        }
    }
}

#[cfg(test)]
mod position_range_tests {
    use super::PositionRange;

    #[test]
    fn position_range_is_compact() {
        assert_eq!(std::mem::size_of::<PositionRange>(), 16);
        assert!(std::mem::size_of::<PositionRange>() <= std::mem::size_of::<Option<(u64, u64)>>());
    }
}

/// A known/existing variant from the variation cache.
#[derive(Debug, Clone)]
pub struct KnownVariant {
    pub name: String,
    pub allele_string: Option<String>,
    pub minor_allele: Option<String>,
    pub minor_allele_freq: Option<f64>,
    pub clinical_significance: Option<String>,
    pub somatic: bool,
    pub phenotype_or_disease: bool,
    pub pubmed: Vec<String>,
    pub frequencies: std::collections::HashMap<String, f64>,
}

impl VariationFeature {
    /// Map internal alleles back to uploaded VCF keys for allele-matched sources.
    pub fn supplementary_query_alleles(&self) -> Vec<(String, u64, String, String)> {
        if let Some(vcf) = &self.vcf_fields {
            let uploaded_alts: Vec<&str> = vcf.alt.split(',').collect();
            return self
                .alt_alleles
                .iter()
                .enumerate()
                .map(|(idx, allele)| {
                    let allele_string = allele.to_string();
                    (
                        allele_string.clone(),
                        vcf.pos,
                        vcf.ref_allele.clone(),
                        uploaded_alts
                            .get(idx)
                            .copied()
                            .unwrap_or(&allele_string)
                            .to_string(),
                    )
                })
                .collect();
        }

        self.alt_alleles
            .iter()
            .map(|allele| {
                (
                    allele.to_string(),
                    self.position.start,
                    self.ref_allele.to_string(),
                    allele.to_string(),
                )
            })
            .collect()
    }

    /// Compute the most severe consequence across all transcript annotations.
    pub fn compute_most_severe(&mut self) {
        let all_consequences: Vec<Consequence> = self
            .transcript_variations
            .iter()
            .flat_map(|tv| {
                tv.allele_annotations
                    .iter()
                    .flat_map(|aa| aa.consequences.iter().copied())
            })
            .collect();
        self.most_severe_consequence = Consequence::most_severe(&all_consequences);
    }

    /// Check if this is an insertion.
    pub fn is_insertion(&self) -> bool {
        self.ref_allele == Allele::Deletion
    }

    /// Check if this is a deletion.
    pub fn is_deletion(&self) -> bool {
        self.alt_alleles.iter().any(|a| *a == Allele::Deletion)
    }

    /// Check if this is an indel.
    pub fn is_indel(&self) -> bool {
        self.ref_allele == Allele::Deletion
            || self.alt_alleles.iter().any(|a| *a == Allele::Deletion)
            || self
                .alt_alleles
                .iter()
                .any(|a| a.len() != self.ref_allele.len())
    }
}
