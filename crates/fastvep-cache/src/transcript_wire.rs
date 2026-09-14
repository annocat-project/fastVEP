//! Frozen bincode layouts. Change the format version, never these field orders.
use fastvep_core::Strand;
use fastvep_genome as model;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gene {
    pub stable_id: Arc<str>,
    pub symbol: Option<Arc<str>>,
    pub symbol_source: Option<String>,
    pub hgnc_id: Option<String>,
    pub biotype: Arc<str>,
    pub chromosome: Arc<str>,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
}

impl From<model::Gene> for Gene {
    fn from(v: model::Gene) -> Self {
        Self {
            stable_id: v.stable_id,
            symbol: v.symbol,
            symbol_source: v.symbol_source,
            hgnc_id: v.hgnc_id,
            biotype: v.biotype,
            chromosome: v.chromosome,
            start: v.start,
            end: v.end,
            strand: v.strand,
        }
    }
}

impl From<Gene> for model::Gene {
    fn from(v: Gene) -> Self {
        Self {
            stable_id: v.stable_id,
            symbol: v.symbol,
            symbol_source: v.symbol_source,
            hgnc_id: v.hgnc_id,
            biotype: v.biotype,
            chromosome: v.chromosome,
            start: v.start,
            end: v.end,
            strand: v.strand,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Exon {
    pub stable_id: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub phase: i8,
    pub end_phase: i8,
    pub rank: u32,
}

impl From<model::Exon> for Exon {
    fn from(v: model::Exon) -> Self {
        Self {
            stable_id: v.stable_id,
            start: v.start,
            end: v.end,
            strand: v.strand,
            phase: v.phase,
            end_phase: v.end_phase,
            rank: v.rank,
        }
    }
}

impl From<Exon> for model::Exon {
    fn from(v: Exon) -> Self {
        Self {
            stable_id: v.stable_id,
            start: v.start,
            end: v.end,
            strand: v.strand,
            phase: v.phase,
            end_phase: v.end_phase,
            rank: v.rank,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Translation {
    pub stable_id: String,
    pub genomic_start: u64,
    pub genomic_end: u64,
    pub start_exon_rank: u32,
    pub start_exon_offset: u64,
    pub end_exon_rank: u32,
    pub end_exon_offset: u64,
}

impl From<model::Translation> for Translation {
    fn from(v: model::Translation) -> Self {
        Self {
            stable_id: v.stable_id,
            genomic_start: v.genomic_start,
            genomic_end: v.genomic_end,
            start_exon_rank: v.start_exon_rank,
            start_exon_offset: v.start_exon_offset,
            end_exon_rank: v.end_exon_rank,
            end_exon_offset: v.end_exon_offset,
        }
    }
}

impl From<Translation> for model::Translation {
    fn from(v: Translation) -> Self {
        Self {
            stable_id: v.stable_id,
            genomic_start: v.genomic_start,
            genomic_end: v.genomic_end,
            start_exon_rank: v.start_exon_rank,
            start_exon_offset: v.start_exon_offset,
            end_exon_rank: v.end_exon_rank,
            end_exon_offset: v.end_exon_offset,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Transcript {
    pub stable_id: Arc<str>,
    pub version: Option<u32>,
    pub gene: Gene,
    pub biotype: Arc<str>,
    pub chromosome: Arc<str>,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub exons: Vec<Exon>,
    pub translation: Option<Translation>,
    pub cdna_coding_start: Option<u64>,
    pub cdna_coding_end: Option<u64>,
    pub coding_region_start: Option<u64>,
    pub coding_region_end: Option<u64>,
    pub spliced_seq: Option<String>,
    pub translateable_seq: Option<String>,
    pub peptide: Option<String>,
    pub canonical: bool,
    pub mane_select: Option<String>,
    pub mane_plus_clinical: Option<String>,
    pub tsl: Option<u8>,
    pub appris: Option<String>,
    pub ccds: Option<String>,
    pub protein_id: Option<String>,
    pub protein_version: Option<u32>,
    pub swissprot: Vec<String>,
    pub trembl: Vec<String>,
    pub uniparc: Vec<String>,
    pub refseq_id: Option<String>,
    pub source: Option<String>,
    pub gencode_primary: bool,
    pub flags: Vec<String>,
    pub codon_table_start_phase: u64,
}

impl From<model::Transcript> for Transcript {
    fn from(v: model::Transcript) -> Self {
        Self {
            stable_id: v.stable_id,
            version: v.version,
            gene: v.gene.into(),
            biotype: v.biotype,
            chromosome: v.chromosome,
            start: v.start,
            end: v.end,
            strand: v.strand,
            exons: v.exons.into_iter().map(Into::into).collect(),
            translation: v.translation.map(Into::into),
            cdna_coding_start: v.cdna_coding_start,
            cdna_coding_end: v.cdna_coding_end,
            coding_region_start: v.coding_region_start,
            coding_region_end: v.coding_region_end,
            spliced_seq: v.spliced_seq,
            translateable_seq: v.translateable_seq,
            peptide: v.peptide,
            canonical: v.canonical,
            mane_select: v.mane_select,
            mane_plus_clinical: v.mane_plus_clinical,
            tsl: v.tsl,
            appris: v.appris,
            ccds: v.ccds,
            protein_id: v.protein_id,
            protein_version: v.protein_version,
            swissprot: v.swissprot,
            trembl: v.trembl,
            uniparc: v.uniparc,
            refseq_id: v.refseq_id,
            source: v.source,
            gencode_primary: v.gencode_primary,
            flags: v.flags,
            codon_table_start_phase: v.codon_table_start_phase,
        }
    }
}

impl From<Transcript> for model::Transcript {
    fn from(v: Transcript) -> Self {
        Self {
            stable_id: v.stable_id,
            version: v.version,
            gene: v.gene.into(),
            biotype: v.biotype,
            chromosome: v.chromosome,
            start: v.start,
            end: v.end,
            strand: v.strand,
            exons: v.exons.into_iter().map(Into::into).collect(),
            translation: v.translation.map(Into::into),
            cdna_coding_start: v.cdna_coding_start,
            cdna_coding_end: v.cdna_coding_end,
            coding_region_start: v.coding_region_start,
            coding_region_end: v.coding_region_end,
            spliced_seq: v.spliced_seq,
            translateable_seq: v.translateable_seq,
            peptide: v.peptide,
            canonical: v.canonical,
            mane_select: v.mane_select,
            mane_plus_clinical: v.mane_plus_clinical,
            tsl: v.tsl,
            appris: v.appris,
            ccds: v.ccds,
            protein_id: v.protein_id,
            protein_version: v.protein_version,
            swissprot: v.swissprot,
            trembl: v.trembl,
            uniparc: v.uniparc,
            refseq_id: v.refseq_id,
            source: v.source,
            gencode_primary: v.gencode_primary,
            flags: v.flags,
            codon_table_start_phase: v.codon_table_start_phase,
            reference_peptide: None,
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct TranscriptWithoutPhase {
    pub stable_id: Arc<str>,
    pub version: Option<u32>,
    pub gene: Gene,
    pub biotype: Arc<str>,
    pub chromosome: Arc<str>,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub exons: Vec<Exon>,
    pub translation: Option<Translation>,
    pub cdna_coding_start: Option<u64>,
    pub cdna_coding_end: Option<u64>,
    pub coding_region_start: Option<u64>,
    pub coding_region_end: Option<u64>,
    pub spliced_seq: Option<String>,
    pub translateable_seq: Option<String>,
    pub peptide: Option<String>,
    pub canonical: bool,
    pub mane_select: Option<String>,
    pub mane_plus_clinical: Option<String>,
    pub tsl: Option<u8>,
    pub appris: Option<String>,
    pub ccds: Option<String>,
    pub protein_id: Option<String>,
    pub protein_version: Option<u32>,
    pub swissprot: Vec<String>,
    pub trembl: Vec<String>,
    pub uniparc: Vec<String>,
    pub refseq_id: Option<String>,
    pub source: Option<String>,
    pub gencode_primary: bool,
    pub flags: Vec<String>,
}

impl From<TranscriptWithoutPhase> for model::Transcript {
    fn from(v: TranscriptWithoutPhase) -> Self {
        Self {
            stable_id: v.stable_id,
            version: v.version,
            gene: v.gene.into(),
            biotype: v.biotype,
            chromosome: v.chromosome,
            start: v.start,
            end: v.end,
            strand: v.strand,
            exons: v.exons.into_iter().map(Into::into).collect(),
            translation: v.translation.map(Into::into),
            cdna_coding_start: v.cdna_coding_start,
            cdna_coding_end: v.cdna_coding_end,
            coding_region_start: v.coding_region_start,
            coding_region_end: v.coding_region_end,
            spliced_seq: v.spliced_seq,
            translateable_seq: v.translateable_seq,
            peptide: v.peptide,
            canonical: v.canonical,
            mane_select: v.mane_select,
            mane_plus_clinical: v.mane_plus_clinical,
            tsl: v.tsl,
            appris: v.appris,
            ccds: v.ccds,
            protein_id: v.protein_id,
            protein_version: v.protein_version,
            swissprot: v.swissprot,
            trembl: v.trembl,
            uniparc: v.uniparc,
            refseq_id: v.refseq_id,
            source: v.source,
            gencode_primary: v.gencode_primary,
            flags: v.flags,
            codon_table_start_phase: 0,
            reference_peptide: None,
        }
    }
}
