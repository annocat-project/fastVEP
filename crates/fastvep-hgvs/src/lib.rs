mod coding;
mod genomic;
mod protein;

pub use coding::{
    hgvsc, hgvsc_intronic, hgvsc_intronic_range, hgvsc_noncoding, hgvsc_noncoding_intronic,
    hgvsc_noncoding_intronic_range, hgvsc_with_seq,
};
pub use genomic::hgvsg;
pub use protein::{
    cds_insertion_point,
    hgvsp, hgvsp_frameshift, hgvsp_frameshift_from_cds, hgvsp_frameshift_from_cds_with_tables,
    hgvsp_frameshift_from_cds_with_tables_and_ref_peptide, hgvsp_frameshift_from_cds_with_context,
    hgvsp_inframe_deletion_from_cds,
    hgvsp_inframe_indel, hgvsp_inframe_insertion_from_cds, hgvsp_shifted_stop_retained_insertion,
    hgvsp_inframe_indel_with_context, hgvsp_inframe_insertion_from_cds_with_start_lost,
    hgvsp_start_lost, hgvsp_stop_lost_from_cds, hgvsp_stop_lost_suffix_from_cds,
};

/// Full HGVS annotation result.
#[derive(Debug, Clone, Default)]
pub struct HgvsAnnotation {
    pub hgvsc: Option<String>,
    pub hgvsp: Option<String>,
    pub hgvsg: Option<String>,
}
