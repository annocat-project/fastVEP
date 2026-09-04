use fastvep_cli::pipeline::{run_annotate, AnnotateConfig};
use std::fs;

#[test]
fn reverse_strand_hgvsc_uses_each_alts_minimal_representation() {
    let tmp = tempfile::tempdir().unwrap();
    let gff3 = tmp.path().join("transcripts.gff3");
    let fasta = tmp.path().join("reference.fa");
    let input = tmp.path().join("input.vcf");
    let output = tmp.path().join("output.vcf");

    fs::write(
        &gff3,
        "##gff-version 3\n\
         1\ttest\tgene\t1\t12\t.\t-\t.\tID=gene:ENSG_REV1;Name=REVTEST;gene_name=REVTEST;biotype=protein_coding\n\
         1\ttest\tmRNA\t1\t12\t.\t-\t.\tID=transcript:ENST_REV1;Parent=gene:ENSG_REV1;biotype=protein_coding\n\
         1\ttest\texon\t1\t12\t.\t-\t.\tID=exon:ENSE_REV1;Parent=transcript:ENST_REV1;rank=1\n\
         1\ttest\tCDS\t1\t12\t.\t-\t0\tID=CDS:ENSP_REV1;Parent=transcript:ENST_REV1\n",
    )
    .unwrap();
    fs::write(&fasta, ">1\nACGTCTGACGTA\n").unwrap();
    fs::write(
        &input,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         1\t5\t.\tC\tA,CAA\t.\tPASS\t.\n",
    )
    .unwrap();

    run_annotate(AnnotateConfig {
        input: input.to_string_lossy().into_owned(),
        output: output.to_string_lossy().into_owned(),
        gff3: vec![gff3.to_string_lossy().into_owned()],
        fasta: Some(fasta.to_string_lossy().into_owned()),
        output_format: "vcf".into(),
        buffer_size: 5000,
        pick: false,
        hgvs: true,
        distance: 0,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: Vec::new(),
        sa_only: false,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
        structured_output: None,
        omit_supplementary_vcf: false,
        show_progress: false,
        profile_output: None,
    })
    .unwrap();

    let output = fs::read_to_string(output).unwrap();
    assert!(output.contains("ENST_REV1:c.7_8insTT"), "{output}");
    assert!(!output.contains("delins"), "{output}");
}
