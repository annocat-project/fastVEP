use fastvep_cli::pipeline::{run_annotate, AnnotateConfig};
use std::fs;
use std::path::Path;

fn gff3(strand: char) -> String {
    format!(
        "##gff-version 3\n\
         1\ttest\tgene\t1\t24\t.\t{strand}\t.\tID=gene:ENSG_REV1;Name=REVTEST;gene_name=REVTEST;biotype=protein_coding\n\
         1\ttest\tmRNA\t1\t24\t.\t{strand}\t.\tID=transcript:ENST_REV1;Parent=gene:ENSG_REV1;biotype=protein_coding;tag=Ensembl_canonical\n\
         1\ttest\texon\t1\t24\t.\t{strand}\t.\tID=exon:ENSE_REV1;Parent=transcript:ENST_REV1;rank=1\n\
         1\ttest\tCDS\t1\t24\t.\t{strand}\t0\tID=CDS:ENSP_REV1;Parent=transcript:ENST_REV1\n"
    )
}

fn annotate(dir: &Path, strand: char) -> String {
    let gff3_path = dir.join("transcripts.gff3");
    let fasta_path = dir.join("reference.fa");
    let input_path = dir.join("input.vcf");
    let output_path = dir.join("output.vcf");
    let (fasta, vcf) = if strand == '-' {
        (
            ">1\nTTAAGCTTCACCTTCACCTTCCAT\n",
            "1\t12\t.\tCTTCACCTTC\tC\t.\tPASS\t.\n",
        )
    } else {
        (
            ">1\nATGGAAGGTGAAGGTGAAGCTTAA\n",
            "1\t3\t.\tGGAAGGTGAA\tG\t.\tPASS\t.\n",
        )
    };
    fs::write(&gff3_path, gff3(strand)).unwrap();
    fs::write(&fasta_path, fasta).unwrap();
    fs::write(
        &input_path,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{vcf}"),
    )
    .unwrap();

    run_annotate(AnnotateConfig {
        input: input_path.to_string_lossy().into_owned(),
        output: output_path.to_string_lossy().into_owned(),
        gff3: vec![gff3_path.to_string_lossy().into_owned()],
        fasta: Some(fasta_path.to_string_lossy().into_owned()),
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

    fs::read_to_string(output_path).unwrap()
}

#[test]
fn periodic_deletion_uses_the_affected_residues_on_both_strands() {
    for strand in ['+', '-'] {
        let tmp = tempfile::tempdir().unwrap();
        let output = annotate(tmp.path(), strand);
        assert!(output.contains("inframe_deletion"), "{output}");
        assert!(output.contains("|2-4|"), "{output}");
        assert!(output.contains("p.Glu2_Glu4del"), "{output}");
        assert!(!output.contains("p.Glu4_Glu6del"), "{output}");
    }
}
