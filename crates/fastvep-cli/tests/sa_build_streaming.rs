//! Regression tests for the streaming `.osa` builder (issue #55 / PR #54).
//!
//! These guard behaviors that the streaming path changed relative to the old
//! buffer-then-sort builder:
//!   * gzip/bgzip inputs are detected by magic bytes, not just by extension;
//!   * a build that fails mid-stream (unsorted input) errors clearly and leaves
//!     no corrupt partial `.osa`/`.osa.idx` behind.

use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_cli::pipeline::{run_sa_build, run_sa_build_inputs};
use fastvep_sa::reader::SaReader;
use std::fs;
use std::io::Write;
use std::process::{Command, Stdio};

const GNOMAD_SORTED: &str = "\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=\"af\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"an\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t.\tPASS\tAF=0.01;AN=1000
chr1\t200\t.\tC\tT\t.\tPASS\tAF=0.02;AN=1000
chr2\t150\t.\tG\tA\t.\tPASS\tAF=0.3;AN=1000
";

const GNOMAD_UNSORTED: &str = "\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=\"af\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"an\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t500\t.\tA\tG\t.\tPASS\tAF=0.01;AN=1000
chr1\t100\t.\tC\tT\t.\tPASS\tAF=0.02;AN=1000
";

fn dbnsfp_fixture(positions: [u32; 2]) -> String {
    let mut header = vec!["#chr", "pos(1-based)", "ref", "alt"];
    header.extend_from_slice(fastvep_sa::sources::dbnsfp::CURATED_FIELDS);
    let field_index = |name: &str| header.iter().position(|field| *field == name).unwrap();
    let mut rows = Vec::new();
    for (index, position) in positions.into_iter().enumerate() {
        let mut row = vec![".".to_string(); header.len()];
        row[0] = "1".into();
        row[1] = position.to_string();
        row[2] = if index == 0 { "A" } else { "C" }.into();
        row[3] = if index == 0 { "G" } else { "T" }.into();
        row[field_index("SIFT_score")] = if index == 0 { "0.032" } else { "0.450" }.into();
        row[field_index("SIFT_pred")] = if index == 0 { "D" } else { "T" }.into();
        row[field_index("Polyphen2_HDIV_score")] =
            if index == 0 { "0.998" } else { "0.100" }.into();
        row[field_index("Polyphen2_HDIV_pred")] = if index == 0 { "D" } else { "B" }.into();
        rows.push(row.join("\t"));
    }
    format!("{}\n{}\n", header.join("\t"), rows.join("\n"))
}

fn gzip(data: &str) -> Vec<u8> {
    use flate2::{write::GzEncoder, Compression};
    use std::io::Write;
    let mut enc = GzEncoder::new(Vec::new(), Compression::default());
    enc.write_all(data.as_bytes()).unwrap();
    enc.finish().unwrap()
}

#[test]
fn streaming_build_rejects_unsorted_input_and_leaves_no_partial_files() {
    let tmp = tempfile::tempdir().unwrap();
    let src = tmp.path().join("gnomad.vcf");
    let out = tmp.path().join("db");
    fs::write(&src, GNOMAD_UNSORTED).unwrap();

    let err = run_sa_build(
        "gnomad",
        src.to_str().unwrap(),
        out.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .expect_err("unsorted input must fail the streaming build");

    let msg = err.to_string();
    assert!(msg.contains("not sorted"), "unexpected error: {msg}");
    assert!(
        msg.contains("bcftools sort") || msg.contains("chr1..chr22"),
        "sort error should tell the user how to fix it: {msg}"
    );

    // The streaming writer fills the .osa incrementally, so a mid-stream failure
    // must not leave a truncated database (or a stale index) on disk.
    assert!(
        !out.with_extension("osa").exists(),
        "partial .osa must be removed when the build fails"
    );
    assert!(
        !out.with_extension("osa.idx").exists(),
        "no .osa.idx should be left behind when the build fails"
    );
}

#[test]
fn streaming_build_detects_gzip_by_magic_bytes_not_extension() {
    let tmp = tempfile::tempdir().unwrap();

    // Identical content: one plain .vcf, one gzip-compressed but deliberately
    // named .vcf (no .gz/.bgz extension) — mimicking a locally renamed release.
    let plain = tmp.path().join("plain.vcf");
    let misnamed = tmp.path().join("misnamed.vcf");
    fs::write(&plain, GNOMAD_SORTED).unwrap();
    fs::write(&misnamed, gzip(GNOMAD_SORTED)).unwrap();

    let out_plain = tmp.path().join("plain_db");
    let out_gz = tmp.path().join("gz_db");
    run_sa_build("gnomad", plain.to_str().unwrap(), out_plain.to_str().unwrap(), "GRCh38", None, &[], false).unwrap();
    run_sa_build("gnomad", misnamed.to_str().unwrap(), out_gz.to_str().unwrap(), "GRCh38", None, &[], false).unwrap();

    // Magic-byte detection must transparently decompress the misnamed gzip, so
    // both builds produce byte-identical databases. Before the fix the gzip
    // bytes were read raw and the build silently produced an empty database.
    let plain_osa = fs::read(out_plain.with_extension("osa")).unwrap();
    let gz_osa = fs::read(out_gz.with_extension("osa")).unwrap();
    assert!(
        plain_osa.len() > 10,
        "build should contain records beyond the 10-byte header"
    );
    assert_eq!(
        plain_osa, gz_osa,
        "gzip-by-magic build must match the plain-text build"
    );
}

#[test]
fn custom_vcf_build_accepts_plain_and_gzip_stdin() {
    let tmp = tempfile::tempdir().unwrap();
    for (label, input) in [
        ("plain", GNOMAD_SORTED.as_bytes().to_vec()),
        ("gzip", gzip(GNOMAD_SORTED)),
    ] {
        let output = tmp.path().join(label);
        let mut child = Command::new(env!("CARGO_BIN_EXE_fastvep"))
            .args([
                "sa-build", "--source", "custom_vcf", "--input", "-", "--output",
                output.to_str().unwrap(), "--assembly", "GRCh38", "--name", "fixture",
                "--no-progress",
            ])
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        child.stdin.take().unwrap().write_all(&input).unwrap();
        let result = child.wait_with_output().unwrap();
        assert!(
            result.status.success(),
            "{label} stdin failed: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert!(output.with_extension("osa").is_file());
        assert!(output.with_extension("osa.idx").is_file());
    }
}

#[test]
fn dbnsfp_build_streams_to_a_loadable_osa() {
    let tmp = tempfile::tempdir().unwrap();
    let src = tmp.path().join("dbnsfp.tsv");
    let out = tmp.path().join("dbnsfp");
    fs::write(&src, dbnsfp_fixture([100, 200])).unwrap();

    run_sa_build(
        "dbnsfp",
        src.to_str().unwrap(),
        out.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();

    let reader = SaReader::open(&out.with_extension("osa")).unwrap();
    let hit = reader
        .annotate_position("chr1", 100, "A", "G")
        .unwrap()
        .expect("dbNSFP allele should be present");
    match hit {
        AnnotationValue::Json(json) => {
            assert!(json.contains("\"SIFT_score\":\"0.032\""));
            assert!(json.contains("\"SIFT_pred\":\"D\""));
            assert!(json.contains("\"Polyphen2_HDIV_score\":\"0.998\""));
            assert!(json.contains("\"Polyphen2_HDIV_pred\":\"D\""));
        }
        other => panic!("expected allele-specific dbNSFP JSON, got {other:?}"),
    }
}

#[test]
fn dbnsfp_stream_rejects_unsorted_input_without_partial_files() {
    let tmp = tempfile::tempdir().unwrap();
    let src = tmp.path().join("dbnsfp.tsv");
    let out = tmp.path().join("dbnsfp");
    fs::write(&src, dbnsfp_fixture([500, 100])).unwrap();

    let error = run_sa_build(
        "dbnsfp",
        src.to_str().unwrap(),
        out.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .expect_err("unsorted dbNSFP input must fail");

    assert!(error.to_string().contains("not sorted"));
    assert!(!out.with_extension("osa").exists());
    assert!(!out.with_extension("osa.idx").exists());
}

/// PhyloP/GERP are per-base genome-wide sources and must build via the
/// streaming path (not buffer every position as a `Vec`) — see
/// `iter_phylop_auto`/`iter_score_tsv`/`iter_wigfix` in fastvep-sa. This
/// exercises both on-disk formats (wigFix and plain TSV) end to end
/// through `run_sa_build` and confirms the resulting `.osa` answers
/// positional queries correctly.
#[test]
fn streaming_build_handles_phylop_wigfix_and_gerp_tsv() {
    let tmp = tempfile::tempdir().unwrap();

    let wigfix = "\
fixedStep chrom=chr1 start=100 step=1
0.5
0.75
-1.25
fixedStep chrom=chr2 start=50 step=1
2.0
";
    let phylop_src = tmp.path().join("hg38.phyloP100way.wigFix");
    fs::write(&phylop_src, wigfix).unwrap();
    let phylop_out = tmp.path().join("phylop");
    run_sa_build(
        "phylop",
        phylop_src.to_str().unwrap(),
        phylop_out.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();

    let phylop_reader = SaReader::open(&phylop_out.with_extension("osa")).unwrap();
    let hit = phylop_reader
        .annotate_position("chr1", 101, "", "")
        .unwrap()
        .expect("position 101 on chr1 should have a phyloP score");
    match hit {
        AnnotationValue::Positional(v) => assert!(v.contains("0.75"), "unexpected phyloP value: {v}"),
        other => panic!("expected a positional phyloP value, got {other:?}"),
    }
    assert!(
        phylop_reader
            .annotate_position("chr2", 50, "", "")
            .unwrap()
            .is_some(),
        "second fixedStep block (chr2) should also be present"
    );

    let gerp_tsv = "\
chr1\t100\t1.234
chr1\t200\t-0.5
";
    let gerp_src = tmp.path().join("gerp_scores.tsv");
    fs::write(&gerp_src, gerp_tsv).unwrap();
    let gerp_out = tmp.path().join("gerp");
    run_sa_build(
        "gerp",
        gerp_src.to_str().unwrap(),
        gerp_out.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();

    let gerp_reader = SaReader::open(&gerp_out.with_extension("osa")).unwrap();
    let hit = gerp_reader
        .annotate_position("chr1", 100, "", "")
        .unwrap()
        .expect("position 100 on chr1 should have a GERP score");
    match hit {
        AnnotationValue::Positional(v) => assert!(v.contains("1.234"), "unexpected GERP value: {v}"),
        other => panic!("expected a positional GERP value, got {other:?}"),
    }
}

#[test]
fn fastvep_merges_multiple_raw_cadd_artifacts() {
    let tmp = tempfile::tempdir().unwrap();
    let snv = tmp.path().join("snv.part");
    let indel = tmp.path().join("indel.part");
    let out = tmp.path().join("cadd");
    fs::write(
        &snv,
        "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n1\t100\tA\tG\t0.1\t10.0\n1\t200\tC\tT\t0.3\t30.0\n",
    )
    .unwrap();
    fs::write(
        &indel,
        "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n1\t150\tAT\tA\t0.2\t20.0\n",
    )
    .unwrap();

    run_sa_build_inputs(
        "cadd",
        &[
            snv.to_string_lossy().into_owned(),
            indel.to_string_lossy().into_owned(),
        ],
        &[0, 0],
        Some("1"),
        out.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();

    let reader = SaReader::open(&out.with_extension("osa")).unwrap();
    for (position, reference, alternate) in
        [(100, "A", "G"), (150, "AT", "A"), (200, "C", "T")]
    {
        assert!(
            reader
                .annotate_position("chr1", position, reference, alternate)
                .unwrap()
                .is_some(),
            "merged CADD cache is missing {position}:{reference}>{alternate}"
        );
    }
}

#[test]
fn fastvep_owns_range_skip_and_chromosome_filtering() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("range.part");
    let output = tmp.path().join("gnomad_chr1");
    let prefix = "discard this decompressed prefix\n";
    fs::write(&input, gzip(&format!("{prefix}{GNOMAD_SORTED}"))).unwrap();

    run_sa_build_inputs(
        "gnomad",
        &[input.to_string_lossy().into_owned()],
        &[prefix.len() as u64],
        Some("chr1"),
        output.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();

    let reader = SaReader::open(&output.with_extension("osa")).unwrap();
    assert!(
        reader
            .annotate_position("chr1", 100, "A", "G")
            .unwrap()
            .is_some()
    );
    assert!(
        reader
            .annotate_position("chr2", 150, "G", "A")
            .unwrap()
            .is_none()
    );
}
