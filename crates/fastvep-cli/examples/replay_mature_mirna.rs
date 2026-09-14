//! Test-only replay of inventoried mature-miRNA ranges through the CLI pipeline.
use anyhow::{ensure, Context, Result};
use fastvep_cache::transcript_cache::load_cache;
use fastvep_cli::pipeline::{run_annotate_with_mature_mirna_ranges, AnnotateConfig};
use std::{collections::HashMap, env, fs::File, path::Path};

fn main() -> Result<()> {
    let args: Vec<_> = env::args().skip(1).collect();
    ensure!(
        args.len() == 6,
        "usage: replay_mature_mirna INPUT OUTPUT CACHE FASTA EXCLUSIONS RANGES_JSON"
    );
    let data: serde_json::Value = serde_json::from_reader(File::open(&args[5])?)?;
    ensure!(
        data["schemaVersion"] == 1 && data["ensemblRelease"] == 115 && data["assembly"] == "GRCh38",
        "wrong overlay identity"
    );
    let entries = data["matureMirnaRanges"]
        .as_object()
        .context("missing ranges")?;
    let mut ranges = HashMap::new();
    let transcripts = load_cache(Path::new(&args[2]))?;
    for tr in &transcripts {
        let Some(entry) = entries.get(tr.stable_id.as_ref()) else {
            continue;
        };
        ensure!(
            tr.biotype.as_ref() == "miRNA" && serde_json::to_value(tr.version)? == entry["version"],
            "transcript mismatch: {}",
            tr.stable_id
        );
        let parsed: Vec<(u64, u64)> = serde_json::from_value(entry["ranges"].clone())?;
        let length: u64 = tr.exons.iter().map(|e| e.end - e.start + 1).sum();
        ensure!(
            parsed.iter().all(|&(s, e)| s > 0 && s <= e && e <= length),
            "invalid ranges: {}",
            tr.stable_id
        );
        ensure!(
            ranges.insert(tr.stable_id.to_string(), parsed).is_none(),
            "duplicate transcript"
        );
    }
    ensure!(
        ranges.len() == entries.len(),
        "overlay contains missing transcripts"
    );
    drop(transcripts);
    let config = AnnotateConfig {
        input: args[0].clone(),
        output: args[1].clone(),
        transcript_cache: Some(args[2].clone()),
        fasta: Some(args[3].clone()),
        gff3: vec![],
        output_format: "vcf".into(),
        buffer_size: 5000,
        pick: false,
        hgvs: true,
        distance: 5000,
        cache_dir: None,
        sa_dir: vec![],
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
    };
    run_annotate_with_mature_mirna_ranges(config, Some(args[4].clone()), ranges)
}
