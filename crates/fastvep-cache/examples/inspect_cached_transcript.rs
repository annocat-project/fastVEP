use anyhow::Result;

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    anyhow::ensure!(args.len() >= 3, "usage: inspect_cached_transcript CACHE TRANSCRIPT...");
    let cache = fastvep_cache::annocat_cache::load(std::path::Path::new(&args[1]))?;
    for tr in cache.transcripts.iter().filter(|tr| args[2..].iter().any(|id| id == tr.stable_id.as_ref())) {
        println!("{}", serde_json::json!({"transcript": tr.stable_id, "peptide": tr.peptide,
            "cds": tr.translateable_seq, "phase": tr.codon_table_start_phase, "flags": tr.flags,
            "enrichment": cache.enrichment.get(&format!("{}:{}", tr.chromosome, tr.stable_id))}));
    }
    Ok(())
}
