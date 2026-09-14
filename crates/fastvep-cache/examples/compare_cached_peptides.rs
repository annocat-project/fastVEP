use anyhow::Result;
use std::collections::HashMap;
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let old = fastvep_cache::annocat_cache::load(std::path::Path::new(&args[1]))?;
    let new = fastvep_cache::annocat_cache::load(std::path::Path::new(&args[2]))?;
    let by_id: HashMap<_, _> = old
        .transcripts
        .iter()
        .map(|tr| {
            (
                (
                    tr.chromosome.trim_start_matches("chr"),
                    tr.stable_id.as_ref(),
                ),
                tr,
            )
        })
        .collect();
    for tr in &new.transcripts {
        if let Some(before) = by_id.get(&(
            tr.chromosome.trim_start_matches("chr"),
            tr.stable_id.as_ref(),
        )) {
            let a = before
                .peptide
                .as_deref()
                .unwrap_or("")
                .trim_end_matches('*');
            let b = tr.peptide.as_deref().unwrap_or("").trim_end_matches('*');
            if a != b {
                let edits: Vec<_> = a
                    .bytes()
                    .zip(b.bytes())
                    .enumerate()
                    .filter(|(_, (x, y))| x != y)
                    .map(|(i, (x, y))| (i + 1, (x as char).to_string(), (y as char).to_string()))
                    .collect();
                println!(
                    "{}",
                    serde_json::json!({"transcript":tr.stable_id,"chromosome":tr.chromosome,"start":tr.start,"flags":tr.flags,"lengths":[a.len(),b.len()],"differences":edits,"sourceEdits":new.enrichment.get(&format!("{}:{}",tr.chromosome,tr.stable_id)).map(|e|&e.translation_edits)})
                );
            }
        }
    }
    Ok(())
}
