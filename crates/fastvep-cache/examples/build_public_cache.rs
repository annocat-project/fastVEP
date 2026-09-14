use anyhow::{ensure, Result};
use std::path::Path;
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    ensure!(
        args.len() == 5,
        "usage: build_public_cache GFF3 FASTA CORE_DIRECTORY MANIFEST OUTPUT"
    );
    let header = fastvep_cache::ensembl_core::build(
        Path::new(&args[0]),
        Path::new(&args[1]),
        Path::new(&args[2]),
        Path::new(&args[3]),
        Path::new(&args[4]),
    )?;
    println!("{}", serde_json::to_string_pretty(&header)?);
    Ok(())
}
