use anyhow::{bail, Context, Result};
use fastvep_cache::transcript_cache::{load_cache, save_cache, verify_cache};
use fastvep_genome::Transcript;
use serde::Deserialize;
use serde_json::json;
use sha2::{Digest, Sha256};
use std::collections::BTreeSet;
use std::env;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

#[derive(Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
struct Overlay {
    schema_version: u32,
    test_only: bool,
    input_cache_sha256: String,
    operations: Vec<Operation>,
}

#[derive(Deserialize)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase",
    deny_unknown_fields
)]
enum Operation {
    AddFlag {
        transcript_id: String,
        transcript_version: u32,
        flag: String,
    },
    SetPeptideResidue {
        transcript_id: String,
        transcript_version: u32,
        position: usize,
        expected: String,
        replacement: String,
    },
}

impl Operation {
    fn transcript(&self) -> (&str, u32) {
        match self {
            Self::AddFlag {
                transcript_id,
                transcript_version,
                ..
            }
            | Self::SetPeptideResidue {
                transcript_id,
                transcript_version,
                ..
            } => (transcript_id, *transcript_version),
        }
    }

    fn apply(&self, transcript: &mut Transcript) -> Result<()> {
        match self {
            Self::AddFlag { flag, .. } => {
                if flag.is_empty() || transcript.flags.iter().any(|value| value == flag) {
                    bail!(
                        "{} already has or was given an empty flag {:?}",
                        transcript.stable_id,
                        flag
                    );
                }
                transcript.flags.push(flag.clone());
                transcript.flags.sort();
                transcript.flags.dedup();
            }
            Self::SetPeptideResidue {
                position,
                expected,
                replacement,
                ..
            } => {
                let expected = single_ascii(expected, "expected residue")?;
                let replacement = single_ascii(replacement, "replacement residue")?;
                let index = position
                    .checked_sub(1)
                    .context("peptide position must be one-based")?;
                let peptide = transcript
                    .peptide
                    .as_mut()
                    .with_context(|| format!("{} has no peptide", transcript.stable_id))?;
                let mut bytes = peptide.as_bytes().to_vec();
                let found = bytes.get(index).copied().with_context(|| {
                    format!(
                        "{} peptide position {} is out of range",
                        transcript.stable_id, position
                    )
                })?;
                if found != expected {
                    bail!(
                        "{} peptide position {} expected {}, found {}",
                        transcript.stable_id,
                        position,
                        expected as char,
                        found as char
                    );
                }
                bytes[index] = replacement;
                *peptide = String::from_utf8(bytes).context("patched peptide is not UTF-8")?;
            }
        }
        Ok(())
    }
}

fn single_ascii(value: &str, label: &str) -> Result<u8> {
    let bytes = value.as_bytes();
    if bytes.len() != 1 || !bytes[0].is_ascii() {
        bail!("{label} must be one ASCII character");
    }
    Ok(bytes[0])
}

fn sha256(path: &Path) -> Result<String> {
    let mut file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut digest = Sha256::new();
    let mut buffer = vec![0u8; 1024 * 1024];
    loop {
        let read = file.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        digest.update(&buffer[..read]);
    }
    Ok(format!("{:x}", digest.finalize()))
}

fn parse_args() -> Result<(PathBuf, PathBuf, Vec<PathBuf>)> {
    let mut args = env::args_os().skip(1);
    let input = args.next().map(PathBuf::from).context("missing CACHE")?;
    let output = args.next().map(PathBuf::from).context("missing OUTPUT")?;
    let overlays: Vec<_> = args.map(PathBuf::from).collect();
    if overlays.is_empty() {
        bail!("usage: apply_causal_overlay CACHE OUTPUT OVERLAY [OVERLAY ...]");
    }
    if input == output || output.exists() {
        bail!("OUTPUT must be a new path distinct from CACHE");
    }
    Ok((input, output, overlays))
}

fn main() -> Result<()> {
    let (input, output, overlay_paths) = parse_args()?;
    let input_sha256 = sha256(&input)?;
    let mut overlays = Vec::with_capacity(overlay_paths.len());
    for overlay_path in &overlay_paths {
        let overlay: Overlay = serde_json::from_reader(BufReader::new(
            File::open(overlay_path)
                .with_context(|| format!("opening {}", overlay_path.display()))?,
        ))
        .with_context(|| format!("parsing {}", overlay_path.display()))?;
        if overlay.schema_version != 1 || !overlay.test_only || overlay.operations.is_empty() {
            bail!("overlay must be non-empty, schemaVersion 1, and explicitly testOnly");
        }
        if !input_sha256.eq_ignore_ascii_case(&overlay.input_cache_sha256) {
            bail!(
                "input cache SHA-256 mismatch for {}: expected {}, found {}",
                overlay_path.display(),
                overlay.input_cache_sha256,
                input_sha256
            );
        }
        overlays.push(overlay);
    }
    let verification = verify_cache(&input, false)?;
    if verification.cache_format != "FSTVEP02" {
        bail!("causal replay requires an FSTVEP02 cache");
    }

    let mut transcripts = load_cache(&input)?;
    let mut touched = BTreeSet::new();
    let mut operations_applied = 0;
    for overlay in &overlays {
        for operation in &overlay.operations {
            let (id, version) = operation.transcript();
            let matches: Vec<_> = transcripts
                .iter()
                .enumerate()
                .filter(|(_, transcript)| {
                    transcript.stable_id.as_ref() == id && transcript.version == Some(version)
                })
                .map(|(index, _)| index)
                .collect();
            if matches.len() != 1 {
                bail!(
                    "expected one {}.{} transcript, found {}",
                    id,
                    version,
                    matches.len()
                );
            }
            operation.apply(&mut transcripts[matches[0]])?;
            touched.insert(format!("{id}.{version}"));
            operations_applied += 1;
        }
    }

    save_cache(&transcripts, &output)?;
    let output_verification = verify_cache(&output, false)?;
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "schemaVersion": 1,
            "testOnly": true,
            "inputCache": input,
            "inputCacheSha256": input_sha256,
            "overlays": overlay_paths,
            "operationsApplied": operations_applied,
            "transcriptsTouched": touched,
            "outputCache": output,
            "outputCacheSha256": sha256(&output)?,
            "outputCacheBytes": output_verification.cache_bytes,
            "transcriptCount": output_verification.transcript_count,
        }))?
    );
    Ok(())
}
