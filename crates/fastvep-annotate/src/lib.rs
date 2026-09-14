//! Shared annotation engine for fastVEP.
//!
//! Provides [`AnnotationContext`] which loads transcript models, reference sequences,
//! and supplementary annotation sources, then annotates VCF variants against them.
//!
//! Used by both `fastvep-web` (production HTTP server) and `fastvep-cli` (embedded web server).
//! The CLI batch pipeline (`run_annotate`) has its own streaming implementation but shares
//! the same underlying crates.

mod hgvs_normalize;

pub use hgvs_normalize::{
    convert_ins_to_dup_range, convert_ins_to_dup_range_noncoding, hgvsc_intronic_shifted,
    intronic_dup_span, intronic_ins_as_dup, three_prime_shift_intronic,
};

use anyhow::{Context, Result};
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_cache::fasta::FastaReader;
use fastvep_cache::gff::parse_gff3;
use fastvep_cache::providers::{
    FastaSequenceProvider, IndexedTranscriptProvider, PrefetchedSequenceProvider, SequenceProvider, TranscriptProvider,
};
use fastvep_consequence::{AlleleConsequenceResult, ConsequencePredictor};
use fastvep_core::{Allele, Consequence};
use fastvep_genome::{Exon, Transcript};
use fastvep_io::output;
use fastvep_io::variant::{AlleleAnnotation, PositionRange, TranscriptVariation, VariationFeature};
use fastvep_io::vcf::VcfParser;
use rayon::prelude::*;
use std::fs::File;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

/// Pre-loaded annotation context shared by web and CLI.
///
/// Holds transcript models, a reference sequence provider, a consequence predictor,
/// and supplementary annotation providers (ClinVar, gnomAD, etc.).
pub struct AnnotationContext {
    pub transcript_provider: IndexedTranscriptProvider,
    pub seq_provider: Option<Box<dyn SequenceProvider + Send + Sync>>,
    pub predictor: ConsequencePredictor,
    pub gff3_source: Option<String>,
    pub distance: u64,
    pub hgvs: bool,
    /// Supplementary annotation providers (ClinVar, gnomAD, etc.)
    /// Wrapped in Mutex because SA readers use internal caches that need &mut.
    pub sa_providers: Vec<Mutex<Box<dyn AnnotationProvider>>>,
    /// Gene-level annotation providers (OMIM, gnomAD gene constraints, ClinVar protein index).
    pub gene_providers: Vec<fastvep_sa::gene::GeneIndex>,
    /// ACMG-AMP classification configuration (None = disabled).
    pub acmg_config: Option<fastvep_classification::AcmgConfig>,
}

impl AnnotationContext {
    /// Build a context from GFF3, optional FASTA, and optional SA directory.
    pub fn new(
        gff3: Option<&str>,
        fasta: Option<&str>,
        sa_dir: Option<&str>,
        distance: u64,
    ) -> Result<Self> {
        let gff3_source: Option<String> = gff3.map(|p| {
            Path::new(p)
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_else(|| p.to_string())
        });

        let mut transcripts = if let Some(gff3_path) = gff3 {
            let cache_path =
                fastvep_cache::transcript_cache::default_cache_path(Path::new(gff3_path));
            let from_cache = if cache_path.exists() {
                let is_fresh = fastvep_cache::transcript_cache::cache_is_fresh(
                    &cache_path,
                    Path::new(gff3_path),
                );
                if is_fresh {
                    fastvep_cache::transcript_cache::load_cache(&cache_path).ok()
                } else {
                    None
                }
            } else {
                None
            };
            if let Some(trs) = from_cache {
                tracing::info!("Loaded {} transcripts from cache", trs.len());
                trs
            } else {
                let gff_file = File::open(gff3_path)
                    .with_context(|| format!("Opening GFF3 file: {}", gff3_path))?;
                // Auto-decompress gzipped GFF3. Without this, parse_gff3
                // reads binary gz bytes as text, yields zero transcripts,
                // and downstream silently produces empty annotations.
                let trs = if gff3_path.ends_with(".gz") || gff3_path.ends_with(".bgz") {
                    parse_gff3(flate2::read::MultiGzDecoder::new(gff_file))?
                } else {
                    parse_gff3(gff_file)?
                };
                if trs.is_empty() {
                    return Err(anyhow::anyhow!(
                        "GFF3 file {} produced 0 transcripts — likely malformed, truncated, or unrecognized format. Refusing to continue with empty transcript set.",
                        gff3_path
                    ));
                }
                tracing::info!("Loaded {} transcripts from {}", trs.len(), gff3_path);
                if let Err(e) = fastvep_cache::transcript_cache::save_cache(&trs, &cache_path) {
                    tracing::warn!("Could not save cache: {}", e);
                }
                trs
            }
        } else {
            Vec::new()
        };

        let seq_provider: Option<Box<dyn SequenceProvider + Send + Sync>> =
            if let Some(fasta_path) = fasta {
                let fai_path = format!("{}.fai", fasta_path);
                if Path::new(&fai_path).exists() {
                    let reader =
                        fastvep_cache::fasta::MmapFastaReader::open(Path::new(fasta_path))?;
                    tracing::info!("Memory-mapped FASTA from {}", fasta_path);
                    Some(Box::new(
                        fastvep_cache::providers::MmapFastaSequenceProvider::new(reader),
                    ))
                } else {
                    let fasta_file = File::open(fasta_path)
                        .with_context(|| format!("Opening FASTA: {}", fasta_path))?;
                    let reader = FastaReader::from_reader(fasta_file)?;
                    tracing::info!("Loaded FASTA from {}", fasta_path);
                    Some(Box::new(FastaSequenceProvider::new(reader)))
                }
            } else {
                None
            };

        // Build sequences for coding transcripts (parallel via rayon)
        if let Some(ref sp) = seq_provider {
            let built = AtomicUsize::new(0);
            transcripts.par_iter_mut().for_each(|tr| {
                if tr.is_coding() && tr.spliced_seq.is_none() {
                    if tr
                        .build_sequences(|chrom, start, end| {
                            sp.fetch_sequence(chrom, start, end)
                                .map_err(|e| e.to_string())
                        })
                        .is_ok()
                    {
                        built.fetch_add(1, Ordering::Relaxed);
                    }
                }
            });
            let built = built.load(Ordering::Relaxed);
            tracing::info!("Built sequences for {} coding transcripts", built);

            // Re-save cache with pre-built sequences so future startups skip this step
            if built > 0 {
                if let Some(gff3_path) = gff3 {
                    let cache_path =
                        fastvep_cache::transcript_cache::default_cache_path(Path::new(gff3_path));
                    match fastvep_cache::transcript_cache::save_cache(&transcripts, &cache_path) {
                        Ok(()) => tracing::info!("Updated cache with pre-built sequences"),
                        Err(e) => tracing::warn!("Could not update cache with sequences: {}", e),
                    }
                }
            }
        }

        let transcript_provider = IndexedTranscriptProvider::new(transcripts);
        let predictor = ConsequencePredictor::new(distance, distance);

        // Load supplementary annotation providers (.osa, .osa2 files)
        let sa_providers = if let Some(dir) = sa_dir {
            load_sa_providers(Path::new(dir))?
        } else {
            Vec::new()
        };

        let gene_providers = if let Some(dir) = sa_dir {
            load_gene_providers(Path::new(dir))?
        } else {
            Vec::new()
        };

        Ok(Self {
            transcript_provider,
            seq_provider,
            predictor,
            gff3_source,
            distance,
            hgvs: true,
            sa_providers,
            gene_providers,
            acmg_config: None,
        })
    }

    pub fn transcript_count(&self) -> usize {
        self.transcript_provider.transcript_count()
    }

    /// Names of loaded supplementary annotation sources.
    pub fn sa_source_names(&self) -> Vec<String> {
        self.sa_providers
            .iter()
            .filter_map(|m| {
                let guard = m.lock().ok()?;
                Some(guard.name().to_string())
            })
            .collect()
    }

    /// Load a genome from GFF3 path (+ optional FASTA + optional SA directory).
    /// Replaces transcripts, sequence provider, and SA providers.
    pub fn load_genome(
        &mut self,
        gff3_path: &str,
        fasta_path: Option<&str>,
        sa_dir: Option<&str>,
    ) -> Result<usize> {
        let new_ctx = Self::new(Some(gff3_path), fasta_path, sa_dir, self.distance)?;
        let tr_count = new_ctx.transcript_provider.transcript_count();
        self.transcript_provider = new_ctx.transcript_provider;
        self.seq_provider = new_ctx.seq_provider;
        self.predictor = new_ctx.predictor;
        self.gff3_source = new_ctx.gff3_source;
        self.sa_providers = new_ctx.sa_providers;
        self.gene_providers = new_ctx.gene_providers;
        Ok(tr_count)
    }

    /// Replace the transcript models by parsing GFF3 text uploaded from the browser.
    pub fn update_gff3_text(&mut self, gff3_text: &str) -> Result<(usize, usize)> {
        let mut transcripts = parse_gff3(gff3_text.as_bytes())?;
        let gene_count = {
            let mut genes = std::collections::HashSet::new();
            for t in &transcripts {
                genes.insert(t.gene.stable_id.clone());
            }
            genes.len()
        };

        if let Some(ref sp) = self.seq_provider {
            let mut built = 0usize;
            for tr in &mut transcripts {
                if tr.is_coding() && tr.spliced_seq.is_none() {
                    if tr
                        .build_sequences(|chrom, start, end| {
                            sp.fetch_sequence(chrom, start, end)
                                .map_err(|e| e.to_string())
                        })
                        .is_ok()
                    {
                        built += 1;
                    }
                }
            }
            if built > 0 {
                tracing::info!("Built sequences for {} coding transcripts", built);
            }
        }

        let tr_count = transcripts.len();
        self.transcript_provider = IndexedTranscriptProvider::new(transcripts);
        self.gff3_source = Some("user-upload".to_string());
        tracing::info!(
            "Updated GFF3: {} genes, {} transcripts",
            gene_count,
            tr_count
        );
        Ok((gene_count, tr_count))
    }

    /// Annotate VCF text and return JSON results, using `self.acmg_config`.
    pub fn annotate_vcf_text(&self, vcf_text: &str, pick: bool) -> Result<Vec<serde_json::Value>> {
        self.annotate_vcf_text_with_acmg(vcf_text, pick, self.acmg_config.as_ref())
    }

    /// Annotate VCF text and return JSON results, using an explicit ACMG
    /// config for this call instead of `self.acmg_config`.
    ///
    /// This only needs `&self` (no mutation of shared context state), so
    /// callers serving concurrent requests over a shared `AnnotationContext`
    /// can take a read lock instead of a write lock — see fastvep-web's
    /// `annotate` handler, which used to mutate `self.acmg_config` per
    /// request under a write lock, serializing all concurrent traffic
    /// (including unrelated reads) behind whichever annotation was running.
    pub fn annotate_vcf_text_with_acmg(
        &self,
        vcf_text: &str,
        pick: bool,
        acmg_config: Option<&fastvep_classification::AcmgConfig>,
    ) -> Result<Vec<serde_json::Value>> {
        let mut vcf_parser = VcfParser::new(vcf_text.as_bytes())?;

        // Extract sample names from VCF #CHROM header
        let sample_names: Vec<String> = vcf_parser
            .header_lines()
            .last()
            .filter(|l| l.starts_with("#CHROM"))
            .map(|l| l.split('\t').skip(9).map(|s| s.to_string()).collect())
            .unwrap_or_default();

        let mut variants = vcf_parser.read_all()?;

        for vf in &mut variants {
            let chrom = &vf.position.chromosome;
            // VEP AnnotationType::Transcript selects overlap after parser
            // minimization, before creating transcript/allele annotations.
            let query_position = vf.alt_alleles.first().map(|alternate| {
                fastvep_consequence::vep_input_position(
                    &vf.position, &vf.ref_allele, alternate, vf.minimised || vf.alt_alleles.len() > 1,
                )
            }).unwrap_or_else(|| vf.position.clone());
            let query_start = query_position.start.saturating_sub(self.distance).max(1);
            let query_end = query_position.end.saturating_add(self.distance);
            let overlapping = self
                .transcript_provider
                .get_transcripts(chrom, query_start, query_end)
                .unwrap_or_default();

            if overlapping.is_empty() {
                annotate_intergenic(vf);
            } else {
                let ref_seq = self
                    .seq_provider
                    .as_ref()
                    .and_then(|sp| sp.fetch_sequence(chrom, query_start, query_end).ok());
                let hgvs_provider = self.seq_provider.as_deref().map(|inner| PrefetchedSequenceProvider {
                    inner, chrom, start: query_start, bases: ref_seq.as_deref(),
                });

                let result = self.predictor.predict_with_parsed_input(
                    &vf.position,
                    &vf.ref_allele,
                    &vf.alt_alleles,
                    &overlapping,
                    ref_seq.as_deref(),
                    vf.minimised,
                );

                for (i, tc) in result.transcript_consequences.iter().enumerate() {
                    let transcript: Option<&Transcript> = overlapping
                        .get(i)
                        .copied()
                        .filter(|t| t.stable_id == tc.transcript_id)
                        .or_else(|| {
                            overlapping
                                .iter()
                                .copied()
                                .find(|t| t.stable_id == tc.transcript_id)
                        });

                    let allele_annotations: Vec<AlleleAnnotation> = tc
                        .allele_consequences
                        .iter()
                        .map(|ac| {
                            let (cdna_position, cds_position, protein_position) = transcript
                                .map(|tr| vep_position_ranges(tr, ac, vf))
                                .unwrap_or_else(|| {
                                    (
                                        zip_positions(ac.cdna_start, ac.cdna_end),
                                        zip_positions(ac.cds_start, ac.cds_end),
                                        ac.protein_range()
                                            .map(|(start, end)| PositionRange::complete(start, end))
                                            .unwrap_or_default(),
                                    )
                                });
                            let mut ann = AlleleAnnotation {
                                allele: ac.allele.clone(),
                                consequences: ac.consequences.clone(),
                                impact: ac.impact,
                                cdna_position,
                                cds_position,
                                protein_position,
                                amino_acids: ac.amino_acids.clone(),
                                codons: ac.codons.clone(),
                                exon: ac.exon,
                                intron: ac.intron,
                                distance: ac.distance,
                                hgvsc: None,
                                hgvsp: None,
                                hgvsg: None,
                                hgvs_offset: None,
                                existing_variation: vec![],
                                sift: None,
                                polyphen: None,
                                supplementary: Vec::new(),
                                acmg_classification: None,
                            };

                            if let Some(change) = transcript
                                .filter(|_| ac.amino_acids.is_some() || ac.codons.is_some())
                                .and_then(|tr| {
                                    self.predictor.display_coding_change(
                                        &vf.position,
                                        &vf.ref_allele,
                                        &ac.allele,
                                        tr,
                                        vf.minimised || vf.alt_alleles.len() > 1,
                                    )
                                })
                            {
                                ann.amino_acids = change.amino_acids;
                                ann.codons = change.codons;
                            }
                            if self.hgvs {
                                ann.hgvsg = Some(fastvep_hgvs::hgvsg(
                                    chrom,
                                    vf.position.start,
                                    vf.position.end,
                                    &vf.ref_allele,
                                    &ac.allele,
                                ));
                                if let Some(tr) = transcript {
                                    let versioned_tid = match tr.version {
                                        Some(v) => format!("{}.{}", tc.transcript_id, v),
                                        None => tc.transcript_id.to_string(),
                                    };
                                    ann.hgvsc = hgvsc_for_variation_allele(
                                        hgvs_provider
                                            .as_ref()
                                            .map(|sp| sp as &dyn SequenceProvider),
                                        chrom,
                                        tr,
                                        &versioned_tid,
                                        vf,
                                        ac,
                                    );
                                    // Calculate the shift before protein HGVS so both notations
                                    // can use the same transcript-relative normalization.
                                    ann.hgvs_offset = hgvs_offset_for_allele(
                                        hgvs_provider
                                            .as_ref()
                                            .map(|sp| sp as &dyn SequenceProvider),
                                        chrom,
                                        tr,
                                        vf,
                                        ac,
                                    );

                                    if let Some(ref pid) = tr.protein_id {
                                        let versioned_pid: String = match tr.protein_version {
                                            Some(v) => {
                                                let suffix = format!(".{}", v);
                                                if pid.ends_with(suffix.as_str()) {
                                                    pid.clone()
                                                } else {
                                                    format!("{}.{}", pid, v)
                                                }
                                            }
                                            None => pid.clone(),
                                        };
                                        ann.hgvsp = hgvsp_for_variation_allele_with_offset(
                                            hgvs_provider
                                                .as_ref()
                                                .map(|sp| sp as &dyn SequenceProvider),
                                            chrom,
                                            tr,
                                            &versioned_pid,
                                            ac,
                                            vf,
                                            ann.hgvsc.as_deref(),
                                            ann.hgvs_offset,
                                        );
                                    }
                                    // HGVS_OFFSET describes a shift applied to a reported
                                    // transcript or protein HGVS notation. Do not expose an
                                    // otherwise internal shift when neither notation rendered.
                                    if ann.hgvsc.is_none() && ann.hgvsp.is_none() {
                                        ann.hgvs_offset = None;
                                    }
                                }
                            }
                            ann
                        })
                        .collect();

                    let should_include =
                        !pick || tc.canonical || vf.transcript_variations.is_empty();
                    if should_include {
                        vf.transcript_variations.push(TranscriptVariation {
                            transcript_id: tc.transcript_id.clone(),
                            gene_id: tc.gene_id.clone(),
                            gene_symbol: tc.gene_symbol.clone(),
                            biotype: tc.biotype.clone(),
                            allele_annotations,
                            canonical: tc.canonical,
                            strand: tc.strand,
                            source: self.gff3_source.clone(),
                            protein_id: transcript.and_then(|t| t.protein_id.clone()),
                            mane_select: transcript.and_then(|t| t.mane_select.clone()),
                            mane_plus_clinical: transcript
                                .and_then(|t| t.mane_plus_clinical.clone()),
                            tsl: transcript.and_then(|t| t.tsl),
                            appris: transcript.and_then(|t| t.appris.clone()),
                            ccds: transcript.and_then(|t| t.ccds.clone()),
                            gencode_primary: transcript.map(|t| t.gencode_primary).unwrap_or(false),
                            symbol_source: transcript.and_then(|t| t.gene.symbol_source.clone()),
                            hgnc_id: transcript.and_then(|t| t.gene.hgnc_id.clone()),
                            flags: transcript.map(|t| t.flags.clone()).unwrap_or_default(),
                        });
                    }
                }
            }

            // Supplementary annotation: query SA providers once per unique
            // allele, then attach the result to each (transcript, allele)
            // slot. Payload lookup is allele-level; transcript-specific vectors
            // remain intact for downstream transcript selection.
            if !self.sa_providers.is_empty() {
                let chrom = &vf.position.chromosome;
                let sa_queries = vf.supplementary_query_alleles();
                // gnomAD stores left-aligned, parsimonious alleles; normalize the
                // query to its minimal representation so indels (especially in
                // repeats) match instead of silently missing — which otherwise
                // makes PM2 misfire on common variants. Only when a reference and
                // a gnomAD provider are present; applied only to the gnomAD lookup
                // (every other source keeps the existing query unchanged).
                let has_gnomad = self.seq_provider.is_some()
                    && self
                        .sa_providers
                        .iter()
                        .any(|sa| sa.lock().unwrap().json_key() == "gnomad");
                let mut allele_results: std::collections::HashMap<String, Vec<(String, String)>> =
                    std::collections::HashMap::new();
                for tv in &vf.transcript_variations {
                    for aa in &tv.allele_annotations {
                        let alt_str = aa.allele.to_string();
                        if allele_results.contains_key(&alt_str) {
                            continue;
                        }
                        let (_, query_pos, ref_str, query_alt) = sa_queries
                            .iter()
                            .find(|(allele, _, _, _)| allele == &alt_str)
                            .expect("annotation allele must belong to the input variant");
                        let gnomad_norm = if has_gnomad {
                            self.seq_provider.as_ref().map(|sp| {
                                fastvep_cache::normalize::normalize_variant(
                                    &**sp,
                                    chrom,
                                    *query_pos,
                                    ref_str,
                                    query_alt,
                                )
                            })
                        } else {
                            None
                        };
                        let mut results: Vec<(String, String)> = Vec::new();
                        for sa in &self.sa_providers {
                            let sa_guard = sa.lock().unwrap();
                            let (q_pos, q_ref, q_alt) = if !sa_guard.metadata().match_by_allele {
                                (vf.position.start, "", "")
                            } else if sa_guard.json_key() == "gnomad" {
                                match &gnomad_norm {
                                    Some(n) => {
                                        (n.pos, n.ref_allele.as_str(), n.alt_allele.as_str())
                                    }
                                    None => (*query_pos, ref_str.as_str(), query_alt.as_str()),
                                }
                            } else {
                                (*query_pos, ref_str.as_str(), query_alt.as_str())
                            };
                            match sa_guard.annotate_position(chrom, q_pos, q_ref, q_alt) {
                                Ok(Some(ann)) => {
                                    let json_str = match ann {
                                        AnnotationValue::Json(j) => j,
                                        AnnotationValue::Positional(j) => j,
                                        AnnotationValue::Interval(v) => {
                                            format!("[{}]", v.join(","))
                                        }
                                    };
                                    results.push((sa_guard.json_key().to_string(), json_str));
                                }
                                Ok(None) => {}
                                Err(error) => {
                                    sa_lookup_errors().record(sa_guard.json_key(), &error)
                                }
                            }
                        }
                        allele_results.insert(alt_str, results);
                    }
                }
                for tv in &mut vf.transcript_variations {
                    for aa in &mut tv.allele_annotations {
                        if let Some(results) = allele_results.get(&aa.allele.to_string()) {
                            aa.supplementary.extend(results.iter().cloned());
                        }
                    }
                }
            }

            // Gene-level annotation pass (OMIM, gnomAD gene constraints, etc.)
            if !self.gene_providers.is_empty() {
                use fastvep_cache::annotation::GeneAnnotationProvider;
                let mut seen_genes = std::collections::HashSet::new();
                for tv in &vf.transcript_variations {
                    if let Some(gene_sym) = tv.gene_symbol.as_deref() {
                        if seen_genes.insert(gene_sym.to_string()) {
                            for gp in &self.gene_providers {
                                if let Ok(Some(json)) = gp.annotate_gene(gene_sym) {
                                    vf.gene_annotations.push(fastvep_core::GeneAnnotation {
                                        gene_symbol: gene_sym.to_string(),
                                        json_key: gp.json_key().to_string(),
                                        json_string: json,
                                    });
                                }
                            }
                        }
                    }
                }
            }

            // ACMG-AMP classification pass (after all SA annotations are attached)
            if let Some(acmg_cfg) = acmg_config {
                // Parse sample genotypes if trio config is present
                let trio_genotypes = extract_trio_genotypes(vf, acmg_cfg, &sample_names);

                for tv in &mut vf.transcript_variations {
                    let gene_sym = tv.gene_symbol.as_deref().unwrap_or("");
                    let gene_anns: Vec<&fastvep_core::GeneAnnotation> = vf
                        .gene_annotations
                        .iter()
                        .filter(|ga| ga.gene_symbol == gene_sym)
                        .collect();
                    for aa in &mut tv.allele_annotations {
                        let input = fastvep_classification::extract_classification_input(
                            &aa.consequences,
                            aa.impact,
                            tv.gene_symbol.as_deref(),
                            tv.canonical,
                            aa.amino_acids.as_ref(),
                            aa.protein_position.first_known(),
                            aa.hgvsc.as_deref(),
                            aa.exon.map(|(first, _, total)| (first, total)),
                            &aa.supplementary,
                            &gene_anns,
                            &vf.supplementary_annotations,
                            trio_genotypes.0.clone(),
                            trio_genotypes.1.clone(),
                            trio_genotypes.2.clone(),
                            vec![], // companion_variants populated in second pass
                        );
                        let result = fastvep_classification::classify(&input, acmg_cfg);
                        aa.acmg_classification = serde_json::to_value(&result).ok();
                    }
                }
            }

            vf.compute_most_severe();
        }

        // Compound-het enrichment pass: re-evaluate PM3/BP2 with companion variant data
        if let Some(acmg_cfg) = acmg_config {
            if acmg_cfg.trio.is_some() {
                enrich_compound_het(&mut variants, acmg_cfg, &sample_names);
            }
        }

        Ok(variants
            .iter()
            .map(|vf| output::format_json(vf, false))
            .collect())
    }
}

/// Per-allele scaffold for `--sa-only` mode: creates a TranscriptVariation
/// per alt allele with empty consequences so the SA attachment loop has a
/// slot to populate while emitting no default-CSQ annotation.
pub fn annotate_sa_only_scaffold(vf: &mut VariationFeature) {
    for alt in &vf.alt_alleles {
        vf.transcript_variations.push(TranscriptVariation {
            transcript_id: "-".into(),
            gene_id: "-".into(),
            gene_symbol: None,
            biotype: "-".into(),
            allele_annotations: vec![AlleleAnnotation {
                allele: alt.clone(),
                consequences: vec![],
                impact: fastvep_core::Impact::Modifier,
                cdna_position: PositionRange::default(),
                cds_position: PositionRange::default(),
                protein_position: PositionRange::default(),
                amino_acids: None,
                codons: None,
                exon: None,
                intron: None,
                distance: None,
                hgvsc: None,
                hgvsp: None,
                hgvsg: None,
                hgvs_offset: None,
                existing_variation: vec![],
                sift: None,
                polyphen: None,
                supplementary: Vec::new(),
                acmg_classification: None,
            }],
            canonical: false,
            strand: fastvep_core::Strand::Forward,
            source: None,
            protein_id: None,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            gencode_primary: false,
            symbol_source: None,
            hgnc_id: None,
            flags: Vec::new(),
        });
    }
}

pub fn annotate_intergenic(vf: &mut VariationFeature) {
    for alt in vf.alt_alleles.iter().filter(|alt| *alt != &vf.ref_allele) {
        vf.transcript_variations.push(TranscriptVariation {
            transcript_id: "-".into(),
            gene_id: "-".into(),
            gene_symbol: None,
            biotype: "-".into(),
            allele_annotations: vec![AlleleAnnotation {
                allele: alt.clone(),
                consequences: vec![Consequence::IntergenicVariant],
                impact: fastvep_core::Impact::Modifier,
                cdna_position: PositionRange::default(),
                cds_position: PositionRange::default(),
                protein_position: PositionRange::default(),
                amino_acids: None,
                codons: None,
                exon: None,
                intron: None,
                distance: None,
                hgvsc: None,
                hgvsp: None,
                hgvsg: None,
                hgvs_offset: None,
                existing_variation: vec![],
                sift: None,
                polyphen: None,
                supplementary: Vec::new(),
                acmg_classification: None,
            }],
            canonical: false,
            strand: fastvep_core::Strand::Forward,
            source: None,
            protein_id: None,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            gencode_primary: false,
            symbol_source: None,
            hgnc_id: None,
            flags: Vec::new(),
        });
    }
    vf.most_severe_consequence = (!vf.transcript_variations.is_empty())
        .then_some(Consequence::IntergenicVariant);
}

pub fn zip_positions(start: Option<u64>, end: Option<u64>) -> PositionRange {
    match (start, end) {
        (Some(s), Some(e)) => PositionRange::complete(s.min(e), s.max(e)),
        pair => PositionRange::new(pair.0, pair.1),
    }
}

/// Project the predictor's genomic-endpoint coordinates the way VEP exposes
/// cDNA, CDS and protein positions.
pub fn vep_position_ranges(
    transcript: &Transcript,
    allele: &AlleleConsequenceResult,
    variation: &VariationFeature,
) -> (PositionRange, PositionRange, PositionRange) {
    // Display coordinates describe VEP's parsed interval. Clipping a retained
    // multiallelic replacement into an insertion must not invoke map_insert's
    // flank adjustment on its original replacement coordinates.
    let position = fastvep_consequence::vep_input_position(
        &variation.position, &variation.ref_allele, &allele.allele,
        variation.minimised || variation.alt_alleles.len() > 1,
    );
    let start = position.start;
    let end = position.end;
    let insertion = end.checked_add(1) == Some(start);
    let mut cdna = (allele.cdna_start, allele.cdna_end);
    let insertion_between_transcript_bases = insertion
        && start >= transcript.start
        && start <= transcript.end
        && end >= transcript.start
        && end <= transcript.end;

    if insertion {
        cdna = match cdna {
            (Some(a), Some(b)) => (Some(a.min(b)), Some(a.max(b))),
            (Some(value), None) if insertion_between_transcript_bases => match transcript.strand {
                fastvep_core::Strand::Forward => (value.checked_sub(1), Some(value)),
                fastvep_core::Strand::Reverse => (Some(value), value.checked_add(1)),
            },
            (None, Some(value)) if insertion_between_transcript_bases => match transcript.strand {
                fastvep_core::Strand::Forward => (Some(value), value.checked_add(1)),
                fastvep_core::Strand::Reverse => (value.checked_sub(1), Some(value)),
            },
            (Some(value), None) | (None, Some(value)) => (Some(value), Some(value)),
            pair => pair,
        };
    } else if transcript.strand == fastvep_core::Strand::Reverse {
        cdna = (cdna.1, cdna.0);
    }

    let cds = (
        cdna.0.and_then(|value| transcript.cdna_to_cds(value)),
        cdna.1.and_then(|value| transcript.cdna_to_cds(value)),
    );
    let vep_range = |(first, last): (Option<u64>, Option<u64>)| match (first, last) {
        (Some(first), Some(last)) => PositionRange::complete(first.min(last), first.max(last)),
        pair => PositionRange::new(pair.0, pair.1),
    };

    let cds_range = if insertion_between_transcript_bases && (cds.0.is_none() || cds.1.is_none()) {
        PositionRange::default()
    } else {
        vep_range(cds)
    };
    let protein_range = PositionRange::new(
        cds_range.start().map(Transcript::cds_to_protein),
        cds_range.end().map(Transcript::cds_to_protein),
    );

    (vep_range(cdna), cds_range, protein_range)
}

pub fn intronic_or_exonic_cdna(transcript: &Transcript, genomic: u64) -> Option<(u64, i64)> {
    transcript
        .genomic_to_intronic_cdna(genomic)
        .or_else(|| transcript.genomic_to_cdna(genomic).map(|cdna| (cdna, 0)))
}

pub fn cds_and_downstream(
    transcript: &Transcript,
    spliced: &str,
    coding_start: u64,
) -> Option<Vec<u8>> {
    let start = usize::try_from(coding_start.checked_sub(1)?).ok()?;
    let bases = spliced.as_bytes();
    if start > bases.len() {
        return None;
    }
    let phase = transcript.codon_table_start_phase as usize;
    let mut cds = vec![b'N'; phase];
    cds.extend_from_slice(&bases[start..]);
    Some(cds)
}

pub fn frameshift_codon_table(_transcript: &Transcript) -> fastvep_genome::CodonTable {
    // VEP 115's protein-HGVS extension calculation calls BioPerl translate
    // without a mitochondrial codon-table override. Consequence prediction
    // still uses the transcript's mitochondrial table; this is only the HGVS
    // rendering compatibility path.
    fastvep_genome::CodonTable::standard()
}

/// Generate HGVSp through the one path shared by batch and in-process annotation.
pub fn hgvsp_for_allele(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &Transcript,
    protein_id: &str,
    allele: &AlleleConsequenceResult,
    ref_allele: &Allele,
    hgvsc: Option<&str>,
) -> Option<String> {
    hgvsp_for_allele_with_offset(
        seq_provider,
        chrom,
        transcript,
        protein_id,
        allele,
        ref_allele,
        hgvsc,
        None,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn hgvsp_for_allele_with_offset(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &Transcript,
    protein_id: &str,
    allele: &AlleleConsequenceResult,
    ref_allele: &Allele,
    hgvsc: Option<&str>,
    hgvs_offset: Option<i64>,
) -> Option<String> {
    hgvsp_with_record_class(seq_provider, chrom, transcript, protein_id, allele,
        ref_allele, hgvsc, hgvs_offset,
        *ref_allele == Allele::Deletion || allele.allele == Allele::Deletion)
}

/// Retain the whole record's class after shared-anchor trimming. An ALT `-`
/// in a mixed replacement record does not take VEP's standalone indel shift.
pub fn hgvsp_for_variation_allele_with_offset(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &Transcript,
    protein_id: &str,
    allele: &AlleleConsequenceResult,
    variation: &VariationFeature,
    hgvsc: Option<&str>,
    hgvs_offset: Option<i64>,
) -> Option<String> {
    hgvsp_with_record_class(seq_provider, chrom, transcript, protein_id, allele,
        &variation.ref_allele, hgvsc, hgvs_offset,
        vep_input_is_shiftable_indel(&variation.ref_allele, &variation.alt_alleles))
}

fn hgvsp_with_record_class(
    seq_provider: Option<&dyn SequenceProvider>, chrom: &str,
    transcript: &Transcript, protein_id: &str, allele: &AlleleConsequenceResult,
    ref_allele: &Allele, hgvsc: Option<&str>, hgvs_offset: Option<i64>,
    shiftable_indel: bool,
) -> Option<String> {
    let mut description = hgvsp_for_allele_with_offset_impl(
        seq_provider, chrom, transcript, protein_id, allele, ref_allele, hgvsc, hgvs_offset,
        shiftable_indel,
    )?;
    // VEP's final formatter adds stop-loss extension to the entire clipped
    // deletion span, not only when its first residue is the reference stop.
    // Transcript ablation suppresses displayed coding terms, but VEP still
    // evaluates stop_lost from the peptide pair in its HGVS formatter.
    let stop_lost = allele.consequences.contains(&Consequence::StopLost)
        || (allele.consequences.contains(&Consequence::TranscriptAblation)
            && allele.amino_acids.as_ref().is_some_and(|(reference, alternate)| {
                reference.contains('*') && !alternate.contains('*')
            }));
    if stop_lost
        && description.ends_with("del")
    {
        // Read the first coordinate of our own generated protein description:
        // shifting and clipping can move it beyond the predictor endpoint.
        let position = description.split(":p.").nth(1)?
            .split(|c: char| !c.is_ascii_digit()).find(|part| !part.is_empty())?
            .parse().ok()?;
        let suffix = transcript.spliced_seq.as_deref()
            .zip(transcript.cdna_coding_start)
            .and_then(|(spliced, start)| cds_and_downstream(transcript, spliced, start))
            .and_then(|cds| fastvep_hgvs::hgvsp_stop_lost_suffix_from_cds(
                position, transcript.peptide.as_deref()?, &cds,
                hgvs_shifted_cds_positions(transcript, allele, hgvs_offset).0,
                hgvs_shifted_cds_positions(transcript, allele, hgvs_offset).1,
                ref_allele, &allele.allele, transcript.strand,
                &frameshift_codon_table(transcript),
            ))
            .unwrap_or_else(|| "extTer?".into());
        description.push_str(&suffix);
    }
    Some(description)
}

// VEP BaseTranscriptVariation::cds_coords applies the offset in genomic
// space before exon mapping. Genomic and CDS distances differ across introns.
fn hgvs_shifted_cds_positions(
    transcript: &Transcript,
    allele: &AlleleConsequenceResult,
    genomic_offset: Option<i64>,
) -> (Option<u64>, Option<u64>) {
    let offset = genomic_offset.unwrap_or(0);
    if offset == 0 { return (allele.cds_start, allele.cds_end); }
    let map = |position: u64| {
        let cdna = transcript.genomic_to_cdna(position.checked_add_signed(offset)?)?;
        transcript.cdna_to_cds(cdna)
    };
    (map(allele.normalized_position.start), map(allele.normalized_position.end))
}

fn hgvsp_for_allele_with_offset_impl(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &Transcript,
    protein_id: &str,
    allele: &AlleleConsequenceResult,
    ref_allele: &Allele,
    hgvsc: Option<&str>,
    hgvs_offset: Option<i64>,
    shiftable_indel: bool,
) -> Option<String> {
    // TVA::hgvs_protein rejects ambiguous ALT and requires a defined reference
    // peptide; TVA::peptide cannot produce one from an ambiguous REF allele.
    if !allele.allele.is_unambiguous_dna() || !ref_allele.is_unambiguous_dna() {
        return None;
    }
    // Retained replacements do not take VEP's standalone indel shift. Its
    // hgvs_protein guard requires both original translation endpoints; an
    // intronic gap must not be revived by the clipped-deletion fallback.
    // A parsed insertion can be non-shiftable when ALT also contains REF;
    // Mapper::map_insert still removes its intronic gap before protein mapping.
    if !shiftable_indel
        && *ref_allele != Allele::Deletion
        && (allele.protein_start.is_none() || allele.protein_end.is_none())
    {
        return None;
    }
    if shiftable_indel
        && hgvs_normalize::vep_hgvs_shift_exceeds_transcript(
        seq_provider, chrom, transcript,
        allele.normalized_position.start, allele.normalized_position.end,
        &allele.normalized_ref_allele, &allele.normalized_alt_allele,
    ) {
        return None;
    }
    // VEP hgvs_protein also requires the unshifted `coding` pre-predicate.
    // Short introns stretch VEP's coding/exonic preclassification by 12 bases.
    // Such an insertion may lack original CDS coordinates and still gain HGVSp.
    if allele.normalized_ref_allele == Allele::Deletion
        && allele.cds_start.is_none() && allele.cds_end.is_none()
        && !vep_frameshift_intron_stretched_coding_overlap(
            transcript, allele.normalized_position.start, allele.normalized_position.end,
        )
    {
        return None;
    }
    let one_sided_reverse_insertion = allele.normalized_ref_allele == Allele::Deletion
        && transcript.strand == fastvep_core::Strand::Reverse
        && allele.protein_start.is_none()
        && allele.protein_end.is_some();
    // VEP 115.2 `translation_coords` retains its boundary-shift fallback when
    // the first mapped segment is a gap, but `translation_end` is still absent
    // when the last segment is a gap. The already-normalized HGVSc carries that
    // distinction: its transcript-order final endpoint must map to the CDS.
    // Ensembl Mapper::map_insert removes the intronic gap at an exon edge.
    // An unshifted insertion can therefore have mapped CDS endpoints even
    // when its HGVS flanking-base notation ends with an intronic coordinate.
    let mapped_unshifted_insertion = allele.normalized_ref_allele == Allele::Deletion
        && matches!(&allele.normalized_alt_allele, Allele::Sequence(bases) if !bases.is_empty())
        && (allele.cds_start.is_some() || allele.cds_end.is_some())
        && hgvs_offset.unwrap_or(0) == 0;
    // Mapper::map_insert removes an intronic gap at an exon edge after the
    // HGVS shift too. Intronic flanking-base HGVSc does not imply an absent
    // protein endpoint when the other flank maps into a coding exon.
    let shifted_insertion_abuts_cds = allele.normalized_ref_allele == Allele::Deletion
        && allele.normalized_position.start.checked_add_signed(hgvs_offset.unwrap_or(0))
            .zip(allele.normalized_position.end.checked_add_signed(hgvs_offset.unwrap_or(0)))
            .is_some_and(|(start, end)| {
                if start < transcript.start || start > transcript.end
                    || end < transcript.start || end > transcript.end { return false; }
                match (transcript.genomic_to_cdna(start), transcript.genomic_to_cdna(end)) {
                    (Some(cdna), None) | (None, Some(cdna)) => transcript.cdna_to_cds(cdna).is_some(),
                    _ => false,
                }
            });
    // A retained replacement keeps its original mapped protein interval even
    // when HGVS clipping moves a flank into the intron. Both original endpoints
    // were required above; the clipped HGVSc must not veto that interval.
    let has_vep_protein_endpoint =
        !shiftable_indel
        || hgvsc.is_some_and(hgvsc_has_vep_protein_endpoint) || mapped_unshifted_insertion || shifted_insertion_abuts_cds;
    if !has_vep_protein_endpoint && !one_sided_reverse_insertion {
        // VEP can still rebuild a deletion from a shifted exonic interval when
        // the original boundary-spanning edit has no reportable HGVSc. It does
        // so only on the peptide-reconstruction path represented by a missing
        // predictor amino-acid window.
        if allele.amino_acids.is_none() {
            return hgvsp_splice_boundary_deletion(
                seq_provider,
                chrom,
                transcript,
                protein_id,
                allele,
            );
        }
        return None;
    }

    // VEP shift_feature_seqs rebuilds insertion peptides after the HGVS shift,
    // including insertions that move from the UTR into the CDS. Consequences
    // and their displayed coordinates still describe the original allele.
    // Mixed multiallelic records retain replacement coordinates and must use
    // the full parsed REF/ALT in the general reconstruction below.
    if *ref_allele == Allele::Deletion
        && matches!(&allele.normalized_alt_allele, Allele::Sequence(bases) if !bases.is_empty())
        && (allele.frameshift
            || (allele.amino_acids.is_none() && has_vep_protein_endpoint))
    {
        let shift = match transcript.strand {
            fastvep_core::Strand::Forward => hgvs_offset.unwrap_or(0),
            fastvep_core::Strand::Reverse => hgvs_offset.unwrap_or(0).checked_neg()?,
        };
        let (start, end) = hgvs_shifted_cds_positions(transcript, allele, hgvs_offset);
        let insertion_point = fastvep_hgvs::cds_insertion_point(start, end, transcript.strand);
        if allele.amino_acids.is_none() && (shift <= 0 || insertion_point.is_none()) {
            return None;
        }
        if shift > 0 && !hgvs_insertion_end_is_coding(transcript, start, end) {
            return None;
        }
        let reference_table = if fastvep_genome::is_mitochondrial(&transcript.chromosome) {
            fastvep_genome::mitochondrial_codon_table()
        } else {
            fastvep_genome::CodonTable::standard()
        };
        let cds = cds_and_downstream(transcript, transcript.spliced_seq.as_deref()?, transcript.cdna_coding_start?)?;
        // The same shifted CDS reconstruction also supplies an in-frame
        // insertion's peptide when its original interval straddles the UTR.
        if allele.normalized_alt_allele.len().is_multiple_of(3) {
            // Mapper::map_insert drops an intronic flank at an exon edge.
            // Reconstruct the zero-length CDS interval from the surviving flank.
            let point = insertion_point?;
            return fastvep_hgvs::hgvsp_inframe_insertion_from_cds_with_start_lost(
                protein_id, &cds, point, point.checked_add(1)?, &allele.normalized_alt_allele,
                transcript.strand, shift.try_into().ok()?, &reference_table,
                transcript.peptide.as_deref().map(str::as_bytes), false,
                transcript.translateable_seq.as_ref().map(String::len), transcript.reference_peptide.as_deref().map(str::as_bytes),
);
        }
        let mut alternate = allele.normalized_alt_allele.clone();
        if let Allele::Sequence(bases) = &mut alternate {
            let rotation = shift.rem_euclid(bases.len() as i64) as usize;
            match transcript.strand {
                fastvep_core::Strand::Forward => bases.rotate_left(rotation),
                // VEP shift_feature_seqs loops seq_length - shift_length
                // times on reverse transcripts. A negative count does not
                // rotate, rather than wrapping modulo the allele length.
                fastvep_core::Strand::Reverse if shift <= bases.len() as i64 => bases.rotate_right(rotation),
                fastvep_core::Strand::Reverse => {},
            }
        }
        return fastvep_hgvs::hgvsp_frameshift_from_cds_with_context(
            protein_id, &cds, start, end, &Allele::Deletion, &alternate,
            transcript.strand, &reference_table, &frameshift_codon_table(transcript),
            transcript.peptide.as_deref().map(str::as_bytes),
            allele.consequences.contains(&Consequence::StopLost),
            allele.consequences.contains(&Consequence::StartLost),
            transcript.translateable_seq.as_ref().map(String::len),
        );
    }

    let Some(amino_acids) = allele.amino_acids.as_ref() else {
        return hgvsp_splice_boundary_deletion(seq_provider, chrom, transcript, protein_id, allele);
    };
    let (amino_acids, mut protein_start) = (amino_acids, allele.protein_start.or(allele.protein_end)?);
    let mut protein_end = allele.protein_end.unwrap_or(protein_start);
    if *ref_allele == Allele::Deletion && amino_acids.0 == "-" {
        if let Some(point) = fastvep_hgvs::cds_insertion_point(allele.cds_start, allele.cds_end, transcript.strand) {
            // A surviving exon endpoint represents a zero-length insertion,
            // not replacement of the one mapped residue (F1 ASCC2).
            protein_start = point / 3 + 1;
            protein_end = point.div_ceil(3);
        }
    }
    // For a replacement, VEP's hgvs_protein starts
    // from TranscriptVariation::translation_start, the lower peptide coordinate
    // returned by genomic2pep. Predictor endpoints retain genomic order, so the
    // reverse-strand pair must be sorted for this class. Pure indels
    // keep their inverted pair because the formatter uses it to locate the
    // deletion span or insertion flanks.
    let (protein_start, protein_end) = if (amino_acids.0.len() == amino_acids.1.len() && amino_acids.0.len() > 1)
        || (*ref_allele != Allele::Deletion && allele.allele != Allele::Deletion) {
        (protein_start.min(protein_end), protein_start.max(protein_end))
    } else {
        (protein_start, protein_end)
    };


    let mut is_frameshift = allele.frameshift;
    let pure_deletion = allele.normalized_alt_allele == Allele::Deletion
        && matches!(&allele.normalized_ref_allele, Allele::Sequence(bases) if !bases.is_empty());
    let pure_insertion = allele.normalized_ref_allele == Allele::Deletion
        && matches!(&allele.normalized_alt_allele, Allele::Sequence(bases) if !bases.is_empty());
    let transcript_shift = hgvs_offset.and_then(|offset| {
        let offset = match transcript.strand {
            fastvep_core::Strand::Forward => offset,
            fastvep_core::Strand::Reverse => offset.checked_neg()?,
        };
        u64::try_from(offset).ok()
    });
    // VEP 115.2 recalculates both translation endpoints after applying the
    // HGVS shift and returns no HGVSp unless both still map into the CDS. A
    // terminal insertion can shift its far flank beyond the translated region
    // even though the unshifted consequence window still contains a peptide.
    if pure_insertion && transcript_shift.is_some_and(|shift| shift > 0) {
        let (start, end) = hgvs_shifted_cds_positions(transcript, allele, hgvs_offset);
        if !hgvs_insertion_end_is_coding(transcript, start, end) { return None; }
    }
    if pure_deletion && transcript_shift.is_some_and(|shift| shift > 0) {
        // VEP hgvs_protein requires both shifted translation endpoints.
        // Do not fall back to the unshifted peptide if an endpoint is intronic.
        let (start, end) = hgvs_shifted_cds_positions(transcript, allele, hgvs_offset);
        let (Some(start), Some(end)) = (start, end) else { return None; };
        // VariationEffect::frameshift reads the shifted CDS span during HGVS,
        // while partial_codon and stop_retained retain their predicate cache.
        is_frameshift = !(start.abs_diff(end) + 1).is_multiple_of(3)
            && !allele.consequences.contains(&Consequence::IncompleteTerminalCodonVariant)
            && !allele.consequences.contains(&Consequence::StopRetainedVariant);
    }

    if is_frameshift
        && allele.consequences.contains(&Consequence::StopLost)
        && amino_acids.0 == "*"
        && allele.normalized_alt_allele == Allele::Deletion
        && fastvep_genome::is_mitochondrial(&transcript.chromosome)
    {
        return Some(format!("{}:p.Ter{}delextTer?", protein_id, protein_start));
    }

    if is_frameshift {
        let deletion_shift = if pure_deletion && ref_allele == &allele.normalized_ref_allele {
            transcript_shift.unwrap_or(0)
        } else {
            0
        };
        let reference_codon_table = if fastvep_genome::is_mitochondrial(&transcript.chromosome) {
            fastvep_genome::mitochondrial_codon_table()
        } else {
            fastvep_genome::CodonTable::standard()
        };
        let alternate_codon_table = frameshift_codon_table(transcript);
        return transcript
            .spliced_seq
            .as_deref()
            .zip(transcript.cdna_coding_start)
            .and_then(|(spliced, coding_start)| {
                cds_and_downstream(transcript, spliced, coding_start)
            })
            .and_then(|cds| {
                fastvep_hgvs::hgvsp_frameshift_from_cds_with_context(
                    protein_id,
                    &cds,
                    hgvs_shifted_cds_positions(transcript, allele, (deletion_shift > 0).then_some(hgvs_offset.unwrap_or(0))).0,
                    hgvs_shifted_cds_positions(transcript, allele, (deletion_shift > 0).then_some(hgvs_offset.unwrap_or(0))).1,
                    ref_allele,
                    &allele.allele,
                    transcript.strand,
                    &reference_codon_table,
                    &alternate_codon_table,
                    transcript.peptide.as_deref().map(str::as_bytes),
                    allele.consequences.contains(&Consequence::StopLost),
                    allele.consequences.contains(&Consequence::StartLost),
            transcript.translateable_seq.as_ref().map(String::len),
                )
            });
    }

    // HGVS shifts and rebuilds the peptide independently of consequence
    // labels, including start_retained and stop_gained whole-codon insertions.
    let shifted_inframe = (pure_deletion && !is_frameshift)
        || (pure_insertion
            && (allele.normalized_alt_allele.len().is_multiple_of(3)
                || (allele.consequences.contains(&Consequence::InframeInsertion)
                    && !allele.consequences.contains(&Consequence::StopRetainedVariant))));

    if shifted_inframe {
        if let Some(exact) = transcript_shift
            .filter(|shift| *shift > 0)
            .and_then(|shift| {
                let (start, end) = hgvs_shifted_cds_positions(transcript, allele, hgvs_offset);
                let (shifted_start, shifted_end) = if pure_insertion {
                    let point = fastvep_hgvs::cds_insertion_point(start, end, transcript.strand)?;
                    (point, point.checked_add(1)?)
                } else { (start?, end?) };
                let table = if fastvep_genome::is_mitochondrial(&transcript.chromosome) {
                    fastvep_genome::mitochondrial_codon_table()
                } else {
                    fastvep_genome::CodonTable::standard()
                };
                transcript
                    .spliced_seq
                    .as_deref()
                    .zip(transcript.cdna_coding_start)
                    .and_then(|(spliced, coding_start)| {
                        cds_and_downstream(transcript, spliced, coding_start)
                    })
                    .map(|cds| {
                        if pure_deletion {
                            fastvep_hgvs::hgvsp_inframe_deletion_from_cds(
                                protein_id,
                                &cds,
                                shifted_start,
                                shifted_end,
                                &allele.normalized_ref_allele,
                                transcript.strand,
                                &table,
                                transcript.peptide.as_deref().map(str::as_bytes),
                                allele.consequences.contains(&Consequence::StartLost),
                                transcript.translateable_seq.as_ref().map(String::len), transcript.reference_peptide.as_deref().map(str::as_bytes),
)
                        } else {
                            fastvep_hgvs::hgvsp_inframe_insertion_from_cds_with_start_lost(
                                protein_id,
                                &cds,
                                shifted_start,
                                shifted_end,
                                &allele.normalized_alt_allele,
                                transcript.strand,
                                shift,
                                &table,
                                transcript.peptide.as_deref().map(str::as_bytes),
                                allele.consequences.contains(&Consequence::StartLost),
                                transcript.translateable_seq.as_ref().map(String::len), transcript.reference_peptide.as_deref().map(str::as_bytes),
)
                        }
                    })
            })
        {
            if exact.is_some() || pure_insertion {
                return exact;
            }
        }
    }

    // VEP applies start_lost after shifting and rebuilding the peptide window.
    if allele.consequences.contains(&Consequence::StartLost)
        && !is_frameshift
    {
        // VEP hgvs_protein clips from translation_start, in peptide order.
        // Reverse-strand predictor endpoints still retain genomic order.
        let protein_start = protein_start.min(protein_end);
        // Duplication uses ordinary CDS translation, while insertion flanks
        // use the full protein. An initiator-normalized Met is not a raw Leu.
        let duplication_peptide = transcript.translateable_seq.as_deref()
            .filter(|_| amino_acids.0 == "-" || amino_acids.1.len() > amino_acids.0.len())
            .map(|cds| fastvep_genome::CodonTable::standard().translate_seq(cds.as_bytes()));
        return fastvep_hgvs::hgvsp_inframe_indel_with_context(
            protein_id,
            protein_start,
            protein_start + amino_acids.0.len().saturating_sub(1) as u64,
            &amino_acids.0,
            &amino_acids.1,
            transcript.peptide.as_deref().map(str::as_bytes),
            fastvep_core::Strand::Forward,
            true,
            duplication_peptide.as_deref(),
            // Clipping to an insertion asks VEP for full-protein flanks
            // before it formats start_lost; those can have a normalized Met.
            transcript.reference_peptide.as_deref().map(str::as_bytes),
        )
        .or_else(|| {
            amino_acids.0.as_bytes().first().map(|&reference| {
                format!(
                    "{}:p.{}{}?",
                    protein_id,
                    fastvep_genome::codon::aa_one_to_three(reference),
                    protein_start
                )
            })
        });
    }


    let shifted_stop_retained_insertion = pure_insertion
        && allele.consequences.contains(&Consequence::InframeInsertion)
        && allele
            .consequences
            .contains(&Consequence::StopRetainedVariant);
    if shifted_stop_retained_insertion {
        let exact = transcript_shift
            .filter(|shift| *shift > 0)
            .and_then(|shift| {
                transcript
                    .spliced_seq
                    .as_deref()
                    .zip(transcript.cdna_coding_start)
                    .and_then(|(spliced, coding_start)| {
                        cds_and_downstream(transcript, spliced, coding_start)
                    })
                    .and_then(|cds| {
                        let table = if fastvep_genome::is_mitochondrial(&transcript.chromosome) {
                            fastvep_genome::mitochondrial_codon_table()
                        } else {
                            fastvep_genome::CodonTable::standard()
                        };
                        fastvep_hgvs::hgvsp_shifted_stop_retained_insertion(
                            protein_id,
                            &cds,
                            allele.cds_start,
                            allele.cds_end,
                            &allele.normalized_alt_allele,
                            transcript.strand,
                            shift,
                            &table,
                            transcript.peptide.as_deref().map(str::as_bytes),
                            transcript.translateable_seq.as_ref().map(String::len),
                            Some(hgvs_shifted_cds_positions(transcript, allele, hgvs_offset)), transcript.reference_peptide.as_deref().map(str::as_bytes),
)
                    })
            });
        // This helper handles shifted partial-codon insertions. Its None
        // also means VEP deliberately has no flanking peptide; do not revive
        // the unshifted annotation after that terminal suppression.
        if exact.is_some()
            || (transcript_shift.is_some_and(|shift| shift > 0)
                && allele.normalized_alt_allele.len() % 3 != 0)
        {
            return exact;
        }

        // _get_hgvs_peptides checks duplication before requiring flanks.
        // Let the shared formatter do both; a terminal duplicate needs no
        // following residue (F1 STYK1 Leu422dup).

    }

    if amino_acids.0 == amino_acids.1 && matches!(amino_acids.0.as_str(), "X" | "*") {
        return fastvep_hgvs::hgvsp(protein_id, protein_start, b'*', b'*', false);
    }

    // VEP _clip_alleles returns the untouched window when its common prefix
    // reaches an internal stop. A one-residue shortcut must not bypass that.
    let unchanged_prefix_has_stop = amino_acids.0.bytes().zip(amino_acids.1.bytes())
        .take_while(|(reference, alternate)| reference == alternate)
        .any(|(reference, _)| reference == b'*');
    let mut changed_residues = amino_acids
        .0
        .bytes()
        .zip(amino_acids.1.bytes())
        .enumerate()
        .filter(|(_, (reference, alternate))| reference != alternate);
    let single_stop_change = changed_residues
        .next()
        .filter(|(_, (reference, _))| *reference == b'*')
        .filter(|_| {
            amino_acids.0.len() == amino_acids.1.len() && changed_residues.next().is_none()
        });
    if let Some((offset, (_, alt_aa))) =
        single_stop_change.filter(|_| !unchanged_prefix_has_stop
            && allele.consequences.contains(&Consequence::StopLost))
    {
        let stop_position = protein_start.min(protein_end) + offset as u64;
        // VEP converts X and * to Ter before its final synonymous check.
        if alt_aa == b'X' {
            return fastvep_hgvs::hgvsp(protein_id, stop_position, b'*', b'*', false);
        }

        let exact = transcript
            .spliced_seq
            .as_deref()
            .zip(transcript.cdna_coding_start)
            .and_then(|(spliced, coding_start)| {
                cds_and_downstream(transcript, spliced, coding_start)
            })
            .and_then(|cds| {
                fastvep_hgvs::hgvsp_stop_lost_from_cds(
                    protein_id,
                    stop_position,
                    transcript.peptide.as_deref()?,
                    alt_aa,
                    &cds,
                    allele.cds_start,
                    allele.cds_end,
                    ref_allele,
                    &allele.allele,
                    transcript.strand,
                    &frameshift_codon_table(transcript),
                )
            });
        return exact
            .or_else(|| fastvep_hgvs::hgvsp(protein_id, stop_position, b'*', alt_aa, false));
    }

    if amino_acids.1 == "-"
        || amino_acids.0.len() != amino_acids.1.len()
        || allele.consequences.contains(&Consequence::StartLost)
        || allele.consequences.contains(&Consequence::InframeDeletion)
        || allele.consequences.contains(&Consequence::InframeInsertion)
        || amino_acids.0.len() > 1
    {
        if amino_acids.0.len() == amino_acids.1.len() && amino_acids.0.len() > 1
            && !unchanged_prefix_has_stop {
            let mut differences = amino_acids
                .0
                .bytes()
                .zip(amino_acids.1.bytes())
                .enumerate()
                .filter(|(_, (reference, alternate))| reference != alternate);
            if let Some((offset, (reference, alternate))) = differences.next() {
                if differences.next().is_none() {
                    return fastvep_hgvs::hgvsp(
                        protein_id,
                        protein_start.min(protein_end) + offset as u64,
                        reference,
                        alternate,
                        false,
                    );
                }
            }
        }
        let reference = amino_acids.0.clone();
        let alternate = amino_acids.1.clone();
        let duplication_peptide = transcript.translateable_seq.as_deref()
            .filter(|_| reference == "-" || alternate.len() > reference.len())
            .map(|cds| fastvep_genome::CodonTable::standard().translate_seq(cds.as_bytes()));
        fastvep_hgvs::hgvsp_inframe_indel_with_context(
            protein_id,
            protein_start,
            protein_end,
            &reference,
            &alternate,
            transcript.peptide.as_deref().map(str::as_bytes),
            transcript.strand,
            false,
            duplication_peptide.as_deref(), transcript.reference_peptide.as_deref().map(str::as_bytes),
)
    } else {
        fastvep_hgvs::hgvsp(
            protein_id,
            protein_start,
            amino_acids.0.as_bytes().first().copied().unwrap_or(b'X'),
            amino_acids.1.as_bytes().first().copied().unwrap_or(b'X'),
            false,
        )
    }
}

fn hgvs_insertion_end_is_coding(
    transcript: &Transcript, start: Option<u64>, end: Option<u64>,
) -> bool {
    let point = fastvep_hgvs::cds_insertion_point(start, end, transcript.strand);
    let length = transcript.cdna_coding_end.zip(transcript.cdna_coding_start)
        .and_then(|(end, start)| end.checked_sub(start))
        .and_then(|n| n.checked_add(1 + transcript.codon_table_start_phase));
    // An absent UTR flank must not hide an insertion's far translation end.
    shifted_translation_end_is_coding(length, point, point.and_then(|p| p.checked_add(1)), 0)
}

fn shifted_translation_end_is_coding(
    cds_length: Option<u64>,
    cds_start: Option<u64>,
    cds_end: Option<u64>,
    transcript_shift: u64,
) -> bool {
    let Some(cds_length) = cds_length else {
        return false;
    };
    cds_start
        .into_iter()
        .chain(cds_end)
        .max()
        .and_then(|end| end.checked_add(transcript_shift))
        .is_some_and(|shifted_end| shifted_end <= cds_length)
}

fn hgvsc_has_vep_protein_endpoint(hgvsc: &str) -> bool {
    let Some(coding) = hgvsc.split_once(":c.").map(|(_, coding)| coding) else {
        return false;
    };
    let coordinate_end = coding
        .find(|character: char| character.is_ascii_alphabetic() || character == '=')
        .unwrap_or(coding.len());
    coding[..coordinate_end]
        .split('_')
        .next_back()
        .is_some_and(|position| {
            !position.is_empty() && position.bytes().all(|byte| byte.is_ascii_digit())
        })
}

fn hgvsp_splice_boundary_deletion(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &Transcript,
    protein_id: &str,
    allele: &AlleleConsequenceResult,
) -> Option<String> {
    if (!allele
        .consequences
        .contains(&Consequence::CodingSequenceVariant)
        && !allele
            .consequences
            .contains(&Consequence::StartRetainedVariant)
        && !allele.consequences.contains(&Consequence::StartLost)
        && !vep_frameshift_intron_stretched_coding_overlap(
            transcript,
            allele.normalized_position.start,
            allele.normalized_position.end,
        ))
        || !matches!(
            (&allele.normalized_ref_allele, &allele.normalized_alt_allele),
            (Allele::Sequence(reference), Allele::Deletion) if !reference.is_empty()
        )
    {
        return None;
    }

    let (cdna_start, cdna_end) = hgvs_normalize::exonic_deletion_cdna_span(
        seq_provider,
        chrom,
        transcript,
        allele.normalized_position.start,
        allele.normalized_position.end,
        &allele.normalized_ref_allele,
        &allele.normalized_alt_allele,
    )?;
    let (cds_start, cds_end) = (
        transcript.cdna_to_cds(cdna_start)?,
        transcript.cdna_to_cds(cdna_end)?,
    );
    let (cds_start, cds_end) = (cds_start.min(cds_end), cds_start.max(cds_end));
    if (cds_end - cds_start + 1) % 3 != 0 {
        let reference_codon_table = if fastvep_genome::is_mitochondrial(&transcript.chromosome) {
            fastvep_genome::mitochondrial_codon_table()
        } else {
            fastvep_genome::CodonTable::standard()
        };
        return transcript
            .spliced_seq
            .as_deref()
            .zip(transcript.cdna_coding_start)
            .and_then(|(spliced, coding_start)| {
                cds_and_downstream(transcript, spliced, coding_start)
            })
            .and_then(|cds| {
                fastvep_hgvs::hgvsp_frameshift_from_cds_with_context(
                    protein_id,
                    &cds,
                    Some(cds_start),
                    Some(cds_end),
                    &allele.normalized_ref_allele,
                    &allele.normalized_alt_allele,
                    transcript.strand,
                    &reference_codon_table,
                    &frameshift_codon_table(transcript),
                    transcript.peptide.as_deref().map(str::as_bytes),
                    allele.consequences.contains(&Consequence::StopLost),
                    allele.consequences.contains(&Consequence::StartLost),
            transcript.translateable_seq.as_ref().map(String::len),
                )
            });
    }
    let table = if fastvep_genome::is_mitochondrial(&transcript.chromosome) {
        fastvep_genome::mitochondrial_codon_table()
    } else {
        fastvep_genome::CodonTable::standard()
    };
    transcript
        .spliced_seq
        .as_deref()
        .zip(transcript.cdna_coding_start)
        .and_then(|(spliced, coding_start)| cds_and_downstream(transcript, spliced, coding_start))
        .and_then(|cds| {
            fastvep_hgvs::hgvsp_inframe_deletion_from_cds(
                protein_id,
                &cds,
                cds_start,
                cds_end,
                &allele.normalized_ref_allele,
                transcript.strand,
                &table,
                transcript.peptide.as_deref().map(str::as_bytes),
                allele.consequences.contains(&Consequence::StartLost),
                transcript.translateable_seq.as_ref().map(String::len), transcript.reference_peptide.as_deref().map(str::as_bytes),
)
        })
}

fn vep_frameshift_intron_stretched_coding_overlap(
    transcript: &Transcript,
    start: u64,
    end: u64,
) -> bool {
    // VEP _bvfo_preds sorts both insertion flanks before its exon query.
    let (start, end) = (start.min(end), start.max(end));
    let Some((coding_start, coding_end)) = transcript
        .coding_region_start
        .zip(transcript.coding_region_end)
    else {
        return false;
    };
    if end < coding_start || start > coding_end || !vep_has_frameshift_intron(&transcript.exons) {
        return false;
    }

    // VEP 115.2 marks a transcript with any very short intron and then expands
    // every exon by 12 bases for its exonic/coding predicate. Its source tests
    // `abs(intron_end - intron_start) <= 12`, which includes a 13-base intron.
    transcript
        .exons
        .iter()
        .any(|exon| exon.start.saturating_sub(12) <= end && start <= exon.end.saturating_add(12))
}

fn vep_has_frameshift_intron(exons: &[Exon]) -> bool {
    let mut exons: Vec<&Exon> = exons.iter().collect();
    exons.sort_unstable_by_key(|exon| exon.start);
    exons.windows(2).any(|pair| {
        pair[1]
            .start
            .checked_sub(pair[0].end)
            .is_some_and(|delta| (2..=14).contains(&delta))
    })
}

/// Supplementary-source lookup failures collected without turning them into
/// false annotation misses.
#[derive(Default)]
pub struct SaLookupErrors {
    inner: std::sync::Mutex<std::collections::HashMap<String, (u64, String)>>,
}

impl SaLookupErrors {
    pub fn record(&self, source: &str, error: &anyhow::Error) {
        let Ok(mut errors) = self.inner.lock() else {
            return;
        };
        if let Some((count, _)) = errors.get_mut(source) {
            *count += 1;
        } else {
            errors.insert(source.to_string(), (1, error.to_string()));
        }
    }

    pub fn report_lines(&self) -> Vec<String> {
        let Ok(errors) = self.inner.lock() else {
            return Vec::new();
        };
        let mut rows: Vec<_> = errors.iter().collect();
        rows.sort_by_key(|(source, _)| *source);
        rows.into_iter()
            .map(|(source, (count, first))| {
                format!(
                    "warning: {source} could not be read for {count} variant lookup(s); those variants carry no {source} annotation. First error: {first}"
                )
            })
            .collect()
    }
}

pub fn sa_lookup_errors() -> &'static SaLookupErrors {
    static ERRORS: std::sync::OnceLock<SaLookupErrors> = std::sync::OnceLock::new();
    ERRORS.get_or_init(SaLookupErrors::default)
}

pub fn report_sa_lookup_errors() {
    for line in sa_lookup_errors().report_lines() {
        eprintln!("{line}");
    }
}

pub fn complement_allele(allele: &Allele) -> Allele {
    match allele {
        Allele::Sequence(bases) => {
            let comp: Vec<u8> = bases
                .iter()
                .rev()
                .map(|&b| match b {
                    b'A' | b'a' => b'T',
                    b'T' | b't' => b'A',
                    b'C' | b'c' => b'G',
                    b'G' | b'g' => b'C',
                    other => other,
                })
                .collect();
            Allele::Sequence(comp)
        }
        other => other.clone(),
    }
}

/// Generate transcript HGVS from the allele-local minimal representation.
///
/// Consequence prediction intentionally uses the uploaded VCF allele, but HGVS
/// must use the minimal representation of each ALT independently. This matters
/// for multiallelic records where the parser cannot trim a shared anchor.
pub fn hgvsc_for_allele(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &Transcript,
    versioned_tid: &str,
    allele: &AlleleConsequenceResult,
) -> Option<String> {
    hgvsc_for_allele_impl(seq_provider, chrom, transcript, versioned_tid, None, allele)
}

/// Generate transcript HGVS while retaining VEP's record-level variant class
/// and full parsed bounds. The compatibility wrapper above remains available
/// to callers that only have an allele consequence.
pub fn hgvsc_for_variation_allele(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &Transcript,
    versioned_tid: &str,
    variation: &VariationFeature,
    allele: &AlleleConsequenceResult,
) -> Option<String> {
    hgvsc_for_allele_impl(
        seq_provider,
        chrom,
        transcript,
        versioned_tid,
        Some(variation),
        allele,
    )
}

fn vep_input_is_snp(reference: &Allele, alternates: &[Allele]) -> bool {
    matches!(reference, Allele::Sequence(bases) if bases.len() == 1)
        && !alternates.is_empty()
        && alternates
            .iter()
            .all(|alternate| matches!(alternate, Allele::Sequence(bases) if bases.len() == 1))
}

// VEP hgvs_variant_notation checks a full-REF duplication before _clip_alleles
// trims the transcript-oriented prefix, then suffix. A retained replacement
// never takes _genomic_shift, even when clipping leaves an insertion/deletion.
fn hgvsc_retained_replacement(
    transcript: &Transcript, tid: &str, position: &fastvep_core::GenomicPosition,
    reference: &[u8], alternate: &[u8],
) -> Option<String> {
    let reverse = transcript.strand == fastvep_core::Strand::Reverse;
    let reference = if reverse { fastvep_genome::codon::reverse_complement(reference) } else { reference.to_vec() };
    let alternate = if reverse { fastvep_genome::codon::reverse_complement(alternate) } else { alternate.to_vec() };
    let duplicate = !reference.is_empty() && alternate == reference.repeat(2);
    // VEP classifies the complete replacement before _clip_alleles. Clipping
    // can expose reverse complements, but does not promote delins to inv.
    let inversion = reference.len() > 1
        && alternate == fastvep_genome::codon::reverse_complement(&reference);
    let (mut start, mut end) = (position.start, position.end);
    let (reference, alternate) = if duplicate {
        (reference, alternate)
    } else {
        let prefix = reference.iter().zip(&alternate).take_while(|(a, b)| a == b).count();
        let mut ref_end = reference.len();
        let mut alt_end = alternate.len();
        while ref_end > prefix && alt_end > prefix && reference[ref_end - 1] == alternate[alt_end - 1] {
            ref_end -= 1; alt_end -= 1;
        }
        let suffix = reference.len() - ref_end;
        start = start.checked_add(if reverse { suffix } else { prefix } as u64)?;
        end = end.checked_sub(if reverse { prefix } else { suffix } as u64)?;
        (reference[prefix..ref_end].to_vec(), alternate[prefix..alt_end].to_vec())
    };
    let map = |p| transcript.genomic_to_cdna(p).map(|p| (p, 0))
        .or_else(|| transcript.genomic_to_intronic_cdna(p));
    let (a, b) = (map(start)?, map(end)?);
    let (a, b) = (a.min(b), a.max(b));
    if duplicate {
        return match transcript.cdna_coding_start {
            Some(coding_start) => hgvs_normalize::convert_ins_to_dup_range(
                &format!("{tid}:c."), a, b, coding_start, transcript.cdna_coding_end),
            None => hgvs_normalize::convert_ins_to_dup_range_noncoding(&format!("{tid}:n."), a, b),
        };
    }
    let as_allele = |bases: Vec<u8>| if bases.is_empty() { Allele::Deletion } else { Allele::Sequence(bases) };
    let (reference, alternate) = (as_allele(reference), as_allele(alternate));
    let notation = match transcript.cdna_coding_start {
        Some(coding_start) => fastvep_hgvs::hgvsc_intronic_range(tid, a.0, a.1,
            Some(b.0), Some(b.1), &reference, &alternate, coding_start, transcript.cdna_coding_end),
        None => fastvep_hgvs::hgvsc_noncoding_intronic_range(tid, a.0, a.1,
            Some(b.0), Some(b.1), &reference, &alternate),
    }?;
    if !inversion {
        if let (Some(span), Allele::Sequence(bases)) = (notation.strip_suffix("inv"), &alternate) {
            return Some(format!("{span}delins{}", std::str::from_utf8(bases).ok()?));
        }
    }
    Some(notation)
}

fn hgvsc_for_allele_impl(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &Transcript,
    versioned_tid: &str,
    variation: Option<&VariationFeature>,
    allele: &AlleleConsequenceResult,
) -> Option<String> {
    // VEP TVA::hgvs_transcript rejects ambiguous alternate DNA before clipping.
    if !allele.allele.is_unambiguous_dna() {
        return None;
    }
    let start = allele.normalized_position.start;
    let end = allele.normalized_position.end;
    // VEP checks the complete parsed VariationFeature against the transcript
    // slice before clipping common bases for HGVS. This excludes a padded
    // substitution whose changed base is inside but whose padding is outside.
    let input_position = variation
        .map(|variation| &variation.position)
        .unwrap_or(&allele.normalized_position);
    if input_position.start < transcript.start
        || input_position.start > transcript.end
        || input_position.end < transcript.start
        || input_position.end > transcript.end
    {
        return None;
    }
    // TVA::hgvs_transcript obtains REF from the genomic slice, including when
    // a concrete input REF disagrees with it. Input alleles still determine
    // parsing, variant class and shifting. The batch callers reuse their
    // already-fetched prediction window for this lookup.
    let input_reference = variation.map(|v| &v.ref_allele)
        .unwrap_or(&allele.normalized_ref_allele);
    if seq_provider.is_none() && !input_reference.is_unambiguous_dna() {
        return None;
    }
    let genomic_reference = if input_position.start <= input_position.end {
        seq_provider.map(|provider| provider.fetch_sequence(
            chrom, input_position.start, input_position.end,
        )).transpose().ok()?.map(Allele::Sequence)
    } else { None };
    if let Some(variation) = variation {
        let reference = genomic_reference.as_ref().unwrap_or(&variation.ref_allele);
        if allele.allele == Allele::Deletion
            && !vep_input_is_shiftable_indel(&variation.ref_allele, &variation.alt_alleles)
        {
            return hgvsc_retained_replacement(transcript, versioned_tid,
                &variation.position, reference.as_bytes(), &[]);
        }
        if let (Allele::Sequence(reference), Allele::Sequence(alternate)) =
            (reference, &allele.allele)
        {
            if !reference.is_empty() && !alternate.is_empty()
                && (reference.len() != alternate.len() || reference.len() > 1)
            {
                return hgvsc_retained_replacement(transcript, versioned_tid,
                    &variation.position, reference, alternate);
            }
        }
    }
    let normalized_reference;
    let reference = if allele.normalized_ref_allele == Allele::Deletion || seq_provider.is_none() {
        &allele.normalized_ref_allele
    } else if start == input_position.start && end == input_position.end {
        genomic_reference.as_ref()?
    } else {
        normalized_reference = Allele::Sequence(seq_provider?.fetch_sequence(chrom, start, end).ok()?);
        &normalized_reference
    };
    if reference == &allele.normalized_alt_allele {
        return None;
    }
    let (hgvs_ref, hgvs_alt) = if transcript.strand == fastvep_core::Strand::Reverse {
        (
            complement_allele(reference),
            complement_allele(&allele.normalized_alt_allele),
        )
    } else {
        (
            reference.clone(),
            allele.normalized_alt_allele.clone(),
        )
    };

    // VEP keeps an equal-length padded substitution intact while calculating
    // consequence/display positions, then `_clip_alleles` removes its common
    // prefix and advances the HGVS start. `normalized_position` is that
    // allele-local clipped interval; `allele.cdna_start/end` intentionally
    // remain the parser/display coordinates and must not be reused here.
    let cdna_span = (
        transcript.genomic_to_cdna(start),
        transcript.genomic_to_cdna(end),
    );
    // `hgvs_transcript` uses phase-adjusted TranscriptVariation CDS
    // coordinates only when the parsed record is a one-base SNP. A padded
    // equal-length substitution has class `substitution`; after `_clip_alleles`
    // it is rendered through `_get_cDNA_position`, which has no phase term.
    let is_input_snp = variation
        .is_none_or(|variation| vep_input_is_snp(&variation.ref_allele, &variation.alt_alleles));
    let substitution_phase = if is_input_snp {
        transcript.codon_table_start_phase
    } else {
        0
    };

    if hgvs_normalize::vep_hgvs_shift_exceeds_transcript(
        seq_provider,
        chrom,
        transcript,
        start,
        end,
        &allele.normalized_ref_allele,
        &allele.normalized_alt_allele,
    ) {
        return None;
    }
    match cdna_span {
        (Some(cdna_start), Some(cdna_end)) => {
            let render_exonic = |cdna_start: u64, cdna_end: u64, spliced_seq: Option<&str>| {
                let (cdna_start, cdna_end) = (cdna_start.min(cdna_end), cdna_start.max(cdna_end));
                match transcript.cdna_coding_start {
                    Some(coding_start) => fastvep_hgvs::hgvsc_with_seq(
                        versioned_tid,
                        cdna_start,
                        cdna_end,
                        &hgvs_ref,
                        &hgvs_alt,
                        coding_start,
                        transcript.cdna_coding_end,
                        spliced_seq,
                        substitution_phase,
                    ),
                    None => fastvep_hgvs::hgvsc_noncoding(
                        versioned_tid,
                        cdna_start,
                        cdna_end,
                        &hgvs_ref,
                        &hgvs_alt,
                        spliced_seq,
                    ),
                }
            };
            if seq_provider.is_some()
                && matches!(
                    (&hgvs_ref, &hgvs_alt),
                    (Allele::Sequence(bases), Allele::Deletion) if !bases.is_empty()
                )
            {
                if let Some((shifted_cdna_start, shifted_cdna_end)) =
                    hgvs_normalize::exonic_deletion_cdna_span(
                        seq_provider,
                        chrom,
                        transcript,
                        start,
                        end,
                        &allele.normalized_ref_allele,
                        &allele.normalized_alt_allele,
                    )
                {
                    // The genomic 3'-shift is complete; passing the spliced
                    // sequence would shift a second time and can jump an intron.
                    return render_exonic(shifted_cdna_start, shifted_cdna_end, None);
                }
                return hgvsc_intronic_shifted(
                    seq_provider,
                    chrom,
                    transcript,
                    versioned_tid,
                    start,
                    end,
                    &allele.normalized_ref_allele,
                    &allele.normalized_alt_allele,
                    &hgvs_ref,
                    &hgvs_alt,
                    transcript.cdna_coding_start,
                    transcript.cdna_coding_end,
                );
            }
            if seq_provider.is_some()
                && matches!(
                    (&hgvs_ref, &hgvs_alt),
                    (Allele::Deletion, Allele::Sequence(bases)) if !bases.is_empty()
                )
            {
                let provider = seq_provider?;
                let (shifted_start, shifted_end) = hgvs_normalize::three_prime_shift_intronic(
                    provider,
                    chrom,
                    start,
                    end,
                    &allele.normalized_ref_allele,
                    &allele.normalized_alt_allele,
                    transcript.strand,
                    transcript.start,
                    transcript.end,
                );
                let shifted_crosses_boundary = transcript.genomic_to_cdna(shifted_start).is_none()
                    || transcript.genomic_to_cdna(shifted_end).is_none();
                let shift_crosses_intron = transcript.genomic_to_cdna(start)
                    .zip(transcript.genomic_to_cdna(shifted_start))
                    .is_some_and(|(a, b)| a.abs_diff(b) != start.abs_diff(shifted_start));
                // Inspect the prospective duplicated block even when it does
                // not match genomic sequence. Splicing can otherwise create
                // a false duplication across an exon junction (F1 USP13).
                let duplicated_block_crosses_boundary = match &allele.normalized_alt_allele {
                    Allele::Sequence(inserted) if !inserted.is_empty() => match transcript.strand {
                        fastvep_core::Strand::Forward => shifted_start.checked_sub(inserted.len() as u64)
                            .zip(shifted_start.checked_sub(1)),
                        fastvep_core::Strand::Reverse => shifted_start.checked_add(inserted.len() as u64 - 1)
                            .map(|end| (shifted_start, end)),
                    }
                    .is_some_and(|(lo, hi)| {
                        transcript.genomic_to_cdna(lo).zip(transcript.genomic_to_cdna(hi))
                            .is_none_or(|(a, b)| a.abs_diff(b) != hi - lo)
                    }),
                    _ => false,
                };
                if shifted_crosses_boundary || shift_crosses_intron || duplicated_block_crosses_boundary {
                    // Ensembl transcripts use `_genomic_shift`, even when both
                    // unshifted insertion flanks are exonic. If either the
                    // shift, shifted point or duplicated block reaches an intron,
                    // the spliced transcript names another event.
                    return hgvsc_intronic_shifted(
                        seq_provider,
                        chrom,
                        transcript,
                        versioned_tid,
                        start,
                        end,
                        &allele.normalized_ref_allele,
                        &allele.normalized_alt_allele,
                        &hgvs_ref,
                        &hgvs_alt,
                        transcript.cdna_coding_start,
                        transcript.cdna_coding_end,
                    );
                }
            }
            // Non-coding sequences are intentionally absent from the cache.
            // Fetch one only when an insertion needs the HGVS 3'-rule and
            // duplication detection, then drop it after this annotation.
            let transient_spliced = if transcript.cdna_coding_start.is_none()
                && transcript.spliced_seq.is_none()
                && matches!(
                    (&hgvs_ref, &hgvs_alt),
                    (Allele::Deletion, Allele::Sequence(bases)) if !bases.is_empty()
                ) {
                hgvs_normalize::transient_spliced_sequence(seq_provider, transcript)
            } else {
                None
            };
            let spliced_seq = transcript
                .spliced_seq
                .as_deref()
                .or(transient_spliced.as_deref());
            render_exonic(cdna_start, cdna_end, spliced_seq)
        }
        _ => hgvsc_intronic_shifted(
            seq_provider,
            chrom,
            transcript,
            versioned_tid,
            start,
            end,
            &allele.normalized_ref_allele,
            &allele.normalized_alt_allele,
            &hgvs_ref,
            &hgvs_alt,
            transcript.cdna_coding_start,
            transcript.cdna_coding_end,
        ),
    }
}

/// Signed genomic displacement applied while rendering the transcript HGVS.
/// VEP reports positive values for shifts toward increasing genomic
/// coordinates and negative values for shifts toward decreasing coordinates.
pub fn hgvs_offset_for_allele(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &Transcript,
    variation: &VariationFeature,
    allele: &AlleleConsequenceResult,
) -> Option<i64> {
    // VEP's HGVS block runs only for `within_feature` rows. FastVEP records a
    // distance only for upstream/downstream consequences, so this is the same
    // predicate without coupling HGVS_OFFSET to whether HGVSc was renderable.
    // TranscriptVariationAllele::_genomic_shift accepts only whole-record
    // insertion/deletion classes, not an indel extracted from a mixed record.
    if allele.distance.is_some()
        || !vep_input_is_shiftable_indel(&variation.ref_allele, &variation.alt_alleles)
    {
        return None;
    }
    let provider = seq_provider?;
    let is_nonempty_indel = match (&allele.normalized_ref_allele, &allele.normalized_alt_allele) {
        (Allele::Sequence(bases), Allele::Deletion)
        | (Allele::Deletion, Allele::Sequence(bases)) => !bases.is_empty(),
        _ => false,
    };
    if !is_nonempty_indel {
        return None;
    }
    let start = allele.normalized_position.start;
    let (shifted_start, _) = hgvs_normalize::three_prime_shift_intronic(
        provider,
        chrom,
        start,
        allele.normalized_position.end,
        &allele.normalized_ref_allele,
        &allele.normalized_alt_allele,
        transcript.strand,
        transcript.start,
        transcript.end,
    );
    // VEP 115 minimizes the parsed REF/ALT pair from the left before the right.
    // The consequence engine uses the opposite order because it directly gives
    // the transcript-most-3' representation. Both reach the same final HGVS,
    // but the reported offset is measured from VEP's left-first baseline.
    vep_hgvs_offset_from_shifted_start(
        shifted_start,
        variation.position.start,
        &variation.ref_allele,
        &allele.allele,
    )
}

fn vep_input_is_shiftable_indel(reference: &Allele, alternates: &[Allele]) -> bool {
    !alternates.is_empty() && match reference {
        Allele::Deletion => alternates.iter().all(|alt|
            matches!(alt, Allele::Sequence(bases) if !bases.is_empty())),
        Allele::Sequence(bases) if !bases.is_empty() =>
            alternates.iter().all(|alt| *alt == Allele::Deletion),
        _ => false,
    }
}

fn vep_hgvs_baseline_start(start: u64, reference: &Allele, alternate: &Allele) -> u64 {
    let (Allele::Sequence(reference), Allele::Sequence(alternate)) = (reference, alternate) else {
        return start;
    };
    if reference.len() == alternate.len() {
        return start;
    }

    let prefix = reference
        .iter()
        .zip(alternate)
        .take_while(|(left, right)| left.eq_ignore_ascii_case(right))
        .count();
    start.saturating_add(prefix as u64)
}

fn vep_hgvs_offset_from_shifted_start(
    shifted_start: u64,
    start: u64,
    reference: &Allele,
    alternate: &Allele,
) -> Option<i64> {
    let baseline_start = vep_hgvs_baseline_start(start, reference, alternate);
    let offset = i128::from(shifted_start) - i128::from(baseline_start);
    i64::try_from(offset).ok().filter(|offset| *offset != 0)
}

#[cfg(test)]
mod allele_tests {
    use super::*;

    #[test]
    fn ambiguous_input_reference_uses_genome_for_transcript_hgvs() {
        struct Genome;
        impl SequenceProvider for Genome {
            fn fetch_sequence(&self, _: &str, start: u64, end: u64) -> anyhow::Result<Vec<u8>> {
                Ok(b"ATGAAATAA"[(start - 1) as usize..end as usize].to_vec())
            }
        }
        let gff = "chr1\ttest\tgene\t1\t9\t.\t+\t.\tID=gene:G;biotype=protein_coding\n\
chr1\ttest\tmRNA\t1\t9\t.\t+\t.\tID=transcript:T;Parent=gene:G;biotype=protein_coding\n\
chr1\ttest\texon\t1\t9\t.\t+\t.\tParent=transcript:T;rank=1\n\
chr1\ttest\tCDS\t1\t9\t.\t+\t0\tParent=transcript:T;protein_id=P\n";
        let mut transcripts = parse_gff3(gff.as_bytes()).unwrap();
        let tr = &mut transcripts[0];
        tr.build_sequences(|chrom, start, end| Genome.fetch_sequence(chrom, start, end)
            .map_err(|error| error.to_string())).unwrap();
        for (reference, alternate, expected) in [
            ("AN", "AT", None), ("AN", "AC", Some("T:c.2T>C")),
            ("N", "T", Some("T:c.1A>T")),
            ("AC", "AT", None), ("AC", "AG", Some("T:c.2T>G")),
            ("C", "T", Some("T:c.1A>T")),
        ] {
            let variation = fastvep_io::vcf::parse_vcf_line(
                &format!("chr1\t1\t.\t{reference}\t{alternate}\t.\tPASS\t."),
            ).unwrap();
            let prediction = ConsequencePredictor::default().predict_with_parsed_input(
                &variation.position, &variation.ref_allele, &variation.alt_alleles,
                &[tr], None, variation.minimised,
            );
            let allele = &prediction.transcript_consequences[0].allele_consequences[0];
            assert_eq!(hgvsc_for_variation_allele(Some(&Genome), "chr1", tr, "T",
                &variation, allele).as_deref(), expected);
        }
        let variation = fastvep_io::vcf::parse_vcf_line(
            "chr1\t1\t.\tAN\tA\t.\tPASS\t.",
        ).unwrap();
        let prediction = ConsequencePredictor::default().predict_with_parsed_input(
            &variation.position, &variation.ref_allele, &variation.alt_alleles,
            &[tr], None, variation.minimised,
        );
        let allele = &prediction.transcript_consequences[0].allele_consequences[0];
        assert_eq!(hgvsp_for_variation_allele_with_offset(None, "chr1", tr, "P",
            allele, &variation, None, None), None);
    }

    #[test]
    fn intergenic_annotation_discards_reference_alt_without_rewriting_the_record() {
        let mut vf = fastvep_io::vcf::parse_vcf_line("chr1\t10\t.\tG\tG,N\t.\tPASS\t.").unwrap();
        annotate_intergenic(&mut vf);
        assert_eq!(vf.alt_alleles.len(), 2);
        assert_eq!(vf.transcript_variations.len(), 1);
        assert_eq!(vf.transcript_variations[0].allele_annotations[0].allele, Allele::from_str("N"));
    }

    #[test]
    fn clipped_start_loss_uses_full_reference_flanks() {
        for sign in ["+", "-"] {
            let (coding_start, coding_end) = if sign == "+" { (11, 22) } else { (1, 12) };
            let gff = format!("chr1\ttest\tgene\t1\t22\t.\t{sign}\t.\tID=gene:G;biotype=protein_coding\n\
chr1\ttest\tmRNA\t1\t22\t.\t{sign}\t.\tID=transcript:T;Parent=gene:G;biotype=protein_coding\n\
chr1\ttest\texon\t1\t22\t.\t{sign}\t.\tParent=transcript:T;rank=1\n\
chr1\ttest\tCDS\t{coding_start}\t{coding_end}\t.\t{sign}\t0\tParent=transcript:T;protein_id=P\n");
            let mut transcripts = parse_gff3(gff.as_bytes()).unwrap();
            let tr = &mut transcripts[0];
            let sequence = if sign == "+" { b"CCCCCCCCCCCTGGACCAATAA".to_vec() }
                else { fastvep_genome::codon::reverse_complement(b"CCCCCCCCCCCTGGACCAATAA") };
            tr.build_sequences(|_, start, end| Ok(sequence[start as usize - 1..end as usize].to_vec())).unwrap();
            tr.peptide = Some("LDQ*".into());
            tr.reference_peptide = Some("MDQ".into());
            for (forward, reverse, alternate) in [
                ("chr1\t12\t.\tT\tA,TACG\t.\tPASS\t.", "chr1\t11\t.\tA\tT,CGTA\t.\tPASS\t.", "LR"),
                ("chr1\t11\t.\tC\tA,CTGA\t.\tPASS\t.", "chr1\t12\t.\tG\tT,TCAG\t.\tPASS\t.", "LM"),
            ] {
            let row = if sign == "+" { forward } else { reverse };
            let variation = fastvep_io::vcf::parse_vcf_line(row).unwrap();
            let prediction = ConsequencePredictor::default().predict_with_parsed_input(
                &variation.position, &variation.ref_allele, &variation.alt_alleles,
                &[tr], None, variation.minimised,
            );
            let allele = &prediction.transcript_consequences[0].allele_consequences[1];
            assert_eq!(allele.amino_acids, Some(("L".into(), alternate.into())));
            assert!(allele.consequences.contains(&Consequence::StartLost));
            assert_eq!(hgvsp_for_variation_allele_with_offset(None, "chr1", tr, "P",
                allele, &variation, None, None).as_deref(), Some("P:p.MetAsp2_?1"));
            }
        }
    }


    #[test]
    fn retained_replacement_hgvs_clips_without_a_standalone_indel_shift() {
        use fastvep_core::{GenomicPosition, Strand};
        for (sign, strand) in [("+", Strand::Forward), ("-", Strand::Reverse)] {
            let gff = format!("chr1\ttest\tgene\t1\t30\t.\t{sign}\t.\tID=gene:G;biotype=protein_coding\n\
chr1\ttest\tmRNA\t1\t30\t.\t{sign}\t.\tID=transcript:T;Parent=gene:G;biotype=protein_coding\n\
chr1\ttest\texon\t1\t30\t.\t{sign}\t.\tParent=transcript:T;rank=1\n\
chr1\ttest\tCDS\t1\t30\t.\t{sign}\t0\tParent=transcript:T;protein_id=P\n");
            let transcripts = parse_gff3(gff.as_bytes()).unwrap();
            let tr = &transcripts[0];
            let variation = fastvep_io::vcf::parse_vcf_line(
                "chr1\t29\t.\tCC\tGG,CCCC\t.\t.\t.",
            ).unwrap();
            let prediction = fastvep_consequence::ConsequencePredictor::default().predict(
                &variation.position, &variation.ref_allele, &variation.alt_alleles, &[tr], None,
            );
            let ranges = vep_position_ranges(tr,
                &prediction.transcript_consequences[0].allele_consequences[1], &variation);
            let expected = if strand == Strand::Forward { (29, 30) } else { (1, 2) };
            assert_eq!(ranges.0.start(), Some(expected.0));
            assert_eq!(ranges.0.end(), Some(expected.1));
            // A shared anchor can leave ALT '-' while the other ALT keeps
            // the whole record in VEP's non-shiftable sequence-alteration class.
            let mixed = fastvep_io::vcf::parse_vcf_line(
                "chr1\t10\t.\tCC\tC,CT\t.\t.\t.",
            ).unwrap();
            let predicted = ConsequencePredictor::default().predict(
                &mixed.position, &mixed.ref_allele, &mixed.alt_alleles, &[tr], None,
            );
            let deletion = &predicted.transcript_consequences[0].allele_consequences[0];
            assert_eq!(deletion.allele, Allele::Deletion);
            assert_eq!(hgvsc_for_variation_allele(None, "chr1", tr, "T", &mixed, deletion),
                Some(if strand == Strand::Forward { "T:c.11del" } else { "T:c.20del" }.into()));
            for missing_start in [true, false] {
                let mut gap = deletion.clone();
                if missing_start { gap.protein_start = None; }
                else { gap.protein_end = None; }
                assert_eq!(hgvsp_for_variation_allele_with_offset(None, "chr1", tr,
                    "P", &gap, &mixed, Some("T:c.11del"), None), None);
            }
            for (reference, alternate, forward, reverse) in [
                ("C", "CC", "T:c.10dup", "T:c.21dup"),
                ("C", "CCC", "T:c.10_11insCC", "T:c.21_22insGG"),
                ("C", "CAA", "T:c.10_11insAA", "T:c.20_21insTT"),
                ("CC", "C", "T:c.11del", "T:c.21del"),
                ("AGCT", "GCTT", "T:c.10_12delinsGCT", "T:c.19_21delinsAGC"),
                ("AGC", "GCT", "T:c.10_12inv", "T:c.19_21inv"),
            ] {
                let position = GenomicPosition::new("chr1", 10, 9 + reference.len() as u64, Strand::Forward);
                assert_eq!(hgvsc_retained_replacement(tr, "T", &position, reference.as_bytes(), alternate.as_bytes()),
                    Some(if strand == Strand::Forward { forward } else { reverse }.to_string()));
                if reference.len() == alternate.len() && reference.len() > 1 {
                    let variation = fastvep_io::vcf::parse_vcf_line(
                        &format!("chr1\t10\t.\t{reference}\t{alternate}\t.\t.\t."),
                    ).unwrap();
                    let prediction = ConsequencePredictor::default().predict(
                        &variation.position, &variation.ref_allele, &variation.alt_alleles, &[tr], None,
                    );
                    let allele = &prediction.transcript_consequences[0].allele_consequences[0];
                    assert_eq!(hgvsc_for_variation_allele(None, "chr1", tr, "T", &variation, allele),
                        Some(if strand == Strand::Forward { forward } else { reverse }.to_string()));
                }
            }
        }
    }

    #[test]
    fn hgvs_shift_requires_a_whole_record_insertion_or_deletion() {
        for (reference, alternates, expected) in [
            ("-", vec!["A", "AA"], true),
            ("AA", vec!["-"], true),
            ("C", vec!["A", "CAA"], false),
            ("C", vec!["CAA", "A"], false),
            ("CA", vec!["C", "T"], false),
            ("CA", vec!["C", "CAA"], false),
            ("-", vec![], false),
        ] {
            assert_eq!(vep_input_is_shiftable_indel(&Allele::from_str(reference),
                &alternates.iter().map(|alt| Allele::from_str(alt)).collect::<Vec<_>>()), expected);
        }
    }

    #[test]
    fn mixed_multiallelic_frameshift_uses_retained_replacement_interval() {
        use fastvep_core::{GenomicPosition, Strand};
        // ARVCF CDS sequence from the VEP source trace, with the affected Val
        // renumbered to residue 2. The same replacement is tested on both strands.
        let cds = b"ATGGTGGACCCCGTGAAGGCCAATGCGGCCGCCTACCTGCAGCATCTGTGCTTTGAGAACGAG";
        let length = cds.len() as u64;
        for strand in [Strand::Forward, Strand::Reverse] {
            let sign = if strand == Strand::Forward { "+" } else { "-" };
            let gff = format!("chr1\ttest\tgene\t1\t{length}\t.\t{sign}\t.\tID=gene:G;biotype=protein_coding\n\
chr1\ttest\tmRNA\t1\t{length}\t.\t{sign}\t.\tID=transcript:T;Parent=gene:G;biotype=protein_coding\n\
chr1\ttest\texon\t1\t{length}\t.\t{sign}\t.\tParent=transcript:T;rank=1\n\
chr1\ttest\tCDS\t1\t{length}\t.\t{sign}\t0\tParent=transcript:T;protein_id=P\n");
            let genomic = if strand == Strand::Forward { cds.to_vec() }
                else { fastvep_genome::codon::reverse_complement(cds) };
            let mut transcripts = fastvep_cache::gff::parse_gff3(gff.as_bytes()).unwrap();
            let tr = &mut transcripts[0];
            tr.build_sequences(|_, start, end| Ok(genomic[start as usize - 1..end as usize].to_vec())).unwrap();
            let (position, reference, alternate) = if strand == Strand::Forward {
                (4, Allele::from_str("G"), Allele::from_str("TTG"))
            } else {
                (length - 3, Allele::from_str("C"), Allele::from_str("CAA"))
            };
            let result = ConsequencePredictor::default().predict(
                &GenomicPosition::new("chr1", position, position, Strand::Forward),
                &reference, &[alternate, Allele::from_str("A")], &[tr], None,
            );
            let allele = &result.transcript_consequences[0].allele_consequences[0];
            assert_eq!((allele.cds_start, allele.cds_end), (Some(4), Some(4)));
            assert_eq!(allele.normalized_ref_allele, Allele::Deletion);
            assert_eq!(hgvsp_for_allele_with_offset(None, "chr1", tr, "P", allele,
                &reference, Some("T:c.3_4insTT"), None), Some("P:p.Val2LeufsTer5".into()));
            assert_eq!(hgvsp_for_allele_with_offset(None, "chr1", tr, "P", allele,
                &reference, Some("T:c.3_4+1del"), None), Some("P:p.Val2LeufsTer5".into()));
            for missing_start in [true, false] {
                let mut gap = allele.clone();
                if missing_start { gap.protein_start = None; }
                else { gap.protein_end = None; }
                assert_eq!(hgvsp_for_allele_with_offset(None, "chr1", tr, "P", &gap,
                    &reference, Some("T:c.3_4insTT"), None), None);
            }
        }
    }


    #[test]
    fn protein_formatting_keeps_hidden_stop_predicates_and_internal_stops() {
        use fastvep_core::Strand;
        for strand in [Strand::Forward, Strand::Reverse] {
            for (sequence, cds_length, lo, reference, alternate, ablation, expected) in [
                ("ATGGGATGA", 9, 1, "ATGGGATGA", "", true, "P:p.Met1_Ter3delextTer?"),
                ("ATGGCTTAAGCTGAA", 14, 4, "GCTTAAGCTGA", "GCTTAAGCTGT", false, "P:p.Ala2_Ter5delinsAlaTer"),
            ] {
                let length = sequence.len() as u64;
                let (symbol, coding_lo, coding_hi, start, end) = if strand == Strand::Forward {
                    ("+", 1, cds_length, lo, cds_length)
                } else {
                    ("-", length - cds_length + 1, length, length - cds_length + 1, length - lo + 1)
                };
                let gff = format!("chr1\ttest\tgene\t1\t{length}\t.\t{symbol}\t.\tID=gene:G;biotype=protein_coding\nchr1\ttest\tmRNA\t1\t{length}\t.\t{symbol}\t.\tID=transcript:T;Parent=gene:G;biotype=protein_coding\nchr1\ttest\texon\t1\t{length}\t.\t{symbol}\t.\tParent=transcript:T;rank=1\nchr1\ttest\tCDS\t{coding_lo}\t{coding_hi}\t.\t{symbol}\t0\tParent=transcript:T;protein_id=P\n");
                let orient = |bases: &[u8]| if strand == Strand::Reverse {
                    fastvep_genome::codon::reverse_complement(bases)
                } else { bases.to_vec() };
                let genomic = orient(sequence.as_bytes());
                let mut transcripts = fastvep_cache::gff::parse_gff3(gff.as_bytes()).unwrap();
                let tr = &mut transcripts[0];
                tr.build_sequences(|_, start, end| Ok(genomic[start as usize - 1..end as usize].to_vec())).unwrap();
                let reference = Allele::Sequence(orient(reference.as_bytes()));
                let alternate = if alternate.is_empty() { Allele::Deletion } else { Allele::Sequence(orient(alternate.as_bytes())) };
                let position = fastvep_core::GenomicPosition::new("chr1", start, end, Strand::Forward);
                let result = ConsequencePredictor::default().predict(&position, &reference, &[alternate], &[tr], None);
                let allele = &result.transcript_consequences[0].allele_consequences[0];
                if ablation { assert_eq!(allele.consequences, vec![Consequence::TranscriptAblation]); }
                assert_eq!(hgvsp_with_record_class(None, "chr1", tr, "P", allele, &reference, None, None, false).as_deref(), Some(expected), "{strand:?}, ablation {ablation}");
            }
        }
    }

    #[test]
    fn stop_loss_after_an_unchanged_residue_finds_the_extension_stop() {
        let gff = "chr1\ttest\tgene\t1\t24\t.\t+\t.\tID=gene:G;biotype=protein_coding\n\
chr1\ttest\tmRNA\t1\t24\t.\t+\t.\tID=transcript:T;Parent=gene:G;biotype=protein_coding\n\
chr1\ttest\texon\t1\t24\t.\t+\t.\tParent=transcript:T;rank=1\n\
chr1\ttest\tCDS\t1\t9\t.\t+\t0\tParent=transcript:T;protein_id=P\n";
        let mut transcripts = fastvep_cache::gff::parse_gff3(gff.as_bytes()).unwrap();
        let tr = &mut transcripts[0];
        let sequence = b"ATGGGATGAGCTGCTGCTGCTTAA";
        tr.build_sequences(|_, start, end| Ok(sequence[start as usize - 1..end as usize].to_vec()))
            .unwrap();
        let reference = Allele::from_str("AT");
        let position =
            fastvep_core::GenomicPosition::new("chr1", 6, 7, fastvep_core::Strand::Forward);
        let result = ConsequencePredictor::default().predict(
            &position,
            &reference,
            &[Allele::from_str("CA")],
            &[tr],
            None,
        );
        let allele = &result.transcript_consequences[0].allele_consequences[0];
        assert_eq!(allele.amino_acids, Some(("G*".into(), "GR".into())));
        assert_eq!(
            hgvsp_for_allele_with_offset(
                None,
                "chr1",
                tr,
                "P",
                allele,
                &reference,
                Some("T:c.6_7delinsCA"),
                None
            ),
            Some("P:p.Ter3ArgextTer5".into())
        );
    }

    #[test]
    fn internal_exon_edge_insertion_retains_protein_hgvs() {
        let gff = "chr1\ttest\tgene\t1\t27\t.\t+\t.\tID=gene:G;biotype=protein_coding\n\
chr1\ttest\tmRNA\t1\t27\t.\t+\t.\tID=transcript:T;Parent=gene:G;biotype=protein_coding\n\
chr1\ttest\texon\t1\t9\t.\t+\t.\tParent=transcript:T;rank=1\n\
chr1\ttest\texon\t19\t27\t.\t+\t.\tParent=transcript:T;rank=2\n\
chr1\ttest\tCDS\t1\t9\t.\t+\t0\tParent=transcript:T;protein_id=P\n\
chr1\ttest\tCDS\t19\t27\t.\t+\t0\tParent=transcript:T;protein_id=P\n";
        let mut transcripts = fastvep_cache::gff::parse_gff3(gff.as_bytes()).unwrap();
        let tr = &mut transcripts[0];
        let sequence = b"ATGGCTGCTNNNNNNNNNGCTGCTTAA";
        tr.build_sequences(|_, start, end| Ok(sequence[start as usize - 1..end as usize].to_vec()))
            .unwrap();
        let position =
            fastvep_core::GenomicPosition::new("chr1", 10, 9, fastvep_core::Strand::Forward);
        let result = ConsequencePredictor::new(5000, 5000).predict(
            &position,
            &Allele::Deletion,
            &[Allele::from_str("T")],
            &[tr],
            None,
        );
        let allele = &result.transcript_consequences[0].allele_consequences[0];
        // Ensembl 115 Mapper::map_insert drops the intronic flank's gap.
        assert!(
            hgvsp_for_allele_with_offset(
                None,
                "chr1",
                tr,
                "P",
                allele,
                &Allele::Deletion,
                Some("T:c.9_9+1insT"),
                None
            )
            .is_some(),
            "{allele:?}; {tr:?}"
        );
        // A repeated REF ALT makes VEP classify -/-/T as sequence_alteration,
        // but Mapper::map_insert still drops the intronic flank's gap.
        assert!(hgvsp_with_record_class(None, "chr1", tr, "P", allele,
            &Allele::Deletion, Some("T:c.9_9+1insT"), None, false).is_some());
    }

    #[test]
    fn reverse_complements_multibase_alleles() {
        assert_eq!(
            complement_allele(&Allele::from_str("AGC")),
            Allele::from_str("GCT")
        );
    }

    #[test]
    fn vep_snp_class_requires_one_base_for_every_record_allele() {
        assert!(vep_input_is_snp(
            &Allele::from_str("A"),
            &[Allele::from_str("C"), Allele::from_str("G")]
        ));
        assert!(!vep_input_is_snp(
            &Allele::from_str("CA"),
            &[Allele::from_str("CC")]
        ));
        assert!(!vep_input_is_snp(
            &Allele::from_str("A"),
            &[Allele::from_str("C"), Allele::from_str("AT")]
        ));
        assert!(!vep_input_is_snp(&Allele::from_str("A"), &[]));
    }

    #[test]
    fn vep_hgvs_offset_baseline_trims_repeated_indels_from_the_left() {
        assert_eq!(
            vep_hgvs_baseline_start(100, &Allele::from_str("TT"), &Allele::from_str("T")),
            101
        );
        assert_eq!(
            vep_hgvs_baseline_start(100, &Allele::from_str("AAA"), &Allele::from_str("AAAA")),
            103
        );
        assert_eq!(
            vep_hgvs_baseline_start(100, &Allele::from_str("AC"), &Allele::from_str("AG")),
            100
        );
        assert_eq!(
            vep_hgvs_offset_from_shifted_start(
                100,
                100,
                &Allele::from_str("TT"),
                &Allele::from_str("T")
            ),
            Some(-1)
        );
        assert_eq!(
            vep_hgvs_offset_from_shifted_start(
                118,
                100,
                &Allele::from_str("AAA"),
                &Allele::from_str("AAAA")
            ),
            Some(15)
        );
    }

    #[test]
    fn protein_hgvs_follows_veps_asymmetric_boundary_shift_fallback() {
        assert!(hgvsc_has_vep_protein_endpoint("ENST1:c.243_244del"));
        assert!(!hgvsc_has_vep_protein_endpoint("ENST1:c.243_243+1insT"));
        assert!(hgvsc_has_vep_protein_endpoint("ENST1:c.243-1_243del"));
        assert!(!hgvsc_has_vep_protein_endpoint("ENST1:c.1097+2_1097+5dup"));
        assert!(!hgvsc_has_vep_protein_endpoint("ENST1:c.*1del"));
        assert!(!hgvsc_has_vep_protein_endpoint("ENST1:n.42A>G"));
    }

    #[test]
    fn shifted_translation_endpoint_must_remain_inside_the_cds() {
        assert!(!shifted_translation_end_is_coding(
            Some(144),
            Some(142),
            Some(143),
            2
        ));
        assert!(shifted_translation_end_is_coding(
            Some(144),
            Some(142),
            Some(143),
            1
        ));
    }

    #[test]
    fn vep_frameshift_intron_definition_includes_thirteen_bases() {
        let exon = |start, end| Exon {
            stable_id: String::new(),
            start,
            end,
            strand: fastvep_core::Strand::Forward,
            phase: 0,
            end_phase: 0,
            rank: 1,
        };
        assert!(vep_has_frameshift_intron(
            &[exon(100, 199), exon(213, 250),]
        ));
        assert!(!vep_has_frameshift_intron(&[
            exon(100, 199),
            exon(214, 250),
        ]));
    }
}

/// Extract trio genotype information from a VariationFeature's VCF sample columns.
///
/// Returns (proband, mother, father) GenotypeInfo tuples.
fn extract_trio_genotypes(
    vf: &VariationFeature,
    acmg_cfg: &fastvep_classification::AcmgConfig,
    sample_names: &[String],
) -> (
    Option<fastvep_classification::GenotypeInfo>,
    Option<fastvep_classification::GenotypeInfo>,
    Option<fastvep_classification::GenotypeInfo>,
) {
    let trio = match &acmg_cfg.trio {
        Some(t) => t,
        None => return (None, None, None),
    };

    let vcf_fields = match &vf.vcf_fields {
        Some(f) => f,
        None => return (None, None, None),
    };

    // rest[0] is FORMAT, rest[1..] are sample columns
    if vcf_fields.rest.is_empty() {
        return (None, None, None);
    }

    let format_str = &vcf_fields.rest[0];
    let sample_strs: Vec<&str> = vcf_fields.rest[1..].iter().map(|s| s.as_str()).collect();

    let samples = fastvep_io::sample::parse_samples(format_str, &sample_strs, sample_names);

    let proband_gt = samples
        .iter()
        .find(|s| s.name == trio.proband)
        .map(|s| sample_data_to_genotype_info(s));

    let mother_gt = trio.mother.as_ref().and_then(|name| {
        samples
            .iter()
            .find(|s| &s.name == name)
            .map(|s| sample_data_to_genotype_info(s))
    });

    let father_gt = trio.father.as_ref().and_then(|name| {
        samples
            .iter()
            .find(|s| &s.name == name)
            .map(|s| sample_data_to_genotype_info(s))
    });

    (proband_gt, mother_gt, father_gt)
}

/// Convert a SampleData to GenotypeInfo.
fn sample_data_to_genotype_info(
    sample: &fastvep_io::sample::SampleData,
) -> fastvep_classification::GenotypeInfo {
    let gt = sample.genotype.as_ref();
    let is_het = gt.map_or(false, |g| g.is_het());
    let is_hom_ref = gt.map_or(false, |g| g.is_hom_ref());
    let is_hom_alt = gt.map_or(false, |g| g.is_hom_alt());
    let is_missing = gt.map_or(true, |g| g.is_missing());
    let is_phased = gt.map_or(false, |g| g.phased);

    // Determine which alt allele index is carried
    let alt_allele_index = gt.and_then(|g| {
        g.alleles
            .iter()
            .filter_map(|a| *a)
            .find(|&a| a > 0)
            .map(|a| a)
    });

    fastvep_classification::GenotypeInfo {
        is_het,
        is_hom_ref,
        is_hom_alt,
        is_missing,
        is_phased,
        depth: sample.depth,
        quality: sample.quality,
        alt_allele_index,
    }
}

/// Compound-het enrichment pass: after all variants are annotated,
/// group by gene and identify companion variant relationships,
/// then re-evaluate PM3/BP2 with companion data.
fn enrich_compound_het(
    variants: &mut [VariationFeature],
    acmg_cfg: &fastvep_classification::AcmgConfig,
    sample_names: &[String],
) {
    use std::collections::HashMap;

    // Collect per-gene variant info: (variant_index, gene_symbol, ClinVar P/LP flags, proband_het, is_phased, hgvsc, allele_indices for phase)
    struct VariantGeneInfo {
        vf_idx: usize,
        tv_idx: usize,
        aa_idx: usize,
        is_clinvar_pathogenic: bool,
        is_clinvar_likely_pathogenic: bool,
        proband_het: bool,
        is_phased: bool,
        /// Proband's allele indices for phase comparison
        proband_alleles: Vec<Option<u32>>,
        hgvsc: Option<String>,
    }

    let mut gene_variants: HashMap<String, Vec<VariantGeneInfo>> = HashMap::new();

    for (vf_idx, vf) in variants.iter().enumerate() {
        let trio_genotypes = extract_trio_genotypes(vf, acmg_cfg, sample_names);
        let proband_gt = &trio_genotypes.0;

        for (tv_idx, tv) in vf.transcript_variations.iter().enumerate() {
            let gene_sym = match tv.gene_symbol.as_deref() {
                Some(g) if !g.is_empty() && g != "-" => g.to_string(),
                _ => continue,
            };

            for (aa_idx, aa) in tv.allele_annotations.iter().enumerate() {
                let is_clinvar_pathogenic = aa
                    .acmg_classification
                    .as_ref()
                    .and_then(|v| v.get("criteria"))
                    .and_then(|c| c.as_array())
                    .map_or(false, |criteria| {
                        // Check if this variant has ClinVar pathogenic data
                        criteria.iter().any(|c| {
                            c.get("code")
                                .and_then(|v| v.as_str())
                                .map_or(false, |code| code == "PP5" || code == "PS4")
                                && c.get("met").and_then(|v| v.as_bool()).unwrap_or(false)
                        })
                    });

                // Classify ClinVar supplementary as Pathogenic / Likely pathogenic
                // separately so PM3 v1.0 can score them at their proper point
                // values. A bare substring match on "pathogenic" matches both
                // "Pathogenic" and "Likely pathogenic" — this would over-score
                // LP companions as P. Strip "Likely pathogenic" first, then
                // see if any "pathogenic" remains: that residual signals true P.
                let (clinvar_p_from_sa, clinvar_lp_from_sa) = aa
                    .supplementary
                    .iter()
                    .filter(|(key, json)| {
                        key == "clinvar"
                            && !json.contains("Conflicting")
                            && !json.contains("conflicting")
                    })
                    .map(|(_, json)| {
                        let lower = json.to_lowercase();
                        let has_lp = lower.contains("likely pathogenic");
                        let stripped = lower.replace("likely pathogenic", "");
                        let has_p = stripped.contains("pathogenic");
                        (has_p, has_lp && !has_p)
                    })
                    .fold((false, false), |(p_acc, lp_acc), (p, lp)| {
                        (p_acc || p, lp_acc || lp)
                    });

                let proband_het = proband_gt.as_ref().map_or(false, |g| g.is_het);
                let is_phased = proband_gt.as_ref().map_or(false, |g| g.is_phased);
                let proband_alleles = if let Some(ref vcf_fields) = vf.vcf_fields {
                    if !vcf_fields.rest.is_empty() && !sample_names.is_empty() {
                        let format_str = &vcf_fields.rest[0];
                        let sample_strs: Vec<&str> =
                            vcf_fields.rest[1..].iter().map(|s| s.as_str()).collect();
                        let samples = fastvep_io::sample::parse_samples(
                            format_str,
                            &sample_strs,
                            sample_names,
                        );
                        if let Some(trio) = &acmg_cfg.trio {
                            samples
                                .iter()
                                .find(|s| s.name == trio.proband)
                                .and_then(|s| s.genotype.as_ref())
                                .map(|g| g.alleles.clone())
                                .unwrap_or_default()
                        } else {
                            vec![]
                        }
                    } else {
                        vec![]
                    }
                } else {
                    vec![]
                };

                gene_variants
                    .entry(gene_sym.clone())
                    .or_default()
                    .push(VariantGeneInfo {
                        vf_idx,
                        tv_idx,
                        aa_idx,
                        is_clinvar_pathogenic: is_clinvar_pathogenic || clinvar_p_from_sa,
                        is_clinvar_likely_pathogenic: clinvar_lp_from_sa,
                        proband_het,
                        is_phased,
                        proband_alleles,
                        hgvsc: aa.hgvsc.clone(),
                    });
            }
        }
    }

    // For each gene with multiple het variants, build companion relationships and re-classify
    for (_gene, gene_infos) in &gene_variants {
        let het_variants: Vec<&VariantGeneInfo> =
            gene_infos.iter().filter(|v| v.proband_het).collect();
        if het_variants.len() < 2 {
            continue;
        }

        // For each het variant, build companion list from other het variants in the gene
        for info in &het_variants {
            let companions: Vec<fastvep_classification::CompanionVariant> = het_variants
                .iter()
                .filter(|other| {
                    other.vf_idx != info.vf_idx
                        || other.tv_idx != info.tv_idx
                        || other.aa_idx != info.aa_idx
                })
                .map(|other| {
                    // Determine trans/cis from phase information
                    let is_in_trans = if info.is_phased && other.is_phased {
                        // Both phased: check if they're on different haplotypes
                        // In a phased genotype like 0|1 vs 1|0, alleles at same index
                        // come from the same parent. So het 0|1 and 1|0 means they're
                        // on different haplotypes (trans).
                        if info.proband_alleles.len() >= 2 && other.proband_alleles.len() >= 2 {
                            let info_alt_on_first = info
                                .proband_alleles
                                .first()
                                .map_or(false, |a| a.map_or(false, |v| v > 0));
                            let other_alt_on_first = other
                                .proband_alleles
                                .first()
                                .map_or(false, |a| a.map_or(false, |v| v > 0));
                            // If alt alleles are on different haplotypes, they're in trans
                            Some(info_alt_on_first != other_alt_on_first)
                        } else {
                            None
                        }
                    } else {
                        None
                    };

                    fastvep_classification::CompanionVariant {
                        is_clinvar_pathogenic: other.is_clinvar_pathogenic,
                        is_clinvar_likely_pathogenic: other.is_clinvar_likely_pathogenic,
                        is_in_trans,
                        proband_het: other.proband_het,
                        hgvsc: other.hgvsc.clone(),
                    }
                })
                .collect();

            if companions.is_empty() {
                continue;
            }

            // Re-extract classification input with companion data and re-classify
            let vf = &variants[info.vf_idx];
            let tv = &vf.transcript_variations[info.tv_idx];
            let aa = &tv.allele_annotations[info.aa_idx];
            let gene_sym = tv.gene_symbol.as_deref().unwrap_or("");
            let gene_anns: Vec<&fastvep_core::GeneAnnotation> = vf
                .gene_annotations
                .iter()
                .filter(|ga| ga.gene_symbol == gene_sym)
                .collect();

            let trio_genotypes = extract_trio_genotypes(vf, acmg_cfg, sample_names);

            let input = fastvep_classification::extract_classification_input(
                &aa.consequences,
                aa.impact,
                tv.gene_symbol.as_deref(),
                tv.canonical,
                aa.amino_acids.as_ref(),
                aa.protein_position.first_known(),
                aa.hgvsc.as_deref(),
                aa.exon.map(|(first, _, total)| (first, total)),
                &aa.supplementary,
                &gene_anns,
                &vf.supplementary_annotations,
                trio_genotypes.0,
                trio_genotypes.1,
                trio_genotypes.2,
                companions,
            );
            let result = fastvep_classification::classify(&input, acmg_cfg);
            variants[info.vf_idx].transcript_variations[info.tv_idx].allele_annotations
                [info.aa_idx]
                .acmg_classification = serde_json::to_value(&result).ok();
        }
    }
}

/// Load supplementary annotation providers (.osa, .osa2, .osi files) from a
/// directory.
///
/// Opening is done in parallel and the results are ordered by path (issue
/// #78). A per-chromosome deployment puts 200+ shards in one `--sa-dir`, and
/// opening a shard is latency-bound rather than CPU-bound, so serial opens made
/// startup scale linearly with shard count for no reason. Sorting by path also
/// makes the provider order - and therefore the output column order - depend on
/// the file names rather than on directory iteration order.
pub fn load_sa_providers(sa_dir: &Path) -> Result<Vec<Mutex<Box<dyn AnnotationProvider>>>> {
    use fastvep_sa::interval::OsiReader;
    use fastvep_sa::reader::AnySaReader;
    use rayon::prelude::*;

    if !sa_dir.is_dir() {
        tracing::warn!(
            "SA directory does not exist: {} (skipping)",
            sa_dir.display()
        );
        return Ok(Vec::new());
    }

    let paths = sorted_sa_paths(sa_dir)?;
    let opened: Result<Vec<Option<Mutex<Box<dyn AnnotationProvider>>>>> = paths
        .par_iter()
        .map(|path| {
            if path
                .file_name()
                .and_then(|name| name.to_str())
                .is_some_and(|name| name.ends_with(fastvep_sa::sharded::SHARD_MANIFEST_SUFFIX))
            {
                let reader =
                    fastvep_sa::sharded::ShardedSaReader::open(path).with_context(|| {
                        format!("Loading required OSA shard manifest {}", path.display())
                    })?;
                tracing::info!("Loaded sharded SA: {} ({})", reader.name(), path.display());
                return Ok(Some(Mutex::new(boxed(reader))));
            }

            let ext = path.extension().and_then(|e| e.to_str());
            match ext {
                Some("osa" | "osa2") => {
                    let reader = AnySaReader::open(path).with_context(|| {
                        format!("Loading required supplementary cache {}", path.display())
                    })?;
                    tracing::info!("Loaded SA: {} ({})", reader.name(), path.display());
                    Ok(Some(Mutex::new(boxed(reader))))
                }
                Some("osi") => match OsiReader::open(path) {
                    Ok(reader) => {
                        tracing::info!(
                            "Loaded SA interval: {} ({})",
                            reader.name(),
                            path.display()
                        );
                        Ok(Some(Mutex::new(boxed(reader))))
                    }
                    Err(error) => {
                        tracing::warn!("Could not load {}: {}", path.display(), error);
                        Ok(None)
                    }
                },
                _ => Ok(None),
            }
        })
        .collect();
    let providers = opened?.into_iter().flatten().collect();

    Ok(providers)
}

fn sorted_sa_paths(dir: &Path) -> Result<Vec<std::path::PathBuf>> {
    let mut paths = Vec::new();
    for entry in std::fs::read_dir(dir)
        .with_context(|| format!("Listing annotation directory {}", dir.display()))?
    {
        let path = entry
            .with_context(|| format!("Reading an entry of {}", dir.display()))?
            .path();
        let is_cache = path
            .extension()
            .and_then(|extension| extension.to_str())
            .is_some_and(|extension| matches!(extension, "osa" | "osa2" | "osi"));
        let is_manifest = path
            .file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| name.ends_with(fastvep_sa::sharded::SHARD_MANIFEST_SUFFIX));
        if is_cache || is_manifest {
            paths.push(path);
        }
    }
    paths.sort();
    Ok(paths)
}

fn boxed<P: AnnotationProvider + 'static>(provider: P) -> Box<dyn AnnotationProvider> {
    Box::new(provider)
}

/// Every file in `dir` whose extension is one of `exts`, sorted by path so the
/// caller's provider order is reproducible across machines and filesystems.
///
/// A failed directory entry is propagated rather than skipped: on a network
/// `--sa-dir` a transient error would otherwise drop a shard, and the run would
/// finish successfully with that source's annotations simply missing.
fn sorted_paths_with_extensions(dir: &Path, exts: &[&str]) -> Result<Vec<std::path::PathBuf>> {
    let mut paths: Vec<std::path::PathBuf> = Vec::new();
    for entry in std::fs::read_dir(dir)
        .with_context(|| format!("Listing annotation directory {}", dir.display()))?
    {
        let path = entry
            .with_context(|| format!("Reading an entry of {}", dir.display()))?
            .path();
        if path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| exts.contains(&e))
        {
            paths.push(path);
        }
    }
    paths.sort();
    Ok(paths)
}

/// Load gene-level annotation providers (.oga files) from a directory.
///
/// Same ordering and parallelism contract as [`load_sa_providers`].
pub fn load_gene_providers(sa_dir: &Path) -> Result<Vec<fastvep_sa::gene::GeneIndex>> {
    use rayon::prelude::*;

    if !sa_dir.is_dir() {
        return Ok(Vec::new());
    }

    let paths = sorted_paths_with_extensions(sa_dir, &["oga"])?;

    let providers: Vec<fastvep_sa::gene::GeneIndex> = paths
        .par_iter()
        .filter_map(|path| {
            match std::fs::File::open(path)
                .map_err(anyhow::Error::from)
                .and_then(|mut f| fastvep_sa::gene::GeneIndex::read_from(&mut f))
            {
                Ok(index) => {
                    tracing::info!(
                        "Loaded gene annotations: {} ({}, {} genes)",
                        index.header.name,
                        path.display(),
                        index.gene_count()
                    );
                    Some(index)
                }
                Err(e) => {
                    tracing::warn!("Could not load {}: {}", path.display(), e);
                    None
                }
            }
        })
        .collect();

    Ok(providers)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// fastvep-annotate is the shared engine both fastvep-cli and
    /// fastvep-web sit on top of, and had zero unit tests before this. This
    /// exercises the full annotate_vcf_text path without needing a real
    /// genome: an empty transcript set (gff3 = None) makes every variant
    /// intergenic, which is enough to verify the pipeline runs end-to-end
    /// and produces well-formed output.
    fn empty_context() -> AnnotationContext {
        AnnotationContext::new(None, None, None, 0).expect("empty context should build")
    }

    #[test]
    fn position_ranges_keep_unknown_endpoints() {
        assert_eq!(
            zip_positions(Some(10), None),
            PositionRange::new(Some(10), None)
        );
        assert_eq!(
            zip_positions(None, Some(20)),
            PositionRange::new(None, Some(20))
        );
        assert_eq!(
            zip_positions(Some(20), Some(10)),
            PositionRange::complete(10, 20)
        );
    }

    #[test]
    fn annotate_vcf_text_returns_one_result_per_variant() {
        let ctx = empty_context();
        let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                   1\t100\t.\tA\tG\t.\tPASS\t.\n\
                   1\t200\t.\tC\tT\t.\tPASS\t.\n";
        let results = ctx.annotate_vcf_text(vcf, false).unwrap();
        assert_eq!(results.len(), 2);
        assert_eq!(
            results[0]["most_severe_consequence"],
            serde_json::json!("intergenic_variant")
        );
    }

    #[test]
    fn annotate_vcf_text_with_acmg_none_disables_classification_regardless_of_self_config() {
        let mut ctx = empty_context();
        // self.acmg_config is enabled...
        ctx.acmg_config = Some(fastvep_classification::AcmgConfig::default());
        let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                   1\t100\t.\tA\tG\t.\tPASS\t.\n";

        // ...but an explicit None override (as fastvep-web now passes for
        // acmg_requested=false) must take precedence over self.acmg_config,
        // since concurrent requests share one AnnotationContext and must not
        // leak each other's ACMG preference.
        let results = ctx.annotate_vcf_text_with_acmg(vcf, false, None).unwrap();
        assert_eq!(results.len(), 1);
    }

    #[test]
    fn hgvsp_frameshift_skips_gracefully_on_inconsistent_spliced_seq_len() {
        // Regression: the HGVSp-frameshift branch slices
        // `spliced.as_bytes()[coding_start_idx..]` where `coding_start_idx`
        // comes from `tr.cdna_coding_start`. If `spliced_seq` is ever shorter
        // than `coding_start_idx` implies -- e.g. malformed/truncated GFF3-
        // or cache-derived transcript data -- this must not panic ("range
        // start index out of range"); it should just skip HGVSp generation.
        //
        // `cdna_coding_start`/`cdna_coding_end` (used by the predictor, via
        // `cdna_to_cds`, to classify the variant and via `translateable_seq`
        // to compute amino acids) are left internally consistent, so the
        // variant is still correctly classified as a frameshift with real
        // amino acids computed -- only `spliced_seq` itself (used solely by
        // fastvep-annotate's separate HGVSp-slicing code, not by the
        // predictor) is corrupted, isolating the guard under test.
        use fastvep_genome::{Exon, Gene, Transcript, Translation};

        // Two-exon transcript on chr1: exon1 is a 50bp 5' UTR (genomic
        // 1-50), exon2 is the fully-coding CDS (genomic 51-62): ATG AAA CCC
        // TAA (Met Lys Pro Stop).
        let mut transcript = Transcript {
            stable_id: "ENST_FS_TEST".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_FS_TEST".into(),
                symbol: Some("FS-TEST".into()),
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "1".into(),
                start: 1,
                end: 62,
                strand: fastvep_core::Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "1".into(),
            start: 1,
            end: 62,
            strand: fastvep_core::Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "ENSE_FS_TEST_1".into(),
                    start: 1,
                    end: 50,
                    strand: fastvep_core::Strand::Forward,
                    phase: -1,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "ENSE_FS_TEST_2".into(),
                    start: 51,
                    end: 62,
                    strand: fastvep_core::Strand::Forward,
                    phase: 0,
                    end_phase: -1,
                    rank: 2,
                },
            ],
            translation: Some(Translation {
                stable_id: "ENSP_FS_TEST".into(),
                genomic_start: 51,
                genomic_end: 62,
                start_exon_rank: 2,
                start_exon_offset: 0,
                end_exon_rank: 2,
                end_exon_offset: 11,
            }),
            cdna_coding_start: Some(51),
            cdna_coding_end: Some(62),
            coding_region_start: Some(51),
            coding_region_end: Some(62),
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: Some("ENSP_FS_TEST".into()),
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
            reference_peptide: None,
        };

        transcript
            .build_sequences(|_chrom, start, _end| {
                if start == 1 {
                    Ok(b"N".repeat(50))
                } else {
                    Ok(b"ATGAAACCCTAA".to_vec())
                }
            })
            .expect("build_sequences should succeed for a well-formed test transcript");
        assert_eq!(
            transcript.translateable_seq.as_deref(),
            Some("ATGAAACCCTAA")
        );
        assert_eq!(transcript.spliced_seq.as_deref().map(|s| s.len()), Some(62));

        // Simulate `spliced_seq` becoming shorter than `cdna_coding_start`
        // implies (e.g. a corrupted/truncated cache reload of just this
        // field) -- `cdna_coding_start` (51) is left untouched, so
        // coding_start_idx (50) now exceeds the truncated spliced_seq's
        // length (20).
        transcript.spliced_seq = transcript.spliced_seq.map(|s| s[..20].to_string());

        let mut ctx = empty_context();
        ctx.transcript_provider = IndexedTranscriptProvider::new(vec![transcript]);

        // 1bp deletion of the third "A" at genomic position 56 (VCF-anchored
        // at 55): this is a frameshift but does not remove the start codon.
        let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                   1\t55\t.\tAA\tA\t.\tPASS\t.\n";

        // Must not panic despite the truncated spliced_seq.
        let results = ctx
            .annotate_vcf_text_with_acmg(vcf, false, None)
            .expect("annotation should succeed even with a truncated spliced_seq");
        assert_eq!(results.len(), 1);

        // Confirm this actually exercised the guarded path: the variant must
        // still be classified as a frameshift with real amino acids (proving
        // the truncated spliced_seq didn't derail classification, since that
        // uses translateable_seq/coding_region_start/coding_region_end
        // instead), but with no hgvsp emitted (proving the slicing guard
        // skipped it rather than panicking).
        let tc = &results[0]["transcript_consequences"][0];
        let terms = tc["consequence_terms"]
            .as_array()
            .expect("consequence_terms should be present");
        assert!(
            terms
                .iter()
                .any(|t| t.as_str() == Some("frameshift_variant")),
            "expected frameshift_variant, got: {:?}",
            terms
        );
        assert!(
            tc.get("amino_acids").is_some(),
            "amino acids should still be computed from the intact translateable_seq: {:?}",
            tc
        );
        assert!(
            tc.get("hgvsp").is_none(),
            "hgvsp should be omitted (skipped), not present, when spliced_seq is \
             shorter than coding_start_idx implies: {:?}",
            tc
        );
    }

    #[test]
    fn hgvsp_reconstructs_a_frameshift_when_shifted_cds_coordinates_are_complete() {
        use fastvep_core::{GenomicPosition, Impact, Strand};
        use fastvep_genome::{Exon, Gene, Transcript, Translation};

        let transcript = Transcript {
            stable_id: "ENST_BOUNDARY_TEST".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_BOUNDARY_TEST".into(),
                symbol: None,
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "1".into(),
                start: 1,
                end: 18,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "1".into(),
            start: 1,
            end: 18,
            strand: Strand::Forward,
            exons: vec![Exon {
                stable_id: "ENSE_BOUNDARY_TEST".into(),
                start: 1,
                end: 18,
                strand: Strand::Forward,
                phase: 0,
                end_phase: -1,
                rank: 1,
            }],
            translation: Some(Translation {
                stable_id: "ENSP_BOUNDARY_TEST".into(),
                genomic_start: 1,
                genomic_end: 12,
                start_exon_rank: 1,
                start_exon_offset: 0,
                end_exon_rank: 1,
                end_exon_offset: 11,
            }),
            cdna_coding_start: Some(1),
            cdna_coding_end: Some(12),
            coding_region_start: Some(1),
            coding_region_end: Some(12),
            spliced_seq: Some("ATGGAAGAATAACCCTAA".into()),
            translateable_seq: Some("ATGGAAGAATAA".into()),
            peptide: Some("MEE*".into()),
            canonical: false,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: Some("ENSP_BOUNDARY_TEST".into()),
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
            reference_peptide: None,
        };
        let allele = AlleleConsequenceResult {
            allele: Allele::Deletion,
            normalized_position: GenomicPosition::new("1", 4, 7, Strand::Forward),
            normalized_ref_allele: Allele::from_str("GAAG"),
            normalized_alt_allele: Allele::Deletion,
            consequences: vec![Consequence::CodingSequenceVariant],
            frameshift: false,
            impact: Impact::High,
            cdna_start: Some(4),
            cdna_end: None,
            cds_start: Some(4),
            cds_end: None,
            protein_start: None,
            protein_end: None,
            amino_acids: None,
            codons: None,
            exon: None,
            intron: None,
            distance: None,
        };

        assert_eq!(
            hgvsp_for_allele(
                None,
                "1",
                &transcript,
                "ENSP_BOUNDARY_TEST",
                &allele,
                &allele.normalized_ref_allele,
                Some("ENST_BOUNDARY_TEST:c.4_7del"),
            ),
            Some("ENSP_BOUNDARY_TEST:p.Glu2AsnfsTer?".into())
        );

        // VEP 115 recomputes the reading frame from shifted CDS endpoints.
        // A genomic three-base deletion can remove only two coding bases
        // after it moves across a one-base intron (P15 SLC37A4).
        let mut transcript = transcript;
        transcript.end = 19;
        transcript.gene.end = 19;
        transcript.coding_region_end = Some(13);
        transcript.exons[0].end = 7;
        let mut right = transcript.exons[0].clone();
        right.start = 9;
        right.end = 19;
        right.rank = 2;
        transcript.exons.push(right);
        let mut allele = allele;
        allele.normalized_position = GenomicPosition::new("1", 4, 6, Strand::Forward);
        allele.normalized_ref_allele = Allele::from_str("GAA");
        allele.cds_start = Some(4);
        allele.cds_end = Some(6);
        allele.protein_start = Some(2);
        allele.protein_end = Some(2);
        allele.amino_acids = Some(("E".into(), "-".into()));
        allele.consequences = vec![Consequence::InframeDeletion];
        assert_eq!(hgvsp_for_allele_with_offset(None, "1", &transcript,
            "ENSP_BOUNDARY_TEST", &allele, &allele.normalized_ref_allele,
            Some("ENST_BOUNDARY_TEST:c.7_8del"), Some(3)),
            Some("ENSP_BOUNDARY_TEST:p.Glu3IlefsTer?".into()));

        // A short intron satisfies VEP's coding preclassification even when
        // both original insertion flanks lack CDS coordinates (P15 ORAI1).
        transcript.end = 21;
        transcript.gene.end = 21;
        transcript.coding_region_end = Some(15);
        transcript.exons[0].end = 6;
        transcript.exons[1].start = 10;
        transcript.exons[1].end = 21;
        allele.normalized_position = GenomicPosition::new("1", 8, 7, Strand::Forward);
        allele.normalized_ref_allele = Allele::Deletion;
        allele.normalized_alt_allele = Allele::from_str("C");
        allele.allele = allele.normalized_alt_allele.clone();
        allele.cds_start = None;
        allele.cds_end = None;
        allele.protein_start = None;
        allele.protein_end = None;
        allele.amino_acids = None;
        allele.consequences = vec![Consequence::CodingSequenceVariant];
        assert_eq!(hgvsp_for_allele_with_offset(None, "1", &transcript,
            "ENSP_BOUNDARY_TEST", &allele, &Allele::Deletion,
            Some("ENST_BOUNDARY_TEST:c.7-1dup"), Some(2)),
            Some("ENSP_BOUNDARY_TEST:p.Glu3ArgfsTer?".into()));

        // Shifting into an exon cannot grant this exception to a long intron.
        transcript.exons[1].start = 100;
        transcript.exons[1].end = 111;
        transcript.end = 111;
        transcript.coding_region_end = Some(105);
        assert_eq!(hgvsp_for_allele_with_offset(None, "1", &transcript,
            "ENSP_BOUNDARY_TEST", &allele, &Allele::Deletion,
            Some("ENST_BOUNDARY_TEST:c.7-1dup"), Some(92)), None);
    }

    #[test]
    fn hgvsp_terminal_insertion_without_genomic_shift_keeps_flanks() {
        use fastvep_genome::{Exon, Gene, Transcript, Translation};

        let mut transcript = Transcript {
            stable_id: "ENST_INS_TEST".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_INS_TEST".into(),
                symbol: Some("INSTEST".into()),
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "1".into(),
                start: 1,
                end: 62,
                strand: fastvep_core::Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "1".into(),
            start: 1,
            end: 62,
            strand: fastvep_core::Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "ENSE_INS_TEST_1".into(),
                    start: 1,
                    end: 50,
                    strand: fastvep_core::Strand::Forward,
                    phase: -1,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "ENSE_INS_TEST_2".into(),
                    start: 51,
                    end: 62,
                    strand: fastvep_core::Strand::Forward,
                    phase: 0,
                    end_phase: -1,
                    rank: 2,
                },
            ],
            translation: Some(Translation {
                stable_id: "ENSP_INS_TEST".into(),
                genomic_start: 51,
                genomic_end: 62,
                start_exon_rank: 2,
                start_exon_offset: 0,
                end_exon_rank: 2,
                end_exon_offset: 11,
            }),
            cdna_coding_start: Some(51),
            cdna_coding_end: Some(62),
            coding_region_start: Some(51),
            coding_region_end: Some(62),
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: Some("ENSP_INS_TEST".into()),
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
            reference_peptide: None,
        };
        transcript
            .build_sequences(|_chrom, start, _end| {
                if start == 1 {
                    Ok(b"N".repeat(50))
                } else {
                    Ok(b"ATGTGGCGGTAA".to_vec())
                }
            })
            .unwrap();

        let mut ctx = empty_context();
        ctx.transcript_provider = IndexedTranscriptProvider::new(vec![transcript]);
        let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                   1\t56\t.\tG\tGCGG\t.\tPASS\t.\n";
        let results = ctx.annotate_vcf_text(vcf, false).unwrap();
        let consequence = &results[0]["transcript_consequences"][0];

        assert_eq!(consequence["amino_acids"], serde_json::json!("-/R"));
        assert_eq!(
            consequence["hgvsp"],
            serde_json::json!("ENSP_INS_TEST:p.Trp2_Arg3insArg")
        );
    }
}
