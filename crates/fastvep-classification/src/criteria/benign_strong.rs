use crate::config::AcmgConfig;
use crate::sa_extract::ClassificationInput;
use crate::types::{EvidenceCriterion, EvidenceDirection, EvidenceStrength};

/// Evaluate all benign strong criteria: BS1, BS2, BS3, BS4.
pub fn evaluate_all(input: &ClassificationInput, config: &AcmgConfig) -> Vec<EvidenceCriterion> {
    vec![
        evaluate_bs1(input, config),
        evaluate_bs2(input, config),
        evaluate_bs3(input, config),
        evaluate_bs4(input, config),
    ]
}

/// BS1: Allele frequency is greater than expected for disorder.
fn evaluate_bs1(input: &ClassificationInput, config: &AcmgConfig) -> EvidenceCriterion {
    let threshold = config.effective_bs1_threshold(input.gene_symbol.as_deref());

    let mut details = serde_json::Map::new();
    details.insert("af_threshold".into(), serde_json::json!(threshold));

    let (met, summary) = if let Some(ref gnomad) = input.gnomad {
        // ClinGen SVI gnomAD v4 guidance (March 2024): require minimum AN
        // before BS1 fires, same as BA1. Treat missing AN as NotEvaluated
        // (the SVI guidance is a requirement, not an opt-in).
        match gnomad.all_an {
            Some(an) if an >= config.min_an_for_frequency_criteria => {}
            other => {
                details.insert(
                    "min_an_for_frequency_criteria".into(),
                    serde_json::json!(config.min_an_for_frequency_criteria),
                );
                let summary = match other {
                    Some(an) => {
                        details.insert("an_below_minimum".into(), serde_json::json!(an));
                        format!(
                            "BS1 not evaluated: gnomAD AN={} below minimum {} (gnomAD v4 guidance)",
                            an, config.min_an_for_frequency_criteria
                        )
                    }
                    None => {
                        details.insert("an_missing".into(), serde_json::json!(true));
                        format!(
                            "BS1 not evaluated: gnomAD AN unavailable; minimum {} required (gnomAD v4 guidance)",
                            config.min_an_for_frequency_criteria
                        )
                    }
                };
                return EvidenceCriterion {
                    code: "BS1".to_string(),
                    direction: EvidenceDirection::Benign,
                    strength: EvidenceStrength::Strong,
                    default_strength: EvidenceStrength::Strong,
                    met: false,
                    evaluated: false,
                    summary,
                    details: serde_json::Value::Object(details),
                };
            }
        }
        // ClinGen SVI guidance applies BS1 against the **max-population
        // AF** (mirroring BA1), not the cohort-wide allAf. Using the
        // cohort AF would let a 5% variant in a single subpopulation slip
        // under a 1 % BS1 threshold whenever the global cohort happens
        // to dilute it. METHODS.md aligns: "max population AF" for both
        // BA1 and BS1.
        let max_pop_af = gnomad.max_pop_af().unwrap_or(0.0);
        let cohort_af = gnomad.all_af.unwrap_or(0.0);
        details.insert("gnomad_allAf".into(), serde_json::json!(cohort_af));
        details.insert("gnomad_max_pop_af".into(), serde_json::json!(max_pop_af));

        // BS1 should not fire if BA1 would fire (BA1 takes precedence)
        if max_pop_af > config.ba1_af_threshold {
            (
                false,
                format!(
                    "BA1 takes precedence (max pop AF={:.4} > BA1 threshold {:.2})",
                    max_pop_af, config.ba1_af_threshold
                ),
            )
        } else if max_pop_af > threshold {
            (
                true,
                format!(
                    "Max-pop AF ({:.6}) exceeds expected for disorder (threshold={:.4})",
                    max_pop_af, threshold
                ),
            )
        } else {
            (
                false,
                format!(
                    "Max-pop AF ({:.6}) within expected range (threshold={:.4})",
                    max_pop_af, threshold
                ),
            )
        }
    } else {
        (false, "No gnomAD data available".to_string())
    };

    EvidenceCriterion {
        code: "BS1".to_string(),
        direction: EvidenceDirection::Benign,
        strength: EvidenceStrength::Strong,
        default_strength: EvidenceStrength::Strong,
        met,
        evaluated: input.gnomad.is_some(),
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// BS2: Observed in a healthy adult individual for a recessive (homozygous),
/// dominant (heterozygous), or X-linked (hemizygous) disorder with full
/// penetrance expected at an early age (Richards 2015).
///
/// - **Recessive** (or unknown inheritance): require ≥1 homozygote in gnomAD.
/// - **Dominant**: require AC ≥ `bs2_ad_min_ac` (default 5) — Richards 2015
///   says "observed in unaffected adult", which is not the same as a single
///   carrier of a novel allele in a 100K cohort. Singletons / doubletons
///   are sequencing-noise plausibility, not evidence the variant is
///   tolerated. ClinGen VCEPs commonly use AC ≥ 5 (Hereditary Cancer
///   VCEP, Lynch Syndrome curation guide).
///
/// Inheritance is inferred from the disease-gene `.oga` (ClinGen GDV
/// preferred, OMIM accepted as legacy).
fn evaluate_bs2(input: &ClassificationInput, config: &AcmgConfig) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();

    let is_dominant = input
        .omim
        .as_ref()
        .map_or(false, |o| o.has_dominant_inheritance());
    let is_recessive = input
        .omim
        .as_ref()
        .map_or(false, |o| o.has_recessive_inheritance());
    details.insert("omim_dominant".into(), serde_json::json!(is_dominant));
    details.insert("omim_recessive".into(), serde_json::json!(is_recessive));

    let (met, evaluated, summary) = if let Some(ref gnomad) = input.gnomad {
        let hc = gnomad.all_hc.unwrap_or(0);
        let an = gnomad.all_an.unwrap_or(0);
        let ac = gnomad.all_ac.unwrap_or(0);
        details.insert("gnomad_allHc".into(), serde_json::json!(hc));
        details.insert("gnomad_allAn".into(), serde_json::json!(an));
        details.insert("gnomad_allAc".into(), serde_json::json!(ac));

        if is_dominant && !is_recessive && ac >= config.bs2_ad_min_ac {
            // For AD-only genes: ≥`bs2_ad_min_ac` heterozygote observations
            // in healthy adults (gnomAD). Singletons / doubletons of a
            // novel allele are not BS2 evidence.
            (
                true,
                true,
                format!(
                    "Observed in gnomAD as ≥{} unaffected heterozygotes (AC={}) for autosomal-dominant disorder",
                    config.bs2_ad_min_ac, ac
                ),
            )
        } else if hc > 0 {
            // Recessive / X-linked / unknown inheritance: ≥1 homozygote is
            // BS2 evidence regardless of inheritance label (a homozygous
            // healthy adult disproves recessive lethality and challenges
            // dominant haploinsufficiency-or-not).
            (
                true,
                true,
                format!(
                    "Observed as homozygous in gnomAD ({} homozygotes), suggesting tolerated in healthy adults",
                    hc
                ),
            )
        } else if is_dominant && !is_recessive {
            (
                false,
                true,
                format!(
                    "AC={} below BS2 threshold ({} required for AD)",
                    ac, config.bs2_ad_min_ac
                ),
            )
        } else {
            (false, true, "No homozygotes observed in gnomAD".to_string())
        }
    } else {
        (false, false, "No gnomAD data available".to_string())
    };

    EvidenceCriterion {
        code: "BS2".to_string(),
        direction: EvidenceDirection::Benign,
        strength: EvidenceStrength::Strong,
        default_strength: EvidenceStrength::Strong,
        met,
        evaluated,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// BS3: Well-established in vitro or in vivo functional studies show no damaging effect.
fn evaluate_bs3(_input: &ClassificationInput, _config: &AcmgConfig) -> EvidenceCriterion {
    EvidenceCriterion {
        code: "BS3".to_string(),
        direction: EvidenceDirection::Benign,
        strength: EvidenceStrength::Strong,
        default_strength: EvidenceStrength::Strong,
        met: false,
        evaluated: false,
        summary: "Requires curated functional study evidence showing no damaging effect — not automatable from variant data".to_string(),
        details: serde_json::Value::Null,
    }
}

/// BS4: Lack of segregation in affected members of a family.
fn evaluate_bs4(_input: &ClassificationInput, _config: &AcmgConfig) -> EvidenceCriterion {
    EvidenceCriterion {
        code: "BS4".to_string(),
        direction: EvidenceDirection::Benign,
        strength: EvidenceStrength::Strong,
        default_strength: EvidenceStrength::Strong,
        met: false,
        evaluated: false,
        summary:
            "Requires multi-generation pedigree with affection status to assess lack of segregation"
                .to_string(),
        details: serde_json::Value::Null,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sa_extract::GnomadData;
    use fastvep_core::Impact;

    fn make_input(gnomad: Option<GnomadData>) -> ClassificationInput {
        ClassificationInput {
            consequences: vec![],
            impact: Impact::Modifier,
            gene_symbol: Some("TEST".to_string()),
            is_canonical: true,
            amino_acids: None,
            protein_position: None,
            gnomad,
            clinvar: None,
            revel: None,
            splice_ai: None,
            dbnsfp: None,
            phylop: None,
            gerp: None,
            gene_constraints: None,
            omim: None,
            clinvar_protein: None,
            hgvs_c: None,
            predicted_nmd: None,
            protein_truncation_pct: None,
            is_last_exon: None,
            in_critical_region: None,
            alt_start_codon_distance: None,
            same_splice_position_pathogenic: None,
            in_repeat_region: None,
            at_exon_edge: None,
            intronic_offset: None,
            proband_genotype: None,
            mother_genotype: None,
            father_genotype: None,
            companion_variants: vec![],
        }
    }

    #[test]
    fn test_bs1_above_threshold() {
        let input = make_input(Some(GnomadData {
            all_af: Some(0.02),
            all_an: Some(100_000),
            ..Default::default()
        }));
        let result = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(result.met);
    }

    #[test]
    fn test_bs1_below_threshold() {
        let input = make_input(Some(GnomadData {
            all_af: Some(0.001),
            all_an: Some(100_000),
            ..Default::default()
        }));
        let result = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_bs1_ba1_takes_precedence() {
        let input = make_input(Some(GnomadData {
            all_af: Some(0.10),
            afr_af: Some(0.10),
            all_an: Some(100_000),
            ..Default::default()
        }));
        let result = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(!result.met); // BA1 would fire, so BS1 should not
    }

    #[test]
    fn test_bs1_low_an_not_evaluated() {
        // gnomAD v4 guidance: AN below 2000 → NotEvaluated, even at high AF.
        let input = make_input(Some(GnomadData {
            all_af: Some(0.02),
            all_an: Some(500),
            ..Default::default()
        }));
        let result = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(!result.met);
        assert!(!result.evaluated);
        assert!(result.summary.contains("below minimum"));
    }

    #[test]
    fn test_bs1_missing_an_not_evaluated() {
        // gnomAD record present but AN is None → NotEvaluated, never fires.
        let input = make_input(Some(GnomadData {
            all_af: Some(0.02),
            all_an: None,
            ..Default::default()
        }));
        let result = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(!result.met);
        assert!(!result.evaluated);
        assert!(result.summary.contains("AN unavailable"));
    }

    #[test]
    fn test_bs2_homozygotes_present() {
        let input = make_input(Some(GnomadData {
            all_hc: Some(5),
            ..Default::default()
        }));
        let result = evaluate_bs2(&input, &AcmgConfig::default());
        assert!(result.met);
    }

    #[test]
    fn test_bs2_no_homozygotes() {
        let input = make_input(Some(GnomadData {
            all_hc: Some(0),
            ..Default::default()
        }));
        let result = evaluate_bs2(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    fn make_input_omim(
        gnomad: Option<GnomadData>,
        omim: Option<crate::sa_extract::OmimData>,
    ) -> ClassificationInput {
        let mut i = make_input(gnomad);
        i.omim = omim;
        i
    }

    #[test]
    fn test_bs2_ad_gene_singleton_does_not_fire() {
        // AD gene + AC=1: a single heterozygote is not "observed in
        // healthy adult" per Richards 2015 BS2; default threshold AC≥5.
        use crate::sa_extract::OmimData;
        let input = make_input_omim(
            Some(GnomadData {
                all_ac: Some(1),
                all_hc: Some(0),
                all_af: Some(0.000001),
                all_an: Some(1_000_000),
                ..Default::default()
            }),
            Some(OmimData {
                mim_number: None,
                phenotypes: Some(vec!["dominant disorder".into()]),
            }),
        );
        let r = evaluate_bs2(&input, &AcmgConfig::default());
        assert!(!r.met, "AC=1 should not fire AD BS2; got {}", r.summary);
        assert!(r.summary.contains("below BS2 threshold"));
    }

    #[test]
    fn test_bs2_ad_gene_meets_threshold_fires() {
        // AD gene + AC≥5 (default `bs2_ad_min_ac`) → BS2 fires.
        use crate::sa_extract::OmimData;
        let input = make_input_omim(
            Some(GnomadData {
                all_ac: Some(7),
                all_hc: Some(0),
                all_af: Some(7e-6),
                all_an: Some(1_000_000),
                ..Default::default()
            }),
            Some(OmimData {
                mim_number: None,
                phenotypes: Some(vec!["dominant disorder".into()]),
            }),
        );
        let r = evaluate_bs2(&input, &AcmgConfig::default());
        assert!(r.met, "AC=7 ≥ default 5 should fire AD BS2");
        assert!(r.summary.contains("autosomal-dominant"));
    }

    #[test]
    fn test_bs2_ad_gene_min_ac_configurable() {
        // Config knob lets a stricter VCEP raise the threshold.
        use crate::sa_extract::OmimData;
        let mut cfg = AcmgConfig::default();
        cfg.bs2_ad_min_ac = 20;
        let input = make_input_omim(
            Some(GnomadData {
                all_ac: Some(7),
                all_hc: Some(0),
                ..Default::default()
            }),
            Some(OmimData {
                mim_number: None,
                phenotypes: Some(vec!["dominant disorder".into()]),
            }),
        );
        let r = evaluate_bs2(&input, &cfg);
        assert!(!r.met, "AC=7 < raised threshold 20 should not fire");
    }

    #[test]
    fn test_bs2_ar_gene_homozygote_fires_regardless_of_ac() {
        // Recessive: ≥1 hom is BS2 evidence even when AC is low.
        use crate::sa_extract::OmimData;
        let input = make_input_omim(
            Some(GnomadData {
                all_ac: Some(2),
                all_hc: Some(1),
                ..Default::default()
            }),
            Some(OmimData {
                mim_number: None,
                phenotypes: Some(vec!["recessive disorder".into()]),
            }),
        );
        let r = evaluate_bs2(&input, &AcmgConfig::default());
        assert!(r.met, "AR + 1 hom should fire BS2");
        assert!(r.summary.contains("homozygous"));
    }

    // ── BS1 (max-pop AF) ──

    #[test]
    fn test_bs1_uses_max_pop_af_not_cohort_af() {
        // Cohort AF below threshold but max-pop AF above: BS1 fires.
        // (Pre-fix this would have used `all_af` and missed the
        // single-population enrichment.)
        let input = make_input(Some(GnomadData {
            all_af: Some(0.001),
            // 5 % in EAS — well above default BS1 threshold of 1 %.
            eas_af: Some(0.05),
            all_an: Some(2_000_000),
            ..Default::default()
        }));
        let r = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(r.met, "max-pop AF should drive BS1, not cohort AF");
        assert!(r.summary.contains("Max-pop"));
    }
}
