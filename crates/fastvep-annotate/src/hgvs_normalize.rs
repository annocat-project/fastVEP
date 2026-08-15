//! HGVS normalization utilities: insertion-to-duplication conversion and 3' shifting.
//!
//! These functions post-process HGVSc strings to convert insertion notation (`ins`)
//! to duplication notation (`dup`) when the inserted bases match the adjacent reference,
//! and to 3'-shift intronic variants per HGVS conventions.

use fastvep_cache::providers::SequenceProvider;
use fastvep_core::{Allele, Strand};
use fastvep_genome::Transcript;

/// Transcript coordinates for a deletion that becomes fully exonic after
/// applying the HGVS 3' rule across an exon boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ShiftedExonicDeletion {
    pub cdna_start: u64,
    pub cdna_end: u64,
    pub cds_start: Option<u64>,
    pub cds_end: Option<u64>,
}

/// Shift a genomic deletion in transcript 3' direction and return transcript
/// coordinates when the equivalent description becomes fully exonic.
///
/// This covers splice-boundary deletions such as a deletion that begins in an
/// intron but can be represented completely in the following exon. It changes
/// HGVS representation only; the original genomic allele remains canonical.
pub fn shifted_exonic_deletion(
    seq_provider: &dyn SequenceProvider,
    chrom: &str,
    start: u64,
    end: u64,
    ref_allele: &Allele,
    alt_allele: &Allele,
    transcript: &Transcript,
) -> Option<ShiftedExonicDeletion> {
    if !matches!(
        (ref_allele, alt_allele),
        (Allele::Sequence(bases), Allele::Deletion) if !bases.is_empty()
    ) {
        return None;
    }

    // Fully exonic deletions already use sequence-aware cDNA normalization.
    if transcript.genomic_to_cdna(start).is_some() && transcript.genomic_to_cdna(end).is_some() {
        return None;
    }

    let (shifted_start, shifted_end) = three_prime_shift_intronic(
        seq_provider,
        chrom,
        start,
        end,
        ref_allele,
        alt_allele,
        transcript.strand,
        transcript.start,
        transcript.end,
    );
    if (shifted_start, shifted_end) == (start, end) {
        return None;
    }

    let first = transcript.genomic_to_cdna(shifted_start)?;
    let last = transcript.genomic_to_cdna(shifted_end)?;
    let (cdna_start, cdna_end) = (first.min(last), first.max(last));
    Some(ShiftedExonicDeletion {
        cdna_start,
        cdna_end,
        cds_start: transcript.cdna_to_cds(cdna_start),
        cds_end: transcript.cdna_to_cds(cdna_end),
    })
}

/// Return the minimal protein replacement implied by an in-frame coding
/// deletion after transcript-level HGVS shifting.
pub fn shifted_inframe_deletion_protein_change(
    transcript: &Transcript,
    cds_start: u64,
    cds_end: u64,
) -> Option<(u64, String, String)> {
    let (cds_start, cds_end) = (cds_start.min(cds_end), cds_start.max(cds_end));
    let deleted_len = cds_end.checked_sub(cds_start)?.checked_add(1)? as usize;
    if deleted_len == 0 || deleted_len % 3 != 0 {
        return None;
    }

    let sequence = transcript.translateable_seq.as_deref()?.as_bytes();
    let deletion_start = cds_start.checked_sub(1)? as usize;
    let deletion_end = deletion_start.checked_add(deleted_len)?;
    if deletion_end > sequence.len() {
        return None;
    }

    let codon_start = deletion_start / 3 * 3;
    let codon_end = deletion_end.checked_add(2)? / 3 * 3;
    if codon_end > sequence.len() {
        return None;
    }

    let mut alternate = sequence.to_vec();
    alternate.drain(deletion_start..deletion_end);
    let alternate_end = codon_end.checked_sub(deleted_len)?;
    let table = if fastvep_genome::is_mitochondrial(&transcript.chromosome) {
        fastvep_genome::mitochondrial_codon_table()
    } else {
        fastvep_genome::CodonTable::standard()
    };
    let translate = |bases: &[u8]| -> String {
        bases
            .chunks_exact(3)
            .map(|codon| table.translate(&[codon[0], codon[1], codon[2]]) as char)
            .collect()
    };
    let reference_aas = translate(&sequence[codon_start..codon_end]);
    let alternate_aas = translate(&alternate[codon_start..alternate_end]);
    let protein_start = (codon_start / 3 + 1) as u64;
    Some((
        protein_start,
        reference_aas,
        if alternate_aas.is_empty() {
            "-".to_string()
        } else {
            alternate_aas
        },
    ))
}

/// Generate non-coding exonic HGVSc, reading adjacent reference bases only
/// when an older transcript cache does not contain the spliced sequence.
pub fn hgvsc_noncoding_exonic(
    seq_provider: Option<&dyn SequenceProvider>,
    transcript: &Transcript,
    transcript_id: &str,
    cdna_start: u64,
    cdna_end: u64,
    ref_allele: &Allele,
    alt_allele: &Allele,
) -> Option<String> {
    if let Some(seq) = transcript.spliced_seq.as_deref() {
        return fastvep_hgvs::hgvsc_noncoding_with_seq(
            transcript_id,
            cdna_start,
            cdna_end,
            ref_allele,
            alt_allele,
            Some(seq),
        );
    }

    let (mut start, mut end) = if cdna_start <= cdna_end {
        (cdna_start, cdna_end)
    } else {
        (cdna_end, cdna_start)
    };
    if matches!((ref_allele, alt_allele), (Allele::Sequence(bases), Allele::Deletion) if !bases.is_empty())
    {
        if let Some(provider) = seq_provider {
            while end < transcript.cdna_length() {
                let Some(first) = transcript_base(provider, transcript, start) else {
                    break;
                };
                let Some(next) = transcript_base(provider, transcript, end + 1) else {
                    break;
                };
                if first != next {
                    break;
                }
                start += 1;
                end += 1;
            }
        }
    }

    fastvep_hgvs::hgvsc_noncoding(transcript_id, start, end, ref_allele, alt_allele)
}

fn transcript_base(
    seq_provider: &dyn SequenceProvider,
    transcript: &Transcript,
    cdna_pos: u64,
) -> Option<u8> {
    let genomic_pos = transcript.cdna_to_genomic(cdna_pos)?;
    let base = *seq_provider
        .fetch_sequence(&transcript.chromosome, genomic_pos, genomic_pos)
        .ok()?
        .first()?;
    Some(match transcript.strand {
        Strand::Forward => base.to_ascii_uppercase(),
        Strand::Reverse => fastvep_genome::codon::reverse_complement(&[base])[0],
    })
}

/// Return the adjacent reference range duplicated by an intronic insertion.
/// If both sides match, prefer the side furthest 3' in transcript direction.
pub fn intronic_duplication_range(
    seq_provider: &dyn SequenceProvider,
    chrom: &str,
    insertion_start: u64,
    insertion_end: u64,
    inserted_bases: &[u8],
    transcript: &Transcript,
) -> Option<((u64, i64), (u64, i64))> {
    if inserted_bases.is_empty() {
        return None;
    }

    let len = inserted_bases.len() as u64;
    let before = insertion_end
        .checked_sub(len - 1)
        .map(|start| (start, insertion_end));
    let after = insertion_start
        .checked_add(len - 1)
        .map(|end| (insertion_start, end));
    let candidates = match transcript.strand {
        Strand::Forward => [after, before],
        Strand::Reverse => [before, after],
    };

    for (start, end) in candidates.into_iter().flatten() {
        let Some(bounds) = transcript.intron_bounds_at(start) else {
            continue;
        };
        if end > bounds.1 || transcript.intron_bounds_at(end) != Some(bounds) {
            continue;
        }
        let Ok(reference) = seq_provider.fetch_sequence_slice(chrom, start, end) else {
            continue;
        };
        if reference.len() != inserted_bases.len()
            || !reference
                .iter()
                .zip(inserted_bases)
                .all(|(a, b)| a.eq_ignore_ascii_case(b))
        {
            continue;
        }

        let (first, last) = match transcript.strand {
            Strand::Forward => (start, end),
            Strand::Reverse => (end, start),
        };
        return Some((
            transcript.genomic_to_intronic_cdna(first)?,
            transcript.genomic_to_intronic_cdna(last)?,
        ));
    }

    None
}

/// Convert an intronic insertion HGVSc to duplication notation (coding transcript).
pub fn convert_ins_to_dup(
    hgvsc: &str,
    intron_offset: i64,
    ins_len: u64,
    nearest_exon_cdna_pos: u64,
    coding_start: u64,
    coding_end: Option<u64>,
) -> Option<String> {
    let prefix_end = hgvsc
        .find(":c.")
        .map(|i| i + 3)
        .or_else(|| hgvsc.find(":n.").map(|i| i + 3))?;
    let prefix = &hgvsc[..prefix_end];

    let build_pos = |cdna: u64, off: i64| -> String {
        let raw = cdna as i64 - coding_start as i64 + 1;
        let cp = if raw <= 0 { raw - 1 } else { raw };
        if cp < 0 {
            if off > 0 {
                format!("{}+{}", cp, off)
            } else {
                format!("{}{}", cp, off)
            }
        } else if coding_end.is_some() && cdna > coding_end.unwrap() {
            let u = cdna - coding_end.unwrap();
            if off > 0 {
                format!("*{}+{}", u, off)
            } else {
                format!("*{}{}", u, off)
            }
        } else if off > 0 {
            format!("{}+{}", cp, off)
        } else {
            format!("{}{}", cp, off)
        }
    };

    if ins_len == 1 {
        let pos = build_pos(nearest_exon_cdna_pos, intron_offset);
        Some(format!("{}{}dup", prefix, pos))
    } else {
        let start_offset = intron_offset - ins_len as i64 + 1;
        let start_pos = build_pos(nearest_exon_cdna_pos, start_offset);
        let end_pos = build_pos(nearest_exon_cdna_pos, intron_offset);
        Some(format!("{}{}_{}dup", prefix, start_pos, end_pos))
    }
}

/// Convert intronic insertion to dup notation with explicit start/end offsets (coding).
pub fn convert_ins_to_dup_range(
    hgvsc: &str,
    start_cdna_pos: u64,
    start_offset: i64,
    end_cdna_pos: u64,
    end_offset: i64,
    coding_start: u64,
    coding_end: Option<u64>,
) -> Option<String> {
    let prefix_end = hgvsc
        .find(":c.")
        .map(|i| i + 3)
        .or_else(|| hgvsc.find(":n.").map(|i| i + 3))?;
    let prefix = &hgvsc[..prefix_end];

    let build_pos = |cdna: u64, off: i64| -> String {
        let raw = cdna as i64 - coding_start as i64 + 1;
        let cp = if raw <= 0 { raw - 1 } else { raw };
        if cp < 0 {
            if off > 0 {
                format!("{}+{}", cp, off)
            } else {
                format!("{}{}", cp, off)
            }
        } else if coding_end.is_some() && cdna > coding_end.unwrap() {
            let u = cdna - coding_end.unwrap();
            if off > 0 {
                format!("*{}+{}", u, off)
            } else {
                format!("*{}{}", u, off)
            }
        } else if off > 0 {
            format!("{}+{}", cp, off)
        } else {
            format!("{}{}", cp, off)
        }
    };

    if start_cdna_pos == end_cdna_pos && start_offset == end_offset {
        let pos = build_pos(start_cdna_pos, start_offset);
        Some(format!("{}{}dup", prefix, pos))
    } else {
        let start_pos = build_pos(start_cdna_pos, start_offset);
        let end_pos = build_pos(end_cdna_pos, end_offset);
        Some(format!("{}{}_{}dup", prefix, start_pos, end_pos))
    }
}

/// Convert intronic insertion to dup notation with explicit start/end offsets (non-coding).
pub fn convert_ins_to_dup_range_noncoding(
    hgvsc: &str,
    start_cdna_pos: u64,
    start_offset: i64,
    end_cdna_pos: u64,
    end_offset: i64,
) -> Option<String> {
    let prefix_end = hgvsc
        .find(":n.")
        .map(|i| i + 3)
        .or_else(|| hgvsc.find(":c.").map(|i| i + 3))?;
    let prefix = &hgvsc[..prefix_end];

    let build_pos = |cdna: u64, off: i64| -> String {
        if off > 0 {
            format!("{}+{}", cdna, off)
        } else {
            format!("{}{}", cdna, off)
        }
    };

    if start_cdna_pos == end_cdna_pos && start_offset == end_offset {
        let pos = build_pos(start_cdna_pos, start_offset);
        Some(format!("{}{}dup", prefix, pos))
    } else {
        let start_pos = build_pos(start_cdna_pos, start_offset);
        let end_pos = build_pos(end_cdna_pos, end_offset);
        Some(format!("{}{}_{}dup", prefix, start_pos, end_pos))
    }
}

/// Convert an intronic insertion HGVSc to duplication notation (non-coding transcript).
pub fn convert_ins_to_dup_noncoding(
    hgvsc: &str,
    intron_offset: i64,
    ins_len: u64,
    nearest_exon_cdna_pos: u64,
) -> Option<String> {
    let prefix_end = hgvsc
        .find(":n.")
        .map(|i| i + 3)
        .or_else(|| hgvsc.find(":c.").map(|i| i + 3))?;
    let prefix = &hgvsc[..prefix_end];

    let build_pos = |off: i64| -> String {
        if off > 0 {
            format!("{}+{}", nearest_exon_cdna_pos, off)
        } else {
            format!("{}{}", nearest_exon_cdna_pos, off)
        }
    };

    if ins_len == 1 {
        let pos = build_pos(intron_offset);
        Some(format!("{}{}dup", prefix, pos))
    } else {
        let start_offset = intron_offset - ins_len as i64 + 1;
        let start_pos = build_pos(start_offset);
        let end_pos = build_pos(intron_offset);
        Some(format!("{}{}_{}dup", prefix, start_pos, end_pos))
    }
}

/// 3' shift an intronic indel along the transcript direction.
///
/// HGVS requires variants to be described at the most 3' position.
/// For intronic deletions and insertions/dups in repetitive regions,
/// the position must be shifted toward the 3' end of the transcript.
///
/// Returns the shifted genomic start and end positions.
pub fn three_prime_shift_intronic(
    seq_provider: &dyn SequenceProvider,
    chrom: &str,
    start: u64,
    end: u64,
    ref_allele: &fastvep_core::Allele,
    alt_allele: &fastvep_core::Allele,
    strand: fastvep_core::Strand,
    intron_genomic_start: u64,
    intron_genomic_end: u64,
) -> (u64, u64) {
    use fastvep_core::Allele;

    match (ref_allele, alt_allele) {
        // Deletion: shift the deleted bases toward 3' end
        (Allele::Sequence(ref_bases), Allele::Deletion) if !ref_bases.is_empty() => {
            let mut s = start;
            let mut e = end;

            match strand {
                fastvep_core::Strand::Forward => loop {
                    let next_pos = e + 1;
                    if next_pos > intron_genomic_end {
                        break;
                    }
                    let next_base = match seq_provider.fetch_sequence(chrom, next_pos, next_pos) {
                        Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                        _ => break,
                    };
                    let first_base = match seq_provider.fetch_sequence(chrom, s, s) {
                        Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                        _ => break,
                    };
                    if next_base == first_base {
                        s += 1;
                        e += 1;
                    } else {
                        break;
                    }
                },
                fastvep_core::Strand::Reverse => loop {
                    if s == 0 || s - 1 < intron_genomic_start {
                        break;
                    }
                    let prev_pos = s - 1;
                    let prev_base = match seq_provider.fetch_sequence(chrom, prev_pos, prev_pos) {
                        Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                        _ => break,
                    };
                    let last_base = match seq_provider.fetch_sequence(chrom, e, e) {
                        Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                        _ => break,
                    };
                    if prev_base == last_base {
                        s -= 1;
                        e -= 1;
                    } else {
                        break;
                    }
                },
            }
            (s, e)
        }
        // Insertion/dup: shift toward 3' end using the actual inserted bases
        (Allele::Deletion, Allele::Sequence(ins_bases)) if !ins_bases.is_empty() => {
            let ins_len = ins_bases.len();
            let mut pos = start;
            let genomic_ins: Vec<u8> = ins_bases.iter().map(|b| b.to_ascii_uppercase()).collect();

            match strand {
                fastvep_core::Strand::Forward => {
                    let mut shift_count = 0u64;
                    loop {
                        if pos > intron_genomic_end {
                            break;
                        }
                        let check_base = match seq_provider.fetch_sequence(chrom, pos, pos) {
                            Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                            _ => break,
                        };
                        let idx = (shift_count as usize) % ins_len;
                        if check_base == genomic_ins[idx] {
                            pos += 1;
                            shift_count += 1;
                        } else {
                            break;
                        }
                    }
                }
                fastvep_core::Strand::Reverse => {
                    let mut shift_count = 0u64;
                    loop {
                        if pos == 0 || pos - 1 < intron_genomic_start {
                            break;
                        }
                        let check_pos = pos - 1;
                        let check_base =
                            match seq_provider.fetch_sequence(chrom, check_pos, check_pos) {
                                Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                                _ => break,
                            };
                        let idx = ins_len - 1 - (shift_count as usize % ins_len);
                        if check_base == genomic_ins[idx] {
                            pos -= 1;
                            shift_count += 1;
                        } else {
                            break;
                        }
                    }
                }
            }
            (pos, pos.saturating_sub(1))
        }
        _ => (start, end),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::{anyhow, Result};
    use fastvep_genome::{Exon, Gene, Translation};

    struct TestSequence(Vec<u8>);

    impl SequenceProvider for TestSequence {
        fn fetch_sequence(&self, _chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
            let first = start.checked_sub(1).ok_or_else(|| anyhow!("position 0"))? as usize;
            let last = end as usize;
            self.0
                .get(first..last)
                .map(|bases| bases.to_vec())
                .ok_or_else(|| anyhow!("out of range"))
        }
    }

    fn coding_transcript() -> Transcript {
        Transcript {
            stable_id: "ENST_TEST".into(),
            version: Some(1),
            gene: Gene {
                stable_id: "ENSG_TEST".into(),
                symbol: Some("TEST".into()),
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start: 1,
                end: 20,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start: 1,
            end: 20,
            strand: Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "E1".into(),
                    start: 1,
                    end: 2,
                    strand: Strand::Forward,
                    phase: -1,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "E2".into(),
                    start: 10,
                    end: 20,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: -1,
                    rank: 2,
                },
            ],
            translation: Some(Translation {
                stable_id: "ENSP_TEST".into(),
                genomic_start: 1,
                genomic_end: 20,
                start_exon_rank: 1,
                start_exon_offset: 0,
                end_exon_rank: 2,
                end_exon_offset: 10,
            }),
            cdna_coding_start: Some(1),
            cdna_coding_end: Some(13),
            coding_region_start: Some(1),
            coding_region_end: Some(20),
            spliced_seq: Some("CCGAGGTAAAAAA".to_string()),
            translateable_seq: Some("ATGGAAGAAGAA".to_string()),
            peptide: Some("MEEE".to_string()),
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: Some(1),
            appris: None,
            ccds: None,
            protein_id: Some("ENSP_TEST".to_string()),
            protein_version: Some(1),
            swissprot: Vec::new(),
            trembl: Vec::new(),
            uniparc: Vec::new(),
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: Vec::new(),
            codon_table_start_phase: 0,
        }
    }

    #[test]
    fn shifts_splice_boundary_deletion_into_exon() {
        // Positions 8..10 delete AGG. The repeating AGG permits a three-base
        // 3' shift to positions 11..13, which are fully exonic.
        let provider = TestSequence(b"CCTTTTTAGGAGGTAAAAAA".to_vec());
        let transcript = coding_transcript();
        let shifted = shifted_exonic_deletion(
            &provider,
            "chr1",
            8,
            10,
            &Allele::Sequence(b"AGG".to_vec()),
            &Allele::Deletion,
            &transcript,
        )
        .unwrap();

        assert_eq!(shifted.cdna_start, 4);
        assert_eq!(shifted.cdna_end, 6);
        assert_eq!(shifted.cds_start, Some(4));
        assert_eq!(shifted.cds_end, Some(6));
    }

    #[test]
    fn derives_inframe_protein_replacement_after_shift() {
        let transcript = coding_transcript();
        let change = shifted_inframe_deletion_protein_change(&transcript, 4, 6).unwrap();
        assert_eq!(change, (2, "E".to_string(), "-".to_string()));
    }
}
