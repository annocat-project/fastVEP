use fastvep_genome::codon::{aa_one_to_three, CodonTable};

/// Generate HGVSp (protein) notation.
///
/// Format: ENSP00000001:p.Arg41Lys (missense)
///         ENSP00000001:p.Arg41Ter (stop gained)
///         ENSP00000001:p.Arg41= (synonymous)
///         ENSP00000001:p.Arg41fs (frameshift)
pub fn hgvsp(
    protein_id: &str,
    protein_pos: u64,
    ref_aa: u8,
    alt_aa: u8,
    is_frameshift: bool,
) -> Option<String> {
    let prefix = format!("{}:p.", protein_id);
    let ref_aa3 = aa_one_to_three(ref_aa);

    if is_frameshift {
        return Some(format!("{}{}{}fs", prefix, ref_aa3, protein_pos));
    }

    if protein_pos == 1 && ref_aa == b'M' && alt_aa != ref_aa {
        return Some(format!("{}Met1?", prefix));
    }

    if ref_aa == alt_aa {
        // Synonymous
        return Some(format!("{}{}{}=", prefix, ref_aa3, protein_pos));
    }

    let alt_aa3 = aa_one_to_three(alt_aa);

    if alt_aa == b'*' {
        // Stop gained
        return Some(format!("{}{}{}{}", prefix, ref_aa3, protein_pos, alt_aa3));
    }

    if ref_aa == b'*' {
        // Stop lost - extension
        return Some(format!("{}{}{}ext*?", prefix, alt_aa3, protein_pos));
    }

    // Missense
    Some(format!("{}{}{}{}", prefix, ref_aa3, protein_pos, alt_aa3))
}

/// Generate HGVSp for a translated substitution spanning one or more codons.
pub fn hgvsp_substitution(
    protein_id: &str,
    protein_start: u64,
    ref_aas: &str,
    alt_aas: &str,
) -> Option<String> {
    let mut reference = ref_aas.as_bytes();
    let mut alternate = alt_aas.as_bytes();
    let mut start = protein_start;

    while !reference.is_empty() && !alternate.is_empty() && reference[0] == alternate[0] {
        reference = &reference[1..];
        alternate = &alternate[1..];
        start += 1;
    }
    while !reference.is_empty() && !alternate.is_empty() && reference.last() == alternate.last() {
        reference = &reference[..reference.len() - 1];
        alternate = &alternate[..alternate.len() - 1];
    }

    if reference.is_empty() && alternate.is_empty() {
        let aa = ref_aas.as_bytes().first().copied()?;
        return hgvsp(protein_id, protein_start, aa, aa, false);
    }
    if reference.len() == 1 && alternate.len() == 1 {
        return hgvsp(protein_id, start, reference[0], alternate[0], false);
    }

    unshifted_inframe_change(&format!("{}:p.", protein_id), start, reference, alternate)
}

fn residue_span(first_pos: u64, residues: &[u8]) -> String {
    let first = aa_one_to_three(residues[0]);
    if residues.len() == 1 {
        format!("{}{}", first, first_pos)
    } else {
        let last = aa_one_to_three(residues[residues.len() - 1]);
        format!(
            "{}{}_{}{}",
            first,
            first_pos,
            last,
            first_pos + residues.len() as u64 - 1
        )
    }
}

fn three_letter(residues: &[u8]) -> String {
    residues.iter().map(|&aa| aa_one_to_three(aa)).collect()
}

fn unshifted_inframe_change(
    prefix: &str,
    start: u64,
    reference: &[u8],
    alternate: &[u8],
) -> Option<String> {
    if reference.is_empty() {
        return None;
    }
    let range = residue_span(start, reference);
    if alternate.is_empty() {
        Some(format!("{}{}del", prefix, range))
    } else {
        Some(format!(
            "{}{}delins{}",
            prefix,
            range,
            three_letter(alternate)
        ))
    }
}

fn peptide_carries(peptide: &[u8], protein_start: u64, reference: &[u8]) -> bool {
    let Some(start) = protein_start.checked_sub(1).map(|value| value as usize) else {
        return false;
    };
    if reference.is_empty() {
        return start <= peptide.len();
    }
    let Some(end) = start.checked_add(reference.len()) else {
        return false;
    };
    peptide.get(start..end) == Some(reference)
}

/// Generate HGVSp notation for an in-frame insertion, deletion, or delins.
///
/// Insertions and deletions use the HGVS 3' rule when the transcript peptide
/// agrees with the reported reference residues. A missing, truncated, or
/// inconsistent peptide falls back to a valid unshifted description. Delins
/// changes are reduced to their minimal changed region but are not shifted.
pub fn hgvsp_inframe_indel(
    protein_id: &str,
    protein_start: u64,
    ref_aas: &str,
    alt_aas: &str,
    peptide: Option<&str>,
) -> Option<String> {
    let strip = |value: &str| {
        if value == "-" {
            Vec::new()
        } else {
            value.as_bytes().to_vec()
        }
    };
    let original_ref = strip(ref_aas);
    let original_alt = strip(alt_aas);
    let prefix = format!("{}:p.", protein_id);
    let fallback =
        || unshifted_inframe_change(&prefix, protein_start, &original_ref, &original_alt);

    // Transcript peptides include the terminator. It is not a protein residue
    // and must not participate in shifting or flank an insertion.
    let peptide = peptide
        .map(str::as_bytes)
        .map(|value| match value.iter().position(|&aa| aa == b'*') {
            Some(stop) => &value[..stop],
            None => value,
        })
        .filter(|value| peptide_carries(value, protein_start, &original_ref));

    // Reduce the caller's replacement to the minimal changed region.
    let mut reference = original_ref.clone();
    let mut alternate = original_alt.clone();
    let mut start = protein_start;
    while !reference.is_empty() && !alternate.is_empty() && reference[0] == alternate[0] {
        reference.remove(0);
        alternate.remove(0);
        start += 1;
    }
    while !reference.is_empty() && !alternate.is_empty() && reference.last() == alternate.last() {
        reference.pop();
        alternate.pop();
    }

    if reference.is_empty() && alternate.is_empty() {
        return None;
    }

    if reference.is_empty() {
        let (Some(peptide), Some(mut at)) = (peptide, start.checked_sub(1).map(|v| v as usize))
        else {
            return fallback();
        };
        let mut inserted = alternate;
        while at < peptide.len() && peptide[at] == inserted[0] {
            inserted.rotate_left(1);
            at += 1;
        }

        let preceding = at
            .checked_sub(inserted.len())
            .and_then(|begin| peptide.get(begin..at));
        if preceding == Some(inserted.as_slice()) {
            let dup_start = (at - inserted.len() + 1) as u64;
            return Some(format!(
                "{}{}dup",
                prefix,
                residue_span(dup_start, &inserted)
            ));
        }

        match (
            at.checked_sub(1).and_then(|index| peptide.get(index)),
            peptide.get(at),
        ) {
            (Some(&before), Some(&after)) => Some(format!(
                "{}{}{}_{}{}ins{}",
                prefix,
                aa_one_to_three(before),
                at,
                aa_one_to_three(after),
                at + 1,
                three_letter(&inserted)
            )),
            _ => fallback(),
        }
    } else if alternate.is_empty() {
        let mut at = start;
        let mut residues = reference;
        if let Some(peptide) = peptide {
            let len = residues.len();
            loop {
                let first = at
                    .checked_sub(1)
                    .and_then(|index| peptide.get(index as usize));
                let following = peptide.get(at as usize + len - 1);
                match (first, following) {
                    (Some(left), Some(right)) if left == right => at += 1,
                    _ => break,
                }
            }
            match at
                .checked_sub(1)
                .and_then(|begin| peptide.get(begin as usize..begin as usize + len))
            {
                Some(block) => residues = block.to_vec(),
                None => return fallback(),
            }
        }
        Some(format!("{}{}del", prefix, residue_span(at, &residues)))
    } else {
        Some(format!(
            "{}{}delins{}",
            prefix,
            residue_span(start, &reference),
            three_letter(&alternate)
        ))
    }
}

/// Generate HGVSp notation for a frameshift variant.
///
/// Scans the frameshifted sequence to find the first changed amino acid and
/// the position of the new stop codon.
///
/// Format: ENSP00000001:p.Ala498ProfsTer28
///   - Ala498 = first amino acid that changes (ref)
///   - Pro = new amino acid at that position
///   - Ter28 = new stop codon 28 positions downstream
///
/// `codon_table` lets the caller select the genetic code to translate with —
/// pass the vertebrate mitochondrial table (NCBI table 2) for MT transcripts
/// so AGA/AGG/ATA/TGA are read correctly instead of with the standard code.
pub fn hgvsp_frameshift(
    protein_id: &str,
    ref_translateable: &[u8],
    alt_translateable: &[u8],
    affected_codon_start: usize, // 0-based codon index where the frameshift starts
    codon_table: &CodonTable,
) -> Option<String> {
    let prefix = format!("{}:p.", protein_id);

    // Translate both sequences from the affected codon onwards
    let ref_start = affected_codon_start * 3;
    if ref_start + 3 > ref_translateable.len() {
        return None;
    }
    if ref_start > alt_translateable.len() {
        return None;
    }

    let ref_peptide: Vec<u8> = ref_translateable[ref_start..]
        .chunks(3)
        .filter(|c| c.len() == 3)
        .map(|c| codon_table.translate(&[c[0], c[1], c[2]]))
        .collect();

    let alt_peptide: Vec<u8> = alt_translateable[ref_start..]
        .chunks(3)
        .filter(|c| c.len() == 3)
        .map(|c| codon_table.translate(&[c[0], c[1], c[2]]))
        .collect();

    // Find the first position where amino acids differ
    let mut first_changed_offset = 0;
    for i in 0..ref_peptide.len().min(alt_peptide.len()) {
        if ref_peptide[i] != alt_peptide[i] {
            first_changed_offset = i;
            break;
        }
        // If we reach a stop codon in ref before finding a change,
        // the change starts at this position
        if ref_peptide[i] == b'*' {
            first_changed_offset = i;
            break;
        }
        first_changed_offset = i + 1;
    }

    if first_changed_offset >= ref_peptide.len() && first_changed_offset >= alt_peptide.len() {
        return None;
    }

    let first_changed_pos = affected_codon_start + first_changed_offset + 1; // 1-based
    let ref_aa = if first_changed_offset < ref_peptide.len() {
        ref_peptide[first_changed_offset]
    } else {
        b'X'
    };
    let alt_aa = if first_changed_offset < alt_peptide.len() {
        alt_peptide[first_changed_offset]
    } else {
        b'X'
    };

    let ref_aa3 = aa_one_to_three(ref_aa);
    let alt_aa3 = aa_one_to_three(alt_aa);

    // Find the new stop codon position in the alt sequence.
    // If the sequence contains unresolved (X) amino acids, use Ter? to indicate uncertainty.
    let mut stop_dist = None;
    let mut hit_unresolved = false;
    let unresolved_count = alt_peptide[first_changed_offset..]
        .iter()
        .take(10)
        .filter(|&&aa| aa == b'X')
        .count();
    let mostly_unresolved = unresolved_count > 5;
    if !mostly_unresolved {
        for (i, &aa) in alt_peptide.iter().enumerate().skip(first_changed_offset) {
            if aa == b'*' {
                stop_dist = Some(i - first_changed_offset + 1);
                break;
            }
            if aa == b'X' {
                hit_unresolved = true;
            }
        }
    } else {
        hit_unresolved = true;
    }

    if let Some(d) = stop_dist {
        Some(format!(
            "{}{}{}{}fsTer{}",
            prefix, ref_aa3, first_changed_pos, alt_aa3, d
        ))
    } else if hit_unresolved || mostly_unresolved {
        // Sequence has unresolved regions - can't determine stop position
        Some(format!(
            "{}{}{}{}fsTer?",
            prefix, ref_aa3, first_changed_pos, alt_aa3
        ))
    } else {
        // No stop found and sequence is clean - true extension
        Some(format!(
            "{}{}{}{}fsTer?",
            prefix, ref_aa3, first_changed_pos, alt_aa3
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fastvep_genome::mitochondrial_codon_table;

    #[test]
    fn test_hgvsp_frameshift_mitochondrial_table_differs() {
        // Same ref/alt translateable sequences, only the codon table differs.
        // Codon 0 changes (Arg CGT -> Pro CCC, same under both tables), so
        // the frameshift starts there regardless of table. Codon 1 is TGA:
        // a stop under the standard table but Trp under the vertebrate
        // mitochondrial table (NCBI table 2), so the two tables must find
        // the new stop codon (Ter) at different downstream distances.
        let ref_translateable = b"CGTCGTCGTCGT"; // Arg Arg Arg Arg
        let alt_translateable = b"CCCTGAAAATAA"; // Pro TGA(*/W) Lys TAA(*)

        let standard = CodonTable::standard();
        let mitochondrial = mitochondrial_codon_table();

        let standard_result =
            hgvsp_frameshift("ENSP1", ref_translateable, alt_translateable, 0, &standard);
        let mito_result = hgvsp_frameshift(
            "ENSP1",
            ref_translateable,
            alt_translateable,
            0,
            &mitochondrial,
        );

        // Standard table: TGA is a stop, so the new terminator is 2 codons in.
        assert_eq!(standard_result, Some("ENSP1:p.Arg1ProfsTer2".to_string()));
        // Mitochondrial table: TGA reads as Trp, so translation continues
        // past it to the real stop (TAA) 4 codons in.
        assert_eq!(mito_result, Some("ENSP1:p.Arg1ProfsTer4".to_string()));
        assert_ne!(standard_result, mito_result);
    }

    #[test]
    fn test_hgvsp_frameshift_short_alt_translateable_returns_none() {
        // Regression: there's a bounds check guarding `ref_translateable`
        // (`ref_start + 3 > ref_translateable.len()`) but nothing equivalent
        // guarded `alt_translateable[ref_start..]` on the next line. If the
        // alt sequence is shorter than `ref_start`, that slice must not
        // panic ("start index out of range") -- it should return None, same
        // as the existing ref-side guard.
        let ref_translateable = b"CGTCGTCGTCGT"; // 12 bases, ref_start=3 is in-bounds
        let alt_translateable = b"CC"; // only 2 bases -- shorter than ref_start (3)

        let standard = CodonTable::standard();
        let result = hgvsp_frameshift("ENSP1", ref_translateable, alt_translateable, 1, &standard);
        assert_eq!(result, None);
    }

    #[test]
    fn test_hgvsp_missense() {
        let result = hgvsp("ENSP00000001", 41, b'R', b'K', false);
        assert_eq!(result, Some("ENSP00000001:p.Arg41Lys".to_string()));
    }

    #[test]
    fn test_hgvsp_start_lost_is_unknown() {
        let result = hgvsp("ENSP00000001", 1, b'M', b'V', false);
        assert_eq!(result, Some("ENSP00000001:p.Met1?".to_string()));
    }

    #[test]
    fn test_hgvsp_substitution_uses_first_changed_residue() {
        let result = hgvsp_substitution("ENSP00000001", 366, "ML", "IL");
        assert_eq!(result, Some("ENSP00000001:p.Met366Ile".to_string()));
    }

    #[test]
    fn test_hgvsp_synonymous() {
        let result = hgvsp("ENSP00000001", 41, b'R', b'R', false);
        assert_eq!(result, Some("ENSP00000001:p.Arg41=".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_deletion_single() {
        // single-residue in-frame deletion
        let r = hgvsp_inframe_indel("ENSP00000001", 157, "F", "-", None);
        assert_eq!(r, Some("ENSP00000001:p.Phe157del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_deletion_range() {
        // multi-residue in-frame deletion (regression for the p.Tyr43??? bug)
        let r = hgvsp_inframe_indel("ENSP00000001", 43, "YXQ", "-", None);
        assert_eq!(r, Some("ENSP00000001:p.Tyr43_Gln45del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_delins() {
        // in-frame deletion-insertion
        let r = hgvsp_inframe_indel("ENSP00000001", 2173, "NL", "K", None);
        assert_eq!(
            r,
            Some("ENSP00000001:p.Asn2173_Leu2174delinsLys".to_string())
        );
    }

    fn peptide_with(at: u64, residues: &str, length: usize) -> String {
        let mut peptide = vec![b'M'; length];
        for (offset, residue) in residues.bytes().enumerate() {
            peptide[at as usize - 1 + offset] = residue;
        }
        String::from_utf8(peptide).unwrap()
    }

    #[test]
    fn test_hgvsp_inframe_insertion_without_peptide_is_delins() {
        let r = hgvsp_inframe_indel("ENSP1", 3, "R", "RR", None);
        assert_eq!(r, Some("ENSP1:p.Arg3delinsArgArg".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_insertion_normalizes_to_duplication() {
        let r = hgvsp_inframe_indel("ENSP1", 3, "R", "RR", Some("MWR*"));
        assert_eq!(r, Some("ENSP1:p.Arg3dup".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_multi_residue_duplication() {
        let duplicated = "NEYFYVDFREYEYD";
        let peptide = peptide_with(587, &format!("{}K", duplicated), 700);
        let alternate = format!("D{}", duplicated);
        let r = hgvsp_inframe_indel("ENSP1", 600, "D", &alternate, Some(&peptide));
        assert_eq!(r, Some("ENSP1:p.Asn587_Asp600dup".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_change_uses_rightmost_repeat() {
        let r = hgvsp_inframe_indel("ENSP1", 2, "A", "-", Some("MAAA*"));
        assert_eq!(r, Some("ENSP1:p.Ala4del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_insertion_uses_flanking_residues() {
        let r = hgvsp_inframe_indel("ENSP1", 2, "W", "WQ", Some("MWR*"));
        assert_eq!(r, Some("ENSP1:p.Trp2_Arg3insGln".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_true_insertion_uses_ins_form() {
        let peptide = peptide_with(92, "SSK", 200);
        let r = hgvsp_inframe_indel("ENSP1", 92, "S", "SG", Some(&peptide));
        assert_eq!(r, Some("ENSP1:p.Ser92_Ser93insGly".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_normalizes_peptides_with_leading_x_and_stop() {
        let r = hgvsp_inframe_indel("ENSP1", 3, "R", "RR", Some("XWR*"));
        assert_eq!(r, Some("ENSP1:p.Arg3dup".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_delins_is_not_shifted() {
        let r = hgvsp_inframe_indel("ENSP1", 2, "AP", "PA", Some("MAPAPQ*"));
        assert_eq!(r, Some("ENSP1:p.Ala2_Pro3delinsProAla".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_unusable_peptides_fall_back() {
        let cases = [
            (6, "KX", "-", "MAAAGK", Some("ENSP1:p.Lys6_Xaa7del")),
            (100, "AK", "-", "MAAAGK", Some("ENSP1:p.Ala100_Lys101del")),
            (500, "-", "R", "MAAAGK", None),
            (1, "F", "-", "", Some("ENSP1:p.Phe1del")),
            (7, "K", "KG", "MAAAGK", Some("ENSP1:p.Lys7delinsLysGly")),
            (0, "F", "-", "MAAAGK", Some("ENSP1:p.Phe0del")),
        ];
        for (start, reference, alternate, peptide, expected) in cases {
            let got = hgvsp_inframe_indel("ENSP1", start, reference, alternate, Some(peptide));
            assert_eq!(got.as_deref(), expected);
        }
    }

    #[test]
    fn test_hgvsp_inframe_does_not_use_the_terminator() {
        let deletion = hgvsp_inframe_indel("ENSP1", 2, "K", "-", Some("MKKG*"));
        assert_eq!(deletion, Some("ENSP1:p.Lys3del".to_string()));

        let insertion = hgvsp_inframe_indel("ENSP1", 4, "G", "GS", Some("MKKG*"));
        assert_eq!(insertion, Some("ENSP1:p.Gly4delinsGlySer".to_string()));
        assert!(!insertion.unwrap().contains("Ter"));
    }

    #[test]
    fn test_hgvsp_inframe_ignores_inconsistent_peptide() {
        let r = hgvsp_inframe_indel("ENSP1", 2, "F", "-", Some("MQQQQQK"));
        assert_eq!(r, Some("ENSP1:p.Phe2del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_terminal_insertion_is_not_dropped() {
        for start in [1, 8] {
            assert!(hgvsp_inframe_indel("ENSP1", start, "G", "GG", Some("MKKRSTV")).is_some());
        }
    }

    #[test]
    fn test_hgvsp_inframe_observed_insertions_never_use_substitution_shape() {
        let insertions = [
            (41, "G", "AG", "p.Gly41Ala"),
            (185, "W", "WR", "p.Trp185="),
            (92, "S", "SG", "p.Ser92="),
            (2927, "E", "DE", "p.Glu2927Asp"),
            (375, "Q", "PLGPAKPPAQQ", "p.Gln375Pro"),
            (1829, "G", "GSSG", "p.Gly1829="),
            (510, "P", "QP", "p.Pro510Gln"),
            (510, "P", "QQP", "p.Pro510Gln"),
            (498, "Q", "QQ", "p.Gln498="),
            (600, "D", "DNEYFYVDFREYEYD", "p.Asp600="),
            (598, "E", "DVDFREYE", "p.Glu598Asp"),
            (596, "E", "VPSDNEYFYVDFRE", "p.Glu596Val"),
            (439, "I", "IKKK", "p.Ile439="),
            (
                1688,
                "S",
                "CSKDLEAFNPESKELLDLVEFTNEIQTLLGSS",
                "p.Ser1688Cys",
            ),
            (188, "S", "SD", "p.Ser188="),
        ];

        for (start, reference, alternate, previous) in insertions {
            let peptide = peptide_with(start, reference, start as usize + reference.len() + 64);
            for context in [None, Some(peptide.as_str())] {
                let output = hgvsp_inframe_indel("ENSP1", start, reference, alternate, context)
                    .expect("in-frame insertion must produce HGVSp");
                let change = output.split(":p.").nth(1).unwrap();
                assert!(
                    change.ends_with("dup") || change.contains("delins") || change.contains("ins"),
                    "{reference}/{alternate} at {start} rendered as {output}"
                );
                assert!(!output.ends_with('='));
                assert_ne!(output.split(':').nth(1).unwrap(), previous);
            }
        }
    }

    #[test]
    fn test_hgvsp_stop_gained() {
        let result = hgvsp("ENSP00000001", 41, b'R', b'*', false);
        assert_eq!(result, Some("ENSP00000001:p.Arg41Ter".to_string()));
    }

    #[test]
    fn test_hgvsp_frameshift() {
        let result = hgvsp("ENSP00000001", 41, b'R', b'X', true);
        assert_eq!(result, Some("ENSP00000001:p.Arg41fs".to_string()));
    }

    #[test]
    fn test_hgvsp_stop_lost() {
        let result = hgvsp("ENSP00000001", 100, b'*', b'R', false);
        assert_eq!(result, Some("ENSP00000001:p.Arg100ext*?".to_string()));
    }
}
