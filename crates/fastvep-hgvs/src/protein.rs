use fastvep_core::{Allele, Strand};
use fastvep_genome::codon::{aa_one_to_three, CodonTable};

// VEP 115.2 converts BioPerl's Xaa spelling to Ter immediately before it
// formats protein HGVS. Keep the conversion local to HGVSp: Amino_acids still
// uses X to represent the translated cache value.
fn hgvs_aa_one_to_three(aa: u8) -> &'static str {
    if aa == b'X' {
        "Ter"
    } else {
        aa_one_to_three(aa)
    }
}

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
    let ref_aa3 = hgvs_aa_one_to_three(ref_aa);

    if is_frameshift {
        return Some(format!("{}{}{}fs", prefix, ref_aa3, protein_pos));
    }

    if ref_aa3 == hgvs_aa_one_to_three(alt_aa) {
        // Synonymous
        return Some(format!("{}{}{}=", prefix, ref_aa3, protein_pos));
    }

    let alt_aa3 = hgvs_aa_one_to_three(alt_aa);

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

/// Format the protein uncertainty for a consequence already classified as
/// `start_lost`. VEP applies that predicate before its protein formatter; the
/// formatter clips common peptide ends and does not infer the consequence from
/// placeholder residues.
pub fn hgvsp_start_lost(
    protein_id: &str,
    protein_start: u64,
    ref_aas: &str,
    alt_aas: &str,
    ref_peptide: Option<&[u8]>,
) -> Option<String> {
    hgvsp_inframe_indel_with_context(
        protein_id,
        protein_start,
        protein_start + ref_aas.len().saturating_sub(1) as u64,
        ref_aas,
        alt_aas,
        ref_peptide,
        Strand::Forward,
        true, None, None,
)
}

/// Render a 1-based inclusive residue range, e.g. `Gly41` or `Asn587_Asp600`.
fn residue_span(first_pos: u64, residues: &[u8]) -> String {
    let first = hgvs_aa_one_to_three(residues[0]);
    if residues.len() == 1 {
        format!("{}{}", first, first_pos)
    } else {
        let last = hgvs_aa_one_to_three(residues[residues.len() - 1]);
        format!(
            "{}{}_{}{}",
            first,
            first_pos,
            last,
            first_pos + residues.len() as u64 - 1
        )
    }
}

/// Render a change that takes out the initiation codon, whose downstream
/// consequence is unknown: `Glu1?` for a single residue, `MetAsnIle1_?3` for a
/// run. The residues are still named - the annotation does say what was
/// affected - and the `?` marks the part that cannot be resolved from sequence.
///
/// `start` is the position of the first residue named. It is not always 1: a
/// change that begins at the initiator can share its leading residue with the
/// replacement, and Ensembl trims that off first. `XP/XAKSTVGA` on
/// ENSP00000500558 is `p.Pro2?` in VEP 115.1, not `p.XaaPro1_?2`.
///
/// The run's second coordinate is written as the last residue named, which is
/// what `p.Met1_?N` reads as. Not verified: every multi-residue case in the
/// ClinVar sample starts at residue 1, where the last position and the residue
/// count are the same number, so nothing here distinguishes them.
fn uncertain_from_initiator(start: u64, residues: &[u8]) -> String {
    if residues.len() == 1 {
        format!("{}{}?", hgvs_aa_one_to_three(residues[0]), start)
    } else {
        format!(
            "{}{}_?{}",
            three_letter(residues),
            start,
            start + residues.len() as u64 - 1
        )
    }
}

fn three_letter(residues: &[u8]) -> String {
    residues.iter().map(|&b| hgvs_aa_one_to_three(b)).collect()
}

// VEP clips the complete peptide pair first, then its final ins/delins
// formatter removes anything after the first translated terminator.
fn through_terminator(residues: &[u8]) -> &[u8] {
    match residues
        .iter()
        .position(|&residue| matches!(residue, b'*' | b'X'))
    {
        Some(stop) => &residues[..=stop],
        None => residues,
    }
}

fn clip_residues(
    mut reference: Vec<u8>,
    mut alternate: Vec<u8>,
    mut start: u64,
) -> (Vec<u8>, Vec<u8>, u64) {
    // VEP _clip_alleles returns the untouched notation if prefix scanning
    // reaches a recreated stop, before committing any preceding prefix trim.
    if reference.iter().zip(&alternate)
        .take_while(|(r, a)| r == a)
        .any(|(r, _)| *r == b'*')
    {
        return (reference, alternate, start);
    }
    while !reference.is_empty() && !alternate.is_empty() && reference[0] == alternate[0] {
        reference.remove(0);
        alternate.remove(0);
        start += 1;
    }
    while !reference.is_empty() && !alternate.is_empty() && reference.last() == alternate.last() {
        reference.pop();
        alternate.pop();
    }
    (reference, alternate, start)
}

/// Describe the change using only the residues the caller supplied, without
/// consulting the peptide. Always available when there is at least one
/// reference residue, and the fallback whenever the peptide is missing or
/// cannot be trusted — a valid, unshifted description beats emitting nothing.
fn unshifted_description(
    prefix: &str,
    start: u64,
    end: u64,
    reference: &[u8],
    alternate: &[u8],
) -> Option<String> {
    if reference.is_empty() {
        return None;
    }
    // VEP peptide() suppresses a partial X after a sole stop, but keeps the
    // translation endpoints. Its delins formatter still names both endpoints.
    let range = if reference == b"*" && end > start && !alternate.is_empty() {
        format!("Ter{}_Ter{}", start, end)
    } else {
        residue_span(start, reference)
    };
    if alternate.is_empty() {
        Some(format!("{}{}del", prefix, range))
    } else {
        Some(format!(
            "{}{}delins{}",
            prefix,
            range,
            three_letter(through_terminator(alternate))
        ))
    }
}

/// Whether `peptide` actually carries `reference` at `protein_start`.
///
/// A transcript whose sequence disagrees with its own coordinates would
/// otherwise produce a confident, well-formed, wrong description — worse than
/// no normalisation at all.
fn peptide_carries(peptide: &[u8], protein_start: u64, reference: &[u8]) -> bool {
    let Some(lo) = protein_start.checked_sub(1).map(|v| v as usize) else {
        return false;
    };
    if reference.is_empty() {
        // Pure insertion: only the flanking position is read.
        return lo <= peptide.len();
    }
    peptide.get(lo..lo + reference.len()) == Some(reference)
}

/// Whether `protein_start` names the *end* of the affected span rather than its
/// start.
///
/// `protein_start` is derived from the genomic left edge. On a transcript that
/// runs left to right that edge is the first affected residue; on one that runs
/// right to left it is the last. The residues themselves are built from the
/// lower of the two CDS coordinates for a shrinking change and from `cds_start`
/// otherwise (see `predict_coding_consequence`), so it is exactly the shrinking
/// change on the reverse strand whose span arrives anchored at its end. An
/// insertion is anchored at `protein_start` on either strand.
///
/// This is a property of the coordinate convention, not of the sequence, which
/// is why it can be decided here rather than guessed from the peptide.
fn anchored_at_span_end(strand: Strand, reference_len: usize, alternate_len: usize) -> bool {
    strand == Strand::Reverse && alternate_len < reference_len
}

/// The positions at which `reference` may sit, given the caller's anchor, most
/// likely first.
///
/// `protein_start` is the first affected residue on the forward strand, but for
/// a shrinking change it can be the *last*: the field is derived from the
/// genomic left edge, which is the end of the affected range when the transcript
/// runs right to left (see #89). Both readings describe the same span, so a
/// reference of length n is anchored either at `protein_start` or at
/// `protein_start - (n - 1)`.
///
/// Those two and no others. The pair is one span read from either end, so any
/// further candidate would be a scan for a coincidental repeat rather than a
/// consequence of the coordinate convention. An empty or single-residue
/// reference has one candidate, which leaves insertions - where `protein_start`
/// is already the far end of the pair - on exactly the path they had before.
///
/// `end_first` puts them in the order [`anchored_at_span_end`] determines, which
/// is what makes the choice between them a reading of the convention rather than
/// a coincidence. Taking the first *corroborated* candidate is not enough on its
/// own: where the reference is periodic with period n-1 both ends corroborate,
/// and before #96 the wrong one won whenever it was tried first. `MEGEGEA` with
/// `EGE` arriving at `protein_start` 4 is the smallest case - residues 2-4 and
/// 4-6 both read `EGE`, and deleting them leaves different proteins (`MGEA` and
/// `MEGA`), so the 3'-rule cannot reconcile the two descriptions afterwards.
///
/// The other end stays as a fallback rather than being dropped. The residues are
/// not guaranteed to sit at either scalar - over 37,122 in-frame ClinVar rows
/// they sat at `protein_start` in 71.1% of cases, one residue earlier in 11.9%,
/// and at neither in 17.0% - so ordering is all the strand can honestly buy. The
/// peptide still decides.
fn anchor_candidates(
    protein_start: u64,
    reference_len: usize,
    end_first: bool,
) -> [Option<u64>; 2] {
    let from_end = reference_len
        .checked_sub(1)
        .filter(|&back| back > 0)
        .and_then(|back| protein_start.checked_sub(back as u64))
        // Residues are numbered from 1, so 0 is not a position to try.
        .filter(|&anchor| anchor > 0);
    match from_end {
        Some(other) if end_first => [Some(other), Some(protein_start)],
        _ => [Some(protein_start), from_end],
    }
}

/// Generate HGVSp notation for an in-frame indel — deletion, delins, insertion
/// or duplication.
///
/// `ref_aas` are the affected reference residues (one-letter, in order);
/// `alt_aas` is the replacement ("-" or empty for a pure deletion, longer than
/// `ref_aas` for an insertion).
///
///   ENSP0:p.Phe157del            single-residue deletion
///   ENSP0:p.Tyr43_Gln45del       multi-residue deletion
///   ENSP0:p.Asn2173_Leu2174delinsLys   delins
///   ENSP0:p.Ser92_Ser93insGly    insertion
///   ENSP0:p.Asn587_Asp600dup     insertion that repeats what it follows
///
/// This is the only correct rendering for an in-frame indel: `hgvsp` compares
/// just the first residue of each side, so it would describe these as a
/// substitution — synonymous when that residue is unchanged, missense when it
/// differs, and neither is true of a variant that changes the protein's length.
///
/// `ref_peptide` is the reference protein in one-letter codes, residue 1 at
/// index 0 — pass `Transcript::peptide`, which is translated in frame from
/// `codon_table_start_phase` and with the right codon table for the contig. Any
/// terminator it carries is trimmed here. Given one, insertions and deletions
/// are normalised per the HGVS 3'-rule - shifted as far C-terminal as they can
/// go, with duplications collapsed to `dup`. `delins` is not shifted.
///
/// The shift bound intentionally follows Ensembl VEP 115.2's `_shift_3prime`
/// (`Bio::EnsEMBL::Variation::TranscriptVariationAllele`): for a change of `n`
/// residues, VEP only scans through `length(post_seq) - n`. This can halt
/// `n - 1` residues before the protein terminus, but reproducing that bound is
/// required for VEP-compatible HGVSp output.
///
/// `strand` says which end of `ref_aas` the caller's `protein_start` names,
/// which the coordinate convention fixes rather than leaves open: see
/// [`anchored_at_span_end`]. Pass the strand of the transcript the residues were
/// derived from - not the strand the caller finds convenient - because a wrong
/// one here reorders the anchor candidates and, where the reference is periodic,
/// selects a span the variant does not touch.
///
/// Every peptide-dependent step degrades to the unshifted description rather
/// than failing: a peptide that is absent, too short, or inconsistent with
/// `protein_start` still yields valid HGVS, just unnormalised.
pub fn hgvsp_inframe_indel(
    protein_id: &str,
    protein_start: u64,
    protein_end: u64,
    ref_aas: &str,
    alt_aas: &str,
    ref_peptide: Option<&[u8]>,
    strand: Strand,
) -> Option<String> {
    hgvsp_inframe_indel_with_context(
        protein_id, protein_start, protein_end, ref_aas, alt_aas, ref_peptide, strand, false, None, None,
)
}

#[allow(clippy::too_many_arguments)]
pub fn hgvsp_inframe_indel_with_context(
    protein_id: &str,
    protein_start: u64,
    protein_end: u64,
    ref_aas: &str,
    alt_aas: &str,
    ref_peptide: Option<&[u8]>,
    strand: Strand,
    start_lost: bool,
    duplication_peptide: Option<&[u8]>,
    full_reference_peptide: Option<&[u8]>,
) -> Option<String> {
    let strip = |s: &str| -> Vec<u8> {
        if s == "-" {
            Vec::new()
        } else {
            s.bytes().collect()
        }
    };
    let original_ref = strip(ref_aas);
    // VEP _get_surrounding_peptides appends original_ref when it begins
    // with a stop, restoring the terminal residue omitted by _peptide().
    let extended_peptide = ref_peptide
        .filter(|_| original_ref.starts_with(b"*"))
        .map(|peptide| [peptide.strip_suffix(b"*").unwrap_or(peptide), &original_ref].concat());
    let ref_peptide = extended_peptide.as_deref().or(ref_peptide);
    // Nothing past a terminator the change introduces is translated, so those
    // residues are not part of the protein and must not be named. Ensembl's
    // `Amino_acids` column keeps the whole translated window - `SL/MEP*S` - but
    // its HGVSp for the same row is `p.Ser269_Leu270delinsMetGluProTer`, and
    // naming the `S` after the `Ter` would describe a residue that does not
    // exist. The terminator itself is kept: it is the last thing the protein has.
    let original_alt = strip(alt_aas);
    // A pure insertion replaces no residue, so it is written between two of
    // them - and the caller's pair names both. `protein_start` comes from the
    // genomic left edge, which is the *upper* residue on the forward strand and
    // the lower one on the reverse; the insertion sits in front of the higher of
    // the two whichever way the transcript runs. Taking `protein_start` alone
    // put every reverse-strand duplication one residue early, which reads as an
    // ordinary insertion rather than a `dup`: `p.Gly559_Asp560insHisGluAsnLys...`
    // where VEP writes `p.His553_Asp560dup`.
    let protein_start = if original_ref.is_empty() {
        protein_start.max(protein_end)
    } else {
        protein_start
    };
    let prefix = format!("{}:p.", protein_id);

    // VEP writes a synonymous multi-residue window as all three-letter residues
    // followed by the first position, for example `p.SerTer22=`.
    if !original_ref.is_empty() && original_ref == original_alt {
        let lo = protein_start.min(protein_end);
        // VEP skips clipping equal peptides, but still applies start_lost
        // before the final synonymous-format check.
        if start_lost {
            return Some(format!("{}{}", prefix, uncertain_from_initiator(lo, &original_ref)));
        }
        return Some(format!("{}{}{}=", prefix, three_letter(&original_ref), lo));
    }

    // VEP sorts the translation endpoints for a delins range. Even when a
    // terminal X is absent from the full peptide, the supplied endpoints still
    // locate the window; the genomic-left endpoint alone is wrong in reverse.
    let fallback = || {
        if start_lost {
            let (reference, _, start) = clip_residues(
                original_ref.clone(), original_alt.clone(), protein_start.min(protein_end),
            );
            return (!reference.is_empty())
                .then(|| format!("{}{}", prefix, uncertain_from_initiator(start, &reference)));
        }
        // A terminal stop is intentionally absent from VEP's full peptide.
        // Clip its shrinking local window even without full-peptide support;
        // retain the existing fallback for unrelated unavailable sequences.
        let (reference, alternate, start) = if matches!(original_ref.last(), Some(b'*' | b'X'))
            && original_alt.len() < original_ref.len()
            && !(original_ref.first() == Some(&b'*') && original_alt.first() == Some(&b'*')) {
            clip_residues(original_ref.clone(), original_alt.clone(), protein_start.min(protein_end))
        } else {
            (original_ref.clone(), original_alt.clone(), protein_start.min(protein_end))
        };
        let prefix_trim = start - protein_start.min(protein_end);
        let suffix_trim = original_ref.len().saturating_sub(prefix_trim as usize + reference.len());
        unshifted_description(
            &prefix, start, protein_start.max(protein_end).saturating_sub(suffix_trim as u64),
            &reference, &alternate,
        )
    };
    if original_ref == b"*" && protein_start != protein_end {
        return fallback();
    }
    // VEP `_clip_alleles` does not clip a recreated leading stop. The final
    // delins formatter still truncates the alternate at that terminator.
    if original_ref.first() == Some(&b'*') && original_alt.first() == Some(&b'*') {
        return fallback();
    }
    // Everything below is positioned relative to an anchor, so none of it is
    // safe unless the peptide corroborates that the caller's residues really sit
    // there. They do not always: for a shrinking change like `FF/F` the call
    // sites pass the end of the affected range rather than its start, so try
    // reading the span from its other end before giving up.
    //
    // Trusting only `protein_start` was not merely a missed normalisation. The
    // un-normalised description takes its residue letters from `ref_aas` and its
    // numbers from the anchor, so a wrong anchor emits a position the protein
    // does not have - `p.Ser1092_Phe1094delinsPhe` on ENSP00000247087 where
    // residues 1092-1094 are FQP and the SSF is at 1090-1092. The letters and
    // the numbers contradict each other. Over 47,013 ClinVar indels that was
    // 17,654 descriptions.
    let Some((peptide, protein_start)) = ref_peptide.and_then(|p| {
        let end_first = protein_start >= protein_end
            && anchored_at_span_end(strand, original_ref.len(), original_alt.len());
        anchor_candidates(protein_start, original_ref.len(), end_first)
            .into_iter()
            .flatten()
            .find(|&anchor| peptide_carries(p, anchor, &original_ref))
            .map(|anchor| (p, anchor))
            .or_else(|| {
                let anchor = protein_start.min(protein_end);
                let (reference, alternate, clipped_start) =
                    clip_residues(original_ref.clone(), original_alt.clone(), anchor);
                // A partial X is absent from the full peptide, but clipping
                // it can reveal an insertion. Still check duplication before
                // requiring the two flanking residues, as VEP does.
                (original_ref.last() == Some(&b'X')
                    && ((!reference.is_empty() && peptide_carries(p, clipped_start, &reference))
                        || (reference.is_empty() && !alternate.is_empty())))
                    .then_some((p, anchor))
            })
    }) else {
        return fallback();
    };

    // Corroborate the allele window before dropping the final CDS terminator:
    // a terminal S*/* deletion must still clip to S/-. VEP's surrounding
    // reference peptide retains internal stops but omits the final one.
    // Corroboration and clipping use the allele/consequence peptide. VEP's
    // post-sequence lookup and displayed flanks use Transcript::_peptide,
    // whose initiator normalization precedes source sequence edits.
    let extended_reference = full_reference_peptide
        .filter(|_| original_ref.starts_with(b"*"))
        .map(|p| [p.strip_suffix(b"*").unwrap_or(p), &original_ref].concat());
    let peptide = extended_reference.as_deref().or(full_reference_peptide).unwrap_or(peptide);
    let peptide = if original_ref.starts_with(b"*") {
        peptide
    } else {
        peptide.strip_suffix(b"*").unwrap_or(peptide)
    };

    // Reduce to the minimal changed region: residues shared at either end are
    // not part of the description. Trimming the prefix moves the start right.
    let (reference, alternate, start) =
        clip_residues(original_ref.clone(), original_alt.clone(), protein_start);

    if reference.is_empty() && alternate.is_empty() {
        return None;
    }

    if reference.is_empty() {
        // Pure insertion, sitting between residues (start - 1) and start.
        let Some(mut at) = start.checked_sub(1).map(|v| v as usize) else {
            return fallback();
        };
        let mut inserted = alternate;
        // _get_hgvs_protein_type replaces the first stop with X before
        // post-sequence shifting and duplication checks.
        if let Some(stop) = inserted.iter_mut().find(|aa| **aa == b'*') {
            *stop = b'X';
        }
        // 3'-rule: slide right while the residue the insertion sits in front of
        // is the one it would place there.
        // VEP requires a residue beyond the requested post-sequence start.
        let can_shift = at + 1 < peptide.len();
        while can_shift && at
            .checked_add(inserted.len())
            .is_some_and(|end| end <= peptide.len())
            && peptide[at] == inserted[0]
        {
            inserted.rotate_left(1);
            at += 1;
        }
        // A duplication is an insertion whose residues repeat those immediately
        // before it.
        // VEP _check_for_peptide_duplication translates CDS afresh, whereas
        // surrounding residues use the cached, sequence-edited peptide.
        let preceding = at
            .checked_sub(inserted.len())
            .and_then(|lo| duplication_peptide.unwrap_or(peptide).get(lo..at));
        if !inserted
            .iter()
            .any(|&residue| residue == b'*')
            && preceding == Some(inserted.as_slice())
        {
            let dup_start = (at - inserted.len() + 1) as u64;
            return Some(format!(
                "{}{}dup",
                prefix,
                if inserted.len() == 1 {
                    format!("{}{}", aa_one_to_three(inserted[0]), dup_start)
                } else {
                    format!("{}{}_{}{}", aa_one_to_three(inserted[0]), dup_start,
                        aa_one_to_three(*inserted.last()?), dup_start + inserted.len() as u64 - 1)
                }
            ));
        }
        if start_lost {
            // VEP checks duplication before start_lost, then obtains the
            // flanking peptide with substr(ref, min(start,end)-1, 2).
            // At residue zero Perl's negative offset selects the final residue.
            if at >= peptide.len() {
                return None;
            }
            let flank_start = at.checked_sub(1).unwrap_or(peptide.len() - 1);
            let flanks = &peptide[flank_start..(flank_start + 2).min(peptide.len())];
            return Some(format!("{}{}{}_?{}", prefix, three_letter(flanks), at + 1, at));
        }
        // Clipping a reference window down to an insertion before residue one
        // gives VEP substr(peptide, -1, 2): one final residue, used for both
        // names in the 0_1 insertion form (F1 RER1 Phe0_Phe1insLys).
        if at == 0 && !original_ref.is_empty() {
            let flank = hgvs_aa_one_to_three(*peptide.last()?);
            return Some(format!("{}{flank}0_{flank}1ins{}", prefix, three_letter(through_terminator(&inserted))));
        }
        // Otherwise name the flanking pair. VEP emits no HGVSp at a terminus:
        // `_get_surrounding_peptides` cannot return the required two residues.
        match (
            at.checked_sub(1).and_then(|i| peptide.get(i)),
            peptide.get(at),
        ) {
            (Some(&before), Some(&after)) => Some(format!(
                "{}{}{}_{}{}ins{}",
                prefix,
                hgvs_aa_one_to_three(before),
                at,
                hgvs_aa_one_to_three(after),
                at + 1,
                three_letter(through_terminator(&inserted))
            )),
            _ => None,
        }
    } else if alternate.is_empty() {
        // Pure deletion of `reference` starting at `start`.
        let mut at = start;
        let mut residues = reference;
        {
            let len = residues.len();
            // 3'-rule: slide right while the residue following the deleted block
            // repeats the first deleted residue.
            // `_get_surrounding_peptides(end + 1)` refuses a post-sequence
            // starting at the final residue. This is an initial bound only;
            // once accepted, `_shift_3prime` can consume that final residue.
            let can_read_post_sequence = at as usize + len < peptide.len();
            let mut rotation = 0;
            while can_read_post_sequence {
                // VEP _shift_3prime rotates the local allele, which can differ
                // from the full protein's normalized initiator or sequence edits.
                let first = residues.get(rotation);
                let following_index = at as usize + len - 1;
                // VEP 115.2 passes the sequence after the deleted block to
                // `_shift_3prime`, then stops when fewer than `len` residues
                // remain in that post-sequence. Preserve that exact bound.
                let following = following_index
                    .checked_add(len)
                    .filter(|&end| end <= peptide.len())
                    .and_then(|_| peptide.get(following_index));
                match (first, following) {
                    (Some(a), Some(b)) if a == b => {
                        at += 1;
                        rotation = (rotation + 1) % len;
                    }
                    _ => break,
                }
            }
            residues.rotate_left(rotation);
        }
        if start_lost {
            return Some(format!(
                "{}{}",
                prefix,
                uncertain_from_initiator(at, &residues)
            ));
        }
        Some(format!("{}{}del", prefix, residue_span(at, &residues)))
    } else if start_lost {
        Some(format!("{}{}", prefix, uncertain_from_initiator(start, &reference)))
    } else if reference.len() == 1 && alternate.len() == 1 {
        // One residue for one residue is a substitution, whatever the window it
        // came from. A change spanning two codons that alters only the second
        // one - `EP/ET` - is `p.Pro154Thr` to Ensembl, not a two-residue delins
        // and not `p.Glu153=`, which is what reading the first residue of each
        // side gave. About 3,000 HGVSp rows per 6,600 ClinVar variants.
        let (r, a) = (reference[0], alternate[0]);
        Some(match (r, a) {
            // A terminator the change removes extends the protein by an unknown
            // amount; one it introduces ends it here.
            (b'*', _) => format!("{}{}{}ext*?", prefix, hgvs_aa_one_to_three(a), start),
            _ => format!(
                "{}{}{}{}",
                prefix,
                hgvs_aa_one_to_three(r),
                start,
                hgvs_aa_one_to_three(a)
            ),
        })
    } else {
        // Replacement of one residue run by another.
        Some(format!(
            "{}{}delins{}",
            prefix,
            residue_span(start, &reference),
            three_letter(through_terminator(&alternate))
        ))
    }
}

/// HGVSp for a frameshift, from the transcript's own sequence and the variant's
/// CDS coordinates.
///
/// The edit is the same one the codon window makes: replace the CDS bases the
/// reference allele covers with the alternate allele's, in transcript
/// orientation. Both per-variant loops used to open-code it, and both got it
/// wrong in the same three ways - they read `cds_start` as the low end of the
/// span (it is the *high* end on the reverse strand), they complemented the
/// inserted bases in place instead of reverse-complementing them, and they had
/// no case at all for a delins, so a replacement was inserted without removing
/// what it replaced. Over a 6,600-variant ClinVar sample that was 3,200 of
/// 3,794 frameshift-delins rows disagreeing with real VEP 115.1, plus 2,900 of
/// 11,596 frameshift deletions and 2,400 of 11,212 frameshift insertions.
///
/// `cds_and_downstream` must be CDS-indexed - its byte `n - 1` is CDS position
/// `n` - and run past the annotated terminator, because a frameshift's new stop
/// is often in what was the 3' UTR.
#[allow(clippy::too_many_arguments)] // each argument is an independent coordinate or allele
pub fn hgvsp_frameshift_from_cds(
    protein_id: &str,
    cds_and_downstream: &[u8],
    cds_start: Option<u64>,
    cds_end: Option<u64>,
    ref_allele: &Allele,
    alt_allele: &Allele,
    strand: Strand,
    codon_table: &CodonTable,
) -> Option<String> {
    hgvsp_frameshift_from_cds_with_tables(
        protein_id,
        cds_and_downstream,
        cds_start,
        cds_end,
        ref_allele,
        alt_allele,
        strand,
        codon_table,
        codon_table,
    )
}

/// Generate frameshift HGVSp while allowing VEP-compatible reference and
/// alternate translation tables to differ.
#[allow(clippy::too_many_arguments)]
pub fn hgvsp_frameshift_from_cds_with_tables(
    protein_id: &str,
    cds_and_downstream: &[u8],
    cds_start: Option<u64>,
    cds_end: Option<u64>,
    ref_allele: &Allele,
    alt_allele: &Allele,
    strand: Strand,
    reference_codon_table: &CodonTable,
    alternate_codon_table: &CodonTable,
) -> Option<String> {
    hgvsp_frameshift_from_cds_with_tables_and_ref_peptide(
        protein_id,
        cds_and_downstream,
        cds_start,
        cds_end,
        ref_allele,
        alt_allele,
        strand,
        reference_codon_table,
        alternate_codon_table,
        None,
        false,
        false,
    )
}

/// Generate VEP-compatible frameshift HGVSp using the transcript's annotated
/// reference peptide when it is available. The existing entry points retain
/// their signatures and fall back to translating the CDS.
#[allow(clippy::too_many_arguments)]
pub fn hgvsp_frameshift_from_cds_with_tables_and_ref_peptide(
    protein_id: &str,
    cds_and_downstream: &[u8],
    cds_start: Option<u64>,
    cds_end: Option<u64>,
    ref_allele: &Allele,
    alt_allele: &Allele,
    strand: Strand,
    reference_codon_table: &CodonTable,
    alternate_codon_table: &CodonTable,
    reference_peptide: Option<&[u8]>,
    stop_lost: bool,
    start_lost: bool,
) -> Option<String> {
    hgvsp_frameshift_from_cds_with_context(
        protein_id, cds_and_downstream, cds_start, cds_end, ref_allele, alt_allele,
        strand, reference_codon_table, alternate_codon_table, reference_peptide,
        stop_lost, start_lost, None,
    )
}

/// Rebuild the alternate CDS before appending its 3' UTR, using the annotated
/// CDS length rather than inferring it from a potentially incomplete peptide.
#[allow(clippy::too_many_arguments)]
pub fn hgvsp_frameshift_from_cds_with_context(
    protein_id: &str,
    cds_and_downstream: &[u8],
    cds_start: Option<u64>,
    cds_end: Option<u64>,
    ref_allele: &Allele,
    alt_allele: &Allele,
    strand: Strand,
    reference_codon_table: &CodonTable,
    alternate_codon_table: &CodonTable,
    reference_peptide: Option<&[u8]>,
    stop_lost: bool,
    start_lost: bool,
    reference_cds_length: Option<usize>,
) -> Option<String> {
    let (mut edited, first) = edited_cds(
        cds_and_downstream,
        cds_start,
        cds_end,
        ref_allele,
        alt_allele,
        strand,
    )?;

    let removed = cds_and_downstream.len() + alt_allele.len() - edited.len();
    // _get_alternate_cds trims a sub-codon CDS BEFORE appending the UTR.
    // Keep the original replacement span for the local peptide window below.
    if let Some(length) = reference_cds_length {
        let remaining = length.checked_sub(removed)?.checked_add(alt_allele.len())?;
        if remaining < 3 { edited.drain(..remaining); }
    }
    if edited.is_empty() || (reference_cds_length.is_none() && edited.len() < 3) {
        return None;
    }

    // VEP _get_fs_peptides returns the initial clipped peptide window when
    // the alternate translation ends before translation_start. Keep its full
    // deletion span instead of reducing it to the first reference residue.
    let codon_start = first / 3 * 3;
    let translate_window = |bases: &[u8], table: &CodonTable| {
        let mut peptide = table.translate_seq(bases);
        if bases.len() % 3 != 0 && peptide != b"*" { peptide.push(b'X'); }
        peptide
    };
    if codon_start >= edited.len() / 3 * 3 {
        let reference_end = (first + removed).div_ceil(3) * 3;
        let alternate_end = reference_end.checked_sub(removed)?.checked_add(alt_allele.len())?;
        let bounded_end = reference_end.min(reference_cds_length.unwrap_or(cds_and_downstream.len()));
        let reference = translate_window(cds_and_downstream.get(codon_start..bounded_end)?, reference_codon_table);
        let alternate = translate_window(edited.get(codon_start..alternate_end.min(edited.len()))?, alternate_codon_table);
        let start = first as u64 / 3 + 1;
        // hgvs_protein skips clipping identical windows, including X/X from
        // two incomplete terminal codons. The retained X becomes Ter ... del.
        let (reference, _, clipped_start) = if reference == alternate {
            (reference, alternate, start)
        } else {
            clip_residues(reference, alternate, start)
        };
        // _get_fs_peptides switches to deletion before start_lost formatting,
        // preserving the complete local reference and its clipped end.
        if start_lost {
            let end = start + reference.len() as u64 - 1 + (clipped_start - start);
            return Some(if start == end {
                format!("{}:p.{}{}?", protein_id, three_letter(&reference), start)
            } else {
                format!("{}:p.{}{}_?{}", protein_id, three_letter(&reference), start, end)
            });
        }
        // A codon-boundary insertion has no reference peptide. VEP's
        // _get_del_peptides clips the two empty terminal tails to synonymy.
        if reference.is_empty() {
            return Some(format!("{}:p.{}=", protein_id, start));
        }
        let end = start + reference.len() as u64 - 1;
        // _get_hgvs_protein_format retains both deletion endpoints when
        // stop_lost applies, even if the first deleted residue is not a stop.
        let extension = if stop_lost { "extTer?" } else { "" };
        return Some(if reference.len() == 1 {
            format!("{}:p.{}{}del{}", protein_id, hgvs_aa_one_to_three(reference[0]), start, extension)
        } else {
            format!("{}:p.{}{}_{}{}del{}", protein_id, hgvs_aa_one_to_three(reference[0]), start, hgvs_aa_one_to_three(*reference.last()?), end, extension)
        });
    }

    hgvsp_frameshift_with_tables(
        protein_id,
        cds_and_downstream,
        &edited,
        first / 3,
        reference_codon_table,
        alternate_codon_table,
        reference_peptide,
        stop_lost,
        start_lost.then(|| {
            let end = if *ref_allele == Allele::Deletion {
                cds_start.into_iter().chain(cds_end).min()
            } else {
                cds_start.into_iter().chain(cds_end).max()
            }.unwrap_or(1);
            let end = end.div_ceil(3);
            // hgvs_protein clips the local codon window before _get_fs_peptides
            // replaces its residues. That suffix trim still determines `end`.
            let reference_end = end as usize * 3;
            let alternate_end = reference_end.checked_sub(removed)?.checked_add(alt_allele.len())?;
            let bounded_end = reference_end.min(reference_cds_length.unwrap_or(cds_and_downstream.len()));
            let reference = translate_window(cds_and_downstream.get(codon_start..bounded_end)?, reference_codon_table);
            let alternate = translate_window(edited.get(codon_start..alternate_end.min(edited.len()))?, alternate_codon_table);
            if reference == alternate { return Some(end); }
            let length = reference.len();
            let start = first as u64 / 3 + 1;
            let (clipped, _, clipped_start) = clip_residues(reference, alternate, start);
            Some(end.saturating_sub(length.saturating_sub(clipped.len() + (clipped_start - start) as usize) as u64))
        }).flatten(),
    )
}

/// Rebuild an in-frame deletion whose shifted coding interval did not produce
/// a predictor peptide window. This follows VEP 115.2's `_get_del_peptides`:
/// translate the alternate CDS, compare both peptide tails from the shifted
/// translation start, and let the common-prefix/suffix clipping identify the
/// deleted residues.
#[allow(clippy::too_many_arguments)]
pub fn hgvsp_inframe_deletion_from_cds(
    protein_id: &str,
    cds_and_downstream: &[u8],
    cds_start: u64,
    cds_end: u64,
    ref_allele: &Allele,
    strand: Strand,
    codon_table: &CodonTable,
    reference_peptide: Option<&[u8]>,
    start_lost: bool,
    reference_cds_length: Option<usize>,
    full_reference_peptide: Option<&[u8]>,
) -> Option<String> {
    let (mut edited, _) = edited_cds(
        cds_and_downstream,
        Some(cds_start),
        Some(cds_end),
        ref_allele,
        &Allele::Deletion,
        strand,
    )?;
    let first = usize::try_from(cds_start.min(cds_end).checked_sub(1)?).ok()?;
    let codon_start = first / 3 * 3;
    let reference_end = usize::try_from(cds_start.max(cds_end).checked_add(2)? / 3 * 3).ok()?;
    let removed = usize::try_from(cds_start.abs_diff(cds_end).checked_add(1)?).ok()?;
    // VEP _trim_incomplete_codon empties a CDS shorter than one codon before
    // appending the UTR; longer partial codons remain in the alternate window.
    let remaining_cds = reference_cds_length.unwrap_or(cds_and_downstream.len()).checked_sub(removed)?;
    if remaining_cds < 3 { edited.drain(..remaining_cds); }
    let alternate_end = reference_end.checked_sub(removed)?;
    let protein_start = codon_start as u64 / 3 + 1;
    let bounded_end = reference_end.min(reference_cds_length.unwrap_or(cds_and_downstream.len()));
    let window = cds_and_downstream.get(codon_start..bounded_end)?;
    let mut reference = codon_table.translate_seq(window);
    // VEP peptide() retains the incomplete terminal CDS codon as X, while
    // the cached full peptide omits it. Preserve cached edits in complete codons.
    if let Some(peptide) = reference_peptide {
        for (offset, residue) in reference.iter_mut().enumerate() {
            if let Some(&cached) = peptide.get(codon_start / 3 + offset) {
                *residue = cached;
            }
        }
    }
    if window.len() % 3 != 0 && reference != b"*" {
        reference.push(b'X');
    }
    let alternate_window = edited.get(codon_start..alternate_end.min(edited.len()))?;
    let mut translated_alternate = codon_table.translate_seq(alternate_window);
    if alternate_window.len() % 3 != 0 && translated_alternate != b"*" {
        translated_alternate.push(b'X');
    }
    hgvsp_inframe_indel_with_context(
        protein_id,
        protein_start,
        protein_start + reference.len() as u64 - 1,
        std::str::from_utf8(&reference).ok()?,
        if translated_alternate.is_empty() {
            "-"
        } else {
            std::str::from_utf8(&translated_alternate).ok()?
        },
        reference_peptide,
        strand,
        start_lost, None, full_reference_peptide,
)
}

/// Rebuild a 3'-shifted in-frame insertion from the shifted CDS position and
/// rotated inserted sequence, matching VEP 115.2's shifted peptide window.
#[allow(clippy::too_many_arguments)]
pub fn hgvsp_inframe_insertion_from_cds(
    protein_id: &str,
    cds_and_downstream: &[u8],
    cds_start: u64,
    cds_end: u64,
    alt_allele: &Allele,
    strand: Strand,
    transcript_shift: u64,
    codon_table: &CodonTable,
    reference_peptide: Option<&[u8]>,
) -> Option<String> {
    hgvsp_inframe_insertion_from_cds_with_start_lost(
        protein_id, cds_and_downstream, cds_start, cds_end, alt_allele, strand,
        transcript_shift, codon_table, reference_peptide, false, None, None,
)
}

#[allow(clippy::too_many_arguments)]
pub fn hgvsp_inframe_insertion_from_cds_with_start_lost(
    protein_id: &str,
    cds_and_downstream: &[u8],
    cds_start: u64,
    cds_end: u64,
    alt_allele: &Allele,
    strand: Strand,
    transcript_shift: u64,
    codon_table: &CodonTable,
    reference_peptide: Option<&[u8]>,
    start_lost: bool,
    reference_cds_length: Option<usize>,
    full_reference_peptide: Option<&[u8]>,
) -> Option<String> {
    let duplication_peptide = CodonTable::standard().translate_seq(cds_and_downstream);
    let point = cds_start.min(cds_end);
    if cds_start.max(cds_end) != point.checked_add(1)? {
        return None;
    }
    let mut inserted = match alt_allele {
        Allele::Sequence(bases) if !bases.is_empty() => match strand {
            Strand::Forward => bases.clone(),
            Strand::Reverse => bases.iter().rev().map(|&base| complement(base)).collect(),
        },
        _ => return None,
    };
    let rotation = usize::try_from(transcript_shift % inserted.len() as u64).ok()?;
    // VEP's reverse shift_feature_seqs loop has a negative bound when the
    // shift exceeds seq_length; it then leaves the inserted sequence intact.
    if strand == Strand::Forward || transcript_shift <= inserted.len() as u64 {
        inserted.rotate_left(rotation);
    }
    let point = usize::try_from(point).ok()?;
    if point > cds_and_downstream.len() {
        return None;
    }
    // VEP's codon-boundary insertion has an empty reference peptide window.
    // Do not borrow a following codon from the UTR when the cached peptide
    // ends here (F1 ENST00000416415: His96_Ser98dup).
    if point % 3 == 0 && inserted.len() % 3 == 0 {
        let alternate = codon_table.translate_seq(&inserted);
        return hgvsp_inframe_indel_with_context(
            protein_id, point as u64 / 3 + 1, point as u64 / 3,
            "-", std::str::from_utf8(&alternate).ok()?, reference_peptide, strand, start_lost, Some(&duplication_peptide), full_reference_peptide,
);
    }
    let mut edited = Vec::with_capacity(cds_and_downstream.len() + inserted.len());
    edited.extend_from_slice(&cds_and_downstream[..point]);
    edited.extend_from_slice(&inserted);
    edited.extend_from_slice(&cds_and_downstream[point..]);
    let codon_start = point / 3 * 3;
    // VEP's reference codon is bounded by CDS, although its alternate
    // window may use downstream sequence. Do not complete a partial reference
    // codon with UTR bases (F1 KCTD21/FAM246C).
    let reference_end = reference_cds_length.unwrap_or(cds_and_downstream.len()).min(cds_and_downstream.len());
    let reference_window = cds_and_downstream.get(codon_start..(codon_start + 3).min(reference_end))?;
    let reference = reference_peptide
        .and_then(|peptide| peptide.get(codon_start / 3))
        .copied()
        .or_else(|| {
            codon_table
                .translate_seq(reference_window)
                .first()
                .copied()
        })
        .or_else(|| (!reference_window.is_empty() && reference_window.len() < 3).then_some(b'X'))?;
    let alternate_window = edited.get(codon_start..codon_start.checked_add(3 + inserted.len())?.min(edited.len()))?;
    let mut alternate = codon_table.translate_seq(alternate_window);
    // VEP peptide() retains a terminal partial codon as X, except after a
    // sole stop. A shifted insertion can complete that codon (X -> Gly).
    if alternate_window.len() % 3 != 0 && alternate != b"*" {
        alternate.push(b'X');
    }
    if reference == b'X' && alternate.len() == 1 {
        return hgvsp(protein_id, codon_start as u64 / 3 + 1, reference, alternate[0], false);
    }
    hgvsp_inframe_indel_with_context(
        protein_id,
        codon_start as u64 / 3 + 1,
        codon_start as u64 / 3 + 1,
        std::str::from_utf8(&[reference]).ok()?,
        std::str::from_utf8(&alternate).ok()?,
        reference_peptide,
        strand,
        start_lost, Some(&duplication_peptide), full_reference_peptide,
)
}

/// Generate the VEP stop-loss form and count to the next translated stop.
///
/// VEP 115 `_stop_loss_extra_AA` measures non-frameshift extensions relative to
/// the full reference peptide length, even for an internal stop. Nonpositive
/// distances and residue-one events are reported as `extTer?`.
#[allow(clippy::too_many_arguments)]
pub fn hgvsp_stop_lost_from_cds(
    protein_id: &str,
    protein_pos: u64,
    reference_peptide: &str,
    alt_aa: u8,
    cds_and_downstream: &[u8],
    cds_start: Option<u64>,
    cds_end: Option<u64>,
    ref_allele: &Allele,
    alt_allele: &Allele,
    strand: Strand,
    codon_table: &CodonTable,
) -> Option<String> {
    let suffix = hgvsp_stop_lost_suffix_from_cds(
        protein_pos, reference_peptide, cds_and_downstream, cds_start, cds_end,
        ref_allele, alt_allele, strand, codon_table,
    )?;
    Some(format!("{}:p.Ter{}{}{}", protein_id, protein_pos,
        if alt_aa == b'-' { "del" } else { hgvs_aa_one_to_three(alt_aa) }, suffix))
}

/// VEP stop-loss suffix, shared by substitutions and complete deletion spans.
#[allow(clippy::too_many_arguments)]
pub fn hgvsp_stop_lost_suffix_from_cds(
    protein_pos: u64,
    reference_peptide: &str,
    cds_and_downstream: &[u8],
    cds_start: Option<u64>,
    cds_end: Option<u64>,
    ref_allele: &Allele,
    alt_allele: &Allele,
    strand: Strand,
    codon_table: &CodonTable,
) -> Option<String> {
    let (edited, _) = edited_cds(
        cds_and_downstream,
        cds_start,
        cds_end,
        ref_allele,
        alt_allele,
        strand,
    )?;
    let reference_length = reference_peptide
        .strip_suffix('*')
        .unwrap_or(reference_peptide)
        .len() as u64;
    let next_stop = codon_table
        .translate_seq(&edited)
        .iter()
        .position(|&aa| aa == b'*')
        .and_then(|zero_based| (zero_based as u64).checked_sub(reference_length))
        .filter(|&distance| distance > 0 && protein_pos > 1);
    let distance = next_stop
        .map(|value| value.to_string())
        .unwrap_or_else(|| "?".to_string());

    Some(format!("extTer{}", distance))
}

/// Recalculate a shifted in-frame insertion when HGVS normalization moves it
/// into or beyond the reference stop codon.
///
/// VEP 115.2 shifts the transcript allele before `hgvs_protein`, invalidates
/// the cached translation coordinates, and obtains new reference/alternate
/// peptides from the shifted allele. This helper mirrors that sequence only
/// for the narrow stop-retained case; ordinary in-frame indels continue to use
/// the consequence engine's peptide window.
#[allow(clippy::too_many_arguments)]
pub fn hgvsp_shifted_stop_retained_insertion(
    protein_id: &str,
    cds_and_downstream: &[u8],
    cds_start: Option<u64>,
    cds_end: Option<u64>,
    alt_allele: &Allele,
    strand: Strand,
    transcript_shift: u64,
    codon_table: &CodonTable,
    reference_peptide: Option<&[u8]>,
    reference_cds_length: Option<usize>,
    mapped_cds: Option<(Option<u64>, Option<u64>)>,
    full_reference_peptide: Option<&[u8]>,
) -> Option<String> {
    if transcript_shift == 0 {
        return None;
    }
    let (shifted_start, shifted_end) = mapped_cds.unwrap_or_else(|| (
        cds_start.and_then(|value| value.checked_add(transcript_shift)),
        cds_end.and_then(|value| value.checked_add(transcript_shift)),
    ));
    let mut shifted_alt = alt_allele.clone();
    if let Allele::Sequence(bases) = &mut shifted_alt {
        if bases.is_empty() {
            return None;
        }
        // VEP `shift_feature_seqs` rotates the allele with its 3' shift.
        let rotation = (transcript_shift % bases.len() as u64) as usize;
        match strand {
            Strand::Forward => bases.rotate_left(rotation),
            Strand::Reverse if transcript_shift <= bases.len() as u64 => bases.rotate_right(rotation),
            Strand::Reverse => {},
        }
    }
    let (edited, first) = edited_cds(
        cds_and_downstream,
        shifted_start,
        shifted_end,
        &Allele::Deletion,
        &shifted_alt,
        strand,
    )?;
    // At a codon boundary the insertion is between the final residue and the
    // stop. VEP keeps that as a pure insertion; it only recomputes a stop-codon
    // delins when the shifted insertion point lies inside the stop codon.
    if first % 3 == 0 {
        let inserted = match alt_allele {
            Allele::Sequence(bases) if bases.len() % 3 != 0 => bases.len(),
            _ => return None,
        };
        let mut alternate = codon_table.translate_seq(edited.get(first..first + inserted)?);
        if alternate != b"*" { alternate.push(b'X'); }
        let duplication_peptide = CodonTable::standard().translate_seq(cds_and_downstream);
        return hgvsp_inframe_indel_with_context(
            protein_id,
            first as u64 / 3 + 1,
            first as u64 / 3 + 1,
            "",
            std::str::from_utf8(&alternate).ok()?,
            reference_peptide,
            strand, false, Some(&duplication_peptide), full_reference_peptide,
);
    }
    let codon_start = first / 3 * 3;
    let reference_end = reference_cds_length.unwrap_or(cds_and_downstream.len()).min(cds_and_downstream.len());
    let reference_window = cds_and_downstream.get(codon_start..(codon_start + 3).min(reference_end))?;
    let reference = codon_table.translate_seq(reference_window);
    let reference = reference
        .first()
        .copied()
        .or_else(|| (!reference_window.is_empty() && reference_window.len() < 3).then_some(b'X'))?;
    // An internal stop can shift into a later ordinary codon. VEP 115
    // `hgvs_protein` recalculates both peptides before testing synonymy.
    let inserted = match alt_allele {
        Allele::Sequence(bases) => bases.len(),
        _ => return None,
    };
    // A whole-codon insertion immediately before the stop belongs to the
    // preceding peptide position and may be a duplication. Only a partial
    // codon makes the shifted allele's local peptide window begin at the stop.
    if inserted % 3 == 0 {
        return None;
    }
    // `shift_feature_seqs` translates the allele-local codon window, not the
    // remainder of the transcript. BioPerl represents its trailing incomplete
    // codon as X, which VEP's protein formatter subsequently spells `Ter`.
    let alternate_end = codon_start.checked_add(3 + inserted)?.min(edited.len());
    let alternate_window = edited.get(codon_start..alternate_end)?;
    let mut alternate = codon_table.translate_seq(alternate_window);
    if alternate_window.len() % 3 != 0 && alternate != b"*" {
        alternate.push(b'X');
    }
    // _clip_alleles retains these windows, then _get_hgvs_protein_type
    // reclassifies a longer alternate as delins. Formatting truncates that
    // alternate at its leading stop; only a sole stop remains synonymous.
    if reference == b'*' && alternate.first() == Some(&b'*') {
        let position = codon_start as u64 / 3 + 1;
        return if alternate.len() == 1 {
            hgvsp(protein_id, position, b'*', b'*', false)
        } else {
            Some(format!("{}:p.Ter{}delinsTer", protein_id, position))
        };
    }
    if reference == b'X' && alternate.len() == 1 {
        return hgvsp(
            protein_id,
            codon_start as u64 / 3 + 1,
            reference,
            alternate[0],
            false,
        );
    }
    let alternate = std::str::from_utf8(&alternate).ok()?;
    hgvsp_inframe_indel_with_context(
        protein_id,
        codon_start as u64 / 3 + 1,
        codon_start as u64 / 3 + 1,
        std::str::from_utf8(&[reference]).ok()?,
        alternate,
        reference_peptide,
        strand,
        false,
        None,
        full_reference_peptide,
    )
}

/// Transcript-oriented insertion offset, retaining Mapper::map_insert's
/// coding endpoint when the other endpoint is an intronic gap.
pub use fastvep_genome::insertion_point as cds_insertion_point;

fn edited_cds(
    cds_and_downstream: &[u8],
    cds_start: Option<u64>,
    cds_end: Option<u64>,
    ref_allele: &Allele,
    alt_allele: &Allele,
    strand: Strand,
) -> Option<(Vec<u8>, usize)> {
    let (first, ref_len) = if *ref_allele == Allele::Deletion {
        // The reference covers no bases, and Ensembl's zero-length interval puts
        // the insertion point just after the lower coordinate. An insertion on
        // an exon's edge has one end in the intron and so only one coordinate:
        // it still abuts the exonic base, on whichever side the strand puts the
        // intron. `cds_start` comes from the genomic left edge and `cds_end`
        // from the right, so the surviving one says which.
        let point = cds_insertion_point(cds_start, cds_end, strand)?;
        (point as usize, 0usize)
    } else {
        let (s, e) = (cds_start?, cds_end?);
        let (lo, hi) = (s.min(e), s.max(e));
        if lo < 1 || hi - lo + 1 > ref_allele.len() as u64 {
            return None;
        }
        // VEP _get_alternate_cds splices the mapped span, excluding introns.
        let ref_len = usize::try_from(hi - lo + 1).ok()?;
        ((lo - 1) as usize, ref_len)
    };
    if first + ref_len > cds_and_downstream.len() {
        return None;
    }
    let alt_cds: Vec<u8> = match alt_allele {
        Allele::Sequence(bases) => match strand {
            Strand::Forward => bases.clone(),
            Strand::Reverse => bases.iter().rev().map(|&b| complement(b)).collect(),
        },
        _ => Vec::new(),
    };

    let mut edited = Vec::with_capacity(cds_and_downstream.len() - ref_len + alt_cds.len());
    edited.extend_from_slice(&cds_and_downstream[..first]);
    edited.extend_from_slice(&alt_cds);
    edited.extend_from_slice(&cds_and_downstream[first + ref_len..]);
    Some((edited, first))
}

fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
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
/// A determinable stop distance is retained when the first changed residue is
/// 1. VEP 115.2 instead passes residue index 0 to a helper that rejects zero,
/// producing `fsTer?`; AnnoCAT records that behavior as a reviewed divergence.
pub fn hgvsp_frameshift(
    protein_id: &str,
    ref_translateable: &[u8],
    alt_translateable: &[u8],
    affected_codon_start: usize, // 0-based codon index where the frameshift starts
    codon_table: &CodonTable,
) -> Option<String> {
    hgvsp_frameshift_with_tables(
        protein_id,
        ref_translateable,
        alt_translateable,
        affected_codon_start,
        codon_table,
        codon_table,
        None,
        false,
        None,
    )
}

fn hgvsp_frameshift_with_tables(
    protein_id: &str,
    ref_translateable: &[u8],
    alt_translateable: &[u8],
    affected_codon_start: usize,
    reference_codon_table: &CodonTable,
    alternate_codon_table: &CodonTable,
    reference_peptide: Option<&[u8]>,
    stop_lost: bool,
    start_lost_end: Option<u64>,
) -> Option<String> {
    let prefix = format!("{}:p.", protein_id);

    // Translate both sequences from the affected codon onwards
    let ref_start = affected_codon_start * 3;
    // VEP appends a reference stop to the annotated peptide, including a
    // partial terminal codon. It can report deletion of that final residue.
    if ref_start + 3 > ref_translateable.len() && reference_peptide.is_none() {
        return None;
    }
    if ref_start > alt_translateable.len() {
        return None;
    }

    let mut ref_peptide: Vec<u8> = match reference_peptide {
        Some(peptide) => {
            let mut annotated = peptide.get(affected_codon_start..)?.to_vec();
            if !annotated.contains(&b'*') {
                annotated.push(b'*');
            }
            annotated
        }
        None => ref_translateable[ref_start..]
            .chunks(3)
            .filter(|c| c.len() == 3)
            .map(|c| reference_codon_table.translate(&[c[0], c[1], c[2]]))
            .collect(),
    };
    if affected_codon_start == 0 {
        reference_codon_table.normalize_reference_initiator(&mut ref_peptide, ref_translateable);
    }

    let alt_peptide: Vec<u8> = alt_translateable[ref_start..]
        .chunks(3)
        .filter(|c| c.len() == 3)
        .map(|c| alternate_codon_table.translate(&[c[0], c[1], c[2]]))
        .collect();
    // VEP 115.2's `_stop_loss_extra_AA` searches the complete alternate
    // translation and uses its first terminator, even when that terminator is
    // upstream of the frameshift. This matters for annotated selenoproteins:
    // their reference peptide contains U, while BioPerl translates the same
    // TGA as `*` in the alternate CDS. The resulting non-positive distance is
    // rendered `fsTer?` by VEP.
    let first_vep_stop = alt_translateable
        .chunks(3)
        .filter(|codon| codon.len() == 3)
        .map(|codon| alternate_codon_table.translate(&[codon[0], codon[1], codon[2]]))
        .position(|residue| residue == b'*')
        .map(|zero_based| zero_based + 1);

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

    // VEP applies start_lost after selecting the first changed full-peptide residue.
    if let Some(end) = start_lost_end {
        let reference = hgvs_aa_one_to_three(ref_aa);
        return Some(if first_changed_pos as u64 == end {
            format!("{}{}{}?", prefix, reference, first_changed_pos)
        } else {
            format!("{}{}{}_?{}", prefix, reference, first_changed_pos, end)
        });
    }

    // VEP's _get_fs_peptides changes the event to a deletion when the altered
    // translation ends before its starting residue. Its ordinary protein
    // formatter then turns deletion of the reference stop into delextTer?.
    if first_changed_offset >= alt_peptide.len() {
        // When the alternate translation contains the starting residue but
        // ends after a shared suffix, VEP's loop increments the position and
        // retains the last equal peptide pair. The `fs` type is unchanged, so
        // the formatter emits that pair at the incremented position with an
        // unknown stop distance.
        if let (Some(&reference), Some(&alternate)) = (
            ref_peptide.get(first_changed_offset.saturating_sub(1)),
            alt_peptide.last(),
        ) {
            if first_changed_offset > 0 && reference == alternate {
                return Some(format!(
                    "{}{}{}{}fsTer?",
                    prefix,
                    hgvs_aa_one_to_three(reference),
                    first_changed_pos,
                    hgvs_aa_one_to_three(alternate)
                ));
            }
        }
        return if ref_aa == b'*' && stop_lost {
            Some(format!("{}Ter{}delextTer?", prefix, first_changed_pos))
        } else {
            Some(format!(
                "{}{}{}del",
                prefix,
                hgvs_aa_one_to_three(ref_aa),
                first_changed_pos
            ))
        };
    }
    let alt_aa = if first_changed_offset < alt_peptide.len() {
        alt_peptide[first_changed_offset]
    } else {
        b'X'
    };

    let ref_aa3 = hgvs_aa_one_to_three(ref_aa);
    let alt_aa3 = hgvs_aa_one_to_three(alt_aa);

    if ref_aa == b'*' && alt_aa == b'*' {
        return Some(format!("{}Ter{}=", prefix, first_changed_pos));
    }

    // A frameshift whose *first* changed residue is already a terminator is
    // described as the nonsense variant it is: `p.Leu1545Ter`, not
    // `p.Leu1545TerfsTer1`. There is no shifted reading frame to describe -
    // translation stops at the residue the change lands on. 1,611 rows over a
    // 6,600-variant ClinVar sample.
    if alt_aa == b'*' {
        return Some(format!("{}{}{}Ter", prefix, ref_aa3, first_changed_pos));
    }

    let stop_dist = first_vep_stop
        .filter(|&stop| stop >= first_changed_pos)
        .map(|stop| stop - first_changed_pos + 1);

    if let Some(d) = stop_dist {
        Some(format!(
            "{}{}{}{}fsTer{}",
            prefix, ref_aa3, first_changed_pos, alt_aa3, d
        ))
    } else {
        // VEP uses `?` when no positive distance can be calculated, including
        // no stop, an upstream stop, or an immediate stop handled above.
        Some(format!(
            "{}{}{}{}fsTer?",
            prefix, ref_aa3, first_changed_pos, alt_aa3
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn depleted_cds_keeps_short_utr_and_partial_reference_windows_distinct() {
        let table = CodonTable::standard();
        for strand in [Strand::Forward, Strand::Reverse] {
            for (sequence, cds_length, lo, reference, peptide, expected) in [
                (b"ATGCGCTAACG".as_slice(), 9, 3, "GCGCTAA", b"MR*".as_slice(), "P:p.MetArgTer1_?3"),
                (b"ATGGCCGGTGAAATAA".as_slice(), 7, 1, "ATGGCCG", b"MA".as_slice(), "P:p.Met1_?2"),
            ] {
                let reference = Allele::Sequence(if strand == Strand::Reverse {
                    fastvep_genome::codon::reverse_complement(reference.as_bytes())
                } else { reference.as_bytes().to_vec() });
                assert_eq!(hgvsp_frameshift_from_cds_with_context(
                    "P", sequence, Some(lo), Some(cds_length), &reference,
                    &Allele::Deletion, strand, &table, &table, Some(peptide),
                    true, true, Some(cds_length as usize),
                ).as_deref(), Some(expected), "{strand:?}, CDS length {cds_length}");
            }
        }
    }

    #[test]
    fn depleted_cds_is_trimmed_before_utr_translation() {
        let table = CodonTable::standard();
        for strand in [Strand::Forward, Strand::Reverse] {
            let reference = Allele::from_str(if strand == Strand::Forward { "TGCGCTAA" } else { "TTAGCGCA" });
            assert_eq!(hgvsp_frameshift_from_cds_with_context(
                "P", b"ATGCGCTAATGCCAACTGA", Some(2), Some(9), &reference,
                &Allele::Deletion, strand, &table, &table, Some(b"MR*"),
                true, true, Some(9),
            ).as_deref(), Some("P:p.Met1_?3"));
        }
    }

    #[test]
    fn start_loss_keeps_the_clipped_partial_codon_window() {
        let table = CodonTable::standard();
        let cds = b"ATGGCCCACCGGA";
        assert_eq!(hgvsp_inframe_deletion_from_cds(
            "P", cds, 2, 13, &Allele::from_str("TGGCCCACCGGA"), Strand::Forward,
            &table, Some(b"MAHR"), true, Some(13), Some(b"MAHR"),
        ).as_deref(), Some("P:p.MetAlaHisArgTer1_?5"));
        assert_eq!(hgvsp_frameshift_from_cds_with_tables_and_ref_peptide(
            "P", cds, Some(3), Some(13), &Allele::from_str("GGCCCACCGGA"),
            &Allele::from_str("CCGGGTGGCCTT"), Strand::Forward,
            &table, &table, Some(b"MAHR"), false, true,
        ).as_deref(), Some("P:p.Met1_?4"));
    }

    #[test]
    fn insertion_clips_raw_peptide_but_names_full_reference_flanks() {
        assert_eq!(hgvsp_inframe_indel_with_context(
            "P", 1, 1, "L", "LGL", Some(b"LTQQ"), Strand::Forward,
            false, Some(b"LTQQ"), Some(b"MTQQ"),
        ), Some("P:p.Met1_Thr2insGlyLeu".into()));
    }

    #[test]
    fn deletion_rotates_local_residues_without_replacing_the_initiator() {
        for strand in [Strand::Forward, Strand::Reverse] {
            for (reference, local, full, expected) in [
                ("L", b"LTQQ".as_slice(), b"MTQQ".as_slice(), "P:p.Leu1del"),
                ("L", b"LLLQ".as_slice(), b"MLLQ".as_slice(), "P:p.Leu3del"),
                ("LM", b"LMLMQR".as_slice(), b"MMLMQR".as_slice(), "P:p.Leu3_Met4del"),
            ] {
                assert_eq!(hgvsp_inframe_indel_with_context(
                    "P", 1, reference.len() as u64, reference, "-", Some(local),
                    strand, false, None, Some(full),
                ).as_deref(), Some(expected));
            }
        }
    }
    use fastvep_genome::mitochondrial_codon_table;

    #[test]
    fn reference_terminal_codon_does_not_borrow_utr_bases() {
        let table = CodonTable::standard();
        assert_eq!(hgvsp_inframe_insertion_from_cds_with_start_lost(
            "P", b"ATGGAG", 4, 5, &Allele::from_str("C"), Strand::Reverse,
            1, &table, Some(b"M"), false, Some(5), None,
), None);
        assert_eq!(hgvsp_inframe_insertion_from_cds_with_start_lost(
            "P", b"ATGGAG", 4, 5, &Allele::from_str("GCGGGGACT"), Strand::Forward,
            10, &table, Some(b"M"), false, Some(5), None,
).as_deref(), Some("P:p.Ter2delinsAlaGlyThrGlu"));
    }

    #[test]
    fn stop_insertion_can_duplicate_a_partial_reference_residue() {
        assert_eq!(hgvsp_inframe_indel_with_context(
            "P", 2, 2, "W", "*W", Some(b"XW"), Strand::Forward,
            false, Some(b"XW"), None,
).as_deref(), Some("P:p.Xaa1dup"));
    }

    #[test]
    fn initiation_edits_do_not_replace_the_duplication_reference() {
        assert_eq!(hgvsp_inframe_indel_with_context(
            "P", 4, 4, "A", "AVESA", Some(b"MESAIA"), Strand::Forward,
            false, Some(b"VESAIA"), None,
).as_deref(), Some("P:p.Val1_Ala4dup"));
        assert_eq!(hgvsp_inframe_insertion_from_cds_with_start_lost(
            "P", b"GTGGAGAGTGCGATT", 4, 5, &Allele::from_str("CCC"), Strand::Reverse,
            2, &CodonTable::standard(), Some(b"MESAI"), true, None, None,
).as_deref(), Some("P:p.MetGlu2_?1"));
    }

    #[test]
    fn terminal_insertion_post_sequence_requires_a_following_residue() {
        assert_eq!(hgvsp_inframe_indel("P", 2, 2, "L", "LP", Some(b"MLP"), Strand::Forward).as_deref(), Some("P:p.Leu2_Pro3insPro"));
        assert_eq!(hgvsp_inframe_indel("P", 2, 2, "L", "LL", Some(b"MLL"), Strand::Forward).as_deref(), Some("P:p.Leu2dup"));
    }

    #[test]
    fn explicit_peptide_range_keeps_its_start_in_a_reverse_repeat() {
        assert_eq!(hgvsp_inframe_indel("P", 3, 4, "PP", "P", Some(b"MPPP"), Strand::Reverse).as_deref(), Some("P:p.Pro4del"));
    }

    #[test]
    fn terminal_windows_follow_vep_clipping_and_stop_loss_formatting() {
        let table = CodonTable::standard();
        assert_eq!(hgvsp_frameshift_from_cds_with_tables_and_ref_peptide(
            "P", b"ATGATTTAA", Some(6), Some(7), &Allele::from_str("TT"),
            &Allele::Deletion, Strand::Forward, &table, &table, Some(b"MI*"), true, false,
        ).as_deref(), Some("P:p.Ile3IlefsTer?"));
        assert_eq!(hgvsp_inframe_deletion_from_cds(
            "P", b"ATGAGTAGT", 6, 8, &Allele::from_str("TAG"), Strand::Forward,
            &table, Some(b"MS"), false, Some(8), None,
).as_deref(), Some("P:p.Ter3del"));
        assert_eq!(hgvsp_inframe_indel("P", 2, 3, "L*", "L", Some(b"ML"), Strand::Forward).as_deref(), Some("P:p.Ter3del"));
        assert_eq!(hgvsp_inframe_indel("P", 2, 3, "LX", "X", Some(b"ML"), Strand::Forward).as_deref(), Some("P:p.Leu2del"));
        assert_eq!(hgvsp_inframe_insertion_from_cds(
            "P", b"ATGCTTTGA", 7, 8, &Allele::from_str("AAA"), Strand::Reverse,
            1, &table, Some(b"ML*"),
        ).as_deref(), Some("P:p.Leu2_Ter3insPhe"));
        assert_eq!(hgvsp_frameshift_from_cds_with_tables_and_ref_peptide(
            "P", b"ATGATTTAA", Some(6), Some(9), &Allele::from_str("TTAA"),
            &Allele::Deletion, Strand::Forward, &table, &table, Some(b"MI*"), true, false,
        ).as_deref(), Some("P:p.Ile2_Ter3delextTer?"));
    }

    #[test]
    fn partial_stop_window_keeps_both_translation_endpoints() {
        for strand in [Strand::Forward, Strand::Reverse] {
            assert_eq!(
                hgvsp_inframe_indel("P", 331, 332, "*", "CX", Some(b"M"), strand).as_deref(),
                Some("P:p.Ter331_Ter332delinsCysTer")
            );
        }
    }

    #[test]
    fn recreated_stop_keeps_the_original_window_before_clipping() {
        for strand in [Strand::Forward, Strand::Reverse] {
            assert_eq!(
                hgvsp_inframe_indel("P", 2, 3, "A*", "A*X", Some(b"MA*"), strand).as_deref(),
                Some("P:p.Ala2_Ter3delinsAlaTer")
            );
        }
    }

    #[test]
    fn shifted_insertion_completes_a_partial_terminal_codon() {
        for (strand, inserted) in [(Strand::Forward, "G"), (Strand::Reverse, "C")] {
            assert_eq!(
                hgvsp_inframe_insertion_from_cds(
                    "P", b"ATGGC", 4, 5, &Allele::from_str(inserted), strand,
                    1, &CodonTable::standard(), Some(b"M"),
                ).as_deref(),
                Some("P:p.Ter2Gly")
            );
        }
    }

    #[test]
    fn deletion_shift_does_not_start_at_the_final_reference_residue() {
        for (strand, start, end) in [(Strand::Forward, 3, 4), (Strand::Reverse, 4, 3)] {
            assert_eq!(
                hgvsp_inframe_indel("P", start, end, "RR", "R", Some(b"MARRR"), strand).as_deref(),
                Some("P:p.Arg4del")
            );
        }
    }

    #[test]
    fn terminal_frameshift_keeps_the_initial_clipped_window() {
        let table = CodonTable::standard();
        for (cds, start, end, reference, alternate, peptide, expected) in [
            ("ATGAAGCAGTATTTCT", 6, 16, "GCAGTATTTCT", "-", "MKQYF", "P:p.Lys2_Phe5del"),
            ("ATGTTTTTTG", 10, 9, "-", "T", "MFF", "P:p.4="),
            ("ATGCAAAG", 7, 7, "A", "-", "MQ", "P:p.Ter3del"),
        ] {
            assert_eq!(hgvsp_frameshift_from_cds_with_tables_and_ref_peptide(
                "P", cds.as_bytes(), Some(start), Some(end),
                &Allele::from_str(reference), &Allele::from_str(alternate),
                Strand::Forward, &table, &table, Some(peptide.as_bytes()), false, false,
            ).as_deref(), Some(expected));
        }
    }

    #[test]
    fn terminal_frameshift_deletion_uses_the_shifted_start() {
        let table = CodonTable::standard();
        for (cds, peptide, expected) in [
            (b"ATGAGTAGT".as_slice(), b"MSS".as_slice(), "P:p.Ser3del"),
            (b"ATGAGTC".as_slice(), b"MS".as_slice(), "P:p.Ter3del"),
        ] {
            assert_eq!(
                hgvsp_frameshift_from_cds_with_tables_and_ref_peptide(
                    "P",
                    cds,
                    Some(7),
                    Some(7),
                    &Allele::from_str("A"),
                    &Allele::Deletion,
                    Strand::Forward,
                    &table,
                    &table,
                    Some(peptide),
                    false,
        false,
    )
                .as_deref(),
                Some(expected)
            );
        }
    }

    #[test]
    fn residue_one_frameshift_keeps_a_determinable_stop_distance() {
        let result = hgvsp_frameshift(
            "ENSP1",
            b"CGTCGTCGT",
            b"CCCTGATAA",
            0,
            &CodonTable::standard(),
        );
        assert_eq!(result, Some("ENSP1:p.Arg1ProfsTer2".to_string()));
    }

    #[test]
    fn vep_uses_an_upstream_alternate_stop_for_frameshift_distance() {
        let table = CodonTable::standard();
        let result = hgvsp_frameshift_with_tables(
            "ENSP1",
            b"ATGTGAAAAAAATAA", // annotated peptide M U K K *
            b"ATGTGAAAACCTTAA", // BioPerl-style translation M * K P *
            3,
            &table,
            &table,
            Some(b"MUKK*"),
            false,
        None,
    );
        assert_eq!(result, Some("ENSP1:p.Lys4ProfsTer?".to_string()));
    }

    #[test]
    fn vep_retains_the_last_equal_residue_when_the_alternate_translation_ends() {
        let table = CodonTable::standard();
        let result = hgvsp_frameshift_with_tables(
            "ENSP1",
            b"ATGAAAAAA", // M K K
            b"ATGAAA",    // M K
            1,
            &table,
            &table,
            Some(b"MKK"),
            false,
        None,
    );
        assert_eq!(result, Some("ENSP1:p.Lys3LysfsTer?".to_string()));
    }

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

        let mixed_result = hgvsp_frameshift_with_tables(
            "ENSP1",
            b"AGACCCCCCCCC",
            b"ACATGAAAATAA",
            0,
            &mitochondrial,
            &standard,
            None,
            false,
        None,
    );
        assert_eq!(mixed_result, Some("ENSP1:p.Ter1ThrfsTer2".to_string()));
    }

    #[test]
    fn frameshift_uses_annotated_reference_peptide_and_vep_terminal_deletions() {
        let table = CodonTable::standard();

        assert_eq!(
            hgvsp_frameshift_with_tables(
                "ENSP1",
                b"CTGAAAAAA",
                b"CCCTGATAA",
                0,
                &table,
                &table,
                Some(b"MKK"),
                false,
        None,
    ),
            Some("ENSP1:p.Met1ProfsTer2".to_string())
        );
        assert_eq!(
            hgvsp_frameshift_with_tables(
                "ENSP1",
                b"GGT",
                b"",
                0,
                &table,
                &table,
                Some(b"G"),
                false,
        None,
    ),
            Some("ENSP1:p.Gly1del".to_string())
        );
        assert_eq!(
            hgvsp_frameshift_with_tables("ENSP1", b"TAA", b"", 0, &table, &table, Some(b"*"), true,
        None,
    ),
            Some("ENSP1:p.Ter1delextTer?".to_string())
        );
    }

    #[test]
    fn shifted_insertion_recomputes_the_terminal_codon_window() {
        assert_eq!(hgvsp_shifted_stop_retained_insertion(
            "P", b"ATGTAAAAATGA", Some(4), Some(3), &Allele::from_str("AAAAA"),
            Strand::Forward, 1, &CodonTable::standard(), Some(b"M*"),
         None,  None, None,
), Some("P:p.Ter2delinsTer".into()));
        assert_eq!(
            hgvsp_shifted_stop_retained_insertion(
                "P",
                b"TAAGCT",
                Some(1),
                Some(2),
                &Allele::Sequence(b"A".to_vec()),
                Strand::Forward,
                1,
                &CodonTable::standard(),
                Some(b"*A"),
             None,  None, None,
),
            Some("P:p.Ter1=".to_string()),
        );
        for (strand, inserted) in [
            (Strand::Forward, b"AGAGTTAGAT".as_slice()),
            (Strand::Reverse, b"ATCTAACTCT".as_slice()),
        ] {
            assert_eq!(
                hgvsp_shifted_stop_retained_insertion(
                    "P",
                    b"GAACGT",
                    Some(1),
                    Some(2),
                    &Allele::Sequence(inserted.to_vec()),
                    strand,
                    1,
                    &CodonTable::standard(),
                    Some(b"ER"),
                 None,  None, None,
),
                Some("P:p.Glu1_Arg2insSerTer".to_string()),
                "the peptide window uses the rotated insertion on both strands"
            );
        }
        assert_eq!(
            hgvsp_shifted_stop_retained_insertion(
                "P",
                b"TGATCT",
                Some(2),
                Some(3),
                &Allele::Sequence(b"A".to_vec()),
                Strand::Forward,
                1,
                &CodonTable::standard(),
                Some(b"*S"),
             None,  None, None,
),
            Some("P:p.Ter1_Ser2insTer".to_string()),
        );
        assert_eq!(
            hgvsp_shifted_stop_retained_insertion(
                "P",
                b"TGAAAACAT",
                Some(2),
                Some(3),
                &Allele::Sequence(b"A".to_vec()),
                Strand::Forward,
                3,
                &CodonTable::standard(),
                Some(b"*KH"),
             None,  None, None,
),
            Some("P:p.Lys2_His3insTer".to_string()),
        );
        assert_eq!(
            hgvsp_shifted_stop_retained_insertion(
                "P",
                b"TGAAGA",
                Some(2),
                Some(3),
                &Allele::Sequence(b"A".to_vec()),
                Strand::Forward,
                2,
                &CodonTable::standard(),
                None,
             None,  None, None,
),
            Some("P:p.Arg2delinsLysTer".to_string()),
            "HGVS shifts beyond an internal stop before comparing peptides"
        );
        let mut cds = vec![b'G'; 138];
        cds.extend_from_slice(b"TAGCTA");
        assert_eq!(
            hgvsp_shifted_stop_retained_insertion(
                "ENSP1",
                &cds,
                Some(138),
                Some(137),
                &Allele::Sequence(b"T".to_vec()),
                Strand::Forward,
                2,
                &CodonTable::standard(),
                None,
             None,  None, None,
),
            Some("ENSP1:p.Ter47delinsLeuTer".to_string())
        );
        assert_eq!(
            hgvsp_shifted_stop_retained_insertion(
                "ENSP1",
                &cds,
                Some(135),
                Some(136),
                &Allele::Sequence(b"TAG".to_vec()),
                Strand::Forward,
                3,
                &CodonTable::standard(),
                None,
             None,  None, None,
),
            None,
            "a whole-codon insertion before the stop remains eligible for duplication"
        );
        assert_eq!(
            hgvsp_shifted_stop_retained_insertion(
                "ENSP1",
                b"GGGTAG",
                Some(2),
                Some(3),
                &Allele::Sequence(b"T".to_vec()),
                Strand::Forward,
                1,
                &CodonTable::standard(),
                None,
             None,  None, None,
),
            None,
            "an insertion between the last residue and stop remains a pure insertion"
        );
        assert_eq!(
            hgvsp_shifted_stop_retained_insertion(
                "ENSP1",
                b"GGGGG",
                Some(2),
                Some(3),
                &Allele::Sequence(b"G".to_vec()),
                Strand::Forward,
                2,
                &CodonTable::standard(),
                None,
             None,  None, None,
),
            Some("ENSP1:p.Ter2Gly".to_string()),
            "VEP spells a completed cds_end_NF codon as a Ter substitution"
        );
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
    fn test_hgvsp_synonymous() {
        let result = hgvsp("ENSP00000001", 41, b'R', b'R', false);
        assert_eq!(result, Some("ENSP00000001:p.Arg41=".to_string()));
    }

    #[test]
    fn test_hgvsp_uses_vep_ter_spelling_for_x() {
        let result = hgvsp("ENSP00000001", 326, b'X', b'V', false);
        assert_eq!(result, Some("ENSP00000001:p.Ter326Val".to_string()));
        assert_eq!(hgvsp("P", 3, b'X', b'*', false), Some("P:p.Ter3=".into()));
        for strand in [Strand::Forward, Strand::Reverse] {
            assert_eq!(hgvsp_inframe_indel("P", 3, 3, "X", "IX", Some(b"MY"), strand), None);
            assert_eq!(hgvsp_inframe_indel("P", 3, 3, "X", "VX", Some(b"MV"), strand), Some("P:p.Val2dup".into()));
        }
    }

    #[test]
    fn shifted_frameshift_start_loss_uses_the_changed_full_peptide_residue() {
        let table = CodonTable::standard();
        for strand in [Strand::Forward, Strand::Reverse] {
            let reference = Allele::from_str(if strand == Strand::Forward { "G" } else { "C" });
            assert_eq!(
                hgvsp_frameshift_from_cds_with_tables_and_ref_peptide(
                    "P", b"ATGGAAAAATAA", Some(4), Some(4), &reference,
                    &Allele::Deletion, strand, &table, &table, Some(b"MEK*"), false, true,
                ),
                Some("P:p.Glu2?".into()),
            );
        }
    }

    #[test]
    fn test_hgvsp_inframe_deletion_single() {
        // single-residue in-frame deletion
        let r = hgvsp_inframe_indel("ENSP00000001", 157, 157, "F", "-", None, Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Phe157del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_deletion_range() {
        // multi-residue in-frame deletion (regression for the p.Tyr43??? bug)
        let r = hgvsp_inframe_indel("ENSP00000001", 43, 43, "YXQ", "-", None, Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Tyr43_Gln45del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_delins() {
        // in-frame deletion-insertion
        let r = hgvsp_inframe_indel("ENSP00000001", 2173, 2173, "NL", "K", None, Strand::Forward);
        assert_eq!(
            r,
            Some("ENSP00000001:p.Asn2173_Leu2174delinsLys".to_string())
        );
    }

    #[test]
    fn protein_terminator_is_trimmed_after_clipping_like_vep() {
        let peptide = b"MVWQ";
        assert_eq!(
            hgvsp_inframe_indel(
                "ENSP00000001",
                3,
                3,
                "W",
                "*W",
                Some(peptide),
                Strand::Forward,
            ),
            Some("ENSP00000001:p.Val2_Trp3insTer".to_string())
        );

        let peptide = b"MAYQ";
        assert_eq!(
            hgvsp_inframe_indel(
                "ENSP00000001",
                3,
                3,
                "Y",
                "*H",
                Some(peptide),
                Strand::Forward,
            ),
            Some("ENSP00000001:p.Tyr3delinsTer".to_string())
        );

        let peptide = b"MALA";
        assert_eq!(
            hgvsp_inframe_indel(
                "ENSP00000001",
                3,
                3,
                "L",
                "LF*Q",
                Some(peptide),
                Strand::Forward,
            ),
            Some("ENSP00000001:p.Leu3_Ala4insPheTer".to_string())
        );
    }

    #[test]
    fn shifted_inframe_deletion_is_clipped_from_translated_peptide_tails() {
        let cds = b"ATGGGTGGTGGTGCTGCTTAA";
        assert_eq!(
            hgvsp_inframe_deletion_from_cds(
                "ENSP1",
                cds,
                5,
                10,
                &Allele::Sequence(b"GTGGTG".to_vec()),
                Strand::Forward,
                &CodonTable::standard(),
                Some(b"MGGGAA*"),
                false,
             None, None,
),
            Some("ENSP1:p.Gly3_Gly4del".to_string())
        );
    }

    #[test]
    fn shifted_inframe_insertion_rotates_before_translating_the_peptide_tail() {
        for (inserted, expected) in [("AAA", "Phe0_Phe1insLys"), ("ATG", "Met1dup")] {
            assert_eq!(hgvsp_inframe_insertion_from_cds(
                "ENSP1", b"ATGTCTTTC", 2, 1, &Allele::from_str(inserted),
                Strand::Forward, 1, &CodonTable::standard(), Some(b"MSF"),
            ), Some(format!("ENSP1:p.{expected}")));
        }
        assert_eq!(hgvsp_inframe_insertion_from_cds(
            "ENSP1", b"ATGGTCAAATAA", 7, 6, &Allele::Sequence(b"GAC".to_vec()),
            Strand::Reverse, 7, &CodonTable::standard(), Some(b"MVK*"),
        ), Some("ENSP1:p.Val2dup".into()));
        for sequence in [b"ATGCATGGAAGTAGCTAA".as_slice(), b"ATGCATGGAAGT"] {
            assert_eq!(hgvsp_inframe_insertion_from_cds(
                "ENSP1", sequence, 12, 13, &Allele::Sequence(b"CATGGAAGT".to_vec()),
                Strand::Forward, 9, &CodonTable::standard(), Some(b"MHGS"),
            ), Some("ENSP1:p.His2_Ser4dup".into()));
        }
        assert_eq!(
            hgvsp_inframe_insertion_from_cds(
                "ENSP1",
                b"TGGTGTAAATAA", // W C K *
                3,
                4,
                &Allele::Sequence(b"TTCTGGTCT".to_vec()),
                Strand::Forward,
                6,
                &CodonTable::standard(),
                Some(b"WCK*"),
            ),
            Some("ENSP1:p.Trp1_Cys2insSerPheTrp".to_string())
        );
    }

    // In-frame insertions previously fell through to the substitution branch,
    // which compares only the first residue of each side. Every case below is
    // real fastVEP output from a clinical panel, checked against Ensembl VEP.

    /// Build a reference peptide with `residues` placed at 1-based `at`, padded
    /// with a filler residue that cannot be confused with the payload.
    fn peptide_with(at: u64, residues: &str, length: usize) -> Vec<u8> {
        let mut pep = vec![b'M'; length];
        for (i, b) in residues.bytes().enumerate() {
            pep[at as usize - 1 + i] = b;
        }
        pep
    }

    #[test]
    fn test_hgvsp_inframe_insertion_collapses_to_duplication() {
        // C8A c.553_554insGGA, amino_acids "W/WR" at 185 — previously p.Trp185=.
        // Residue 186 is already Arg, so inserting Arg duplicates it.
        let pep = peptide_with(185, "WRQ", 200);
        let r = hgvsp_inframe_indel(
            "ENSP00000001",
            185,
            185,
            "W",
            "WR",
            Some(&pep),
            Strand::Forward,
        );
        assert_eq!(r, Some("ENSP00000001:p.Arg186dup".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_insertion_multi_residue_duplication() {
        // FLT3 c.1759_1800dup, a 14-codon ITD — previously p.Asp600=. The delins
        // spelling runs to ~55 characters; the duplication form is 18.
        let dup = "NEYFYVDFREYEYD";
        let pep = peptide_with(587, &format!("{}K", dup), 700);
        let alt = format!("D{}", dup);
        let r = hgvsp_inframe_indel(
            "ENSP00000001",
            600,
            600,
            "D",
            &alt,
            Some(&pep),
            Strand::Forward,
        );
        assert_eq!(r, Some("ENSP00000001:p.Asn587_Asp600dup".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_insertion_true_insertion_uses_ins_form() {
        // ITPKB c.275_276insGGT, amino_acids "S/SG" at 92 — previously p.Ser92=.
        // Gly does not repeat the preceding residues, so it stays an insertion.
        let pep = peptide_with(92, "SSK", 200);
        let r = hgvsp_inframe_indel(
            "ENSP00000001",
            92,
            92,
            "S",
            "SG",
            Some(&pep),
            Strand::Forward,
        );
        assert_eq!(r, Some("ENSP00000001:p.Ser92_Ser93insGly".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_deletion_is_three_prime_shifted() {
        // A deletion inside a homopolymer run is reported at the most C-terminal
        // position it can occupy: deleting one Ala from AAA at 2..4 is p.Ala4del.
        let pep: Vec<u8> = "MAAAGK".bytes().collect();
        let r = hgvsp_inframe_indel("ENSP00000001", 2, 2, "A", "-", Some(&pep), Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Ala4del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_deletion_uses_vep_115_terminal_shift_bound() {
        // VEP 115.2's `_shift_3prime` leaves a two-residue deletion one residue
        // before the end because it requires a full two-residue comparison
        // window after the original deletion. Real case: NT5C2 c.1674_1679del.
        let pep: Vec<u8> = "MKEEEEE*".bytes().collect();
        let r = hgvsp_inframe_indel("ENSP00000001", 3, 3, "EE", "-", Some(&pep), Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Glu5_Glu6del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_indel_without_peptide_stays_valid() {
        for (start, end, strand) in [(2, 3, Strand::Forward), (3, 2, Strand::Reverse)] {
            assert_eq!(
                hgvsp_inframe_indel("P", start, end, "YX", "*", Some(b"MY"), strand),
                Some("P:p.Tyr2_Ter3delinsTer".into()),
            );
        }
        // No sequence context: emit an unshifted but well-formed description
        // rather than nothing, and never a substitution shape.
        let r = hgvsp_inframe_indel("ENSP00000001", 185, 185, "W", "WR", None, Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Trp185delinsTrpArg".to_string()));
        let d = hgvsp_inframe_indel("ENSP00000001", 157, 157, "F", "-", None, Strand::Forward);
        assert_eq!(d, Some("ENSP00000001:p.Phe157del".to_string()));
    }

    // The peptide is caller-derived, so every index into it has to survive a
    // transcript whose sequence disagrees with its own coordinates. Each case
    // below aborted the process before the bounds checks were added.

    /// `(case name, protein_start, ref_aas, alt_aas, peptide, expected)`.
    type UnusablePeptideCase<'a> = (&'a str, u64, &'a str, &'a str, &'a [u8], Option<&'a str>);

    #[test]
    fn test_hgvsp_inframe_indel_survives_unusable_peptides() {
        assert_eq!(
            hgvsp_inframe_indel("P", 2, 3, "S*", "*", Some(b"MS*"), Strand::Forward),
            Some("P:p.Ser2del".to_string()),
        );
        assert_eq!(
            hgvsp_inframe_indel("P", 3, 2, "P*", "R", Some(b"MP*"), Strand::Reverse),
            Some("P:p.Pro2_Ter3delinsArg".to_string()),
        );
        assert_eq!(
            hgvsp_inframe_indel("P", 2, 2, "*", "*CKX", Some(b"M*"), Strand::Forward),
            Some("P:p.Ter2delinsTer".to_string()),
        );
        let short: Vec<u8> = "MAAAGK".bytes().collect();
        let cases: Vec<UnusablePeptideCase<'_>> = vec![
            // Deleted block overruns the peptide end.
            (
                "block overruns end",
                6,
                "KX",
                "-",
                &short,
                Some("ENSP00000001:p.Lys6_Ter7del"),
            ),
            // protein_start past the peptide entirely.
            (
                "start beyond peptide",
                100,
                "AK",
                "-",
                &short,
                Some("ENSP00000001:p.Ala100_Lys101del"),
            ),
            // Insertion anchored past the peptide.
            ("insertion beyond peptide", 500, "-", "R", &short, None),
            // Empty peptide (truncated transcript).
            (
                "empty peptide",
                1,
                "F",
                "-",
                &[],
                Some("ENSP00000001:p.Phe1del"),
            ),
            // Insertion one residue past the C-terminus.
            (
                "insertion at C-terminus",
                7,
                "K",
                "KG",
                &short,
                Some("ENSP00000001:p.Lys7delinsLysGly"),
            ),
            // protein_start of 0 would underflow a 1-based conversion.
            (
                "zero protein_start",
                0,
                "F",
                "-",
                &short,
                Some("ENSP00000001:p.Phe0del"),
            ),
        ];
        for (name, start, reference, alternate, pep, expected) in cases {
            let got = hgvsp_inframe_indel(
                "ENSP00000001",
                start,
                start,
                reference,
                alternate,
                Some(pep),
                Strand::Forward,
            );
            assert_eq!(got.as_deref(), expected, "case: {name}");
        }
    }

    #[test]
    fn test_hgvsp_inframe_indel_does_not_shift_onto_the_terminator() {
        // Transcript::peptide ends with `*`. A change abutting the stop must not
        // shift onto it or name it as a flanking residue: Ter is not a residue
        // of the protein, and a position at or past it is not a real position.
        let pep: Vec<u8> = "MKKG*".bytes().collect();

        // Deleting one Lys from the KK run shifts to the 3'-most Lys (3), not
        // onto Gly4 or the terminator.
        let deletion =
            hgvsp_inframe_indel("ENSP00000001", 2, 2, "K", "-", Some(&pep), Strand::Forward);
        assert_eq!(deletion, Some("ENSP00000001:p.Lys3del".to_string()));

        // An insertion immediately before the terminator has no residue on its
        // 3' side once the stop is excluded, so VEP emits no HGVSp rather than
        // a Ter-flanked range.
        let insertion =
            hgvsp_inframe_indel("ENSP00000001", 4, 4, "G", "GS", Some(&pep), Strand::Forward);
        assert_eq!(insertion, None);
        assert!(!deletion.unwrap().contains("Ter"));
    }

    #[test]
    fn test_hgvsp_inframe_indel_ignores_a_peptide_that_disagrees() {
        // The peptide says Gln at 2; the caller says Phe. Trusting the peptide
        // would emit a confident, well-formed, wrong description (p.Gln6del).
        // Fall back to the caller's residues instead.
        let pep: Vec<u8> = "MQQQQQK".bytes().collect();
        let r = hgvsp_inframe_indel("ENSP00000001", 2, 2, "F", "-", Some(&pep), Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Phe2del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_indel_reads_the_span_from_either_end() {
        // For a shrinking change the call sites do not always pass the start of
        // `ref_aas`: real ClinVar output has `Protein_position` 328-329 with
        // `Amino_acids` FF/F arriving here as protein_start 329, the end of the
        // pair rather than its start. Both numbers name the same span, so the
        // description must not depend on which end the caller happened to send.
        //
        // The peptide is M R I F F A S M, so the F pair is at 4-5. Which end
        // the caller sends follows from the strand (#96), so the two readings
        // are paired with the strand that produces them - but the answer has to
        // be the same either way, which is what this test is for.
        let pep: Vec<u8> = "MRIFFASM".bytes().collect();
        for (anchor, strand) in [(4u64, Strand::Forward), (5, Strand::Reverse)] {
            assert_eq!(
                hgvsp_inframe_indel(
                    "ENSP00000001",
                    anchor,
                    anchor,
                    "FF",
                    "F",
                    Some(&pep),
                    strand
                ),
                Some("ENSP00000001:p.Phe5del".to_string()),
                "anchor {anchor} on {strand:?} should describe the same deletion"
            );
        }

        // And the pairing is a preference, not a requirement: a reverse-strand
        // caller whose residues really do sit at `protein_start` still
        // normalises, because the other end stays as a fallback.
        assert_eq!(
            hgvsp_inframe_indel("ENSP00000001", 4, 4, "FF", "F", Some(&pep), Strand::Reverse),
            Some("ENSP00000001:p.Phe5del".to_string()),
            "the unpreferred end must still be tried"
        );

        // This case used to take the un-normalised path and emit
        // `p.Phe5_Phe6delinsPhe`, naming residue 6 as Phe when the peptide has
        // Ala there. Reading the span from its other end is what fixes it: the
        // letters and the numbers now agree with the protein.
    }

    #[test]
    fn start_loss_uses_vep_peptide_preparation_and_explicit_predicate() {
        assert_eq!(
            hgvsp_start_lost("P", 1, "L", "L", Some(b"LKGN")),
            Some("P:p.Leu1?".into())
        );
        assert_eq!(
            hgvsp_start_lost("P", 1, "M", "NM", Some(b"MKGN")),
            Some("P:p.Asn1_?0".into())
        );
        assert_eq!(
            hgvsp_start_lost("P", 1, "M", "MM", Some(b"MKGN")),
            Some("P:p.Met1dup".into())
        );
        for start_lost in [false, true] {
            assert_eq!(
                hgvsp_inframe_deletion_from_cds(
                    "P", b"ATGGCCGGGGCCATCAAATAA", 3, 14,
                    &Allele::Sequence(b"GGCCGGGGCCAT".to_vec()), Strand::Forward,
                    &CodonTable::standard(), Some(b"MAGAIK*"), start_lost,
                 None, None,
),
                Some(if start_lost { "P:p.MetAlaGlyAla1_?4" } else { "P:p.Met1_Ala4del" }.into())
            );
        }
    }

    #[test]
    fn a_deletion_of_the_initiation_codon_is_described_as_unresolvable() {
        // Removing the start codon leaves the protein's fate undetermined: the
        // ribosome may initiate downstream, or not at all, and the sequence does
        // not say which. Describing it as an ordinary deletion asserts a protein
        // that begins where nothing says it begins.
        //
        // Both shapes come from real ClinVar rows checked against Ensembl VEP:
        // KCNA2 `MNII/I` at residues 1-4, and POLE `EA/A` at 1-2.
        assert_eq!(
            hgvsp_start_lost("ENSP00000491354", 1, "MNII", "I", None),
            Some("ENSP00000491354:p.MetAsnIle1_?3".to_string())
        );

        assert_eq!(
            hgvsp_start_lost("ENSP00000500921", 1, "EA", "A", None),
            Some("ENSP00000500921:p.Glu1?".to_string())
        );
    }

    #[test]
    fn a_change_that_loses_the_initiator_is_unresolvable() {
        // The marker is specific to the initiation codon. A deletion anywhere
        // else is an ordinary `del`, including one that starts at residue 2.
        let pep: Vec<u8> = "MKFFASM".bytes().collect();
        assert_eq!(
            hgvsp_inframe_indel("ENSP00000001", 2, 2, "K", "-", Some(&pep), Strand::Forward),
            Some("ENSP00000001:p.Lys2del".to_string())
        );

        // A delins that takes out the initiator is *not* an ordinary delins.
        // Replacing `MK` with `W` removes the ATG, so where translation begins
        // is exactly what is no longer known, and Ensembl marks it `?` for every
        // shape of change that loses the start - `p.Met1?` for `M/T`,
        // `p.MetAla1_?2` for `MA/IS`, `p.Pro2?` for `XP/XAKSTVGA`, all from real
        // VEP 115.1. Naming the replacement instead described a protein that may
        // never be made.
        assert_eq!(
            hgvsp_start_lost("ENSP00000001", 1, "MK", "W", None),
            Some("ENSP00000001:p.MetLys1_?2".to_string())
        );

        // VEP's start_lost formatter does not require the cached reference
        // peptide to corroborate residues already supplied by the consequence
        // calculation.
        assert_eq!(
            hgvsp_start_lost("ENSP00000001", 1, "MA", "IP", None),
            Some("ENSP00000001:p.MetAla1_?2".to_string())
        );

        // A replacement that keeps the reference residues at one end of itself
        // has not lost the start. `MK` -> `MWK` is an insertion between them,
        // and it is described as one.
        assert_eq!(
            hgvsp_inframe_indel(
                "ENSP00000001",
                1,
                1,
                "MK",
                "MWK",
                Some(&pep),
                Strand::Forward
            ),
            Some("ENSP00000001:p.Met1_Lys2insTrp".to_string())
        );
    }

    #[test]
    fn test_hgvsp_inframe_indel_still_declines_when_neither_end_corroborates() {
        // Reading from the other end is a consequence of the coordinate
        // convention, not a search. When the peptide carries the reference at
        // neither end of the span, there is no evidence for any anchor and the
        // un-normalised description is still the honest answer.
        let pep: Vec<u8> = "MRIFFASM".bytes().collect();
        for strand in [Strand::Forward, Strand::Reverse] {
            let r = hgvsp_inframe_indel("ENSP00000001", 5, 5, "KK", "K", Some(&pep), strand);
            assert_eq!(
                r,
                Some("ENSP00000001:p.Lys5_Lys6delinsLys".to_string()),
                "on {strand:?} an uncorroborated span must stay unshifted"
            );
        }

        // A single residue has only one candidate: there is no other end to try,
        // so an uncorroborated anchor cannot be rescued and must not be guessed.
        let r = hgvsp_inframe_indel("ENSP00000001", 2, 2, "W", "-", Some(&pep), Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Trp2del".to_string()));
    }

    #[test]
    fn anchor_candidates_offers_the_other_end_only_when_there_is_one() {
        // An empty or single reference has one reading, which keeps insertions -
        // where protein_start is already the far end of the pair - unchanged.
        assert_eq!(anchor_candidates(10, 0, false), [Some(10), None]);
        assert_eq!(anchor_candidates(10, 1, false), [Some(10), None]);
        assert_eq!(anchor_candidates(10, 2, false), [Some(10), Some(9)]);
        assert_eq!(anchor_candidates(10, 5, false), [Some(10), Some(6)]);
        // Never underflows past residue 1.
        assert_eq!(anchor_candidates(2, 5, false), [Some(2), None]);
        assert_eq!(anchor_candidates(1, 2, false), [Some(1), None]);
    }

    #[test]
    fn anchor_candidates_puts_the_determined_end_first() {
        // Same pair, opposite order: which one is tried first is what decides a
        // periodic reference, where both corroborate (#96).
        assert_eq!(anchor_candidates(10, 5, true), [Some(6), Some(10)]);
        assert_eq!(anchor_candidates(10, 5, false), [Some(10), Some(6)]);

        // With only one candidate there is nothing to order, and the request for
        // the other end must not invent one or lose the one there is.
        assert_eq!(anchor_candidates(10, 1, true), [Some(10), None]);
        assert_eq!(anchor_candidates(10, 0, true), [Some(10), None]);
        assert_eq!(anchor_candidates(2, 5, true), [Some(2), None]);
        assert_eq!(anchor_candidates(1, 2, true), [Some(1), None]);
    }

    #[test]
    fn only_a_shrinking_change_on_the_reverse_strand_is_anchored_at_its_end() {
        // The four combinations, because the rule is a reading of how the
        // coordinates were built (see `predict_coding_consequence`): the
        // residues come from the lower CDS coordinate for a shrinking change and
        // from `cds_start` otherwise, and only on the reverse strand are those
        // two different ends of the span.
        assert!(anchored_at_span_end(Strand::Reverse, 3, 0), "deletion");
        assert!(
            anchored_at_span_end(Strand::Reverse, 3, 1),
            "shrinking delins"
        );
        assert!(
            !anchored_at_span_end(Strand::Forward, 3, 0),
            "forward strand"
        );
        assert!(!anchored_at_span_end(Strand::Reverse, 1, 3), "insertion");
        assert!(
            !anchored_at_span_end(Strand::Reverse, 3, 3),
            "equal-length replacement does not shrink"
        );
    }

    #[test]
    fn a_periodic_reference_on_the_reverse_strand_names_the_span_the_caller_meant() {
        // Issue #96. `EGE` sits at residues 2-4 *and* at 4-6, so both anchors
        // are corroborated and taking whichever came first picked the wrong one.
        // These are not two spellings of one variant: deleting 2-4 leaves MGEA
        // and deleting 4-6 leaves MEGA, so the 3'-rule cannot merge them
        // afterwards - the wrong answer named residues the variant never touched.
        let pep: Vec<u8> = "MEGEGEA".bytes().collect();
        assert_eq!(
            hgvsp_inframe_indel(
                "ENSP00000001",
                4,
                4,
                "EGE",
                "-",
                Some(&pep),
                Strand::Reverse
            ),
            Some("ENSP00000001:p.Glu2_Glu4del".to_string())
        );

        // Period 3 at n = 4, the same shape one residue longer.
        let pep: Vec<u8> = "MABCABCA".bytes().collect();
        assert_eq!(
            hgvsp_inframe_indel(
                "ENSP00000001",
                5,
                5,
                "ABCA",
                "-",
                Some(&pep),
                Strand::Reverse
            ),
            Some("ENSP00000001:p.Ala2_Ala5del".to_string())
        );
    }

    #[test]
    fn a_periodic_reference_on_the_forward_strand_is_read_from_the_start() {
        // The mirror of the case above, and the reason the fix is an ordering
        // rather than a preference for the earlier residue: on the forward strand
        // `protein_start` *is* the first affected residue, so the same peptide
        // and the same reference must resolve to the other span.
        let pep: Vec<u8> = "MEGEGEA".bytes().collect();
        assert_eq!(
            hgvsp_inframe_indel(
                "ENSP00000001",
                4,
                4,
                "EGE",
                "-",
                Some(&pep),
                Strand::Forward
            ),
            Some("ENSP00000001:p.Glu4_Glu6del".to_string())
        );
        // And the caller that means residues 2-4 on the forward strand says so.
        assert_eq!(
            hgvsp_inframe_indel(
                "ENSP00000001",
                2,
                2,
                "EGE",
                "-",
                Some(&pep),
                Strand::Forward
            ),
            Some("ENSP00000001:p.Glu2_Glu4del".to_string())
        );
    }

    /// Deterministic pseudo-random source. A seeded LCG rather than a `rand`
    /// dependency: the sweep below has to fail the same way twice, or a failure
    /// cannot be investigated.
    struct Lcg(u64);

    impl Lcg {
        fn below(&mut self, n: usize) -> usize {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((self.0 >> 33) as usize) % n
        }
    }

    /// Apply an emitted HGVSp description back to the reference peptide, and
    /// check on the way that every residue it names is the residue the peptide
    /// has at that position.
    ///
    /// This is the inverse of the description, which is what makes it a check
    /// worth having: a span that is well-formed, corroborated and *wrong* still
    /// reconstructs the wrong protein.
    fn apply_to_peptide(description: &str, peptide: &str) -> String {
        let one = |aa3: &str| -> char {
            match aa3 {
                "Ala" => 'A',
                "Glu" => 'E',
                "Gly" => 'G',
                "Lys" => 'K',
                "Met" => 'M',
                other => panic!("unexpected residue {other} in {description}"),
            }
        };
        let body = description
            .split(":p.")
            .nth(1)
            .unwrap_or_else(|| panic!("no :p. in {description}"));
        let (span, replacement) = match body.split_once("delins") {
            Some((span, inserted)) => (
                span,
                inserted
                    .as_bytes()
                    .chunks(3)
                    .map(|c| one(std::str::from_utf8(c).unwrap()))
                    .collect::<String>(),
            ),
            None => (
                body.strip_suffix("del")
                    .unwrap_or_else(|| panic!("neither del nor delins: {description}")),
                String::new(),
            ),
        };

        // `Glu2` or `Glu2_Glu4`, and both ends are checked against the peptide.
        let mut bounds = Vec::new();
        for part in span.split('_') {
            let (aa3, digits) = part.split_at(3);
            let pos: usize = digits
                .parse()
                .unwrap_or_else(|_| panic!("unparsable position in {description}"));
            assert_eq!(
                peptide.chars().nth(pos - 1),
                Some(one(aa3)),
                "{description} names {aa3} at {pos}, peptide has {:?} ({peptide})",
                peptide.chars().nth(pos - 1)
            );
            bounds.push(pos);
        }
        let lo = bounds[0];
        let hi = *bounds.last().unwrap();
        assert!(lo <= hi, "inverted span in {description}");
        format!("{}{}{}", &peptide[..lo - 1], replacement, &peptide[hi..])
    }

    #[test]
    fn every_description_reconstructs_the_protein_the_variant_produces() {
        // The property that matters, over 4,000 shrinking in-frame changes on a
        // four-residue alphabet chosen to make periodic references common: apply
        // the description back to the reference peptide and you must get the
        // protein the variant actually produces. A description can be
        // well-formed, name real residues, and still fail this - that is exactly
        // what #96 was, and what #91 was before it.
        //
        // Both strands, each with the anchor its coordinate convention produces:
        // the start of the span on the forward strand, the end of it on the
        // reverse. Spans start at residue 2 or later, because a deletion reaching
        // the initiation codon is deliberately described as unresolvable (`?`)
        // rather than as a reconstructible event.
        const ALPHABET: &[u8] = b"AEGK";
        let mut rng = Lcg(0x5EED_1234_9ABC_DEF0);
        let mut checked = 0;

        for _ in 0..4_000 {
            let len = 6 + rng.below(14);
            let peptide: String = (0..len)
                .map(|_| ALPHABET[rng.below(ALPHABET.len())] as char)
                .collect();
            let n = 1 + rng.below(4);
            if len < n + 2 {
                continue;
            }
            // 1-based, and never residue 1.
            let lo = 2 + rng.below(len - n);
            let hi = lo + n - 1;
            let reference = &peptide[lo - 1..hi];
            // Shrinking: a pure deletion, or a replacement by fewer residues.
            let m = rng.below(n);
            let replacement: String = (0..m)
                .map(|_| ALPHABET[rng.below(ALPHABET.len())] as char)
                .collect();
            let alt_aas = if m == 0 { "-".to_string() } else { replacement };
            let expected = format!(
                "{}{}{}",
                &peptide[..lo - 1],
                alt_aas.replace('-', ""),
                &peptide[hi..]
            );

            for (anchor, strand) in [(lo as u64, Strand::Forward), (hi as u64, Strand::Reverse)] {
                let got = hgvsp_inframe_indel(
                    "P",
                    anchor,
                    anchor,
                    reference,
                    &alt_aas,
                    Some(peptide.as_bytes()),
                    strand,
                )
                .unwrap_or_else(|| {
                    panic!("no description for {reference}/{alt_aas} at {lo}-{hi} in {peptide}")
                });
                assert_eq!(
                    apply_to_peptide(&got, &peptide),
                    expected,
                    "{got} does not reconstruct {expected} from {peptide} \
                     ({reference}/{alt_aas} at {lo}-{hi}, anchor {anchor} on {strand:?})"
                );
                checked += 1;
            }
        }
        assert!(checked > 5_000, "sweep covered only {checked} cases");
    }

    #[test]
    fn a_two_residue_homopolymer_reads_the_same_from_either_end() {
        // Where the reference is a homopolymer inside a longer run both anchors
        // are corroborated at n = 2, but the 3'-shift converges on one answer, so
        // the strand cannot change it. Worth pinning: it bounds the exposure the
        // ordering fix was needed for to n >= 3.
        let pep: Vec<u8> = "MAKKKA".bytes().collect();
        for (anchor, strand) in [
            (3u64, Strand::Forward),
            (4, Strand::Forward),
            (3, Strand::Reverse),
            (4, Strand::Reverse),
        ] {
            assert_eq!(
                hgvsp_inframe_indel(
                    "ENSP00000001",
                    anchor,
                    anchor,
                    "KK",
                    "-",
                    Some(&pep),
                    strand
                ),
                Some("ENSP00000001:p.Lys4_Lys5del".to_string()),
                "anchor {anchor} on {strand:?}"
            );
        }
    }

    #[test]
    fn test_hgvsp_inframe_indel_matches_vep_terminal_insertion_suppression() {
        for peptide in [b"MKE".as_slice(), b"MKE*".as_slice()] {
            for strand in [Strand::Forward, Strand::Reverse] {
                assert_eq!(
                    hgvsp_inframe_indel("P", 4, 4, "*", "F*", Some(peptide), strand),
                    Some("P:p.Glu3_Ter4insPhe".into())
                );
            }
        }
        // VEP requires both surrounding reference residues for an insertion.
        // At the protein terminus the second one does not exist, so HGVSp is
        // absent rather than an invented delins.
        let pep: Vec<u8> = "MKKRSTV".bytes().collect();
        for strand in [Strand::Forward, Strand::Reverse] {
            assert_eq!(
                hgvsp_inframe_indel("ENSP00000001", 7, 7, "V", "VX", Some(&pep), strand,),
                None,
                "terminal insertion on {strand:?}"
            );
        }
    }

    #[test]
    fn test_hgvsp_inframe_indel_never_emits_substitution_shape() {
        // Every in-frame insertion observed in a clinical panel run, as
        // (protein_start, amino_acids ref, amino_acids alt, prior output).
        // All fifteen previously rendered in a substitution shape: eight as
        // synonymous, seven as a plausible missense. Two (FLT3 at 598 and 596)
        // are ITDs, where a missense reading is clinically misleading.
        let insertions: &[(u64, &str, &str, &str)] = &[
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
        let deletions: &[(u64, &str, &str)] = &[(157, "F", "-"), (43, "YXQ", "-")];

        // Run each case twice: with no peptide, and with a peptide that really
        // carries the stated reference residues at the stated position.
        for &(start, ref_aas, alt_aas, prior) in insertions {
            let pep = peptide_with(start, ref_aas, start as usize + ref_aas.len() + 64);
            for context in [None, Some(pep.as_slice())] {
                // Both strands. An insertion is anchored at `protein_start`
                // regardless (#96), so every assertion below has to hold either
                // way; a strand-dependent answer here would be a bug.
                for strand in [Strand::Forward, Strand::Reverse] {
                    let out = hgvsp_inframe_indel(
                        "ENSP00000001",
                        start,
                        start,
                        ref_aas,
                        alt_aas,
                        context,
                        strand,
                    )
                    .expect("in-frame indel must produce a protein description");
                    let change = out.split(":p.").nth(1).unwrap();
                    // "delins" contains both "del" and "ins", so require one of the
                    // whole forms rather than a substring of another.
                    assert!(
                        change.ends_with("del")
                            || change.ends_with("dup")
                            || change.contains("delins")
                            || change.contains("ins"),
                        "{ref_aas}/{alt_aas} at {start} (was {prior}) rendered \
                     without an indel form: {out}"
                    );
                    assert!(!out.contains('?'), "placeholder residue in {out}");
                    assert!(!out.ends_with('='), "synonymous shape: {out} (was {prior})");
                    assert_ne!(
                        out.split(':').nth(1).unwrap(),
                        prior,
                        "still emitting {prior}"
                    );
                }
            }
        }

        for &(start, ref_aas, alt_aas) in deletions {
            let pep = peptide_with(start, ref_aas, start as usize + ref_aas.len() + 64);
            for context in [None, Some(pep.as_slice())] {
                for strand in [Strand::Forward, Strand::Reverse] {
                    let out = hgvsp_inframe_indel(
                        "ENSP00000001",
                        start,
                        start,
                        ref_aas,
                        alt_aas,
                        context,
                        strand,
                    )
                    .expect("deletion must produce a protein description");
                    assert!(
                        out.ends_with("del"),
                        "deletion regressed on {strand:?}: {out}"
                    );
                }
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

#[cfg(test)]
mod window_tests {
    use super::*;

    /// A change spanning two codons that alters only one of them is a
    /// substitution, not a two-residue delins and not "unchanged".
    ///
    /// `hgvsp()` compares one residue per side, so `EP/ET` read as unchanged and
    /// rendered `p.Glu153=` for a change real VEP 115.1 calls `p.Pro154Thr`.
    /// About 3,000 HGVSp rows per 6,600 ClinVar variants.
    #[test]
    fn a_two_residue_window_that_changes_one_residue_is_a_substitution() {
        let pep = b"MKEPQR".to_vec();
        let call = |start: u64, r: &str, a: &str| {
            hgvsp_inframe_indel("P", start, start, r, a, Some(&pep), Strand::Forward).unwrap()
        };
        // Residues 3 and 4 are E and P; only the second changes.
        assert_eq!(call(3, "EP", "ET"), "P:p.Pro4Thr");
        // Both change: the delins form names the whole run.
        assert_eq!(call(3, "EP", "MG"), "P:p.Glu3_Pro4delinsMetGly");
        // A terminator the change introduces ends the description there.
        assert_eq!(call(3, "EP", "E*"), "P:p.Pro4Ter");
    }

    /// A frameshift's new stop is usually *past* the annotated terminator, and
    /// the Ter distance counts to that stop, not to the end of the reference
    /// protein.
    ///
    /// This is the largest place fastVEP and Ensembl disagree on HGVSp - 373
    /// rows over a 6,600-variant ClinVar sample - and Ensembl is the one that is
    /// wrong. Ten were checked by rebuilding the CDS from the Ensembl 115 GFF3
    /// and FASTA, applying the variant and translating, independently of either
    /// tool; the reference protein lengths came back matching UniProt (TP53 393,
    /// FLCN 579, SLC17A5 495, MPV17 176, NRIP1 1158), and fastVEP matched the
    /// computed distance on all ten, over deltas of 1 to 89 in both directions:
    ///
    /// | first changed residue | computed | fastVEP | VEP 115.1 |
    /// |---|---:|---:|---:|
    /// | MPV17 `p.Leu151Profs` | 39 | 39 | 50 |
    /// | NRIP1 `p.Lys1155Asnfs` | 15 | 15 | 6 |
    /// | TP53 `p.Asp393Thrfs` | 29 | 29 | 89 |
    /// | ARSA `p.Arg498Profs` | 76 | 76 | 21 |
    /// | FLCN `p.Ala541Cysfs` | 61 | 61 | 60 |
    ///
    /// In every one of the ten the new stop lay beyond the reference protein's
    /// end, which is why `cds_and_downstream` has to run past the annotated
    /// terminator: a translation stopping there would have no stop to find.
    #[test]
    fn a_frameshift_stop_past_the_annotated_terminator_is_counted_to() {
        let table = CodonTable::standard();
        //          M   K   L   F   *  | what was the 3' UTR
        let reference = b"ATGAAACTTTTTTAA\
                          CCCGTAACCCTAG"
            .iter()
            .filter(|b| !b.is_ascii_whitespace())
            .copied()
            .collect::<Vec<u8>>();
        // Drop one base from codon 2, shifting the frame from there.
        let mut edited = reference.clone();
        edited.remove(4);

        let out = hgvsp_frameshift("P", &reference, &edited, 1, &table).unwrap();

        // Reference reads M K L F *; edited reads M N F F N P *. The new stop is
        // residue 7 and the first changed residue is 2, so Ter counts 6 - and
        // residue 7 is two past the terminator the reference had at residue 5.
        assert_eq!(out, "P:p.Lys2AsnfsTer6");
    }

    /// A frameshift whose first changed residue is already a terminator is the
    /// nonsense variant it is: `p.Leu1545Ter`, not `p.Leu1545TerfsTer1`.
    ///
    /// 1,611 rows over a 6,600-variant ClinVar sample.
    #[test]
    fn a_frameshift_landing_on_a_terminator_is_written_as_nonsense() {
        let table = CodonTable::standard();
        // ATG AAA CTT TTT TAA: M K L F *
        let cds = b"ATGAAACTTTTTTAA";
        // Deleting the C of codon 3 shifts the frame: ATG AAA TTT TTT AA ->
        // M K F F, no terminator where the reference had one.
        let shifted = hgvsp_frameshift_from_cds(
            "P",
            cds,
            Some(7),
            Some(7),
            &Allele::Sequence(b"C".to_vec()),
            &Allele::Deletion,
            Strand::Forward,
            &table,
        )
        .unwrap();
        assert!(shifted.starts_with("P:p.Leu3Phefs"), "got {shifted}");

        // ATG TAC AAA CTT TAA: deleting the C of codon 2 leaves `TA` in front of
        // the A that follows, so the first changed residue is itself a
        // terminator - which HGVS writes without the `fs`.
        let cds = b"ATGTACAAACTTTAA";
        let nonsense = hgvsp_frameshift_from_cds(
            "P",
            cds,
            Some(6),
            Some(6),
            &Allele::Sequence(b"C".to_vec()),
            &Allele::Deletion,
            Strand::Forward,
            &table,
        )
        .unwrap();
        assert_eq!(nonsense, "P:p.Tyr2Ter");
    }

    /// The edit is "replace the CDS bases the reference covers", which is not
    /// what either per-variant loop used to do: both read `cds_start` as the low
    /// end of the span, complemented an insertion in place instead of
    /// reverse-complementing it, and had no case for a delins at all.
    #[test]
    fn the_edited_cds_replaces_the_reference_bases_on_either_strand() {
        let table = CodonTable::standard();
        // ATG AAA CTT TTT TAA on the transcript's own strand.
        let cds = b"ATGAAACTTTTTTAA";
        // A delins replacing CDS 7-8 (`TT`) with one base: the reference bases
        // have to come out, not just the replacement go in.
        let forward = hgvsp_frameshift_from_cds(
            "P",
            cds,
            Some(7),
            Some(8),
            &Allele::Sequence(b"TT".to_vec()),
            &Allele::Sequence(b"G".to_vec()),
            Strand::Forward,
            &table,
        );
        // The same edit on a reverse-strand transcript arrives with its CDS
        // coordinates the other way round and its alternate reverse-complemented.
        let reverse = hgvsp_frameshift_from_cds(
            "P",
            cds,
            Some(8),
            Some(7),
            &Allele::Sequence(b"AA".to_vec()),
            &Allele::Sequence(b"C".to_vec()),
            Strand::Reverse,
            &table,
        );
        assert_eq!(forward, reverse, "the two strands describe the same edit");
        assert!(forward.is_some());
    }

    #[test]
    fn stop_loss_counts_to_the_next_stop_on_either_strand() {
        let table = CodonTable::standard();
        // M K * Q *; changing the first base of the first TAA to C makes Gln,
        // and the next terminator is two residues later.
        let cds = b"ATGAAATAACAATAA";
        let forward = hgvsp_stop_lost_from_cds(
            "P",
            3,
            "MK*",
            b'Q',
            cds,
            Some(7),
            Some(7),
            &Allele::Sequence(b"T".to_vec()),
            &Allele::Sequence(b"C".to_vec()),
            Strand::Forward,
            &table,
        );
        let reverse = hgvsp_stop_lost_from_cds(
            "P",
            3,
            "MK*",
            b'Q',
            cds,
            Some(7),
            Some(7),
            &Allele::Sequence(b"A".to_vec()),
            &Allele::Sequence(b"G".to_vec()),
            Strand::Reverse,
            &table,
        );

        assert_eq!(forward, Some("P:p.Ter3GlnextTer2".to_string()));
        assert_eq!(reverse, forward);
        assert_eq!(hgvsp_stop_lost_from_cds(
            "P", 3, "MK*", b'-', cds, Some(7), Some(9),
            &Allele::Sequence(b"TAA".to_vec()), &Allele::Deletion,
            Strand::Forward, &table,
        ), Some("P:p.Ter3delextTer1".into()));

        // VEP measures a non-frameshift extension from the full reference
        // peptide length, even when the affected stop is internal.
        let internal_stop = hgvsp_stop_lost_from_cds(
            "P",
            3,
            "MK*Q*",
            b'Q',
            cds,
            Some(7),
            Some(7),
            &Allele::Sequence(b"T".to_vec()),
            &Allele::Sequence(b"C".to_vec()),
            Strand::Forward,
            &table,
        );
        assert_eq!(internal_stop, Some("P:p.Ter3GlnextTer?".to_string()));

        let no_later_stop = hgvsp_stop_lost_from_cds(
            "P",
            3,
            "MK*",
            b'Q',
            b"ATGAAATAACAA",
            Some(7),
            Some(7),
            &Allele::Sequence(b"T".to_vec()),
            &Allele::Sequence(b"C".to_vec()),
            Strand::Forward,
            &table,
        );
        assert_eq!(no_later_stop, Some("P:p.Ter3GlnextTer?".to_string()));
    }
}
