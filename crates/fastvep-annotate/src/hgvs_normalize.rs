//! The intronic HGVSc path: 3'-shifting, and naming a shifted insertion as the
//! duplication it is.
//!
//! [`hgvsc_intronic_shifted`] is the entry point, and both annotation loops call
//! it. They used to carry a copy each - the CLI's shifted, the library's did
//! not - so the same intronic duplication came out normalised from
//! `fastvep annotate` and unnormalised from the server.

use fastvep_cache::providers::SequenceProvider;

/// Convert intronic insertion to dup notation with explicit start/end positions
/// (coding).
///
/// Each end carries its own `(exon anchor, offset)` pair. A span deep inside one
/// intron shares an anchor, but one running past the intron's midpoint is
/// written from the exon on either side - `c.5044+27_5045-47dup` - so the two
/// ends are not interchangeable.
pub fn convert_ins_to_dup_range(
    hgvsc: &str,
    start: (u64, i64),
    end: (u64, i64),
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
        // VEP's `_get_cDNA_position` names an intronic offset anchored on the
        // last coding base from `*`, including the deliberately unusual
        // negative form `*-11` on partial reverse-strand transcripts.
        let terminal_coding_anchor = coding_end == Some(cdna) && off != 0;
        // An offset of 0 is an exonic end, written as the anchor alone. Folding
        // it into the negative arm renders `c.21` as `c.210`.
        let anchor = match coding_end.filter(|&ce| !terminal_coding_anchor && cp >= 0 && cdna > ce)
        {
            Some(ce) => format!("*{}", cdna - ce),
            None if terminal_coding_anchor => "*".to_string(),
            None => cp.to_string(),
        };
        match off.cmp(&0) {
            std::cmp::Ordering::Greater if terminal_coding_anchor => format!("{}{}", anchor, off),
            std::cmp::Ordering::Greater => format!("{}+{}", anchor, off),
            std::cmp::Ordering::Less => format!("{}{}", anchor, off),
            std::cmp::Ordering::Equal => anchor,
        }
    };

    if start == end {
        Some(format!("{}{}dup", prefix, build_pos(start.0, start.1)))
    } else {
        Some(format!(
            "{}{}_{}dup",
            prefix,
            build_pos(start.0, start.1),
            build_pos(end.0, end.1)
        ))
    }
}

/// Convert intronic insertion to dup notation with explicit start/end positions
/// (non-coding). See [`convert_ins_to_dup_range`] for why each end carries its
/// own anchor.
pub fn convert_ins_to_dup_range_noncoding(
    hgvsc: &str,
    start: (u64, i64),
    end: (u64, i64),
) -> Option<String> {
    let prefix_end = hgvsc
        .find(":n.")
        .map(|i| i + 3)
        .or_else(|| hgvsc.find(":c.").map(|i| i + 3))?;
    let prefix = &hgvsc[..prefix_end];

    let build_pos = |cdna: u64, off: i64| -> String {
        match off.cmp(&0) {
            std::cmp::Ordering::Greater => format!("{}+{}", cdna, off),
            std::cmp::Ordering::Less => format!("{}{}", cdna, off),
            std::cmp::Ordering::Equal => format!("{}", cdna),
        }
    };

    if start == end {
        Some(format!("{}{}dup", prefix, build_pos(start.0, start.1)))
    } else {
        Some(format!(
            "{}{}_{}dup",
            prefix,
            build_pos(start.0, start.1),
            build_pos(end.0, end.1)
        ))
    }
}

/// A window of reference the shift walks through, refilled in blocks rather than
/// fetched a base at a time.
///
/// The walk used to call `fetch_sequence` once per base, and twice per base for
/// a deletion, each returning an owned `Vec` - so a variant sliding 4,000 bases
/// down a repeat cost 4,000 heap allocations, or 8,000 as a deletion, once per
/// (variant x transcript x allele).
///
/// The block **grows**, and that is the point: almost every shift stops within a
/// base or two, so a window that opened at its full size would read and copy
/// hundreds of bases to answer two questions - work done before knowing it is
/// needed, which measured 11 % slower over a real indel-only callset than the
/// per-base reads it replaced. Starting small and doubling keeps the common case
/// at one short read and still collapses a long walk to a handful.
struct RefWindow<'a> {
    provider: &'a dyn SequenceProvider,
    chrom: &'a str,
    /// 1-based genomic position of `bases[0]`. Zero until the first fill.
    origin: u64,
    bases: Vec<u8>,
    next_block: u64,
}

impl<'a> RefWindow<'a> {
    /// Enough for a shift that stops immediately, which is nearly all of them.
    const FIRST_BLOCK: u64 = 16;
    /// Past this a walk is in a long repeat and the reads are already amortised.
    const MAX_BLOCK: u64 = 1024;

    fn new(provider: &'a dyn SequenceProvider, chrom: &'a str) -> Self {
        Self {
            provider,
            chrom,
            origin: 0,
            bases: Vec::new(),
            next_block: Self::FIRST_BLOCK,
        }
    }

    /// The base at `pos`, uppercased, or `None` past the contig.
    fn base(&mut self, pos: u64) -> Option<u8> {
        if pos == 0 {
            return None;
        }
        if self.origin == 0 || pos < self.origin || pos >= self.origin + self.bases.len() as u64 {
            let block = self.next_block;
            self.next_block = (block * 2).min(Self::MAX_BLOCK);
            // Centre the block on the request so a walk in either direction has
            // room; the caller's direction is not known here.
            let origin = pos.saturating_sub(block / 2).max(1);
            self.bases = self
                .provider
                .fetch_sequence_slice(self.chrom, origin, origin + block - 1)
                .ok()?;
            self.origin = origin;
            if self.bases.is_empty() {
                return None;
            }
        }
        self.bases
            .get((pos - self.origin) as usize)
            .map(|b| b.to_ascii_uppercase())
    }
}

/// 3' shift an intronic indel along the transcript direction.
///
/// HGVS requires variants to be described at the most 3' position.
/// For intronic deletions and insertions/dups in repetitive regions,
/// the position must be shifted toward the 3' end of the transcript.
///
/// Returns the shifted genomic start and end positions.
// Each argument is an independent coordinate, allele or flag with no
// natural grouping; bundling them into a struct would only move the
// argument list to the call site.
#[allow(clippy::too_many_arguments)]
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

    // VEP's contig-clipped mitochondrial window affects the shift itself,
    // not only the later check for whether HGVS fits inside the transcript.
    if fastvep_genome::is_mitochondrial(chrom) {
        return vep_mitochondrial_genomic_shift(
            seq_provider, chrom, start, end, ref_allele, alt_allele, strand,
        ).unwrap_or((start, end));
    }

    match (ref_allele, alt_allele) {
        // TVA::perform_shift rotates the supplied deletion allele, including
        // ambiguous bases. Reading it back from the genome changes that input.
        (Allele::Sequence(ref_bases), Allele::Deletion) if !ref_bases.is_empty() => {
            let limit = vep_genomic_shift_steps(1000, ref_bases.len(), strand) as u64;
            let (mut s, mut e) = (start, end);
            let mut ahead = RefWindow::new(seq_provider, chrom);
            match strand {
                fastvep_core::Strand::Forward => loop {
                    let next = e + 1;
                    if next > intron_genomic_end || s.abs_diff(start) >= limit {
                        break;
                    }
                    let expected = ref_bases[((s - start) % ref_bases.len() as u64) as usize].to_ascii_uppercase();
                    match ahead.base(next) {
                        Some(base) if base == expected => {
                            s += 1;
                            e += 1;
                        }
                        _ => break,
                    }
                },
                fastvep_core::Strand::Reverse => loop {
                    if s == 0 || s - 1 < intron_genomic_start || s.abs_diff(start) >= limit {
                        break;
                    }
                    let expected = ref_bases[ref_bases.len() - 1 - ((start - s) % ref_bases.len() as u64) as usize].to_ascii_uppercase();
                    match ahead.base(s - 1) {
                        Some(base) if base == expected => {
                            s -= 1;
                            e -= 1;
                        }
                        _ => break,
                    }
                },
            }
            (s, e)
        }
        // An insertion slides over a base that repeats the next base of the
        // inserted sequence, which rotates the sequence by one each step.
        (Allele::Deletion, Allele::Sequence(ins_bases)) if !ins_bases.is_empty() => {
            let ins_len = ins_bases.len();
            let limit = vep_genomic_shift_steps(1000, ins_len, strand);
            let mut pos = start;
            let mut shift = 0usize;
            let mut window = RefWindow::new(seq_provider, chrom);
            match strand {
                fastvep_core::Strand::Forward => loop {
                    if pos > intron_genomic_end || shift >= limit {
                        break;
                    }
                    let expected = ins_bases[shift % ins_len].to_ascii_uppercase();
                    match window.base(pos) {
                        Some(b) if b == expected => {
                            pos += 1;
                            shift += 1;
                        }
                        _ => break,
                    }
                },
                fastvep_core::Strand::Reverse => loop {
                    if pos == 0 || pos - 1 < intron_genomic_start || shift >= limit {
                        break;
                    }
                    let expected = ins_bases[ins_len - 1 - (shift % ins_len)].to_ascii_uppercase();
                    match window.base(pos - 1) {
                        Some(b) if b == expected => {
                            pos -= 1;
                            shift += 1;
                        }
                        _ => break,
                    }
                },
            }
            (pos, pos.saturating_sub(1))
        }
        _ => (start, end),
    }
}

/// Apply VEP's genomic HGVS 3'-shift to an exonic deletion, then map the result
/// back to cDNA. `None` means the shifted deletion no longer has two exonic
/// endpoints and must be rendered by the intronic path or suppressed.
pub(crate) fn exonic_deletion_cdna_span(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &fastvep_genome::Transcript,
    start: u64,
    end: u64,
    ref_allele: &fastvep_core::Allele,
    alt_allele: &fastvep_core::Allele,
) -> Option<(u64, u64)> {
    use fastvep_core::Allele;

    let original = || {
        Some((
            transcript.genomic_to_cdna(start)?,
            transcript.genomic_to_cdna(end)?,
        ))
    };
    let Some(provider) = seq_provider else {
        return original();
    };
    if !matches!(
        (ref_allele, alt_allele),
        (Allele::Sequence(bases), Allele::Deletion) if !bases.is_empty()
    ) {
        return original();
    }

    let (shifted_start, shifted_end) = three_prime_shift_intronic(
        provider,
        chrom,
        start,
        end,
        ref_allele,
        alt_allele,
        transcript.strand,
        transcript.start,
        transcript.end,
    );
    Some((
        transcript.genomic_to_cdna(shifted_start)?,
        transcript.genomic_to_cdna(shifted_end)?,
    ))
}

/// Build a transcript's spliced sequence for one HGVS normalization without
/// retaining it in the transcript cache.
pub(crate) fn transient_spliced_sequence(
    seq_provider: Option<&dyn SequenceProvider>,
    transcript: &fastvep_genome::Transcript,
) -> Option<String> {
    let provider = seq_provider?;
    let mut transcript = transcript.clone();
    transcript
        .build_sequences(|chrom, start, end| {
            provider
                .fetch_sequence(chrom, start, end)
                .map_err(|error| error.to_string())
        })
        .ok()?;
    transcript.spliced_seq
}

/// Whether VEP's genomic 3'-shift makes transcript HGVS unavailable.
///
/// VEP maps the unshifted variation to the transcript slice, adds the genomic
/// shift, and returns no transcript HGVS when the shifted slice end is past the
/// transcript (`TranscriptVariationAllele.pm::hgvs_transcript`). An insertion
/// carries the coordinate before it as its slice end, so a final insertion
/// point one base past the genomic transcript end still fits. Treating both
/// insertion coordinates as ordinary bases incorrectly suppresses a valid
/// terminal duplication.
pub(crate) fn vep_hgvs_shift_exceeds_transcript(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &fastvep_genome::Transcript,
    start: u64,
    end: u64,
    ref_allele: &fastvep_core::Allele,
    alt_allele: &fastvep_core::Allele,
) -> bool {
    use fastvep_core::{Allele, Strand};

    let Some(provider) = seq_provider else {
        return false;
    };
    let is_insertion = matches!(
        (ref_allele, alt_allele),
        (Allele::Deletion, Allele::Sequence(bases)) if !bases.is_empty()
    );
    let is_deletion = matches!(
        (ref_allele, alt_allele),
        (Allele::Sequence(bases), Allele::Deletion) if !bases.is_empty()
    );
    if !is_insertion && !is_deletion {
        return false;
    }

    // One outside base is enough to prove that the shifted representation no
    // longer fits. It also avoids scanning a long repeat past a transcript.
    let (shifted_start, shifted_end) = if fastvep_genome::is_mitochondrial(chrom) {
        vep_mitochondrial_genomic_shift(
            provider,
            chrom,
            start,
            end,
            ref_allele,
            alt_allele,
            transcript.strand,
        )
        .unwrap_or((start, end))
    } else {
        three_prime_shift_intronic(
            provider,
            chrom,
            start,
            end,
            ref_allele,
            alt_allele,
            transcript.strand,
            transcript.start.saturating_sub(1).max(1),
            transcript.end.saturating_add(1),
        )
    };

    match (is_insertion, transcript.strand) {
        // For an insertion, `start` is the base after the insertion and `end`
        // is the base before it. VEP checks the latter transcript-slice
        // coordinate after adding the shift.
        (true, Strand::Forward) => shifted_start > transcript.end.saturating_add(1),
        (true, Strand::Reverse) => shifted_start < transcript.start,
        (false, Strand::Forward) => shifted_end > transcript.end,
        (false, Strand::Reverse) => shifted_start < transcript.start,
    }
}

fn vep_genomic_shift_steps(flank: usize, motif: usize, strand: fastvep_core::Strand) -> usize {
    // TVA perform_shift: inclusive forward loop starts at zero, reverse at
    // one; a negative loop limit resets to the flank length on either strand.
    match strand {
        fastvep_core::Strand::Forward => flank.checked_sub(motif).map(|n| n + 1).unwrap_or(flank),
        fastvep_core::Strand::Reverse => (flank + 1).checked_sub(motif).unwrap_or(flank),
    }.min(flank)
}

/// Reproduce VEP's contig-clipped mitochondrial window. It takes the final
/// 1,000 slice bases as post_seq even near the contig end, where that window
/// begins before the variant. This affects whether terminal HGVS is emitted.
fn vep_mitochondrial_genomic_shift(
    provider: &dyn SequenceProvider,
    chrom: &str,
    start: u64,
    end: u64,
    ref_allele: &fastvep_core::Allele,
    alt_allele: &fastvep_core::Allele,
    strand: fastvep_core::Strand,
) -> Option<(u64, u64)> {
    use fastvep_core::{Allele, Strand};

    let mut motif = match (ref_allele, alt_allele) {
        (Allele::Sequence(bases), Allele::Deletion) if !bases.is_empty() => bases.clone(),
        (Allele::Deletion, Allele::Sequence(bases)) if !bases.is_empty() => bases.clone(),
        _ => return Some((start, end)),
    };
    motif.make_ascii_uppercase();

    const FLANK: u64 = 1_000;
    let slice_start = start.saturating_sub(FLANK).max(1);
    let slice_end = end.saturating_add(FLANK).min(fastvep_genome::MT_LENGTH);
    let sequence = provider
        .fetch_sequence_slice(chrom, slice_start, slice_end)
        .ok()?;
    let flank = FLANK as usize;
    let mut shift = 0u64;

    match strand {
        Strand::Forward => {
            let post = &sequence[sequence.len().saturating_sub(flank)..];
            let limit = vep_genomic_shift_steps(post.len(), motif.len(), strand);
            for &base in post.iter().take(limit) {
                if motif.first().copied()? != base.to_ascii_uppercase() {
                    break;
                }
                motif.rotate_left(1);
                shift += 1;
            }
            Some((start.saturating_add(shift), end.saturating_add(shift)))
        }
        Strand::Reverse => {
            let pre = &sequence[..sequence.len().min(flank)];
            let limit = vep_genomic_shift_steps(pre.len(), motif.len(), strand);
            for &base in pre.iter().rev().take(limit) {
                if motif.last().copied()? != base.to_ascii_uppercase() {
                    break;
                }
                motif.rotate_right(1);
                shift += 1;
            }
            Some((start.saturating_sub(shift), end.saturating_sub(shift)))
        }
    }
}

/// The block a 3'-shifted intronic insertion duplicates, in genomic coordinates.
///
/// After a maximal 3'-shift the duplicated copy can only sit immediately 5' of
/// the insertion point *in transcript orientation*; anything further 3' would
/// have been shifted over. The inserted string rotates one base per position
/// shifted, so the block is compared against that rotated form and not against
/// the string the VCF carried.
///
/// `shifted_start` follows the same convention [`three_prime_shift_intronic`]
/// returns: the insertion sits between genomic `shifted_start - 1` and
/// `shifted_start`.
///
/// The dup anchor used to be derived by re-shifting the *unshifted* insertion
/// point through `three_prime_shift_intronic` over a single position, which
/// walks while the next base repeats the current one - a homopolymer test. A
/// `TG` insertion in a `TGTGTG…` repeat therefore never moved at all and named
/// the copy 16 bases 5' of the one HGVS asks for. That accounted for 1,390 of
/// the 1,538 HGVSc rows disagreeing with real Ensembl VEP 115.1 over a
/// genome-wide HG002 sample, and it is why the two tools split by strand: this
/// walk stays where the VCF put the variant, so fastVEP's anchor was always the
/// genomically-left one whichever way the transcript ran.
pub fn intronic_dup_span(
    seq_provider: &dyn SequenceProvider,
    chrom: &str,
    shifted_start: u64,
    ins_bases: &[u8],
    shift: u64,
    strand: fastvep_core::Strand,
) -> Option<(u64, u64)> {
    let len = ins_bases.len();
    if len == 0 {
        return None;
    }
    // One base of rotation per position travelled, in the direction of travel:
    // shifting 3' on the forward strand walks the first base to the end, and on
    // the reverse strand - where 3' runs towards lower coordinates - the last
    // base walks to the front.
    let rot = (shift % len as u64) as usize;
    let mut rotated = ins_bases.to_vec();
    match strand {
        fastvep_core::Strand::Forward => rotated.rotate_left(rot),
        fastvep_core::Strand::Reverse => rotated.rotate_right(rot),
    }

    // 5' of the insertion point is genomically before it on the forward strand
    // and after it on the reverse.
    let (lo, hi) = match strand {
        fastvep_core::Strand::Forward => {
            let hi = shifted_start.checked_sub(1)?;
            (hi.checked_sub(len as u64 - 1)?, hi)
        }
        fastvep_core::Strand::Reverse => (shifted_start, shifted_start + len as u64 - 1),
    };
    if lo == 0 {
        return None;
    }

    let block = seq_provider.fetch_sequence_slice(chrom, lo, hi).ok()?;
    if block.len() != len {
        return None;
    }
    block
        .iter()
        .zip(rotated.iter())
        .all(|(a, b)| a.eq_ignore_ascii_case(b))
        .then_some((lo, hi))
}

/// Rewrite a 3'-shifted intronic insertion as a duplication, when it is one.
///
/// `hgvsc` is the insertion notation already built for the shifted position, and
/// is returned rewritten. `None` means the insertion does not duplicate the
/// adjacent block.
///
/// `coding_start` is `None` for a transcript numbered from its first base, which
/// selects `n.` numbering.
// The arguments are independent coordinates, alleles and transcript state with
// no natural grouping; a struct would only move the list to the call site.
#[allow(clippy::too_many_arguments)]
pub fn intronic_ins_as_dup(
    seq_provider: &dyn SequenceProvider,
    chrom: &str,
    transcript: &fastvep_genome::Transcript,
    hgvsc: &str,
    shifted_start: u64,
    ins_bases: &[u8],
    shift: u64,
    coding_start: Option<u64>,
    coding_end: Option<u64>,
) -> Option<String> {
    let (lo, hi) = intronic_dup_span(
        seq_provider,
        chrom,
        shifted_start,
        ins_bases,
        shift,
        transcript.strand,
    )?;
    // Transcript order, not genomic: on the reverse strand the block's 5' end is
    // its higher coordinate.
    let (first, last) = match transcript.strand {
        fastvep_core::Strand::Forward => (lo, hi),
        fastvep_core::Strand::Reverse => (hi, lo),
    };
    // VEP's `_genomic_shift` is not bounded by the intron containing the input.
    // It maps each shifted endpoint afterwards, so a duplicated block may be
    // intronic, exonic, or cross a splice boundary.
    let start = crate::intronic_or_exonic_cdna(transcript, first)?;
    let end = crate::intronic_or_exonic_cdna(transcript, last)?;
    match coding_start {
        Some(cs) => convert_ins_to_dup_range(hgvsc, start, end, cs, coding_end),
        None => convert_ins_to_dup_range_noncoding(hgvsc, start, end),
    }
}

/// Build the HGVSc for a variant reaching into an intron: 3'-shifted, and
/// written as a duplication where the shifted insertion sits against the block
/// it copies.
///
/// Both annotation loops call this. They had drifted: the CLI's copy shifted and
/// converted to `dup`, the library's did neither, so the same intronic
/// duplication came out normalised from `fastvep annotate` and unnormalised from
/// the server and the web UI.
///
/// `genomic_ref` and `genomic_alt` are the alleles as the VCF carried them;
/// `hgvs_ref` and `hgvs_alt` are the same pair in transcript orientation. The
/// shift reads the reference, so it needs the genomic pair; the notation is
/// written from the transcript pair.
///
/// `coding_start` is `None` for a transcript numbered from its first base, which
/// selects `n.` numbering.
// The arguments are independent coordinates, alleles and transcript state with
// no natural grouping; a struct would only move the list to the call site.
#[allow(clippy::too_many_arguments)]
pub fn hgvsc_intronic_shifted(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &fastvep_genome::Transcript,
    versioned_tid: &str,
    var_start: u64,
    var_end: u64,
    genomic_ref: &fastvep_core::Allele,
    genomic_alt: &fastvep_core::Allele,
    hgvs_ref: &fastvep_core::Allele,
    hgvs_alt: &fastvep_core::Allele,
    coding_start: Option<u64>,
    coding_end: Option<u64>,
) -> Option<String> {
    use fastvep_core::{Allele, Strand};

    let is_insertion = matches!(
        (hgvs_ref, hgvs_alt),
        (Allele::Deletion, Allele::Sequence(_))
    );
    let is_deletion = matches!(
        (hgvs_ref, hgvs_alt),
        (Allele::Sequence(_), Allele::Deletion)
    );
    let is_indel = is_insertion || is_deletion;

    // Decline malformed cache models whose declared transcript span extends
    // beyond every exon and intron. Real VEP transcript slices do not have such
    // holes, and shifting an unmapped variant into a later exon invents HGVS.
    let maps_to_transcript = |position| {
        transcript.genomic_to_cdna(position).is_some()
            || transcript.intron_bounds_at(position).is_some()
    };
    if !maps_to_transcript(var_start) && !maps_to_transcript(var_end) {
        return None;
    }

    // VEP 115's `_genomic_shift` walks on the genomic reference before mapping
    // the result back to the transcript. It can therefore cross splice
    // boundaries; limiting this to the input intron changes both the final
    // coordinate and whether an insertion is recognized as a duplication.
    let bounds = is_indel.then_some((transcript.start, transcript.end));
    let (shifted_start, shifted_end) = match seq_provider.filter(|_| is_indel).zip(bounds) {
        Some((sp, (intron_start, intron_end))) => three_prime_shift_intronic(
            sp,
            chrom,
            var_start,
            var_end,
            genomic_ref,
            genomic_alt,
            transcript.strand,
            intron_start,
            intron_end,
        ),
        None => (var_start, var_end),
    };
    // Distance travelled, which is also how far the inserted string rotated.
    let shift = match transcript.strand {
        Strand::Forward => shifted_start.saturating_sub(var_start),
        Strand::Reverse => var_start.saturating_sub(shifted_start),
    };
    // `hgvs_alt` is already in transcript orientation, where the shift always
    // travels 3' whichever way the transcript runs, so the string always rotates
    // left - the strand is spent before this point, on reverse-complementing it.
    //
    // Rotating right on the reverse strand instead, which is what this did,
    // named the right position with the wrong bases: `c.21-21_21-20insGTC` for a
    // variant the same transcript calls `insCGT` when the VCF spells it one base
    // over. Genomically that rotation is correct and [`intronic_dup_span`] keeps
    // it, because that function reads the reference rather than the transcript.
    let shifted_alt = match hgvs_alt {
        Allele::Sequence(ins) if is_insertion && shift > 0 && !ins.is_empty() => {
            let mut rotated = ins.clone();
            let k = (shift % rotated.len() as u64) as usize;
            rotated.rotate_left(k);
            Allele::Sequence(rotated)
        }
        other => other.clone(),
    };

    // An insertion is written over the two bases it sits between, and both are
    // mapped: `shifted_end` and `shifted_start` are the pair, in transcript
    // order on the forward strand and reversed on the reverse.
    //
    // Letting the renderer infer the second coordinate as `offset + 1` instead
    // breaks across the middle of an intron, where `+n` counts from one exon and
    // `-m` from the next: an insertion between `c.20+30` and `c.21-30` came out
    // `c.20+30_20+31ins…`, naming a base past the half the donor-side offsets
    // reach.
    let (anchor_pos, second_pos) = if is_insertion {
        match transcript.strand {
            Strand::Forward => (shifted_end, Some(shifted_end + 1)),
            Strand::Reverse => (shifted_end + 1, Some(shifted_end)),
        }
    } else {
        (
            shifted_start,
            (shifted_start != shifted_end).then_some(shifted_end),
        )
    };
    let (cdna_pos, offset) = crate::intronic_or_exonic_cdna(transcript, anchor_pos)?;
    let (end_cdna, end_offset) = second_pos
        .and_then(|p| crate::intronic_or_exonic_cdna(transcript, p))
        .map(|(c, o)| (Some(c), Some(o)))
        .unwrap_or((None, None));

    let hgvsc = match coding_start {
        Some(cs) => fastvep_hgvs::hgvsc_intronic_range(
            versioned_tid,
            cdna_pos,
            offset,
            end_cdna,
            end_offset,
            hgvs_ref,
            &shifted_alt,
            cs,
            coding_end,
        ),
        None => fastvep_hgvs::hgvsc_noncoding_intronic_range(
            versioned_tid,
            cdna_pos,
            offset,
            end_cdna,
            end_offset,
            hgvs_ref,
            &shifted_alt,
        ),
    }?;

    if let (true, Allele::Sequence(ins), Some(sp)) = (is_insertion, genomic_alt, seq_provider) {
        if hgvsc.contains("ins") && !ins.is_empty() {
            if let Some(dup) = intronic_ins_as_dup(
                sp,
                chrom,
                transcript,
                &hgvsc,
                shifted_start,
                ins,
                shift,
                coding_start,
                coding_end,
            ) {
                return Some(dup);
            }
        }
    }
    // VEP `_get_cDNA_position` cannot map a base past the final genomic exon.
    // A terminal duplication can still be valid because its two coordinates
    // describe the preceding duplicated block, handled above.
    if anchor_pos > transcript.end || second_pos.is_some_and(|pos| pos > transcript.end) {
        return None;
    }
    Some(hgvsc)
}

#[cfg(test)]
mod tests {
    #[test]
    fn deletion_shift_rotates_input_reference_including_ambiguity() {
        use fastvep_core::{Allele, Strand};
        let genome = StrRef("AAAAAAAAAA");
        for strand in [Strand::Forward, Strand::Reverse] {
            assert_eq!(three_prime_shift_intronic(&genome, "1", 5, 5,
                &Allele::from_str("N"), &Allele::Deletion, strand, 1, 10), (5, 5));
        }
    }
    use super::*;
    use anyhow::{anyhow, Result};
    use fastvep_core::Strand;
    use fastvep_genome::{Exon, Gene, Transcript};

    /// Minimal `SequenceProvider` over a 1-based reference string for one contig,
    /// mirroring the real readers' contract: 1-based inclusive, `Err` past the end.
    struct StrRef(&'static str);
    impl SequenceProvider for StrRef {
        fn fetch_sequence(&self, _chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
            if start < 1 || end < start {
                return Err(anyhow!("bad range"));
            }
            let b = self.0.as_bytes();
            let s0 = (start - 1) as usize;
            if s0 >= b.len() {
                return Err(anyhow!("past contig end"));
            }
            Ok(b[s0..(end as usize).min(b.len())].to_vec())
        }
    }

    struct HomopolymerRef;
    impl SequenceProvider for HomopolymerRef {
        fn fetch_sequence(&self, _chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
            if start < 1 || end < start || end > 1_000 {
                return Err(anyhow!("bad range"));
            }
            Ok(vec![b'C'; (end - start + 1) as usize])
        }
    }

    struct TerminalInsertionRef;
    impl SequenceProvider for TerminalInsertionRef {
        fn fetch_sequence(&self, _chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
            if start < 1 || end < start || end > 1_000 {
                return Err(anyhow!("bad range"));
            }
            Ok((start..=end)
                .map(|position| if position == 100 { b'C' } else { b'G' })
                .collect())
        }
    }

    struct MitoBoundaryRef;
    impl SequenceProvider for MitoBoundaryRef {
        fn fetch_sequence(&self, _chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
            if start < 1 || end < start || end > fastvep_genome::MT_LENGTH {
                return Err(anyhow!("bad range"));
            }
            Ok((start..=end)
                .map(|position| if position == 10_059 || position == 15_954 { b'A' } else { b'G' })
                .collect())
        }
    }

    struct SpliceBoundaryRef;
    impl SequenceProvider for SpliceBoundaryRef {
        fn fetch_sequence(&self, _chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
            if start < 1 || end < start || end > 1_000 {
                return Err(anyhow!("bad range"));
            }
            Ok((start..=end)
                .map(|pos| {
                    if matches!(pos, 20 | 21 | 80 | 81) {
                        b'T'
                    } else {
                        b'G'
                    }
                })
                .collect())
        }
    }

    struct SpliceJumpRef;
    impl SequenceProvider for SpliceJumpRef {
        fn fetch_sequence(&self, _chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
            if start < 1 || end < start || end > 1_000 {
                return Err(anyhow!("bad range"));
            }
            Ok((start..=end)
                .map(|pos| if matches!(pos, 20 | 81) { b'T' } else { b'G' })
                .collect())
        }
    }

    struct ExonicInsertionRef;
    impl SequenceProvider for ExonicInsertionRef {
        fn fetch_sequence(&self, _chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
            if start < 1 || end < start || end > 1_000 {
                return Err(anyhow!("bad range"));
            }
            Ok((start..=end)
                .map(|position| match position {
                    20 | 22 => b'T',
                    21 => b'G',
                    _ => b'A',
                })
                .collect())
        }
    }

    /// Two exons on `strand` with one intron between them, so an intronic
    /// position has an anchor on either side. Exon 1 is 1..=20, exon 2 is
    /// 81..=100, and the intron is 21..=80.
    fn transcript(strand: Strand) -> Transcript {
        let exon = |start: u64, end: u64, rank: u32| Exon {
            stable_id: format!("ENSE{}", rank),
            start,
            end,
            strand,
            phase: 0,
            end_phase: 0,
            rank,
        };
        Transcript {
            stable_id: "ENST00000000001".into(),
            version: Some(1),
            gene: Gene {
                stable_id: "ENSG00000000001".into(),
                symbol: Some("TEST".into()),
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "1".into(),
                start: 1,
                end: 100,
                strand,
            },
            biotype: "protein_coding".into(),
            chromosome: "1".into(),
            start: 1,
            end: 100,
            strand,
            exons: vec![exon(1, 20, 1), exon(81, 100, 2)],
            translation: None,
            cdna_coding_start: Some(1),
            cdna_coding_end: Some(40),
            coding_region_start: None,
            coding_region_end: None,
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
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
        }
    }

    #[test]
    fn transient_spliced_sequence_is_not_cached_on_the_transcript() {
        let tr = transcript(Strand::Forward);
        let sequence = transient_spliced_sequence(Some(&HomopolymerRef), &tr).unwrap();

        assert_eq!(sequence, "C".repeat(40));
        assert!(tr.spliced_seq.is_none());
    }

    #[test]
    fn duplication_at_a_terminal_coding_anchor_uses_vep_star_notation() {
        assert_eq!(
            convert_ins_to_dup_range(
                "ENST00000000001.1:c.318-11dup",
                (400, -11),
                (400, -11),
                83,
                Some(400),
            )
            .as_deref(),
            Some("ENST00000000001.1:c.*-11dup")
        );
    }

    #[test]
    fn vep_boundary_check_allows_a_terminal_insertion_but_not_a_shift_past_it() {
        let tr = transcript(Strand::Forward);
        let inserted = fastvep_core::Allele::from_str("C");

        assert!(!vep_hgvs_shift_exceeds_transcript(
            Some(&TerminalInsertionRef),
            "1",
            &tr,
            100,
            99,
            &fastvep_core::Allele::Deletion,
            &inserted,
        ));
        assert!(vep_hgvs_shift_exceeds_transcript(
            Some(&HomopolymerRef),
            "1",
            &tr,
            100,
            99,
            &fastvep_core::Allele::Deletion,
            &inserted,
        ));
        for (inserted, expected) in [("C", Some("T:c.40dup")), ("CA", None)] {
            let inserted = fastvep_core::Allele::from_str(inserted);
            assert_eq!(
                hgvsc_intronic_shifted(
                    Some(&TerminalInsertionRef),
                    "1",
                    &tr,
                    "T",
                    100,
                    99,
                    &fastvep_core::Allele::Deletion,
                    &inserted,
                    &fastvep_core::Allele::Deletion,
                    &inserted,
                    Some(1),
                    Some(40),
                )
                .as_deref(),
                expected,
            );
        }
    }

    #[test]
    fn mitochondrial_shift_handles_alleles_longer_than_the_search_flank() {
        use fastvep_core::Allele;
        for (length, forward, reverse) in [(999, 2, 2), (1000, 1, 1), (1001, 1000, 0), (1002, 1000, 1000), (1039, 1000, 1000)] {
            let sequence = Allele::Sequence(vec![b'G'; length]);
            for (reference, alternate, end) in [(&sequence, &Allele::Deletion, 4000 + length as u64 - 1), (&Allele::Deletion, &sequence, 3999)] {
                for (strand, shift) in [(Strand::Forward, forward), (Strand::Reverse, -reverse)] {
                    assert_eq!(vep_mitochondrial_genomic_shift(
                        &MitoBoundaryRef, "MT", 4000, end, reference, alternate, strand,
                    ), Some((4000u64.checked_add_signed(shift).unwrap(), end.checked_add_signed(shift).unwrap())),
                        "length={length}, strand={strand:?}");
                    assert_eq!(three_prime_shift_intronic(
                        &MitoBoundaryRef, "1", 4000, end, reference, alternate, strand, 1, 16000,
                    ), (4000u64.checked_add_signed(shift).unwrap(), end.checked_add_signed(shift).unwrap()),
                        "nuclear length={length}, strand={strand:?}");
                }
            }
        }
    }

    #[test]
    fn vep_boundary_check_matches_linear_and_clipped_mitochondrial_windows() {
        let tr = transcript(Strand::Forward);
        let reference = fastvep_core::Allele::from_str("C");

        assert!(vep_hgvs_shift_exceeds_transcript(
            Some(&HomopolymerRef),
            "1",
            &tr,
            100,
            100,
            &reference,
            &fastvep_core::Allele::Deletion,
        ));

        let mut early_mt = transcript(Strand::Forward);
        early_mt.start = 9_991;
        early_mt.end = 10_058;
        early_mt.gene.start = early_mt.start;
        early_mt.gene.end = early_mt.end;
        early_mt.exons[0].start = early_mt.start;
        early_mt.exons[0].end = early_mt.end;
        let mt_reference = fastvep_core::Allele::from_str("A");
        assert!(vep_hgvs_shift_exceeds_transcript(
            Some(&MitoBoundaryRef),
            "MT",
            &early_mt,
            early_mt.end,
            early_mt.end,
            &mt_reference,
            &fastvep_core::Allele::Deletion,
        ));

        let mut late_mt = early_mt.clone();
        late_mt.start = 15_888;
        late_mt.end = 15_953;
        late_mt.gene.start = late_mt.start;
        late_mt.gene.end = late_mt.end;
        late_mt.exons[0].start = late_mt.start;
        late_mt.exons[0].end = late_mt.end;
        // Immediate genomic sequence permits an A insertion shift here, but
        // VEP's contig-clipped post window begins earlier with G.
        assert_eq!(three_prime_shift_intronic(
            &MitoBoundaryRef, "MT", 15_954, 15_953,
            &fastvep_core::Allele::Deletion, &mt_reference, Strand::Forward,
            1, fastvep_genome::MT_LENGTH,
        ), (15_954, 15_953));
        assert!(!vep_hgvs_shift_exceeds_transcript(
            Some(&MitoBoundaryRef),
            "MT",
            &late_mt,
            late_mt.end,
            late_mt.end,
            &mt_reference,
            &fastvep_core::Allele::Deletion,
        ));
    }

    #[test]
    fn exonic_deletion_can_shift_to_the_first_intronic_base() {
        for (strand, position) in [(Strand::Forward, 20), (Strand::Reverse, 81)] {
            let tr = transcript(strand);
            let reference = fastvep_core::Allele::from_str("T");
            assert_eq!(
                exonic_deletion_cdna_span(
                    Some(&SpliceBoundaryRef),
                    "1",
                    &tr,
                    position,
                    position,
                    &reference,
                    &fastvep_core::Allele::Deletion,
                ),
                None
            );

            let out = hgvsc_intronic_shifted(
                Some(&SpliceBoundaryRef),
                "1",
                &tr,
                "ENST00000000001.1",
                position,
                position,
                &reference,
                &fastvep_core::Allele::Deletion,
                &reference,
                &fastvep_core::Allele::Deletion,
                Some(1),
                Some(40),
            );
            assert_eq!(
                out.as_deref(),
                Some("ENST00000000001.1:c.20+1del"),
                "{strand:?}"
            );
        }
    }

    #[test]
    fn exonic_duplication_retains_rotation_when_the_shift_crosses_an_intron() {
        use fastvep_core::{Allele, GenomicPosition};
        // OAZ3 source-traced window: the genomic shift crosses a one-base
        // intron, but both final flanks and the duplicated block are exonic.
        let reference = StrRef("GCCTCCAGTGCTCCTGAGTCCCTAGTAGGCCTCCAGGAGGGCAAAAGCAC");
        let length = reference.0.len();
        let gff = format!("1\ttest\tgene\t1\t{length}\t.\t+\t.\tID=gene:G;biotype=protein_coding\n\
1\ttest\tmRNA\t1\t{length}\t.\t+\t.\tID=transcript:T;Parent=gene:G;biotype=protein_coding\n\
1\ttest\texon\t1\t14\t.\t+\t.\tParent=transcript:T;rank=1\n\
1\ttest\texon\t16\t{length}\t.\t+\t.\tParent=transcript:T;rank=2\n\
1\ttest\tCDS\t1\t14\t.\t+\t0\tParent=transcript:T;protein_id=P\n\
1\ttest\tCDS\t16\t{length}\t.\t+\t1\tParent=transcript:T;protein_id=P\n");
        let mut transcripts = fastvep_cache::gff::parse_gff3(gff.as_bytes()).unwrap();
        let tr = &mut transcripts[0];
        tr.build_sequences(|chrom, start, end| reference.fetch_sequence(chrom, start, end)
            .map_err(|error| error.to_string())).unwrap();
        let result = fastvep_consequence::ConsequencePredictor::default().predict(
            &GenomicPosition::new("1", 13, 12, Strand::Forward),
            &Allele::Deletion, &[Allele::from_str("CCTGAGTC")], &[tr], None,
        );
        let allele = &result.transcript_consequences[0].allele_consequences[0];
        assert_eq!(crate::hgvsc_for_allele(Some(&reference), "1", tr, "T", allele),
            Some("T:c.15_22dup".into()));
    }

    #[test]
    fn a_genomic_deletion_shift_does_not_jump_an_intron() {
        let tr = transcript(Strand::Forward);
        let span = exonic_deletion_cdna_span(
            Some(&SpliceJumpRef),
            "1",
            &tr,
            20,
            20,
            &fastvep_core::Allele::from_str("T"),
            &fastvep_core::Allele::Deletion,
        );

        // cDNA bases 20 and 21 are both T, but genomic base 21 is G. VEP
        // therefore leaves the deletion at cDNA 20 instead of skipping the
        // intron and shifting it to the next exon.
        assert_eq!(span, Some((20, 20)));
    }

    #[test]
    fn exonic_insertion_can_shift_to_an_intronic_duplication() {
        let tr = transcript(Strand::Forward);
        let out = hgvsc_intronic_shifted(
            Some(&ExonicInsertionRef),
            "1",
            &tr,
            "ENST00000000001.1",
            20,
            19,
            &fastvep_core::Allele::Deletion,
            &fastvep_core::Allele::from_str("TG"),
            &fastvep_core::Allele::Deletion,
            &fastvep_core::Allele::from_str("TG"),
            Some(1),
            Some(40),
        );

        assert_eq!(out.as_deref(), Some("ENST00000000001.1:c.20+1_20+2dup"));
    }

    #[test]
    fn an_unmapped_deletion_is_not_promoted_into_a_later_exon() {
        let mut tr = transcript(Strand::Forward);
        tr.exons.remove(0);
        tr.cdna_coding_start = None;
        tr.cdna_coding_end = None;

        let out = hgvsc_intronic_shifted(
            Some(&HomopolymerRef),
            "1",
            &tr,
            "ENST00000000001.1",
            80,
            80,
            &fastvep_core::Allele::from_str("C"),
            &fastvep_core::Allele::Deletion,
            &fastvep_core::Allele::from_str("C"),
            &fastvep_core::Allele::Deletion,
            None,
            None,
        );

        assert_eq!(out, None);
    }

    #[test]
    fn reverse_exonic_deletion_shifts_without_a_cached_transcript_sequence() {
        let tr = transcript(Strand::Reverse);
        let span = exonic_deletion_cdna_span(
            Some(&HomopolymerRef),
            "1",
            &tr,
            91,
            91,
            &fastvep_core::Allele::from_str("C"),
            &fastvep_core::Allele::Deletion,
        );
        assert_eq!(span, Some((40, 40)));
    }

    #[test]
    fn deletion_can_shift_from_an_intron_into_an_exon() {
        let tr = transcript(Strand::Forward);
        let out = hgvsc_intronic_shifted(
            Some(&HomopolymerRef),
            "1",
            &tr,
            "ENST00000000001.1",
            80,
            81,
            &fastvep_core::Allele::from_str("CC"),
            &fastvep_core::Allele::Deletion,
            &fastvep_core::Allele::from_str("CC"),
            &fastvep_core::Allele::Deletion,
            Some(1),
            Some(40),
        );

        assert_eq!(out.as_deref(), Some("ENST00000000001.1:c.39_40del"));
    }

    /// 20 exonic bases, then a `TG` repeat filling the intron, then exon 2. An
    /// insertion of `TG` anywhere in that repeat is the same variant.
    const TG_REPEAT: &str = "AAAAAAAAAAAAAAAAAAAA\
                             TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG\
                             CCCCCCCCCCCCCCCCCCCC";

    /// The duplicated block sits immediately 5' of the *shifted* insertion
    /// point, so a `TG` insertion that travelled the whole repeat names the last
    /// two bases of it, not the two it was written against.
    #[test]
    fn dup_span_follows_the_shifted_insertion_not_the_vcf_position() {
        let r = StrRef(TG_REPEAT);
        // Insertion point right after exon 1; the repeat runs 21..=80, so the
        // maximal 3' shift on the forward strand travels its full 60 bases.
        let span = intronic_dup_span(&r, "1", 81, b"TG", 60, Strand::Forward);
        assert_eq!(span, Some((79, 80)));
    }

    /// The inserted string rotates one base per position travelled. After an odd
    /// shift of a two-base insert the block is `GT`, and reading the unrotated
    /// `TG` against the reference would reject a duplication that is real.
    #[test]
    fn dup_span_compares_against_the_rotated_insert() {
        let r = StrRef(TG_REPEAT);
        // An odd shift lands the insertion point between a G and a T.
        let span = intronic_dup_span(&r, "1", 80, b"TG", 59, Strand::Forward);
        assert_eq!(span, Some((78, 79)));
        let block = r.fetch_sequence("1", 78, 79).unwrap();
        assert_eq!(&block, b"GT", "the block is the rotated form, not `TG`");
    }

    /// 3' on the reverse strand runs towards lower coordinates, so the block a
    /// reverse-strand insertion duplicates lies *after* the insertion point.
    #[test]
    fn dup_span_reads_the_other_side_on_the_reverse_strand() {
        let r = StrRef(TG_REPEAT);
        let span = intronic_dup_span(&r, "1", 21, b"TG", 0, Strand::Reverse);
        assert_eq!(span, Some((21, 22)));
    }

    /// A one-base walk through a homopolymer is not the repeat test: this is the
    /// shape the old dup anchor used, and it stopped after one base of `TGTG…`.
    #[test]
    fn dup_span_is_none_when_the_block_does_not_repeat() {
        let r = StrRef(TG_REPEAT);
        // Inside the exon's poly-A run, a `TG` insert duplicates nothing.
        assert_eq!(
            intronic_dup_span(&r, "1", 15, b"TG", 0, Strand::Forward),
            None
        );
    }

    #[test]
    fn dup_span_declines_rather_than_reading_past_the_contig() {
        let r = StrRef(TG_REPEAT);
        assert_eq!(
            intronic_dup_span(&r, "1", 1, b"TG", 0, Strand::Forward),
            None
        );
        assert_eq!(
            intronic_dup_span(&r, "1", 0, b"TG", 0, Strand::Forward),
            None
        );
        assert_eq!(
            intronic_dup_span(&r, "1", 81, b"", 0, Strand::Forward),
            None
        );
    }

    /// The whole rewrite: a shifted intronic insertion becomes a `dup` naming the
    /// block at its 3' end, written from the anchor its own side of the intron.
    #[test]
    fn intronic_insertion_becomes_a_dup_at_its_shifted_position() {
        let r = StrRef(TG_REPEAT);
        let tr = transcript(Strand::Forward);
        let out = intronic_ins_as_dup(
            &r,
            "1",
            &tr,
            "ENST00000000001.1:c.20+59_20+60insTG",
            81,
            b"TG",
            60,
            Some(1),
            Some(40),
        );
        // 79 and 80 are the last two intronic bases, one and two before exon 2,
        // so they are written from the *downstream* anchor: c.21-2_21-1.
        assert_eq!(out.as_deref(), Some("ENST00000000001.1:c.21-2_21-1dup"));
    }

    /// VEP shifts across splice boundaries and maps the resulting duplication
    /// endpoints independently.
    #[test]
    fn a_block_leaving_the_intron_is_mapped_as_a_duplication() {
        // Exon 2 begins with `C`s, so a `CC` insert at the intron's 3' edge
        // duplicates a block that starts inside the exon.
        let r = StrRef(TG_REPEAT);
        let tr = transcript(Strand::Reverse);
        let out = intronic_ins_as_dup(
            &r,
            "1",
            &tr,
            "ENST00000000001.1:c.21-1_21insCC",
            81,
            b"CC",
            0,
            Some(1),
            Some(40),
        );
        assert_eq!(out.as_deref(), Some("ENST00000000001.1:c.19_20dup"));
    }

    /// The property the whole shift exists to provide, and the one that broke:
    /// every spelling of the same insertion inside a repeat must name the same
    /// block. Walking every position of the `TG` repeat, each with the rotation
    /// that spelling carries, must land on one answer.
    ///
    /// Against real Ensembl VEP 115.1 this was 0 of 168 transcript rows agreeing
    /// across the two spellings of six genome-wide variants; VEP agreed on all
    /// 168.
    #[test]
    fn every_spelling_of_one_insertion_names_the_same_block() {
        let r = StrRef(TG_REPEAT);
        // The repeat runs 21..=80. An insertion written at offset `k` into it
        // carries the insert rotated by `k`, and has `60 - k` left to travel.
        let answers: std::collections::HashSet<_> = (0..60)
            .map(|k| {
                let mut ins = b"TG".to_vec();
                ins.rotate_left(k % 2);
                intronic_dup_span(&r, "1", 81, &ins, (60 - k) as u64, Strand::Forward)
            })
            .collect();
        assert_eq!(
            answers,
            [Some((79, 80))].into_iter().collect(),
            "every spelling must name the last two bases of the repeat"
        );
    }

    /// A span running past the intron's own midpoint is still one intron, and is
    /// written from the exon on either side. Requiring one shared anchor instead
    /// left 160 such rows as `ins` where Ensembl writes `c.5044+27_5045-47dup`.
    #[test]
    fn a_span_crossing_the_intron_midpoint_is_still_one_dup() {
        // A 60-base insert filling the whole intron duplicates all of it.
        let r = StrRef(TG_REPEAT);
        let tr = transcript(Strand::Forward);
        let ins: Vec<u8> = TG_REPEAT.as_bytes()[20..80].to_vec();
        let out = intronic_ins_as_dup(
            &r,
            "1",
            &tr,
            "ENST00000000001.1:c.20+60_21-0ins…",
            81,
            &ins,
            0,
            Some(1),
            Some(40),
        );
        assert_eq!(out.as_deref(), Some("ENST00000000001.1:c.20+1_21-1dup"));
    }
}
