use fastvep_core::{Allele, Consequence, GenomicPosition, Impact, Strand};
use fastvep_genome::{is_mitochondrial, mitochondrial_codon_table, CodonTable, Transcript};
use std::collections::HashMap;
use std::sync::Arc;

use crate::splice;

/// Result of consequence prediction for a variant against a transcript.
#[derive(Debug, Clone)]
pub struct TranscriptConsequence {
    pub transcript_id: Arc<str>,
    pub gene_id: Arc<str>,
    pub gene_symbol: Option<Arc<str>>,
    pub biotype: Arc<str>,
    pub allele_consequences: Vec<AlleleConsequenceResult>,
    pub canonical: bool,
    pub strand: Strand,
}

/// Consequence result for a single allele against a transcript.
#[derive(Debug, Clone)]
pub struct AlleleConsequenceResult {
    /// ALT as represented by the parsed VCF record. Keep this for output and
    /// supplementary-source matching.
    pub allele: Allele,
    /// Allele-local minimal representation retained for lookup and output consumers.
    /// Consequence prediction uses the original anchored VCF allele.
    pub normalized_position: GenomicPosition,
    pub normalized_ref_allele: Allele,
    pub normalized_alt_allele: Allele,
    pub consequences: Vec<Consequence>,
    /// Coding frameshift predicate before genomic allele-shape output gates.
    /// VEP calls this predicate directly when selecting the HGVSp path.
    pub frameshift: bool,
    pub impact: Impact,
    /// Transcript-relative coordinate pairs, in *strand order*: `start` is the
    /// position of the genomically leftmost affected base and `end` the
    /// rightmost. On the reverse strand a transcript runs right-to-left, so
    /// `start > end` there and the pair is not `(lo, hi)`.
    ///
    /// Use [`AlleleConsequenceResult::protein_range`] (and sort the others the
    /// same way) whenever a span is wanted. A bare `protein_start` is only the
    /// first affected residue on the forward strand; the exported
    /// `Protein_position` matches Ensembl VEP because the output layer sorts the
    /// pair before printing it.
    pub cdna_start: Option<u64>,
    pub cdna_end: Option<u64>,
    pub cds_start: Option<u64>,
    pub cds_end: Option<u64>,
    pub protein_start: Option<u64>,
    pub protein_end: Option<u64>,
    /// `(replaced residues, replacement residues)`.
    ///
    /// The replaced residues are **not** guaranteed to begin at
    /// `protein_start`. They are built from the lower of the two CDS
    /// coordinates for deletions and from `cds_start` otherwise, so for a
    /// shrinking in-frame change on the reverse strand `protein_start` is the
    /// *end* of the affected range and the residues begin one codon earlier.
    /// Over 37,122 in-frame ClinVar rows the replaced residues sat at
    /// `protein_start` in 71.1% of cases, at `protein_start - 1` in 11.9%, and
    /// at neither in 17.0%.
    ///
    /// Anything that needs the true anchor - HGVSp, in particular - has to
    /// corroborate it against the transcript peptide rather than trust either
    /// scalar. See `fastvep_hgvs::protein` and issue #89.
    pub amino_acids: Option<(String, String)>,
    pub codons: Option<(String, String)>,
    pub exon: Option<(u32, u32, u32)>,
    pub intron: Option<(u32, u32, u32)>,
    pub distance: Option<i64>,
}

/// Whether the transcript's annotation claims a complete initiator codon.
///
/// A `cds_start_NF` transcript does not, and neither does one whose CDS is
/// annotated as beginning part-way through a codon - the three bases at
/// `cdna_coding_start` are then the tail of a codon that starts before them.
/// Ensembl checks the flag before every start-codon predicate; the phase is the
/// same fact in a different field, and without it every length-changing variant
/// reaching that end would be `start_lost`, because those three bases never
/// read as ATG.
fn start_codon_known(transcript: &Transcript) -> bool {
    transcript.codon_table_start_phase == 0 && !transcript.flags.iter().any(|f| f == "cds_start_NF")
}

/// Resolve a reference codon against the transcript's annotated peptide.
///
/// The codon table cannot see an `initial_met` edit or selenocysteine
/// readthrough. The cached peptide can. Keep the correction deliberately
/// limited to those two source annotations; other peptide/translation
/// disagreements require their own causal replay before they may affect output.
///
/// VEP applies these edits only to the reference allele, including when an
/// alternate window contains an unchanged copy of the edited codon.
fn resolve_annotated_residue(transcript: &Transcript, codon_index: usize, translated: u8) -> u8 {
    let annotated = transcript
        .peptide
        .as_deref()
        .and_then(|peptide| peptide.as_bytes().get(codon_index).copied());
    match annotated {
        Some(b'M') if codon_index == 0 => b'M',
        Some(residue) if translated == b'*' && residue != b'*' => residue,
        _ => translated,
    }
}

/// Drop the common prefix and the common suffix of two sequences, Ensembl's
/// `Bio::EnsEMBL::Variation::Utils::Sequence::trim_sequences`.
///
/// Used to ask whether one sequence is an interior slice of the other: if
/// nothing is left of the shorter side, everything it contained matched in
/// place.
fn trim_common_ends<'a>(a: &'a [u8], b: &'a [u8]) -> (&'a [u8], &'a [u8]) {
    let mut lo = 0;
    while lo < a.len() && lo < b.len() && a[lo] == b[lo] {
        lo += 1;
    }
    let mut hi = 0;
    while hi < a.len() - lo && hi < b.len() - lo && a[a.len() - 1 - hi] == b[b.len() - 1 - hi] {
        hi += 1;
    }
    (&a[lo..a.len() - hi], &b[lo..b.len() - hi])
}

impl AlleleConsequenceResult {
    /// The affected residue span as `(lo, hi)`, whichever way the transcript
    /// runs.
    ///
    /// `protein_start`/`protein_end` are strand-ordered, so on the reverse
    /// strand the raw pair arrives reversed - a 3 bp deletion spanning residues
    /// 3 and 4 comes back as `protein_start = 4`, `protein_end = 3`. Every
    /// consumer that wants a span therefore has to sort the pair, and each one
    /// that open-coded it was one `min`/`max` away from an off-by-one that
    /// nothing downstream could detect.
    ///
    /// A single present coordinate describes a single residue.
    pub fn protein_range(&self) -> Option<(u64, u64)> {
        match (self.protein_start, self.protein_end) {
            (Some(s), Some(e)) => Some((s.min(e), s.max(e))),
            (Some(s), None) => Some((s, s)),
            (None, Some(e)) => Some((e, e)),
            (None, None) => None,
        }
    }
}

/// Full prediction result for a variant.
#[derive(Debug, Clone)]
pub struct PredictionResult {
    pub transcript_consequences: Vec<TranscriptConsequence>,
    pub most_severe: Option<Consequence>,
}

/// The consequence prediction engine.
/// A ref/alt pair rendered for display: amino acids (`("T", "M")`) or codons
/// (`("aCg", "aTg")`). `None` when the change does not produce one.
pub type DisplayPair = Option<(String, String)>;

/// What a change inside the CDS resolves to.
///
/// `additional` carries every SO term beyond the first. Ensembl evaluates each
/// predicate in `VariationEffect.pm` independently and keeps all that hold, so
/// one delins can be both `stop_gained` and `protein_altering_variant` - 164 of
/// the 1,803 coding rows the ClinVar 2-star in-frame delins produce are exactly
/// that pair.
///
/// This was one `Option` slot, on the stated grounds that no shape yields more
/// than two coding terms. Three do: a change that replaces the initiator and
/// both introduces a stop and shifts the frame fires `stop_gained`,
/// `frameshift_variant` *and* `start_lost`, and `start_lost` - the lowest-ranked
/// of the three - was dropped without a diagnostic. The `Vec` costs nothing
/// extra on the hot path because `terms_for_window` already builds one and this
/// now takes ownership of it rather than copying out of it.
pub struct CodingChange {
    pub consequence: Consequence,
    pub additional: Vec<Consequence>,
    pub amino_acids: DisplayPair,
    pub codons: DisplayPair,
}

/// The codon window a CDS-internal change edits, translated on both sides.
///
/// This is the whole of Ensembl's coding model. Every term in
/// `VariationEffect.pm` from `start_lost` to `synonymous_variant` is decided
/// by comparing these two peptides and these two codon strings, and every
/// `Amino_acids` / `Codons` value VEP prints is one of them. Building it once
/// and reading the terms off it is not a refactor of convenience: the four
/// shapes that used to have their own arithmetic (SNV/MNV, pure insertion,
/// pure deletion, delins) each got a different subset of the model right,
/// and the disagreements were exactly where the subsets differed.
struct CodonWindow {
    /// The window's residues before and after the change. `X` stands for a
    /// trailing codon the window does not complete, which is how a frameshift's
    /// alternate peptide ends. Either side is empty when the window is - an
    /// insertion between two codons replaces no residues - and stays empty here
    /// rather than becoming `-`, because every predicate below tests it as a
    /// string and `-` would be a residue to them. Ensembl draws the same line:
    /// `_get_peptide_alleles` maps its `-` to `''` before any predicate sees it.
    ref_aas: String,
    alt_aas: String,
    /// Ambiguous input alleles have undefined peptides, unlike translated X residues.
    peptides_defined: bool,
    alt_dna_unambiguous: bool,
    ref_window: Vec<u8>,
    alt_window: Vec<u8>,
    /// Rendered with the replaced span uppercase; `-` for an empty side.
    ref_codons: String,
    alt_codons: String,
    /// CDS bases the reference covers, and bases the alternate supplies.
    ref_len: usize,
    alt_len: usize,
    /// 1-based protein positions of the window's first and last codons.
    /// Ensembl's `translation_start` / `translation_end`, so an empty window -
    /// an insertion between two codons - has `tl_end == tl_start - 1`.
    tl_start: usize,
    tl_end: usize,
    /// Whether the codon the change starts in runs off the end of the CDS -
    /// Ensembl's `partial_codon`. It is a property of the window's *first*
    /// codon, not its last: a window that merely ends on an incomplete codon is
    /// still a codon edit, and Ensembl reports `inframe_deletion` for one.
    partial_codon: bool,
}

impl CodingChange {
    /// The common case: one term, with whatever ref/alt pairs go with it.
    fn single(consequence: Consequence, amino_acids: DisplayPair, codons: DisplayPair) -> Self {
        Self {
            consequence,
            additional: Vec::new(),
            amino_acids,
            codons,
        }
    }
}

pub struct ConsequencePredictor {
    pub upstream_distance: u64,
    pub downstream_distance: u64,
    codon_table: CodonTable,
    /// Vertebrate mitochondrial codon table (NCBI translation table 2), used
    /// instead of `codon_table` whenever the transcript being predicted is on
    /// the mitochondrial chromosome (see [`is_mitochondrial`]). Built once
    /// up front so per-allele translation doesn't reconstruct the table.
    mt_codon_table: CodonTable,
    mature_mirna_ranges: Arc<HashMap<String, Vec<(u64, u64)>>>,
}

impl ConsequencePredictor {
    pub fn new(upstream_distance: u64, downstream_distance: u64) -> Self {
        Self {
            upstream_distance,
            downstream_distance,
            codon_table: CodonTable::standard(),
            mt_codon_table: mitochondrial_codon_table(),
            mature_mirna_ranges: Arc::new(HashMap::new()),
        }
    }

    pub fn with_mature_mirna_ranges(mut self, ranges: HashMap<String, Vec<(u64, u64)>>) -> Self {
        self.mature_mirna_ranges = Arc::new(ranges);
        self
    }

    /// Output-only codon window for a non-minimal biallelic replacement.
    /// Do not feed this window back into consequence or HGVS construction.
    pub fn display_coding_change(
        &self,
        position: &GenomicPosition,
        reference: &Allele,
        alternate: &Allele,
        transcript: &Transcript,
        preserve_input: bool,
    ) -> Option<CodingChange> {
        if reference.len() <= 1 || alternate.len() <= 1 {
            return None;
        }
        let (pos, r, a) = vep_input_variant(position, reference, alternate, preserve_input);
        if r == *reference && a == *alternate {
            return None;
        }
        let start = transcript.genomic_to_cdna(pos.start);
        let end = transcript.genomic_to_cdna(pos.end);
        self.predict_coding_consequence(
            &r,
            &a,
            transcript,
            start.and_then(|p| transcript.cdna_to_cds(p)),
            end.and_then(|p| transcript.cdna_to_cds(p)),
            start,
            end,
        )
    }

    /// Select the codon table to use for a given transcript: the vertebrate
    /// mitochondrial table (NCBI table 2) for MT transcripts, the standard
    /// nuclear table otherwise. AGA/AGG (Arg->Stop), ATA (Ile->Met), and TGA
    /// (Stop->Trp) all differ between the two, so using the wrong table on an
    /// MT variant silently mis-predicts stop_gained/missense/synonymous.
    fn codon_table_for(&self, transcript: &Transcript) -> &CodonTable {
        if is_mitochondrial(&transcript.chromosome) {
            &self.mt_codon_table
        } else {
            &self.codon_table
        }
    }

    /// Predict consequences of a variant against a set of transcripts.
    pub fn predict(
        &self,
        position: &GenomicPosition,
        ref_allele: &Allele,
        alt_alleles: &[Allele],
        transcripts: &[&Transcript],
        ref_seq: Option<&[u8]>,
    ) -> PredictionResult {
        self.predict_with_parsed_input(position, ref_allele, alt_alleles, transcripts, ref_seq, false)
    }

    /// Preserve input coordinates when a parser has already completed its
    /// minimization pass. Repeating it after case conversion is not idempotent.
    pub fn predict_with_parsed_input(
        &self,
        position: &GenomicPosition,
        ref_allele: &Allele,
        alt_alleles: &[Allele],
        transcripts: &[&Transcript],
        ref_seq: Option<&[u8]>,
        input_parsed: bool,
    ) -> PredictionResult {
        let mut transcript_consequences = Vec::new();

        for transcript in transcripts {
            let tc =
                self.predict_transcript(position, ref_allele, alt_alleles, transcript, ref_seq, input_parsed);
            transcript_consequences.push(tc);
        }

        let all_consequences: Vec<Consequence> = transcript_consequences
            .iter()
            .flat_map(|tc| {
                tc.allele_consequences
                    .iter()
                    .flat_map(|ac| ac.consequences.iter().copied())
            })
            .collect();

        let most_severe = Consequence::most_severe(&all_consequences);

        PredictionResult {
            transcript_consequences,
            most_severe,
        }
    }

    fn predict_transcript(
        &self,
        position: &GenomicPosition,
        ref_allele: &Allele,
        alt_alleles: &[Allele],
        transcript: &Transcript,
        ref_seq: Option<&[u8]>,
        input_parsed: bool,
    ) -> TranscriptConsequence {
        let allele_consequences: Vec<AlleleConsequenceResult> = alt_alleles
            .iter()
            // VEP discards repeated REF after parsing; retain raw ALT count for normalization.
            .filter(|alt| *alt != ref_allele)
            .map(|alt| {
                self.predict_allele(
                    position,
                    ref_allele,
                    alt,
                    input_parsed || alt_alleles.len() > 1,
                    transcript,
                    ref_seq,
                )
            })
            .collect();

        TranscriptConsequence {
            transcript_id: transcript.stable_id.clone(),
            gene_id: transcript.gene.stable_id.clone(),
            gene_symbol: transcript.gene.symbol.clone(),
            biotype: transcript.biotype.clone(),
            allele_consequences,
            canonical: transcript.canonical,
            strand: transcript.strand,
        }
    }

    fn predict_allele(
        &self,
        position: &GenomicPosition,
        ref_allele: &Allele,
        alt_allele: &Allele,
        preserve_input: bool,
        transcript: &Transcript,
        _ref_seq: Option<&[u8]>,
    ) -> AlleleConsequenceResult {
        // `Allele::Missing` covers the VCF placeholder alleles that assert no
        // alternate sequence at this site: `*` (spanning upstream deletion)
        // and the non-variant forms `.` / `<NON_REF>` / `<*>` that the VCF
        // reader maps here. There is no alternate sequence to translate, so
        // emitting any consequence for them is meaningless. VEP likewise
        // gives `*` no consequence.
        if *alt_allele == Allele::Missing {
            return AlleleConsequenceResult {
                allele: alt_allele.clone(),
                normalized_position: position.clone(),
                normalized_ref_allele: ref_allele.clone(),
                normalized_alt_allele: alt_allele.clone(),
                consequences: Vec::new(),
                frameshift: false,
                impact: Impact::Modifier,
                cdna_start: None,
                cdna_end: None,
                cds_start: None,
                cds_end: None,
                protein_start: None,
                protein_end: None,
                amino_acids: None,
                codons: None,
                exon: None,
                intron: None,
                distance: None,
            };
        }

        let uploaded_allele = alt_allele.clone();
        let (normalized_position, normalized_ref_allele, normalized_alt_allele) =
            minimize_allele(position, ref_allele, alt_allele);
        let var_start = position.start;
        let var_end = position.end;
        // VEP 115 minimises a biallelic non-minimal indel before it compares
        // the variant with transcript and splice boundaries. Keep FastVEP's
        // established consequence/HGVS coordinates below, but use the VEP
        // parser's left-first representation for location predicates and
        // displayed positions.
        let (predicate_position, predicate_ref, predicate_alt) =
            vep_input_variant(position, ref_allele, alt_allele, preserve_input);
        let distance_start = predicate_position.start;
        let distance_end = predicate_position.end;
        let tr_start = transcript.start;
        let tr_end = transcript.end;

        let mut consequences = Vec::new();
        let mut cds_start = None;
        let mut cds_end = None;
        let mut protein_start = None;
        let mut protein_end = None;
        let mut amino_acids = None;
        let mut frameshift = false;
        let mut codons = None;
        let mut distance = None;

        // 1. Check if variant overlaps the transcript at all
        let overlaps = distance_start <= tr_end && distance_end >= tr_start;

        if !overlaps {
            // Check upstream/downstream
            let dist = self.distance_to_transcript(distance_start, distance_end, transcript);
            if let Some(d) = dist {
                distance = Some(d);
                let abs_dist = d.unsigned_abs();
                if abs_dist <= self.upstream_distance {
                    match transcript.strand {
                        Strand::Forward => {
                            if distance_end < tr_start {
                                consequences.push(Consequence::UpstreamGeneVariant);
                            } else {
                                consequences.push(Consequence::DownstreamGeneVariant);
                            }
                        }
                        Strand::Reverse => {
                            if distance_start > tr_end {
                                consequences.push(Consequence::UpstreamGeneVariant);
                            } else {
                                consequences.push(Consequence::DownstreamGeneVariant);
                            }
                        }
                    }
                }
            }

            if consequences.is_empty() {
                consequences.push(Consequence::IntergenicVariant);
            }

            let impact = Consequence::worst_impact(&consequences).unwrap_or(Impact::Modifier);
            return AlleleConsequenceResult {
                allele: uploaded_allele,
                normalized_position,
                normalized_ref_allele,
                normalized_alt_allele,
                consequences,
                frameshift: false,
                impact,
                cdna_start: None,
                cdna_end: None,
                cds_start: None,
                cds_end: None,
                protein_start: None,
                protein_end: None,
                amino_acids: None,
                codons: None,
                exon: None,
                intron: None,
                distance,
            };
        }

        // 2. Map to cDNA coordinates
        let cdna_start = transcript.genomic_to_cdna(var_start);
        let cdna_end = transcript.genomic_to_cdna(var_end);

        // 3. Determine exon/intron location. An insertion is the zero-length
        // interval between `var_end` and `var_start`; it belongs to a region
        // only when both flanking bases belong to that same region. Treating
        // the flanks as an inclusive range made exon-boundary insertions both
        // exonic and intronic, unlike VEP.
        let predicate_start = predicate_position.start;
        let predicate_end = predicate_position.end;
        let is_insertion = var_end.checked_add(1) == Some(var_start);
        // Coding prediction still needs the adjacent exon at a boundary even
        // though VEP leaves the insertion's EXON field blank.
        // Use VEP's input coordinates for exon overlap, just as for introns
        // and splice predicates; padding can straddle an exon boundary.
        let reaches_exon = transcript
            .exon_at(predicate_start)
            .or_else(|| transcript.exon_overlapping(predicate_start, predicate_end));
        let exon_info = if is_insertion {
            match (transcript.exon_at(predicate_end), transcript.exon_at(predicate_start)) {
                (Some((rank, total)), Some(right)) if (rank, total) == right => Some((rank, rank, total)),
                _ => None,
            }
        } else {
            // VEP exon_number reports the full range in transcript order.
            transcript.exon_range_overlapping(predicate_start, predicate_end)
        }
        .map(|(first, last, t)| (first as u32 + 1, last as u32 + 1, t as u32));
        let intron_info = if is_insertion {
            match (
                transcript.intron_at(predicate_end),
                transcript.intron_at(predicate_start),
            ) {
                (Some((i, total)), Some(right)) if (i, total) == right => Some((i, i, total)),
                _ => None,
            }
        } else {
            // VEP BaseTranscriptVariation::intron_number reports every
            // overlapped intron, sorted in transcript order.
            transcript.intron_overlapping(predicate_start, predicate_end)
        }
        .map(|(first, last, t)| (first as u32 + 1, last as u32 + 1, t as u32));

        let preclassified_exon = splice::overlaps_exon_for_consequence_predicates(
            transcript,
            predicate_start,
            predicate_end,
        );
        let within_frameshift_intron =
            splice::within_frameshift_intron(transcript, predicate_start, predicate_end);

        // 4. Check splice sites (always check regardless of coding status)
        //
        // Every term comes from one pass over the introns, because Ensembl's
        // suppression rules are not a fallback chain: a donor variant that also
        // hits the 5th base earns both terms, and so does an acceptor variant in
        // the polypyrimidine tract. Only `splice_region_variant` is displaced by
        // the more specific terms, and `splice_effects` applies that itself.
        //
        // Measured over a 6,600-variant ClinVar sample against real VEP 115.1,
        // this took the splice-term error from 1,503 rows of 150,725 to none.
        let splice = splice::splice_effects_preclassified(
            transcript,
            predicate_start,
            predicate_end,
            &predicate_ref,
            &predicate_alt,
            preclassified_exon,
        );
        if splice.donor {
            consequences.push(Consequence::SpliceDonorVariant);
        }
        if splice.acceptor {
            consequences.push(Consequence::SpliceAcceptorVariant);
        }
        if splice.donor_5th_base {
            consequences.push(Consequence::SpliceDonorFifthBaseVariant);
        }
        if splice.donor_region {
            consequences.push(Consequence::SpliceDonorRegionVariant);
        }
        if splice.polypyrimidine_tract {
            consequences.push(Consequence::SplicePolypyrimidineTractVariant);
        }
        if splice.region {
            consequences.push(Consequence::SpliceRegionVariant);
        }

        // 5. Coding vs non-coding transcript
        if transcript.is_coding() {
            let coding_start = transcript.coding_region_start.unwrap_or(0);
            let coding_end = transcript.coding_region_end.unwrap_or(0);

            // Map to CDS coordinates
            if let Some(cs) = cdna_start {
                cds_start = transcript.cdna_to_cds(cs);
            }
            if let Some(ce) = cdna_end {
                cds_end = transcript.cdna_to_cds(ce);
            }
            if let Some(cds_s) = cds_start {
                protein_start = Some(Transcript::cds_to_protein(cds_s));
            }
            if let Some(cds_e) = cds_end {
                protein_end = Some(Transcript::cds_to_protein(cds_e));
            }

            // Which regions the variant *reaches*, not which region it starts
            // in. A variant may span the CDS/UTR boundary, and then it belongs
            // to both: VEP reports the pair, `stop_lost&3_prime_UTR_variant`.
            //
            // Testing only `var_start` - and choosing one branch from the
            // result - made the answer depend on the strand, because the
            // genomic start of a reverse-strand variant is its *last* base in
            // transcript order. A delins over FOXL2's stop codon came back as
            // `3_prime_UTR_variant` alone, MODIFIER where VEP says HIGH, and
            // the mirror case on the forward strand lost `start_lost` the same
            // way (#100). Each region is now its own additive test.
            let reaches_coding_exon = reaches_exon.is_some() || within_frameshift_intron;
            let hits_coding = reaches_coding_exon
                && self.is_in_coding_region(var_start, var_end, coding_start, coding_end);
            let hits_5_utr =
                reaches_coding_exon && self.overlaps_5_utr(var_start, var_end, transcript);
            let hits_3_utr =
                reaches_coding_exon && self.overlaps_3_utr(var_start, var_end, transcript);

            if hits_5_utr {
                consequences.push(Consequence::FivePrimeUtrVariant);
            }
            if hits_3_utr {
                consequences.push(Consequence::ThreePrimeUtrVariant);
            }
            if hits_coding {
                // Coding exonic variant - determine coding consequence
                let coding_conseq = self.predict_coding_consequence(
                    ref_allele, alt_allele, transcript, cds_start, cds_end, cdna_start, cdna_end,
                );
                if let Some(change) = coding_conseq {
                    // Constants.pm filters output terms by genomic allele
                    // shape; predicates called by other predicates still use
                    // mapped CDS lengths, which exclude intervening introns.
                    let genomic_ref_len = ref_allele.len();
                    let genomic_alt_len = alt_allele.len();
                    frameshift = change.consequence == Consequence::FrameshiftVariant
                        || change.additional.contains(&Consequence::FrameshiftVariant);
                    consequences.extend(std::iter::once(change.consequence)
                        .chain(change.additional)
                        .filter(|term| match term {
                            Consequence::FrameshiftVariant => genomic_ref_len != genomic_alt_len,
                            Consequence::InframeInsertion => genomic_alt_len > genomic_ref_len,
                            Consequence::InframeDeletion => genomic_alt_len < genomic_ref_len,
                            Consequence::MissenseVariant => genomic_alt_len == genomic_ref_len,
                            // A resolved coding window can fire no predicate.
                            // Apply VEP's default only after splice/NMD terms.
                            Consequence::IntergenicVariant => false,
                            _ => true,
                        }));
                    amino_acids = change.amino_acids;
                    codons = change.codons;
                } else {
                    consequences.push(Consequence::CodingSequenceVariant);
                }
            }
            // `intron_variant` is not a fallback. Ensembl reports it whenever
            // the variant reaches the intron's interior, beside whatever exonic
            // or splice term it also earned - `splice_donor_variant,intron_variant`
            // is its usual answer for a deletion over a donor site. Reserving it
            // for a variant that earned nothing else dropped it from 1,294 rows
            // of a 6,600-variant ClinVar sample.
            if splice.intronic {
                consequences.push(Consequence::IntronVariant);
            }
        } else {
            // Non-coding transcript
            if self.overlaps_mature_mirna(transcript, predicate_start, predicate_end) {
                consequences.push(Consequence::MatureMirnaVariant);
            // VEP non_coding_exon_variant rechecks actual exon overlap;
            // frameshift-intron preclassification alone is insufficient.
            } else if exon_info.is_some() {
                consequences.push(Consequence::NonCodingTranscriptExonVariant);
            } else {
                consequences.push(Consequence::NonCodingTranscriptVariant);
            }
            if splice.intronic {
                consequences.push(Consequence::IntronVariant);
            }
        }

        // Add NMD_transcript_variant modifier for nonsense_mediated_decay transcripts
        if &*transcript.biotype == "nonsense_mediated_decay" {
            consequences.push(Consequence::NmdTranscriptVariant);
        }

        // VEP feature_ablation and feature_amplification are tier 1: whole
        // transcript losses or exact tandem gains exclude lower consequence
        // tiers. Keep independently calculated position fields.
        if predicate_start <= tr_start && predicate_end >= tr_end {
            if predicate_alt.len() < predicate_ref.len() {
                consequences.clear();
                consequences.push(Consequence::TranscriptAblation);
            } else if !predicate_ref.is_empty()
                && predicate_alt.len() > predicate_ref.len()
                && predicate_alt.len() % predicate_ref.len() == 0
                && predicate_alt.as_bytes().chunks(predicate_ref.len())
                    .all(|copy| copy == predicate_ref.as_bytes())
            {
                // VariationEffect::feature_amplification -> tandem_repeat.
                consequences.clear();
                consequences.push(Consequence::TranscriptAmplification);
            } else {
                // coding_unknown and non_coding_exon_variant both reject
                // complete overlap, independently of genomic allele shape.
                consequences.retain(|term| !matches!(term,
                    Consequence::CodingSequenceVariant | Consequence::NonCodingTranscriptExonVariant));
                if transcript.is_coding() {
                    if &*transcript.biotype == "protein_coding" {
                        consequences.push(Consequence::CodingTranscriptVariant);
                    } else {
                        consequences.retain(|term| *term != Consequence::CodingTranscriptVariant);
                    }
                } else if !consequences.contains(&Consequence::MatureMirnaVariant) {
                    consequences.push(Consequence::NonCodingTranscriptVariant);
                }
            }
        }

        // BaseVariationFeatureOverlapAllele::get_all_OverlapConsequences uses
        // DEFAULT_OVERLAP_CONSEQUENCE only after every predicate, including NMD.
        // Even a transcript-overlapping variant can have no applicable term.
        if consequences.is_empty() {
            consequences.push(Consequence::IntergenicVariant);
        }

        // Deduplicate
        consequences.sort_by_key(|c| c.rank());
        consequences.dedup();

        let impact = Consequence::worst_impact(&consequences).unwrap_or(Impact::Modifier);
        // VEP exposes cDNA/CDS/protein position columns from the input
        // parser's left-first coordinates, independently of the later
        // transcript-most-3' HGVS normalization. Keep the local values above
        // for consequence and HGVSp prediction, then project cDNA for the
        // output-layer position calculation.
        let reported_cdna_start = transcript.genomic_to_cdna(distance_start);
        let reported_cdna_end = transcript.genomic_to_cdna(distance_end);

        AlleleConsequenceResult {
            allele: uploaded_allele,
            normalized_position,
            normalized_ref_allele,
            normalized_alt_allele,
            consequences,
            frameshift,
            impact,
            cdna_start: reported_cdna_start,
            cdna_end: reported_cdna_end,
            cds_start,
            cds_end,
            protein_start,
            protein_end,
            amino_acids,
            codons,
            exon: exon_info,
            intron: intron_info,
            distance,
        }
    }

    /// Match VEP's `within_mature_miRNA`: stored ranges are transcript-relative
    /// cDNA coordinates and may cross exon boundaries, so compare each mapped
    /// exon segment rather than treating the complete range as genomic.
    fn overlaps_mature_mirna(
        &self,
        transcript: &Transcript,
        variant_start: u64,
        variant_end: u64,
    ) -> bool {
        if transcript.biotype.as_ref() != "miRNA" {
            return false;
        }
        let Some(ranges) = self.mature_mirna_ranges.get(transcript.stable_id.as_ref()) else {
            return false;
        };

        let mut exons = transcript.exons.iter().collect::<Vec<_>>();
        match transcript.strand {
            Strand::Forward => exons.sort_by_key(|exon| exon.start),
            Strand::Reverse => exons.sort_by(|left, right| right.start.cmp(&left.start)),
        }
        let mut exon_cdna_start = 1;
        for exon in exons {
            let exon_length = exon.end - exon.start + 1;
            let exon_cdna_end = exon_cdna_start + exon_length - 1;
            for &(range_start, range_end) in ranges {
                let segment_start = range_start.max(exon_cdna_start);
                let segment_end = range_end.min(exon_cdna_end);
                if segment_start > segment_end {
                    continue;
                }
                let (genomic_start, genomic_end) = match transcript.strand {
                    Strand::Forward => (
                        exon.start + segment_start - exon_cdna_start,
                        exon.start + segment_end - exon_cdna_start,
                    ),
                    Strand::Reverse => (
                        exon.end - (segment_end - exon_cdna_start),
                        exon.end - (segment_start - exon_cdna_start),
                    ),
                };
                if variant_start <= genomic_end && variant_end >= genomic_start {
                    return true;
                }
            }
            exon_cdna_start = exon_cdna_end + 1;
        }
        false
    }

    /// Predict the coding consequence (missense, synonymous, frameshift, etc.)
    // Each argument is an independent coordinate or allele; grouping them into
    // a struct would only move the argument list to the call site.
    #[allow(clippy::too_many_arguments)]
    fn predict_coding_consequence(
        &self,
        ref_allele: &Allele,
        alt_allele: &Allele,
        transcript: &Transcript,
        cds_start: Option<u64>,
        cds_end: Option<u64>,
        cdna_start: Option<u64>,
        cdna_end: Option<u64>,
    ) -> Option<CodingChange> {
        // A variant reaching past either end of the CDS is not a codon edit:
        // the bases falling in the UTR belong to no codon, and the frame
        // arithmetic below would count them anyway - which is how a delins
        // over FOXL2's stop codon came back `frameshift_variant` on the
        // forward strand and never reached this function at all on the reverse
        // one. VEP reaches the same conclusion structurally: `cds_start` /
        // `cds_end` are undefined when that end of the variant maps outside
        // the coding region, and `frameshift` returns 0 unless both are
        // defined, so the term comes from its start/stop codon predicates.
        //
        // Test the cDNA spans rather than `cds_start.is_none()`: a CDS
        // coordinate is also None for a flank that falls in an intron, and an
        // indel at an exon edge is an ordinary indel.
        let reaches_past_cds = match (
            cdna_start,
            cdna_end,
            transcript.cdna_coding_start,
            transcript.cdna_coding_end,
        ) {
            (Some(s), Some(e), Some(coding_s), Some(coding_e)) => {
                s.min(e) < coding_s || s.max(e) > coding_e
            }
            _ => false,
        };
        // TranscriptMapper preserves a leading intronic gap, so VEP's
        // translation_start is undefined there. A UTR endpoint still has
        // cDNA coordinates and is handled within its mapped exon instead.
        let mapped_translation_start = match transcript.strand {
            Strand::Forward => cdna_start,
            Strand::Reverse => cdna_end,
        }.and_then(|cdna| transcript.cdna_to_cds(cdna));
        let partial_codon = mapped_translation_start.is_some() && transcript
            .translateable_seq
            .as_deref()
            .zip(cds_start.into_iter().chain(cds_end).min())
            .is_some_and(|(seq, cds_s)| {
                // VEP `partial_codon` uses transcript-order translation
                // start, including a sole mapped reverse-strand endpoint.
                let codon_start = (cds_s.saturating_sub(1) as usize / 3) * 3;
                codon_start < seq.len() && codon_start + 3 > seq.len()
            });
        if reaches_past_cds {
            let mut change = self.predict_cds_boundary_consequence(
                ref_allele, alt_allele, transcript, cdna_start, cdna_end,
            );
            // VEP partial_codon is independent of the UTR boundary predicate.
            if partial_codon {
                // stop_retained rejects partial_codon before its CDS/UTR
                // fallback; stop_lost remains an independent predicate.
                if let Some(change) = &mut change {
                    if change.consequence == Consequence::StopRetainedVariant {
                        change.consequence = Consequence::CodingSequenceVariant;
                    }
                }
                change
                    .get_or_insert_with(|| {
                        CodingChange::single(Consequence::CodingSequenceVariant, None, None)
                    })
                    .additional
                    .push(Consequence::IncompleteTerminalCodonVariant);
            }
            return change;
        }
        // One codon window answers every shape. A CDS-internal change is a
        // codon edit whatever its allele string looks like: an SNV, a
        // multi-nucleotide substitution, a pure insertion, a pure deletion and a
        // delins all differ only in how many bases each side of the window
        // contributes. The four separate code paths this replaced disagreed with
        // real VEP 115.1 on 16,000 `Amino_acids` values and 18,000 `Codons`
        // values over a 6,600-variant ClinVar sample, each in its own way.
        let Some(window) =
            self.build_codon_window(transcript, ref_allele, alt_allele, cds_start, cds_end)
        else {
            // The window ran off the end of an incomplete CDS: Ensembl's
            // `partial_codon`, which suppresses every codon term and reports
            // `incomplete_terminal_codon_variant` beside
            // `coding_sequence_variant`. Anything else that has no window - a
            // change not contiguous in CDS space - is `coding_sequence_variant`
            // alone, which is what `None` gives the caller.
            return partial_codon.then_some(CodingChange {
                consequence: Consequence::CodingSequenceVariant,
                additional: vec![Consequence::IncompleteTerminalCodonVariant],
                amino_acids: None,
                codons: None,
            });
        };

        // Ensembl's `_overlaps_start_codon`: a cDNA-space overlap with the three
        // bases at `cdna_coding_start`. A `cds_start_NF` transcript, or one whose
        // CDS begins part-way through a codon, does not claim to carry an
        // initiator, so no term about one is available for it.
        //
        // Ensembl's CDS coordinates run in transcript order on both strands;
        // ours are strand-ordered, so a reverse-strand deletion arrives with its
        // ends swapped and the overlap has to be asked of the sorted pair. An
        // insertion is the exception both ways round: it is the zero-length
        // interval `end = start - 1`, and keeping it inverted is what makes an
        // insertion count only where it sits *inside* the initiator.
        // The window already resolves Mapper::map_insert's surviving exon
        // endpoint. Requiring two directly mapped endpoints loses split starts.
        let overlaps_initiator =
            start_codon_known(transcript) && window.tl_start == 1 && window.tl_end >= 1;

        // Ensembl's `_peptide` is the protein without its terminator; ours
        // carries it, because that is what the annotation's own translation
        // ends with.
        let peptide = transcript.reference_peptide.as_deref().or(transcript.peptide.as_deref());
        // Retained-reference deletions can preserve the full CDS suffix too.
        // Reuse the UTR + CDS edit used by VEP's boundary predicates instead
        // of restricting this check to an empty alternate allele.
        let start_deletion_retains_cds = overlaps_initiator
            && window.alt_len < window.ref_len
            && self.predict_cds_boundary_consequence(
                ref_allele, alt_allele, transcript, cdna_start, cdna_end,
            ).is_some_and(|change| change.consequence == Consequence::StartRetainedVariant
                || change.additional.contains(&Consequence::StartRetainedVariant));
        // VEP tests the edited 5' UTR + CDS before declaring a shorter CDS
        // altered. A following G can replace deleted c.3 while preserving ATG.
        let utr_start_retained = overlaps_initiator
            && window.ref_len != window.alt_len
            && transcript
                .cdna_coding_start
                .filter(|&start| start > 1)
                .zip(if *ref_allele == Allele::Deletion {
                    fastvep_genome::insertion_point(cdna_start, cdna_end, transcript.strand)
                        .and_then(|point| Some((point.checked_add(1)?, point)))
                } else {
                    cdna_start.zip(cdna_end).map(|(a, b)| (a.min(b), a.max(b)))
                })
                .is_some_and(|(start, (lo, hi))| {
                    lo >= start
                        && transcript.cdna_coding_end.is_some_and(|end| {
                            i128::from(end) + alt_allele.len() as i128
                                - (i128::from(hi) - i128::from(lo) + 1)
                                >= i128::from(start) + 2
                        })
                        && self.edited_codon(transcript, lo, hi, alt_allele, start) == Some(*b"ATG")
                });
        // _get_peptide_alleles returns neither peptide if either is undefined.
        // stop_lost/stop_retained then use the independent edited-DNA helper.
        let stop_fallback = (!window.peptides_defined && window.alt_dna_unambiguous)
            .then(|| self.predict_cds_boundary_consequence(
                ref_allele, alt_allele, transcript, cdna_start, cdna_end,
            ))
            .flatten()
            .and_then(|change| std::iter::once(change.consequence).chain(change.additional)
                .find(|term| matches!(term, Consequence::StopLost | Consequence::StopRetainedVariant)));
        Some(self.terms_for_window(
            &window,
            overlaps_initiator,
            peptide,
            start_deletion_retains_cds,
            utr_start_retained,
            transcript.spliced_seq.is_some()
                && transcript.cdna_coding_start.is_some_and(|start| start > 1),
            stop_fallback,
        ))
    }

    /// The SO terms Ensembl derives from a codon window.
    ///
    /// A port of the predicate graph in `VariationEffect.pm` release/115. The
    /// order matters and is not a severity order: `stop_retained` gates both
    /// `frameshift` (l. 1435) and `stop_gained` (l. 1208), `start_lost` gates
    /// `inframe_insertion` (l. 1100) and `protein_altering_variant` (l. 375),
    /// and `protein_altering_variant` defers to almost everything. Reading them
    /// as a severity chain - pick the worst and stop - loses the second term on
    /// the ~1,100 rows per 6,600 ClinVar variants that earn two.
    ///
    /// `overlaps_initiator` is Ensembl's `_overlaps_start_codon`, which the
    /// caller computes because it is a coordinate question rather than a
    /// peptide one.
    fn terms_for_window(
        &self,
        w: &CodonWindow,
        overlaps_initiator: bool,
        peptide: Option<&str>,
        start_deletion_retains_cds: bool,
        utr_start_retained: bool,
        has_five_prime_utr: bool,
        stop_fallback: Option<Consequence>,
    ) -> CodingChange {
        let (ref_pep, alt_pep) = (w.ref_aas.as_str(), w.alt_aas.as_str());
        let extends = |pep: &str| pep.starts_with(ref_pep) || pep.ends_with(ref_pep);

        // `stop_lost` and `stop_retained` are asked first because everything
        // below defers to one or the other.
        let stop_lost = (w.peptides_defined && ref_pep.contains('*') && !alt_pep.contains('*'))
            || stop_fallback == Some(Consequence::StopLost);
        // `ref_eq_alt_sequence` (l. 1321). It matters because `frameshift` and
        // `stop_gained` both defer to `stop_retained`: an insertion that leaves
        // the terminator where it was is an `inframe_insertion` to Ensembl
        // however many bases it adds.
        //
        // The first clause asks whether the replacement keeps the single
        // reference residue at its start and introduces a terminator anywhere.
        // It can suppress an otherwise apparent frameshift, but reproducing it
        // is required for VEP 115 compatibility.
        let ref_matches_alt_start_with_stop =
            ref_pep.len() == 1 && alt_pep.starts_with(ref_pep) && alt_pep.contains('*');

        // The remaining clauses test the annotated terminator. One asks whether
        // it sits at the same residue on both sides. The other is VEP's exact
        // `ref_eq_alt_sequence` test: edit the complete reference peptide, then
        // ask whether its original-length prefix is unchanged and only one or
        // two residues were appended. This can hold before the last residue in
        // a repeat (for example inserting Glu into a poly-Glu tail), so position
        // alone is insufficient. Compare iterators to avoid allocating a full
        // mutant peptide for every coding indel.
        let preserves_reference_with_short_tail = |annotated: &str| {
            let reference = annotated.strip_suffix('*').unwrap_or(annotated).as_bytes();
            let Some(start) = w.tl_start.checked_sub(1) else {
                return false;
            };
            let replaced = if w.tl_end >= w.tl_start {
                w.tl_end - w.tl_start + 1
            } else {
                0
            };
            let Some(suffix_start) = start.checked_add(replaced) else {
                return false;
            };
            if start > reference.len() {
                return false;
            }
            // Perl's `substr($sequence, $start, $length) = $replacement`
            // removes only the available suffix when `$length` runs past the
            // end. VEP relies on that for an incomplete terminal codon.
            let suffix_start = suffix_start.min(reference.len());
            let mutated_len = start + alt_pep.len() + reference.len() - suffix_start;
            mutated_len > reference.len()
                && mutated_len - reference.len() < 3
                && reference[..start]
                    .iter()
                    .chain(alt_pep.as_bytes())
                    .chain(&reference[suffix_start..])
                    .take(reference.len())
                    .eq(reference.iter())
        };
        // `stop_retained` also declines on an incomplete terminal codon.
        let stop_retained = !stop_lost && !w.partial_codon
            && (stop_fallback == Some(Consequence::StopRetainedVariant) || (w.peptides_defined
            && !alt_pep.is_empty()
            && (ref_matches_alt_start_with_stop
                // VEP ref_eq_alt_sequence accepts a leading inserted stop
                // after the full peptide, even with more than two new residues.
                || (alt_pep.starts_with('*') && peptide.is_some_and(|sequence| {
                    w.tl_start > sequence.strip_suffix('*').unwrap_or(sequence).len()
                }))
                || (ref_pep.contains('*') && ref_pep.find('*') == alt_pep.find('*'))
                || peptide.is_some_and(preserves_reference_with_short_tail))));

        // `frameshift` and `inframe_deletion` both decline when the codon the
        // change starts in is incomplete; `protein_altering_variant` does not.
        let frameshift =
            !stop_retained && !w.partial_codon && !w.alt_len.abs_diff(w.ref_len).is_multiple_of(3);
        let stop_gained = w.peptides_defined && !stop_retained && alt_pep.contains('*') && !ref_pep.contains('*');

        // Ensembl withholds the codon strings from a frameshift
        // (`_get_codon_alleles` returns nothing for one), so neither in-frame
        // term can hold there.
        let (mut inframe_insertion, inframe_deletion) = if frameshift {
            (false, false)
        } else if w.alt_len > w.ref_len && w.alt_window.len() > w.ref_window.len() {
            // Constants.pm also requires the DNA-level insertion predicate:
            // UTR extension can enlarge a partial codon after a deletion.
            // `inframe_insertion` cuts the alternate peptide back to its first
            // terminator before the prefix/suffix test, so an insertion that
            // also introduces a stop is still an insertion.
            let alt_to_stop = match alt_pep.find('*') {
                Some(i) => &alt_pep[..=i],
                None => alt_pep,
            };
            (w.peptides_defined && extends(alt_to_stop), false)
        } else if w.alt_window.len() < w.ref_window.len() {
            // `inframe_deletion` tests the codon strings rather than the
            // peptides: the replacement counts as a deletion only when it
            // reproduces a prefix, a suffix, or a whole-codon interior slice of
            // what it replaced.
            let interior = {
                let (r, a) = trim_common_ends(&w.ref_window, &w.alt_window);
                a.is_empty() && r.len().is_multiple_of(3)
            };
            let matched = !w.partial_codon
                && (w.ref_window.starts_with(&w.alt_window)
                    || w.ref_window.ends_with(&w.alt_window)
                    || interior);
            (false, matched)
        } else {
            (false, false)
        };

        // VEP _ins_del_start_altered first accepts an unchanged UTR and ATG.
        // Without that shortcut it compares the complete CDS-length suffix;
        // retaining ATG alone does not preserve an inserted transcript's CDS.
        let start_altered = !utr_start_retained
            && (w.alt_len < w.ref_len || !w.alt_window.ends_with(&w.ref_window));
        let length_changed = w.ref_len != w.alt_len;
        // `start_retained_variant` is independent of `start_lost`: an insertion
        // can displace the active initiator while preserving the original `ATG`
        // at one end of the enlarged window, and Ensembl reports both terms.
        // VEP 115 also reports both for `c.1del` when the preceding UTR base
        // repeats CDS base one, so its CDS-length suffix remains unchanged.
        // Unlike VEP's negated ambiguity guard, these comparisons positively
        // establish an unchanged ATG or CDS suffix, even with N elsewhere.
        let start_retained = overlaps_initiator
            // VEP's pre-predicate `snp` means equal allele lengths, including MNVs.
            && ((!length_changed && w.alt_window.starts_with(b"ATG"))
                || (length_changed && ((w.alt_len > w.ref_len
                && w.alt_window.ends_with(&w.ref_window))
                || start_deletion_retains_cds
                || utr_start_retained)));
        // Ensembl does not call an insertion before the retained start codon an
        // in-frame insertion; the inserted bases are translated before it.
        if start_retained && alt_pep.ends_with(ref_pep) {
            inframe_insertion = false;
        }
        let start_lost = w.alt_dna_unambiguous && overlaps_initiator
            && ((length_changed && start_altered && !(inframe_insertion || inframe_deletion))
                // VEP continues to the inversion and peptide predicates even
                // when the earlier indel predicate did not establish start loss.
                || (has_five_prime_utr && !utr_start_retained && !w.alt_window.starts_with(b"ATG"))
                || (w.peptides_defined && !ref_pep.is_empty() && !alt_pep.is_empty() && alt_pep != "X" && !extends(alt_pep)));
        inframe_insertion &= !start_lost;

        let protein_altering = w.peptides_defined && ref_pep.len() != alt_pep.len()
            && !ref_pep.starts_with('*')
            && !alt_pep.starts_with('*')
            && !extends(alt_pep)
            && !inframe_deletion
            && !start_lost
            && !frameshift;

        let mut terms: Vec<Consequence> = [
            (stop_gained, Consequence::StopGained),
            (frameshift, Consequence::FrameshiftVariant),
            (stop_lost, Consequence::StopLost),
            (start_lost, Consequence::StartLost),
            (inframe_insertion, Consequence::InframeInsertion),
            (inframe_deletion, Consequence::InframeDeletion),
            (protein_altering, Consequence::ProteinAlteringVariant),
            (start_retained, Consequence::StartRetainedVariant),
            (stop_retained, Consequence::StopRetainedVariant),
        ]
        .into_iter()
        .filter_map(|(fired, term)| fired.then_some(term))
        .collect();

        // A same-length replacement that earned nothing above changed the
        // protein or it did not. Ensembl evaluates `missense_variant` and
        // `coding_unknown` independently: `FX/LX` is both missense (F became L)
        // and unresolved coding (the other residue is X).
        let unresolved = !w.peptides_defined || ref_pep.contains('X') || alt_pep.contains('X');
        if w.peptides_defined && terms.is_empty() && !length_changed {
            if ref_pep != alt_pep && !w.partial_codon {
                terms.push(Consequence::MissenseVariant);
            } else if !unresolved {
                terms.push(Consequence::SynonymousVariant);
            }
        }
        // Synonymous is independent of start_lost (e.g. CTG to CTC).
        if ref_pep == alt_pep && !unresolved && !stop_retained && !start_retained
            && !terms.contains(&Consequence::SynonymousVariant)
        {
            terms.push(Consequence::SynonymousVariant);
        }
        // VEP's `coding_unknown` deliberately does not defer to an in-frame
        // insertion or a missense term, so those can carry this term too.
        let coding_unknown = unresolved
            && !(frameshift
                || inframe_deletion
                || protein_altering
                || start_retained
                || start_lost
                || stop_retained
                || stop_lost);
        if coding_unknown {
            terms.push(Consequence::CodingSequenceVariant);
        }
        // `incomplete_terminal_codon_variant` sits beside whatever else held,
        // and beside `coding_sequence_variant` when nothing else did.
        if w.partial_codon {
            if terms.is_empty() {
                terms.push(Consequence::CodingSequenceVariant);
            }
            terms.push(Consequence::IncompleteTerminalCodonVariant);
        }
        terms.sort_by_key(|c| c.rank());

        // `-` is how an absent side is written, for residues as for codons.
        let shown = |pep: &str| {
            if pep.is_empty() {
                "-".to_string()
            } else {
                pep.to_string()
            }
        };
        let amino_acids = w.peptides_defined.then(|| (shown(&w.ref_aas), shown(&w.alt_aas)));
        let codons = Some((w.ref_codons.clone(), w.alt_codons.clone()));
        if terms.is_empty() {
            // coding_unknown requires missing/ambiguous peptides, not merely
            // an unclassified edit. Preserve the resolved display window and
            // defer the default term until all transcript predicates finish.
            return CodingChange::single(Consequence::IntergenicVariant, amino_acids, codons);
        }
        let first = terms.remove(0);
        CodingChange {
            consequence: first,
            additional: terms,
            amino_acids,
            codons,
        }
    }

    /// Build the codon window for a change whose reference bases are contiguous
    /// in CDS space.
    ///
    /// The window is every codon the reference allele touches. For an insertion
    /// the reference covers no bases at all, so the window is the single codon
    /// it falls inside - and *nothing*, when it falls exactly on a codon
    /// boundary, which is why VEP writes `-/HENKTKGD` for a codon-aligned
    /// insertion rather than repeating the flanking residue on both sides.
    ///
    /// VEP edits the mapped CDS span, including when the genomic reference
    /// includes an intron. Returns `None` when coordinates or sequence are absent.
    fn build_codon_window(
        &self,
        transcript: &Transcript,
        ref_allele: &Allele,
        alt_allele: &Allele,
        cds_start: Option<u64>,
        cds_end: Option<u64>,
    ) -> Option<CodonWindow> {
        let seq = transcript.translateable_seq.as_deref()?.as_bytes();
        // One end may be missing for an insertion at an exon edge; the branch
        // below reads the pair itself rather than through these.
        let (cds_lo, cds_hi) = match (cds_start, cds_end) {
            (Some(s), Some(e)) => (s.min(e), s.max(e)),
            _ => (0, 0),
        };

        // `first` is the 0-based CDS offset the change starts at; for an
        // insertion it is the offset it is inserted *before*, so a change of no
        // reference bases still has a position.
        let (first, ref_len) = if *ref_allele == Allele::Deletion {
            // Ensembl's zero-length interval: `cds_start == cds_end + 1`, so the
            // insertion goes after the lower of the two.
            //
            // An insertion on an exon's edge has one end in the intron, and so
            // only one CDS coordinate. The insertion point is still determined:
            // it abuts the exonic base, on whichever side the strand puts the
            // intron. `cds_start` comes from the genomic left edge and `cds_end`
            // from the right, so the surviving coordinate says which. Without
            // this an insertion at an exon's first base was
            // `coding_sequence_variant` - 39 rows of a 6,600-variant ClinVar
            // sample, LOW where VEP says HIGH, every one a frameshift.
            let point = fastvep_genome::insertion_point(cds_start, cds_end, transcript.strand)?;
            (point as usize, 0usize)
        } else {
            if cds_lo < 1 || cds_hi - cds_lo + 1 > ref_allele.len() as u64 {
                return None;
            }
            let ref_len = usize::try_from(cds_hi - cds_lo + 1).ok()?;
            ((cds_lo - 1) as usize, ref_len)
        };

        // Alt bases in transcript orientation. On the reverse strand the VCF
        // gives them in genomic order, so they are reversed *as well as*
        // complemented: complementing in place reported HPS4
        // c.1060_1061delTCinsAG as `stop_gained` when it is synonymous.
        let alt_cds: Vec<u8> = match alt_allele {
            Allele::Sequence(bases) => match transcript.strand {
                Strand::Forward => bases.clone(),
                Strand::Reverse => bases.iter().rev().map(|&b| complement(b)).collect(),
            },
            _ => Vec::new(),
        };

        let win_start = (first / 3) * 3;
        let win_end = if ref_len == 0 {
            // An insertion on a codon boundary belongs to no codon.
            if first.is_multiple_of(3) {
                win_start
            } else {
                win_start + 3
            }
        } else {
            ((first + ref_len - 1) / 3 + 1) * 3
        };
        // The last codon of an incomplete CDS is short, and the window ends
        // wherever the sequence does: Ensembl reports `aAGA/a` for a deletion
        // there, translating the leftover as `X` rather than declining to
        // describe the change. Refusing the window instead left 3 rows per 230
        // IMPACT-changing ClinVar variants as `coding_sequence_variant` where
        // VEP says `inframe_deletion` or `frameshift_variant`.
        let nominal_win_end = win_end;
        let mut utr_bases = win_end.saturating_sub(seq.len());
        let win_end = win_end.min(seq.len());
        if win_start > win_end {
            return None;
        }

        let ref_window = seq[win_start..win_end].to_vec();
        let lead = first - win_start; // unchanged bases before the change
        let trail = win_end.saturating_sub(first + ref_len); // unchanged bases after it
        let mut alt_window = Vec::with_capacity(lead + alt_cds.len() + trail);
        if lead > ref_window.len() {
            return None;
        }
        alt_window.extend_from_slice(&ref_window[..lead]);
        alt_window.extend_from_slice(&alt_cds);
        alt_window.extend_from_slice(&ref_window[ref_window.len() - trail..]);
        // VEP's _trim_incomplete_codon returns an empty CDS when the
        // edited whole CDS is shorter than three bases, before appending UTR.
        if seq.len().saturating_sub(ref_len) + alt_cds.len() < 3 {
            alt_window.clear();
            utr_bases = (nominal_win_end - win_start + alt_cds.len()).saturating_sub(ref_len);
        }
        // VEP _get_alternate_cds appends the 3' UTR before codon() slices
        // the alternate window. The reference window stays CDS-only.
        if utr_bases > 0 {
            if let Some(utr) = transcript
                .spliced_seq
                .as_deref()
                .zip(transcript.cdna_coding_end)
                .and_then(|(spliced, end)| spliced.as_bytes().get(end as usize..))
            {
                alt_window.extend_from_slice(&utr[..utr_bases.min(utr.len())]);
            }
        }

        let table = self.codon_table_for(transcript);
        // VEP peptide() applies sequence edits only inside its is_reference
        // branch. Alternate windows always use ordinary codon translation.
        let translate = |window: &[u8], anchored: bool| -> String {
            let mut pep: String = window
                .as_chunks::<3>()
                .0
                .iter()
                .enumerate()
                .map(|(i, codon)| {
                    let index = win_start / 3 + i;
                    // At residue one, an annotated mitochondrial initiator is
                    // methionine even when its ordinary table-2 translation is
                    // Ile (for example human MT-ND2 starts with ATT).
                    let translated = if anchored
                        && index == 0
                        && is_mitochondrial(&transcript.chromosome)
                        && start_codon_known(transcript)
                    {
                        b'M'
                    } else {
                        table.translate(codon)
                    };
                    if anchored {
                        resolve_annotated_residue(transcript, index, translated) as char
                    } else {
                        translated as char
                    }
                })
                .collect();
            // Ensembl marks a codon the window does not complete with `X`, but
            // only when the codons before it did not already end translation:
            // `if($partial_codon && $pep ne '*') { $pep .= 'X' }`
            // (`TranscriptVariationAllele.pm` release/115). So `tAac` is `*` and
            // `gttg` is `VX`, and a peptide that merely *contains* a terminator
            // still gets its `X` - `SIFNYIITLFQ*YSFIPYX`.
            if !window.len().is_multiple_of(3) && pep != "*" {
                pep.push('X');
            }
            pep
        };

        // VEP's codon rendering for an edit of any length: the unchanged flanks
        // of the window stay lowercase and exactly the replaced or inserted
        // bases are uppercase, on both sides - `gGg/gTCCCg` for one base
        // replaced by four, and `TAT/CAC` for a whole codon replaced even though
        // its middle base did not change. That is not the same as "uppercase
        // wherever the two differ", which is what produced `TaT/CaC`.
        let render = |window: &[u8], upper_len: usize| -> String {
            if window.is_empty() {
                return "-".to_string();
            }
            window
                .iter()
                .enumerate()
                .map(|(i, &b)| {
                    if i >= lead && i < lead + upper_len {
                        (b as char).to_ascii_uppercase()
                    } else {
                        (b as char).to_ascii_lowercase()
                    }
                })
                .collect()
        };

        Some(CodonWindow {
            peptides_defined: ref_allele.is_unambiguous_dna() && alt_allele.is_unambiguous_dna(),
            alt_dna_unambiguous: alt_allele.is_unambiguous_dna(),
            ref_aas: translate(&ref_window, true),
            alt_aas: translate(&alt_window, false),
            // VEP display_codon uppercases feature_seq length (genomic),
            // independently of the mapped CDS length used for reconstruction.
            ref_codons: render(&ref_window, ref_allele.len()),
            alt_codons: render(&alt_window, alt_cds.len()),
            tl_start: win_start / 3 + 1,
            tl_end: win_end.div_ceil(3),
            partial_codon: win_start + 3 > seq.len(),
            ref_window,
            alt_window,
            ref_len,
            alt_len: alt_cds.len(),
        })
    }

    /// Consequence for a coding variant that reaches past one end of the CDS.
    ///
    /// Such a variant almost always touches the initiator or the terminator:
    /// reaching past `cdna_coding_start` while still overlapping the CDS means
    /// covering part of `[cdna_coding_start, cdna_coding_start + 2]`, and the
    /// same holds at the other end. So the question is usually only whether
    /// that codon survives the change.
    ///
    /// It is not a guarantee, because the overlap that put us here is measured
    /// in genomic space and can be satisfied by intronic bases alone. Returning
    /// `None` then is correct rather than merely safe: the caller falls back to
    /// `coding_sequence_variant`, which is what a coding change nobody can
    /// resolve any further is.
    ///
    /// The amino-acid and codon pair is deliberately left unset. A variant of
    /// this shape has no single replaced codon to name - part of it is not in
    /// the CDS at all - and the alternative is inventing one: the forward-strand
    /// case used to reach the frameshift formatter and report
    /// `p.Tyr175TerfsTer1` for a change that does not shift any frame.
    fn predict_cds_boundary_consequence(
        &self,
        ref_allele: &Allele,
        alt_allele: &Allele,
        transcript: &Transcript,
        cdna_start: Option<u64>,
        cdna_end: Option<u64>,
    ) -> Option<CodingChange> {
        // VEP's DNA alteration helpers also reject ambiguous alleles. Their
        // negation-based "retained" fallback is not evidence of retention:
        // leave this unresolved instead of asserting either loss or retention.
        if !alt_allele.is_unambiguous_dna() {
            return None;
        }
        let coding_start = transcript.cdna_coding_start?;
        let coding_end = transcript.cdna_coding_end?;
        let (s, e) = (cdna_start?, cdna_end?);
        let (cdna_lo, cdna_hi) = (s.min(e), s.max(e));
        if coding_end < coding_start + 2 {
            return None; // no room for a codon at either end
        }

        // VEP's `_overlaps_stop_codon` / `_overlaps_start_codon`: a cDNA-space
        // overlap with the three bases of the terminator or the initiator. A
        // `cds_end_NF` / `cds_start_NF` transcript is one whose annotation does
        // not claim to carry that codon, so neither term is available for it -
        // VEP checks the same two flags before either predicate.
        let flagged = |flag: &str| transcript.flags.iter().any(|f| f == flag);
        let overlaps_stop =
            cdna_lo <= coding_end && cdna_hi >= coding_end - 2 && !flagged("cds_end_NF");
        let overlaps_start =
            cdna_lo <= coding_start + 2 && cdna_hi >= coding_start && start_codon_known(transcript);

        // VEP evaluates start and stop predicates independently. Its stop
        // fallback declines to call an alteration when the leading cDNA
        // endpoint has no CDS coordinate, even if the edit spans the CDS.
        // An equal-length replacement declines only the stop fallback; it
        // must still evaluate the independent start-retention predicate.
        let stop_consequence = (overlaps_stop && ref_allele.len() != alt_allele.len()).then(|| {
            let still_a_stop = transcript.cdna_to_cds(cdna_lo).is_none()
                || self.edited_codon(transcript, cdna_lo, cdna_hi, alt_allele, coding_end - 2)
                    .is_some_and(|codon| self.codon_table_for(transcript).translate(&codon) == b'*');
            if still_a_stop { Consequence::StopRetainedVariant } else { Consequence::StopLost }
        });

        let start_change = (|| {
        if overlaps_start {
            if ref_allele.len() != alt_allele.len() {
                // VEP `_ins_del_start_altered` edits UTR + CDS. Repeats can
                // preserve both the original UTR and ATG even when the edit
                // crosses their boundary. An unchanged full CDS suffix is
                // also retained. This allocation is confined to that rare
                // boundary path; ordinary coding variants use codon windows.
                let (retained, fixed_start_altered) = transcript
                    .spliced_seq
                    .as_deref()
                    .zip(transcript.translateable_seq.as_deref())
                    .and_then(|(spliced, cds)| {
                        let utr = spliced
                            .as_bytes()
                            .get(..coding_start.checked_sub(1)? as usize)?;
                        let mut edited = [utr, cds.as_bytes()].concat();
                        let range = if *ref_allele == Allele::Deletion {
                            cdna_lo as usize..cdna_lo as usize
                        } else {
                            cdna_lo.checked_sub(1)? as usize..cdna_hi as usize
                        };
                        let range = range.start..range.end.min(edited.len());
                        edited.get(range.clone())?;
                        let alt = match alt_allele {
                            Allele::Sequence(bases) => bases.as_slice(),
                            Allele::Deletion => &[],
                            _ => return None,
                        };
                        match transcript.strand {
                            Strand::Forward => {
                                edited.splice(range, alt.iter().copied());
                            }
                            Strand::Reverse => {
                                edited
                                    .splice(range, alt.iter().rev().map(|&base| complement(base)));
                            }
                        }
                        let atg = edited.get(utr.len()..utr.len() + 3) == Some(b"ATG");
                        Some((
                            (!utr.is_empty() && edited.get(..utr.len()) == Some(utr) && atg)
                                || edited.ends_with(cds.as_bytes()),
                            !utr.is_empty() && !atg,
                        ))
                    })
                    .unwrap_or((false, false));
                let mut change = CodingChange::single(
                    if retained {
                        Consequence::StartRetainedVariant
                    } else {
                        Consequence::StartLost
                    },
                    None,
                    None,
                );
                // VEP start_lost also calls `_inv_start_altered` for indels.
                // Its fixed-offset test is independent of the preserved-CDS
                // suffix above, so both start terms can be reported.
                if retained && fixed_start_altered {
                    change.additional.push(Consequence::StartLost);
                }
                return Some(change);
            }
            // VEP's SNP predicate reads the CDS-length suffix, while its
            // inversion/start-loss predicate reads the fixed original offset.
            // Across an intron, equal genomic REF/ALT lengths can still expand
            // the edited cDNA, so these two codons need not be the same.
            // `_inv_start_altered` requires a 5' UTR and declines edits
            // extending past UTR + CDS. `_snp_start_altered` instead lets
            // substr clip that removal at the CDS end before taking its suffix.
            let fixed_start_altered = coding_start > 1 && cdna_hi <= coding_end && self.edited_codon(
                transcript, cdna_lo, cdna_hi, alt_allele, coding_start,
            )? != *b"ATG";
            let suffix_start = coding_start.checked_add(alt_allele.len() as u64)?
                .checked_sub(cdna_hi.min(coding_end) - cdna_lo + 1)?;
            let retained = self.edited_codon(
                transcript, cdna_lo, cdna_hi, alt_allele, suffix_start,
            )? == *b"ATG";
            if !fixed_start_altered && !retained { return None; }
            let mut change = CodingChange::single(
                if retained { Consequence::StartRetainedVariant }
                else { Consequence::StartLost }, None, None,
            );
            if retained && fixed_start_altered { change.additional.push(Consequence::StartLost); }
            return Some(change);
        }

        None
        })();
        let mut change = start_change.or_else(|| stop_consequence
            .map(|term| CodingChange::single(term, None, None)))?;
        if let Some(stop) = stop_consequence.filter(|term| *term != change.consequence) {
            change.additional.push(stop);
        }
        Some(change)
    }

    /// The three bases at cDNA position `codon_start` *after* `alt_allele`
    /// replaces `[cdna_lo, cdna_hi]`, in transcript orientation.
    ///
    /// VEP decides whether an indel took the terminator or the initiator away
    /// by editing the transcript's own sequence and re-reading one codon at a
    /// fixed offset (`_ins_del_stop_altered`, `_ins_del_start_altered`) rather
    /// than by arithmetic on the variant's length, and it has to: length
    /// arithmetic cannot tell a delins that removes the stop codon from one
    /// that happens to rebuild it, and those are `stop_lost` and
    /// `stop_retained_variant` respectively.
    ///
    /// Read the codon in place instead of materialising the edited sequence.
    /// This runs inside the per-variant loop and the sequence being edited is
    /// the whole transcript, so a copy here would be a copy per variant.
    ///
    /// `None` when the run loaded no sequences, or when the edited cDNA is too
    /// short to reach the codon.
    fn edited_codon(
        &self,
        transcript: &Transcript,
        cdna_lo: u64,
        cdna_hi: u64,
        alt_allele: &Allele,
        codon_start: u64,
    ) -> Option<[u8; 3]> {
        let seq = transcript.spliced_seq.as_deref()?.as_bytes();
        let alt: &[u8] = match alt_allele {
            Allele::Sequence(bases) => bases,
            Allele::Deletion => &[],
            // `*` and the non-variant placeholders assert no alternate
            // sequence, and a symbolic `<DEL>`/`<DUP>` names one without
            // spelling it, so neither has an edited codon to read. Structural
            // alleles are the SV predictor's job.
            Allele::Missing | Allele::Symbolic(_) => return None,
        };
        let alt_len = alt.len() as u64;

        let mut codon = [0u8; 3];
        for (i, slot) in codon.iter_mut().enumerate() {
            let pos = codon_start + i as u64; // 1-based, in the *edited* cDNA
            let base = if pos < cdna_lo {
                *seq.get((pos - 1) as usize)?
            } else if pos < cdna_lo + alt_len {
                // Inside the replacement. The alternate bases arrive in
                // genomic order, so a reverse-strand transcript reads them
                // reverse-complemented, not merely complemented in place.
                let offset = (pos - cdna_lo) as usize;
                match transcript.strand {
                    Strand::Forward => alt[offset],
                    Strand::Reverse => complement(alt[alt.len() - 1 - offset]),
                }
            } else {
                // Past the replacement, so the next original base is the one
                // after the replaced span rather than the one at `pos`.
                let original = cdna_hi + (pos - (cdna_lo + alt_len)) + 1;
                *seq.get((original - 1) as usize)?
            };
            *slot = base.to_ascii_uppercase();
        }
        Some(codon)
    }

    /// Compute amino acids and codons affected by an indel variant.
    /// Returns (amino_acids, codons) tuples.
    /// For frameshifts: ref codon with VEP-style case formatting, truncated alt codon.
    fn distance_to_transcript(
        &self,
        var_start: u64,
        var_end: u64,
        transcript: &Transcript,
    ) -> Option<i64> {
        // VEP uses the shortest absolute distance among all four variant and
        // transcript endpoints. Do not reorder insertion coordinates: their
        // Ensembl representation intentionally has start == end + 1.
        [
            var_start.abs_diff(transcript.start),
            var_start.abs_diff(transcript.end),
            var_end.abs_diff(transcript.start),
            var_end.abs_diff(transcript.end),
        ]
        .into_iter()
        .min()
        .and_then(|distance| i64::try_from(distance).ok())
    }

    fn is_in_coding_region(
        &self,
        var_start: u64,
        var_end: u64,
        coding_start: u64,
        coding_end: u64,
    ) -> bool {
        var_start <= coding_end && var_end >= coding_start
    }

    /// Does the variant's span reach the untranslated region on the low
    /// genomic side of the coding region - `[transcript.start, cds_start - 1]`?
    ///
    /// Insertions widen the span by the base on either side of the insertion
    /// point, which is what makes an insertion sitting exactly on the coding
    /// boundary a UTR variant. VEP carries the same rule as the two special
    /// cases at the top of `_before_coding` / `_after_coding`.
    ///
    /// The span tested is genomic, so it includes any intron lying between the
    /// transcript's edge and the coding region, and a variant whose only bases
    /// there are intronic still counts. That is VEP's rule rather than an
    /// oversight - `_before_coding` and `_after_coding` are plain `overlap`
    /// calls against the same genomic interval, gated only on the variant
    /// touching an exon somewhere - and matching it keeps the UTR terms
    /// comparable with VEP's. Its inclusive overlap also reports a span that
    /// straddles a coding boundary at the transcript edge, even though that
    /// side has no annotated UTR. It can only add a MODIFIER term next to a
    /// correctly-derived one, never change the reported impact.
    fn reaches_low_side_utr(&self, var_start: u64, var_end: u64, transcript: &Transcript) -> bool {
        let Some(coding_start) = transcript.coding_region_start else {
            return false;
        };
        let (lo, hi) = (var_start.min(var_end), var_start.max(var_end));
        lo < coding_start && hi >= transcript.start
    }

    /// The same test on the high genomic side - `[cds_end + 1, transcript.end]`.
    fn reaches_high_side_utr(&self, var_start: u64, var_end: u64, transcript: &Transcript) -> bool {
        let Some(coding_end) = transcript.coding_region_end else {
            return false;
        };
        let (lo, hi) = (var_start.min(var_end), var_start.max(var_end));
        hi > coding_end && lo <= transcript.end
    }

    /// Does the variant's span reach the 5' UTR? Which genomic side that is
    /// depends on the strand; the overlap test itself does not.
    fn overlaps_5_utr(&self, var_start: u64, var_end: u64, transcript: &Transcript) -> bool {
        if transcript.coding_region_start.is_none() || transcript.coding_region_end.is_none() {
            return false;
        }
        match transcript.strand {
            Strand::Forward => self.reaches_low_side_utr(var_start, var_end, transcript),
            Strand::Reverse => self.reaches_high_side_utr(var_start, var_end, transcript),
        }
    }

    /// Does the variant's span reach the 3' UTR?
    fn overlaps_3_utr(&self, var_start: u64, var_end: u64, transcript: &Transcript) -> bool {
        if transcript.coding_region_start.is_none() || transcript.coding_region_end.is_none() {
            return false;
        }
        match transcript.strand {
            Strand::Forward => self.reaches_high_side_utr(var_start, var_end, transcript),
            Strand::Reverse => self.reaches_low_side_utr(var_start, var_end, transcript),
        }
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        other => other,
    }
}

fn minimize_allele(
    position: &GenomicPosition,
    ref_allele: &Allele,
    alt_allele: &Allele,
) -> (GenomicPosition, Allele, Allele) {
    let (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases)) = (ref_allele, alt_allele)
    else {
        return (position.clone(), ref_allele.clone(), alt_allele.clone());
    };

    let mut ref_end = ref_bases.len();
    let mut alt_end = alt_bases.len();
    while ref_end > 0
        && alt_end > 0
        && ref_bases[ref_end - 1].eq_ignore_ascii_case(&alt_bases[alt_end - 1])
    {
        ref_end -= 1;
        alt_end -= 1;
    }

    let mut prefix = 0;
    while prefix < ref_end
        && prefix < alt_end
        && ref_bases[prefix].eq_ignore_ascii_case(&alt_bases[prefix])
    {
        prefix += 1;
    }

    if prefix == 0 && ref_end == ref_bases.len() && alt_end == alt_bases.len() {
        return (position.clone(), ref_allele.clone(), alt_allele.clone());
    }

    let start = position.start + prefix as u64;
    let normalized_ref = if ref_end == prefix {
        Allele::Deletion
    } else {
        Allele::Sequence(ref_bases[prefix..ref_end].to_vec())
    };
    let normalized_alt = if alt_end == prefix {
        Allele::Deletion
    } else {
        Allele::Sequence(alt_bases[prefix..alt_end].to_vec())
    };
    let end = if normalized_ref == Allele::Deletion {
        start.saturating_sub(1)
    } else {
        start + normalized_ref.len() as u64 - 1
    };

    (
        GenomicPosition::new(position.chromosome.clone(), start, end, position.strand),
        normalized_ref,
        normalized_alt,
    )
}

/// Return the coordinates VEP 115's input parser uses for transcript overlap
/// and distance calculations.
///
/// Ordinary biallelic VCF indels are minimised automatically, with the common
/// prefix removed before the common suffix. VEP's rejoined multi-allelic CSQ
/// output retains the uploaded record's shared interval for these fields, so
/// it is intentionally left unchanged here. FastVEP keeps its existing
/// end-first minimisation for consequence and HGVS calculation, where it is
/// needed for transcript-most-3' representation.
fn vep_input_variant(
    position: &GenomicPosition,
    ref_allele: &Allele,
    alt_allele: &Allele,
    preserve_input: bool,
) -> (GenomicPosition, Allele, Allele) {
    let (Allele::Sequence(reference), Allele::Sequence(alternate)) = (ref_allele, alt_allele)
    else {
        return (position.clone(), ref_allele.clone(), alt_allele.clone());
    };

    if preserve_input || reference.len() == alternate.len() {
        return (position.clone(), ref_allele.clone(), alt_allele.clone());
    }

    let mut prefix = 0;
    while prefix < reference.len()
        && prefix < alternate.len()
        && reference[prefix].eq_ignore_ascii_case(&alternate[prefix])
    {
        prefix += 1;
    }

    let mut ref_end = reference.len();
    let mut alt_end = alternate.len();
    while ref_end > prefix
        && alt_end > prefix
        && reference[ref_end - 1].eq_ignore_ascii_case(&alternate[alt_end - 1])
    {
        ref_end -= 1;
        alt_end -= 1;
    }

    if prefix == 0 && ref_end == reference.len() && alt_end == alternate.len() {
        return (position.clone(), ref_allele.clone(), alt_allele.clone());
    }

    let position = GenomicPosition::new(
        position.chromosome.clone(),
        position.start + prefix as u64,
        position
            .end
            .saturating_sub((reference.len() - ref_end) as u64),
        position.strand,
    );
    let reference = if ref_end == prefix {
        Allele::Deletion
    } else {
        Allele::Sequence(reference[prefix..ref_end].to_vec())
    };
    let alternate = if alt_end == prefix {
        Allele::Deletion
    } else {
        Allele::Sequence(alternate[prefix..alt_end].to_vec())
    };

    (position, reference, alternate)
}

/// VEP 115 parser coordinates used for transcript lookup and location predicates.
pub fn vep_input_position(
    position: &GenomicPosition,
    ref_allele: &Allele,
    alt_allele: &Allele,
    preserve_input: bool,
) -> GenomicPosition {
    vep_input_variant(position, ref_allele, alt_allele, preserve_input).0
}

impl Default for ConsequencePredictor {
    fn default() -> Self {
        Self::new(5000, 5000)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn protein_range_sorts_a_reverse_strand_pair() {
        // The trap this exists to remove: on the reverse strand the raw pair
        // arrives with start past end, so a consumer reading `protein_start` as
        // "the first affected residue" is one codon off, and `Protein_position`
        // would print backwards.
        let mut ac = bare_result();

        ac.protein_start = Some(4);
        ac.protein_end = Some(3);
        assert_eq!(ac.protein_range(), Some((3, 4)));

        // Forward strand: already ordered, unchanged.
        ac.protein_start = Some(3);
        ac.protein_end = Some(4);
        assert_eq!(ac.protein_range(), Some((3, 4)));

        // A single coordinate describes a single residue, from either side.
        ac.protein_start = Some(7);
        ac.protein_end = None;
        assert_eq!(ac.protein_range(), Some((7, 7)));
        ac.protein_start = None;
        ac.protein_end = Some(7);
        assert_eq!(ac.protein_range(), Some((7, 7)));

        ac.protein_end = None;
        assert_eq!(ac.protein_range(), None);
    }

    fn bare_result() -> AlleleConsequenceResult {
        AlleleConsequenceResult {
            allele: Allele::from_str("A"),
            normalized_position: GenomicPosition::new("chr1", 1, 1, Strand::Forward),
            normalized_ref_allele: Allele::from_str("C"),
            normalized_alt_allele: Allele::from_str("A"),
            consequences: vec![],
            frameshift: false,
            impact: Impact::Modifier,
            cdna_start: None,
            cdna_end: None,
            cds_start: None,
            cds_end: None,
            protein_start: None,
            protein_end: None,
            amino_acids: None,
            codons: None,
            exon: None,
            intron: None,
            distance: None,
        }
    }

    use fastvep_genome::{Exon, Gene, Translation};

    fn make_coding_transcript() -> Transcript {
        // A simple protein-coding transcript on forward strand:
        // Exon1: 1000-1200 (UTR: 1000-1049, CDS: 1050-1200)
        // Intron: 1201-1999
        // Exon2: 2000-2300 (all CDS)
        // Intron: 2301-3999
        // Exon3: 4000-5000 (CDS: 4000-4500, UTR: 4501-5000)
        //
        // CDS length: 151 + 301 + 501 = 953 bases
        // translateable_seq: from cDNA pos 51 to 953+50=1003
        let translateable = "ATGGCTTCAAAGCCC".to_string() + &"A".repeat(938); // starts with ATG

        Transcript {
            stable_id: "ENST00000001".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG00000001".into(),
                symbol: Some("TESTGENE".into()),
                symbol_source: Some("HGNC".into()),
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start: 1000,
                end: 5000,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start: 1000,
            end: 5000,
            strand: Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "E1".into(),
                    start: 1000,
                    end: 1200,
                    strand: Strand::Forward,
                    phase: -1,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "E2".into(),
                    start: 2000,
                    end: 2300,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 1,
                    rank: 2,
                },
                Exon {
                    stable_id: "E3".into(),
                    start: 4000,
                    end: 5000,
                    strand: Strand::Forward,
                    phase: 1,
                    end_phase: -1,
                    rank: 3,
                },
            ],
            translation: Some(Translation {
                stable_id: "ENSP00000001".into(),
                genomic_start: 1050,
                genomic_end: 4500,
                start_exon_rank: 1,
                start_exon_offset: 50,
                end_exon_rank: 3,
                end_exon_offset: 500,
            }),
            cdna_coding_start: Some(51),
            cdna_coding_end: Some(1003),
            coding_region_start: Some(1050),
            coding_region_end: Some(4500),
            spliced_seq: None,
            translateable_seq: Some(translateable),
            peptide: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: Some(1),
            appris: None,
            ccds: None,
            protein_id: Some("ENSP00000001".into()),
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

    fn make_noncoding_transcript() -> Transcript {
        Transcript {
            stable_id: "ENST_NC".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_NC".into(),
                symbol: Some("NCRNA1".into()),
                symbol_source: None,
                hgnc_id: None,
                biotype: "lncRNA".into(),
                chromosome: "chr1".into(),
                start: 10000,
                end: 12000,
                strand: Strand::Forward,
            },
            biotype: "lncRNA".into(),
            chromosome: "chr1".into(),
            start: 10000,
            end: 12000,
            strand: Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "E1".into(),
                    start: 10000,
                    end: 10500,
                    strand: Strand::Forward,
                    phase: -1,
                    end_phase: -1,
                    rank: 1,
                },
                Exon {
                    stable_id: "E2".into(),
                    start: 11500,
                    end: 12000,
                    strand: Strand::Forward,
                    phase: -1,
                    end_phase: -1,
                    rank: 2,
                },
            ],
            translation: None,
            cdna_coding_start: None,
            cdna_coding_end: None,
            coding_region_start: None,
            coding_region_end: None,
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: false,
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
    fn stop_retention_compares_the_full_annotated_reference_protein() {
        let mut transcript = make_coding_transcript();
        let reference = "CTGAGAGGCGCTGCTGGGCGTCTGGGCGGAGGACTCCTGGTTCTGG";
        transcript.translateable_seq = Some(reference.into());
        transcript.peptide = Some("LRGAAGRLGGGLLVL".into());
        // Transcript::translate normalizes a CTG initiator even with cds_start_NF.
        transcript.reference_peptide = Some("MRGAAGRLGGGLLVL".into());
        transcript.flags = vec!["cds_start_NF".into(), "cds_end_NF".into()];
        transcript.cdna_coding_end = Some(96);
        transcript.coding_region_end = Some(1095);
        let translation = transcript.translation.as_mut().unwrap();
        translation.genomic_end = 1095;
        translation.end_exon_rank = 1;
        translation.end_exon_offset = 95;
        let result = ConsequencePredictor::default().predict_with_parsed_input(
            &GenomicPosition::new("chr1", 1050, 1095, Strand::Forward),
            &Allele::from_str(reference), &[Allele::from_str(&format!("{reference}C"))],
            &[&transcript], None, true,
        );
        let allele = &result.transcript_consequences[0].allele_consequences[0];
        assert_eq!(allele.consequences, vec![Consequence::FrameshiftVariant]);
        assert!(allele.frameshift);
        assert_eq!(allele.amino_acids, Some(("LRGAAGRLGGGLLVLX".into(), "LRGAAGRLGGGLLVLX".into())));
    }

    #[test]
    fn complete_transcript_tandem_gain_uses_unminimized_input_span() {
        let mut transcript = make_noncoding_transcript();
        transcript.end = transcript.start + 3;
        transcript.exons.truncate(1);
        transcript.exons[0].end = transcript.end;
        let predictor = ConsequencePredictor::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            transcript.strand = strand;
            transcript.exons[0].strand = strand;
            for (alternate, parsed, expected) in [
                ("ACGTACGT", true, true),
                ("ACGTACGTACGT", true, true),
                ("ACGTACG", true, false),
                ("ACGTACGA", true, false),
                ("ACGTACGT", false, false),
            ] {
                let result = predictor.predict_with_parsed_input(
                    &GenomicPosition::new("chr1", transcript.start, transcript.end, Strand::Forward),
                    &Allele::from_str("ACGT"), &[Allele::from_str(alternate)],
                    &[&transcript], None, parsed,
                );
                let allele = &result.transcript_consequences[0].allele_consequences[0];
                assert_eq!(allele.consequences.contains(&Consequence::TranscriptAmplification), expected,
                    "{strand:?} {alternate} parsed={parsed}");
                if expected {
                    assert_eq!(allele.consequences, vec![Consequence::TranscriptAmplification]);
                    assert_eq!(allele.impact, Impact::High);
                }
            }
        }
    }

    #[test]
    fn mature_mirna_ranges_replace_the_broader_noncoding_exon_term() {
        let mut transcript = make_noncoding_transcript();
        transcript.stable_id = "ENST_MIRNA".into();
        transcript.biotype = "miRNA".into();
        let predictor = ConsequencePredictor::default()
            .with_mature_mirna_ranges(HashMap::from([("ENST_MIRNA".to_owned(), vec![(490, 520)])]));
        let consequences_at = |position| {
            predictor
                .predict(
                    &GenomicPosition::new("chr1", position, position, Strand::Forward),
                    &Allele::from_str("A"),
                    &[Allele::from_str("G")],
                    &[&transcript],
                    None,
                )
                .transcript_consequences[0]
                .allele_consequences[0]
                .consequences
                .clone()
        };

        // The cDNA interval crosses the intron: 490-501 maps to exon 1 and
        // 502-520 maps to exon 2. Both exonic segments are mature miRNA.
        for position in [10495, 11510] {
            let consequences = consequences_at(position);
            assert!(consequences.contains(&Consequence::MatureMirnaVariant));
            assert!(!consequences.contains(&Consequence::NonCodingTranscriptExonVariant));
        }
        let consequences = consequences_at(11530);
        assert!(!consequences.contains(&Consequence::MatureMirnaVariant));
        assert!(consequences.contains(&Consequence::NonCodingTranscriptExonVariant));
    }

    #[test]
    fn mature_mirna_ranges_follow_reverse_strand_cdna_order() {
        let mut transcript = make_noncoding_transcript();
        transcript.stable_id = "ENST_MIRNA_REVERSE".into();
        transcript.biotype = "miRNA".into();
        transcript.strand = Strand::Reverse;
        transcript.gene.strand = Strand::Reverse;
        for exon in &mut transcript.exons {
            exon.strand = Strand::Reverse;
        }
        let predictor =
            ConsequencePredictor::default().with_mature_mirna_ranges(HashMap::from([(
                "ENST_MIRNA_REVERSE".to_owned(),
                vec![(10, 20)],
            )]));
        let result = predictor.predict(
            &GenomicPosition::new("chr1", 11985, 11985, Strand::Forward),
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&transcript],
            None,
        );
        let consequences = &result.transcript_consequences[0].allele_consequences[0].consequences;
        assert!(consequences.contains(&Consequence::MatureMirnaVariant));
        assert!(!consequences.contains(&Consequence::NonCodingTranscriptExonVariant));
    }

    /// Nine predicates run over the codon window and every one that holds is a
    /// term Ensembl would report. Keeping two dropped the third silently.
    ///
    /// A change that replaces the initiator, introduces a stop and shifts the
    /// frame holds `stop_gained`, `frameshift_variant` and `start_lost` at once.
    /// `start_lost` ranks below the other two, so it was the one lost.
    ///
    /// PVS1 happens to be unaffected here - `NonsenseOrFrameshift` outranks
    /// `StartLost` when both are present, so the criterion takes the same branch
    /// either way - but the reported consequence set was still short a term
    /// Ensembl's model holds.
    #[test]
    fn every_coding_term_that_holds_is_reported() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // CDS 1-3 is the initiator.
        let pos = GenomicPosition::new("chr1", 1050, 1052, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("ATG"),
            &[Allele::from_str("TAAA")],
            &[&tr],
            None,
        );
        let got = &result.transcript_consequences[0].allele_consequences[0].consequences;
        for term in [
            Consequence::StopGained,
            Consequence::FrameshiftVariant,
            Consequence::StartLost,
        ] {
            assert!(got.contains(&term), "{term:?} missing from {got:?}");
        }
    }

    #[test]
    fn test_upstream_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        let pos = GenomicPosition::new("chr1", 500, 500, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        assert_eq!(result.transcript_consequences.len(), 1);
        let tc = &result.transcript_consequences[0];
        assert_eq!(tc.allele_consequences.len(), 1);
        assert!(tc.allele_consequences[0]
            .consequences
            .contains(&Consequence::UpstreamGeneVariant));
        assert_eq!(tc.allele_consequences[0].distance, Some(500));
    }

    #[test]
    fn test_downstream_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        let pos = GenomicPosition::new("chr1", 5500, 5500, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac
            .consequences
            .contains(&Consequence::DownstreamGeneVariant));
        assert_eq!(ac.distance, Some(500));
    }

    #[test]
    fn padded_downstream_indel_uses_vep_left_first_distance() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // This is the post-VCF-parser form of a padded insertion. VEP removes
        // the remaining three-base prefix, leaving insertion coordinates
        // 5504-5503. The closest endpoint is therefore 503 bases from the
        // transcript end, not 501 bases from the still-padded start.
        let pos = GenomicPosition::new("chr1", 5501, 5503, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("AAA"),
            &[Allele::from_str("AAAA")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac
            .consequences
            .contains(&Consequence::DownstreamGeneVariant));
        assert_eq!(ac.distance, Some(503));
    }

    #[test]
    fn padded_reverse_insertion_uses_vep_coordinates_for_direction() {
        let mut tr = make_coding_transcript();
        tr.strand = Strand::Reverse;
        let result = ConsequencePredictor::default().predict(
            &GenomicPosition::new("chr1", 4999, 5002, Strand::Forward),
            &Allele::from_str("ACAC"), &[Allele::from_str("ACACAC")], &[&tr], None,
        );
        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert_eq!(ac.consequences, vec![Consequence::UpstreamGeneVariant]);
        assert_eq!(ac.distance, Some(2));
    }

    #[test]
    fn padded_intronic_deletion_does_not_reach_the_exon() {
        for strand in [Strand::Forward, Strand::Reverse] {
            for (mut tr, start) in [(make_coding_transcript(), 1200), (make_noncoding_transcript(), 10500)] {
                tr.strand = strand;
                let result = ConsequencePredictor::default().predict(
                    &GenomicPosition::new("chr1", start, start + 11, Strand::Forward),
                    &Allele::from_str("ATTTATTTATTT"), &[Allele::from_str("ATTTATTT")], &[&tr], None,
                );
                let ac = &result.transcript_consequences[0].allele_consequences[0];
                assert_eq!(ac.exon, None, "{strand:?}: {:?}", ac.consequences);
                assert!(!ac.consequences.contains(&Consequence::NonCodingTranscriptExonVariant));
                assert!(!ac.consequences.contains(&Consequence::CodingSequenceVariant));
                assert!(ac.amino_acids.is_none());
            }
        }
    }

    #[test]
    fn padded_indel_crossing_distance_limit_is_intergenic() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // The padded interval begins within 5 kb, but VEP's left-first input
        // coordinates are 10002-10001, whose nearest endpoint is 5001 bases
        // away. It must not be called a downstream-gene variant.
        let pos = GenomicPosition::new("chr1", 9999, 10001, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("AAA"),
            &[Allele::from_str("AAAA")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert_eq!(ac.consequences, vec![Consequence::IntergenicVariant]);
        assert_eq!(ac.distance, Some(5001));
    }

    #[test]
    fn vep_input_position_trims_suffix_after_prefix() {
        let pos = GenomicPosition::new("chr1", 500, 502, Strand::Forward);
        let (got, reference, alternate) = vep_input_variant(
            &pos,
            &Allele::from_str("ACA"),
            &Allele::from_str("AA"),
            false,
        );
        assert_eq!((got.start, got.end), (501, 501));
        assert_eq!(reference, Allele::from_str("C"));
        assert_eq!(alternate, Allele::Deletion);
    }

    #[test]
    fn vep_input_position_preserves_equal_length_and_multi_alt_intervals() {
        let pos = GenomicPosition::new("chr1", 5500, 5501, Strand::Forward);
        let reference = Allele::from_str("AC");
        let alternate = Allele::from_str("AT");

        let biallelic = vep_input_position(&pos, &reference, &alternate, false);
        assert_eq!((biallelic.start, biallelic.end), (5500, 5501));

        let multi_alt =
            vep_input_position(&pos, &Allele::from_str("C"), &Allele::from_str("CAA"), true);
        assert_eq!((multi_alt.start, multi_alt.end), (5500, 5501));
    }

    #[test]
    fn test_intergenic_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Very far away
        let pos = GenomicPosition::new("chr1", 100000, 100000, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::IntergenicVariant));
    }

    #[test]
    fn test_5_prime_utr_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1020 is in exon1 (1000-1200), before CDS start (1050)
        let pos = GenomicPosition::new("chr1", 1020, 1020, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::FivePrimeUtrVariant));
    }

    #[test]
    fn test_3_prime_utr_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 4600 is in exon3 (4000-5000), after CDS end (4500)
        let pos = GenomicPosition::new("chr1", 4600, 4600, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::ThreePrimeUtrVariant));
    }

    #[test]
    fn a_span_straddling_a_coding_transcript_edge_matches_vep_utr_overlap() {
        let predictor = ConsequencePredictor::default();
        let mut tr = make_boundary_transcript(Strand::Forward);

        tr.coding_region_start = Some(tr.start);
        assert!(predictor.reaches_low_side_utr(tr.start - 1, tr.start, &tr));

        tr.coding_region_end = Some(tr.end);
        assert!(predictor.reaches_high_side_utr(tr.end, tr.end + 1, &tr));
    }

    #[test]
    fn test_intron_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1500 is in intron1 (1201-1999), away from splice sites
        let pos = GenomicPosition::new("chr1", 1500, 1500, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_splice_donor_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1201 is first base of intron1 → splice donor
        let pos = GenomicPosition::new("chr1", 1201, 1201, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("A")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::SpliceDonorVariant));
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn test_splice_acceptor_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1999 is last base of intron1 → splice acceptor
        let pos = GenomicPosition::new("chr1", 1999, 1999, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("A")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac
            .consequences
            .contains(&Consequence::SpliceAcceptorVariant));
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn resolved_repeat_without_a_coding_predicate_defers_to_the_transcript_default() {
        let predictor = ConsequencePredictor::default();
        let mut tr = make_coding_transcript();
        tr.translateable_seq.as_mut().unwrap().replace_range(3..12, "TATGAATAA");
        for (biotype, expected) in [
            ("protein_coding", Consequence::IntergenicVariant),
            ("nonsense_mediated_decay", Consequence::NmdTranscriptVariant),
        ] {
            tr.biotype = biotype.into();
            let result = predictor.predict(
                &GenomicPosition::new("chr1", 1053, 1059, Strand::Forward),
                &Allele::from_str("TATGAAT"),
                &[Allele::from_str("ATACTTA"), Allele::from_str(&"TATGAAT".repeat(4))],
                &[&tr], None,
            );
            let allele = &result.transcript_consequences[0].allele_consequences[1];
            assert_eq!(allele.consequences, vec![expected]);
            assert_eq!(allele.amino_acids, Some(("YE*".into(), "YEL*IMNYE*".into())));
        }
    }

    #[test]
    fn test_synonymous_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // First CDS position is at genomic 1050, which is cDNA pos 51, CDS pos 1
        // translateable_seq starts with "ATG" (Met)
        // CDS pos 3 (third base of first codon) - change G to A: ATA still codes for... wait
        // ATG -> Met. Let's change position 3 of ATG from G to something that's still Met: not possible
        // Let's use a different codon. CDS pos 4-6 is "GCT" (Ala). GCC also codes for Ala.
        // Genomic pos for CDS pos 4 = 1050 + 3 = 1053
        // Change T at CDS pos 6 to C: GCT -> GCC both = Ala
        let pos = GenomicPosition::new("chr1", 1055, 1055, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("T"),
            &[Allele::from_str("C")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::SynonymousVariant),
            "Expected synonymous, got: {:?}",
            ac.consequences
        );
        assert_eq!(ac.impact, Impact::Low);
    }

    #[test]
    fn ambiguous_coding_boundary_does_not_assert_loss_or_retention() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        for start in [tr.cdna_coding_start.unwrap(), tr.cdna_coding_end.unwrap()-2] {
            assert!(predictor.predict_cds_boundary_consequence(
                &Allele::from_str("ATG"), &Allele::from_str("N"), &tr,
                Some(start), Some(start+2)).is_none());
        }
    }

    #[test]
    fn ambiguous_insertion_can_preserve_a_known_initiator() {
        let tr = make_coding_transcript();
        let result = ConsequencePredictor::default().predict(
            &GenomicPosition::new("chr1", 1051, 1050, Strand::Forward),
            &Allele::Deletion, &[Allele::from_str("NA"), Allele::from_str("NN")],
            &[&tr], None,
        );
        let alleles = &result.transcript_consequences[0].allele_consequences;
        // A[NA]TG retains the original ATG suffix; A[NN]TG cannot prove it.
        assert!(alleles[0].consequences.contains(&Consequence::StartRetainedVariant));
        assert!(!alleles[1].consequences.contains(&Consequence::StartRetainedVariant));
        for allele in alleles {
            assert!(allele.consequences.contains(&Consequence::FrameshiftVariant));
            assert!(!allele.consequences.contains(&Consequence::StartLost));
            assert!(allele.amino_acids.is_none());
        }
    }

    #[test]
    fn reference_repeated_as_an_alternate_has_no_annotation() {
        let tr = make_coding_transcript();
        let result = ConsequencePredictor::default().predict(
            &GenomicPosition::new("chr1", 1053, 1053, Strand::Forward),
            &Allele::from_str("G"), &[Allele::from_str("G"), Allele::from_str("N")],
            &[&tr], None);
        let alleles = &result.transcript_consequences[0].allele_consequences;
        assert_eq!(alleles.len(), 1);
        assert_eq!(alleles[0].allele, Allele::from_str("N"));
        assert_eq!(alleles[0].consequences, vec![Consequence::CodingSequenceVariant]);
    }

    #[test]
    fn parsed_input_does_not_repeat_minimization_after_case_conversion() {
        let tr = make_coding_transcript();
        // VEP parses raw a>AT without trimming, then validates it as A>AT.
        let result = ConsequencePredictor::default().predict_with_parsed_input(
            &GenomicPosition::new("chr1", 1050, 1050, Strand::Forward),
            &Allele::from_str("A"), &[Allele::from_str("AT")], &[&tr], None, true,
        );
        let allele = &result.transcript_consequences[0].allele_consequences[0];
        assert_eq!((allele.cdna_start, allele.cdna_end), (Some(51), Some(51)));
        assert_eq!(allele.allele, Allele::from_str("AT"));
    }

    #[test]
    fn ambiguous_alternate_has_no_peptide_but_retains_dna_consequences() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        let pos = GenomicPosition::new("chr1", 1053, 1053, Strand::Forward);
        for (alternate, expected) in [("N", Consequence::CodingSequenceVariant),
            ("NN", Consequence::FrameshiftVariant), ("NNNN", Consequence::CodingSequenceVariant)] {
            let result = predictor.predict(&pos, &Allele::from_str("G"),
                &[Allele::from_str(alternate)], &[&tr], None);
            let allele = &result.transcript_consequences[0].allele_consequences[0];
            assert_eq!(allele.consequences, vec![expected], "{alternate}");
            assert!(allele.amino_acids.is_none(), "{alternate}");
            assert!(allele.codons.is_some(), "DNA codons remain available");
        }
    }

    #[test]
    fn test_missense_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // CDS pos 4-6 is "GCT" (Ala). Change first base G to T: TCT = Ser (different!)
        // Genomic pos for CDS pos 4 = 1050 + 3 = 1053
        let pos = GenomicPosition::new("chr1", 1053, 1053, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("T")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::MissenseVariant),
            "Expected missense, got: {:?}",
            ac.consequences
        );
        assert_eq!(ac.impact, Impact::Moderate);
    }

    #[test]
    fn source_edited_reference_residue_does_not_rewrite_changed_alt_residue() {
        let predictor = ConsequencePredictor::default();
        let mut tr = make_coding_transcript();
        tr.translateable_seq
            .as_mut()
            .unwrap()
            .replace_range(..3, "CTG");
        tr.peptide = Some("M".to_string());

        let result = predictor.predict(
            &GenomicPosition::new("chr1", 1051, 1051, Strand::Forward),
            &Allele::from_str("T"),
            &[Allele::from_str("A")],
            &[&tr],
            None,
        );
        let ac = &result.transcript_consequences[0].allele_consequences[0];

        assert!(ac.consequences.contains(&Consequence::StartLost));
        assert_eq!(ac.amino_acids, Some(("M".into(), "Q".into())));
        assert_eq!(ac.codons, Some(("cTg".into(), "cAg".into())));
    }

    #[test]
    fn annotated_residue_resolution_is_limited_to_known_source_edits() {
        let mut tr = make_coding_transcript();
        tr.peptide = Some("MWG".to_string());

        assert_eq!(resolve_annotated_residue(&tr, 0, b'L'), b'M');
        assert_eq!(resolve_annotated_residue(&tr, 1, b'*'), b'W');
        assert_eq!(resolve_annotated_residue(&tr, 2, b'A'), b'A');
    }

    #[test]
    fn unchanged_stop_codon_in_an_alt_window_does_not_inherit_reference_edits() {
        let predictor = ConsequencePredictor::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            for edited in ['U', 'X', 'W'] {
                let mut tr = make_boundary_transcript(strand);
                tr.translateable_seq.as_mut().unwrap().replace_range(3..6, "TGA");
                tr.spliced_seq.as_mut().unwrap().replace_range(13..16, "TGA");
                tr.peptide = Some(format!("M{edited}AAAAAAA*"));
                let (start, end, alternate) = if strand == Strand::Forward {
                    (6, 5, "A")
                } else {
                    (5, 6, "T")
                };
                let change = predictor.predict_coding_consequence(
                    &Allele::Deletion, &Allele::from_str(alternate), &tr,
                    Some(start), Some(end), Some(start + 10), Some(end + 10),
                ).unwrap();
                assert_eq!(change.amino_acids, Some((edited.to_string(), "*".into())));
                assert!(change.consequence == Consequence::StopGained
                    || change.additional.contains(&Consequence::StopGained));
            }
        }
    }

    #[test]
    fn cds_start_nf_suppresses_a_start_lost_claim() {
        let predictor = ConsequencePredictor::default();
        let mut tr = make_coding_transcript();
        tr.flags.push("cds_start_NF".into());

        let result = predictor.predict(
            &GenomicPosition::new("chr1", 1050, 1050, Strand::Forward),
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );
        let ac = &result.transcript_consequences[0].allele_consequences[0];

        assert_eq!(ac.consequences, vec![Consequence::MissenseVariant]);
        assert_eq!(ac.impact, Impact::Moderate);
        assert_eq!(ac.amino_acids, Some(("M".into(), "V".into())));
    }

    #[test]
    fn missense_and_unknown_coding_are_independent() {
        let predictor = ConsequencePredictor::default();
        let change = predictor.terms_for_window(
            &CodonWindow {
                peptides_defined: true,
                alt_dna_unambiguous: true,
                ref_aas: "FX".into(),
                alt_aas: "LX".into(),
                ref_window: b"TTCN".to_vec(),
                alt_window: b"CTCN".to_vec(),
                ref_codons: "ttCN".into(),
                alt_codons: "ttAC".into(),
                ref_len: 2,
                alt_len: 2,
                tl_start: 1,
                tl_end: 2,
                partial_codon: false,
            },
            false,
            None,
            false,
            false,
            false,
            None,
        );

        assert_eq!(change.consequence, Consequence::MissenseVariant);
        assert_eq!(change.additional, vec![Consequence::CodingSequenceVariant]);
    }

    #[test]
    fn inframe_insertion_and_unknown_coding_are_independent() {
        let predictor = ConsequencePredictor::default();
        let change = predictor.terms_for_window(
            &CodonWindow {
                peptides_defined: true,
                alt_dna_unambiguous: true,
                ref_aas: "X".into(),
                alt_aas: "XL".into(),
                ref_window: b"N".to_vec(),
                alt_window: b"NTTC".to_vec(),
                ref_codons: "-".into(),
                alt_codons: "TTC".into(),
                ref_len: 0,
                alt_len: 3,
                tl_start: 1,
                tl_end: 1,
                partial_codon: true,
            },
            false,
            None,
            false,
            false,
            false,
            None,
        );

        assert_eq!(change.consequence, Consequence::InframeInsertion);
        assert_eq!(
            change.additional,
            vec![
                Consequence::IncompleteTerminalCodonVariant,
                Consequence::CodingSequenceVariant,
            ]
        );
    }

    #[test]
    fn partial_codon_survives_one_unmapped_endpoint_on_either_strand() {
        let predictor = ConsequencePredictor::default();
        for (strand, start, end) in [
            (Strand::Forward, Some(2), None),
            (Strand::Reverse, None, Some(2)),
        ] {
            let mut tr = make_coding_transcript();
            tr.strand = strand;
            tr.translateable_seq = Some("GC".into());
            tr.cdna_coding_start = Some(1);
            tr.cdna_coding_end = Some(2);
            let change = predictor
                .predict_coding_consequence(
                    &Allele::from_str("CT"),
                    &Allele::from_str("GA"),
                    &tr,
                    start,
                    end,
                    start.or(Some(3)),
                    end.or(Some(3)),
                )
                .unwrap();
            assert_eq!(change.consequence, Consequence::CodingSequenceVariant);
            assert_eq!(
                change.additional,
                vec![Consequence::IncompleteTerminalCodonVariant]
            );
            let (gap_start, gap_end) = if strand == Strand::Forward {
                (None, Some(2))
            } else {
                (Some(2), None)
            };
            let gap = predictor.predict_coding_consequence(
                &Allele::from_str("CT"), &Allele::from_str("GA"), &tr,
                gap_start, gap_end, gap_start, gap_end,
            );
            assert!(gap.is_none());
        }
    }

    #[test]
    fn alternate_partial_codon_can_read_into_the_utr() {
        let predictor = ConsequencePredictor::default();
        let mut tr = make_boundary_transcript(Strand::Forward);
        tr.translateable_seq = Some("ATGAA".into());
        tr.cdna_coding_end = Some(15);
        tr.spliced_seq = Some("CCCCCCCCCCATGAAA".into());
        let window = predictor
            .build_codon_window(
                &tr,
                &Allele::from_str("A"),
                &Allele::from_str("C"),
                Some(4),
                Some(4),
            )
            .unwrap();
        assert_eq!(
            (window.ref_aas.as_str(), window.alt_aas.as_str()),
            ("X", "Q")
        );
        assert_eq!(
            (window.ref_codons.as_str(), window.alt_codons.as_str()),
            ("Aa", "Caa")
        );
        assert!(window.partial_codon);
        let snv = predictor
            .predict_coding_consequence(
                &Allele::from_str("A"),
                &Allele::from_str("C"),
                &tr,
                Some(4),
                Some(4),
                Some(14),
                Some(14),
            )
            .unwrap();
        assert_eq!(snv.consequence, Consequence::IncompleteTerminalCodonVariant);
        assert!(!snv.additional.contains(&Consequence::MissenseVariant));
        tr.translateable_seq = Some("ATGA".into());
        tr.cdna_coding_end = Some(14);
        tr.spliced_seq = Some("CCCCCCCCCCATGATG".into());
        let deletion = predictor.predict_coding_consequence(
            &Allele::from_str("A"), &Allele::Deletion, &tr,
            Some(4), Some(4), Some(14), Some(14),
        ).unwrap();
        assert_eq!(deletion.consequence, Consequence::IncompleteTerminalCodonVariant);
        assert!(!deletion.additional.contains(&Consequence::InframeInsertion));
        tr.translateable_seq = Some("ATGAA".into());
        tr.cdna_coding_end = Some(15);
        tr.spliced_seq = Some("CCCCCCCCCCATGAAA".into());
        let change = predictor
            .predict_coding_consequence(
                &Allele::from_str("AA"),
                &Allele::from_str("CC"),
                &tr,
                Some(5),
                None,
                Some(15),
                Some(16),
            )
            .unwrap();
        assert!(change
            .additional
            .contains(&Consequence::IncompleteTerminalCodonVariant));
    }

    #[test]
    fn test_stop_gained() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // CDS pos 4-6 is "GCT". Change to "TAA" (stop) → need to change pos 4,5,6
        // For simplicity, change CDS pos 4: G->T, pos 5: C->A, pos 6: T->A
        // But our predictor works one SNV at a time. Let's pick a codon that's one base
        // away from a stop. "TCA" (Ser) → change C→A: TAA (stop). But we'd need that codon.
        // Actually, translateable_seq[6..9] = "TCA" (positions 7-9 in 1-based)
        // CDS pos 7 is at genomic 1050+6 = 1056
        // Change T to T (no), we need C at pos 8 to become something.
        // Let's just use translateable[3..6] = "GCT" and change pos 4 (G) to T: "TCT" = Ser
        // That's missense, not stop. Let's try another approach.
        // translateable[9..12] = "AAG" (Lys). Change A at pos 10 to T: TAG = stop!
        // CDS pos 10 is at genomic 1050+9 = 1059
        let pos = GenomicPosition::new("chr1", 1059, 1059, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("T")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::StopGained),
            "Expected stop_gained, got: {:?}. translateable[9..12]={:?}",
            ac.consequences,
            &tr.translateable_seq.as_ref().unwrap()[9..12]
        );
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn test_frameshift_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Deletion of 1 base in CDS → frameshift
        // CDS pos 4 at genomic 1053
        let pos = GenomicPosition::new("chr1", 1053, 1053, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::Deletion],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::FrameshiftVariant),
            "Expected frameshift, got: {:?}",
            ac.consequences
        );
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn test_inframe_deletion() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Deletion of 3 bases in CDS → inframe deletion
        let pos = GenomicPosition::new("chr1", 1053, 1055, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("GCT"),
            &[Allele::Deletion],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::InframeDeletion),
            "Expected inframe_deletion, got: {:?}",
            ac.consequences
        );
        assert_eq!(ac.impact, Impact::Moderate);
    }

    #[test]
    fn test_noncoding_exon_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_noncoding_transcript();
        let pos = GenomicPosition::new("chr1", 10100, 10100, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac
            .consequences
            .contains(&Consequence::NonCodingTranscriptExonVariant));
    }

    #[test]
    fn reported_cdna_positions_use_vep_left_first_input_coordinates() {
        let predictor = ConsequencePredictor::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            let mut tr = make_noncoding_transcript();
            tr.strand = strand;
            tr.gene.strand = strand;
            for exon in &mut tr.exons {
                exon.strand = strand;
            }

            // Post-parser form of GCA/GCACA. VEP trims the remaining CA
            // prefix for position reporting, producing 10102-10101.
            let pos = GenomicPosition::new("chr1", 10100, 10101, Strand::Forward);
            let result = predictor.predict(
                &pos,
                &Allele::from_str("CA"),
                &[Allele::from_str("CACA")],
                &[&tr],
                None,
            );
            let ac = &result.transcript_consequences[0].allele_consequences[0];

            assert_eq!(ac.cdna_start, tr.genomic_to_cdna(10102));
            assert_eq!(ac.cdna_end, tr.genomic_to_cdna(10101));
        }
    }

    #[test]
    fn an_insertion_between_a_noncoding_exon_and_intron_is_not_in_either() {
        let predictor = ConsequencePredictor::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            let mut tr = make_noncoding_transcript();
            tr.strand = strand;
            tr.gene.strand = strand;
            for exon in &mut tr.exons {
                exon.strand = strand;
            }

            // A VCF insertion after genomic base 10500 becomes Ensembl's
            // zero-length interval between the exon and intron: 10501-10500.
            let pos = GenomicPosition::new("chr1", 10501, 10500, Strand::Forward);
            let result = predictor.predict(
                &pos,
                &Allele::Deletion,
                &[Allele::from_str("T")],
                &[&tr],
                None,
            );
            let ac = &result.transcript_consequences[0].allele_consequences[0];

            assert!(ac.consequences.contains(&Consequence::SpliceRegionVariant));
            assert!(ac
                .consequences
                .contains(&Consequence::NonCodingTranscriptVariant));
            assert!(!ac
                .consequences
                .contains(&Consequence::NonCodingTranscriptExonVariant));
            assert!(!ac.consequences.contains(&Consequence::IntronVariant));
            assert_eq!(ac.exon, None);
            assert_eq!(ac.intron, None);
        }
    }

    #[test]
    fn left_minimal_splice_predicates_remove_a_false_region_term() {
        let predictor = ConsequencePredictor::default();
        let mut tr = make_noncoding_transcript();
        tr.strand = Strand::Reverse;
        tr.gene.strand = Strand::Reverse;
        for exon in &mut tr.exons {
            exon.strand = Strand::Reverse;
        }
        // Intron 10501..11499. This post-parser A/AA insertion begins at +7
        // from the genomic start, where the padded representation appears to
        // be in splice_region on a reverse-strand transcript.
        // VEP removes the remaining shared prefix first, placing the actual
        // insertion at 10509..10508: still intronic and in the
        // polypyrimidine tract, but no longer in splice_region.
        let pos = GenomicPosition::new("chr1", 10508, 10508, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("AA")],
            &[&tr],
            None,
        );
        let ac = &result.transcript_consequences[0].allele_consequences[0];

        assert!(ac
            .consequences
            .contains(&Consequence::SplicePolypyrimidineTractVariant));
        assert!(ac.consequences.contains(&Consequence::IntronVariant));
        assert!(!ac.consequences.contains(&Consequence::SpliceRegionVariant));
    }

    #[test]
    fn a_frameshift_intron_is_preclassified_as_coding_without_an_exon_label() {
        let predictor = ConsequencePredictor::default();
        let mut tr = make_coding_transcript();
        tr.exons[1].start = 1202; // one-base intron at 1201

        let pos = GenomicPosition::new("chr1", 1201, 1201, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("T"),
            &[Allele::Deletion],
            &[&tr],
            None,
        );
        let ac = &result.transcript_consequences[0].allele_consequences[0];

        assert_eq!(ac.consequences, vec![Consequence::CodingSequenceVariant]);
        assert_eq!(ac.exon, None);
        assert_eq!(ac.intron, Some((1, 1, 2)));

        // The short intron also stretches the last exon for preclassification,
        // suppressing a polypyrimidine candidate at -10. It must not make this
        // ordinary long-intron position coding.
        let ordinary_intron = GenomicPosition::new("chr1", 3990, 3990, Strand::Forward);
        let result = predictor.predict(
            &ordinary_intron,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );
        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert_eq!(ac.consequences, vec![Consequence::IntronVariant]);

        tr.translation = None;
        tr.biotype = "unprocessed_pseudogene".into();
        let result = predictor.predict(
            &pos,
            &Allele::from_str("T"),
            &[Allele::from_str("A")],
            &[&tr],
            None,
        );
        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert_eq!(
            ac.consequences,
            vec![Consequence::NonCodingTranscriptVariant]
        );
        assert_eq!(ac.exon, None);
    }

    #[test]
    fn test_start_lost() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // CDS pos 1 (first base of ATG) is at genomic 1050
        // Change A to G: GTG is not a standard start codon
        let pos = GenomicPosition::new("chr1", 1050, 1050, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::StartLost),
            "Expected start_lost, got: {:?}",
            ac.consequences
        );
    }

    #[test]
    fn an_empty_predicate_result_uses_veps_default_after_nmd() {
        let predictor = ConsequencePredictor::default();
        let mut tr = make_coding_transcript();
        tr.exons[1].start = 1202;
        // Genomic MNV, but unequal CDS lengths across the one-base intron.
        // VEP's genomic-shape filters reject the resulting frameshift term.
        let position = GenomicPosition::new("chr1", 1198, 1206, Strand::Forward);
        for (biotype, expected) in [
            ("protein_coding", Consequence::IntergenicVariant),
            ("nonsense_mediated_decay", Consequence::NmdTranscriptVariant),
        ] {
            tr.biotype = biotype.into();
            let result = predictor.predict(&position, &Allele::from_str("AAAAAAAAA"),
                &[Allele::from_str("TTTTTTTTT")], &[&tr], None);
            assert_eq!(result.transcript_consequences[0].allele_consequences[0].consequences,
                vec![expected]);
        }
    }

    #[test]
    fn insertion_retaining_atg_needs_a_utr_for_the_start_shortcut() {
        let predictor = ConsequencePredictor::default();
        for has_utr in [false, true] {
            let mut tr = make_boundary_transcript(Strand::Forward);
            if !has_utr {
                tr.spliced_seq = tr.translateable_seq.clone();
                tr.cdna_coding_start = Some(1);
                tr.cdna_coding_end = Some(tr.translateable_seq.as_ref().unwrap().len() as u64);
            }
            let coding_start = tr.cdna_coding_start.unwrap();
            let change = predictor.predict_coding_consequence(
                &Allele::Deletion, &Allele::from_str("GGGCC"), &tr,
                Some(3), Some(2), Some(coding_start + 2), Some(coding_start + 1),
            ).unwrap();
            assert_eq!(change.additional.contains(&Consequence::StartLost), !has_utr);
            assert_eq!(change.consequence, Consequence::FrameshiftVariant);
        }
    }

    #[test]
    fn non_atg_start_predicates_are_independent_of_peptide_identity() {
        let predictor = ConsequencePredictor::default();
        let mut tr = make_boundary_transcript(Strand::Forward);
        tr.translateable_seq.as_mut().unwrap().replace_range(..3, "CTG");
        tr.spliced_seq.as_mut().unwrap().replace_range(10..13, "CTG");
        tr.peptide = Some("L".into());
        let synonymous = predictor.predict_coding_consequence(
            &Allele::from_str("G"), &Allele::from_str("C"), &tr,
            Some(3), Some(3), Some(13), Some(13),
        ).unwrap();
        assert_eq!(synonymous.consequence, Consequence::StartLost);
        assert!(synonymous.additional.contains(&Consequence::SynonymousVariant));
        let retained = predictor.predict_coding_consequence(
            &Allele::from_str("C"), &Allele::from_str("A"), &tr,
            Some(1), Some(1), Some(11), Some(11),
        ).unwrap();
        assert_eq!(retained.consequence, Consequence::StartLost);
        assert!(retained.additional.contains(&Consequence::StartRetainedVariant));
        let mnv = predictor.predict_coding_consequence(
            &Allele::from_str("CT"), &Allele::from_str("CA"), &tr,
            Some(1), Some(2), Some(11), Some(12),
        ).unwrap();
        assert!(!mnv.additional.contains(&Consequence::StartRetainedVariant));
        tr.peptide = Some("M".into());
        let expanded_start = predictor.predict_coding_consequence(
            &Allele::Deletion, &Allele::from_str("G"), &tr,
            Some(3), Some(2), Some(13), Some(12),
        ).unwrap();
        assert_eq!(expanded_start.amino_acids, Some(("M".into(), "LX".into())));
        tr.translateable_seq.as_mut().unwrap().replace_range(..3, "GTT");
        tr.spliced_seq.as_mut().unwrap().replace_range(10..13, "GTT");
        tr.peptide = Some("V".into());
        let insertion = predictor.predict_coding_consequence(
            &Allele::Deletion, &Allele::from_str("T"), &tr,
            Some(2), Some(1), Some(12), Some(11),
        ).unwrap();
        assert!(!insertion.additional.contains(&Consequence::StartRetainedVariant));
    }

    #[test]
    fn an_annotated_mitochondrial_initiator_is_the_reference_residue() {
        let mut tr = make_boundary_transcript(Strand::Forward);
        tr.chromosome = "chrM".into();
        tr.gene.chromosome = "chrM".into();
        tr.spliced_seq
            .as_mut()
            .unwrap()
            .replace_range(10..13, "ATT");
        tr.translateable_seq
            .as_mut()
            .unwrap()
            .replace_range(0..3, "ATT");
        // Human mitochondrial CDSs have no 5' UTR. Keep this fixture's
        // artificial flank out of VEP's separate fixed-offset UTR predicate.
        tr.spliced_seq = None;

        // Vertebrate mitochondrial ATT initiates as Met, and ATA is also Met.
        let result = ConsequencePredictor::default().predict(
            &GenomicPosition::new("chrM", 1012, 1012, Strand::Forward),
            &Allele::from_str("T"),
            &[Allele::from_str("A"), Allele::from_str("C")],
            &[&tr],
            None,
        );
        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::SynonymousVariant));
        assert!(!ac.consequences.contains(&Consequence::StartLost));
        assert_eq!(ac.amino_acids, Some(("M".into(), "M".into())));

        // ATC is Ile in a translated human mitochondrial CDS, so this change
        // still loses the annotated methionine initiator.
        let lost = &result.transcript_consequences[0].allele_consequences[1];
        assert!(lost.consequences.contains(&Consequence::StartLost));
        assert_eq!(lost.amino_acids, Some(("M".into(), "I".into())));

        // VEP peptide() applies the annotated initial-Met edit only to the
        // reference allele. An insertion retaining ATT at the start of its
        // longer codon window must translate that alternate ATT as Ile.
        let inserted = ConsequencePredictor::default().predict(
            &GenomicPosition::new("chrM", 1011, 1010, Strand::Forward),
            &Allele::Deletion,
            &[Allele::from_str("TTA")],
            &[&tr],
            None,
        );
        let ac = &inserted.transcript_consequences[0].allele_consequences[0];
        assert_eq!(ac.amino_acids, Some(("M".into(), "II".into())));
        assert!(ac.consequences.contains(&Consequence::StartLost));
    }

    #[test]
    fn test_multiple_transcripts() {
        let predictor = ConsequencePredictor::default();
        let tr1 = make_coding_transcript();
        let tr2 = make_noncoding_transcript();

        // Position in tr1's intron, not overlapping tr2
        let pos = GenomicPosition::new("chr1", 1500, 1500, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr1, &tr2],
            None,
        );

        assert_eq!(result.transcript_consequences.len(), 2);
        // tr1: intron variant
        assert!(result.transcript_consequences[0].allele_consequences[0]
            .consequences
            .contains(&Consequence::IntronVariant));
        // tr2: 8500bp away (>5000), so intergenic
        assert!(result.transcript_consequences[1].allele_consequences[0]
            .consequences
            .contains(&Consequence::IntergenicVariant));
    }

    #[test]
    fn test_most_severe_across_transcripts() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();

        // Splice donor is more severe than intron variant
        let pos = GenomicPosition::new("chr1", 1201, 1201, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("A")],
            &[&tr],
            None,
        );

        assert_eq!(result.most_severe, Some(Consequence::SpliceDonorVariant));
    }

    #[test]
    fn test_multi_allelic() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();

        // Two alt alleles at a coding position
        let pos = GenomicPosition::new("chr1", 1053, 1053, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("T"), Allele::from_str("C")],
            &[&tr],
            None,
        );

        let tc = &result.transcript_consequences[0];
        assert_eq!(tc.allele_consequences.len(), 2);
    }

    // ---- variants spanning either end of the CDS (#100) ----

    /// A single-exon transcript carrying its own sequence, so the terminator
    /// and initiator tests can be decided rather than guessed.
    ///
    /// cDNA position `p` sits at genomic `999 + p` on the forward strand and
    /// `1201 - p` on the reverse one. The layout either way:
    ///
    ///   cDNA   1..10    5' UTR
    ///   cDNA  11..40    CDS - `ATG`, eight `GCT`, `TAA`
    ///   cDNA  41..201   3' UTR - `CC`, then `TGA`, then filler
    ///
    /// The `TGA` three bases into the 3' UTR is deliberate: it lands on the old
    /// terminator's offset after exactly five cDNA bases are removed, which is
    /// what separates `stop_lost` from `stop_retained_variant` below.
    #[test]
    fn ambiguous_reference_keeps_dna_stop_alteration_fallback() {
        for (strand, position) in [(Strand::Forward, 1037), (Strand::Reverse, 1163)] {
            let tr = make_boundary_transcript(strand);
            let result = ConsequencePredictor::default().predict(
                &GenomicPosition::new("chr1", position, position, Strand::Forward),
                &Allele::from_str("N"), &[Allele::Deletion], &[&tr], None,
            );
            let allele = &result.transcript_consequences[0].allele_consequences[0];
            assert!(allele.consequences.contains(&Consequence::StopLost), "{:?}", allele.consequences);
            assert!(allele.consequences.contains(&Consequence::FrameshiftVariant));
            assert!(allele.amino_acids.is_none());
        }
    }

    fn make_boundary_transcript(strand: Strand) -> Transcript {
        let five_utr = "C".repeat(10);
        let cds = format!("ATG{}TAA", "GCT".repeat(8));
        let three_utr = format!("CCTGA{}", "G".repeat(156));
        let spliced = format!("{five_utr}{cds}{three_utr}");
        assert_eq!(spliced.len(), 201);

        let (coding_region_start, coding_region_end) = match strand {
            // CDS is cDNA 11..40
            Strand::Forward => (1010, 1039),
            Strand::Reverse => (1161, 1190),
        };

        Transcript {
            stable_id: "ENST_BOUND".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_BOUND".into(),
                symbol: Some("BOUNDGENE".into()),
                symbol_source: Some("HGNC".into()),
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start: 1000,
                end: 1200,
                strand,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start: 1000,
            end: 1200,
            strand,
            exons: vec![Exon {
                stable_id: "E1".into(),
                start: 1000,
                end: 1200,
                strand,
                phase: -1,
                end_phase: -1,
                rank: 1,
            }],
            translation: Some(Translation {
                stable_id: "ENSP_BOUND".into(),
                genomic_start: coding_region_start,
                genomic_end: coding_region_end,
                start_exon_rank: 1,
                start_exon_offset: 10,
                end_exon_rank: 1,
                end_exon_offset: 40,
            }),
            cdna_coding_start: Some(11),
            cdna_coding_end: Some(40),
            coding_region_start: Some(coding_region_start),
            coding_region_end: Some(coding_region_end),
            spliced_seq: Some(spliced),
            translateable_seq: Some(cds),
            peptide: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: Some(1),
            appris: None,
            ccds: None,
            protein_id: Some("ENSP_BOUND".into()),
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

    /// Genomic span of a cDNA range on the given strand, as `(start, end)`.
    fn genomic_span(strand: Strand, cdna_lo: u64, cdna_hi: u64) -> (u64, u64) {
        match strand {
            Strand::Forward => (999 + cdna_lo, 999 + cdna_hi),
            Strand::Reverse => (1201 - cdna_hi, 1201 - cdna_lo),
        }
    }

    fn consequences_for(
        strand: Strand,
        cdna_lo: u64,
        cdna_hi: u64,
        ref_allele: &Allele,
        alt_allele: &Allele,
    ) -> Vec<Consequence> {
        let predictor = ConsequencePredictor::default();
        let tr = make_boundary_transcript(strand);
        let (g_start, g_end) = genomic_span(strand, cdna_lo, cdna_hi);
        let pos = GenomicPosition::new("chr1", g_start, g_end, Strand::Forward);
        let result = predictor.predict(
            &pos,
            ref_allele,
            std::slice::from_ref(alt_allele),
            &[&tr],
            None,
        );
        result.transcript_consequences[0].allele_consequences[0]
            .consequences
            .clone()
    }

    /// The reported case: a delins over the terminator and on into the 3' UTR
    /// is `stop_lost&3_prime_UTR_variant`, on either strand.
    ///
    /// It used to depend on the strand. The region was chosen from the
    /// variant's genomic *start*, which on the reverse strand is its last base
    /// in transcript order, so FOXL2's `c.1127_*4delinsCG` matched the 3' UTR
    /// first and never reached the coding branch at all - MODIFIER where VEP
    /// says HIGH. On the forward strand the same shape matched the coding
    /// branch first and lost the UTR term instead.
    #[test]
    fn delins_over_the_stop_codon_is_stop_lost_on_either_strand() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // cDNA 38..44: the whole terminator plus four 3' UTR bases.
            let got = consequences_for(
                strand,
                38,
                44,
                &Allele::from_str("TAACCTG"),
                &Allele::from_str("AC"),
            );
            assert!(
                got.contains(&Consequence::StopLost),
                "{strand:?}: expected stop_lost, got {got:?}"
            );
            assert!(
                got.contains(&Consequence::ThreePrimeUtrVariant),
                "{strand:?}: expected the 3'UTR term alongside it, got {got:?}"
            );
            assert_eq!(Consequence::worst_impact(&got), Some(Impact::High));
        }
    }

    #[test]
    fn boundary_predicates_use_the_full_deleted_prefix_and_partial_codon_guard() {
        let predictor = ConsequencePredictor::default();
        let mut tr = make_boundary_transcript(Strand::Forward);
        tr.spliced_seq.as_mut().unwrap().replace_range(8..10, "AT");
        let start = predictor.predict_coding_consequence(
            &Allele::from_str("AT"), &Allele::Deletion, &tr,
            Some(1), Some(2), Some(11), Some(12),
        ).unwrap();
        assert!(start.additional.contains(&Consequence::StartRetainedVariant));

        tr.cdna_coding_end = Some(41);
        tr.spliced_seq.as_mut().unwrap().replace_range(38..43, "TAACA");
        tr.translateable_seq = Some(tr.spliced_seq.as_ref().unwrap()[10..41].into());
        let terminal = predictor.predict_coding_consequence(
            &Allele::from_str("AC"), &Allele::Deletion, &tr,
            Some(31), None, Some(41), Some(42),
        ).unwrap();
        assert_eq!(terminal.consequence, Consequence::CodingSequenceVariant);
        assert_eq!(terminal.additional, vec![Consequence::IncompleteTerminalCodonVariant]);
    }

    #[test]
    fn a_same_length_change_across_the_stop_boundary_is_generic_coding() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // cDNA 40 is the last base of the stop codon and 41 is 3' UTR.
            // Alleles are supplied in genomic orientation.
            let (reference, alternate) = match strand {
                Strand::Forward => (Allele::from_str("AC"), Allele::from_str("CT")),
                Strand::Reverse => (Allele::from_str("GT"), Allele::from_str("AG")),
            };
            let got = consequences_for(strand, 40, 41, &reference, &alternate);
            assert!(
                got.contains(&Consequence::CodingSequenceVariant),
                "{strand:?}: expected coding_sequence_variant, got {got:?}"
            );
            assert!(got.contains(&Consequence::ThreePrimeUtrVariant));
            assert!(!got.contains(&Consequence::StopLost));
            assert!(!got.contains(&Consequence::StopRetainedVariant));
        }
    }

    /// A deletion reaching past the terminator is `stop_retained_variant` when
    /// the bases that move up into its place still read as a stop, and
    /// `stop_lost` when they do not. One base of difference decides it, which
    /// is why the codon is re-read from the sequence rather than inferred from
    /// the variant's length.
    #[test]
    fn a_deletion_past_the_stop_codon_distinguishes_lost_from_retained() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // Removing cDNA 38..42 slides the 3' UTR's `TGA` onto the
            // terminator's offset: still a stop.
            let retained = consequences_for(
                strand,
                38,
                42,
                &Allele::from_str("TAACC"),
                &Allele::Deletion,
            );
            assert!(
                retained.contains(&Consequence::StopRetainedVariant),
                "{strand:?}: expected stop_retained_variant, got {retained:?}"
            );
            assert!(!retained.contains(&Consequence::StopLost));

            // One base less, and `CTG` lands there instead: the stop is gone.
            let lost =
                consequences_for(strand, 38, 41, &Allele::from_str("TAAC"), &Allele::Deletion);
            assert!(
                lost.contains(&Consequence::StopLost),
                "{strand:?}: expected stop_lost, got {lost:?}"
            );
            assert!(!lost.contains(&Consequence::StopRetainedVariant));
        }
    }

    /// The mirror case at the other end of the CDS, which had the same bug with
    /// the strands swapped: a length change reaching across the initiator is
    /// `start_lost&5_prime_UTR_variant`.
    #[test]
    fn delins_over_the_start_codon_is_start_lost_on_either_strand() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // cDNA 9..13: two 5' UTR bases and the first three of the CDS.
            let got = consequences_for(
                strand,
                9,
                13,
                &Allele::from_str("CCATG"),
                &Allele::from_str("TT"),
            );
            assert!(
                got.contains(&Consequence::StartLost),
                "{strand:?}: expected start_lost, got {got:?}"
            );
            assert!(
                got.contains(&Consequence::FivePrimeUtrVariant),
                "{strand:?}: expected the 5'UTR term alongside it, got {got:?}"
            );
            assert_eq!(Consequence::worst_impact(&got), Some(Impact::High));
        }
    }

    #[test]
    fn boundary_start_retention_requires_preserved_utr_and_atg() {
        let predictor = ConsequencePredictor::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            for (utr, expected) in [
                ("CA", Consequence::StartRetainedVariant),
                ("CG", Consequence::StartLost),
            ] {
                let mut tr = make_boundary_transcript(strand);
                tr.cdna_coding_start = Some(3);
                tr.cdna_coding_end = Some(14);
                tr.translateable_seq = Some("ATGAATGCTTAA".into());
                tr.spliced_seq = Some(format!("{utr}ATGAATGCTTAA"));
                let (start, end) = if strand == Strand::Forward {
                    (2, 5)
                } else {
                    (5, 2)
                };
                let change = predictor
                    .predict_cds_boundary_consequence(
                        &Allele::from_str("AATG"),
                        &Allele::Deletion,
                        &tr,
                        Some(start),
                        Some(end),
                    )
                    .unwrap();
                assert_eq!(change.consequence, expected, "{strand:?}, {utr}");
            }
        }
    }

    #[test]
    fn depleted_cds_cannot_borrow_a_start_codon_from_three_prime_utr() {
        let predictor = ConsequencePredictor::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            let mut tr = make_boundary_transcript(strand);
            tr.spliced_seq.as_mut().unwrap().replace_range(40..42, "TG");
            let reference = &tr.translateable_seq.as_ref().unwrap().as_bytes()[1..];
            let reference = Allele::Sequence(if strand == Strand::Reverse {
                fastvep_genome::codon::reverse_complement(reference)
            } else { reference.to_vec() });
            let change = predictor.predict_coding_consequence(
                &reference, &Allele::Deletion, &tr,
                Some(2), Some(30), Some(12), Some(40),
            ).unwrap();
            assert!(change.additional.contains(&Consequence::StartLost), "{strand:?}: {:?}", change.additional);
            assert!(!change.additional.contains(&Consequence::StartRetainedVariant), "{strand:?}: {:?}", change.additional);
        }
    }

    #[test]
    fn deleting_c3_can_retain_the_initiator_using_the_next_base() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let got = consequences_for(strand, 13, 13, &Allele::from_str("G"), &Allele::Deletion);
            assert!(got.contains(&Consequence::FrameshiftVariant));
            assert!(
                got.contains(&Consequence::StartRetainedVariant),
                "{strand:?}: {got:?}"
            );
            assert!(
                !got.contains(&Consequence::StartLost),
                "{strand:?}: {got:?}"
            );
        }
    }

    #[test]
    fn a_c3_deletion_cannot_use_the_c1_utr_repeat_exception() {
        let predictor = ConsequencePredictor::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            let mut tr = make_boundary_transcript(strand);
            tr.spliced_seq.as_mut().unwrap().replace_range(9..10, "A");
            tr.spliced_seq.as_mut().unwrap().replace_range(13..14, "C");
            tr.translateable_seq
                .as_mut()
                .unwrap()
                .replace_range(3..4, "C");
            let reference = Allele::from_str(if strand == Strand::Forward { "G" } else { "C" });
            let change = predictor
                .predict_coding_consequence(
                    &reference,
                    &Allele::Deletion,
                    &tr,
                    Some(3),
                    Some(3),
                    Some(13),
                    Some(13),
                )
                .unwrap();
            assert!(change.additional.contains(&Consequence::StartLost));
            assert!(!change
                .additional
                .contains(&Consequence::StartRetainedVariant));
        }
    }

    #[test]
    fn insertion_in_a_split_start_codon_uses_the_surviving_endpoint() {
        let predictor = ConsequencePredictor::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            let tr = make_boundary_transcript(strand);
            let (cds_start, cds_end, cdna_start, cdna_end, alt) = match strand {
                Strand::Forward => (Some(3), None, Some(13), None, "C"),
                Strand::Reverse => (None, Some(3), None, Some(13), "G"),
            };
            let change = predictor
                .predict_coding_consequence(
                    &Allele::Deletion,
                    &Allele::from_str(alt),
                    &tr,
                    cds_start,
                    cds_end,
                    cdna_start,
                    cdna_end,
                )
                .unwrap();
            assert!(
                change.consequence == Consequence::StartLost
                    || change.additional.contains(&Consequence::StartLost),
                "{strand:?}: {:?}, {:?}",
                change.consequence,
                change.additional
            );
        }
    }

    /// VEP retains the start when removing UTR bases leaves the full CDS suffix intact.
    #[test]
    fn boundary_edit_can_preserve_the_entire_cds_suffix() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let got = consequences_for(
                strand,
                9,
                13,
                &Allele::from_str("CCATG"),
                &Allele::from_str(if strand == Strand::Forward {
                    "ATG"
                } else {
                    "CAT"
                }),
            );
            assert!(
                got.contains(&Consequence::StartRetainedVariant),
                "{strand:?}: {got:?}"
            );
            assert!(got.contains(&Consequence::StartLost));
        }
    }

    #[test]
    fn retained_deletion_at_start_can_preserve_the_cds_suffix() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let mut tr = make_boundary_transcript(strand);
            tr.spliced_seq.as_mut().unwrap().replace_range(9..10, "A");
            let (start, end) = genomic_span(strand, 11, 12);
            let reference = Allele::from_str("AT");
            let alternate = Allele::from_str(if strand == Strand::Forward { "T" } else { "A" });
            let result = ConsequencePredictor::default().predict(
                &GenomicPosition::new("chr1", start, end, Strand::Forward),
                &reference, &[alternate, Allele::from_str("TA")], &[&tr], None,
            );
            let terms = &result.transcript_consequences[0].allele_consequences[0].consequences;
            for term in [Consequence::FrameshiftVariant, Consequence::StartLost, Consequence::StartRetainedVariant] {
                assert!(terms.contains(&term), "{strand:?}: {terms:?}");
            }
        }
    }

    #[test]
    fn equal_length_boundary_replacement_evaluates_start_independently_of_stop() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let tr = make_boundary_transcript(strand);
            for (lo, hi, retains_start) in [(2, 40, true), (2, 41, false), (38, 41, false)] {
                let reference = tr.spliced_seq.as_ref().unwrap().as_bytes()
                    [(lo - 1) as usize..hi as usize].to_vec();
                let mut alternate = reference.clone();
                *alternate.last_mut().unwrap() = b'G';
                let orient = |bases: Vec<u8>| Allele::Sequence(if strand == Strand::Forward {
                    bases
                } else { bases.into_iter().rev().map(complement).collect() });
                let terms = consequences_for(strand, lo, hi, &orient(reference), &orient(alternate));
                assert_eq!(terms.contains(&Consequence::StartRetainedVariant), retains_start,
                    "{strand:?} {lo}..{hi}: {terms:?}");
                assert!(!terms.contains(&Consequence::StopRetainedVariant), "{terms:?}");
                assert!(!terms.contains(&Consequence::StopLost), "{terms:?}");
            }
        }
    }

    #[test]
    fn intron_spanning_replacement_tests_start_at_the_cds_suffix() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let tr = make_boundary_transcript(strand);
            // The 21 genomic bases map to five cDNA bases after removing the
            // intron. Keeping genomic REF/ALT lengths equal does not keep the
            // edited cDNA length unchanged.
            for (text, retained) in [
                (format!("CCATG{}", "G".repeat(16)), false),
                (format!("CCCCC{}ATG", "G".repeat(13)), true),
            ] {
                let bases = if strand == Strand::Forward { text.into_bytes() }
                    else { text.bytes().rev().map(complement).collect() };
                let change = ConsequencePredictor::default().predict_cds_boundary_consequence(
                    &Allele::Sequence(vec![b'C'; 21]), &Allele::Sequence(bases),
                    &tr, Some(9), Some(13),
                );
                if retained {
                    let change = change.unwrap();
                    assert_eq!(change.consequence, Consequence::StartRetainedVariant);
                    assert!(change.additional.contains(&Consequence::StartLost));
                } else { assert!(change.is_none()); }
            }
        }
    }

    /// A variant that stays inside one region keeps exactly the term it had, on
    /// both strands - the additive region tests must not add a second one.
    #[test]
    fn variants_inside_a_single_region_gain_no_extra_term() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // Wholly in the 3' UTR.
            let utr3 =
                consequences_for(strand, 60, 62, &Allele::from_str("GGG"), &Allele::Deletion);
            assert_eq!(
                utr3,
                vec![Consequence::ThreePrimeUtrVariant],
                "{strand:?}: 3'UTR deletion"
            );

            // Wholly in the 5' UTR.
            let utr5 = consequences_for(strand, 3, 5, &Allele::from_str("CCC"), &Allele::Deletion);
            assert_eq!(
                utr5,
                vec![Consequence::FivePrimeUtrVariant],
                "{strand:?}: 5'UTR deletion"
            );

            // Wholly inside the CDS, clear of both terminator and initiator.
            let cds = consequences_for(strand, 20, 22, &Allele::from_str("TGC"), &Allele::Deletion);
            assert_eq!(
                cds,
                vec![Consequence::InframeDeletion],
                "{strand:?}: in-CDS deletion"
            );
        }
    }

    /// A transcript whose annotation does not claim to carry the terminator
    /// cannot lose it. VEP checks `cds_end_NF` before its stop predicates and
    /// `cds_start_NF` before its start ones.
    #[test]
    fn a_cds_end_nf_transcript_does_not_report_stop_lost() {
        let predictor = ConsequencePredictor::default();
        let mut tr = make_boundary_transcript(Strand::Forward);
        tr.flags = vec!["cds_end_NF".to_string()];
        let (g_start, g_end) = genomic_span(Strand::Forward, 38, 44);
        let pos = GenomicPosition::new("chr1", g_start, g_end, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("TAACCTG"),
            &[Allele::from_str("AC")],
            &[&tr],
            None,
        );
        let got = &result.transcript_consequences[0].allele_consequences[0].consequences;
        assert!(
            !got.contains(&Consequence::StopLost),
            "cds_end_NF transcript should not report stop_lost, got {got:?}"
        );
        assert!(got.contains(&Consequence::ThreePrimeUtrVariant));
    }

    /// The codon at `codon_start` read the obvious way: materialise the whole
    /// edited cDNA, then slice. `edited_codon` must agree with this in every
    /// case; the only reason it does not work this way is that the sequence
    /// being edited is the whole transcript and the loop runs per variant.
    fn naive_edited_codon(
        seq: &str,
        strand: Strand,
        cdna_lo: u64,
        cdna_hi: u64,
        alt: &Allele,
        codon_start: u64,
    ) -> Option<[u8; 3]> {
        let bytes = seq.as_bytes();
        let alt_tx: Vec<u8> = match alt {
            Allele::Sequence(b) => match strand {
                Strand::Forward => b.clone(),
                Strand::Reverse => b.iter().rev().map(|&x| complement(x)).collect(),
            },
            Allele::Deletion => Vec::new(),
            _ => return None,
        };
        let lo = (cdna_lo - 1) as usize;
        let hi = cdna_hi as usize; // exclusive
        if hi > bytes.len() || lo > hi {
            return None;
        }
        let mut edited: Vec<u8> = Vec::with_capacity(bytes.len());
        edited.extend_from_slice(&bytes[..lo]);
        edited.extend_from_slice(&alt_tx);
        edited.extend_from_slice(&bytes[hi..]);

        let start = (codon_start - 1) as usize;
        if start + 3 > edited.len() {
            return None;
        }
        let mut out = [0u8; 3];
        for (i, slot) in out.iter_mut().enumerate() {
            *slot = edited[start + i].to_ascii_uppercase();
        }
        Some(out)
    }

    /// `edited_codon` reads a codon out of the edited cDNA without building it.
    /// Check it against the build-it-and-slice version over a wide sweep of
    /// shapes, on both strands, including the ones the boundary path never
    /// generates - the arithmetic should not depend on that.
    #[test]
    fn edited_codon_agrees_with_building_the_edited_sequence() {
        let predictor = ConsequencePredictor::default();
        // A deterministic xorshift, so a failure is reproducible.
        let mut state: u64 = 0x2545_F491_4F6C_DD1D;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };

        let bases = *b"ACGT";
        let mut checked = 0usize;
        for strand in [Strand::Forward, Strand::Reverse] {
            let mut tr = make_boundary_transcript(strand);
            let seq = tr.spliced_seq.clone().unwrap();
            let seq_len = seq.len() as u64;

            for _ in 0..4_000 {
                let lo = 1 + next() % (seq_len - 1);
                let span = next() % 12; // 0..11 extra bases
                let hi = (lo + span).min(seq_len);
                let alt_len = (next() % 6) as usize; // 0 => a pure deletion
                let alt = if alt_len == 0 {
                    Allele::Deletion
                } else {
                    Allele::Sequence((0..alt_len).map(|_| bases[(next() % 4) as usize]).collect())
                };
                let codon_start = 1 + next() % (seq_len - 2);

                let got = predictor.edited_codon(&tr, lo, hi, &alt, codon_start);
                let want = naive_edited_codon(&seq, strand, lo, hi, &alt, codon_start);
                assert_eq!(
                    got, want,
                    "{strand:?}: lo={lo} hi={hi} alt={alt:?} codon_start={codon_start}"
                );
                checked += 1;
            }

            // And with no sequence loaded there is nothing to read.
            tr.spliced_seq = None;
            assert_eq!(
                predictor.edited_codon(&tr, 10, 12, &Allele::Deletion, 11),
                None,
                "a transcript with no sequence has no edited codon"
            );
        }
        assert_eq!(checked, 8_000);
    }

    /// A CDS annotated as beginning part-way through a codon has no complete
    /// initiator, so nothing at its start can be `start_lost`. Without the
    /// phase guard every length-changing variant reaching that end would be,
    /// because the three bases there never read as ATG.
    #[test]
    fn a_phase_offset_transcript_does_not_report_start_lost() {
        let predictor = ConsequencePredictor::default();
        let mut tr = make_boundary_transcript(Strand::Forward);
        tr.codon_table_start_phase = 2;
        let (g_start, g_end) = genomic_span(Strand::Forward, 9, 13);
        let pos = GenomicPosition::new("chr1", g_start, g_end, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("CCATG"),
            &[Allele::from_str("TT")],
            &[&tr],
            None,
        );
        let got = &result.transcript_consequences[0].allele_consequences[0].consequences;
        assert!(
            !got.contains(&Consequence::StartLost),
            "a phase-offset CDS has no initiator to lose, got {got:?}"
        );
        // The same transcript with a complete initiator does report it, so the
        // guard is what makes the difference and not the variant.
        let tr0 = make_boundary_transcript(Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("CCATG"),
            &[Allele::from_str("TT")],
            &[&tr0],
            None,
        );
        assert!(result.transcript_consequences[0].allele_consequences[0]
            .consequences
            .contains(&Consequence::StartLost));
    }

    /// An insertion sitting exactly on either coding boundary is a UTR variant
    /// and nothing more. Its span is the two bases it sits between, so it never
    /// deletes a codon, and the boundary path must not claim one was lost.
    #[test]
    fn an_insertion_on_a_coding_boundary_reports_only_the_utr_term() {
        let predictor = ConsequencePredictor::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            let tr = make_boundary_transcript(strand);
            // The initiator is cDNA 11 and the terminator ends at cDNA 40.
            for (cdna_flank, forbidden) in
                [(11u64, Consequence::StartLost), (41, Consequence::StopLost)]
            {
                let genomic = match strand {
                    Strand::Forward => 999 + cdna_flank,
                    // On the reverse strand the preceding transcript base is
                    // one genomic coordinate higher, so the zero-length
                    // interval starts there.
                    Strand::Reverse => 1202 - cdna_flank,
                };
                // Ensembl's zero-length interval: end = start - 1.
                let pos = GenomicPosition::new("chr1", genomic, genomic - 1, Strand::Forward);
                let result = predictor.predict(
                    &pos,
                    &Allele::Deletion,
                    &[Allele::from_str("TTT")],
                    &[&tr],
                    None,
                );
                let got = &result.transcript_consequences[0].allele_consequences[0].consequences;
                assert!(
                    !got.contains(&forbidden),
                    "{strand:?} insertion at cDNA {cdna_flank}: unexpected {forbidden:?} in {got:?}"
                );
            }
        }
    }
    /// A delins in CDS terms, run through the full predictor.
    ///
    /// `cds_lo` is the CDS coordinate of the first replaced base in transcript
    /// order; `reference` and `replacement` are read in transcript order too, so
    /// the reverse-strand case is written the same way and the helper supplies
    /// the reverse complement the VCF would carry.
    fn delins_at(
        strand: Strand,
        cds_lo: u64,
        reference: &str,
        replacement: &str,
    ) -> AlleleConsequenceResult {
        let tr = make_boundary_transcript(strand);
        // CDS n is cDNA n + 10 on this transcript.
        let (lo, hi) = genomic_span(
            strand,
            cds_lo + 10,
            cds_lo + 10 + reference.len() as u64 - 1,
        );
        let orient = |s: &str| -> Allele {
            Allele::Sequence(match strand {
                Strand::Forward => s.as_bytes().to_vec(),
                Strand::Reverse => fastvep_genome::codon::reverse_complement(s.as_bytes()),
            })
        };
        let pos = GenomicPosition::new("chr1", lo, hi, Strand::Forward);
        let alt = orient(replacement);
        let result = ConsequencePredictor::default().predict(
            &pos,
            &orient(reference),
            std::slice::from_ref(&alt),
            &[&tr],
            None,
        );
        result.transcript_consequences[0].allele_consequences[0].clone()
    }

    #[test]
    fn padded_repeat_deletion_displays_only_the_minimized_codon() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let tr = make_boundary_transcript(strand);
            let (lo, hi) = genomic_span(strand, 14, 22);
            let orient = |s: &[u8]| {
                Allele::Sequence(if strand == Strand::Forward {
                    s.to_vec()
                } else {
                    fastvep_genome::codon::reverse_complement(s)
                })
            };
            let ac = ConsequencePredictor::default()
                .display_coding_change(
                    &GenomicPosition::new("chr1", lo, hi, Strand::Forward),
                    &orient(b"GCTGCTGCT"),
                    &orient(b"GCTGCT"),
                    &tr,
                    false,
                )
                .unwrap();
            assert_eq!(ac.codons, Some(("GCT".into(), "-".into())));
            assert_eq!(ac.amino_acids, Some(("A".into(), "-".into())));
        }
    }

    /// A delins that replaces residues without preserving the reference ones at
    /// either end of the replacement is `protein_altering_variant`, not an
    /// in-frame indel.
    ///
    /// Choosing the term from the direction of the length change - which is what
    /// this code did - called all of these `inframe_deletion` or
    /// `inframe_insertion`. Real VEP 115.1 over the 156 in-frame delins in the
    /// ClinVar 2-star set gives `protein_altering_variant` on 1,231 transcript
    /// rows, `inframe_insertion` on 81 and `inframe_deletion` on none.
    #[test]
    fn a_delins_that_replaces_residues_is_protein_altering_on_either_strand() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // CDS 4-9 is Ala2 Ala3 (`GCTGCT`), replaced by `TGG` (Trp).
            let ac = delins_at(strand, 4, "GCTGCT", "TGG");
            assert!(
                ac.consequences
                    .contains(&Consequence::ProteinAlteringVariant),
                "{strand:?}: expected protein_altering_variant, got {:?}",
                ac.consequences
            );
            assert_eq!(
                ac.amino_acids,
                Some(("AA".to_string(), "W".to_string())),
                "{strand:?}"
            );
            // Every base of the window is replaced, so none of it stays lower.
            assert_eq!(
                ac.codons,
                Some(("GCTGCT".to_string(), "TGG".to_string())),
                "{strand:?}"
            );
        }
    }

    /// The replacement preserving the reference residues at one end of itself is
    /// what makes a delins an in-frame indel. The codon rendering keeps the
    /// unchanged flanks of the window lowercase and uppercases exactly the
    /// replaced and inserted bases - `gGg/gTCCCg`, not "uppercase from the first
    /// difference", which would lose the trailing lowercase base.
    #[test]
    fn a_delins_that_extends_the_reference_residues_is_an_inframe_insertion() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // CDS 5 is the middle base of Ala2's `gCt`. Replacing `C` with
            // `CTGC` rebuilds the window as `GCTGCT`, so the alt peptide still
            // starts with the reference residue and gains one.
            let ac = delins_at(strand, 5, "C", "CTGC");
            assert!(
                ac.consequences.contains(&Consequence::InframeInsertion),
                "{strand:?}: expected inframe_insertion, got {:?}",
                ac.consequences
            );
            assert_eq!(
                ac.amino_acids,
                Some(("A".to_string(), "AA".to_string())),
                "{strand:?}"
            );
            assert_eq!(
                ac.codons,
                Some(("gCt".to_string(), "gCTGCt".to_string())),
                "{strand:?}"
            );
        }
    }

    /// VEP's first `ref_eq_alt_sequence` clause treats a length-changing edit
    /// that preserves its one reference residue and introduces a later stop as
    /// an in-frame insertion with the stop retained.
    #[test]
    fn a_preserved_first_residue_and_later_stop_matches_vep() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let ac = delins_at(strand, 4, "GCT", "GCTTAGA");
            assert!(
                ac.consequences.contains(&Consequence::InframeInsertion),
                "{strand:?}: expected inframe_insertion, got {:?}",
                ac.consequences
            );
            assert!(
                ac.consequences.contains(&Consequence::StopRetainedVariant),
                "{strand:?}: expected stop_retained_variant, got {:?}",
                ac.consequences
            );
            assert!(
                !ac.consequences.contains(&Consequence::FrameshiftVariant)
                    && !ac.consequences.contains(&Consequence::StopGained),
                "{strand:?}: VEP suppresses frameshift and stop_gained here, got {:?}",
                ac.consequences
            );
            assert_eq!(ac.impact, Impact::Moderate, "{strand:?}");
        }
    }

    #[test]
    fn an_insertion_in_a_terminal_repeat_can_retain_the_reference_protein() {
        let predictor = ConsequencePredictor::default();
        let mut window = CodonWindow {
                peptides_defined: true,
                alt_dna_unambiguous: true,
            ref_aas: String::new(),
            alt_aas: "E".into(),
            ref_window: Vec::new(),
            alt_window: b"GAG".to_vec(),
            ref_codons: "-".into(),
            alt_codons: "GAG".into(),
            ref_len: 0,
            alt_len: 3,
            tl_start: 4,
            tl_end: 3,
            partial_codon: false,
        };

        let retained = predictor.terms_for_window(&window, false, Some("MEEEE*"), false, false, false, None);
        assert_eq!(retained.consequence, Consequence::InframeInsertion);
        assert_eq!(retained.additional, vec![Consequence::StopRetainedVariant]);

        let displaced = predictor.terms_for_window(&window, false, Some("MEEDE*"), false, false, false, None);
        assert_eq!(displaced.consequence, Consequence::InframeInsertion);
        assert!(displaced.additional.is_empty());
        window.tl_start = 6;
        window.tl_end = 5;
        window.alt_aas = "*CX".into();
        window.alt_window = b"TGATGTA".to_vec();
        window.alt_codons = "TGATGTA".into();
        window.alt_len = 7;
        let terminal = predictor.terms_for_window(&window, false, Some("MEEEE*"), false, false, false, None);
        assert_eq!(terminal.consequence, Consequence::InframeInsertion);
        assert_eq!(terminal.additional, vec![Consequence::StopRetainedVariant]);
    }

    #[test]
    fn initiator_spanning_deletion_keeps_the_later_start_loss_predicate() {
        let window = CodonWindow {
                peptides_defined: true,
                alt_dna_unambiguous: true,
            ref_aas: "MAFEV".into(), alt_aas: "M".into(),
            ref_window: b"ATGGCTTTCGAGGTG".to_vec(), alt_window: b"ATG".to_vec(),
            ref_codons: "aTGGCTTTCGAGGtg".into(), alt_codons: "atg".into(),
            ref_len: 12, alt_len: 0, tl_start: 1, tl_end: 5, partial_codon: false,
        };
        let change = ConsequencePredictor::default().terms_for_window(
            &window, true, None, false, true, true, None,
        );
        assert_eq!(change.consequence, Consequence::StartLost);
        assert!(change.additional.contains(&Consequence::InframeDeletion));
        assert!(change.additional.contains(&Consequence::StartRetainedVariant));
    }

    #[test]
    fn a_replacement_past_an_incomplete_peptide_uses_its_available_suffix() {
        let predictor = ConsequencePredictor::default();
        let window = CodonWindow {
                peptides_defined: true,
                alt_dna_unambiguous: true,
            ref_aas: "SX".into(),
            alt_aas: "SX".into(),
            ref_window: b"AGCN".to_vec(),
            alt_window: b"AGCN".to_vec(),
            ref_codons: "AGCN".into(),
            alt_codons: "AGCN".into(),
            ref_len: 2,
            alt_len: 2,
            tl_start: 3,
            tl_end: 4,
            partial_codon: false,
        };

        let retained = predictor.terms_for_window(&window, false, Some("MSS"), false, false, false, None);
        assert_eq!(retained.consequence, Consequence::StopRetainedVariant);
        assert!(retained.additional.is_empty());
    }

    /// Ensembl's insertion fixture at the third base of the initiator retains
    /// the original `ATG` downstream while changing the codon at position one.
    #[test]
    fn an_insertion_can_both_displace_and_retain_the_start_codon() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let tr = make_boundary_transcript(strand);
            // Insert transcript-oriented `CAT` between CDS bases two and three.
            let (start, end) = genomic_span(strand, 13, 12);
            let inserted = match strand {
                Strand::Forward => Allele::from_str("CAT"),
                Strand::Reverse => Allele::from_str("ATG"),
            };
            let result = ConsequencePredictor::default().predict(
                &GenomicPosition::new("chr1", start, end, Strand::Forward),
                &Allele::Deletion,
                std::slice::from_ref(&inserted),
                &[&tr],
                None,
            );
            let got = &result.transcript_consequences[0].allele_consequences[0].consequences;
            for expected in [Consequence::StartLost, Consequence::StartRetainedVariant] {
                assert!(
                    got.contains(&expected),
                    "{strand:?}: {expected:?} missing from {got:?}"
                );
            }
        }
    }

    #[test]
    fn deleting_cds_base_one_has_both_vep_start_terms() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let mut tr = make_boundary_transcript(strand);
            // VEP's indel predicate compares the edited transcript's CDS-length
            // suffix with the original CDS. If the last 5' UTR base repeats
            // CDS base one, c.1del leaves that suffix unchanged.
            tr.spliced_seq.as_mut().unwrap().replace_range(9..10, "A");
            let reference = match strand {
                Strand::Forward => Allele::from_str("A"),
                Strand::Reverse => Allele::from_str("T"),
            };
            let (start, end) = genomic_span(strand, 11, 11);
            let result = ConsequencePredictor::default().predict(
                &GenomicPosition::new("chr1", start, end, Strand::Forward),
                &reference,
                &[Allele::Deletion],
                &[&tr],
                None,
            );
            let got = &result.transcript_consequences[0].allele_consequences[0].consequences;
            for expected in [
                Consequence::FrameshiftVariant,
                Consequence::StartLost,
                Consequence::StartRetainedVariant,
            ] {
                assert!(
                    got.contains(&expected),
                    "{strand:?}: {expected:?} missing from {got:?}"
                );
            }

            let without_repeat = consequences_for(strand, 11, 11, &reference, &Allele::Deletion);
            assert!(
                !without_repeat.contains(&Consequence::StartRetainedVariant),
                "{strand:?}: a different upstream base must not retain the start: {without_repeat:?}"
            );
        }
    }

    /// A delins that introduces a terminator earns both terms. Ensembl evaluates
    /// each predicate independently and keeps all that hold, so the new stop
    /// gives `stop_gained` and the changed residue count gives
    /// `protein_altering_variant` - 164 of the 1,803 coding rows in the ClinVar
    /// in-frame delins set are that pair in real VEP.
    ///
    /// `Amino_acids` keeps the whole translated window including what follows
    /// the new terminator, which is what VEP reports (`SL/MEP*S`).
    #[test]
    fn a_delins_introducing_a_terminator_reports_both_terms() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let ac = delins_at(strand, 4, "GCTGCT", "GCCTAGGCC");
            for expected in [Consequence::StopGained, Consequence::ProteinAlteringVariant] {
                assert!(
                    ac.consequences.contains(&expected),
                    "{strand:?}: expected {expected:?}, got {:?}",
                    ac.consequences
                );
            }
            assert_eq!(ac.impact, Impact::High, "{strand:?}");
            assert_eq!(
                ac.amino_acids,
                Some(("AA".to_string(), "A*A".to_string())),
                "{strand:?}"
            );
        }
    }

    /// A replacement whose peptide *begins* with the terminator is `stop_gained`
    /// alone: `protein_altering_variant` declines when either peptide starts
    /// with `*`. That is what separates VEP's `HQ/*` rows, which carry one term,
    /// from its `SL/MEP*S` rows, which carry two.
    #[test]
    fn a_delins_whose_peptide_begins_with_a_terminator_is_stop_gained_alone() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let ac = delins_at(strand, 4, "GCTGCT", "TAG");
            assert_eq!(
                ac.amino_acids,
                Some(("AA".to_string(), "*".to_string())),
                "{strand:?}"
            );
            assert!(
                ac.consequences.contains(&Consequence::StopGained),
                "{strand:?}: got {:?}",
                ac.consequences
            );
            assert!(
                !ac.consequences
                    .contains(&Consequence::ProteinAlteringVariant),
                "{strand:?}: protein_altering_variant must decline, got {:?}",
                ac.consequences
            );
        }
    }

    /// A delins over the initiator is `start_lost`, and that outranks the length
    /// change. Ensembl asks whether the reference residues survived at either
    /// end of the replacement, not whether some ATG still sits at the coding
    /// start.
    #[test]
    fn a_delins_over_the_initiator_is_start_lost() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let ac = delins_at(strand, 1, "ATGGCT", "CCC");
            assert!(
                ac.consequences.contains(&Consequence::StartLost),
                "{strand:?}: got {:?}",
                ac.consequences
            );
            assert_eq!(ac.impact, Impact::High, "{strand:?}");
        }
    }

    /// A change whose reference bases are not contiguous in CDS space has no
    /// codon window, and building one from the low CDS coordinate and the allele
    /// length translates codons the variant never touched. VEP reports
    /// `coding_sequence_variant` with no residues for these, which is what
    /// resolving to no coding term produces.
    ///
    /// This is what a change straddling a splice junction looks like from here:
    /// its two ends are both in the CDS, but further apart than it has bases.
    /// Over a 6,600-variant ClinVar sample it was 156 rows of invented residue
    /// change, ten of them a false `stop_gained`.
    #[test]
    fn a_change_that_is_not_contiguous_in_the_cds_names_no_residues() {
        let predictor = ConsequencePredictor::default();
        let tr = make_boundary_transcript(Strand::Forward);
        for (reference, replacement) in [
            ("GCTGCT", "TGG"),    // delins
            ("GCTGCT", "TGGCCC"), // equal-length MNV
        ] {
            let change = predictor.predict_coding_consequence(
                &Allele::Sequence(reference.as_bytes().to_vec()),
                &Allele::Sequence(replacement.as_bytes().to_vec()),
                &tr,
                Some(4),
                Some(20), // 17 CDS positions for six reference bases
                Some(14),
                Some(30),
            );
            assert!(
                change.is_none(),
                "{reference}/{replacement}: a non-contiguous change must resolve to \
                 no coding term, got {:?}",
                change.map(|c| c.consequence)
            );
        }
    }

    /// An insertion that falls exactly between two codons replaces no codon, so
    /// Ensembl's window is empty on the reference side.
    ///
    /// It reports `-/HENKTKGD` and `-/CATGAG...`, not the flanking residue
    /// repeated on both sides. Anchoring to the preceding codon instead was 791
    /// of 1,808 in-frame insertion rows over a 6,600-variant ClinVar sample, and
    /// it also cost the `dup` collapsing downstream: with the flanking residue
    /// in the way, `p.His553_Asp560dup` came out as an eight-residue `ins`.
    #[test]
    fn a_codon_aligned_insertion_names_no_reference_residue() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // Insert three bases after CDS 3, on a codon boundary.
            let ac = insertion_at(strand, 3, "GGG");
            assert_eq!(
                ac.amino_acids,
                Some(("-".to_string(), "G".to_string())),
                "{strand:?}: {:?}",
                ac.amino_acids
            );
            assert_eq!(
                ac.codons,
                Some(("-".to_string(), "GGG".to_string())),
                "{strand:?}"
            );
            assert!(
                ac.consequences.contains(&Consequence::InframeInsertion),
                "{strand:?}: got {:?}",
                ac.consequences
            );

            // One base further in, the insertion sits inside a codon and the
            // window is that codon.
            let ac = insertion_at(strand, 4, "GGG");
            let (r, a) = ac.amino_acids.clone().unwrap();
            assert_eq!(r.len(), 1, "{strand:?}: {r}/{a}");
            assert_eq!(a.len(), 2, "{strand:?}: {r}/{a}");
        }
    }

    /// Ensembl marks a codon the window does not complete with `X`, but only
    /// when the codons before it did not already end translation.
    ///
    /// `if($partial_codon && $pep ne '*') { $pep .= 'X' }`
    /// (`TranscriptVariationAllele.pm` release/115). Appending it
    /// unconditionally reported `Y/*X` where VEP writes `Y/*` - 288 rows of a
    /// 6,600-variant ClinVar sample.
    #[test]
    fn a_partial_codon_after_a_terminator_adds_no_placeholder() {
        // The fixture's CDS begins ATG GCT GCT ...; inserting one base after
        // CDS 3 shifts the frame from codon 2 on.
        let ac = insertion_at(Strand::Forward, 3, "T");
        let (_, alt) = ac.amino_acids.clone().unwrap();
        assert!(
            alt.ends_with('X') || alt == "*",
            "a frameshift window ends in X unless it ended in a terminator: {alt}"
        );
        assert!(
            !alt.ends_with("*X"),
            "nothing is translated past a terminator: {alt}"
        );
    }

    /// An insertion on an exon's edge has one end in the intron, so only one of
    /// its two CDS coordinates exists - and it is still a frameshift.
    ///
    /// Requiring both coordinates left 39 rows of a 6,600-variant ClinVar sample
    /// as `coding_sequence_variant`, LOW where VEP says HIGH, every one of them
    /// a frameshift at the first or last base of an exon.
    #[test]
    fn an_insertion_on_an_exon_edge_is_still_a_codon_edit() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let tr = make_boundary_transcript(strand);
            // The fixture's first exon ends at cDNA 20, so an insertion between
            // cDNA 20 and 21 has one end in the intron.
            let (edge, _) = genomic_span(strand, 20, 20);
            let (next, _) = genomic_span(strand, 21, 21);
            let (lo, hi) = (edge.min(next), edge.max(next));
            let pos = GenomicPosition::new("chr1", hi, lo, Strand::Forward);
            let alt = Allele::Sequence(b"T".to_vec());
            let result = ConsequencePredictor::default().predict(
                &pos,
                &Allele::Deletion,
                std::slice::from_ref(&alt),
                &[&tr],
                None,
            );
            let got = &result.transcript_consequences[0].allele_consequences[0].consequences;
            assert!(
                got.contains(&Consequence::FrameshiftVariant),
                "{strand:?}: got {got:?}"
            );
        }
    }

    /// An insertion of `bases` after CDS position `cds_lo`, on either strand.
    fn insertion_at(strand: Strand, cds_lo: u64, bases: &str) -> AlleleConsequenceResult {
        let tr = make_boundary_transcript(strand);
        // CDS n is cDNA n + 10 on this transcript, and an insertion is the
        // zero-length interval between two adjacent bases.
        let (a, _) = genomic_span(strand, cds_lo + 10, cds_lo + 10);
        let (b, _) = genomic_span(strand, cds_lo + 11, cds_lo + 11);
        let pos = GenomicPosition::new("chr1", a.max(b), a.min(b), Strand::Forward);
        let alt = Allele::Sequence(match strand {
            Strand::Forward => bases.as_bytes().to_vec(),
            Strand::Reverse => fastvep_genome::codon::reverse_complement(bases.as_bytes()),
        });
        let result = ConsequencePredictor::default().predict(
            &pos,
            &Allele::Deletion,
            std::slice::from_ref(&alt),
            &[&tr],
            None,
        );
        result.transcript_consequences[0].allele_consequences[0].clone()
    }
}
