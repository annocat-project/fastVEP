# VEP 115 concordance change inventory

This inventory groups the fork changes into **113 review families**. Several families contain multiple repairs; others cover supporting cache or output work. This is not a count of independent bugs, commits, or failures in current upstream.

The source comparison is `0e13c5bdb92f22a5f780cf314699f8402fb8383e` (recorded upstream base) to `99c427a551ab93506fd427c5e0d2cb0c679ff28b` (reviewed fork). Later upstream changes may already address some of these behaviors.

Start with the [four-case review guide](UPSTREAM-REVIEW.md). Its four groups have verified before/after evidence using retained fork builds. The rest of this inventory links related source tests or implementation locations; it does not claim an isolated before/after demonstration for every row.

## Evidence and limits

- The rationale describes the intended behavior. Related tests are review entry points, not proof that each test was introduced by that repair.
- The [full suite](FULL-SUITE.md) contains 216 input banks, 2,531,230 record observations and 2,451,891 distinct exact record keys. These are not independent statistical trials or a per-family denominator. The contracts define compared fields, normalizers and exceptions.
- This documentation update adds no annotation run. The current validation and before-fix evidence for the four focused groups remain in [their review index](focused/review-index.json). Broader historical bank-to-repair links have not been reconstructed here.

## Source review starting points

The following files were read from the unmodified VEP 115.2 image extraction. They identify routines to inspect; this inventory does not provide a new complete source trace for every family.

- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`: Consequence, start/stop, splice and complete-overlap predicates.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`: Peptide preparation, HGVS protein generation and allele clipping.
- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`: Overlap preclassification and insertion endpoint ordering.
- `Bio/EnsEMBL/Variation/Utils/Sequence.pm`: HGVS notation and inversion classification.
- `Bio/EnsEMBL/Transcript.pm`: Transcript translation and reference initiator adjustment.

## Repair families

### Alleles and record semantics

#### A01 — Prefix-first input minimization

Preserve VEP's parsed interval independently from later HGVS shifting.

Source: [vcf.rs](../../crates/fastvep-io/src/vcf.rs).

Related tests: [`lowercase_alleles_are_uppercased_after_trimming`](../../crates/fastvep-io/src/vcf.rs#L366), [`test_parse_insertion`](../../crates/fastvep-io/src/vcf.rs#L426), [`test_parse_deletion`](../../crates/fastvep-io/src/vcf.rs#L439), [`biallelic_indels_trim_prefix_then_suffix_before_annotation`](../../crates/fastvep-io/src/vcf.rs#L510).

Related source tests; no isolated before/after replay added for this inventory.

#### A02 — Case conversion after trimming

Uppercasing must not trigger a second minimization of an already parsed allele.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`parsed_input_does_not_repeat_minimization_after_case_conversion`](../../crates/fastvep-consequence/src/predictor.rs#L2665).

Related source tests; no isolated before/after replay added for this inventory.

#### A03 — Padded allele direction and distance

Use the left-first input coordinates for transcript direction, intron membership and the distance limit.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`padded_downstream_indel_uses_vep_left_first_distance`](../../crates/fastvep-consequence/src/predictor.rs#L2343), [`padded_reverse_insertion_uses_vep_coordinates_for_direction`](../../crates/fastvep-consequence/src/predictor.rs#L2367), [`padded_intronic_deletion_does_not_reach_the_exon`](../../crates/fastvep-consequence/src/predictor.rs#L2380), [`padded_indel_crossing_distance_limit_is_intergenic`](../../crates/fastvep-consequence/src/predictor.rs#L2398), [`vep_input_position_trims_suffix_after_prefix`](../../crates/fastvep-consequence/src/predictor.rs#L2419), [`vep_input_position_preserves_equal_length_and_multi_alt_intervals`](../../crates/fastvep-consequence/src/predictor.rs#L2433).

Related source tests; no isolated before/after replay added for this inventory.

#### A04 — Record-wide SNP and indel shape

HGVS eligibility depends on the whole record's allele shape; preserve mixed-ALT replacement intervals.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs).

Related tests: [`hgvs_shift_requires_a_whole_record_insertion_or_deletion`](../../crates/fastvep-annotate/src/lib.rs#L2386), [`mixed_multiallelic_frameshift_uses_retained_replacement_interval`](../../crates/fastvep-annotate/src/lib.rs#L2402), [`vep_snp_class_requires_one_base_for_every_record_allele`](../../crates/fastvep-annotate/src/lib.rs#L2570).

Related source tests; no isolated before/after replay added for this inventory.

#### A05 — Reference-equal alternate alleles

Suppress reference-equivalent annotation rows while preserving the original input record.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`reference_repeated_as_an_alternate_has_no_annotation`](../../crates/fastvep-consequence/src/predictor.rs#L2652).

Related source tests; no isolated before/after replay added for this inventory.

#### A06 — Intergenic reference-alternate handling

Apply reference-alternate filtering on the intergenic path without rewriting records.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs).

Related tests: [`intergenic_annotation_discards_reference_alt_without_rewriting_the_record`](../../crates/fastvep-annotate/src/lib.rs#L2276).

Related source tests; no isolated before/after replay added for this inventory.

#### A07 — Ambiguous input reference for HGVS

Use the genomic reference for transcript HGVS when the input reference is ambiguous.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs).

Related tests: [`ambiguous_input_reference_uses_genome_for_transcript_hgvs`](../../crates/fastvep-annotate/src/lib.rs#L2231).

Related source tests; no isolated before/after replay added for this inventory.

#### A08 — CSQ allele and reference projection

Project Allele and REF_ALLELE with VEP input rules while retaining equal-length and multiallelic representations.

Source: [output.rs](../../crates/fastvep-io/src/output.rs).

Related tests: [`csq_allele_is_minimized_without_changing_source_allele`](../../crates/fastvep-io/src/output.rs#L2269), [`csq_allele_preserves_equal_length_padded_substitution`](../../crates/fastvep-io/src/output.rs#L2279), [`csq_allele_uses_vep_left_first_minimization_in_repeats`](../../crates/fastvep-io/src/output.rs#L2287), [`csq_ref_allele_uses_vep_left_first_minimization`](../../crates/fastvep-io/src/output.rs#L2315), [`csq_ref_allele_preserves_multi_alt_and_equal_length_values`](../../crates/fastvep-io/src/output.rs#L2327), [`parsed_csq_alleles_are_not_minimized_again_after_case_conversion`](../../crates/fastvep-io/src/output.rs#L2385).

Related source tests; no isolated before/after replay added for this inventory.

#### A09 — Uploaded-allele spelling

Preserve indel case while using the declared substitution convention.

Source: [output.rs](../../crates/fastvep-io/src/output.rs).

Related tests: [`uploaded_allele_preserves_indel_case_but_uppercases_substitutions`](../../crates/fastvep-io/src/output.rs#L2366).

Related source tests; no isolated before/after replay added for this inventory.

#### A10 — HGVS offset baseline and export

Calculate HGVS_OFFSET from the input minimization baseline and retain it in transcript-specific structured output.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs).

Related tests: [`vep_hgvs_offset_baseline_trims_repeated_indels_from_the_left`](../../crates/fastvep-annotate/src/lib.rs#L2587).

Related source tests; no isolated before/after replay added for this inventory.

#### A11 — Structured HGVS offset retention

Carry each transcript's HGVS offset into JSON rather than dropping it in projection.

Source: [output.rs](../../crates/fastvep-io/src/output.rs).

Related tests: [`json_emits_hgvs_offset_with_its_transcript_consequence`](../../crates/fastvep-io/src/output.rs#L2575).

Related source tests; no isolated before/after replay added for this inventory.

#### A12 — Genomic reference for concrete REF disagreements

Use the genomic slice for HGVSc even when a concrete input REF disagrees with it; retain input alleles for parsing, class and shift eligibility.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs#L1863).

Implementation inspection; no isolated fixture or direct unit-test mapping recorded here.

### Membership and coordinates

#### M01 — Complete indexed GFF transcripts

A regional query must recover the entire overlapping transcript rather than truncated exon/CDS fragments.

Source: [gff.rs](../../crates/fastvep-cache/src/gff.rs).

Related tests: [`indexed_selection_expands_an_overlapping_transcript_to_its_complete_span`](../../crates/fastvep-cache/src/gff.rs#L1128), [`indexed_selection_merges_overlapping_transcript_spans`](../../crates/fastvep-cache/src/gff.rs#L1141).

Related source tests; no isolated before/after replay added for this inventory.

#### M02 — Processed and noncoding transcript membership

Retain source-declared processed/noncoding models, their exons and unconfirmed annotation metadata.

Source: [gff.rs](../../crates/fastvep-cache/src/gff.rs).

Related tests: [`test_parse_gff3_mane_and_metadata`](../../crates/fastvep-cache/src/gff.rs#L993), [`test_parse_gff3_ensembl_ncrna_gene_and_unconfirmed_transcript`](../../crates/fastvep-cache/src/gff.rs#L1031), [`test_parse_gff3_preserves_protein_version`](../../crates/fastvep-cache/src/gff.rs#L1045), [`test_parse_gff3_with_source_tags_every_transcript`](../../crates/fastvep-cache/src/gff.rs#L1106), [`processed_transcript_features_preserve_membership_and_exons`](../../crates/fastvep-cache/src/gff.rs#L1163).

Related source tests; no isolated before/after replay added for this inventory.

#### M03 — Nested transcript interval lookup

An indexed lookup must agree with a linear scan for nested and overlapping transcript intervals.

Source: [providers.rs](../../crates/fastvep-cache/src/providers.rs).

Related tests: [`test_indexed_provider_matches_linear_scan_for_nested_intervals`](../../crates/fastvep-cache/src/providers.rs#L530), [`test_indexed_provider_real_cache_boundaries`](../../crates/fastvep-cache/src/providers.rs#L565).

Related source tests; no isolated before/after replay added for this inventory.

#### M04 — cDNA and intron mapping

Map both strands and complete intron-overlapping spans instead of only a single endpoint.

Source: [transcript.rs](../../crates/fastvep-genome/src/transcript.rs).

Related tests: [`test_cdna_to_genomic_forward_and_reverse`](../../crates/fastvep-genome/src/transcript.rs#L584), [`test_intron_overlapping`](../../crates/fastvep-genome/src/transcript.rs#L617).

Related source tests; no isolated before/after replay added for this inventory.

#### M05 — Insertion boundary region membership

An insertion between an exon and intron must not claim membership in either region merely from one flank.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`a_span_straddling_a_coding_transcript_edge_matches_vep_utr_overlap`](../../crates/fastvep-consequence/src/predictor.rs#L2501), [`an_insertion_between_a_noncoding_exon_and_intron_is_not_in_either`](../../crates/fastvep-consequence/src/predictor.rs#L3100), [`an_insertion_on_a_coding_boundary_reports_only_the_utr_term`](../../crates/fastvep-consequence/src/predictor.rs#L4120).

Related source tests; no isolated before/after replay added for this inventory.

#### M06 — Reported cDNA input coordinates

Use the VEP input interval for reported cDNA positions, distinct from HGVS-normalized coordinates.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`reported_cdna_positions_use_vep_left_first_input_coordinates`](../../crates/fastvep-consequence/src/predictor.rs#L3072).

Related source tests; no isolated before/after replay added for this inventory.

#### M07 — Reverse protein ranges

Order projected protein endpoints correctly on reverse-strand spans.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`protein_range_sorts_a_reverse_strand_pair`](../../crates/fastvep-consequence/src/predictor.rs#L1929).

Related source tests; no isolated before/after replay added for this inventory.

#### M08 — Unknown position endpoints

Preserve partial position ranges rather than discarding known endpoints or inventing missing ones.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs).

Related tests: [`position_ranges_keep_unknown_endpoints`](../../crates/fastvep-annotate/src/lib.rs#L3150).

Related source tests; no isolated before/after replay added for this inventory.

#### M09 — VCF/JSON partial-coordinate parity

Preserve the same partial cDNA/CDS/protein ranges in VCF and structured output.

Source: [output.rs](../../crates/fastvep-io/src/output.rs).

Related tests: [`test_format_position_range`](../../crates/fastvep-io/src/output.rs#L2345), [`partial_positions_match_vcf_and_structured_json`](../../crates/fastvep-io/src/output.rs#L2494).

Related source tests; no isolated before/after replay added for this inventory.

#### M10 — Mature-miRNA ranges

Use mature-miRNA subranges in transcript order rather than only the broad noncoding-exon term.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`mature_mirna_ranges_replace_the_broader_noncoding_exon_term`](../../crates/fastvep-consequence/src/predictor.rs#L2206), [`mature_mirna_ranges_follow_reverse_strand_cdna_order`](../../crates/fastvep-consequence/src/predictor.rs#L2240).

Related source tests; no isolated before/after replay added for this inventory.

#### M11 — Whole-transcript tandem gain

Evaluate full-transcript amplification on the unminimized input span.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`complete_transcript_tandem_gain_uses_unminimized_input_span`](../../crates/fastvep-consequence/src/predictor.rs#L2173).

Related source tests; no isolated before/after replay added for this inventory.

#### M12 — Complete exon and intron rank ranges

Carry first rank, last rank and total count through prediction and output so changes spanning multiple exons or introns are not reduced to one rank.

Source: [transcript.rs](../../crates/fastvep-genome/src/transcript.rs#L276).

Implementation inspection; no isolated fixture or direct unit-test mapping recorded here.

### Splice and overlap predicates

#### S01 — Interval-tree boundary ordering

Match the selected VEP interval-tree boundary ordering when several splice regions overlap.

Source: [splice.rs](../../crates/fastvep-consequence/src/splice.rs).

Related tests: [`boundary_order_matches_vep_interval_tree_012`](../../crates/fastvep-consequence/src/splice.rs#L532).

Related source tests; no isolated before/after replay added for this inventory.

#### S02 — Short-intron exon stretching

Preclassify coding overlap using stretched exons; the VEP threshold includes thirteen-base introns.

Source: [splice.rs](../../crates/fastvep-consequence/src/splice.rs).

Related tests: [`a_short_intron_stretches_every_exon_for_preclassification`](../../crates/fastvep-consequence/src/splice.rs#L666), [`vep_frameshift_intron_threshold_includes_thirteen_bases`](../../crates/fastvep-consequence/src/splice.rs#L697).

Related source tests; no isolated before/after replay added for this inventory.

#### S03 — Coding preclassification without exon label

Preserve short-intron coding eligibility separately from an emitted exon membership label.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`a_frameshift_intron_is_preclassified_as_coding_without_an_exon_label`](../../crates/fastvep-consequence/src/predictor.rs#L3168).

Related source tests; no isolated before/after replay added for this inventory.

#### S04 — Essential splice-site span overlap

Recognize variants reaching into donor/acceptor dinucleotides, including insertion boundaries.

Source: [splice.rs](../../crates/fastvep-consequence/src/splice.rs).

Related tests: [`a_variant_reaching_onto_the_donor_dinucleotide_is_a_donor_variant`](../../crates/fastvep-consequence/src/splice.rs#L932), [`an_insertion_on_a_splice_boundary_is_in_the_splice_region`](../../crates/fastvep-consequence/src/splice.rs#L977).

Related source tests; no isolated before/after replay added for this inventory.

#### S05 — Differing-base splice predicates

Use differing bases and the applicable differing region rather than padded allele spans.

Source: [splice.rs](../../crates/fastvep-consequence/src/splice.rs).

Related tests: [`only_the_differing_bases_are_matched_against_a_site`](../../crates/fastvep-consequence/src/splice.rs#L1004), [`the_last_differing_region_decides_the_splice_region`](../../crates/fastvep-consequence/src/splice.rs#L1100).

Related source tests; no isolated before/after replay added for this inventory.

#### S06 — Polypyrimidine insertion boundary

Retain an insertion abutting the polypyrimidine tract as an overlap.

Source: [splice.rs](../../crates/fastvep-consequence/src/splice.rs).

Related tests: [`an_insertion_abutting_the_polypyrimidine_tract_is_in_it`](../../crates/fastvep-consequence/src/splice.rs#L1038).

Related source tests; no isolated before/after replay added for this inventory.

#### S07 — Independent extended splice terms

An essential splice-site term must not suppress applicable extended splice terms.

Source: [splice.rs](../../crates/fastvep-consequence/src/splice.rs).

Related tests: [`an_essential_site_does_not_swallow_the_extended_terms`](../../crates/fastvep-consequence/src/splice.rs#L1067).

Related source tests; no isolated before/after replay added for this inventory.

#### S08 — Minimized splice-region coordinates

Remove false splice-region terms introduced by the untrimmed allele representation.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`left_minimal_splice_predicates_remove_a_false_region_term`](../../crates/fastvep-consequence/src/predictor.rs#L3136).

Related source tests; no isolated before/after replay added for this inventory.

### Coding consequences

#### C01 — Independent applicable coding terms

Report all applicable coding predicates rather than selecting one early and suppressing the others.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`every_coding_term_that_holds_is_reported`](../../crates/fastvep-consequence/src/predictor.rs#L2278), [`missense_and_unknown_coding_are_independent`](../../crates/fastvep-consequence/src/predictor.rs#L2797), [`inframe_insertion_and_unknown_coding_are_independent`](../../crates/fastvep-consequence/src/predictor.rs#L2828).

Related source tests; no isolated before/after replay added for this inventory.

#### C02 — Transcript-default fallback

When no coding predicate survives, use the appropriate VEP transcript default, including after NMD classification.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`resolved_repeat_without_a_coding_predicate_defers_to_the_transcript_default`](../../crates/fastvep-consequence/src/predictor.rs#L2571), [`an_empty_predicate_result_uses_veps_default_after_nmd`](../../crates/fastvep-consequence/src/predictor.rs#L3242).

Related source tests; no isolated before/after replay added for this inventory.

#### C03 — Ambiguous alternate translation

Suppress unsupported peptide claims while retaining determinable DNA consequences and start/stop predicates.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`ambiguous_coding_boundary_does_not_assert_loss_or_retention`](../../crates/fastvep-consequence/src/predictor.rs#L2622), [`ambiguous_insertion_can_preserve_a_known_initiator`](../../crates/fastvep-consequence/src/predictor.rs#L2633), [`ambiguous_alternate_has_no_peptide_but_retains_dna_consequences`](../../crates/fastvep-consequence/src/predictor.rs#L2678), [`ambiguous_reference_keeps_dna_stop_alteration_fallback`](../../crates/fastvep-consequence/src/predictor.rs#L3452).

Related source tests; no isolated before/after replay added for this inventory.

#### C04 — Source edits restricted to the reference

Apply only corroborated source edits; do not copy an edited reference residue onto a changed alternate residue.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`source_edited_reference_residue_does_not_rewrite_changed_alt_residue`](../../crates/fastvep-consequence/src/predictor.rs#L2718), [`annotated_residue_resolution_is_limited_to_known_source_edits`](../../crates/fastvep-consequence/src/predictor.rs#L2742), [`unchanged_stop_codon_in_an_alt_window_does_not_inherit_reference_edits`](../../crates/fastvep-consequence/src/predictor.rs#L2752).

Related source tests; no isolated before/after replay added for this inventory.

#### C05 — Mitochondrial reference initiator

Keep the annotated reference initiator adjustment separate from ordinary alternate codon translation.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`an_annotated_mitochondrial_initiator_is_the_reference_residue`](../../crates/fastvep-consequence/src/predictor.rs#L3322).

Verified retained pre-fix fork failure and current all-transcript VEP comparison; see focused review index.

[Focused regression and field examples](focused/mitochondrial-initiator/before-fix.json).

#### C06 — Incomplete CDS start and phase

Do not assert start loss when source completeness flags or phase exclude a known coding start.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`cds_start_nf_suppresses_a_start_lost_claim`](../../crates/fastvep-consequence/src/predictor.rs#L2777), [`a_phase_offset_transcript_does_not_report_start_lost`](../../crates/fastvep-consequence/src/predictor.rs#L4083).

Related source tests; no isolated before/after replay added for this inventory.

#### C07 — Incomplete CDS end

Do not assert stop loss when the source marks the CDS end incomplete.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`a_cds_end_nf_transcript_does_not_report_stop_lost`](../../crates/fastvep-consequence/src/predictor.rs#L3960).

Related source tests; no isolated before/after replay added for this inventory.

#### C08 — Partial-codon endpoints and UTR context

Retain one-ended partial codons and the permitted alternate UTR context on either strand.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`partial_codon_survives_one_unmapped_endpoint_on_either_strand`](../../crates/fastvep-consequence/src/predictor.rs#L2865), [`alternate_partial_codon_can_read_into_the_utr`](../../crates/fastvep-consequence/src/predictor.rs#L2906), [`a_replacement_past_an_incomplete_peptide_uses_its_available_suffix`](../../crates/fastvep-consequence/src/predictor.rs#L4362).

Related source tests; no isolated before/after replay added for this inventory.

#### C09 — Start retention requires valid UTR context

Constrain the UTR-based start-retention shortcut; do not borrow a start from the wrong CDS boundary.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`insertion_retaining_atg_needs_a_utr_for_the_start_shortcut`](../../crates/fastvep-consequence/src/predictor.rs#L3262), [`boundary_start_retention_requires_preserved_utr_and_atg`](../../crates/fastvep-consequence/src/predictor.rs#L3709), [`depleted_cds_cannot_borrow_a_start_codon_from_three_prime_utr`](../../crates/fastvep-consequence/src/predictor.rs#L3741), [`a_c3_deletion_cannot_use_the_c1_utr_repeat_exception`](../../crates/fastvep-consequence/src/predictor.rs#L3776).

Related source tests; no isolated before/after replay added for this inventory.

#### C10 — Independent start-loss and retention predicates

A change can displace and retain the initiator; do not collapse the predicates into one decision.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`non_atg_start_predicates_are_independent_of_peptide_identity`](../../crates/fastvep-consequence/src/predictor.rs#L3282), [`initiator_spanning_deletion_keeps_the_later_start_loss_predicate`](../../crates/fastvep-consequence/src/predictor.rs#L4344), [`an_insertion_can_both_displace_and_retain_the_start_codon`](../../crates/fastvep-consequence/src/predictor.rs#L4388), [`deleting_cds_base_one_has_both_vep_start_terms`](../../crates/fastvep-consequence/src/predictor.rs#L4415).

Related source tests; no isolated before/after replay added for this inventory.

#### C11 — Start-spanning replacements on both strands

Evaluate initiator overlap and the preserved CDS suffix for replacements and intron-spanning edits.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`delins_over_the_start_codon_is_start_lost_on_either_strand`](../../crates/fastvep-consequence/src/predictor.rs#L3686), [`equal_length_boundary_replacement_evaluates_start_independently_of_stop`](../../crates/fastvep-consequence/src/predictor.rs#L3878), [`intron_spanning_replacement_tests_start_at_the_cds_suffix`](../../crates/fastvep-consequence/src/predictor.rs#L3899), [`a_delins_over_the_initiator_is_start_lost`](../../crates/fastvep-consequence/src/predictor.rs#L4514).

Related source tests; no isolated before/after replay added for this inventory.

#### C12 — Start-codon deletion suffix reconstruction

Recognize retained CDS suffixes and the next-base reconstruction of a partially deleted initiator.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`deleting_c3_can_retain_the_initiator_using_the_next_base`](../../crates/fastvep-consequence/src/predictor.rs#L3760), [`boundary_edit_can_preserve_the_entire_cds_suffix`](../../crates/fastvep-consequence/src/predictor.rs#L3837), [`retained_deletion_at_start_can_preserve_the_cds_suffix`](../../crates/fastvep-consequence/src/predictor.rs#L3859).

Related source tests; no isolated before/after replay added for this inventory.

#### C13 — Split start-codon insertions

Use the surviving mapped endpoint when an insertion splits a start codon across exons.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`insertion_in_a_split_start_codon_uses_the_surviving_endpoint`](../../crates/fastvep-consequence/src/predictor.rs#L3806).

Related source tests; no isolated before/after replay added for this inventory.

#### C14 — Independent stop-boundary predicates

Distinguish stop loss, retention and generic coding changes across the stop/UTR boundary.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`delins_over_the_stop_codon_is_stop_lost_on_either_strand`](../../crates/fastvep-consequence/src/predictor.rs#L3584), [`boundary_predicates_use_the_full_deleted_prefix_and_partial_codon_guard`](../../crates/fastvep-consequence/src/predictor.rs#L3607), [`a_same_length_change_across_the_stop_boundary_is_generic_coding`](../../crates/fastvep-consequence/src/predictor.rs#L3629), [`a_deletion_past_the_stop_codon_distinguishes_lost_from_retained`](../../crates/fastvep-consequence/src/predictor.rs#L3654).

Related source tests; no isolated before/after replay added for this inventory.

#### C15 — Full-reference protein stop retention

Compare against the full annotated reference protein, including terminal-repeat insertions.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`stop_retention_compares_the_full_annotated_reference_protein`](../../crates/fastvep-consequence/src/predictor.rs#L2147), [`a_preserved_first_residue_and_later_stop_matches_vep`](../../crates/fastvep-consequence/src/predictor.rs#L4283), [`an_insertion_in_a_terminal_repeat_can_retain_the_reference_protein`](../../crates/fastvep-consequence/src/predictor.rs#L4307).

Related source tests; no isolated before/after replay added for this inventory.

#### C16 — Protein-altering versus in-frame replacement

Distinguish replaced residues from extension of the retained reference peptide on both strands.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`a_delins_that_replaces_residues_is_protein_altering_on_either_strand`](../../crates/fastvep-consequence/src/predictor.rs#L4225), [`a_delins_that_extends_the_reference_residues_is_an_inframe_insertion`](../../crates/fastvep-consequence/src/predictor.rs#L4255).

Related source tests; no isolated before/after replay added for this inventory.

#### C17 — Terminator introduction combinations

Preserve independent applicable stop/in-frame terms, including a peptide beginning with a terminator.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`a_delins_introducing_a_terminator_reports_both_terms`](../../crates/fastvep-consequence/src/predictor.rs#L4463), [`a_delins_whose_peptide_begins_with_a_terminator_is_stop_gained_alone`](../../crates/fastvep-consequence/src/predictor.rs#L4487).

Related source tests; no isolated before/after replay added for this inventory.

#### C18 — Codon reconstruction and capitalization

Build the edited codon consistently and display only the minimized codon for padded repeats.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`edited_codon_agrees_with_building_the_edited_sequence`](../../crates/fastvep-consequence/src/predictor.rs#L4028), [`padded_repeat_deletion_displays_only_the_minimized_codon`](../../crates/fastvep-consequence/src/predictor.rs#L4190).

Related source tests; no isolated before/after replay added for this inventory.

#### C19 — Discontinuous CDS and exon-edge codon edits

Do not name unsupported residues for discontinuous CDS changes; preserve valid exon-edge edits.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`a_change_that_is_not_contiguous_in_the_cds_names_no_residues`](../../crates/fastvep-consequence/src/predictor.rs#L4537), [`a_codon_aligned_insertion_names_no_reference_residue`](../../crates/fastvep-consequence/src/predictor.rs#L4571), [`a_partial_codon_after_a_terminator_adds_no_placeholder`](../../crates/fastvep-consequence/src/predictor.rs#L4609), [`an_insertion_on_an_exon_edge_is_still_a_codon_edit`](../../crates/fastvep-consequence/src/predictor.rs#L4631).

Related source tests; no isolated before/after replay added for this inventory.

#### C20 — No extra generic region term

Changes wholly inside one region must not acquire an additional unsupported regional consequence.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs).

Related tests: [`variants_inside_a_single_region_gain_no_extra_term`](../../crates/fastvep-consequence/src/predictor.rs#L3927).

Related source tests; no isolated before/after replay added for this inventory.

#### C21 — Whole-transcript overlap precedence

Apply transcript ablation/amplification before lower consequence tiers, retain independently calculated positions, and exclude inapplicable generic exon/CDS terms on complete overlap.

Source: [predictor.rs](../../crates/fastvep-consequence/src/predictor.rs#L706).

Implementation inspection; no isolated fixture or direct unit-test mapping recorded here.

### Transcript HGVS

#### H01 — Genomic shifting before transcript mapping

Allow legitimate exon/intron transitions without jumping introns through adjacent spliced bases.

Source: [hgvs_normalize.rs](../../crates/fastvep-annotate/src/hgvs_normalize.rs).

Related tests: [`exonic_deletion_can_shift_to_the_first_intronic_base`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1093), [`a_genomic_deletion_shift_does_not_jump_an_intron`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1159), [`an_unmapped_deletion_is_not_promoted_into_a_later_exon`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1199), [`deletion_can_shift_from_an_intron_into_an_exon`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1239).

Related source tests; no isolated before/after replay added for this inventory.

#### H02 — Reverse shifting without persistent cDNA cache

Build temporary transcript sequence when needed and keep it out of the stored transcript cache.

Source: [hgvs_normalize.rs](../../crates/fastvep-annotate/src/hgvs_normalize.rs).

Related tests: [`transient_spliced_sequence_is_not_cached_on_the_transcript`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L945), [`reverse_exonic_deletion_shifts_without_a_cached_transcript_sequence`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1224).

Related source tests; no isolated before/after replay added for this inventory.

#### H03 — Transcript and mitochondrial shift boundaries

Enforce full transcript-slice eligibility, including clipped mitochondrial windows and long alleles.

Source: [hgvs_normalize.rs](../../crates/fastvep-annotate/src/hgvs_normalize.rs).

Related tests: [`vep_boundary_check_allows_a_terminal_insertion_but_not_a_shift_past_it`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L969), [`mitochondrial_shift_handles_alleles_longer_than_the_search_flank`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1015), [`vep_boundary_check_matches_linear_and_clipped_mitochondrial_windows`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1035).

Related source tests; no isolated before/after replay added for this inventory.

#### H04 — Rotate shifted reference and inserted sequence

Retain allele rotation, including ambiguity, across shifting and exon/intron transitions.

Source: [hgvs_normalize.rs](../../crates/fastvep-annotate/src/hgvs_normalize.rs).

Related tests: [`deletion_shift_rotates_input_reference_including_ambiguity`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L773), [`exonic_duplication_retains_rotation_when_the_shift_crosses_an_intron`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1133), [`dup_span_compares_against_the_rotated_insert`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1281).

Related source tests; no isolated before/after replay added for this inventory.

#### H05 — Duplication at the shifted insertion

Choose duplication spans from the shifted insertion and strand-specific reference side.

Source: [hgvs_normalize.rs](../../crates/fastvep-annotate/src/hgvs_normalize.rs).

Related tests: [`exonic_insertion_can_shift_to_an_intronic_duplication`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1178), [`dup_span_follows_the_shifted_insertion_not_the_vcf_position`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1269), [`dup_span_reads_the_other_side_on_the_reverse_strand`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1293), [`intronic_insertion_becomes_a_dup_at_its_shifted_position`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1331).

Related source tests; no isolated before/after replay added for this inventory.

#### H06 — Duplication boundary and representation checks

Reject nonrepeating/out-of-contig spans and retain a single coherent duplicated block across representations.

Source: [hgvs_normalize.rs](../../crates/fastvep-annotate/src/hgvs_normalize.rs).

Related tests: [`dup_span_is_none_when_the_block_does_not_repeat`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1302), [`dup_span_declines_rather_than_reading_past_the_contig`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1312), [`a_block_leaving_the_intron_is_mapped_as_a_duplication`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1353), [`every_spelling_of_one_insertion_names_the_same_block`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1381), [`a_span_crossing_the_intron_midpoint_is_still_one_dup`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L1403).

Related source tests; no isolated before/after replay added for this inventory.

#### H07 — Retained replacements do not take pure-indel shifts

Keep replacement classification and clipping separate from standalone indel normalization.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs).

Related tests: [`retained_replacement_hgvs_clips_without_a_standalone_indel_shift`](../../crates/fastvep-annotate/src/lib.rs#L2320).

Related source tests; no isolated before/after replay added for this inventory.

#### H08 — Inversion before clipping

Classify the complete replacement before clipping; clipping must not promote a delins into an inversion.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs#L1809).

Related tests: [`a_replacement_by_the_reverse_complement_is_an_inversion`](../../crates/fastvep-hgvs/src/coding.rs#L1056).

Verified retained pre-fix fork failure and current all-transcript VEP comparison; see focused review index.

[Focused regression and field examples](focused/inversion-before-clipping/before-fix.json).

The coding test covers inversion recognition; the focused fixture verifies classification before clipping in the annotation path.

#### H09 — CDS/UTR range numbering

Keep both endpoints across the initiator and terminator without c.0 or star-zero coordinates.

Source: [coding.rs](../../crates/fastvep-hgvs/src/coding.rs).

Related tests: [`a_delins_reaching_past_the_terminator_numbers_its_end_from_the_stop`](../../crates/fastvep-hgvs/src/coding.rs#L853), [`a_delins_reaching_out_of_the_five_prime_utr_keeps_both_coordinates`](../../crates/fastvep-hgvs/src/coding.rs#L874), [`a_deletion_from_the_five_prime_utr_into_the_cds_numbers_its_end_in_the_cds`](../../crates/fastvep-hgvs/src/coding.rs#L892), [`a_span_across_the_initiator_skips_the_nonexistent_position_zero`](../../crates/fastvep-hgvs/src/coding.rs#L909), [`the_last_cds_base_is_never_numbered_star_zero`](../../crates/fastvep-hgvs/src/coding.rs#L926).

Related source tests; no isolated before/after replay added for this inventory.

#### H10 — Terminal intronic and duplication anchors

Use the applicable stop-relative anchor at terminal coding positions.

Source: [coding.rs](../../crates/fastvep-hgvs/src/coding.rs).

Related tests: [`intron_at_terminal_coding_base_uses_vep_star_anchor`](../../crates/fastvep-hgvs/src/coding.rs#L774).

Related source tests; no isolated before/after replay added for this inventory.

#### H11 — Terminal duplication star notation

Use the terminal coding anchor convention for a shifted duplication.

Source: [hgvs_normalize.rs](../../crates/fastvep-annotate/src/hgvs_normalize.rs).

Related tests: [`duplication_at_a_terminal_coding_anchor_uses_vep_star_notation`](../../crates/fastvep-annotate/src/hgvs_normalize.rs#L954).

Related source tests; no isolated before/after replay added for this inventory.

#### H12 — Phase offset only for substitutions

Do not apply the CDS substitution phase adjustment to other allele classes.

Source: [coding.rs](../../crates/fastvep-hgvs/src/coding.rs).

Related tests: [`only_a_substitution_carries_the_cds_phase_offset`](../../crates/fastvep-hgvs/src/coding.rs#L954).

Related source tests; no isolated before/after replay added for this inventory.

#### H13 — Intronic ranges in transcript order

Render multibase intronic changes as ordered transcript-relative ranges.

Source: [coding.rs](../../crates/fastvep-hgvs/src/coding.rs).

Related tests: [`a_multi_base_intronic_change_is_written_as_a_range`](../../crates/fastvep-hgvs/src/coding.rs#L1095), [`an_intronic_span_is_written_in_transcript_order`](../../crates/fastvep-hgvs/src/coding.rs#L1122).

Related source tests; no isolated before/after replay added for this inventory.

#### H14 — Shift before insertion/duplication choice

Apply the transcript-oriented repeat shift before choosing the duplication representation.

Source: [coding.rs](../../crates/fastvep-hgvs/src/coding.rs).

Related tests: [`an_insertion_in_a_repeat_shifts_before_it_can_be_a_duplication`](../../crates/fastvep-hgvs/src/coding.rs#L1018).

Related source tests; no isolated before/after replay added for this inventory.

### Protein HGVS

#### P01 — Separate consequence and full-reference peptides

Keep raw consequence windows distinct from the full reference used for flank and duplication context.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`insertion_clips_raw_peptide_but_names_full_reference_flanks`](../../crates/fastvep-hgvs/src/protein.rs#L1565), [`initiation_edits_do_not_replace_the_duplication_reference`](../../crates/fastvep-hgvs/src/protein.rs#L1611), [`frameshift_uses_annotated_reference_peptide_and_vep_terminal_deletions`](../../crates/fastvep-hgvs/src/protein.rs#L1835).

Related source tests; no isolated before/after replay added for this inventory.

#### P02 — Recreated-stop clipping order

Return the original window when prefix scanning reaches a recreated stop, before committing a prefix trim.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`recreated_stop_keeps_the_original_window_before_clipping`](../../crates/fastvep-hgvs/src/protein.rs#L1667).

Verified retained pre-fix fork failure and current all-transcript VEP comparison; see focused review index.

[Focused regression and field examples](focused/recreated-stop-clipping/before-fix.json).

#### P03 — Start-loss uncertainty and peptide preparation

Use the explicit start-loss predicate and prepared peptide window; decline unsupported initiator descriptions.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`start_loss_keeps_the_clipped_partial_codon_window`](../../crates/fastvep-hgvs/src/protein.rs#L1550), [`shifted_frameshift_start_loss_uses_the_changed_full_peptide_residue`](../../crates/fastvep-hgvs/src/protein.rs#L2070), [`start_loss_uses_vep_peptide_preparation_and_explicit_predicate`](../../crates/fastvep-hgvs/src/protein.rs#L2471), [`a_deletion_of_the_initiation_codon_is_described_as_unresolvable`](../../crates/fastvep-hgvs/src/protein.rs#L2498), [`a_change_that_loses_the_initiator_is_unresolvable`](../../crates/fastvep-hgvs/src/protein.rs#L2518).

Related source tests; no isolated before/after replay added for this inventory.

#### P04 — Full-reference flanks after clipped start loss

Use the full-reference peptide for flank naming after start-loss clipping.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs).

Related tests: [`clipped_start_loss_uses_full_reference_flanks`](../../crates/fastvep-annotate/src/lib.rs#L2285).

Related source tests; no isolated before/after replay added for this inventory.

#### P05 — Shifted CDS endpoint eligibility

Keep shifted translation endpoints inside the CDS and preserve VEP's asymmetric boundary fallback.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs).

Related tests: [`protein_hgvs_follows_veps_asymmetric_boundary_shift_fallback`](../../crates/fastvep-annotate/src/lib.rs#L2621), [`shifted_translation_endpoint_must_remain_inside_the_cds`](../../crates/fastvep-annotate/src/lib.rs#L2631).

Related source tests; no isolated before/after replay added for this inventory.

#### P06 — Recover frameshift from mapped CDS coordinates

Reconstruct the protein window when both shifted CDS coordinates are available.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs).

Related tests: [`hgvsp_reconstructs_a_frameshift_when_shifted_cds_coordinates_are_complete`](../../crates/fastvep-annotate/src/lib.rs#L3357).

Related source tests; no isolated before/after replay added for this inventory.

#### P07 — Internal exon-edge insertion eligibility

Retain applicable HGVSp at an internal exon boundary.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs).

Related tests: [`internal_exon_edge_insertion_retains_protein_hgvs`](../../crates/fastvep-annotate/src/lib.rs#L2518).

Related source tests; no isolated before/after replay added for this inventory.

#### P08 — Sort insertion flanks before stretched overlap

Order the insertion endpoints before the short-intron exon-overlap gate.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs#L1659).

Related tests: [`vep_frameshift_intron_definition_includes_thirteen_bases`](../../crates/fastvep-annotate/src/lib.rs#L2647).

Verified retained pre-fix fork failure and current all-transcript VEP comparison; see focused review index.

[Focused regression and field examples](focused/insertion-flank-order/before-fix.json).

The unit test covers the intron threshold, not endpoint sorting. The focused fixture verifies endpoint sorting.

#### P09 — Keep terminal insertion flanks without genomic shift

Do not require a nonzero genomic shift to retain valid terminal protein flanks.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs).

Related tests: [`hgvsp_terminal_insertion_without_genomic_shift_keeps_flanks`](../../crates/fastvep-annotate/src/lib.rs#L3518).

Related source tests; no isolated before/after replay added for this inventory.

#### P10 — Hidden consequence predicates for protein formatting

Retain applicable internal-stop and stop-retention information even when the emitted consequence term is gated.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs).

Related tests: [`protein_formatting_keeps_hidden_stop_predicates_and_internal_stops`](../../crates/fastvep-annotate/src/lib.rs#L2447), [`stop_loss_after_an_unchanged_residue_finds_the_extension_stop`](../../crates/fastvep-annotate/src/lib.rs#L2480).

Related source tests; no isolated before/after replay added for this inventory.

#### P11 — Depleted CDS versus UTR context

Trim depleted CDS before alternate UTR translation and keep short-UTR and partial-reference windows distinct.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`depleted_cds_keeps_short_utr_and_partial_reference_windows_distinct`](../../crates/fastvep-hgvs/src/protein.rs#L1517), [`depleted_cds_is_trimmed_before_utr_translation`](../../crates/fastvep-hgvs/src/protein.rs#L1537).

Related source tests; no isolated before/after replay added for this inventory.

#### P12 — Reference partial codon stays inside CDS

Do not borrow UTR bases to complete a reference terminal codon.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`reference_terminal_codon_does_not_borrow_utr_bases`](../../crates/fastvep-hgvs/src/protein.rs#L1590).

Related source tests; no isolated before/after replay added for this inventory.

#### P13 — Partial terminal insertion reconstruction

Recompute shifted terminal codons and allow duplication of a supported partial reference residue.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`stop_insertion_can_duplicate_a_partial_reference_residue`](../../crates/fastvep-hgvs/src/protein.rs#L1603), [`shifted_insertion_completes_a_partial_terminal_codon`](../../crates/fastvep-hgvs/src/protein.rs#L1677), [`shifted_insertion_recomputes_the_terminal_codon_window`](../../crates/fastvep-hgvs/src/protein.rs#L1875).

Related source tests; no isolated before/after replay added for this inventory.

#### P14 — Terminal post-sequence eligibility

A terminal insertion flank needs an actual following residue; preserve explicit endpoint ranges.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`terminal_insertion_post_sequence_requires_a_following_residue`](../../crates/fastvep-hgvs/src/protein.rs#L1623), [`partial_stop_window_keeps_both_translation_endpoints`](../../crates/fastvep-hgvs/src/protein.rs#L1657).

Related source tests; no isolated before/after replay added for this inventory.

#### P15 — Terminator and unknown-residue clipping

Keep terminators and unknown residues through the applicable clipping stage before trimming the output.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`terminal_windows_follow_vep_clipping_and_stop_loss_formatting`](../../crates/fastvep-hgvs/src/protein.rs#L1634), [`test_hgvsp_uses_vep_ter_spelling_for_x`](../../crates/fastvep-hgvs/src/protein.rs#L2059), [`protein_terminator_is_trimmed_after_clipping_like_vep`](../../crates/fastvep-hgvs/src/protein.rs#L2109).

Related source tests; no isolated before/after replay added for this inventory.

#### P16 — Terminal frameshift windows and shifted starts

Preserve the initial clipped window and use the shifted start for terminal frameshift deletions.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`terminal_frameshift_keeps_the_initial_clipped_window`](../../crates/fastvep-hgvs/src/protein.rs#L1700), [`terminal_frameshift_deletion_uses_the_shifted_start`](../../crates/fastvep-hgvs/src/protein.rs#L1716), [`vep_retains_the_last_equal_residue_when_the_alternate_translation_ends`](../../crates/fastvep-hgvs/src/protein.rs#L1773).

Related source tests; no isolated before/after replay added for this inventory.

#### P17 — Frameshift stop-distance calculation

Count determinable stops, including upstream alternate stops and stops beyond the annotated terminator.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`residue_one_frameshift_keeps_a_determinable_stop_distance`](../../crates/fastvep-hgvs/src/protein.rs#L1744), [`vep_uses_an_upstream_alternate_stop_for_frameshift_distance`](../../crates/fastvep-hgvs/src/protein.rs#L1756), [`a_frameshift_stop_past_the_annotated_terminator_is_counted_to`](../../crates/fastvep-hgvs/src/protein.rs#L3065), [`a_frameshift_landing_on_a_terminator_is_written_as_nonsense`](../../crates/fastvep-hgvs/src/protein.rs#L3091).

Related source tests; no isolated before/after replay added for this inventory.

#### P18 — Mitochondrial frameshift translation context

Preserve the VEP-specific distinction between reference codon translation and frameshift alternate translation.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`test_hgvsp_frameshift_mitochondrial_table_differs`](../../crates/fastvep-hgvs/src/protein.rs#L1790).

Related source tests; no isolated before/after replay added for this inventory.

#### P19 — Rotate protein indels before tail translation

Shift local residue windows without replacing the initiator, then translate the appropriate peptide tails.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`deletion_rotates_local_residues_without_replacing_the_initiator`](../../crates/fastvep-hgvs/src/protein.rs#L1573), [`shifted_inframe_deletion_is_clipped_from_translated_peptide_tails`](../../crates/fastvep-hgvs/src/protein.rs#L2154), [`shifted_inframe_insertion_rotates_before_translating_the_peptide_tail`](../../crates/fastvep-hgvs/src/protein.rs#L2174).

Related source tests; no isolated before/after replay added for this inventory.

#### P20 — Terminal deletion shift bounds

Do not start or end an in-frame deletion shift at an unsupported terminal residue.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`deletion_shift_does_not_start_at_the_final_reference_residue`](../../crates/fastvep-hgvs/src/protein.rs#L1690), [`test_hgvsp_inframe_deletion_is_three_prime_shifted`](../../crates/fastvep-hgvs/src/protein.rs#L2275), [`test_hgvsp_inframe_deletion_uses_vep_115_terminal_shift_bound`](../../crates/fastvep-hgvs/src/protein.rs#L2284), [`test_hgvsp_inframe_indel_does_not_shift_onto_the_terminator`](../../crates/fastvep-hgvs/src/protein.rs#L2395).

Related source tests; no isolated before/after replay added for this inventory.

#### P21 — Insertion versus single/multiple duplication

Distinguish a genuine insertion from a repeated single- or multi-residue block.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`test_hgvsp_inframe_insertion_collapses_to_duplication`](../../crates/fastvep-hgvs/src/protein.rs#L2222), [`test_hgvsp_inframe_insertion_multi_residue_duplication`](../../crates/fastvep-hgvs/src/protein.rs#L2239), [`test_hgvsp_inframe_insertion_true_insertion_uses_ins_form`](../../crates/fastvep-hgvs/src/protein.rs#L2258).

Related source tests; no isolated before/after replay added for this inventory.

#### P22 — Fallback when peptide context is unavailable

Keep valid unshifted descriptions where supported and decline contexts contradicted by both endpoints.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`test_hgvsp_inframe_indel_without_peptide_stays_valid`](../../crates/fastvep-hgvs/src/protein.rs#L2294), [`test_hgvsp_inframe_indel_survives_unusable_peptides`](../../crates/fastvep-hgvs/src/protein.rs#L2317), [`test_hgvsp_inframe_indel_ignores_a_peptide_that_disagrees`](../../crates/fastvep-hgvs/src/protein.rs#L2417), [`test_hgvsp_inframe_indel_still_declines_when_neither_end_corroborates`](../../crates/fastvep-hgvs/src/protein.rs#L2565).

Related source tests; no isolated before/after replay added for this inventory.

#### P23 — Anchor choice for reverse and periodic sequences

Use determined endpoints and strand/shape-aware anchors, preserving the caller's intended peptide span.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`explicit_peptide_range_keeps_its_start_in_a_reverse_repeat`](../../crates/fastvep-hgvs/src/protein.rs#L1629), [`test_hgvsp_inframe_indel_reads_the_span_from_either_end`](../../crates/fastvep-hgvs/src/protein.rs#L2427), [`anchor_candidates_offers_the_other_end_only_when_there_is_one`](../../crates/fastvep-hgvs/src/protein.rs#L2587), [`anchor_candidates_puts_the_determined_end_first`](../../crates/fastvep-hgvs/src/protein.rs#L2600), [`only_a_shrinking_change_on_the_reverse_strand_is_anchored_at_its_end`](../../crates/fastvep-hgvs/src/protein.rs#L2615), [`a_periodic_reference_on_the_reverse_strand_names_the_span_the_caller_meant`](../../crates/fastvep-hgvs/src/protein.rs#L2638), [`a_periodic_reference_on_the_forward_strand_is_read_from_the_start`](../../crates/fastvep-hgvs/src/protein.rs#L2675).

Related source tests; no isolated before/after replay added for this inventory.

#### P24 — Description shape and reconstruction checks

Keep substitution versus delins shape valid and reconstruct the intended altered protein in repeat controls.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`every_description_reconstructs_the_protein_the_variant_produces`](../../crates/fastvep-hgvs/src/protein.rs#L2783), [`a_two_residue_homopolymer_reads_the_same_from_either_end`](../../crates/fastvep-hgvs/src/protein.rs#L2852), [`test_hgvsp_inframe_indel_never_emits_substitution_shape`](../../crates/fastvep-hgvs/src/protein.rs#L2904), [`a_two_residue_window_that_changes_one_residue_is_a_substitution`](../../crates/fastvep-hgvs/src/protein.rs#L3028).

Related source tests; no isolated before/after replay added for this inventory.

#### P25 — Terminal insertion suppression

Decline the terminal insertion shapes VEP suppresses instead of inventing unsupported protein flanks.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`test_hgvsp_inframe_indel_matches_vep_terminal_insertion_suppression`](../../crates/fastvep-hgvs/src/protein.rs#L2881).

Related source tests; no isolated before/after replay added for this inventory.

#### P26 — Edited-CDS stop-loss extension

Apply replacement sequence on either strand and count to the next translated stop.

Source: [protein.rs](../../crates/fastvep-hgvs/src/protein.rs).

Related tests: [`the_edited_cds_replaces_the_reference_bases_on_either_strand`](../../crates/fastvep-hgvs/src/protein.rs#L3133), [`stop_loss_counts_to_the_next_stop_on_either_strand`](../../crates/fastvep-hgvs/src/protein.rs#L3166).

Related source tests; no isolated before/after replay added for this inventory.

#### P27 — Inconsistent spliced-sequence fallback

Decline frameshift reconstruction when the supplied spliced sequence length contradicts its model.

Source: [lib.rs](../../crates/fastvep-annotate/src/lib.rs#L3197).

Implementation inspection; no isolated fixture or direct unit-test mapping recorded here.

#### P28 — Selenocysteine and pyrrolysine residue names

Render U and O as Sec and Pyl rather than unknown residue names; this naming change does not infer a source sequence edit.

Source: [codon.rs](../../crates/fastvep-genome/src/codon.rs#L197).

Related tests: [`test_aa_one_to_three`](../../crates/fastvep-genome/src/codon.rs#L300).

Related source tests; no isolated before/after replay added for this inventory.

### Transcript data and cache correctness

#### K01 — Source-declared peptide edits

Build consequence peptides from codons and explicit source edits rather than inferred amino-acid corrections.

Source: [ensembl_core.rs](../../crates/fastvep-cache/src/ensembl_core.rs).

Related tests: [`codon_peptide_uses_only_source_declared_edits`](../../crates/fastvep-cache/src/ensembl_core.rs#L540).

Related source tests; no isolated before/after replay added for this inventory.

#### K02 — Corroborated selenocysteine only

Resolve internal TGA as selenocysteine only with corroborating source information; leave mitochondrial TGA unchanged.

Source: [transcript.rs](../../crates/fastvep-genome/src/transcript.rs).

Related tests: [`resolves_only_corroborated_internal_tga_as_selenocysteine`](../../crates/fastvep-genome/src/transcript.rs#L454), [`leaves_mitochondrial_tga_translation_unchanged`](../../crates/fastvep-genome/src/transcript.rs#L462).

Related source tests; no isolated before/after replay added for this inventory.

#### K03 — Reference prefetch interval identity

Reuse a prefetched sequence only for the same complete interval.

Source: [providers.rs](../../crates/fastvep-cache/src/providers.rs).

Related tests: [`prefetched_reference_reuses_only_the_matching_complete_interval`](../../crates/fastvep-cache/src/providers.rs#L384).

Related source tests; no isolated before/after replay added for this inventory.

#### K04 — Legacy and enriched cache decoding

Retain legacy cache readability and the enriched payload's ranges; reject damaged payloads.

Source: [transcript_cache.rs](../../crates/fastvep-cache/src/transcript_cache.rs).

Related tests: [`pristine_and_patched_legacy_layouts_remain_readable`](../../crates/fastvep-cache/src/transcript_cache.rs#L410), [`annocat_roundtrip_preserves_ranges_and_rejects_damaged_payloads`](../../crates/fastvep-cache/src/transcript_cache.rs#L436), [`test_legacy_gzip_cache_loads`](../../crates/fastvep-cache/src/transcript_cache.rs#L515).

Related source tests; no isolated before/after replay added for this inventory.

#### K05 — Primary coding-sequence completeness

Require primary-contig coding sequences while permitting documented non-primary omissions.

Source: [transcript_cache.rs](../../crates/fastvep-cache/src/transcript_cache.rs).

Related tests: [`verifies_a_sequence_complete_primary_cache`](../../crates/fastvep-cache/src/transcript_cache.rs#L559), [`rejects_missing_primary_coding_sequences_when_required`](../../crates/fastvep-cache/src/transcript_cache.rs#L572), [`permits_missing_non_primary_sequences`](../../crates/fastvep-cache/src/transcript_cache.rs#L583).

Related source tests; no isolated before/after replay added for this inventory.

#### K06 — Atomic cache replacement

Do not leave a partial cache when replacing an existing one.

Source: [transcript_cache.rs](../../crates/fastvep-cache/src/transcript_cache.rs).

Related tests: [`replaces_an_existing_cache_without_leaving_a_partial_file`](../../crates/fastvep-cache/src/transcript_cache.rs#L593).

Related source tests; no isolated before/after replay added for this inventory.

#### K07 — Public Core identity and version checks

Reject inconsistent GFF/Core transcript and gene identities; preserve gene and protein versions.

Source: [ensembl_core.rs](../../crates/fastvep-cache/src/ensembl_core.rs#L185).

Implementation inspection; no isolated fixture or direct unit-test mapping recorded here.

#### K08 — Source membership exclusions

Match source-declared artifact and readthrough exclusions during cache construction. This is not a blanket exclusion inferred from transcript names.

Source: [ensembl_core.rs](../../crates/fastvep-cache/src/ensembl_core.rs#L397).

Implementation inspection; no isolated fixture or direct unit-test mapping recorded here.

#### K09 — Public metadata enrichment

Populate canonical, display cross-references, HGNC, CCDS, APPRIS, TSL, MANE, GENCODE and completeness flags from Core attributes.

Source: [ensembl_core.rs](../../crates/fastvep-cache/src/ensembl_core.rs#L410).

Implementation inspection; no isolated fixture or direct unit-test mapping recorded here.

#### K10 — Mature-miRNA and explicit translation edits

Retain mature-miRNA ranges and the four source translation-edit codes. Reject RNA coordinate edits until an explicit mapper exists.

Source: [ensembl_core.rs](../../crates/fastvep-cache/src/ensembl_core.rs#L451).

Implementation inspection; no isolated fixture or direct unit-test mapping recorded here.

#### K11 — Enriched cache validation and provenance

Validate enrichment and retain a semantic digest alongside the enriched payload. Cache-format support is supporting infrastructure, not a count of VEP bugs.

Source: [annocat_cache.rs](../../crates/fastvep-cache/src/annocat_cache.rs#L88).

Implementation inspection; no isolated fixture or direct unit-test mapping recorded here.

#### K12 — Indexed reference provider integration

Use the existing indexed reference reader during cache building instead of loading the complete FASTA. Prefetch interval correctness is covered separately by K03. This is caller integration, not a newly invented reader or cache format.

Source: [pipeline.rs](../../crates/fastvep-cli/src/pipeline.rs#L2342).

Implementation inspection; no isolated fixture or direct unit-test mapping recorded here.

#### K13 — Cache build isolation and read-only annotation

Serialize concurrent builds for an output, preserve explicitly supplied caches during annotation and prevent regional selections from being published as complete reusable caches.

Source: [pipeline.rs](../../crates/fastvep-cli/src/pipeline.rs#L2466).

Related tests: [`rejects_a_second_builder_for_the_same_output`](../../crates/fastvep-cli/src/pipeline.rs#L2501), [`explicit_transcript_cache_is_read_only_during_annotation`](../../crates/fastvep-cli/src/pipeline.rs#L2512), [`region_restricted_transcripts_are_not_published_as_a_reusable_cache`](../../crates/fastvep-cli/src/pipeline.rs#L2519).

Related source tests; no isolated before/after replay added for this inventory.

### Output consistency

#### O01 — Intergenic CSQ identity

Use VEP's empty transcript identity in intergenic VCF rows without changing the structured JSON model.

Source: [output.rs](../../crates/fastvep-io/src/output.rs).

Related tests: [`intergenic_vcf_uses_vep_empty_feature_identity_without_changing_json`](../../crates/fastvep-io/src/output.rs#L2515).

Related source tests; no isolated before/after replay added for this inventory.

#### O02 — CSQ column integrity

Resolve each default CSQ column once and preserve the declared field count and formatting contract.

Source: [output.rs](../../crates/fastvep-io/src/output.rs).

Related tests: [`every_default_csq_column_resolves_to_a_distinct_writer`](../../crates/fastvep-io/src/output.rs#L2228).

Related source tests; no isolated before/after replay added for this inventory.

#### O03 — Remove stale owned INFO fields

Remove declared fastVEP-owned INFO fields even when the new annotation has no replacement value; ownership no longer depends only on fields projected for this record.

Source: [output.rs](../../crates/fastvep-io/src/output.rs).

Related tests: [`vcf_info_replaces_existing_fastvep_owned_fields`](../../crates/fastvep-io/src/output.rs#L2654).

Related source tests; no isolated before/after replay added for this inventory.

#### O05 — Shared CLI and library annotation paths

Route CLI and library annotation through the same record-aware HGVSc/HGVSp and position helpers; carry normalized allele and hidden frameshift state through the result types. This is supporting integration, not a separate biological algorithm.

Source: [pipeline.rs](../../crates/fastvep-cli/src/pipeline.rs#L376).

Implementation inspection; no isolated fixture or direct unit-test mapping recorded here.

#### O06 — Explicit transcript-selection and range inputs

Validate explicitly requested transcript exclusions without changing the cache and accept source mature-miRNA ranges. The all-transcript comparison profile does not enable user exclusions.

Source: [pipeline.rs](../../crates/fastvep-cli/src/pipeline.rs#L268).

Related tests: [`transcript_exclusion_is_validated_and_does_not_modify_the_cache`](../../crates/fastvep-cli/src/pipeline.rs#L6190).

Related source tests; no isolated before/after replay added for this inventory.

## Deliberate differences and comparison policy

Known exceptions remain separate from repairs intended to reproduce VEP. The established policies preserve correct CDS-to-UTR coordinate order and omit unsupported HGVSc beyond the terminal exon. Native-only controls retain records and complete where VEP fails; an oracle failure cannot establish output equality. Exact applicability is defined by each comparison contract, not by a blanket exemption for an entire variant class.

FLAGS ordering is a declared normalizer in applicable contracts. Normalized or exception-qualified equality must not be described as raw byte equality. See the contracts distributed with the full suite for all accepted differences.

## Work outside this inventory

Supplementary-source loading and lookup, ACMG/QC settings, allocation reductions, parallel execution, CLI/library integration and direct-Parquet integration are separate changes. Their tests may appear in the source catalog because the same files changed. They are not included as additional VEP-concordance fixes. Structural-variant, RefSeq and other profiles require separate qualification.

The catalog also retains ordinary baseline consequence/HGVS tests, GFF parsing and cache round-trip controls, and diagnostic inventory utilities. Their presence does not establish an additional concordance repair. A future upstream patch submission should select a family, inspect its implementation delta and reduce its regression to that patch; the bundled history does not supply one clean commit per family.
