# Reviewing the VEP 115 fixes

The [concordance change inventory](FIX-INVENTORY.md) covers 113 review families,
including supporting transcript-cache and output changes. The four groups below
are selected examples with verified before/after results, not the full fix list.

Start with these four regression groups: **13 input records**, with all their
VEP transcript annotations. Each group reproduces a failure in a retained
pre-fix fork build and passes with the current fork against the frozen VEP
115.2 oracle. No exception entries are needed for these cases.

These are selected regressions, not an independent accuracy sample. The earlier
builds are fork builds; this does not establish that current upstream still has
each defect.

## Run the focused suite

From the fastVEP repository root, with Python 3.11 or later:

```sh
python scripts/concordance/run.py --data tests/concordance/focused --binary ./fastvep --cache ./transcripts.cache --fasta ./reference.fa --output ./focused-results
```

Supply a built executable, the matching GRCh38 transcript cache and an indexed
reference FASTA. The small fixture includes VEP outputs, so a VEP installation
and the comprehensive oracle download are unnecessary for this check. The
reference and transcript cache are not bundled. Their tested fingerprints and
the annotation profile are recorded in the manifest and validation record.

The command checks all transcript identities, record preservation and declared
fields, including existing normalizers such as FLAGS order. It rejects failed
comparisons. The output directory must be new.

## Fix-to-test index

| Fix and trigger | Previous behavior | Corrected behavior | Implementation | Fixture |
| --- | --- | --- | --- | --- |
| Mitochondrial initiator context | Ordinary codon translation was used where the annotated reference initiator required methionine; amino acids and some consequence/HGVSp fields disagreed. | Preserve the reference initiator adjustment separately from alternate-allele translation. | `fastvep-consequence/src/predictor.rs`, coding peptide construction near line 1387 | [4 records](focused/mitochondrial-initiator/input.vcf), [field examples](focused/mitochondrial-initiator/before-fix.json) |
| Recreated stop during peptide-prefix clipping | A preceding prefix was trimmed before encountering the shared stop, producing different HGVSp. | Retain the original peptide window when prefix scanning reaches the recreated stop. | `fastvep-hgvs/src/protein.rs::clip_residues`, near line 139 | [1 record](focused/recreated-stop-clipping/input.vcf), [field examples](focused/recreated-stop-clipping/before-fix.json) |
| Inversion classification before allele clipping | Clipping a retained replacement could promote the remaining sequence to an inversion, changing HGVSc. | Classify the complete replacement before clipping; preserve that classification. | `fastvep-annotate/src/lib.rs::hgvsc_retained_replacement`, near line 1808 | [4 records](focused/inversion-before-clipping/input.vcf), [field examples](focused/inversion-before-clipping/before-fix.json) |
| Insertion flanks in stretched exon overlap | Reversed insertion endpoints caused applicable HGVSp to be absent. | Sort insertion flanks before the exon-overlap query. | `fastvep-annotate/src/lib.rs::vep_frameshift_intron_stretched_coding_overlap`, near line 1659 | [4 records](focused/insertion-flank-order/input.vcf), [field examples](focused/insertion-flank-order/before-fix.json) |

Implementation paths are under `crates/`. The exact input, expected transcript
rows and comparison contract are adjacent in each fixture directory. The small
before-fix JSON files contain selected field differences, not complete fastVEP
output baselines. Their counts cover all changed transcript identities; examples
are limited to one affected identity per input record.

## VEP source checks

The source index records file hashes from the unmodified VEP image extraction:

- Mitochondrial initiator: `Bio/EnsEMBL/Transcript.pm`, translation and reference
  start-codon adjustment; alternate peptide construction is a separate path.
- Recreated stop: `TranscriptVariationAllele.pm::_clip_alleles`, reached from
  protein HGVS formatting.
- Inversion classification: `Variation/Utils/Sequence.pm::hgvs_variant_notation`,
  followed by `TranscriptVariationAllele.pm::_clip_alleles`.
- Insertion overlap: `BaseVariationFeatureOverlapAllele.pm::_bvfo_preds`, which
  orders endpoints before the transcript overlap queries.

See [source and validation fingerprints](focused/review-index.json). These are
starting points for reviewing the relevant call chains, not a claim that a
single function accounts for every transcript configuration.

## Adopting changes

The listed implementation changes are carried in commit `546adea8`. That commit
also contains other concordance and cache-builder work; it is not four isolated
patches. The changes have not been shown to work as independent cherry-picks;
that would require extraction and testing.

Review annotation fixes separately from cache-format, supplementary-source and
application-integration changes. The focused fixtures provide a starting point
for extracting individual fixes. Before adopting one, run it against upstream's
own tests and appropriate transcript data as well.

The [full collection](FULL-SUITE.md) is supplementary evidence. It includes
documented VEP exceptions and controls that are deliberately absent from this
strict focused suite. Broader source-generation recipes and qualification
limits remain documented there.
