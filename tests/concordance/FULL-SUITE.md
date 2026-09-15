# Full concordance evidence archive — review draft

This archive is intended for the fastVEP fork. It is not published yet.
Start with the [four regression groups](UPSTREAM-REVIEW.md) for a small runnable
example. This archive contains the broader evidence behind the testing counts.

The full inventory contains 216 frozen input files and 2,531,230 record
observations. Across those files there are 2,451,891 distinct exact record keys
and 2,180,583 distinct exact split-ALT keys. These counts preserve chromosome
spelling, position, reference and alternate sequence. They do not collapse
equivalent padded or shifted representations into biological events.

The retained oracle files match the complete input record sets for 206 banks.
The other ten are identified in the manifest as control-only or unavailable.
The completion inventory cites historical comparisons for 140 inputs, output
preservation for 68, native preservation with an oracle failure for five, and
direct comparisons for three. Those historical evidence types are retained in
the manifest; the packaging audit does not replace them with new test results.

The inputs include public ClinVar and GIAB selections, generated transcript
boundary and feature cases, representation pairs, source-path probes, and
single-allele controls. They cover autosomes, X, Y and mitochondrial sequence.
There are no symbolic or breakend alleles in this inventory. Reused loci and
targeted discovery cases are not an independent probability sample.

## What the files prove

The manifest distinguishes these forms of evidence:

- **VEP output with matching record keys:** a retained oracle output whose
  chromosome, position, reference, alternate and record multiplicities match
  the input. Matching keys establish the pairing, not annotation concordance.
- **Controls only or unavailable:** a bank whose complete input has no matched
  oracle in the retained set. Single-allele controls or surviving subsets can
  provide evidence, but cannot be presented as a complete raw VEP comparison.

Comparison contracts retain their declared fields, normalizers and exception
entries. Their presence is historical evidence; do not assume every associated
contract can be applied unchanged to every oracle or mixed-allele control.
For the four regression groups, use the tested `run.py` command in
[the review guide](UPSTREAM-REVIEW.md).

## Verify and compare

The input download contains all frozen inputs, contracts and provenance.
VEP oracle outputs are a separate optional download. Extract downloaded parts
into the same parent directory to populate `vep115-tests`. No fastVEP output
baselines or baseline-replay tool are included in the publication package.
Those remain local development evidence.

Verify the input package, and optionally the downloaded VEP outputs:

```sh
python scripts/concordance/verify-assets.py --data vep115-tests
python scripts/concordance/verify-assets.py --data vep115-tests --oracles
```

Generate candidate annotations with the fastVEP build under test. Use the
manifest to locate the input, VEP output and contract for the selected bank.
Decompress the selected input to `.vcf` before comparison: reviewed exception
entries bind to the uncompressed input checksum. The general comparator accepts
these files:

```sh
python scripts/concordance/compare-vep-concordance.py candidate.vcf oracle.vcf.gz --input input.vcf --contract comparison-contract.json --json comparison.json
```

Controls and documented exceptions require their corresponding comparison
method; a matching record set alone does not establish concordance. The
four regression groups have their own tested `run.py` workflow. Developers
who regenerate VEP outputs can skip the optional oracle download and use the
pinned profile with `run-oracle.py`.

## How the concordance work was performed

Inputs were frozen before candidate outputs were inspected. Comparisons used
the local Ensembl VEP 115.2 oracle with the official human release-115 GRCh38
cache. The work compared all transcript annotations rather than selecting one
representative transcript. Differences were traced through the VEP source,
reduced to targeted cases and checked again after fixes. Historical banks were
replayed to detect regressions, and later banks added fresh discovery cases.

Source-path checks supplemented output comparisons. Unreachable paths in the
selected profile, measurement limitations and explicit exceptions were tracked
separately. Closing that review does not establish exhaustive combinations of
all inputs, full branch coverage, or correctness for other VEP profiles.

Accepted differences include documented mixed-allele state behavior in VEP,
reversed coding-to-UTR HGVS endpoint order, terminal offsets beyond the final
exon, and sequence-checked stop-distance cases. Failed VEP runs are not treated
as expected successful annotations. Preserve the raw discrepancy when reviewing
one of these cases; do not turn it into a general field exemption.

## Generating additional inputs

The existing generation tools are included with their command-line interfaces:

- `select-vep-clinvar-corpus.py`: deterministic reviewed ClinVar selections.
- `select-vep-giab-discovery.py`: GIAB benchmark and difficult-region selections.
- `build-vep-boundary-corpus.py`: transcript boundary and risk-feature probes.
- `build-vep-feature-corpus.py`: additional feature cases; uses the boundary helper.
- `build-vep-metamorphic-corpus.py`: paired sequence representations.
- `check-vep-insertion-allele-independence.py`: single-allele insertion controls.

Run each with `--help` for required inputs. A new source release, salt or
exclusion set produces a different bank. These tools do not, by themselves,
reconstruct all 216 historical files: the frozen inputs are the reproducible
replay inputs. Reconstructing targeted historical probes also requires their
original generation recipes.

The subset manifest records the pinned VEP image, cache URL and checksum,
field list and oracle arguments. Supply explicit local cache, FASTA, input and
output paths when running that oracle. Keep each oracle exit status and input
hash; never substitute native output for missing VEP output.

`run-oracle.py` supplies these arguments to a local VEP installation. Use its
`--help` output for the required paths. Run it inside the pinned runtime with
the recorded cache. The wrapper records file hashes and exit status, but does
not verify the complete Perl installation or the unpacked cache on your behalf.

## Publication transformations

Execution metadata headers are omitted. Functional VCF headers and variant
values are retained, with LF line endings. One deliberately generated `CONTROL`
genotype column is retained; participant sample columns are not accepted by the
packaging script. No application results, logs, transcript caches or reference
genomes are included.

The explicitly synthetic prior-VEP header used by the header-handling control is
also retained.

The manifest records original file hashes, published file hashes and record-body
hashes. Exception entries are inlined from their original ledgers; an input hash
is updated only when it matches an inventoried input whose headers were
sanitized. Original contract and ledger hashes remain in the manifest.

This is an evidence review archive. Publication review and annotation
qualification are separate checks, and both must be described accurately.
