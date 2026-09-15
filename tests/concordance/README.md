# VEP 115 concordance tests

Start with the [upstream review guide](UPSTREAM-REVIEW.md): four regression
groups, 13 input records and one command. The focused fixtures include frozen
VEP outputs, so the larger oracle download is unnecessary for that check.

The [concordance inventory](FIX-INVENTORY.md) explains the broader changes.
The [full-suite guide](FULL-SUITE.md) describes all 216 frozen inputs, comparison
contracts, evidence limits and optional VEP outputs.

Build the fork and supply a matching transcript cache and indexed GRCh38 FASTA.
Python 3.11 or later is required. From the repository root:

```sh
python scripts/concordance/run.py --data tests/concordance/focused --verify-only
python scripts/concordance/run.py --data tests/concordance/focused --binary ./fastvep --cache ./transcripts.cache --fasta ./reference.fa --output ./focused-results
```

The runner checks input hashes, annotates all transcripts and rejects failed
comparisons. The output directory must be new. The contracts define fields and
normalizers; testing a different cache does not inherit previous qualification.

Inputs are generated/public scientific fixtures. The full input collection
retains one synthetic genotype control; these are not user samples. Input and
oracle hashes and annotation settings are recorded in the manifests.

The comparator also has negative controls:

```sh
python scripts/concordance/compare-vep-concordance.py --self-test
```

Local result directories are not part of the publication package. Full fastVEP
output baselines and internal development ledgers remain excluded.
