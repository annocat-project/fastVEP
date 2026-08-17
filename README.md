# fastVEP

[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE.md)
[![Upstream](https://img.shields.io/badge/upstream-Huang--lab%2FfastVEP-4057d6.svg)](https://github.com/Huang-lab/fastVEP)
[![Used by AnnoCAT](https://img.shields.io/badge/used%20by-AnnoCAT-4057d6.svg)](https://github.com/annocat-project/AnnoCAT)

fastVEP is a Rust implementation of the Ensembl Variant Effect Predictor model.
It predicts transcript and protein consequences for VCF variants and can attach
clinical, population, conservation, splicing, and prediction annotations from
local supplementary databases.

The upstream engine handles small variants and structural variants, reports 49
[Sequence Ontology](http://www.sequenceontology.org/) consequence terms,
generates HGVS descriptions, parses multi-sample genotypes, and writes VCF, tab,
or JSON output. It supports merged Ensembl and RefSeq gene models, custom VCF or
BED annotations, expression-based filtering, and non-human genomes with matching
GFF3 and reference data.

This repository is the fastVEP fork maintained for AnnoCAT. It is based on
[Huang-lab/fastVEP](https://github.com/Huang-lab/fastVEP) and remains licensed
under Apache-2.0. Use the upstream repository for the original project and the
AnnoCAT fork for the exact engine bundled with AnnoCAT.

## AnnoCAT fork

AnnoCAT processes large local VCF files, installs annotation sources as verified
chromosome shards, and stores complete structured evidence for its result viewer.
This fork adds the interfaces and correctness rules needed for that workflow.

### Supplementary annotation caches

- Read OSA1 and OSA2 caches through one verified provider interface.
- Build chromosome-sharded databases from plain, gzip, BGZF, or streamed input.
- Read indexed BGZF chromosome ranges without dropping complete boundary records.
- Compose chromosome shards from strict manifests, with deterministic fallback to
  an all-chromosome shard when a source provides one.
- Retain selected source fields for ClinVar, dbSNP, gnomAD, dbNSFP, CADD,
  SpliceAI, REVEL, and conservation sources. Preserve dbNSFP transcript-aligned
  fields instead of collapsing them during cache construction.
- Preserve multiple exact-key records for sources such as dbNSFP, SpliceAI, and
  REVEL when their evidence is transcript- or gene-scoped.
- Preserve ambiguous-reference allele keys without changing the compact path for
  ordinary A, C, G, and T alleles.
- Verify archive structure, checksums, metadata, JSON records, key ordering, and
  lookup parity before a cache is promoted.
- Convert verified CADD and SpliceAI OSA1 shards to OSA2 without replacing the
  source cache.

### Annotation output

- Write newline-delimited structured annotations beside VCF output in the same
  annotation pass.
- Store supplementary evidence once per allele instead of once per transcript.
- Optionally omit duplicate supplementary VCF INFO fields while retaining the
  complete structured evidence used by AnnoCAT.
- Keep annotated VCF output complete when that output is requested.

### Consequences and HGVS

- Incorporate upstream fastVEP v0.3.0 correctness and mitochondrial updates.
- Normalize each alternate allele independently, including complex and
  reverse-strand alleles.
- Preserve transcript metadata needed for deterministic consequence selection.
- Correct consequence and HGVS handling for multiallelic records, noncoding
  genes, splice-boundary deletions, protein versions, and mitochondrial start
  uncertainty.

### Performance and safety

- Memory-map indexed FASTA and OSA2 archives instead of loading complete files.
- Stream dbNSFP, CADD, and REVEL records into cache writers without staging full
  source tables.
- Reuse parser buffers, serialize selected fields directly, and avoid redundant
  JSON parsing, record cloning, and sorting.
- Parse configurable supplementary sources in bounded ordered batches.
- Limit source parsing to one to four workers and overlap ordered OSA compression
  with parsing.
- Query independent sources in parallel while preserving deterministic output.
- Open sharded supplementary providers in deterministic parallel order.
- Share byte-bounded OSA caches across readers and keep decoded records in
  contiguous storage.
- Read OSA2 chunks directly from mapped archives.
- Read OSA2 ZIP central directories sequentially and defer local-header reads
  until the corresponding chunk is queried.
- Keep installed transcript caches read-only during annotation and fully verify
  newly built transcript caches before installation.
- Summarize repeated missing-contig warnings instead of emitting one warning per
  transcript.
- Report aggregate phase timings, source lookups, cache hits and misses, decoded
  bytes, ZIP inflation, JSON decoding, serialization, and blocked output time
  without recording variant values.

The executable keeps the upstream fastVEP version. AnnoCAT identifies a tested
fork build by Git commit, Cargo lockfile hash, and binary checksum. The current
pin and ordered change ledger are in
[`config/fastvep-pin.json`](https://github.com/annocat-project/AnnoCAT/blob/main/config/fastvep-pin.json).

## Build

Install a current Rust toolchain, clone this repository, and build the CLI:

```bash
cargo build --release --locked -p fastvep-cli
target/release/fastvep --version
```

On Windows, the binary is `target\release\fastvep.exe`.

## Usage

Build and verify a transcript cache:

```bash
fastvep cache \
  --gff3 Homo_sapiens.GRCh38.115.gff3 \
  --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --output grch38.fastvep.cache

fastvep cache-verify \
  --input grch38.fastvep.cache \
  --require-primary-coding-sequences
```

Annotate a VCF with a transcript cache and supplementary databases:

```bash
fastvep annotate \
  --input variants.vcf.gz \
  --output annotated.vcf \
  --transcript-cache grch38.fastvep.cache \
  --sa-dir annotation-databases \
  --hgvs --symbol --canonical
```

Write the structured annotations used by AnnoCAT without duplicating
supplementary evidence in the temporary VCF:

```bash
fastvep annotate \
  --input variants.vcf.gz \
  --output annotated.vcf \
  --transcript-cache grch38.fastvep.cache \
  --sa-dir annotation-databases \
  --structured-output annotations.ndjson \
  --omit-supplementary-vcf \
  --hgvs --symbol --canonical
```

## Supplementary databases

`sa-build --format auto` selects OSA2 for supported sources and OSA1 for the
remaining sources. Use `--format osa` or `--format osa2` only when a specific
format is required.

```bash
fastvep sa-build \
  --source clinvar \
  --input clinvar.vcf.gz \
  --output clinvar \
  --assembly GRCh38

fastvep sa-verify --input clinvar.osa2 --assembly GRCh38
```

For chromosome-sharded sources, build one cache per chromosome and place the
verified shards and manifest in the same source directory. `--sa-dir` loads
compatible OSA1, OSA2, interval, and gene-level providers from that directory.

The `sa-convert` command converts verified CADD or SpliceAI OSA1 shards to OSA2.
It never overwrites its input:

```bash
fastvep sa-convert --input chr1.osa --output chr1.osa2
fastvep sa-verify --input chr1.osa2 --chromosome 1 --assembly GRCh38
```

See [Supplementary annotations](docs/SUPPLEMENTARY_ANNOTATIONS.md) for source
schemas and output fields. See [ACMG](docs/ACMG.md) and
[ACMG setup](docs/ACMG_SETUP.md) for the upstream experimental classification
workflow.

## Commands

| Command | Purpose |
|---|---|
| `annotate` | Predict consequences and attach supplementary annotations |
| `cache` | Build a transcript cache from GFF3 and optional FASTA data |
| `cache-verify` | Fully decode and validate a transcript cache |
| `sa-build` | Build an OSA or interval annotation database |
| `sa-convert` | Convert a verified CADD or SpliceAI OSA1 shard to OSA2 |
| `sa-verify` | Fully validate an OSA1 or OSA2 database |
| `filter` | Filter VEP-compatible annotated VCF output |
| `web` | Launch the upstream interactive web interface |

Run `fastvep <command> --help` for the complete command reference.

## Development

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```

Changes to consequence prediction, HGVS, cache encoding, or source parsing must
keep the existing compatibility and parity tests passing. AnnoCAT release builds
also verify the pinned history, dependency lock, complete test suite, and binary
checksum before packaging.

## Citation

If you use fastVEP in research, cite the upstream project:

> Kuan-lin Huang. **fastVEP: A Fast, Comprehensive Variant Effect Predictor
> Written in Rust.** bioRxiv (2026).
> [doi:10.64898/2026.04.14.718452](https://doi.org/10.64898/2026.04.14.718452)

When reproducibility matters, also record the AnnoCAT fork commit and the source
database releases used for annotation.

## License

fastVEP and this fork are licensed under the [Apache License 2.0](LICENSE.md).
fastVEP is inspired by
[Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) and
[Illumina Nirvana](https://github.com/Illumina/Nirvana). The OSA2 format uses
encoding techniques derived from [echtvar](https://github.com/brentp/echtvar).

Report fork-specific problems through
[annocat-project/fastVEP issues](https://github.com/annocat-project/fastVEP/issues).
