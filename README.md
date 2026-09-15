# fastVEP for AnnoCat

This fork of [Huang-lab/fastVEP](https://github.com/Huang-lab/fastVEP) is maintained
for [AnnoCat](https://github.com/annocat-project/AnnoCAT). It extends the annotation
engine, transcript-cache builder, supplementary-source handling and library
interfaces used by the application. It also includes performance and memory-use
changes for large local annotation jobs. Ensembl VEP 115.2 concordance under
the Ensembl human GRCh38 configuration described below is one part of that work.

The standalone fastVEP CLI and AnnoCat's embedded annotation engine use the same
Rust source. AnnoCat owns the application, Parquet writer and viewer.

## Differences from upstream

The comparison is against recorded upstream base `0e13c5bdb92f22a5f780cf314699f8402fb8383e`
(v0.3.0). It does not cover every later upstream change.

- Additional consequence and HGVS repairs for allele minimization, transcript
  membership, splice boundaries, incomplete CDS, start/stop codons and repeats.
- A public Ensembl Core transcript-cache builder with source metadata and
  declared sequence edits. FSTVEP02 reading is retained alongside ANNOCATC1.
- Supplementary-source builders and OSA1/OSA2 readers, verified chromosome
  shards, multiple records per allele key and transcript-aligned dbNSFP fields.
- Shared CLI/library annotation paths and output interfaces for AnnoCat's
  direct-Parquet worker. The Parquet writer and viewer schemas live in AnnoCat.
- Memory-mapped reference and supplementary data, bounded decoded caches,
  packed supplementary records, reduced parsing and cloning, and deterministic
  parallel source loading and querying.

The [concordance inventory](tests/concordance/FIX-INVENTORY.md) describes 113
review families with source and related test links. They include supporting
cache and caller changes; they are not 113 independent bugs or isolated patches.

## AnnoCat integration

AnnoCat manages source downloads and installation. fastVEP builds and verifies
the transcript and supplementary caches, then supplies annotations through its
library or standalone commands. The transcript builder combines matching GFF3,
reference FASTA and Ensembl Core inputs; older FSTVEP02 caches remain readable,
but do not gain metadata that their files lack.

Supplementary-source changes preserve multiple exact-key records and
transcript-specific evidence. Source lookup keys remain distinct from HGVS
strings: providers use uploaded allele mappings, positions or normalized
reference alleles according to their own contracts. Supplementary lookup and
output parity need separate checks from VEP consequence agreement.

The library exposes annotations and CSQ field projection in memory for AnnoCat's
direct-Parquet worker. AnnoCat handles Parquet writing, result schemas, displayed
transcript selection, task recovery and the viewer. The standalone CLI retains
VCF and structured text output; it does not provide an AnnoCat direct-Parquet
command.

The performance changes address data loading, decoding and output preparation.
Their effect depends on the input, enabled sources and hardware. This README
does not claim a general speed or memory advantage over current upstream.

## ANNOCATC1 transcript caches

ANNOCATC1 is this fork's versioned transcript-cache format. GFF3 alone does not
provide all the transcript metadata and sequence adjustments needed to reproduce
the tested VEP 115.2 / Ensembl GRCh38 configuration. The public builder combines GFF3 and indexed reference
FASTA with six Ensembl Core tables: `gene`, `transcript`, `translation`,
`attrib_type`, `transcript_attrib` and `translation_attrib`. Input sizes and
SHA-256 hashes are checked against the build manifest.

The cache contains transcript and gene identifiers, exon structure, coding
coordinates, transcript/CDS/peptide sequences, and annotations such as MANE,
canonical status, APPRIS, TSL, CCDS and completeness flags where available.
An enrichment map adds gene versions, mature-miRNA ranges and declared sequence
edits. The current public builder supports translation edits and rejects RNA
edits it cannot apply; storing an edit field does not imply support for every
edit type.

The file has an `ANNOCATC` signature, a version number and an inspectable JSON
header, followed by a Zstandard-compressed bincode payload. The header records
species, assembly, Ensembl/VEP releases, capabilities, source and builder
provenance, transcript counts and two hashes: one for the compressed payload and
one for its serialized content before compression. Loading checks those hashes,
file length, transcript/contig counts and enrichment references. The serialized
field layout is frozen for this format version.

FSTVEP02 remains readable through the compatibility loader. It does not acquire
ANNOCATC1 provenance or enrichment merely by being opened. ANNOCATC1 is also
separate from Ensembl's own VEP cache format and from AnnoCat result files.
Successful decoding and checksum verification establish file integrity;
annotation concordance still requires comparison with the reference, cache and
options under test.

See the [format reader/writer](crates/fastvep-cache/src/annocat_cache.rs),
[frozen payload layout](crates/fastvep-cache/src/transcript_wire.rs) and
[public builder](crates/fastvep-cache/src/ensembl_core.rs).

## OSA2 supplementary caches

OSA2 is the upstream ZIP-based supplementary-annotation format. It groups
records into genomic chunks with encoded allele keys, field arrays, string
tables and optional structured JSON values. This fork extends its handling for
AnnoCat's sources and large local workloads.

- **Record preservation:** lookup retains multiple records with the same allele
  key. The `record_list` metadata flag distinguishes a list of source records
  from an array-valued annotation. Source readers also preserve transcript-aligned
  fields, such as dbNSFP evidence, through result projection.
- **Alleles outside the compact encoding:** an optional `raw-alleles.enc` entry
  retains exact allele strings that cannot use the compact A/C/G/T encoding,
  rather than dropping those source records.
- **Chromosome shards:** a separate `.osa-shards.json` manifest presents multiple
  cache files as one source. The reader checks source metadata and format
  consistency, chromosome mappings and relative file paths. Sharding does not
  replace the OSA2 container format.
- **Loading and memory:** the reader parses the ZIP directory from a memory map
  and resolves entry offsets on demand. Readers share a byte-budgeted decoded
  chunk cache. Packed JSON storage reduces per-record allocations, and decoding
  limits bound individual structured-value columns.

The fork continues to read OSA1 and OSA2. Runtime improvements such as shared
decoded caches do not require rewriting installed files. The `record_list` and
raw-allele extensions do require reader support; compatibility with an older
upstream reader must not be assumed simply because the file is named `.osa2`.
OSA2 JSON values still use inner Zstandard compression within the ZIP container;
the runtime work does not remove that on-disk compression layer.

See the [OSA2 reader](crates/fastvep-sa/src/reader_v2.rs),
[writer](crates/fastvep-sa/src/writer_v2.rs) and
[shard reader](crates/fastvep-sa/src/sharded.rs). These storage and lookup changes
have their own correctness checks; the VEP concordance totals below do not
qualify every supplementary source or score.

## VEP concordance testing and limits

### Concordance target

The target was to reproduce **Ensembl VEP 115.2's per-allele, per-transcript
annotations for human GRCh38**, using the official Ensembl release-115 indexed
cache and a reference FASTA pinned by checksum. VEP ran offline with HGVS
enabled and a 5,000-base upstream/downstream distance. MANE, canonical status,
TSL, CCDS, gene symbols, biotypes, protein identifiers and exon/intron numbers
were requested as output metadata. All applicable transcript annotations were
retained; MANE and canonical labels did not restrict which transcripts were
compared.

Concordance covers transcript identities and their multiplicity as well as the
declared consequence, impact, HGVSc/HGVSp, cDNA/CDS/protein position, amino-acid,
codon, allele and transcript-metadata fields. Missing or extra transcript
annotations count as differences. The target is agreement under each test's
comparison contract, including its explicit normalizers and narrowly documented
exceptions for VEP behavior we intentionally do not reproduce. It is not a
claim of byte-identical VCF files or agreement with every VEP option, plugin or
database configuration.

The [recorded VEP configuration](tests/concordance/data-manifest.json) supplies
the exact image and cache fingerprints, command-line options and output field
list used by the runnable comparison workflow. Supplementary-source scores and
AnnoCat's displayed transcript selection are separate from this VEP target.

### Evidence and limits

Start with the [four regression groups](tests/concordance/UPSTREAM-REVIEW.md):
13 input records, frozen VEP outputs and one test command. Each group reproduces
a failure in a retained pre-fix fork build and passes with the tested fork.
These tests do not establish that current upstream still has those defects.

The collection contains 216 inputs and 2,531,230 record observations:
2,451,891 distinct exact record keys and 2,180,583 distinct exact split-ALT keys.
Inputs include public ClinVar/GIAB selections and generated cases on autosomes,
X, Y and mitochondrial sequence. They cover coding and noncoding transcripts,
both strands, splice boundaries, UTRs, partial CDS, repeats and multiallelic
representations.

These totals mix oracle comparisons, historical output-preservation checks and
completion controls for inputs on which VEP fails. They are not 2.5 million
independent exact matches or a population-wide accuracy estimate. Contracts
specify fields, normalizers and narrowly defined exceptions. Normalized or
exception-qualified equality is not raw byte equality.

The tested configuration does not qualify every inherited option, structural
variants, RefSeq configuration or non-human dataset. Record the fork commit,
cache, reference and options when reproducing a result.

See the [testing guide](tests/concordance/FULL-SUITE.md) for methods, counts,
accepted differences and commands. Inputs and generators are provided; the
larger frozen VEP outputs are optional downloads. Full fastVEP output baselines
and internal audit ledgers are not included.

## Build

With a current Rust toolchain, run from the repository root:

```sh
cargo build --release --locked -p fastvep-cli
target/release/fastvep --version
```

On Windows, the executable is `target/release/fastvep.exe`.

## Run the regression examples

Provide the matching indexed reference FASTA and transcript cache:

```sh
python scripts/concordance/run.py --data tests/concordance/focused --binary ./target/release/fastvep --cache ./transcripts.cache --fasta ./reference.fa --output ./focused-results
```

Python 3.11 or later is required. The focused fixture includes VEP outputs;
the FASTA and transcript cache are not bundled.

## Citation and license

Cite the [upstream fastVEP paper](https://doi.org/10.64898/2026.04.14.718452)
and record the fork commit and source database releases used for annotation.
This fork retains the [Apache License 2.0](LICENSE.md).
Report fork-specific problems through
[the fork's issue tracker](https://github.com/annocat-project/fastVEP/issues).
