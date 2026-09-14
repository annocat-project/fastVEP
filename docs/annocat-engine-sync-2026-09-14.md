# AnnoCat engine synchronization, 2026-09-14

The maintained source is on `codex/vep115-concordance`, in the checkout
`D:/AnnoCat/tools/fastVEP-concordance`. Commit
`38349c52092ccfd62b6454e1844b22090bce959a` preserves the previously uncommitted
VEP 115 concordance and public transcript-cache work above `43f6bce`.

The subsequent synchronization imports AnnoCat's later direct-output interface,
CSQ cell projection, supplementary-cache decoding and packed-record performance
changes, cache profiling, explicit Cargo workspace references and dependency
lockfile. All 161 files from AnnoCat's embedded engine in checkpoint `c7a4d61`
match byte-for-byte. AnnoCat-only `ANNOCAT.md` and `qualified-snapshot.json` are
excluded from that runtime-source comparison. Existing standalone documentation,
scripts and packaging support are retained. No branch-only crate or web files
were found outside the embedded inventory.

The adjacent JSON receipt records every compared path and SHA-256. Source
synchronization does not change AnnoCat's existing executable, companion pin,
cache formats, biological behavior or previously recorded qualification evidence.
The full frozen VEP corpus and application benchmarks are not rerun merely for
copying identical files; workspace library and integration tests verify the
standalone checkout. Fresh qualification remains required for subsequent code
changes and release artifacts.

The older `tools/fastVEP-src` checkpoint `623a8b9` is preserved history, not the
current qualified engine. Use this concordance branch for the planned selective
upstream ports. No such ports are included in this synchronization.
