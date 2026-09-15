#!/usr/bin/env python3
"""Select a deterministic, candidate-blind GIAB VEP discovery corpus."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import heapq
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


CONTIG_ORDER = {f"chr{value}": value for value in range(1, 23)}
CONTIG_ORDER.update({"chrX": 23, "chrY": 24})
QUOTAS = {
    "cmrg": {
        "snv": 400,
        "insertion_short": 225,
        "deletion_short": 225,
        "insertion_long": 25,
        "deletion_long": 25,
        "complex": 100,
    },
    "difficult": {
        "snv": 400,
        "insertion_short": 225,
        "deletion_short": 225,
        "insertion_long": 25,
        "deletion_long": 25,
        "complex": 100,
    },
    "not_difficult": {
        "snv": 460,
        "insertion_short": 245,
        "deletion_short": 250,
        "insertion_long": 5,
        "deletion_long": 0,
        "complex": 40,
    },
}
COHORTS = ("cmrg", "difficult", "not_difficult")
SOURCE_PINS = {
    "cmrg": {
        "url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz",
        "compressedBytes": 231237,
        "compressedSha256": "80f18efa283d953cfc615d30dd01b1ff80debcc46f52042b942da783acd48da2",
    },
    "benchmark": {
        "url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        "compressedBytes": 156252944,
        "compressedSha256": "adb4d4a50048aa13353a06b84fcfcbca09a5d17525efaa4cea44f8822e81175c",
    },
    "difficultRegions": {
        "url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/Union/GRCh38_alldifficultregions.bed.gz",
        "compressedBytes": 32855595,
        "compressedSha256": "d44a65699cb8056578c629039b0f78a8b34a7e36eb720eadbf45c53e18fbb21f",
        "compressedMd5": "38b6844660d937e3d13d0b32ad3b2358",
    },
}


def file_digest(path, algorithm="sha256"):
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_pin(path, pin):
    if path.stat().st_size != pin["compressedBytes"]:
        raise SystemExit(f"compressed size mismatch for {path}")
    actual = file_digest(path)
    if actual != pin["compressedSha256"]:
        raise SystemExit(f"compressed SHA-256 mismatch for {path}: {actual}")
    if "compressedMd5" in pin:
        actual_md5 = file_digest(path, "md5")
        if actual_md5 != pin["compressedMd5"]:
            raise SystemExit(f"compressed MD5 mismatch for {path}: {actual_md5}")


def allele_bucket(reference, alternate):
    if len(reference) == len(alternate):
        return "snv" if len(reference) == 1 else "complex"
    if reference[0] == alternate[0]:
        if len(reference) == 1:
            return "insertion_long" if len(alternate) > 50 else "insertion_short"
        if len(alternate) == 1:
            return "deletion_long" if len(reference) > 50 else "deletion_short"
    return "complex"


def called_alternates(columns):
    alternates = columns[4].split(",")
    if len(columns) < 10:
        return enumerate(alternates, 1)
    formats = columns[8].split(":")
    values = columns[9].split(":")
    genotype = values[formats.index("GT")] if "GT" in formats and formats.index("GT") < len(values) else ""
    called = {int(value) for value in re.split(r"[|/]", genotype) if value.isdigit() and int(value) > 0}
    return ((index, alt) for index, alt in enumerate(alternates, 1) if not called or index in called)


def selection_score(record, salt):
    identity = f'{record["chromosome"]}\t{record["position"]}\t{record["reference"]}\t{record["alternate"]}'
    return int.from_bytes(hashlib.sha256(f"{salt}\t{identity}".encode()).digest()[:8], "big")


def records(path, cohort, salt):
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 8 or columns[0] not in CONTIG_ORDER or columns[6] not in {"PASS", "."}:
                continue
            reference = columns[3].upper()
            if not re.fullmatch(r"[ACGT]+", reference):
                continue
            for alt_index, alternate in called_alternates(columns):
                alternate = alternate.upper()
                if not re.fullmatch(r"[ACGT]+", alternate) or max(len(reference), len(alternate)) > 200:
                    continue
                record = {
                    "chromosome": columns[0],
                    "position": int(columns[1]),
                    "id": columns[2],
                    "reference": reference,
                    "alternate": alternate,
                    "sourceAltIndex": alt_index,
                    "cohort": cohort,
                    "bucket": allele_bucket(reference, alternate),
                }
                record["score"] = selection_score(record, salt)
                yield record


class DifficultRegions:
    """Stream a coordinate-sorted BED alongside a coordinate-sorted VCF."""

    def __init__(self, path):
        self.handle = gzip.open(path, "rt", encoding="utf-8", errors="strict")
        self.current = None
        self._advance()

    def close(self):
        self.handle.close()

    def _advance(self):
        for line in self.handle:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end, *_ = line.rstrip("\n").split("\t")
            if chrom in CONTIG_ORDER:
                self.current = (chrom, int(start), int(end))
                return
        self.current = None

    def overlaps(self, chrom, one_based_start, reference_length):
        start = one_based_start - 1
        end = start + reference_length
        rank = CONTIG_ORDER[chrom]
        while self.current is not None:
            current_chrom, _, current_end = self.current
            current_rank = CONTIG_ORDER[current_chrom]
            if current_rank < rank or (current_rank == rank and current_end <= start):
                self._advance()
                continue
            break
        if self.current is None:
            return False
        current_chrom, current_start, current_end = self.current
        return current_chrom == chrom and current_start < end and current_end > start


def retain(heap, record, limit):
    item = (-record["score"], f'{record["chromosome"]}:{record["position"]}:{record["reference"]}:{record["alternate"]}', record)
    if len(heap) < limit:
        heapq.heappush(heap, item)
    elif item[0] > heap[0][0]:
        heapq.heapreplace(heap, item)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cmrg", type=Path, required=True)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--difficult-bed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--salt", default="annocat-vep115-giab-discovery-v1")
    args = parser.parse_args()

    for key, path in (("cmrg", args.cmrg), ("benchmark", args.benchmark), ("difficultRegions", args.difficult_bed)):
        verify_pin(path, SOURCE_PINS[key])

    heaps = defaultdict(list)
    eligible = Counter()
    limit_multiplier = 4
    for record in records(args.cmrg, "cmrg", args.salt):
        key = ("cmrg", record["bucket"])
        eligible[key] += 1
        retain(heaps[key], record, QUOTAS["cmrg"][record["bucket"]] * limit_multiplier)

    regions = DifficultRegions(args.difficult_bed)
    try:
        for record in records(args.benchmark, "benchmark", args.salt):
            cohort = "difficult" if regions.overlaps(record["chromosome"], record["position"], len(record["reference"])) else "not_difficult"
            record["cohort"] = cohort
            key = (cohort, record["bucket"])
            eligible[key] += 1
            retain(heaps[key], record, QUOTAS[cohort][record["bucket"]] * limit_multiplier)
    finally:
        regions.close()

    selected = []
    counts = Counter()
    identities = set()
    loci = set()
    for cohort in COHORTS:
        for bucket, quota in QUOTAS[cohort].items():
            if quota == 0:
                continue
            queue = sorted((item[2] for item in heaps[(cohort, bucket)]), key=lambda item: item["score"])
            for record in queue:
                identity = (record["chromosome"], record["position"], record["reference"], record["alternate"])
                locus = identity[:2]
                if identity in identities or locus in loci:
                    continue
                selected.append(record)
                identities.add(identity)
                loci.add(locus)
                counts[(cohort, bucket)] += 1
                if counts[(cohort, bucket)] == quota:
                    break
            if counts[(cohort, bucket)] != quota:
                raise RuntimeError(f"unable to fill {cohort}/{bucket}: {counts[(cohort, bucket)]} of {quota}")

    selected.sort(key=lambda item: (CONTIG_ORDER[item["chromosome"]], item["position"], item["reference"], item["alternate"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n##reference=GRCh38\n")
        handle.write("##source=GIAB-HG002-CMRG-v1.00-and-NIST-v4.2.1-candidate-blind-discovery\n")
        handle.write('##INFO=<ID=ANNOCAT_GIAB_STRATUM,Number=1,Type=String,Description="Candidate-blind GIAB selection stratum">\n')
        handle.write('##INFO=<ID=ANNOCAT_SHAPE,Number=1,Type=String,Description="Input allele-shape stratum">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for index, record in enumerate(selected, 1):
            identifier = record["id"] if record["id"] != "." else f"annocat-giab-{index}"
            handle.write(
                f'{record["chromosome"]}\t{record["position"]}\t{identifier}\t{record["reference"]}\t{record["alternate"]}'
                f'\t.\tPASS\tANNOCAT_GIAB_STRATUM={record["cohort"]};ANNOCAT_SHAPE={record["bucket"]}\n'
            )

    manifest = {
        "schemaVersion": 1,
        "corpusId": "ensembl-115-giab-discovery",
        "role": "candidate-blind-targeted-discovery",
        "assembly": "GRCh38",
        "sources": SOURCE_PINS,
        "generator": {
            "path": "scripts/select-vep-giab-discovery.py",
            "sha256": file_digest(Path(__file__)),
            "command": "python scripts/select-vep-giab-discovery.py --cmrg <pinned-cmrg-vcf.gz> --benchmark <pinned-hg002-v4.2.1-vcf.gz> --difficult-bed <pinned-giab-v3.1-bed.gz> --output fixtures/fastvep/ensembl-115-giab-discovery.vcf --manifest fixtures/fastvep/ensembl-115-giab-discovery-manifest.json",
        },
        "selection": {
            "salt": args.salt,
            "cohorts": list(COHORTS),
            "quotas": QUOTAS,
            "maximumAlleleLength": 200,
            "calledAlternateAllelesOnly": True,
            "locusCap": 1,
            "usesCandidateAnnotations": False,
            "difficultDefinition": "variant reference interval overlaps GIAB v3.1 GRCh38 Union/GRCh38_alldifficultregions.bed.gz",
        },
        "eligibleCounts": {f"{cohort}/{bucket}": eligible[(cohort, bucket)] for cohort in COHORTS for bucket in QUOTAS[cohort]},
        "counts": {f"{cohort}/{bucket}": counts[(cohort, bucket)] for cohort in COHORTS for bucket in QUOTAS[cohort]},
        "variantCount": len(selected),
        "output": {
            "path": str(args.output).replace("\\", "/"),
            "records": len(selected),
            "alleles": len(selected),
            "sha256": file_digest(args.output),
        },
    }
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"selected": len(selected), "counts": manifest["counts"], "outputSha256": manifest["output"]["sha256"]}, indent=2))


if __name__ == "__main__":
    main()
