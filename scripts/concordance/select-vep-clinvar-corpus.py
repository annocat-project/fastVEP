#!/usr/bin/env python3
"""Stream a pinned ClinVar VCF and select a deterministic reviewed edge cohort."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import heapq
import io
import json
import re
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path


EFFECTS = [
    ("splice_acceptor", {"splice_acceptor_variant"}),
    ("splice_donor", {"splice_donor_variant"}),
    ("frameshift", {"frameshift_variant"}),
    ("stop_gained", {"stop_gained", "nonsense"}),
    ("stop_lost", {"stop_lost"}),
    ("start_lost", {"start_lost", "initiator_codon_variant"}),
    ("inframe_deletion", {"inframe_deletion"}),
    ("inframe_insertion", {"inframe_insertion"}),
    ("splice_region", {"splice_region_variant"}),
    ("missense", {"missense_variant"}),
    ("synonymous", {"synonymous_variant"}),
]


class HashingReader(io.RawIOBase):
    def __init__(self, source):
        self.source = source
        self.sha256 = hashlib.sha256()
        self.md5 = hashlib.md5()
        self.bytes_read = 0

    def readable(self):
        return True

    def read(self, size=-1):
        data = self.source.read(size)
        self.sha256.update(data)
        self.md5.update(data)
        self.bytes_read += len(data)
        return data

    def readinto(self, buffer):
        data = self.read(len(buffer))
        buffer[: len(data)] = data
        return len(data)


def info_fields(raw: str) -> dict[str, str]:
    result = {}
    for field in raw.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            result[key] = value
    return result


def reviewed(status: str) -> bool:
    return (
        "practice_guideline" in status
        or "reviewed_by_expert_panel" in status
        or ("multiple_submitters" in status and "no_conflicts" in status)
    )


def significance(value: str) -> str | None:
    lowered = value.lower()
    if "pathogenic" in lowered and "conflicting" not in lowered:
        return "pathogenic"
    if "benign" in lowered and "conflicting" not in lowered:
        return "benign"
    return None


def effect(value: str) -> str | None:
    terms = {item.split("|")[-1] for item in value.split(",")}
    for label, accepted_terms in EFFECTS:
        if terms & accepted_terms:
            return label
    return None


def stable_score(columns: list[str]) -> int:
    key = "\t".join(columns[:5])
    return int.from_bytes(hashlib.sha256(key.encode()).digest()[:8], "big")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--quota", type=int, default=20)
    parser.add_argument("--expected-bytes", type=int)
    parser.add_argument("--expected-md5")
    args = parser.parse_args()

    heaps: dict[tuple[str, str], list] = defaultdict(list)
    reviewed_terms = Counter()
    scanned = eligible = 0
    request = urllib.request.Request(args.url, headers={"User-Agent": "AnnoCAT-oracle-validation/1"})
    with urllib.request.urlopen(request, timeout=120) as response:
        hashing = HashingReader(response)
        with gzip.GzipFile(fileobj=io.BufferedReader(hashing)) as compressed:
            for raw in io.TextIOWrapper(compressed, encoding="utf-8"):
                if raw.startswith("#"):
                    continue
                scanned += 1
                columns = raw.rstrip("\n").split("\t")
                if len(columns) < 8 or "," in columns[4]:
                    continue
                chrom, pos, identifier, reference, alternate = columns[:5]
                if chrom not in {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}:
                    continue
                if not re.fullmatch(r"[ACGT]+", reference) or not re.fullmatch(r"[ACGT]+", alternate):
                    continue
                if max(len(reference), len(alternate)) > 50:
                    continue
                info = info_fields(columns[7])
                review = info.get("CLNREVSTAT", "")
                sig = significance(info.get("CLNSIG", ""))
                molecular = info.get("MC", "")
                if reviewed(review) and sig:
                    reviewed_terms.update(item.split("|")[-1] for item in molecular.split(",") if item)
                category = effect(molecular)
                if not reviewed(review) or not sig or not category:
                    continue
                eligible += 1
                key = (category, sig)
                record = {
                    "chromosome": "chrM" if chrom == "MT" else f"chr{chrom}",
                    "position": int(pos),
                    "id": identifier,
                    "reference": reference,
                    "alternate": alternate,
                    "effect": category,
                    "significance": sig,
                    "reviewStatus": review,
                    "clinicalSignificance": info.get("CLNSIG", ""),
                    "molecularConsequence": info.get("MC", ""),
                }
                item = (-stable_score(columns), "\t".join(columns[:5]), record)
                heap = heaps[key]
                if len(heap) < args.quota:
                    heapq.heappush(heap, item)
                elif item[0] > heap[0][0]:
                    heapq.heapreplace(heap, item)

    source = {
        "url": args.url,
        "compressedBytes": hashing.bytes_read,
        "compressedSha256": hashing.sha256.hexdigest(),
        "compressedMd5": hashing.md5.hexdigest(),
    }
    if args.expected_bytes and hashing.bytes_read != args.expected_bytes:
        raise SystemExit(f"compressed size mismatch: {hashing.bytes_read} != {args.expected_bytes}")
    if args.expected_md5 and hashing.md5.hexdigest().lower() != args.expected_md5.lower():
        raise SystemExit(f"compressed MD5 mismatch: {hashing.md5.hexdigest()} != {args.expected_md5}")

    records = [item[2] for heap in heaps.values() for item in heap]
    records.sort(key=lambda item: (int(item["chromosome"][3:]) if item["chromosome"][3:].isdigit() else 99, item["chromosome"], item["position"], item["reference"], item["alternate"]))
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##reference=GRCh38\n")
        handle.write("##source=ClinVar-20260728-reviewed-deterministic-sample\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for record in records:
            handle.write(
                f'{record["chromosome"]}\t{record["position"]}\t{record["id"]}\t'
                f'{record["reference"]}\t{record["alternate"]}\t.\tPASS\t.\n'
            )

    manifest = {
        "schemaVersion": 1,
        "source": source,
        "selection": {
            "scannedRecords": scanned,
            "eligibleRecords": eligible,
            "quotaPerEffectAndSignificance": args.quota,
            "reviewRequirement": "ClinVar two-star or higher",
            "significanceGroups": ["pathogenic", "benign"],
            "algorithm": "lowest deterministic SHA-256 scores within each effect/significance stratum",
            "observedReviewedMolecularConsequences": dict(reviewed_terms.most_common()),
        },
        "counts": {
            f"{category}/{sig}": len(heap)
            for (category, sig), heap in sorted(heaps.items())
        },
        "variantCount": len(records),
        "variants": records,
    }
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"scanned": scanned, "eligible": eligible, "selected": len(records), "counts": manifest["counts"]}, indent=2))


if __name__ == "__main__":
    main()
