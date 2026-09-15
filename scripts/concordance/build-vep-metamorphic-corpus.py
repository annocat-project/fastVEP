#!/usr/bin/env python3
"""Build deterministic representation-pair inputs for VEP 115.2 qualification."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


FIELDS = (
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "REF_ALLELE",
    "UPLOADED_ALLELE",
    "DISTANCE",
    "STRAND",
    "FLAGS",
    "CANONICAL",
    "MANE",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "TSL",
    "CCDS",
    "ENSP",
    "SOURCE",
    "HGVS_OFFSET",
)
PADDED_INVARIANTS = (
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "STRAND",
    "FLAGS",
    "CANONICAL",
    "MANE",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "TSL",
    "CCDS",
    "ENSP",
    "SOURCE",
)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class IndexedFasta:
    def __init__(self, path):
        self.handle = path.open("rb")
        self.index = {}
        with Path(str(path) + ".fai").open(encoding="utf-8") as handle:
            for line in handle:
                name, length, offset, line_bases, line_bytes = line.rstrip().split("\t")[:5]
                self.index[name] = tuple(map(int, (length, offset, line_bases, line_bytes)))

    def base(self, chromosome, position):
        name = chromosome.removeprefix("chr")
        name = "MT" if name == "M" else name
        length, offset, line_bases, line_bytes = self.index.get(chromosome, self.index.get(name))
        if position < 1 or position > length:
            raise ValueError(f"{chromosome}:{position} is outside the reference")
        zero_based = position - 1
        self.handle.seek(offset + zero_based // line_bases * line_bytes + zero_based % line_bases)
        return self.handle.read(1).decode().upper()

    def close(self):
        self.handle.close()


def source_records(path, salt):
    records = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip().split("\t")
            if len(columns[3]) != 1 or len(columns[4]) != 1 or "," in columns[4]:
                continue
            score = hashlib.sha256((salt + "\t" + "\t".join(columns[:5])).encode()).digest()
            records.append((score, columns))
    return [columns for _, columns in sorted(records)]


def identity(columns):
    return {
        "chromosome": columns[0],
        "position": int(columns[1]),
        "id": columns[2],
        "reference": columns[3],
        "alternate": columns[4],
    }


def write_vcf(path, records, side):
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n##reference=GRCh38\n")
        handle.write("##source=AnnoCAT-VEP-115.2-metamorphic-v1\n")
        handle.write('##INFO=<ID=ANNOCAT_PAIR,Number=1,Type=String,Description="Metamorphic pair identifier">\n')
        handle.write('##INFO=<ID=ANNOCAT_RELATION,Number=1,Type=String,Description="Metamorphic relation">\n')
        handle.write('##INFO=<ID=ANNOCAT_SIDE,Number=1,Type=String,Description="Metamorphic input side">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for columns, pair_id, relation in records:
            handle.write(
                "\t".join(columns[:7])
                + f"\tANNOCAT_PAIR={pair_id};ANNOCAT_RELATION={relation};ANNOCAT_SIDE={side}\n"
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-vcf", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output-a", type=Path, required=True)
    parser.add_argument("--output-b", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--pairs-per-relation", type=int, default=100)
    parser.add_argument("--salt", default="annocat-vep115-metamorphic-v1")
    args = parser.parse_args()

    needed = args.pairs_per_relation * 3
    source = source_records(args.source_vcf, args.salt)
    if len(source) < needed:
        raise ValueError(f"need {needed} eligible source records; found {len(source)}")
    source = source[:needed]
    fasta = IndexedFasta(args.fasta)
    records_a = []
    records_b = []
    pairs = []

    for index, original in enumerate(source[: args.pairs_per_relation], 1):
        pair_id = f"PAD{index:03d}"
        minimal = original.copy()
        minimal[2] = f"{pair_id}_A"
        previous = fasta.base(original[0], int(original[1]) - 1)
        padded = original.copy()
        padded[1] = str(int(original[1]) - 1)
        padded[2] = f"{pair_id}_B"
        padded[3] = previous + original[3]
        padded[4] = previous + original[4]
        records_a.append((minimal, pair_id, "minimal-padded"))
        records_b.append((padded, pair_id, "minimal-padded"))
        pairs.append(
            {
                "id": pair_id,
                "relation": "minimal-padded",
                "a": [identity(minimal)],
                "b": [identity(padded)],
                "invariantFields": list(PADDED_INVARIANTS),
                "mayDiffer": [
                    "vcf_coordinates",
                    *[field for field in FIELDS if field not in PADDED_INVARIANTS],
                ],
            }
        )

    start = args.pairs_per_relation
    end = start + args.pairs_per_relation
    bases = "ACGT"
    for index, original in enumerate(source[start:end], 1):
        pair_id = f"MULTI{index:03d}"
        second_alt = next(base for base in bases if base not in {original[3], original[4]})
        unsplit = original.copy()
        unsplit[2] = f"{pair_id}_A"
        unsplit[4] = f"{original[4]},{second_alt}"
        first = original.copy()
        first[2] = f"{pair_id}_B1"
        second = original.copy()
        second[2] = f"{pair_id}_B2"
        second[4] = second_alt
        records_a.append((unsplit, pair_id, "split-unsplit"))
        records_b.extend(((first, pair_id, "split-unsplit"), (second, pair_id, "split-unsplit")))
        pairs.append(
            {
                "id": pair_id,
                "relation": "split-unsplit",
                "a": [identity(unsplit)],
                "b": [identity(first), identity(second)],
                "invariantFields": [field for field in FIELDS if field != "UPLOADED_ALLELE"],
                "mayDiffer": ["record_id", "alt_cardinality", "UPLOADED_ALLELE"],
            }
        )

    order_a = []
    order_b = []
    for index, original in enumerate(source[end:needed], 1):
        pair_id = f"ORDER{index:03d}"
        record = original.copy()
        record[2] = pair_id
        order_a.append((record, pair_id, "input-order"))
        order_b.append((record.copy(), pair_id, "input-order"))
        pairs.append(
            {
                "id": pair_id,
                "relation": "input-order",
                "a": [identity(record)],
                "b": [identity(record)],
                "invariantFields": list(FIELDS),
                "mayDiffer": ["output_record_order"],
            }
        )
    records_a.extend(order_a)
    records_b.extend(reversed(order_b))
    fasta.close()

    write_vcf(args.output_a, records_a, "A")
    write_vcf(args.output_b, records_b, "B")
    manifest = {
        "schemaVersion": 1,
        "corpusId": "ensembl-115-metamorphic",
        "role": "qualification-metamorphic",
        "assembly": "GRCh38",
        "selection": {
            "algorithm": "lowest deterministic SHA-256 source-record scores",
            "salt": args.salt,
            "pairsPerRelation": args.pairs_per_relation,
            "pairCount": len(pairs),
            "candidateOutputUsed": False,
        },
        "sources": {
            "vcf": {"path": str(args.source_vcf).replace("\\", "/"), "sha256": sha256(args.source_vcf)},
            "reference": {"path": str(args.fasta).replace("\\", "/")},
        },
        "generator": {
            "path": "scripts/build-vep-metamorphic-corpus.py",
            "sha256": sha256(Path(__file__)),
        },
        "outputs": {
            "a": {
                "path": str(args.output_a).replace("\\", "/"),
                "records": len(records_a),
                "sha256": sha256(args.output_a),
            },
            "b": {
                "path": str(args.output_b).replace("\\", "/"),
                "records": len(records_b),
                "sha256": sha256(args.output_b),
            },
        },
        "pairs": pairs,
    }
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"pairs": len(pairs), "recordsA": len(records_a), "recordsB": len(records_b)}, indent=2))


if __name__ == "__main__":
    main()
