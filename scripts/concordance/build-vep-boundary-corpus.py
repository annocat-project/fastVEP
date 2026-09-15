#!/usr/bin/env python3
"""Build a deterministic VCF from Ensembl transcript geometry, not fastVEP calls."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path



def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


PRIMARY = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
TRANSCRIPT_FEATURES = {
    "mRNA",
    "transcript",
    "lnc_RNA",
    "ncRNA",
    "miRNA",
    "snRNA",
    "snoRNA",
    "rRNA",
    "scRNA",
    "tRNA",
    "pseudogenic_transcript",
}


def attributes(raw: str) -> dict[str, str]:
    result = {}
    for item in raw.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def transcript_id(value: str | None) -> str | None:
    if not value:
        return None
    value = value.split(",", 1)[0]
    return value.removeprefix("transcript:") if value.startswith("transcript:") else None


def score(value: str, salt: str = "") -> int:
    return int.from_bytes(hashlib.sha256((f"{salt}\t{value}" if salt else value).encode()).digest()[:8], "big")


@dataclass
class Transcript:
    id: str
    chrom: str
    start: int
    end: int
    strand: str
    biotype: str
    feature_type: str
    tags: tuple[str, ...] = ()
    exons: list[tuple[int, int, int]] = field(default_factory=list)
    cds: list[tuple[int, int, int]] = field(default_factory=list)
    start_codons: list[tuple[int, int]] = field(default_factory=list)
    stop_codons: list[tuple[int, int]] = field(default_factory=list)

    @property
    def kind(self) -> str:
        return "coding" if self.cds else "noncoding"

    @property
    def exon_bin(self) -> str:
        return "0" if not self.exons else "1" if len(self.exons) == 1 else "2-5" if len(self.exons) <= 5 else "6+"

    @property
    def phases(self) -> str:
        values = sorted({phase for _, _, phase in self.cds if phase >= 0})
        return "-".join(map(str, values)) if values else "none"

    @property
    def in_par(self) -> bool:
        intervals = {
            "chrX": ((10001, 2781479), (155701383, 156030895)),
            "chrY": ((10001, 2781479), (56887903, 57217415)),
        }
        return any(self.start <= end and self.end >= start for start, end in intervals.get(self.chrom, ()))


class IndexedFasta:
    def __init__(self, fasta: Path):
        self.handle = fasta.open("rb")
        self.index = {}
        with Path(str(fasta) + ".fai").open(encoding="utf-8") as handle:
            for line in handle:
                name, length, offset, line_bases, line_width, *_ = line.rstrip().split("\t")
                self.index[name] = tuple(map(int, (length, offset, line_bases, line_width)))

    def base(self, chrom: str, position: int) -> str | None:
        item = self.index.get(chrom)
        if not item:
            return None
        length, offset, line_bases, line_width = item
        if position < 1 or position > length:
            return None
        zero = position - 1
        byte = offset + (zero // line_bases) * line_width + zero % line_bases
        self.handle.seek(byte)
        value = self.handle.read(1).decode("ascii").upper()
        return value if value in "ACGT" else None

    def sequence(self, chrom: str, start: int, end: int) -> str | None:
        values = [self.base(chrom, position) for position in range(start, end + 1)]
        return "".join(values) if all(values) else None


def read_candidates(gff: Path, modulus: int, expanded: bool, salt: str = "") -> list[Transcript]:
    selected: dict[str, Transcript] = {}
    opener = gzip.open if gff.suffix == ".gz" else open
    with opener(gff, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip().split("\t")
            if len(columns) != 9 or columns[0] not in PRIMARY:
                continue
            chrom, _, feature, start, end, _, strand, phase, raw_attrs = columns
            attrs = attributes(raw_attrs)
            if feature in TRANSCRIPT_FEATURES and (expanded or feature != "pseudogenic_transcript"):
                tid = transcript_id(attrs.get("ID"))
                if tid and (chrom in {"Y", "MT"} or score(tid, salt) % modulus == 0):
                    selected[tid] = Transcript(
                        tid,
                        "chrM" if chrom == "MT" else f"chr{chrom}",
                        int(start),
                        int(end),
                        strand,
                        attrs.get("biotype", feature),
                        feature,
                        tuple(sorted(filter(None, attrs.get("tag", "").split(",")))),
                    )
                continue
            tid = transcript_id(attrs.get("Parent"))
            transcript = selected.get(tid or "")
            if transcript is None:
                continue
            span = (int(start), int(end))
            if feature == "exon":
                transcript.exons.append((*span, int(attrs.get("rank", "0"))))
            elif feature == "CDS":
                transcript.cds.append((*span, int(phase) if phase in {"0", "1", "2"} else -1))
            elif feature == "start_codon":
                transcript.start_codons.append(span)
            elif feature == "stop_codon":
                transcript.stop_codons.append(span)
    return [item for item in selected.values() if item.exons or expanded]


def choose_transcripts(candidates: list[Transcript], expanded: bool, salt: str = "") -> list[Transcript]:
    if expanded:
        return choose_expanded_transcripts(candidates, salt)
    quotas = {
        ("coding", "+"): 8,
        ("coding", "-"): 8,
        ("noncoding", "+"): 4,
        ("noncoding", "-"): 4,
    }
    chosen = []
    for bucket, quota in quotas.items():
        members = [item for item in candidates if (item.kind, item.strand) == bucket]
        groups: dict[tuple[str, str, str], list[Transcript]] = defaultdict(list)
        for item in members:
            groups[(item.chrom, item.exon_bin, item.phases)].append(item)
        for values in groups.values():
            values.sort(key=lambda item: score(item.id, salt))
        keys = sorted(groups, key=lambda value: score("|".join(value), salt))
        while len([item for item in chosen if (item.kind, item.strand) == bucket]) < quota:
            progressed = False
            for key in keys:
                if groups[key]:
                    chosen.append(groups[key].pop(0))
                    progressed = True
                    if len([item for item in chosen if (item.kind, item.strand) == bucket]) == quota:
                        break
            if not progressed:
                raise RuntimeError(f"not enough transcripts for {bucket}: wanted {quota}")
    for chrom in ("chrY", "chrM"):
        additions = sorted(
            (item for item in candidates if item.chrom == chrom and item not in chosen),
            key=lambda item: (item.kind != "coding", score(item.id, salt)),
        )[:2]
        chosen.extend(additions)
    return chosen


def choose_expanded_transcripts(candidates: list[Transcript], salt: str = "") -> list[Transcript]:
    quotas = {
        ("coding", "+"): 50,
        ("coding", "-"): 50,
        ("noncoding", "+"): 40,
        ("noncoding", "-"): 40,
    }
    selectors = {
        "mane_select": lambda item: "MANE_Select" in item.tags,
        "mane_plus_clinical": lambda item: "MANE_Plus_Clinical" in item.tags,
        "nmd": lambda item: item.biotype == "nonsense_mediated_decay",
        "pseudogene": lambda item: "pseudogene" in item.biotype,
        "lncrna": lambda item: item.biotype == "lncRNA",
        "small_rna": lambda item: item.biotype in {"miRNA", "snRNA", "snoRNA", "rRNA", "tRNA"},
        "partial_cds_proxy": lambda item: bool(item.cds) and (not item.start_codons or not item.stop_codons),
        "par": lambda item: item.in_par,
        "x": lambda item: item.chrom == "chrX",
        "y": lambda item: item.chrom == "chrY",
        "mitochondrial": lambda item: item.chrom == "chrM",
        "single_exon": lambda item: len(item.exons) == 1,
        "many_exons": lambda item: len(item.exons) >= 20,
        "phase_1": lambda item: any(phase == 1 for _, _, phase in item.cds),
        "phase_2": lambda item: any(phase == 2 for _, _, phase in item.cds),
    }
    selectors.update(
        {
            f"chromosome_{chrom}": lambda item, chrom=chrom: item.chrom == chrom
            for chrom in [*(f"chr{index}" for index in range(1, 23)), "chrX", "chrY", "chrM"]
        }
    )
    chosen = []
    for label, matches in selectors.items():
        options = sorted((item for item in candidates if matches(item)), key=lambda item: score(f"{label}|{item.id}", salt))
        if not options:
            raise RuntimeError(f"no transcript satisfies required stratum {label}")
        if options[0] not in chosen:
            chosen.append(options[0])

    for bucket, quota in quotas.items():
        members = [item for item in candidates if (item.kind, item.strand) == bucket and item not in chosen]
        groups: dict[tuple[str, str, str, str], list[Transcript]] = defaultdict(list)
        for item in members:
            groups[(item.chrom, item.exon_bin, item.phases, item.biotype)].append(item)
        for values in groups.values():
            values.sort(key=lambda item: score(item.id, salt))
        keys = sorted(groups, key=lambda value: score("|".join(value), salt))
        while len([item for item in chosen if (item.kind, item.strand) == bucket]) < quota:
            progressed = False
            for key in keys:
                if groups[key]:
                    chosen.append(groups[key].pop(0))
                    progressed = True
                    if len([item for item in chosen if (item.kind, item.strand) == bucket]) == quota:
                        break
            if not progressed:
                raise RuntimeError(f"not enough transcripts for {bucket}: wanted {quota}")
    if len(chosen) != sum(quotas.values()):
        raise RuntimeError(f"expanded selection produced {len(chosen)} transcripts")
    return chosen


def representative_spans(spans: list[tuple[int, ...]]) -> list[tuple[int, ...]]:
    if len(spans) <= 3:
        return spans
    ordered = sorted(spans)
    return [ordered[0], ordered[len(ordered) // 2], ordered[-1]]


def boundary_points(transcript: Transcript) -> dict[int, set[str]]:
    points: dict[int, set[str]] = defaultdict(set)
    points[transcript.start].add("transcript_start")
    points[transcript.end].add("transcript_end")
    for start, end, *_ in representative_spans(transcript.exons):
        points[start].add("exon_start")
        points[end].add("exon_end")
    for start, end, *_ in representative_spans(transcript.cds):
        points[start].add("cds_segment_start")
        points[end].add("cds_segment_end")
    if transcript.cds:
        low = min(start for start, _, _ in transcript.cds)
        high = max(end for _, end, _ in transcript.cds)
        coding_start, coding_stop = (low, high) if transcript.strand == "+" else (high, low)
        points[coding_start].add("coding_start_boundary")
        points[coding_stop].add("coding_stop_boundary")
    for label, spans in (("start_codon", transcript.start_codons), ("stop_codon", transcript.stop_codons)):
        for start, end in spans:
            points[start].add(f"{label}_start")
            points[end].add(f"{label}_end")
    return points


def add_variant(variants: dict, chrom: str, pos: int, ref: str | None, alt: str | None, note: dict):
    if not ref or not alt or ref == alt:
        return
    key = (chrom, pos, ref, alt)
    variants.setdefault(key, []).append(note)


def build_variants(transcripts: list[Transcript], fasta: IndexedFasta, coding_only: bool = False) -> dict:
    variants: dict[tuple[str, int, str, str], list[dict]] = {}
    alternate = {"A": "C", "C": "G", "G": "T", "T": "A"}
    for transcript in transcripts:
        for boundary, labels in boundary_points(transcript).items():
            if coding_only and not labels & {"coding_start_boundary", "coding_stop_boundary"}:
                continue
            common = {
                "transcript": transcript.id,
                "strand": transcript.strand,
                "biotype": transcript.biotype,
                "kind": transcript.kind,
                "exonBin": transcript.exon_bin,
                "phases": transcript.phases,
                "boundaries": sorted(labels),
            }
            for delta in (-2, -1, 0, 1, 2):
                position = boundary + delta
                ref = fasta.base(transcript.chrom, position)
                add_variant(
                    variants,
                    transcript.chrom,
                    position,
                    ref,
                    alternate.get(ref or ""),
                    {**common, "shape": "snv", "offset": delta},
                )
            ref = fasta.base(transcript.chrom, boundary)
            add_variant(
                variants,
                transcript.chrom,
                boundary,
                ref,
                (ref or "") + alternate.get(ref or "", ""),
                {**common, "shape": "insertion", "offset": 0},
            )
            deletion = fasta.sequence(transcript.chrom, boundary - 1, boundary)
            add_variant(
                variants,
                transcript.chrom,
                boundary - 1,
                deletion,
                deletion[:1] if deletion else None,
                {**common, "shape": "deletion", "offset": 0},
            )
            ref2 = fasta.sequence(transcript.chrom, boundary, boundary + 1)
            alt2 = "".join(alternate.get(base, "") for base in ref2 or "")
            add_variant(
                variants,
                transcript.chrom,
                boundary,
                ref2,
                alt2,
                {**common, "shape": "mnv2", "offset": 0},
            )
    return variants


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--sample-modulus", type=int, default=128)
    parser.add_argument("--expanded", action="store_true")
    parser.add_argument("--risk-features", action="store_true",
                        help="All primary GFF transcripts with a partial CDS or non-ATG first codon; coding boundaries only")
    parser.add_argument("--salt", default="")
    parser.add_argument("--exclude-vcf", type=Path, action="append", default=[])
    args = parser.parse_args()
    if args.vcf.exists() or args.manifest.exists():
        parser.error("Refusing to overwrite an existing corpus or manifest")

    if args.risk_features:
        args.sample_modulus = 1
        args.expanded = True
    candidates = read_candidates(args.gff, args.sample_modulus, args.expanded, args.salt)
    fasta = IndexedFasta(args.fasta)
    if args.risk_features:
        transcripts = []
        for item in candidates:
            spans = sorted(item.cds, reverse=item.strand == "-")
            if not spans:
                continue
            positions = []
            for start, end, _ in spans:
                positions.extend(list(range(start, min(end + 1, start + 3))) if item.strand == "+"
                                 else list(range(end, max(start - 1, end - 3), -1)))
                if len(positions) >= 3:
                    break
            bases = [fasta.base(item.chrom, pos) for pos in positions[:3]]
            first = "".join(base or "N" for base in bases)
            if item.strand == "-":
                first = first.translate(str.maketrans("ACGT", "TGCA"))
            if sum(end - start + 1 for start, end, _ in spans) % 3 or first != "ATG":
                transcripts.append(item)
    else:
        transcripts = choose_transcripts(candidates, args.expanded, args.salt)
    variants = build_variants(transcripts, fasta, args.risk_features)
    excluded = set()
    for path in args.exclude_vcf:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith("#"):
                    chrom, pos, _, ref, alts, *_ = line.rstrip().split("\t")
                    chrom = chrom.removeprefix("chr")
                    chrom = "chrM" if chrom in {"M", "MT"} else f"chr{chrom}"
                    excluded.update((chrom, int(pos), ref.upper(), alt.upper()) for alt in alts.split(","))
    excluded_count = len(variants.keys() & excluded)
    variants = {key: value for key, value in variants.items() if key not in excluded}

    ordered = sorted(variants, key=lambda item: (PRIMARY_ORDER.get(item[0], 99), item[1], item[2], item[3]))
    with args.vcf.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##reference=GRCh38\n")
        handle.write("##source=annocat-independent-transcript-geometry-20260904\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for index, (chrom, pos, ref, alt) in enumerate(ordered, 1):
            handle.write(f"{chrom}\t{pos}\tBOUNDARY{index}\t{ref}\t{alt}\t.\tPASS\t.\n")

    manifest = {
        "schemaVersion": 1,
        "corpusId": "ensembl-115-boundary-expanded" if args.expanded else "ensembl-115-boundary",
        "role": "targeted-development-discovery" if args.risk_features else "qualification-discovery" if args.expanded else "release-regression",
        "assembly": "GRCh38",
        "riskFeatures": args.risk_features,
        "sources": {
            "gff3": {
                "url": "https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz",
                "sha256": "1e553efa8496d662e7264061a5cecf3001eb9a1157aaa66d80cd7ac35841509c",
            },
            "reference": {
                "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
                "archiveSha256": "fb4243ebb014caf27111f24dd62b7ce42160f28581da6f8fcd6cba5977778d02",
                "preparedSha256": "9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7",
            },
        },
        "generator": {
            "path": "scripts/build-vep-boundary-corpus.py",
            "sha256": file_sha256(Path(__file__)),
            "arguments": sys.argv[1:],
            "command": "python scripts/build-vep-boundary-corpus.py --gff <pinned-gff3> --fasta <pinned-fasta> --vcf <output-vcf> --manifest <output-manifest> --sample-modulus "
            + str(args.sample_modulus)
            + (" --expanded" if args.expanded else "")
            + (" --risk-features" if args.risk_features else "")
            + (" --salt <selection.salt>" if args.salt else "")
            + (" --exclude-vcf <each-selection.excludedVcfInputs.path>" if args.exclude_vcf else ""),
        },
        "selection": {
            "salt": args.salt,
            "excludedVcfInputs": [{"path": str(path), "sha256": file_sha256(path)} for path in args.exclude_vcf],
            "excludedPriorUseVariants": excluded_count,
            "sampleModulus": args.sample_modulus,
            "candidateTranscripts": len(candidates),
            "selectedTranscripts": len(transcripts),
            "algorithm": ("All primary GFF CDS models with total CDS length modulo 3 nonzero or non-ATG first genomic CDS triplet; coding boundaries only"
                          if args.risk_features else "SHA-256 transcript sample, then balanced coding/strand/geometry strata"),
            "expanded": args.expanded,
            "bucketCounts": dict(
                sorted(
                    Counter(f"{item.kind}/{item.strand}" for item in transcripts).items()
                )
            ),
            "chromosomeCounts": dict(
                sorted(Counter(item.chrom for item in transcripts).items())
            ),
            "biotypeCounts": dict(
                sorted(Counter(item.biotype for item in transcripts).items())
            ),
            "parTranscripts": sum(item.in_par for item in transcripts),
            "maneSelectTranscripts": sum("MANE_Select" in item.tags for item in transcripts),
            "manePlusClinicalTranscripts": sum(
                "MANE_Plus_Clinical" in item.tags for item in transcripts
            ),
            "partialCdsProxyTranscripts": sum(
                bool(item.cds) and (not item.start_codons or not item.stop_codons)
                for item in transcripts
            ),
        },
        "transcripts": [
            {
                "id": item.id,
                "chromosome": item.chrom,
                "start": item.start,
                "end": item.end,
                "strand": item.strand,
                "biotype": item.biotype,
                "featureType": item.feature_type,
                "kind": item.kind,
                "exonCount": len(item.exons),
                "exonBin": item.exon_bin,
                "cdsSegmentCount": len(item.cds),
                "phases": item.phases,
                "tags": list(item.tags),
                "inPar": item.in_par,
                "partialCdsProxy": bool(item.cds)
                and (not item.start_codons or not item.stop_codons),
            }
            for item in transcripts
        ],
        "variantCount": len(ordered),
        "variants": [
            {
                "chromosome": chrom,
                "position": pos,
                "reference": ref,
                "alternate": alt,
                "probes": variants[(chrom, pos, ref, alt)],
            }
            for chrom, pos, ref, alt in ordered
        ],
        "output": {
            "path": str(args.vcf).replace("\\", "/"),
            "records": len(ordered),
            "alleles": len(ordered),
            "sha256": file_sha256(args.vcf),
        },
    }
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"transcripts": len(transcripts), "variants": len(ordered)}, indent=2))


PRIMARY_ORDER = {f"chr{i}": i for i in range(1, 23)} | {"chrX": 23, "chrY": 24, "chrM": 25}


if __name__ == "__main__":
    main()
