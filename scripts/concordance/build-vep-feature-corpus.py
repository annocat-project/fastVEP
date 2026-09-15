#!/usr/bin/env python3
"""Source-guided VEP discovery: transcript feature combinations and dense short edits."""
import argparse
import importlib.util
import itertools
import json
import sys
from collections import Counter
from pathlib import Path

helper_path = Path(__file__).with_name("build-vep-boundary-corpus.py")
spec = importlib.util.spec_from_file_location("vep_boundary", helper_path)
geometry = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = geometry
spec.loader.exec_module(geometry)


def features(transcript, fasta):
    spans = sorted(transcript.cds, reverse=transcript.strand == "-")
    if not spans:
        return None
    low = min(start for start, _, _ in spans)
    high = max(end for _, end, _ in spans)
    positions = []
    for start, end, _ in spans:
        positions.extend(range(start, min(end + 1, start + 3)) if transcript.strand == "+"
                         else range(end, max(start - 1, end - 3), -1))
        if len(positions) >= 3:
            break
    first = "".join(fasta.base(transcript.chrom, pos) or "N" for pos in positions[:3])
    if transcript.strand == "-":
        first = first.translate(str.maketrans("ACGT", "TGCA"))
    utr5, utr3 = (transcript.start < low, high < transcript.end)
    if transcript.strand == "-":
        utr5, utr3 = utr3, utr5
    return (transcript.strand, utr5, utr3, spans[0][2],
            sum(end - start + 1 for start, end, _ in spans) % 3, first == "ATG")


def edits(fasta, chrom, position):
    """All single-base substitutions and 1/2-base insertions; longer repeat edits."""
    base = fasta.base(chrom, position)
    if not base:
        return
    for alt in "ACGT":
        if alt != base:
            yield position, base, alt, "snv"
    for length in (1, 2):
        for letters in itertools.product("ACGT", repeat=length):
            yield position, base, base + "".join(letters), f"ins-{length}"
    for length in range(1, 13):
        deleted = fasta.sequence(chrom, position, position + length)
        if deleted:
            yield position, deleted, base, f"del-{length}"
        for side, start, end in [("left", position - length + 1, position),
                                 ("right", position + 1, position + length)]:
            repeat = fasta.sequence(chrom, start, end)
            if repeat:
                yield position, base, base + repeat, f"repeat-{side}-{length}"
    for length in (3, 5, 8):
        for letter in "ACGT":
            yield position, base, base + letter * length, f"homopolymer-{length}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--salt", required=True)
    parser.add_argument("--exclude-vcf", type=Path, action="append", default=[])
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error("Refusing to overwrite a frozen corpus")
    fasta = geometry.IndexedFasta(args.fasta)
    groups = {}
    population = Counter()
    for transcript in geometry.read_candidates(args.gff, 1, True, args.salt):
        group = features(transcript, fasta)
        if group is None:
            continue
        population[group] += 1
        if group not in groups or geometry.score(transcript.id, args.salt) < geometry.score(groups[group].id, args.salt):
            groups[group] = transcript
    excluded = set()
    for path in args.exclude_vcf:
        with path.open(encoding="utf-8") as source:
            for line in source:
                if line.startswith("#"):
                    continue
                chrom, pos, _, ref, alts, *_ = line.rstrip().split("\t")
                chrom = chrom.removeprefix("chr")
                chrom = "chrM" if chrom in {"M", "MT"} else "chr" + chrom
                excluded.update((chrom, int(pos), ref.upper(), alt.upper()) for alt in alts.split(","))
    variants = {}
    skipped = set()
    selected = []
    for group, transcript in sorted(groups.items()):
        start = min(start for start, _, _ in transcript.cds)
        end = max(end for _, end, _ in transcript.cds)
        if transcript.strand == "-":
            start, end = end, start
        direction = 1 if transcript.strand == "+" else -1
        selected.append({"transcript": transcript.id, "stratum": group,
                         "population": population[group], "chromosome": transcript.chrom,
                         "codingStart": start, "codingEnd": end})
        points = [("coding-start", start + direction * delta, delta)
                  for delta in (-2, -1, 0, 1, 2, 3, 5, 8, 11)]
        points += [("coding-end", end + direction * delta, delta)
                   for delta in (-11, -8, -5, -3, -2, -1, 0, 1, 2)]
        # Query selection is independent of HGVS shifting. Include both sides
        # of the configured 5kb boundary and several padded spellings.
        points += [("query-edge", edge + delta, delta)
                   for edge in (transcript.start - 5000, transcript.end + 5000)
                   for delta in (-2, -1, 0, 1, 2)]
        for endpoint, position, delta in points:
            for pos, ref, alt, shape in edits(fasta, transcript.chrom, position):
                if endpoint == "query-edge" and shape not in {"snv", "del-1", "del-2", "ins-1"}:
                    continue
                forms = [(pos, ref, alt, shape)]
                if len(ref) != len(alt) and (endpoint == "query-edge" or shape in {"ins-1", "del-1"}):
                    for padding in (1, 3):
                        prefix = fasta.sequence(transcript.chrom, pos - padding, pos - 1)
                        if prefix:
                            forms.append((pos - padding, prefix + ref, prefix + alt, shape + f"-pad-{padding}"))
                for pos, ref, alt, shape in forms:
                    key = (transcript.chrom, pos, ref, alt)
                    if key in excluded:
                        skipped.add(key)
                        continue
                    variants.setdefault(key, set()).add((transcript.id, endpoint, delta, shape))
    args.output_dir.mkdir(parents=True)
    header = "##fileformat=VCFv4.2\n##reference=GRCh38\n##source=annocat-vep115-feature-discovery\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    ordered = sorted(variants, key=lambda k: (geometry.PRIMARY_ORDER.get(k[0], 99), *k[1:]))
    outputs = []
    with (args.output_dir / "probes.jsonl").open("w", encoding="utf-8", newline="\n") as probes:
        for offset in range(0, len(ordered), 50000):
            name = f"features-{offset // 50000 + 1:03}"
            path = args.output_dir / (name + ".vcf")
            batch = ordered[offset:offset + 50000]
            with path.open("w", encoding="utf-8", newline="\n") as out:
                out.write(header)
                for index, key in enumerate(batch, offset + 1):
                    chrom, pos, ref, alt = key
                    out.write(f"{chrom}\t{pos}\tFEATURE{index}\t{ref}\t{alt}\t.\tPASS\t.\n")
                    probes.write(json.dumps({"variant": key, "probes": sorted(variants[key])}, separators=(",", ":")) + "\n")
            outputs.append({"corpus": name, "path": path.name, "records": len(batch),
                            "sha256": geometry.file_sha256(path), "priorUseOverlap": 0})
    manifest = {"qualificationEligible": False, "purpose": __doc__, "salt": args.salt,
                "stratumFields": ["strand", "fivePrimeUtr", "threePrimeUtr", "firstCdsPhase", "cdsLengthModulo3", "firstTripletATG"],
                "selection": "Lowest salted transcript hash in every occupied feature combination; GFF/FASTA only, no candidate output.",
                "selected": selected, "records": len(ordered), "excludedPriorUseVariants": len(skipped),
                "pins": {str(p.resolve()): geometry.file_sha256(p) for p in [Path(__file__), helper_path, args.gff, args.fasta, *args.exclude_vcf]},
                "outputs": outputs, "probesSha256": geometry.file_sha256(args.output_dir / "probes.jsonl")}
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"strata": len(groups), "records": len(ordered), "shards": len(outputs)}))


if __name__ == "__main__":
    main()
