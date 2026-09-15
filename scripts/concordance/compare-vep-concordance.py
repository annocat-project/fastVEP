#!/usr/bin/env python3
"""Fail unless fastVEP and an Ensembl VEP oracle emit the same consequence rows."""

import argparse
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
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

PRODUCTION_FIELDS = (
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
    "Existing_variation",
    "REF_ALLELE",
    "UPLOADED_ALLELE",
    "DISTANCE",
    "STRAND",
    "FLAGS",
    "CANONICAL",
    "SYMBOL_SOURCE",
    "HGNC_ID",
    "MANE",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "TSL",
    "APPRIS",
    "CCDS",
    "ENSP",
    "SOURCE",
    "HGVS_OFFSET",
    "SIFT",
    "PolyPhen",
    "AF",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",
    "ACMG",
    "ACMG_CRITERIA",
)

# The archive database adds HGNC provenance absent from Ensembl GFF3, does not
# expose VEP's source label, and fastVEP's GFF3 loader does not retain APPRIS.
# Everything below has a direct representation in both outputs.
REST_FIELDS = (
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
    "DISTANCE",
    "STRAND",
    "FLAGS",
    "CANONICAL",
    "MANE",
    "TSL",
    "CCDS",
    "ENSP",
    "HGVS_OFFSET",
)
ALL_FIELDS = PRODUCTION_FIELDS
NORMALIZERS = {
    "unordered-flags",
    "empty-allele-dash",
    "hgvsp-gff-protein-version",
    "uploaded-allele-delimiter",
}


def open_text(path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalized_info(value):
    return tuple(
        sorted(
            item
            for item in value.split(";")
            if item and item != "." and not item.startswith("CSQ=")
        )
    )


def record_value(columns):
    return tuple(columns[:7]) + (normalized_info(columns[7]),) + tuple(columns[8:])


def parse_input_vcf(path):
    records = Counter()
    samples = None
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("#CHROM"):
                columns = line.rstrip("\n\r").split("\t")
                samples = tuple(columns[9:]) if len(columns) > 8 else ()
                continue
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n\r").split("\t")
            if len(columns) < 8:
                raise ValueError(f"{path}:{line_number}: invalid VCF row")
            records[record_value(columns)] += 1
    if samples is None:
        raise ValueError(f"{path}: #CHROM header is missing")
    if not records:
        raise ValueError(f"{path}: no VCF records")
    return records, samples


def normalize_field(value, normalizer):
    if normalizer == "unordered-flags":
        return "&".join(sorted(value.split("&")))
    if normalizer == "empty-allele-dash":
        return "" if value == "-" else value
    if normalizer == "hgvsp-gff-protein-version":
        return re.sub(r"^([^:]+?)\.\d+(:p\.)", r"\1\2", value)
    if normalizer == "uploaded-allele-delimiter":
        return value.replace("&", "/")
    raise ValueError(f"unsupported field normalizer {normalizer!r}")


def parse_vcf(path, fields=FIELDS, normalizers=None):
    normalizers = normalizers or {}
    csq_fields = None
    variants = Counter()
    annotations = Counter()
    records = Counter()
    samples = None
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("##INFO=<ID=CSQ"):
                if csq_fields is not None:
                    raise ValueError(f"{path}:{line_number}: duplicate CSQ header")
                match = re.search(r'Format: ([^">]+)', line)
                if not match:
                    raise ValueError(f"{path}:{line_number}: CSQ format is missing")
                csq_fields = match.group(1).split("|")
                if len(csq_fields) != len(set(csq_fields)):
                    raise ValueError(f"{path}:{line_number}: duplicate CSQ field name")
                continue
            if line.startswith("#CHROM"):
                columns = line.rstrip("\n\r").split("\t")
                samples = tuple(columns[9:]) if len(columns) > 8 else ()
                continue
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n\r").split("\t")
            if len(columns) < 8:
                raise ValueError(f"{path}:{line_number}: invalid VCF row")
            if csq_fields is None:
                raise ValueError(f"{path}: CSQ header is missing")
            key = tuple(columns[index] for index in (0, 1, 3, 4))
            variants[key] += 1
            records[record_value(columns)] += 1
            csq_values = [
                item[4:] for item in columns[7].split(";") if item.startswith("CSQ=")
            ]
            if len(csq_values) > 1:
                raise ValueError(f"{path}:{line_number}: duplicate CSQ INFO value")
            csq_value = csq_values[0] if csq_values else ""
            for encoded in filter(None, csq_value.split(",")):
                values = encoded.split("|")
                if len(values) != len(csq_fields):
                    raise ValueError(
                        f"{path}:{line_number}: CSQ row has {len(values)} values; "
                        f"header declares {len(csq_fields)}"
                    )
                row = dict(zip(csq_fields, values))
                for field, normalizer in normalizers.items():
                    if field in row:
                        row[field] = normalize_field(row[field], normalizer)
                annotations[(key, tuple(row.get(field, "") for field in fields))] += 1
    if csq_fields is None:
        raise ValueError(f"{path}: CSQ header is missing")
    if samples is None:
        raise ValueError(f"{path}: #CHROM header is missing")
    return variants, annotations, set(csq_fields), records, samples


def csq_escape(value):
    return (
        str(value)
        .replace(",", "&")
        .replace("|", "&")
        .replace(";", "%3B")
        .replace("=", "%3D")
    )


def joined(value):
    if value is None:
        return ""
    if isinstance(value, list):
        return "&".join(csq_escape(item) for item in value)
    return csq_escape(value)


def position(row, prefix):
    start = row.get(f"{prefix}_start")
    end = row.get(f"{prefix}_end")
    if start is None and end is None:
        return ""
    if start is None:
        return f"?-{end}"
    if end is None:
        return f"{start}-?"
    if start == end:
        return str(start)
    return f"{min(start, end)}-{max(start, end)}"


def rest_row(row, feature_type, feature_key, fields=REST_FIELDS):
    values = {
        "Allele": joined(row.get("variant_allele")),
        "Consequence": joined(row.get("consequence_terms")),
        "IMPACT": joined(row.get("impact")),
        "SYMBOL": joined(row.get("gene_symbol")),
        "Gene": joined(row.get("gene_id")),
        "Feature_type": feature_type,
        "Feature": joined(row.get(feature_key)),
        "BIOTYPE": joined(row.get("biotype")),
        "EXON": joined(row.get("exon")),
        "INTRON": joined(row.get("intron")),
        "HGVSc": joined(row.get("hgvsc")),
        "HGVSp": joined(row.get("hgvsp")),
        "cDNA_position": position(row, "cdna"),
        "CDS_position": position(row, "cds"),
        "Protein_position": position(row, "protein"),
        "Amino_acids": joined(row.get("amino_acids")),
        "Codons": joined(row.get("codons")),
        "DISTANCE": joined(row.get("distance")),
        "STRAND": joined(row.get("strand")),
        "FLAGS": joined(row.get("flags")),
        "CANONICAL": "YES" if row.get("canonical") else "",
        "MANE": joined(row.get("mane")),
        "TSL": joined(row.get("tsl")),
        "CCDS": joined(row.get("ccds")),
        "ENSP": joined(row.get("protein_id")),
        "HGVS_OFFSET": joined(row.get("hgvs_offset")),
    }
    return tuple(values[field] for field in fields)


def parse_rest(path, fields=REST_FIELDS):
    document = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(document, list) or not document:
        raise ValueError(f"{path}: VEP REST response must be a non-empty array")
    variants = Counter()
    annotations = Counter()
    collections = (
        ("transcript_consequences", "Transcript", "transcript_id"),
        ("intergenic_consequences", "Intergenic", ""),
    )
    for index, result in enumerate(document, 1):
        if not isinstance(result, dict):
            raise ValueError(f"{path}: response item {index} is not an object")
        if result.get("assembly_name") != "GRCh38":
            raise ValueError(f"{path}: response item {index} is not GRCh38")
        columns = str(result.get("input", "")).split()
        if len(columns) < 5:
            raise ValueError(f"{path}: response item {index} has no VCF input identity")
        key = tuple(columns[item] for item in (0, 1, 3, 4))
        variants[key] += 1
        for collection, feature_type, feature_key in collections:
            rows = result.get(collection, [])
            if not isinstance(rows, list):
                raise ValueError(f"{path}: response item {index} has invalid {collection}")
            for row in rows:
                if not isinstance(row, dict):
                    raise ValueError(f"{path}: {collection} contains a non-object value")
                annotations[(key, rest_row(row, feature_type, feature_key, fields))] += 1
    return variants, annotations


def allowed_extras(document, path):
    allowed = {}
    for item in document.get("allowedExtraIdentities", []):
        variant = item.get("variant")
        if not isinstance(variant, list) or len(variant) != 4 or not item.get("reason"):
            raise ValueError(f"{path}: invalid allowed extra identity")
        identity = (
            tuple(str(value) for value in variant),
            str(item.get("allele", "")),
            str(item.get("featureType", "")),
            str(item.get("feature", "")),
        )
        if not all(identity[1:]) or identity in allowed:
            raise ValueError(f"{path}: duplicate or incomplete allowed extra identity")
        allowed[identity] = int(item.get("count", 1))
        if allowed[identity] < 1:
            raise ValueError(f"{path}: allowed extra count must be positive")
    return allowed


def allowed_field_differences(document, path, fields):
    source = document.get("allowedFieldDifferences", [])
    source_report = None
    default_input_sha = None
    default_reason = None
    if isinstance(source, str):
        source_path = path.parent / source
        source_document = json.loads(source_path.read_text(encoding="utf-8"))
        if source_document.get("schemaVersion") != 1:
            raise ValueError(f"{source_path}: unsupported field-difference schema")
        source = source_document.get("differences")
        default_input_sha = source_document.get("inputSha256")
        default_reason = source_document.get("reason")
        source_report = {"path": str(source_path), "sha256": sha256(source_path)}
    if not isinstance(source, list):
        raise ValueError(f"{path}: allowedFieldDifferences must be a list or path")

    allowed = []
    seen = set()
    for item in source:
        variant = item.get("variant") if isinstance(item, dict) else None
        input_sha = item.get("inputSha256", default_input_sha) if isinstance(item, dict) else None
        field = item.get("field", "") if isinstance(item, dict) else ""
        reason = item.get("reason", default_reason) if isinstance(item, dict) else None
        identity = (
            tuple(str(value) for value in variant) if isinstance(variant, list) else (),
            str(item.get("allele", "")) if isinstance(item, dict) else "",
            str(item.get("featureType", "")) if isinstance(item, dict) else "",
            str(item.get("feature", "")) if isinstance(item, dict) else "",
        )
        if (
            len(identity[0]) != 4
            or not all(identity[1:])
            or field not in fields
            or not isinstance(input_sha, str)
            or not re.fullmatch(r"[0-9a-f]{64}", input_sha)
            or not reason
            or not isinstance(item.get("candidate"), str)
            or not isinstance(item.get("oracle"), str)
            or item["candidate"] == item["oracle"]
        ):
            raise ValueError(f"{path}: invalid allowed field difference")
        key = (input_sha, identity, field)
        if key in seen:
            raise ValueError(f"{path}: duplicate allowed field difference")
        seen.add(key)
        count = int(item.get("count", 1))
        if count < 1:
            raise ValueError(f"{path}: allowed field-difference count must be positive")
        allowed.append(
            {
                "inputSha256": input_sha,
                "identity": identity,
                "field": field,
                "candidate": item["candidate"],
                "oracle": item["oracle"],
                "reason": reason,
                "count": count,
            }
        )
    return allowed, source_report


def load_contract(path, default_fields):
    if path is None:
        return {
            "fields": tuple(default_fields),
            "candidateFields": set(default_fields),
            "oracleFields": set(default_fields),
            "candidateEmptyFields": set(),
            "normalizers": {},
            "allowedExtras": {},
            "allowedFieldDifferences": [],
            "report": None,
        }
    document = json.loads(path.read_text(encoding="utf-8"))
    schema = document.get("schemaVersion")
    if schema == 1:
        ignored = set()
        for item in document.get("ignoredFields", []):
            field = item.get("field", "")
            if field not in default_fields or not item.get("reason"):
                raise ValueError(f"{path}: invalid ignored field {field!r}")
            ignored.add(field)
        fields = tuple(field for field in default_fields if field not in ignored)
        report = {
            "path": str(path),
            "sha256": sha256(path),
            "schemaVersion": schema,
            "ignoredFields": sorted(ignored),
        }
        candidate_fields = set(fields)
        oracle_fields = set(fields)
        empty_fields = set()
        normalizers = {}
    elif schema == 2:
        entries = document.get("fields")
        if not isinstance(entries, list) or not entries:
            raise ValueError(f"{path}: schema 2 contract requires fields")
        names = []
        dispositions = {}
        empty_fields = set()
        normalizers = {}
        for item in entries:
            name = item.get("name", "") if isinstance(item, dict) else ""
            disposition = item.get("disposition", "") if isinstance(item, dict) else ""
            if name not in PRODUCTION_FIELDS or name in dispositions:
                raise ValueError(f"{path}: invalid or duplicate production field {name!r}")
            if disposition not in {"exact", "input-derived", "source-derived", "excluded"}:
                raise ValueError(f"{path}: invalid disposition for {name!r}")
            if disposition != "exact" and not item.get("reason"):
                raise ValueError(f"{path}: {name!r} requires a disposition reason")
            if item.get("requireCandidateEmpty"):
                empty_fields.add(name)
            normalizer = item.get("normalizer")
            if normalizer:
                if normalizer == "unordered-flags" and name != "FLAGS":
                    raise ValueError(f"{path}: unordered-flags is only valid for FLAGS")
                if normalizer not in NORMALIZERS or disposition not in {
                    "exact",
                    "input-derived",
                }:
                    raise ValueError(f"{path}: invalid normalizer for {name!r}")
                if not item.get("reason"):
                    raise ValueError(f"{path}: normalized field {name!r} requires a reason")
                normalizers[name] = normalizer
            names.append(name)
            dispositions[name] = disposition
        if tuple(names) != PRODUCTION_FIELDS:
            raise ValueError(f"{path}: fields must exactly match fastVEP production order")
        fields = tuple(
            name
            for name in names
            if dispositions[name] in {"exact", "input-derived"}
        )
        if not {"Allele", "Feature_type", "Feature"}.issubset(fields):
            raise ValueError(f"{path}: comparison identity fields must be exact")
        candidate_fields = set(names)
        oracle_fields = set(fields)
        report = {
            "path": str(path),
            "sha256": sha256(path),
            "schemaVersion": schema,
            "fieldDispositions": dispositions,
            "fieldNormalizers": normalizers,
        }
    else:
        raise ValueError(f"{path}: unsupported contract schema")

    allowed = allowed_extras(document, path)
    field_differences, field_difference_source = allowed_field_differences(
        document, path, fields
    )
    report["allowedExtraIdentities"] = len(allowed)
    report["allowedFieldDifferences"] = len(field_differences)
    if field_difference_source is not None:
        report["allowedFieldDifferenceSource"] = field_difference_source
    return {
        "fields": fields,
        "candidateFields": candidate_fields,
        "oracleFields": oracle_fields,
        "candidateEmptyFields": empty_fields,
        "normalizers": normalizers,
        "allowedExtras": allowed,
        "allowedFieldDifferences": field_differences,
        "report": report,
    }


def apply_allowed_extras(candidate, oracle, fields, allowed):
    if not allowed:
        return candidate, 0, []
    indexes = tuple(fields.index(field) for field in ("Allele", "Feature_type", "Feature"))
    adjusted = candidate.copy()
    remaining = dict(allowed)
    applied = 0
    for (key, values), count in candidate.items():
        identity = (key, *(values[index] for index in indexes))
        allowance = remaining.get(identity, 0)
        oracle_count = oracle.get((key, values), 0)
        removable = min(max(count - oracle_count, 0), allowance)
        if removable:
            adjusted[(key, values)] -= removable
            if adjusted[(key, values)] == 0:
                del adjusted[(key, values)]
            remaining[identity] -= removable
            applied += removable
    unused = [
        {
            "variant": list(identity[0]),
            "allele": identity[1],
            "featureType": identity[2],
            "feature": identity[3],
            "count": count,
        }
        for identity, count in remaining.items()
        if count
    ]
    return adjusted, applied, unused


def apply_allowed_field_differences(candidate, oracle, fields, allowed, input_sha):
    if not allowed:
        return candidate, 0, 0, []
    if input_sha is None:
        raise ValueError("allowed field differences require the original input VCF")

    applicable = [item for item in allowed if item["inputSha256"] == input_sha]
    adjusted = candidate.copy()
    candidate_by_variant = defaultdict(dict)
    oracle_by_variant = defaultdict(dict)
    for row in adjusted:
        candidate_by_variant[row[0]][row] = None
    for row in oracle:
        oracle_by_variant[row[0]][row] = None
    applied = 0
    identity_indexes = tuple(
        fields.index(field) for field in ("Allele", "Feature_type", "Feature")
    )
    unused = []
    for item in applicable:
        field_index = fields.index(item["field"])
        variant, allele, feature_type, feature = item["identity"]

        def matches(row, expected):
            key, values = row
            return (
                key == variant
                and tuple(values[index] for index in identity_indexes)
                == (allele, feature_type, feature)
                and values[field_index] == expected
            )

        candidate_rows = [row for row in candidate_by_variant[variant] if matches(row, item["candidate"])]
        oracle_rows = [row for row in oracle_by_variant[variant] if matches(row, item["oracle"])]
        remaining = item["count"]
        if len(candidate_rows) == 1 and len(oracle_rows) == 1:
            source = candidate_rows[0]
            usable = min(remaining, adjusted[source], oracle[oracle_rows[0]])
            if usable:
                values = list(source[1])
                values[field_index] = item["oracle"]
                target = (source[0], tuple(values))
                adjusted[source] -= usable
                if adjusted[source] == 0:
                    del adjusted[source]
                    del candidate_by_variant[variant][source]
                adjusted[target] += usable
                candidate_by_variant[variant][target] = None
                applied += usable
                remaining -= usable
        if remaining:
            unused.append(
                {
                    "variant": list(variant),
                    "allele": allele,
                    "featureType": feature_type,
                    "feature": feature,
                    "field": item["field"],
                    "candidate": item["candidate"],
                    "oracle": item["oracle"],
                    "count": remaining,
                }
            )
    return adjusted, len(applicable), applied, unused


def examples(counter, fields, limit=10):
    rows = []
    for (key, values), count in counter.most_common(limit):
        identity = dict(zip(fields, values))
        rows.append(
            {
                "variant": ":".join(key),
                "allele": identity["Allele"],
                "feature": identity["Feature"],
                "count": count,
                "fields": identity,
            }
        )
    return rows


def annotation_rows(counter, fields):
    rows = []
    for (key, values), count in sorted(counter.items()):
        rows.append(
            {
                "variant": list(key),
                "count": count,
                "fields": dict(zip(fields, values)),
            }
        )
    return rows


def record_rows(counter):
    names = ("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    rows = []
    for values, count in sorted(counter.items(), key=lambda item: repr(item[0])):
        fixed = dict(zip(names, values[:8]))
        fixed["INFO"] = list(fixed["INFO"])
        if len(values) > 8:
            fixed["FORMAT"] = values[8]
            fixed["samples"] = list(values[9:])
        rows.append({"count": count, "record": fixed})
    return rows


def field_mismatches(candidate, oracle, fields):
    identity_indexes = tuple(fields.index(field) for field in ("Allele", "Feature_type", "Feature"))

    def indexed(rows):
        result = {}
        for (key, values), count in rows.items():
            identity = (key, *(values[index] for index in identity_indexes))
            result.setdefault(identity, []).extend([values] * count)
        return result

    candidate_rows = indexed(candidate)
    oracle_rows = indexed(oracle)
    shared = candidate_rows.keys() & oracle_rows.keys()
    mismatches = Counter()
    ambiguous = 0
    for identity in shared:
        left = candidate_rows[identity]
        right = oracle_rows[identity]
        if len(left) != 1 or len(right) != 1:
            ambiguous += 1
        if left == right:
            continue
        for index, field in enumerate(fields):
            if Counter(row[index] for row in left) != Counter(row[index] for row in right):
                mismatches[field] += 1
    return {
        "candidateIdentities": len(candidate_rows),
        "oracleIdentities": len(oracle_rows),
        "sharedIdentities": len(shared),
        "comparedIdentities": len(shared),
        "ambiguousMultirowIdentities": ambiguous,
        "missingIdentities": len(oracle_rows.keys() - candidate_rows.keys()),
        "extraIdentities": len(candidate_rows.keys() - oracle_rows.keys()),
        "mismatchesByField": dict(sorted(mismatches.items())),
    }


def nonempty_candidate_fields(annotations, fields, required_empty):
    indexes = {field: fields.index(field) for field in required_empty if field in fields}
    found = {field: Counter() for field in indexes}
    for (_key, values), count in annotations.items():
        for field, index in indexes.items():
            if values[index]:
                found[field][values[index]] += count
    return {
        field: [{"value": value, "count": count} for value, count in sorted(values.items())]
        for field, values in found.items()
        if values
    }


def identity_differences(candidate, oracle, fields):
    indexes = tuple(fields.index(field) for field in ("Allele", "Feature_type", "Feature"))

    def indexed(rows):
        result = Counter()
        for (key, values), count in rows.items():
            result[(key, *(values[index] for index in indexes))] += count
        return result

    def rendered(identities, rows):
        return [
            {
                "variant": list(identity[0]),
                "allele": identity[1],
                "featureType": identity[2],
                "feature": identity[3],
                "rows": rows[identity],
            }
            for identity in sorted(identities)
        ]

    candidate_rows = indexed(candidate)
    oracle_rows = indexed(oracle)
    return {
        "missing": rendered(oracle_rows.keys() - candidate_rows.keys(), oracle_rows),
        "extra": rendered(candidate_rows.keys() - oracle_rows.keys(), candidate_rows),
    }


def compare(candidate, oracle, oracle_format="auto", contract=None, input_path=None):
    if oracle_format == "auto":
        oracle_format = "rest-json" if oracle.suffix.lower() == ".json" else "vcf"
    default_fields = REST_FIELDS if oracle_format == "rest-json" else FIELDS
    contract_data = load_contract(contract, default_fields)
    fields = contract_data["fields"]
    (
        candidate_variants,
        candidate_annotations,
        candidate_fields,
        candidate_records,
        candidate_samples,
    ) = parse_vcf(candidate, fields, contract_data["normalizers"])
    if oracle_format == "rest-json":
        oracle_variants, oracle_annotations = parse_rest(oracle, fields)
        oracle_fields = set(fields)
        oracle_records = None
        oracle_samples = None
    else:
        (
            oracle_variants,
            oracle_annotations,
            oracle_fields,
            oracle_records,
            oracle_samples,
        ) = parse_vcf(oracle, fields, contract_data["normalizers"])
    input_sha = sha256(input_path) if input_path is not None else None
    (
        candidate_annotations,
        applicable_field_differences,
        applied_field_differences,
        unused_field_differences,
    ) = apply_allowed_field_differences(
        candidate_annotations,
        oracle_annotations,
        fields,
        contract_data["allowedFieldDifferences"],
        input_sha,
    )
    candidate_annotations, applied_extras, unused_extras = apply_allowed_extras(
        candidate_annotations,
        oracle_annotations,
        fields,
        contract_data["allowedExtras"],
    )
    missing_variants = oracle_variants - candidate_variants
    extra_variants = candidate_variants - oracle_variants
    missing_annotations = oracle_annotations - candidate_annotations
    extra_annotations = candidate_annotations - oracle_annotations
    missing_candidate_fields = sorted(contract_data["candidateFields"] - candidate_fields)
    missing_oracle_fields = sorted(contract_data["oracleFields"] - oracle_fields)

    nonempty_fields = {}
    if contract_data["candidateEmptyFields"]:
        _, all_candidate_annotations, _, _, _ = parse_vcf(candidate, PRODUCTION_FIELDS)
        nonempty_fields = nonempty_candidate_fields(
            all_candidate_annotations,
            PRODUCTION_FIELDS,
            contract_data["candidateEmptyFields"],
        )

    record_report = None
    record_failures = []
    if oracle_records is not None:
        missing_records = oracle_records - candidate_records
        extra_records = candidate_records - oracle_records
        record_failures.extend((missing_records, extra_records))
        record_report = {
            "candidateVsOracle": {
                "missing": record_rows(missing_records),
                "extra": record_rows(extra_records),
            },
            "sampleNames": {
                "candidate": list(candidate_samples),
                "oracle": list(oracle_samples),
                "match": candidate_samples == oracle_samples,
            },
        }
        record_failures.append(candidate_samples != oracle_samples)
    if input_path is not None:
        input_records, input_samples = parse_input_vcf(input_path)
        candidate_missing_input = input_records - candidate_records
        candidate_extra_input = candidate_records - input_records
        record_failures.extend((candidate_missing_input, candidate_extra_input))
        if record_report is None:
            record_report = {}
        record_report["input"] = {"path": str(input_path), "sha256": sha256(input_path)}
        record_report["candidateVsInput"] = {
            "missing": record_rows(candidate_missing_input),
            "extra": record_rows(candidate_extra_input),
        }
        record_report["candidateSampleNamesMatchInput"] = candidate_samples == input_samples
        record_failures.append(candidate_samples != input_samples)
        if oracle_records is not None:
            oracle_missing_input = input_records - oracle_records
            oracle_extra_input = oracle_records - input_records
            record_failures.extend((oracle_missing_input, oracle_extra_input))
            record_report["oracleVsInput"] = {
                "missing": record_rows(oracle_missing_input),
                "extra": record_rows(oracle_extra_input),
            }
            record_report["oracleSampleNamesMatchInput"] = oracle_samples == input_samples
            record_failures.append(oracle_samples != input_samples)

    passed = not any(
        (
            missing_candidate_fields,
            missing_oracle_fields,
            missing_variants,
            extra_variants,
            missing_annotations,
            extra_annotations,
            unused_extras,
            unused_field_differences,
            nonempty_fields,
            *record_failures,
        )
    ) and bool(candidate_variants and candidate_annotations)
    report = {
        "schemaVersion": 4,
        "oracleFormat": oracle_format,
        "candidate": {"path": str(candidate), "sha256": sha256(candidate)},
        "oracle": {"path": str(oracle), "sha256": sha256(oracle)},
        "variantRecords": {
            "candidate": sum(candidate_variants.values()),
            "oracle": sum(oracle_variants.values()),
            "missing": sum(missing_variants.values()),
            "extra": sum(extra_variants.values()),
        },
        "annotationRows": {
            "candidate": sum(candidate_annotations.values()),
            "oracle": sum(oracle_annotations.values()),
            "missing": sum(missing_annotations.values()),
            "extra": sum(extra_annotations.values()),
        },
        "requiredFields": list(fields),
        "missingCandidateFields": missing_candidate_fields,
        "missingOracleFields": missing_oracle_fields,
        "nonemptyExcludedCandidateFields": nonempty_fields,
        "missingAnnotationExamples": examples(missing_annotations, fields),
        "extraAnnotationExamples": examples(extra_annotations, fields),
        "annotationDifferences": {
            "missing": annotation_rows(missing_annotations, fields),
            "extra": annotation_rows(extra_annotations, fields),
        },
        "identityComparison": field_mismatches(
            candidate_annotations, oracle_annotations, fields
        ),
        "identityDifferences": identity_differences(
            candidate_annotations, oracle_annotations, fields
        ),
        "passed": passed,
    }
    if record_report is not None:
        report["recordIntegrity"] = record_report
    if contract_data["report"] is not None:
        contract_report = dict(contract_data["report"])
        contract_report.update(
            appliedExtraRows=applied_extras,
            unusedAllowedExtraIdentities=unused_extras,
            applicableFieldDifferences=applicable_field_differences,
            appliedFieldDifferenceRows=applied_field_differences,
            unusedAllowedFieldDifferences=unused_field_differences,
        )
        report["contract"] = contract_report
    return report


def self_test():
    assert normalize_field("cds_start_NF&cds_end_NF", "unordered-flags") == normalize_field("cds_end_NF&cds_start_NF", "unordered-flags")
    assert normalize_field("cds_start_NF&cds_start_NF", "unordered-flags") != normalize_field("cds_start_NF", "unordered-flags")
    import tempfile

    assert position({}, "cdna") == ""
    assert position({"cdna_end": 7}, "cdna") == "?-7"
    assert position({"cdna_start": 7}, "cdna") == "7-?"
    assert position({"cdna_start": 9, "cdna_end": 7}, "cdna") == "7-9"

    header = (
        '##fileformat=VCFv4.2\n'
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: '
        + "|".join(ALL_FIELDS)
        + '">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    )
    values = dict.fromkeys(ALL_FIELDS, "")
    values.update(
        Allele="G",
        Consequence="missense_variant",
        Gene="ENSG1",
        Feature_type="Transcript",
        Feature="ENST1",
    )
    row = "|".join(values[field] for field in ALL_FIELDS)
    with tempfile.TemporaryDirectory() as directory:
        left = Path(directory) / "left.vcf"
        right = Path(directory) / "right.vcf"
        text = header + f"1\t10\t.\tA\tG\t.\tPASS\tCSQ={row}\n"
        left.write_text(text, encoding="utf-8")
        right.write_text(text, encoding="utf-8")
        assert compare(left, right)["passed"]

        compressed = Path(directory) / "right.vcf.gz"
        with gzip.open(compressed, mode="wt", encoding="utf-8") as handle:
            handle.write(text)
        assert compare(left, compressed)["passed"]

        input_vcf = Path(directory) / "input.vcf"
        input_vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t10\t.\tA\tG\t.\tPASS\t.\n",
            encoding="utf-8",
        )
        assert compare(left, right, input_path=input_vcf)["passed"]
        mutated_record = text.replace("1\t10\t.\tA\tG", "1\t10\trs-mutated\tA\tG")
        left.write_text(mutated_record, encoding="utf-8")
        failed = compare(left, right, input_path=input_vcf)
        assert not failed["passed"]
        assert failed["recordIntegrity"]["candidateVsInput"]["missing"]
        assert failed["recordIntegrity"]["candidateVsInput"]["extra"]
        left.write_text(text, encoding="utf-8")

        right.write_text(text.replace("missense_variant", "synonymous_variant"), encoding="utf-8")
        failed = compare(left, right)
        assert not failed["passed"]
        assert failed["annotationRows"]["missing"] == 1
        assert failed["annotationRows"]["extra"] == 1
        assert failed["identityDifferences"] == {"missing": [], "extra": []}
        assert len(failed["annotationDifferences"]["missing"]) == 1
        assert len(failed["annotationDifferences"]["extra"]) == 1

        right.write_text(text, encoding="utf-8")
        left.write_text(text.rstrip() + f",{row}\n", encoding="utf-8")
        duplicate = compare(left, right)
        assert not duplicate["passed"]
        assert duplicate["annotationRows"]["extra"] == 1
        assert duplicate["identityComparison"]["ambiguousMultirowIdentities"] == 1
        left.write_text(text, encoding="utf-8")

        rest = Path(directory) / "oracle.json"
        rest.write_text(
            json.dumps(
                [
                    {
                        "assembly_name": "GRCh38",
                        "input": "1 10 . A G . . .",
                        "transcript_consequences": [
                            {
                                "variant_allele": "G",
                                "consequence_terms": ["missense_variant"],
                                "gene_id": "ENSG1",
                                "transcript_id": "ENST1",
                            }
                        ],
                    }
                ]
            ),
            encoding="utf-8",
        )
        assert compare(left, rest)["passed"]

        offset_values = values.copy()
        offset_values["HGVS_OFFSET"] = "2"
        offset_row = "|".join(offset_values[field] for field in ALL_FIELDS)
        left.write_text(
            header + f"1\t10\t.\tA\tG\t.\tPASS\tCSQ={offset_row}\n",
            encoding="utf-8",
        )
        rest_document = json.loads(rest.read_text(encoding="utf-8"))
        rest_document[0]["transcript_consequences"][0]["hgvs_offset"] = 2
        rest.write_text(json.dumps(rest_document), encoding="utf-8")
        assert compare(left, rest)["passed"]
        rest_document[0]["transcript_consequences"][0]["hgvs_offset"] = 3
        rest.write_text(json.dumps(rest_document), encoding="utf-8")
        offset_failure = compare(left, rest)
        assert offset_failure["identityComparison"]["mismatchesByField"] == {
            "HGVS_OFFSET": 1
        }
        left.write_text(text, encoding="utf-8")
        del rest_document[0]["transcript_consequences"][0]["hgvs_offset"]
        rest.write_text(json.dumps(rest_document), encoding="utf-8")

        contract = Path(directory) / "contract.json"
        contract.write_text(
            json.dumps(
                {
                    "schemaVersion": 1,
                    "ignoredFields": [
                        {"field": "FLAGS", "reason": "fixture source difference"}
                    ],
                    "allowedExtraIdentities": [
                        {
                            "variant": ["1", "10", "A", "G"],
                            "allele": "G",
                            "featureType": "Transcript",
                            "feature": "ENST2",
                            "reason": "fixture source difference",
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )
        extra_values = values.copy()
        extra_values["Feature"] = "ENST2"
        extra_row = "|".join(extra_values[field] for field in ALL_FIELDS)
        left.write_text(text.rstrip() + f",{extra_row}\n", encoding="utf-8")
        uncontracted = compare(left, rest)
        assert uncontracted["identityDifferences"]["extra"][0]["feature"] == "ENST2"
        contracted = compare(left, rest, contract=contract)
        assert contracted["passed"]
        assert contracted["contract"]["appliedExtraRows"] == 1
        assert contracted["identityDifferences"] == {"missing": [], "extra": []}
        contract_document = json.loads(contract.read_text(encoding="utf-8"))
        contract_document["allowedExtraIdentities"][0]["feature"] = "ENST3"
        contract.write_text(json.dumps(contract_document), encoding="utf-8")
        assert not compare(left, rest, contract=contract)["passed"]

        exact = set(FIELDS)
        v2_contract = Path(directory) / "field-contract.json"
        v2_contract.write_text(
            json.dumps(
                {
                    "schemaVersion": 2,
                    "fields": [
                        {
                            "name": field,
                            "disposition": "exact" if field in exact else "excluded",
                            **(
                                {}
                                if field in exact
                                else {"reason": "not emitted by the source-matched oracle"}
                            ),
                            **(
                                {"requireCandidateEmpty": True}
                                if field == "Existing_variation"
                                else {}
                            ),
                        }
                        for field in PRODUCTION_FIELDS
                    ],
                    "allowedExtraIdentities": [],
                }
            ),
            encoding="utf-8",
        )
        left.write_text(text, encoding="utf-8")
        right.write_text(text, encoding="utf-8")
        assert compare(left, right, contract=v2_contract)["passed"]

        field_difference_ledger = Path(directory) / "field-differences.json"
        field_difference_ledger.write_text(
            json.dumps(
                {
                    "schemaVersion": 1,
                    "inputSha256": sha256(input_vcf),
                    "differences": [
                        {
                            "variant": ["1", "10", "A", "G"],
                            "allele": "G",
                            "featureType": "Transcript",
                            "feature": "ENST1",
                            "field": "HGVSp",
                            "candidate": "ENSP1:p.Arg1ProfsTer2",
                            "oracle": "ENSP1:p.Arg1ProfsTer?",
                            "reason": "fixture for one reviewed VEP behavior",
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )
        field_contract = json.loads(v2_contract.read_text(encoding="utf-8"))
        field_contract["allowedFieldDifferences"] = field_difference_ledger.name
        v2_contract.write_text(json.dumps(field_contract), encoding="utf-8")
        candidate_values = values.copy()
        candidate_values["HGVSp"] = "ENSP1:p.Arg1ProfsTer2"
        oracle_values = values.copy()
        oracle_values["HGVSp"] = "ENSP1:p.Arg1ProfsTer?"
        left.write_text(
            header
            + "1\t10\t.\tA\tG\t.\tPASS\tCSQ="
            + "|".join(candidate_values[field] for field in ALL_FIELDS)
            + "\n",
            encoding="utf-8",
        )
        right.write_text(
            header
            + "1\t10\t.\tA\tG\t.\tPASS\tCSQ="
            + "|".join(oracle_values[field] for field in ALL_FIELDS)
            + "\n",
            encoding="utf-8",
        )
        allowed = compare(left, right, contract=v2_contract, input_path=input_vcf)
        assert allowed["passed"]
        assert allowed["contract"]["applicableFieldDifferences"] == 1
        assert allowed["contract"]["appliedFieldDifferenceRows"] == 1
        left.write_text(left.read_text(encoding="utf-8").replace("fsTer2", "fsTer3"))
        stale = compare(left, right, contract=v2_contract, input_path=input_vcf)
        assert not stale["passed"]
        assert stale["contract"]["unusedAllowedFieldDifferences"]
        field_contract.pop("allowedFieldDifferences")
        v2_contract.write_text(json.dumps(field_contract), encoding="utf-8")

        normalized_contract = json.loads(v2_contract.read_text(encoding="utf-8"))
        normalizers = {
            "HGVSp": "hgvsp-gff-protein-version",
            "REF_ALLELE": "empty-allele-dash",
            "UPLOADED_ALLELE": "uploaded-allele-delimiter",
        }
        for field in normalized_contract["fields"]:
            if field["name"] in normalizers:
                field["normalizer"] = normalizers[field["name"]]
                field["reason"] = "equivalent source-specific serialization"
        v2_contract.write_text(json.dumps(normalized_contract), encoding="utf-8")
        candidate_values = values.copy()
        candidate_values.update(
            HGVSp="ENSP00000001.8:p.Arg1Gly",
            REF_ALLELE="-",
            UPLOADED_ALLELE="A/G&T",
        )
        oracle_values = values.copy()
        oracle_values.update(
            HGVSp="ENSP00000001.1:p.Arg1Gly",
            REF_ALLELE="",
            UPLOADED_ALLELE="A/G/T",
        )
        left.write_text(
            header
            + "1\t10\t.\tA\tG\t.\tPASS\tCSQ="
            + "|".join(candidate_values[field] for field in ALL_FIELDS)
            + "\n",
            encoding="utf-8",
        )
        right.write_text(
            header
            + "1\t10\t.\tA\tG\t.\tPASS\tCSQ="
            + "|".join(oracle_values[field] for field in ALL_FIELDS)
            + "\n",
            encoding="utf-8",
        )
        assert compare(left, right, contract=v2_contract)["passed"]
        right.write_text(
            right.read_text(encoding="utf-8").replace("p.Arg1Gly", "p.Arg1Val"),
            encoding="utf-8",
        )
        assert not compare(left, right, contract=v2_contract)["passed"]

        left.write_text(text, encoding="utf-8")
        right.write_text(text, encoding="utf-8")
        populated = values.copy()
        populated["Existing_variation"] = "rs1"
        populated_row = "|".join(populated[field] for field in ALL_FIELDS)
        left.write_text(
            header + f"1\t10\t.\tA\tG\t.\tPASS\tCSQ={populated_row}\n",
            encoding="utf-8",
        )
        excluded = compare(left, right, contract=v2_contract)
        assert not excluded["passed"]
        assert "Existing_variation" in excluded["nonemptyExcludedCandidateFields"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate", nargs="?", type=Path)
    parser.add_argument("oracle", nargs="?", type=Path)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--contract", type=Path)
    parser.add_argument("--input", type=Path, help="original VCF for record-integrity checks")
    parser.add_argument(
        "--oracle-format", choices=("auto", "vcf", "rest-json"), default="auto"
    )
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        print("VEP concordance comparator self-test passed")
        return
    if not args.candidate or not args.oracle:
        parser.error("candidate and oracle VCF files are required")
    report = compare(
        args.candidate,
        args.oracle,
        args.oracle_format,
        args.contract,
        args.input,
    )
    rendered = json.dumps(report, indent=2, sort_keys=True)
    if args.json:
        args.json.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
