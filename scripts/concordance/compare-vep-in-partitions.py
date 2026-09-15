"""Run the unchanged raw VEP comparator in bounded, disjoint input partitions."""
import argparse
from collections import Counter
import gc
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def variant_key(line):
    fields = line.split('\t', 5)
    return tuple(fields[i] for i in (0, 1, 3, 4))


def run(candidate, oracle, original, contract, output, partition_size=500):
    comparator = Path(__file__).with_name('compare-vep-concordance.py')
    spec = importlib.util.spec_from_file_location('vep_comparator', comparator)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rules = module.load_contract(contract, module.FIELDS)
    assert not rules['allowedExtras'] and not rules['allowedFieldDifferences'], 'Only raw contracts may be partitioned'
    assert partition_size > 0
    pins = {str(p): sha(p) for p in (candidate, oracle, original, contract, comparator, Path(__file__))}
    output.mkdir(exist_ok=False)
    assignment = {}
    with original.open() as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            key = variant_key(line)
            assert key not in assignment, 'Input variant keys must be unique; duplicates require the full comparator'
            assignment[key] = len(assignment) // partition_size
    assert assignment
    count = max(assignment.values()) + 1
    folders = [output / f'part-{i:04d}' for i in range(count)]
    for folder in folders:
        folder.mkdir()
    for source, name in ((original, 'input.vcf'), (candidate, 'candidate.vcf'), (oracle, 'oracle.vcf')):
        opener = gzip.open if source.suffix == '.gz' else open
        writers = [(folder / name).open('w', encoding='utf-8', newline='') for folder in folders]
        try:
            with opener(source, 'rt', encoding='utf-8', newline='') as stream:
                for line in stream:
                    if line.startswith('#'):
                        for writer in writers:
                            writer.write(line)
                    else:
                        key = variant_key(line)
                        if key not in assignment:
                            error = dict(passed=False, source=str(source), unexpectedVariant=key)
                            (output / 'partition-error.json').write_text(json.dumps(error, indent=2))
                            raise AssertionError(error)
                        writers[assignment[key]].write(line)
        finally:
            for writer in writers:
                writer.close()

    totals = {name: Counter() for name in ('variantRecords', 'annotationRows', 'identityComparison')}
    fields = Counter()
    parts = []
    for folder in folders:
        report = module.compare(folder / 'candidate.vcf', folder / 'oracle.vcf',
                                contract=contract, input_path=folder / 'input.vcf')
        assert report['requiredFields'] == list(rules['fields'])
        path = folder / 'comparison.json'
        with path.open('w') as stream:
            json.dump(report, stream, indent=2)
            stream.write('\n')
        for name, total in totals.items():
            total.update({key: value for key, value in report[name].items() if isinstance(value, int)})
        fields.update(report['identityComparison']['mismatchesByField'])
        parts.append(dict(path=str(path), sha256=sha(path), passed=report['passed']))
        del report
        gc.collect()
    for path, expected in pins.items():
        assert sha(Path(path)) == expected, path
    result = dict(complete=True, passed=all(p['passed'] for p in parts), qualificationEligible=False,
                  scope='Raw comparison of complete pre-existing outputs. Every input variant key belongs to exactly one partition; all its records and annotations stay together. No source or contract changes.',
                  inputRecords=len(assignment), partitionSize=partition_size,
                  parts=parts, pins=pins, **{name: dict(total) for name, total in totals.items()})
    result['identityComparison']['mismatchesByField'] = dict(fields)
    with (output / 'report.json').open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in ('passed', 'inputRecords', 'identityComparison')}), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('candidate', 'oracle', 'input', 'contract', 'output'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--partition-size', type=int, default=500)
    args = parser.parse_args()
    run(args.candidate, args.oracle, args.input, args.contract, args.output, args.partition_size)
