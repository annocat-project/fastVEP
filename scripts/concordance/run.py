"""Run a frozen, sample-free VEP regression bank with an explicit fastVEP build."""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def member(root, name):
    path = (root / name).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError('Manifest member is outside the suite directory')
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--verify-only', action='store_true')
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--cache', type=Path)
    parser.add_argument('--fasta', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--corpus', action='append', help='Run only named banks; repeat as needed')
    args = parser.parse_args()
    manifest = json.loads((args.data / 'manifest.json').read_text(encoding='utf-8'))
    if manifest['schemaVersion'] != 1:
        parser.error('Unsupported suite manifest')
    names = {c['name'] for c in manifest['corpora']}
    if args.corpus and not set(args.corpus) <= names:
        parser.error('Unknown corpus')
    selected = [c for c in manifest['corpora'] if not args.corpus or c['name'] in args.corpus]
    for corpus in selected:
        for key, expected in [('input', corpus['inputProvenance']['publishedSha256']),
                              ('oracle', corpus['oracleProvenance']['publishedSha256']),
                              ('contract', corpus['contractSha256'])]:
            if digest(member(args.data, corpus[key])) != expected:
                raise ValueError(f"Checksum mismatch: {corpus['name']} {key}")
    if args.verify_only:
        print(json.dumps({'verified': True, 'corpora': len(selected)}))
        return
    if any(getattr(args, name) is None for name in ['binary', 'cache', 'fasta', 'output']):
        parser.error('Annotation requires --binary, --cache, --fasta and --output')
    args.output.mkdir(parents=True, exist_ok=False)
    report = {'suite': manifest['suite'], 'binarySha256': digest(args.binary),
              'cacheSha256': digest(args.cache), 'fastaSha256': digest(args.fasta),
              'matchesRecordedCache': digest(args.cache) == manifest['testedCacheSha256'],
              'qualificationEligible': False, 'corpora': [], 'passed': False}
    compare = Path(__file__).with_name('compare-vep-in-partitions.py')
    for corpus in selected:
        destination = args.output / corpus['name']
        destination.mkdir()
        candidate = destination / 'candidate.vcf'
        with (destination / 'run.log').open('w', encoding='utf-8') as log:
            subprocess.run([str(args.binary.resolve()), 'annotate', '--input',
                            str(member(args.data, corpus['input'])), '--output', str(candidate.resolve()),
                            '--transcript-cache', str(args.cache.resolve()), '--fasta', str(args.fasta.resolve()),
                            *manifest['annotationArgs']], stdout=log, stderr=subprocess.STDOUT, check=True)
            result = subprocess.run([sys.executable, '-X', 'utf8', str(compare), str(candidate),
                                     str(member(args.data, corpus['oracle'])), str(member(args.data, corpus['input'])),
                                     str(member(args.data, corpus['contract'])), str(destination / 'comparison')],
                                    stdout=log, stderr=subprocess.STDOUT)
        result_path = destination / 'comparison/report.json'
        passed = result.returncode == 0 and result_path.exists() and json.loads(result_path.read_text())['passed']
        report['corpora'].append({'name': corpus['name'], 'passed': passed})
        (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
        print(json.dumps(report['corpora'][-1]), flush=True)
        if not passed:
            raise SystemExit(1)
    report['passed'] = True
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
