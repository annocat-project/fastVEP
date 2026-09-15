"""Run the recorded local VEP profile on an explicit input; retain failures as failures."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True, type=Path,
                        help='Three-bank manifest containing pinned oracle fields and arguments')
    parser.add_argument('--vep', required=True, type=Path)
    parser.add_argument('--perl', default='perl')
    parser.add_argument('--cache-dir', required=True, type=Path)
    parser.add_argument('--fasta', required=True, type=Path)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path, help='New output directory')
    args = parser.parse_args()
    oracle = json.loads(args.manifest.read_text(encoding='utf-8'))['oracle']
    fields = oracle['fields']
    if not isinstance(fields, list) or not all(isinstance(f, str) for f in fields):
        parser.error('Oracle manifest must contain a field-name list')
    args.output.mkdir(parents=True, exist_ok=False)
    command = [args.perl, str(args.vep.resolve()), *oracle['args'], '--fields', ','.join(fields),
               '--dir_cache', str(args.cache_dir.resolve()), '--fasta', str(args.fasta.resolve()),
               '--input_file', str(args.input.resolve()), '--output_file', str((args.output / 'oracle.vcf').resolve())]
    report = {'complete': False, 'qualificationEligible': False, 'inputSha256': digest(args.input),
              'vepScriptSha256': digest(args.vep), 'fastaSha256': digest(args.fasta),
              'profileManifestSha256': digest(args.manifest),
              'runtimeIdentityVerified': False,
              'note': 'The caller must verify the pinned VEP runtime and cache before treating this output as an oracle.'}
    with (args.output / 'oracle.log').open('w', encoding='utf-8') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    report['exitCode'] = result.returncode
    output = args.output / 'oracle.vcf'
    report['complete'] = result.returncode == 0 and output.is_file()
    if output.is_file():report['outputSha256'] = digest(output)
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'complete': report['complete'], 'exitCode': result.returncode}))
    raise SystemExit(0 if report['complete'] else 1)


if __name__ == '__main__':
    main()
