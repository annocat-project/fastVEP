"""Verify the full input package, optionally including downloaded VEP oracles."""
import argparse
import json
from pathlib import Path
from run import digest, member


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', required=True, type=Path)
    parser.add_argument('--oracles', action='store_true')
    args = parser.parse_args()
    manifest = json.loads((args.data / 'manifest.json').read_text(encoding='utf-8'))
    selected = [a for a in manifest['assets'] if a['role'] == 'input' or args.oracles]
    for asset in selected:
        if digest(member(args.data, asset['path'])) != asset['sha256']:
            raise ValueError('Asset checksum mismatch: ' + asset['path'])
    for corpus in manifest['corpora']:
        if 'contract' in corpus and digest(member(args.data, corpus['contract'])) != corpus['contractSha256']:
            raise ValueError('Contract checksum mismatch: ' + corpus['name'])
    print(json.dumps({'verified': True, 'assets': len(selected), 'oraclesIncluded': args.oracles}))


if __name__ == '__main__':
    main()
