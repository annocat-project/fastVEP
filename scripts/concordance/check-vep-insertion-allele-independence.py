"""Check an insertion matrix against its frozen, single-ALT VEP controls.

Development evidence only: preserves raw multi-ALT oracle disagreements and does
not edit a concordance contract or certify biological correctness of the oracle.
Every input must be a pure anchored insertion, with a single-ALT control for each
ALT at the same locus. CSQ row multiplicity and transcript membership are checked.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path


def records(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('##INFO=<ID=CSQ,'):
                fields = line.split('Format: ', 1)[1].split('"', 1)[0].split('|')
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            alts = cols[4].split(',')
            assert all(a.startswith(cols[3]) and len(a) > len(cols[3]) for a in alts), cols[:5]
            csq = next(x[4:] for x in cols[7].split(';') if x.startswith('CSQ='))
            rows = [dict(zip(fields, x.split('|'))) for x in csq.split(',')]
            yield cols, {a[len(cols[3]):] for a in alts}, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('candidate', type=Path)
    parser.add_argument('oracle', type=Path)
    parser.add_argument('contract', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    contract = json.loads(args.contract.read_text())
    fields = [f['name'] for f in contract['fields'] if f['disposition'] == 'exact']
    # FLAGS list ordering follows the existing matrix contract, not a new waiver.
    flags_policy = next(f for f in contract['fields'] if f['name'] == 'FLAGS')
    assert flags_policy['disposition'] in ('exact', 'normalized'), flags_policy
    if 'FLAGS' not in fields:
        fields.append('FLAGS')

    def value(row, field):
        v = row.get(field, '')
        if field == 'FLAGS' and flags_policy.get('normalizer') == 'unordered-flags':
            return sorted(v.split('&')) if v else []
        return v

    def identity(cols, row):
        return (*cols[:2], cols[3], row['Allele'], row['Feature_type'], row['Feature'])

    def singles(path):
        controls = {}
        count = 0
        for cols, alts, rows in records(path):
            if len(alts) != 1:
                continue
            count += 1
            for row in rows:
                key = identity(cols, row)
                assert key not in controls, ('duplicate control identity', key)
                controls[key] = row
        assert count > 0
        return controls, count

    gold, control_count = singles(args.oracle)
    native, native_count = singles(args.candidate)
    assert gold.keys() == native.keys(), 'single-ALT transcript membership differs'
    assert control_count == native_count
    by_allele = {}
    for key in gold:
        by_allele.setdefault(key[:4], set()).add(key)
    summary = dict(qualificationEligible=False, controlRecords=control_count,
                   controlIdentities=len(gold), fields=fields, records=0, identities=0)
    counts = {k: Counter() for k in ('candidateVsSingleOracle', 'candidateOrderDependence',
                                    'rawOracleDisagreements', 'disagreementsMatchingSingleOracle',
                                    'disagreementsNotMatchingSingleOracle')}
    # Only the complete non-crashing multi-ALT oracle subset has output to compare.
    raw_oracle = {}
    for cols, alts, rows in records(args.oracle):
        for row in rows:
            key = (cols[2], *identity(cols, row))
            assert key not in raw_oracle
            raw_oracle[key] = row
    with (args.output / 'differences.jsonl').open('x') as ledger:
        for cols, alts, rows in records(args.candidate):
            summary['records'] += 1
            expected = set()
            for allele in alts:
                expected.update(by_allele[(*cols[:2], cols[3], allele)])
            observed = Counter(identity(cols, row) for row in rows)
            assert set(observed) == expected and all(n == 1 for n in observed.values()), ('membership', cols[:5])
            for row in rows:
                summary['identities'] += 1
                key = identity(cols, row)
                raw = raw_oracle.get((cols[2], *key))
                for field in fields:
                    a, b, n = value(row, field), value(gold[key], field), value(native[key], field)
                    categories = []
                    if a != b:
                        categories.append('candidateVsSingleOracle')
                    if a != n:
                        categories.append('candidateOrderDependence')
                    if raw is not None and a != value(raw, field):
                        categories.extend(['rawOracleDisagreements',
                            'disagreementsMatchingSingleOracle' if a == b else 'disagreementsNotMatchingSingleOracle'])
                    if categories:
                        for category in categories:
                            counts[category][field] += 1
                        ledger.write(json.dumps(dict(record=cols[2], identity=key, field=field,
                            candidate=a, singleOracle=b, singleCandidate=n,
                            rawOracle=value(raw, field) if raw else None, categories=categories)) + '\n')
    summary.update({k: dict(v) for k, v in counts.items()})
    summary['alleleIndependencePassed'] = not counts['candidateOrderDependence']
    summary['singleOracleComparisonPassed'] = not counts['candidateVsSingleOracle']
    summary['pins'] = {}
    for path in (args.candidate, args.oracle, args.contract, Path(__file__)):
        with path.open('rb') as stream:
            summary['pins'][str(path)] = hashlib.file_digest(stream, 'sha256').hexdigest()
    (args.output / 'report.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    return 0 if summary['alleleIndependencePassed'] and summary['singleOracleComparisonPassed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
