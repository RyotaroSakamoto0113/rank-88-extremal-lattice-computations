#!/usr/bin/env python3
"""Locate the last affine target by its exact Gram, then verify its O-basis.

Defaults are relative to this script, so a fresh proof run may copy this script
and lowspan_seed.json into output/rank4_audit without any old candidate list.
No modules containing earlier proof computations are imported.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time


def add(x, y):
    return x[0]+y[0], x[1]+y[1]


def mul(x, y):
    return x[0]*y[0]-6*x[1]*y[1], x[0]*y[1]+x[1]*y[0]+x[1]*y[1]


def conj(x):
    return x[0]+x[1], -x[1]


def determinant(H):
    n = len(H)
    result = (0, 0)
    for p in itertools.permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i+1, n))
        value = (-1 if inversions % 2 else 1, 0)
        for i in range(n):
            value = mul(value, H[i][p[i]])
        result = add(result, value)
    return result


def inner(H, u, v):
    n = len(H)
    result = (0, 0)
    for i in range(n):
        for j in range(n):
            result = add(result, mul(mul((u[i], u[n+i]), H[i][j]), conj((v[j], v[n+j]))))
    return result


def checked(condition, message):
    # These checks also execute under python -O.
    if not condition:
        raise ValueError(message)


def main():
    start = time.perf_counter()
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--targets', type=Path, default=here/'independent_h3_targets.jsonl')
    parser.add_argument('--seed', type=Path, default=here/'lowspan_seed.json')
    parser.add_argument('--output', type=Path, default=here/'h3_lowspan_rebased.json')
    args = parser.parse_args()
    seed = json.loads(args.seed.read_text())
    H, C = seed['source_H'], seed['basis_rows']
    checked(len(H) == 4 and all(len(row) == 4 for row in H), 'Source Gram is not 4 by 4')
    checked(len(C) == 4 and all(len(row) == 8 for row in C), 'Basis must have four rows of eight integer coefficients')
    checked(all(type(x) is int for row in C for x in row), 'Basis coefficients must be integers')
    checked(all(len(z) == 2 and all(type(x) is int for x in z) for row in H for z in row), 'Gram coefficients must be pairs of integers')
    checked(all(tuple(H[j][i]) == conj(H[i][j]) for i in range(4) for j in range(4)), 'Source Gram is not Hermitian')
    targets_bytes = args.targets.read_bytes()
    matches = [t for line in targets_bytes.splitlines() if (t := json.loads(line))['H'] == H]
    checked(len(matches) == 1, f'Expected exactly one target with the seed Gram; found {len(matches)}')
    source = matches[0]
    detH = determinant(H)
    checked(detH[1] == 0 and detH[0] > 0, 'Source determinant is not positive rational')
    checked(source.get('det') == detH[0], 'Target determinant metadata differs from direct computation')
    for k in range(1, 5):
        minor = determinant([row[:k] for row in H[:k]])
        checked(minor[1] == 0 and minor[0] > 0, 'Source Gram is not positive definite')
    fieldC = [[(row[j], row[j+4]) for j in range(4)] for row in C]
    detC = determinant(fieldC)
    checked(detC == (1, 0), f'Basis determinant is not 1: {detC}')
    checked(list(detC) == seed['basis_determinant'], 'Seed basis determinant metadata differs')
    HC = [[list(inner(H, u, v)) for v in C] for u in C]
    checked(HC == seed['H'], 'Transformed Gram differs from independently supplied expected Gram')
    norms = [HC[i][i] for i in range(4)]
    checked(norms == [[12, 0], [14, 0], [13, 0], [13, 0]], 'Transformed norms differ')
    checked(HC[0][1] == [-3, -4], 'Initial affine binary inner product differs')
    checked(determinant([row[:2] for row in HC[:2]]) == (51, 0), 'Initial binary determinant differs')
    checked(determinant(HC) == detH, 'Unimodular change altered the Gram determinant')
    result = dict(source_id=source['id'], source_H=H, basis_rows=C,
                  basis_determinant=detC, H=HC, source_determinant=detH[0],
                  target_sha256=hashlib.sha256(targets_bytes).hexdigest(),
                  seed_sha256=hashlib.sha256(args.seed.read_bytes()).hexdigest(),
                  seconds=time.perf_counter()-start)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(source_id=source['id'], basis_determinant=detC,
                         norms=[12, 14, 13, 13], source_determinant=detH[0],
                         seconds=result['seconds'], complete=True)))


if __name__ == '__main__':
    main()
