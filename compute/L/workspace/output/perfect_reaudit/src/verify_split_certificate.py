#!/usr/bin/env python3
"""Independently verify a C115 split perfection certificate using Python stdlib.

The certificate proves the supplied norm-8 orbit span has dimension 3916.
Its implication for perfection uses the separately proved minimum min(L)=8.
No search, eigenbasis-generation, or production rank code is imported.
"""
import argparse
import datetime
import hashlib
import json
import math
import resource
import time
from pathlib import Path


def need(condition, message):
    if not condition:
        raise ValueError(message)


def integer(value):
    return type(value) is int


def matrix(value, n, m, name):
    need(isinstance(value, list) and len(value) == n, name + ': row count')
    need(all(isinstance(row, list) and len(row) == m and all(integer(x) for x in row)
             for row in value), name + ': integer dimensions')
    return value


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(x) for x in zip(*a)]


def multiply(a, b, modulus=None):
    """Sparse-aware exact multiplication; reduction only after each output row."""
    n, m = len(a), len(b[0])
    result = []
    for i in range(n):
        row = [0] * m
        for k, entry in enumerate(a[i]):
            if entry:
                bk = b[k]
                for j in range(m):
                    row[j] += entry * bk[j]
        if modulus is not None:
            row = [x % modulus for x in row]
        result.append(row)
    return result


def power(a, exponent):
    r = identity(len(a))
    b = a
    while exponent:
        if exponent & 1:
            r = multiply(r, b)
        exponent //= 2
        if exponent:
            b = multiply(b, b)
    return r


def determinant_mod(a, p):
    """Direct Gaussian determinant, no production elimination dependencies."""
    n = len(a)
    b = [[x % p for x in row] for row in a]
    det = 1
    for k in range(n):
        pivot = next((i for i in range(k, n) if b[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            b[k], b[pivot] = b[pivot], b[k]
            det = -det
        v = b[k][k]
        det = det * v % p
        inv = pow(v, -1, p)
        for i in range(k + 1, n):
            factor = b[i][k] * inv % p
            if factor:
                for j in range(k + 1, n):
                    b[i][j] = (b[i][j] - factor * b[k][j]) % p
            b[i][k] = 0
    return det % p


def positive_bareiss_det(a):
    """Positive leading principal minors verify positive definiteness."""
    b = [row[:] for row in a]
    n = len(b)
    previous = 1
    leading = []
    for k in range(n):
        pivot = b[k][k]
        need(pivot > 0, 'Gram matrix has nonpositive leading principal minor')
        leading.append(pivot)
        if k + 1 == n:
            break
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = pivot * b[i][j] - b[i][k] * b[k][j]
                quotient, remainder = divmod(numerator, previous)
                need(remainder == 0, 'nonexact Bareiss division')
                b[i][j] = quotient
        previous = pivot
        for i in range(k + 1, n):
            b[i][k] = 0
    return leading[-1], leading


def prime(p):
    return p >= 2 and all(p % d for d in range(2, math.isqrt(p) + 1))


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for piece in iter(lambda: f.read(1048576), b''):
            h.update(piece)
    return h.hexdigest()


def snapshot():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return time.perf_counter(), r.ru_utime, r.ru_stime


def elapsed(start):
    stop = snapshot()
    return dict(zip(('wall_seconds', 'user_cpu_seconds', 'system_cpu_seconds'),
                    (stop[i] - start[i] for i in range(3))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--certificate', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path,
                        help='New output directory; existing directories are refused')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = snapshot()
    summary = {
        'status': 'running',
        'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'input': str(args.input.resolve()),
        'certificate': str(args.certificate.resolve()),
        'stages': [],
        'external_prerequisite': 'The minimum of this same lattice is 8.',
    }

    def save():
        (out / 'verification.json').write_text(json.dumps(summary, indent=2) + '\n')

    def stage(label, function):
        t = snapshot()
        result = function()
        summary['stages'].append({'name': label, **elapsed(t)})
        save()
        return result

    save()
    try:
        paths = [args.input.resolve(), args.certificate.resolve(), Path(__file__).resolve()]
        summary['input_hashes'] = {str(path): sha(path) for path in paths}

        def read_inputs():
            ambient = json.loads(args.input.read_text())
            certificate = json.loads(args.certificate.read_text())
            return ambient, certificate

        ambient, cert = stage('read_and_parse_inputs', read_inputs)
        n, order = 88, 115
        gmat = matrix(ambient['G'], n, n, 'G')
        m5 = matrix(ambient['M5'], n, n, 'M5')
        m23 = matrix(ambient['M23'], n, n, 'M23')
        p, zeta = cert['p'], cert['zeta']

        def field_checks():
            need(integer(p) and 2 < p <= 2**31 and prime(p), 'small odd prime required')
            need(p % order == 1, 'prime must split C115')
            need(integer(zeta), 'integer root required')
            need(pow(zeta, order, p) == 1 and pow(zeta, 5, p) != 1
                 and pow(zeta, 23, p) != 1, 'root has exact order 115')
            need(cert['characters'] and all(integer(x) for x in cert['characters']),
                 'integer character labels required')
            need(sorted(cert['characters']) == [a for a in range(order)
                                               if math.gcd(a, order) == 1],
                 'all 88 primitive characters required exactly once')

        stage('prime_root_and_characters', field_checks)

        def ambient_checks():
            need(gmat == transpose(gmat), 'Gram symmetry')
            need(all(gmat[i][i] % 2 == 0 for i in range(n)), 'Gram evenness')
            determinant, leading = positive_bareiss_det(gmat)
            need(determinant == 1, 'Gram determinant must be 1')
            (out / 'gram_leading_minors.json').write_text(json.dumps(leading) + '\n')
            one = identity(n)
            need(power(m5, 5) == one and m5 != one, 'M5 exact order 5')
            need(power(m23, 23) == one and m23 != one, 'M23 exact order 23')
            product = multiply(m5, m23)
            need(product == multiply(m23, m5), 'commuting integral actions')
            need(multiply(transpose(m5), multiply(gmat, m5)) == gmat, 'M5 isometry')
            need(multiply(transpose(m23), multiply(gmat, m23)) == gmat, 'M23 isometry')
            return product

        g = stage('exact_gram_orders_commutation_isometries', ambient_checks)
        transform = matrix(cert['T'], n, n, 'T')
        transform = [[x % p for x in row] for row in transform]

        def eigenbasis_checks():
            determinant = determinant_mod(transform, p)
            need(determinant != 0, 'left eigenbasis invertibility')
            lhs = multiply(transform, g, p)
            rhs = [[pow(zeta, cert['characters'][i], p) * entry % p
                    for entry in transform[i]] for i in range(n)]
            need(lhs == rhs, 'left eigenbasis Tg=DT')
            return determinant

        summary['T_determinant_mod_p'] = stage('left_eigenbasis', eigenbasis_checks)
        seeds = cert['seeds']
        need(isinstance(seeds, list) and len(seeds) >= 44, 'at least 44 seeds required')
        matrix(seeds, len(seeds), n, 'seeds')

        def seed_checks():
            for number, x in enumerate(seeds):
                norm = sum(gmat[i][i] * x[i] * x[i] for i in range(n))
                norm += 2 * sum(gmat[i][j] * x[i] * x[j]
                                for i in range(n) for j in range(i + 1, n))
                need(norm == 8, 'seed {} norm {} instead of 8'.format(number, norm))
            return [[sum(row[i] * x[i] for i in range(n)) % p
                     for row in transform] for x in seeds]

        images = stage('all_seed_norms_and_transformed_coordinates', seed_checks)
        expected = [[] for _ in range(order)]
        for i, a in enumerate(cert['characters']):
            for j in range(i, n):
                expected[(a + cert['characters'][j]) % order].append([i, j])
        need(sum(map(len, expected)) == 3916, 'symmetric-square dimension')
        for k in range(order):
            dimension = 44 if k == 0 else 42 if k % 5 == 0 else 33 if k % 23 == 0 else 32
            need(len(expected[k]) == dimension, 'character dimension formula')

        def block_checks():
            blocks = cert['blocks']
            need(isinstance(blocks, list) and len(blocks) == order, 'all 115 blocks')
            need(all(integer(block['character']) for block in blocks), 'integer block labels')
            need(sorted(block['character'] for block in blocks) == list(range(order)),
                 'each block character appears exactly once')
            results = []
            for block in blocks:
                t = snapshot()
                k = block['character']
                pairs = matrix(block['pairs'], len(expected[k]), 2, 'block pairs')
                need(sorted(pairs) == expected[k], 'exact pair partition for block {}'.format(k))
                ids = block['seed_ids']
                need(isinstance(ids, list) and len(ids) == len(pairs)
                     and all(integer(i) and 0 <= i < len(seeds) for i in ids)
                     and len(set(ids)) == len(ids), 'square minor seed IDs')
                minor = [[images[s][i] * images[s][j] % p for i, j in pairs] for s in ids]
                determinant = determinant_mod(minor, p)
                need(determinant != 0, 'singular block minor k={}'.format(k))
                results.append({'character': k, 'dimension': len(pairs),
                                'minor_determinant_mod_p': determinant,
                                'seed_ids': ids, **elapsed(t)})
            results.sort(key=lambda x: x['character'])
            (out / 'verified_blocks.json').write_text(json.dumps(results, indent=2) + '\n')
            return results

        results = stage('115_independent_square_minor_determinants', block_checks)
        need(sum(result['dimension'] for result in results) == 3916, 'full modular span')
        for path in paths:
            need(sha(path) == summary['input_hashes'][str(path)], 'input changed during verification')
        summary.update(status='verified', p=p, zeta=zeta % p, seed_count=len(seeds),
                       block_count=115, modular_span_dimension=3916,
                       rational_span_dimension=3916,
                       positive_even_unimodular_gram_verified=True,
                       integral_isometries_and_orders_verified=True,
                       all_seed_norms_8=True,
                       conclusion='The norm-8 g-orbit outer products span Sym_88(Q).',
                       proves_perfection_given_external_minimum_8=True)
    except BaseException as error:
        summary.update(status='failed', error=type(error).__name__ + ': ' + str(error))
        raise
    finally:
        summary['total_timing'] = elapsed(started)
        summary['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        save()
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
