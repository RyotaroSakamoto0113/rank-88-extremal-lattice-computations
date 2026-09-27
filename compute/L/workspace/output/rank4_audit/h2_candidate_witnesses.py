"""Find and independently check short binary obstructions in every h=2 target."""
import collections
import hashlib
import json
import pathlib
import subprocess
import time

ROOT = pathlib.Path(__file__).resolve().parent
TYPES = [(6, 9, (2, 2)), (6, 9, (4, -2)), (8, 8, (3, 2)), (8, 8, (5, -2))]


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def mul(x, y):
    a, b = x
    c, d = y
    return a*c-6*b*d, a*d+b*c+b*d


def conj(x):
    return x[0]+x[1], -x[1]


def inner(H, u, v):
    n = len(H)
    value = (0, 0)
    for i in range(n):
        for j in range(n):
            value = add(value, mul(mul((u[i], u[n+i]), H[i][j]), conj((v[j], v[n+j]))))
    return value


def gram_input(H):
    n = len(H)
    G = [[0]*(2*n) for _ in range(2*n)]
    W = [[0]*(2*n) for _ in range(2*n)]
    for i in range(n):
        for j in range(n):
            a, b = H[i][j]
            assert tuple(H[j][i]) == conj((a, b))
            G[i][j] = 2*a+b
            G[i][j+n] = a+12*b
            G[i+n][j] = a-11*b
            G[i+n][j+n] = 6*(2*a+b)
        W[i][i+n] = -6
        W[i+n][i] = W[i+n][i+n] = 1
    assert all(G[i][j] == G[j][i] for i in range(2*n) for j in range(2*n))
    return f'{2*n} 18\n' + ''.join(' '.join(map(str, row))+'\n' for row in G+W)


def main():
    started = time.perf_counter()
    source = ROOT/'independent_h2_targets.jsonl'
    candidates = [json.loads(s) for s in source.read_text().splitlines()]
    inputs = ROOT/'h2_candidate_inputs'
    inputs.mkdir(exist_ok=True)
    records = []
    for target in candidates:
        idx, H = target['id'], target['H']
        path = inputs/f'{idx}.txt'
        path.write_text(gram_input(H))
        result = json.loads(subprocess.check_output([str(ROOT/'h2_binary_check'), str(path)], text=True))
        assert result['complete']
        record = {'id': idx, 'det': target['det'], 'shell': result['shell_pairs'], 'all_binary_counts': result['obstruction_representations']}
        # Any explicit pair with positive determinant below d2(M) is sufficient.
        witness = None
        for i in range(4):
            for j in range(i+1, 4):
                c = tuple(H[i][j])
                delta = H[i][i][0]*H[j][j][0]-mul(c, conj(c))[0]
                if 0 < delta < 16:
                    u, v = [0]*8, [0]*8
                    u[i] = v[j] = 1
                    witness = {'kind': 'd2_violation', 'u': u, 'v': v, 'norms': [H[i][i][0], H[j][j][0]], 'inner': list(c), 'det': delta}
                    break
            if witness:
                break
        if witness is None:
            for typ, w in zip(TYPES, result['first_witnesses']):
                if w is not None:
                    q1, q2, c = typ
                    witness = {'kind': 'absent_binary', **w, 'norms': [q1, q2], 'inner': list(c), 'det': q1*q2-mul(c, conj(c))[0]}
                    break
        assert witness is not None, f'No exclusion found for target {idx}'
        # This verifier uses direct arithmetic in O, not trace-G matrix products.
        u, v = witness['u'], witness['v']
        q1, q2 = witness['norms']
        assert inner(H, u, u) == (q1, 0)
        assert inner(H, v, v) == (q2, 0)
        assert inner(H, u, v) == tuple(witness['inner'])
        assert witness['det'] > 0
        record['witness'] = witness
        records.append(record)
    (ROOT/'h2_candidate_witnesses.json').write_text(json.dumps(records, indent=2)+'\n')
    summary = {
        'candidates': len(records), 'all_excluded': True,
        'kinds': dict(collections.Counter(r['witness']['kind'] for r in records)),
        'determinants': dict(collections.Counter(r['det'] for r in records)),
        'witness_types': dict(collections.Counter(str((r['witness']['norms'],r['witness']['inner'])) for r in records)),
        'input_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'seconds': time.perf_counter()-started,
    }
    (ROOT/'h2_candidate_witnesses_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
