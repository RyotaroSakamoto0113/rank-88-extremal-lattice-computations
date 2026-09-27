"""Exact prescribed-Gram completion with independently verified affine kernels."""
import json
import os
import shutil
import pathlib
import subprocess
import sys
import time
from affine_metadata_verify import mat, num, mv, dot, tr, verify_metadata, enumerate_affine

HERE = pathlib.Path(__file__).resolve().parent
PROJECT = HERE.parents[1]


def gp_H(H):
    return '['+';'.join(','.join(f'({a})+({b})*w' for a, b in row) for row in H)+']'


def gp_states(states):
    return '['+','.join('Mat(['+','.join('['+','.join(map(str, v))+']~' for v in state)+'])' for state in states)+']'


def ambient():
    v = iter((HERE/'h2_M_input.txt').read_text().split())
    n = int(next(v))
    next(v)
    return [[int(next(v)) for _ in range(n)] for _ in range(n)], [[int(next(v)) for _ in range(n)] for _ in range(n)]


def inner(u, v, G, W):
    t, s = dot(u, mv(G, v)), dot(u, mv(G, mv(W, v)))
    assert (12*t-s) % 23 == 0 and (2*s-t) % 23 == 0
    return (12*t-s)//23, (2*s-t)//23


def run_stage(H, states, folder, stage, independent=True):
    G, W = ambient()
    prefix = folder/f'stage{stage}'
    call = prefix.with_suffix('.gp')
    call.write_text('default(parisizemax,1000000000);\nread("output/rank4_audit/affine_stage_generator.gp");\n'
                    f'Target={gp_H(H)};\nStates={gp_states(states)};\n'
                    f'affine_stage(States,Target,"{prefix}");\nquit;\n')
    started = time.perf_counter()
    p = subprocess.run([os.environ.get('AUDIT_GP') or shutil.which('gp') or '/opt/homebrew/bin/gp', '-q', str(call)], cwd=PROJECT, text=True, capture_output=True, check=True)
    (folder/f'stage{stage}_generator.log').write_text(p.stdout+p.stderr)
    assert not any('***' in line and 'Warning' not in line for line in p.stderr.splitlines()), p.stderr
    metadata = [json.loads(s) for s in pathlib.Path(str(prefix)+'_metadata.jsonl').read_text().splitlines()]
    assert len(metadata) == len(states)
    assert {m['state_id'] for m in metadata} == set(range(1, len(states)+1))
    generated = time.perf_counter()
    p = subprocess.run([str(HERE/'exact_affine_enum'), str(prefix)+'_input.txt', '55'], text=True, capture_output=True, check=True)
    pathlib.Path(str(prefix)+'_results.jsonl').write_text(p.stdout)
    results = {row['id']: row for row in map(json.loads, p.stdout.splitlines())}
    assert len(results) == len(p.stdout.splitlines())
    assert set(results) == {m['affine_id'] for m in metadata if m['consistent']}
    enumerated = time.perf_counter()
    output = []
    for meta in metadata:
        state = states[meta['state_id']-1]
        assert meta['k'] == len(state) == stage-1
        for i, u in enumerate(state):
            for j, v in enumerate(state):
                assert inner(u, v, G, W) == tuple(H[i][j])
        expected_A = [mv(G, v) for v in state]+[mv(tr(W), mv(G, v)) for v in state]
        prescribed = [H[i][stage-1] for i in range(stage-1)]
        expected_rhs = [2*a+b for a, b in prescribed]+[a+12*b for a, b in prescribed]
        assert mat(meta['A']) == expected_A
        assert list(map(num, meta['rhs'])) == expected_rhs
        verify_metadata(meta, G)
        if not meta['consistent']:
            continue
        assert num(meta['target_trace_norm']) == 2*H[stage-1][stage-1][0]
        assert H[stage-1][stage-1][1] == 0
        res = results[meta['affine_id']]
        assert res['complete'] and len(res['solutions']) == res['count']
        assert len(set(map(tuple, res['solutions']))) == res['count']
        Q, K = mat(meta['Q']), mat(meta['K'])
        c, x0 = ([num(x) for x in meta[k]] for k in ('center', 'x0'))
        if independent:
            expected = enumerate_affine(Q, c, num(meta['radius']))
            assert expected == set(map(tuple, res['solutions'])), meta['state_id']
        for coordinates in res['solutions']:
            x = [a+b for a, b in zip(x0, mv(K, coordinates))]
            assert all(isinstance(v, int) for v in x)
            assert mv(mat(meta['A']), x) == list(map(num, meta['rhs']))
            q2 = dot(x, mv(G, x))
            assert q2 <= num(meta['target_trace_norm'])
            if q2 != num(meta['target_trace_norm']):
                continue
            new = state+[x]
            for i, u in enumerate(new):
                for j, v in enumerate(new):
                    assert inner(u, v, G, W) == tuple(H[i][j])
            output.append(new)
    finished = time.perf_counter()
    summary = dict(stage=stage, input_states=len(states), inconsistent=sum(not m['consistent'] for m in metadata),
                   affine_problems=len(results), all_sphere_solutions=sum(r['count'] for r in results.values()),
                   output_states=len(output), independent_python=independent,
                   generation_seconds=generated-started, enumeration_seconds=enumerated-generated,
                   verification_seconds=finished-enumerated, seconds=finished-started)
    pathlib.Path(str(prefix)+'_states.json').write_text(json.dumps(output)+'\n')
    pathlib.Path(str(prefix)+'_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary), flush=True)
    return output, summary


def main():
    target_file = pathlib.Path(sys.argv[1]).resolve()
    target = json.loads(target_file.read_text())
    H = target['H']
    label = sys.argv[2] if len(sys.argv) > 2 else 'affine_target'
    folder = HERE/label
    folder.mkdir(exist_ok=True)
    first_norm = H[0][0][0]
    states = []
    for line in (HERE/'M_shell12_orbit_reps.txt').read_text().splitlines():
        q, *v = map(int, line.split())
        if q == first_norm:
            states.append([v])
    assert states and first_norm <= 12
    summaries = []
    for stage in range(2, 5):
        states, summary = run_stage(H, states, folder, stage)
        summaries.append(summary)
        if not states:
            break
    result = dict(target=str(target_file), first_norm=first_norm, completed=True,
                  embeddings_found=len(states), stages=summaries)
    (folder/'summary.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
