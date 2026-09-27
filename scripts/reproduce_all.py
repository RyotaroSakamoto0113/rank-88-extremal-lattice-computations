#!/usr/bin/env python3
"""Run all complete computations in a new directory. Requires PARI, C++/GMP, GAP."""
from pathlib import Path
import argparse, json, shutil, subprocess, sys, time
ROOT = Path(__file__).resolve().parents[1]

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', required=True, type=Path)
    ap.add_argument('--gap', default=shutil.which('gap'))
    args = ap.parse_args()
    if not args.gap or not shutil.which('gp'):
        raise RuntimeError('PARI/GP and GAP (SmallGrp and GAPDoc) must be installed')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter(); stages = []
    commands = [
        ('L', [ROOT/'compute/L/REPRODUCE_FULL.py','--output',out/'L']),
        ('J0', [ROOT/'compute/j0/reproduce.py','--run-dir',out/'J0']),
        ('groups', [ROOT/'compute/groups/reproduce.py','--output',out/'groups','--gap',args.gap]),
        ('factor_claims', [ROOT/'scripts/verify_factor_claims.py','--repository',ROOT,'--groups',out/'groups','--gap',args.gap,'--output',out/'factor_claims']),
        ('L0_input', [ROOT/'scripts/build_L0.py','--output',out/'L0_input']),
        ('L_perfection', [ROOT/'scripts/verify_split_certificate.py','--input',ROOT/'data/L/lattice.json','--certificate',ROOT/'certificates/L_perfection.json','--output',out/'L_perfection']),
        ('L0_perfection', [ROOT/'scripts/verify_split_certificate.py','--input',out/'L0_input/lattice.json','--certificate',ROOT/'certificates/L0_perfection.json','--output',out/'L0_perfection']),
        ('formulas', [ROOT/'scripts/verify_formulas.py']),
    ]
    for name, command in commands:
        t = time.perf_counter(); print('START',name,flush=True)
        with (out/(name+'.stdout')).open('w') as so, (out/(name+'.stderr')).open('w') as se:
            p = subprocess.run([sys.executable]+list(map(str,command)), stdout=so, stderr=se)
        stages.append({'name':name,'returncode':p.returncode,'wall_seconds':time.perf_counter()-t})
        (out/'stages.json').write_text(json.dumps(stages,indent=2)+'\n')
        if p.returncode:
            raise RuntimeError(name+' failed; see logs in '+str(out))
    def read(path): return json.loads((out/path).read_text())
    assert read('L/REPRODUCTION.json')['status'] == 'fully_reproduced'
    assert read('L/FINAL_VERIFICATION.json')['status'] == 'verified'
    j0 = read('J0/output/rank4_audit/FULL_RUN_SUMMARY.json')
    assert j0['verified'] and not j0['rank4_h3_exists']
    assert read('groups/verification.json')['verified']
    assert read('factor_claims/verification.json')['status'] == 'verified'
    assert read('L_perfection/verification.json')['status'] == 'verified'
    assert read('L0_perfection/verification.json')['status'] == 'verified'
    # Exact equality binds the newly rebuilt J0 tensor to the published certificate.
    assert read('L0_input/lattice.json') == json.loads((ROOT/'data/L0/lattice.json').read_text())
    # Bind the reproduced L lattice to the distributed input.
    assert read('L/workspace/output/perfect_reaudit/input/lattice.json') == json.loads((ROOT/'data/L/lattice.json').read_text())
    result = {'status':'fully_reproduced','wall_seconds':time.perf_counter()-start,'stages':stages}
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__ == '__main__':
    main()
