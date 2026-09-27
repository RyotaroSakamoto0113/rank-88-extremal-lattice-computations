#!/usr/bin/env python3
"""Optional randomized witness finder followed by independent exact verification."""
from pathlib import Path
import argparse, json, os, shutil, subprocess, sys
ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--lattice',choices=['L','L0'],required=True)
    ap.add_argument('--trials',type=int,default=160)
    ap.add_argument('--seed',type=int,default=20260926)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    gp=shutil.which('gp');cxx=shutil.which(os.environ.get('CXX','clang++')) or shutil.which('g++')
    if not gp or not cxx: raise RuntimeError('PARI/GP and a C++17 compiler are required')
    inp=ROOT/'data'/a.lattice
    launch='default(nbthreads,1);default(parisizemax,3000000000);default(parisize,256000000);\n'
    launch+='OUT='+json.dumps(str(out))+';TRIALS='+str(a.trials)+';RNGSEED='+str(a.seed)+';\n'
    launch+='read('+json.dumps(str(inp/'lattice.gp'))+');read('+json.dumps(str(ROOT/'scripts/find_seeds.gp'))+');\n'
    (out/'launch.gp').write_text(launch)
    q=subprocess.run([gp,'-fq',out/'launch.gp'],capture_output=True,text=True)
    (out/'search.stdout').write_text(q.stdout);(out/'search.stderr').write_text(q.stderr)
    if q.returncode or 'SEED_SEARCH_COMPLETE' not in q.stdout: raise RuntimeError('Witness search failed')
    rows=[json.loads(s)['x']for s in (out/'candidates.jsonl').read_text().splitlines()]
    def save_candidates(path,rows):path.write_text(str(len(rows))+' 88\n'+''.join(' '.join(map(str,r))+'\n'for r in rows))
    save_candidates(out/'candidates.txt',rows)
    subprocess.run([cxx,'-O3','-std=c++17',ROOT/'scripts/split_rank.cpp','-o',out/'split_rank'],check=True)
    q=subprocess.run([out/'split_rank',inp/'integer_input.txt',out/'candidates.txt',out/'certificate.json','461'],capture_output=True,text=True)
    (out/'rank.stdout').write_text(q.stdout);(out/'rank.stderr').write_text(q.stderr)
    if q.returncode:raise RuntimeError('No complete rank certificate found. Increase trials or change seed; this is not a disproof of perfection.')
    cert=json.loads((out/'certificate.json').read_text())
    block=next(b for b in cert['blocks']if b['character']==0)
    selected=[cert['seeds'][i]for i in block['seed_ids']]
    save_candidates(out/'selected44.txt',selected)
    q=subprocess.run([out/'split_rank',inp/'integer_input.txt',out/'selected44.txt',out/'selected44_certificate.json','461'],capture_output=True,text=True)
    final=out/('selected44_certificate.json' if q.returncode==0 else 'certificate.json')
    subprocess.run([sys.executable,ROOT/'scripts/verify_split_certificate.py','--input',inp/'lattice.json','--certificate',final,'--output',out/'independent_verification'],check=True)
    (out/'result.json').write_text(json.dumps({'status':'verified','certificate':str(final),'trials':a.trials,'seed':a.seed},indent=2)+'\n')
    print('Exact certificate independently verified:',final)

if __name__=='__main__':main()
