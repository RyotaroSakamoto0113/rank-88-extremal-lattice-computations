"""Verify source Gram equivalences and exact coverage, for any candidate count.

Use --directory DIR to select a proof workspace; default: this file's parent.
No abstract-search implementation is imported. Checks remain active under -O.
"""
import argparse,hashlib,itertools,json,pathlib,time

def require(ok,message):
 if not ok:raise ValueError(message)
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def mul(x,y):return(x[0]*y[0]-6*x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def conj(x):return(x[0]+x[1],-x[1])
def scale(n,x):return(n*x[0],n*x[1])
def norm(x):return x[0]*x[0]+x[0]*x[1]+6*x[1]*x[1]
def transform(A,H):
 out=[]
 for u in A:
  row=[]
  for v in A:
   z=(0,0)
   for i in range(4):
    for j in range(4):z=add(z,mul(mul(u[i],H[i][j]),conj(v[j])))
   row.append(z)
  out.append(row)
 return out

def read_rows(path):return [json.loads(s)for s in path.read_text().splitlines()if s.strip()]
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--directory',type=pathlib.Path,default=pathlib.Path(__file__).resolve().parent);args=ap.parse_args();p=args.directory;start=time.monotonic()
 rows=read_rows(p/'independent_h3_targets.jsonl');targets={r['id']:r for r in rows};require(len(targets)==len(rows),'duplicate target ID')
 sieve=read_rows(p/'independent_h3_binary_sieve.jsonl');sids=[r['id']for r in sieve]
 require(len(sids)==len(set(sids)) and set(sids)==set(targets),'binary sieve coverage mismatch')
 surv={r['id']for r in sieve if r['status']=='survives'}
 m=json.loads((p/'h3_source_orbit_map.json').read_text());av=json.loads((p/'h3_source_orbit_automorphisms.json').read_text())
 require(all(len(v)==32 and all(type(x)is int for x in v)for v in av),'nonintegral or malformed basis-change matrix')
 aut=[[[(v[i+4*j],v[i+4*j+16])for j in range(4)]for i in range(4)]for v in av]
 reps=json.loads((p/'h3_source_representative_ids.json').read_text())
 require(len(reps)==len(set(reps)),'duplicate representative ID')
 require(set(reps)<=surv,'representative not in surviving candidates')
 for A in aut:
  det=(0,0)
  for per in itertools.permutations(range(4)):
   sign=(-1)**sum(per[i]>per[j]for i in range(4)for j in range(i+1,4));term=(1,0)
   for i in range(4):term=mul(term,A[i][per[i]])
   det=add(det,scale(sign,term))
  require(norm(det)==1,'source basis-change matrix not O-unimodular')
 require(all(len(v)==3 and all(type(x)is int for x in v)for v in m),'malformed orbit witness')
 mids=[v[0]for v in m];require(len(mids)==len(set(mids)) and set(mids)==surv,'orbit map does not cover all survivors exactly')
 require(set(reps)=={v[1]for v in m},'representative set differs from orbit map')
 for target,rep,ai in m:
  require(1<=ai<=len(aut),'basis-change index out of range')
  H=transform(aut[ai-1],targets[rep]['H'])
  require(H==[[tuple(x)for x in r]for r in targets[target]['H']],'source Gram equivalence fails')
 outs=[targets[i]for i in reps]
 out='\n'.join(json.dumps(t)for t in outs);(p/'independent_h3_representatives.jsonl').write_text(out+('\n'if out else''))
 out='\n'.join(str(t['id'])+' '+' '.join(str(x)for r in t['H']for e in r for x in e)for t in outs);(p/'independent_h3_representatives.flat').write_text(out+('\n'if out else''))
 status=dict(verified=True,coverage_exact=True,unimodular_witness_matrices=len(aut),targets=len(m),representatives=len(reps),seconds=time.monotonic()-start,hashes={f:hashlib.sha256((p/f).read_bytes()).hexdigest()for f in ['independent_h3_targets.jsonl','independent_h3_binary_sieve.jsonl','h3_source_orbit_map.json','h3_source_orbit_automorphisms.json','h3_source_representative_ids.json','independent_h3_representatives.jsonl']})
 (p/'h3_source_orbit_verification.json').write_text(json.dumps(status,indent=2)+'\n');print(json.dumps(status,indent=2))
if __name__=='__main__':main()
