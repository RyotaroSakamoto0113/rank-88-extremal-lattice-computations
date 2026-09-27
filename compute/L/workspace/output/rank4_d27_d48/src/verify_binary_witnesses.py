"""Independently check every explicit unary/binary obstruction and candidate ID.

Use --directory DIR to select a proof workspace.  The default is this file's
parent, allowing the script to replace the legacy verifier in a fresh run.
All proof checks stay active under python -O.
"""
from pathlib import Path
import argparse,collections,json,time

def require(ok,message):
 if not ok:raise ValueError(message)
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def mul(x,y):return(x[0]*y[0]-6*x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def conj(x):return(x[0]+x[1],-x[1])
def inner(u,H,v):
 require(len(u)==len(v)==8,'coefficient vector dimension')
 s=(0,0)
 for i in range(4):
  for j in range(4):s=add(s,mul(mul((u[i],u[i+4]),H[i][j]),conj((v[j],v[j+4]))))
 return s

def norm(x):return x[0]**2+x[0]*x[1]+6*x[1]**2

def read_rows(path):return [json.loads(s)for s in path.read_text().splitlines()if s.strip()]
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--directory',type=Path,default=Path(__file__).resolve().parent);args=ap.parse_args();P=args.directory
 start=time.perf_counter();targets=read_rows(P/'independent_h3_targets.jsonl');H={r['id']:r['H']for r in targets}
 require(len(H)==len(targets),'duplicate target ID')
 rows=read_rows(P/'independent_h3_binary_sieve.jsonl');ids=[r['id']for r in rows]
 require(len(ids)==len(set(ids)),'duplicate sieve ID')
 require(set(ids)==set(H),'sieve IDs do not cover the input candidates exactly')
 types=set();counts=collections.Counter()
 for r in rows:
  status=r['status'];counts[status]+=1
  require(status in {'survives','unary_absent','det_lt16','binary_absent'},'unknown sieve status')
  if status=='survives':continue
  A=H[r['id']];u=r['u'];qu=inner(u,A,u);require(qu==(r['q_u'],0),'unary witness norm mismatch')
  if status=='unary_absent':
   # Norm 7 is excluded by the independently certified full q<=12 shell of M.
   require(0<qu[0]<6 or qu[0]==7,'unsupported unary obstruction');continue
  v=r['v'];qv=inner(v,A,v);h=inner(u,A,v)
  require(qv==(r['q_v'],0) and h==tuple(r['h']),'binary witness Gram mismatch')
  d=qu[0]*qv[0]-norm(h);require(d==r['det'] and d>0,'binary witness determinant mismatch or dependent pair')
  if status=='det_lt16':require(d<16,'binary determinant is not below 16')
  else:
   require(0<qu[0]<=12 and 0<qv[0]<=12,'binary obstruction lies outside certified shell')
   types.add((qu[0],qv[0],*h))
 out='\n'.join(' '.join(map(str,t))for t in sorted(types))
 (P/'binary_obstruction_types.txt').write_text(out+('\n'if out else''))
 result={'verified':True,'coverage_exact':True,'candidates':len(rows),'counts':dict(counts),'distinct_binary_types':len(types),'seconds':time.perf_counter()-start}
 (P/'binary_witness_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
