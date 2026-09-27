"""Read only explicit witnesses; do not trust the candidate search's short lists."""
from pathlib import Path
import json,collections,time
P=Path(__file__).resolve().parent
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def mul(x,y):return(x[0]*y[0]-6*x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def conj(x):return(x[0]+x[1],-x[1])
def inner(u,H,v):
 s=(0,0)
 for i in range(4):
  for j in range(4):s=add(s,mul(mul((u[i],u[i+4]),H[i][j]),conj((v[j],v[j+4]))))
 return s
def norm(x):return x[0]**2+x[0]*x[1]+6*x[1]**2
def main():
 start=time.perf_counter();H={r['id']:r['H']for r in map(json.loads,(P/'independent_h3_targets.jsonl').read_text().splitlines())}
 rows=list(map(json.loads,(P/'independent_h3_binary_sieve.jsonl').read_text().splitlines()))
 assert len(H)==len(rows)==13272 and len({r['id']for r in rows})==13272 and {r['id']for r in rows}==set(H)
 types=set();counts=collections.Counter()
 for r in rows:
  status=r['status'];counts[status]+=1
  if status=='survives':continue
  A=H[r['id']];u=r['u'];qu=inner(u,A,u);assert qu==(r['q_u'],0)
  if status=='unary_absent':assert 0<qu[0]<6 or qu[0]==7;continue
  v=r['v'];qv=inner(v,A,v);h=inner(u,A,v)
  assert qv==(r['q_v'],0) and h==tuple(r['h'])
  d=qu[0]*qv[0]-norm(h);assert d==r['det'] and d>0
  if status=='det_lt16':assert d<16
  else:
   assert status=='binary_absent' and qu[0]<=12 and qv[0]<=12
   types.add((qu[0],qv[0],*h))
 (P/'binary_obstruction_types.txt').write_text('\n'.join(' '.join(map(str,t))for t in sorted(types))+'\n')
 result={'verified':True,'candidates':len(rows),'counts':dict(counts),'distinct_binary_types':len(types),'seconds':time.perf_counter()-start}
 (P/'binary_witness_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
if __name__=='__main__':main()
