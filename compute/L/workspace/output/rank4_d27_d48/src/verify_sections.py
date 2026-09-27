#!/usr/bin/env python3
"""Verify integer Hermitian section witnesses, independently of C++ formulas."""
import argparse,itertools,json,time
from pathlib import Path

def add(x,y):return(x[0]+y[0],x[1]+y[1])
def mul(x,y):
 a,b=x;c,d=y
 return(a*c-6*b*d,a*d+b*c+b*d)
def conj(x):a,b=x;return(a+b,-b)
def neg(x):return(-x[0],-x[1])
def det(H):
 n=len(H);z=(0,0)
 for p in itertools.permutations(range(n)):
  prod=(1,0)
  for i in range(n):prod=mul(prod,H[i][p[i]])
  if sum(p[i]>p[j]for i in range(n)for j in range(i+1,n))%2:prod=neg(prod)
  z=add(z,prod)
 return z

def h(u,v,H):
 z=(0,0)
 for i in range(4):
  for j in range(4):z=add(z,mul(mul((u[i],u[i+4]),H[i][j]),conj((v[j],v[j+4]))))
 return z

def main():
 ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('result');args=ap.parse_args();t=time.perf_counter()
 targets={}
 for l in Path(args.input).read_text().splitlines():
  a=list(map(int,l.split()));assert len(a)==33
  targets[a[0]]=[[tuple(a[1+8*i+2*j:3+8*i+2*j])for j in range(4)]for i in range(4)]
 checked=0;seen=set();counts={}
 for l in Path(args.result).read_text().splitlines():
  a=json.loads(l);i=a['id'];assert i in targets and i not in seen;seen.add(i);s=a['status'];counts[s]=counts.get(s,0)+1
  if 'witness' not in a:
   assert s=='survives';continue
  H=targets[i];w=a['witness'];G=[[h(u,v,H)for v in w]for u in w]
  assert G==[[tuple(z)for z in row]for row in a['witness_Gram']],(i,'Gram')
  d=det(G);assert d==(a['det'],0),(i,d,a['det'])
  assert s=={1:'det1_lt6',2:'det2_lt16',3:'det3_lt27'}[len(w)]
  assert 0<d[0]<{1:6,2:16,3:27}[len(w)];checked+=1
 assert seen==set(targets)
 print(json.dumps({'complete':True,'candidates':len(seen),'independent_witnesses_checked':checked,'counts':counts,'seconds':time.perf_counter()-t},indent=2))
if __name__=='__main__':main()
