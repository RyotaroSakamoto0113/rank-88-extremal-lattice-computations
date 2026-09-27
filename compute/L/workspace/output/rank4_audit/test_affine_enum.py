"""Independent finite-box oracle; covers shifted spheres and exact boundaries."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json, math, random, subprocess, time

HERE=Path(__file__).resolve().parent
def inverse(A):
 n=len(A); B=[list(map(F,r))+[F(i==j) for j in range(n)] for i,r in enumerate(A)]
 for i in range(n):
  p=next(j for j in range(i,n) if B[j][i]); B[i],B[p]=B[p],B[i]
  q=B[i][i]; B[i]=[v/q for v in B[i]]
  for j in range(n):
   if j!=i:
    q=B[j][i]; B[j]=[a-q*b for a,b in zip(B[j],B[i])]
 return [r[n:] for r in B]
def oracle(A,c,r):
 if r<0:return set()
 B=inverse(A); ranges=[]
 for i in range(len(c)):
  # The continuous coordinate projection satisfies |x_i-c_i|²≤r(A^-1)_ii.
  s=r*B[i][i]; k=math.isqrt(s.numerator//s.denominator)+1
  ranges.append(range(math.floor(c[i])-k,math.ceil(c[i])+k+1))
 found=set()
 for x in product(*ranges):
  y=[F(v)-z for v,z in zip(x,c)]
  q=sum(y[i]*A[i][j]*y[j] for i in range(len(x)) for j in range(len(x)))
  if q<=r:found.add(x)
 return found
def main():
 rng=random.Random(20260910); tests=[]
 for n in (1,2,3):
  for k in range(8):
   U=[[F(0) if j<i else F(1) if i==j else F(rng.randint(-1,1),2) for j in range(n)] for i in range(n)]
   D=[F(rng.randint(1,4)) for _ in range(n)]
   A=[[sum(U[k][i]*D[k]*U[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
   c=[F(rng.randint(-3,3),rng.randint(1,4)) for i in range(n)]
   r=F(rng.randint(0,5),2)
   if k==0:r=F(-1)
   if k==1:c=[F(0)]*n;r=F(0)
   if k==2:r=sum(c[i]*A[i][j]*c[j] for i in range(n) for j in range(n))
   tests.append((A,c,r))
 tests.extend([([[F(1)]],[F(10**30,3)],F(1,9)),([[F(i==j) for j in range(4)] for i in range(4)],[F(1,2)]*4,F(1))])
 s=[]
 for A,c,r in tests:s += [f'{len(c)} {r}',' '.join(map(str,c))]+[' '.join(map(str,row)) for row in A]
 path=HERE/'affine_regression_input.txt';path.write_text('\n'.join(s)+'\n')
 started=time.perf_counter();p=subprocess.run([str(HERE/'exact_affine_enum'),str(path)],text=True,capture_output=True,check=True)
 rows=[json.loads(x) for x in p.stdout.splitlines()];assert len(rows)==len(tests)
 for row,(A,c,r) in zip(rows,tests):
  expected=oracle(A,c,r); actual={tuple(x) for x in row['solutions']}
  assert row['complete'] and len(actual)==row['count'] and actual==expected,(row,expected)
 result={'passed':True,'cases':len(rows),'method':'independent finite-box exact rational oracle','seconds':time.perf_counter()-started}
 (HERE/'affine_regression_result.json').write_text(json.dumps(result,indent=2)+'\n'); print(result)
if __name__=='__main__':main()
