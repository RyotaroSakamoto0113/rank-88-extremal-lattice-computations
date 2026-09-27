"""Independent exhaustive h=2 abstract Gram sieve, integer arithmetic only."""
import itertools, json, math, time
from pathlib import Path
OUT=Path(__file__).parent
Z=(0,0)
W=(0,1)
D=(-1,2)
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def sub(x,y): return (x[0]-y[0],x[1]-y[1])
def neg(x): return (-x[0],-x[1])
def scale(n,x):return(n*x[0],n*x[1])
def mul(x,y):return(x[0]*y[0]-6*x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def bar(x):return(x[0]+x[1],-x[1])
def norm(x):return x[0]*x[0]+x[0]*x[1]+6*x[1]*x[1]
def trace(x):return 2*x[0]+x[1]
def dot(x,y):
 r=Z
 for a,b in zip(x,y):r=add(r,mul(a,b))
 return r
def rowmul(r,H):return tuple(dot(r,[H[i][j] for i in range(len(r))]) for j in range(len(H)))
def hgram(r,H,s):return dot(rowmul(r,H),[bar(a) for a in s])
def qgram(r,H):
 q=hgram(r,H,r)
 assert q[1]==0
 return q[0]
def det3(H,ids):
 a,b,c=ids
 return H[a][a][0]*H[b][b][0]*H[c][c][0]-H[a][a][0]*norm(H[b][c])-H[b][b][0]*norm(H[a][c])-H[c][c][0]*norm(H[a][b])+trace(mul(mul(H[a][b],H[b][c]),H[c][a]))
def transform(R,H):return [[hgram(r,H,s) for s in R]for r in R]
def residue(x):return (x[0]+12*x[1])%23
def enum_norm(B):
 ans=[]
 for b in range(-math.isqrt(4*B//23)-1,math.isqrt(4*B//23)+2):
  for a in range(-math.isqrt(B)-abs(b)-1,math.isqrt(B)+abs(b)+2):
   if norm((a,b))<=B:ans.append((a,b))
 return ans
C1=((-4,0),(-1,1));C2=(W,(-4,0))
charts=[]
for v in json.loads((OUT/'charts_scaled23.json').read_text()):
 charts.append([[(v[i+4*j],v[i+4*j+16])for j in range(4)]for i in range(4)])
assert len(charts)==30
triples=list(itertools.combinations(range(4),3));pairs=list(itertools.combinations(range(4),2))
E={}
def normvals(B,r):
 key=B,r
 if key not in E:E[key]=[x for x in enum_norm(B)if residue(x)==r]
 return E[key]
def extensions(K,C,n,delta):
 r0=rowmul(C,K);b=K[0][0][0];a=K[1][1][0]
 out=[]
 for u in normvals(n*b-16,residue(r0[0])):
  for v in normvals(n*a-16,residue(r0[1])):
   # r adj(K) bar(r)^t
   value=a*norm(u)+b*norm(v)-trace(mul(mul(u,K[0][1]),bar(v)))
   ell=n*delta-value
   if ell>=552:out.append(((u,v),ell))
 return out
def checkcharts(H,s):
 # Scaled by 529 because R is stored as 23R.
 for R in charts:
  Q=transform(R,H)
  ns=[Q[i][i][0] for i in range(4)]
  if any(Q[i][i][1] or ns[i]%529 or ns[i]<6*529 or ns[i]>28*529 or ns[i]==7*529 for i in range(4)):return False
  if min(ns[0]+ns[1],ns[2]+ns[3])<s*529:return False
  if any(ns[i]*ns[j]-norm(Q[i][j])<16*529**2 for i,j in pairs):return False
  if any(det3(Q,ids)<552*529**3 for ids in triples):return False
 return True
start=time.monotonic()
counts=dict(cores=0,branches=0,extension_branches=0,extension_pairs=0,k_candidates=0,positive=0,initial_minors=0,charts=0)
results=[]
for a in range(6,12):
 if a==7:continue
 for b in range(a,24-a):
  if b==7:continue
  s=a+b
  for c in enum_norm(a*b-16):
   counts['cores']+=1
   delta=a*b-norm(c)
   if delta*(46-s)<1104:continue
   K=[[(b,0),c],[bar(c),(a,0)]]
   q1=qgram(C1,K);q2=qgram(C2,K)
   for n1 in range(6,47-s-6):
    n2=46-s-n1
    if n1==7 or n2==7 or (n1-q1)%23 or(n2-q2)%23:continue
    counts['branches']+=1
    ex1=extensions(K,C1,n1,delta);ex2=extensions(K,C2,n2,delta)
    if not ex1 or not ex2:continue
    counts['extension_branches']+=1
    for (r1,l1),(r2,l2) in itertools.product(ex1,ex2):
     counts['extension_pairs']+=1
     # k = h(y1,t2)+h(t1,y2)-h(t1,t2) modulo 23 O.
     k0=sub(add(dot(r1,[bar(x)for x in C2]),dot(C1,[bar(x)for x in r2])),hgram(C1,K,C2))
     adjK=[[(a,0),neg(c)],[neg(bar(c)),(b,0)]]
     p12=hgram(r1,adjK,r2)
     for k in enum_norm(n1*n2-16):
      if (k[0]-k0[0])%23 or(k[1]-k0[1])%23:continue
      counts['k_candidates']+=1
      num=l1*l2-norm(sub(scale(delta,k),p12))
      if num<=0:continue
      assert num%delta==0
      detY=num//delta
      assert detY%529==0
      if detY//529>33:continue
      counts['positive']+=1
      H=[[(n1,0),k,*r1],[bar(k),(n2,0),*r2],[bar(r1[0]),bar(r2[0]),(b,0),c],[bar(r1[1]),bar(r2[1]),bar(c),(a,0)]]
      if any(H[i][i][0]*H[j][j][0]-norm(H[i][j])<16 for i,j in pairs)or any(det3(H,ids)<552 for ids in triples):continue
      counts['initial_minors']+=1
      if not checkcharts(H,s):continue
      counts['charts']+=1
      results.append(dict(detX=detY//529,HY=H))
payload=dict(counts=counts,seconds=time.monotonic()-start,candidates=results)
(OUT/'independent_h2_abstract.json').write_text(json.dumps(payload,indent=2))
print(json.dumps({k:v for k,v in payload.items()if k!='candidates'},indent=2))
