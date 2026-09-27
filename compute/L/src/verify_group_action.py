"""Identify the 24-point action explicitly with PSL(2,23), without GAP."""
from pathlib import Path
import json,time
O=Path(__file__).resolve().parents[1];R=O/'results/group';t=time.perf_counter()
data=json.loads((R/'sylow_action.json').read_text());gens=[tuple(g) for g in data['generator_permutations']]
identity=tuple(range(24));group=[identity];seen={identity}
for p in group:
 for g in gens:
  q=tuple(g[p[i]] for i in range(24))
  if q not in seen:seen.add(q);group.append(q)
assert len(group)==6072
def cycles(p):
 seen=set();r=[]
 for i in range(24):
  if i in seen:continue
  c=[];j=i
  while j not in seen:seen.add(j);c.append(j);j=p[j]
  r.append(c)
 return sorted(r,key=len)
translation=next(p for p in group if [len(c) for c in cycles(p)]==[1,23])
c=cycles(translation);labels=[None]*24;labels[c[0][0]]=23
for i,j in enumerate(c[1]):labels[j]=i
inverse=[labels.index(i) for i in range(24)]
def nullvector(A):
 A=[r[:] for r in A];piv=[];i=0
 for j in range(4):
  k=next((k for k in range(i,3) if A[k][j]%23),None)
  if k is None:continue
  A[i],A[k]=A[k],A[i];u=pow(A[i][j]%23,-1,23);A[i]=[x*u%23 for x in A[i]]
  for k in range(3):
   if k!=i:
    u=A[k][j];A[k]=[(x-u*y)%23 for x,y in zip(A[k],A[i])]
  piv.append(j);i+=1
  if i==3:break
 assert len(piv)==3
 free=next(j for j in range(4) if j not in piv);v=[0]*4;v[free]=1
 for i,j in enumerate(piv):v[j]=-A[i][free]%23
 return v
cert=[]
for g in gens:
 q=[labels[g[inverse[i]]] for i in range(24)];A=[]
 for x in [23,0,1]:
  y=q[x]
  A.append(([0,0,1,0] if y==23 else [1,0,-y,0]) if x==23 else ([0,0,x,1] if y==23 else [x,1,-y*x,-y]))
 a,b,c,d=nullvector(A);det=(a*d-b*c)%23;assert det and pow(det,11,23)==1
 for x in range(24):
  num,den=(a,c) if x==23 else ((a*x+b)%23,(c*x+d)%23)
  y=23 if den==0 else num*pow(den,-1,23)%23;assert y==q[x]
 cert.append({'permutation':list(g),'mobius_matrix':[[a,b],[c,d]],'determinant':det})
report={'status':'verified','full_group_order':12144,'sylow_action_order':6072,'kernel_order':2,'projective_line_labels':labels,'translation_permutation':translation,'generators':cert,'identification':'PSL(2,23)','full_identification':'C2 x PSL(2,23)','splitting_argument':'The kernel of the 24-point action is {I,-I}. The F-linear determinant maps the finite group to the roots of unity {1,-1}; det_F(-I)=(-1)^11=-1. The product of this character and the 24-point action is therefore an injection into C2 x PSL(2,23), between groups of the same order.','wall_seconds':time.perf_counter()-t}
(R/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
