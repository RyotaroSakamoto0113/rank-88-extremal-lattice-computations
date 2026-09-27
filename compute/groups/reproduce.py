"""Exact factor reconstruction and group recognition; Python standard library."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import json,subprocess,re,shutil,sys,argparse,time
import cyclotomic_arithmetic as c
P=Path(__file__).resolve().parent
D=json.loads((P/'input.json').read_text())
def transpose(a):return list(map(list,zip(*a)))
def mul(a,b):return [[sum(x*y for x,y in zip(r,v))for v in zip(*b)]for r in a]
def eye(n):return [[int(i==j)for j in range(n)]for i in range(n)]
def det(a):
 a=[list(map(F,r))for r in a];v=F(1)
 for k in range(len(a)):
  j=next((j for j in range(k,len(a))if a[j][k]),None)
  if j is None:return 0
  if j!=k:a[k],a[j]=a[j],a[k];v=-v
  q=a[k][k];v*=q
  for j in range(k+1,len(a)):
   z=a[j][k]/q
   for l in range(k+1,len(a)):a[j][l]-=z*a[k][l]
 return v

def ldl(g):
 n=len(g);l=[[F(i==j)for j in range(n)]for i in range(n)];d=[]
 for i in range(n):
  d.append(F(g[i][i])-sum(l[i][k]**2*d[k]for k in range(i)))
  assert d[-1]>0
  for j in range(i+1,n):l[j][i]=(g[j][i]-sum(l[j][k]*l[i][k]*d[k]for k in range(i)))/d[i]
 return l,d

def roots(g):
 l,d=ldl(g);n=len(g);x=[0]*n;out=[]
 def rec(i,rem):
  if i<0:
   if rem==0:out.append(x[:])
   return
  cc=sum((l[j][i]*x[j]for j in range(i+1,n)),F(0));a,b=cc.numerator,cc.denominator
  r=rem/d[i];k=isqrt(r.numerator*b*b//r.denominator)
  for t in range((-k-a+b-1)//b,(k-a)//b+1):x[i]=t;rec(i-1,rem-d[i]*(t+cc)**2)
 rec(n-1,F(2));return out

def gp(a):return '['+';'.join(','.join(str(x)for x in r)for r in a)+']'
def run(cmd,cwd):
 q=subprocess.run(list(map(str,cmd)),cwd=cwd,capture_output=True,text=True)
 if q.returncode or re.search(r'(^|\n)(Error,|  \*\*\*   (?!Warning))',q.stdout+q.stderr):raise RuntimeError(q.stdout+'\n'+q.stderr)
 return q.stdout

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--gap',default=shutil.which('gap'));ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.perf_counter()
 eta=[F(0),F(1)]+[F(0)]*20;alpha=c.ONE[:]
 for r in c.RESIDUES:alpha=c.add(alpha,c.power(eta,r))
 t=c.add(eta,c.sigma(eta,-1));aa=c.ZERO[:]
 for x in reversed([49,54,-112,-116,66,86,-2,-27,-6,3,1]):aa=c.add(c.mul(aa,t),c.scale(c.ONE,x))
 delta=c.mul(c.power(c.add(c.scale(c.ONE,2),c.scale(t,-1)),5),aa);weight=c.inverse(delta)
 moments=[]
 for k in range(23):
  v=c.mul(weight,c.power(eta,k));moments.append(22*v[0]-sum(v[1:]))
 Gpower=[[moments[(i-j)%23]for j in range(22)]for i in range(22)]
 B=D['J_basis_polynomials_columns'];B0=mul(B,D['J0_basis_in_J']);assert abs(det(B))==47 and abs(det(B0))==47*2**11
 for tag,bas,scale in [('J',B,1),('J0',B0,F(1,2))]:
  g=D[tag+'_G'];w=D[tag+'_W'];assert mul(transpose(bas),mul(Gpower,bas))==[[F(x)/scale for x in r]for r in g]
  assert mul(bas,w)==transpose([c.mul(alpha,list(map(F,v)))for v in transpose(bas)])
  assert det(g)==23**11;ldl(g)
  for v in transpose(bas):assert sum(x*pow(21,i,47)for i,x in enumerate(v))%47==0
 rem=c.gf2_rem;gg=sum(1<<k for k in [0,1,5,6,7,9,11])
 for v in transpose(B0):assert rem(sum((int(x)%2)<<i for i,x in enumerate(v)),gg)==0
 GE=D['G_E'];W=D['W_E'];Id=eye(8);assert det(GE)==1
 R=roots(GE);assert len(R)==240
 T=mul(R,mul(GE,transpose(R)));S=mul(R,mul(GE,mul(W,transpose(R))))
 U=[[23,0,4,12],[0,23,11,4],[-4,12,23,0],[11,-4,0,23]];V=[[0,0,-8,-1],[0,0,1,-8],[8,-1,0,0],[1,8,0,0]]
 HU=[[12*T[i][j]-S[i][j]for j in range(240)]for i in range(240)];HV=[[2*S[i][j]-T[i][j]for j in range(240)]for i in range(240)]
 groups=[]
 for anti in [False,True]:
  aut=[];AW=[[Id[i][j]-W[i][j]if anti else W[i][j]for j in range(8)]for i in range(8)]
  WR=transpose(mul(AW,transpose(R)))
  def extend(ids):
   k=len(ids)
   if k==4:
    m=transpose([R[i]for i in ids]+[WR[i]for i in ids]);assert mul(transpose(m),mul(GE,m))==GE;aut.append(m);return
   for r in range(240):
    if all(HU[s][r]==U[j][k]+(V[j][k]if anti else 0) and HV[s][r]==(-V[j][k]if anti else V[j][k])for j,s in enumerate(ids)):extend(ids+[r])
  extend([]);assert len(aut)==240;groups.append(aut)
 assert len({tuple(r[0] for r in g)for g in groups[0]})==240
 assert all(g in groups[0]for g in D['E_generators']);assert D['E_anti'] in groups[1];assert mul(D['E_anti'],D['E_anti'])==Id
 pair=next((R[i],R[j])for i in range(240)for j in range(240)if HU[i][j]==8 and HV[i][j]==7)
 print('Field matrices, E roots and all 480 semilinear maps verified.',flush=True)
 gs='default(parisizemax,1000000000);\nneed(b,s)={if(!b,error(s))};\n'
 for tag,order in [('J',12144),('J0',24288)]:
  gs+=f'G={gp(D[tag+"_G"])};W={gp(D[tag+"_W"])};AA=qfauto(G);need(AA[1]=={order},"order");for(i=1,#AA[2],g=AA[2][i];need(g~*G*g==G && abs(matdet(g))==1,"isometry");need(g*W==W*g || g*W==(matid(22)-W)*g,"normalizes F"));write("{tag}_aut.gp",AA);\n'
 gs+='quit;\n';(out/'factors.gp').write_text(gs);run([shutil.which('gp'),'-fq','factors.gp'],out)
 for tag in ['J','J0']:
  s=(out/(tag+'_aut.gp')).read_text();D[tag+'_generators']=[[[int(x)for x in row.split(',')]for row in m.strip('[]').split(';')]for m in re.findall(r'\[[^\[\]]*;[^\[\]]*\]',s)]
 tau=D['J0_anti'];assert mul(tau,tau)==eye(22);assert mul(transpose(tau),mul(D['J0_G'],tau))==D['J0_G'];assert mul(tau,D['J0_W'])==mul([[int(i==j)-D['J0_W'][i][j]for j in range(22)]for i in range(22)],tau)
 if not a.gap:raise RuntimeError('GAP executable required; pass --gap /path/to/gap')
 # A local startup file requires only packages actually used here.
 gr=run([a.gap,'--print-gaproot'],out).strip().splitlines()[-1];(out/'gaproot').mkdir();(out/'gaproot/gap.ini').write_text('GAPInfo.Dependencies:=rec(NeededOtherPackages:=[["gapdoc",">= 1.2"],["smallgrp",">= 1.0"]]);\n')
 gap='check:=function(b,s) if not b then Error(s); fi; end;;\n'
 for tag,key in [('C','E_generators'),('J','J_generators'),('J0','J0_generators')]:gap+=tag+':=Group('+repr(D[key])+');;\n'
 gap+='tau:='+repr(D['E_anti'])+';; sigma:='+repr(D['J0_anti'])+';;\n'
 gap+='check(IdGroup(C)=[240,89],"C identification");;check(Size(DerivedSubgroup(C))=120 and IsPerfectGroup(DerivedSubgroup(C)),"C derived");;\n'
 gap+='check(SortedList(List(ConjugacyClasses(DerivedSubgroup(C)),q->[Size(q),Order(Representative(q)),TraceMat(Representative(q))/2]))=SortedList([[1,1,4],[1,2,-4],[20,3,-2],[30,4,0],[12,5,-1],[12,5,-1],[20,6,2],[12,10,1],[12,10,1]]),"character restriction");;\n'
 gap+='H:=Group(Concatenation(GeneratorsOfGroup(C),[tau]));;check(Size(H)=480 and Size(DerivedSubgroup(H))=120,"semilinear E");;\n'
 gap+='K:=DerivedSubgroup(J);;check(Size(J)=12144 and Size(K)=6072 and IsSimpleGroup(K),"J derived");;check(IsomorphismGroups(K,PSL(2,23))<>fail,"PSL recognition");;\n'
 gap+='check(Size(J0)=24288 and Size(Center(J0))=2 and Order(sigma)=2,"J0");;check(Size(Group(Concatenation(GeneratorsOfGroup(DerivedSubgroup(J0)),[sigma])))=12144,"J0 splitting");;\n'
 gap+='K48:=Group('+repr(D['E_order48_subgroup_generators'])+');;check(Size(K48)=48 and IsSubgroup(C,K48),"K48");;HH:=Set(List(Elements(C),g->K48^g));;check(Length(HH)=5,"five 48 subgroups");;act:=ActionHomomorphism(C,HH,OnPoints);;check(Size(Kernel(act))=2 and Size(Image(act))=120,"S5 action");;check(ForAll(Filtered(Elements(C),g->CycleStructurePerm(Image(act,g))=[1]),g->Order(g)=4),"transposition lifts");;\n'
 gap+='Print("GROUPS_VERIFIED\\n");QUIT;\n';(out/'groups.g').write_text(gap)
 result=run([a.gap,'-l',str(out/'gaproot')+';'+gr,'-A','-q','-b','groups.g'],out);(out/'gap.stdout').write_text(result);assert 'GROUPS_VERIFIED' in result
 report={'verified':True,'E_order':240,'E_smallgroup':[240,89],'E_semilinear_order':480,'E_binary_witness_coordinates':pair,'J_trace_order':12144,'J0_trace_order':24288,'seconds':time.perf_counter()-start}
 (out/'verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
