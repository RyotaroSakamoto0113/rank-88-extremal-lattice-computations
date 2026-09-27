#!/usr/bin/env python3
"""Independent field arithmetic audit of source frames and all completed embeddings."""
import argparse,json,math,time
from fractions import Fraction as F
from pathlib import Path

def add(x,y):return(x[0]+y[0],x[1]+y[1])
def neg(x):return(-x[0],-x[1])
def mul(x,y):a,b=x;c,d=y;return(a*c-6*b*d,a*d+b*c+b*d)
def con(x):return(x[0]+x[1],-x[1])
def inv(x):a,b=x;n=a*a+a*b+6*b*b;return((a+b)/n,-b/n)
def field_inverse(c):
 a=[[tuple(map(F,z)) for z in row]+[(F(i==j),F(0)) for j in range(4)] for i,row in enumerate(c)]
 for j in range(4):
  p=next(i for i in range(j,4) if a[i][j]!=(0,0));a[p],a[j]=a[j],a[p];z=inv(a[j][j]);a[j]=[mul(z,v) for v in a[j]]
  for i in range(4):
   if i!=j:
    z=a[i][j];a[i]=[add(v,neg(mul(z,w))) for v,w in zip(a[i],a[j])]
 return [r[4:] for r in a]
def matvec(a,v):return tuple(sum(x*y for x,y in zip(r,v)) for r in a)
def sourcegram(c,h):
 return [[sum_pair(mul(mul(c[i][k],h[k][l]),con(c[j][l])) for k in range(4) for l in range(4)) for j in range(4)] for i in range(4)]
def sum_pair(it):
 a=b=0
 for x,y in it:a+=x;b+=y
 return a,b

def main():
 p=argparse.ArgumentParser();p.add_argument('--targets',required=True);p.add_argument('--results',required=True);p.add_argument('--quadruples',required=True);p.add_argument('--output',required=True);a=p.parse_args();t0=time.monotonic()
 targets={j['id']:j['H'] for j in map(json.loads,open(a.targets))};results={j['id']:j for j in map(json.loads,open(a.results))};quads=[list(map(int,l.split())) for l in open(a.quadruples)]
 data=list(map(int,open('output/rank4_audit/M_integer_data.txt').read().split()));assert data[0]==22;G=[data[1+22*i:1+22*(i+1)] for i in range(22)];W=[data[485+22*i:485+22*(i+1)] for i in range(22)]
 need={abs(k) for row in quads for k in row[2:5]};shell={}
 with open('output/rank4_audit/M_shell12.txt') as sf:
  n,dim=map(int,next(sf).split());assert dim==22
  for i,l in enumerate(sf,1):
   if i in need:shell[i]=tuple(map(int,l.split()))
 assert len(shell)==len(need)
 reps=[tuple(map(int,l.split()))[1:] for l in open('output/rank4_audit/M_shell12_orbit_reps.txt')]
 inverses={};grams={};positive=0;maxinv=0
 for i,j in results.items():
  c=[[(row[k],row[k+4]) for k in range(4)] for row in j['frame']];ci=field_inverse(c);ell=math.lcm(*(z.denominator for row in ci for pair in row for z in pair));P=[[int(z[0]*ell) for z in row] for row in ci];Q=[[int(z[1]*ell) for z in row] for row in ci]
  assert ell==j['ell'] and P==j['inverseP'] and Q==j['inverseQ'],i
  maxinv=max(maxinv,max(abs(z) for m in [P,Q] for row in m for z in row));inverses[i]=(P,Q,ell);grams[i]=sourcegram(c,targets[i])
  # Known integral X: four distinct standard basis vectors. C X and inverse recover X exactly.
  X=[tuple(int(k==r) for k in range(22)) for r in range(4)];WX=[matvec(W,x) for x in X]
  U=[tuple(sum(c[r][s][0]*X[s][k]+c[r][s][1]*WX[s][k] for s in range(4)) for k in range(22)) for r in range(4)];WU=[matvec(W,u) for u in U]
  assert all(sum(P[r][s]*U[s][k]+Q[r][s]*WU[s][k] for s in range(4))==ell*X[r][k] for r in range(4) for k in range(22)),i
  positive+=1
 cov={};omega={}
 def vectors(row):
  vv=[reps[row[1]]]
  for s in row[2:5]:vv.append(tuple((1 if s>0 else -1)*x for x in shell[abs(s)]))
  return vv
 def prepare(v):
  if v not in cov:
   wv=matvec(W,v);omega[v]=wv;cov[v]=(matvec(G,v),matvec(G,wv))
 def h(u,v):
  prepare(u);gu,gwu=cov[u];s=sum(x*y for x,y in zip(gu,v));t=sum(x*y for x,y in zip(gwu,v));assert (11*s+t)%23==0 and (s-2*t)%23==0
  return((11*s+t)//23,(s-2*t)//23)
 count={i:0 for i in results};integral_count={i:0 for i in results};bad_examples=[]
 for row in quads:
  i,r,j,k,l,claimed=row;u=vectors(row);gh=grams[i]
  assert all(h(u[s],u[t])==gh[s][t] for s in range(4) for t in range(4)),row
  for v in u:prepare(v)
  P,Q,ell=inverses[i];failure=None
  for s in range(4):
   for kk in range(22):
    z=sum(P[s][t]*u[t][kk]+Q[s][t]*omega[u[t]][kk] for t in range(4))
    if z%ell:failure=(s,kk,z,ell);break
   if failure:break
  ok=failure is None;assert ok==bool(claimed),row
  count[i]+=1;integral_count[i]+=ok
  if failure and len(bad_examples)<5:bad_examples.append({'id':i,'failure':failure})
 assert all(count[i]==j['quadruples'] and integral_count[i]==j['integral'] for i,j in results.items())
 out={'complete':True,'source_frames':len(results),'field_inverse_matches_restriction_inverse':True,'positive_integral_reconstruction_controls':positive,'quadruples_checked':len(quads),'all_gram_conditions_match':True,'integral_embeddings':sum(integral_count.values()),'maximum_inverse_integer_coefficient':maxinv,'sample_nonintegral_coordinates':bad_examples,'seconds':time.monotonic()-t0}
 Path(a.output).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
