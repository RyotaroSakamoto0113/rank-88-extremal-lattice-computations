#!/usr/bin/env python3
"""Verify actual J upper witnesses using only Python integers and Fraction.

The search uses the right omega action; this checker reconstructs h via
the left omega action. Field determinants are expanded by permutations;
trace determinants are independently computed by rational elimination.
"""
import argparse, datetime, hashlib, itertools, json, resource, time
from fractions import Fraction as Q
from pathlib import Path

def need(ok, text):
    if not ok:
        raise ValueError(text)

def mul(x, y):
    a,b=x; c,d=y
    return (a*c-6*b*d, a*d+b*c+b*d)

def add(x,y): return (x[0]+y[0],x[1]+y[1])
def conj(x): return (x[0]+x[1],-x[1])
def mm(A,B): return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def tr(A): return list(map(list,zip(*A)))
def mv(A,v): return [sum(a*b for a,b in zip(row,v)) for row in A]
def bil(x,A,y): return sum(a*b for a,b in zip(x,mv(A,y)))

def det(A):
    a=[[Q(z) for z in row] for row in A]; n=len(a); answer=Q(1)
    for i in range(n):
        p=next((j for j in range(i,n) if a[j][i]),None)
        if p is None: return Q(0)
        if p!=i: a[i],a[p]=a[p],a[i]; answer=-answer
        pivot=a[i][i]; answer*=pivot
        for j in range(i+1,n):
            c=a[j][i]/pivot
            for k in range(i+1,n): a[j][k]-=c*a[i][k]
    return answer

def det_field(H):
    out=(0,0); n=len(H)
    for perm in itertools.permutations(range(n)):
        inv=sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
        v=((-1)**inv,0)
        for i in range(n): v=mul(v,H[i][perm[i]])
        out=add(out,v)
    return out

def verify(input_path, witness_path):
    data=list(map(int,input_path.read_text().split())); need(data[0]==22,"rank22 input")
    need(len(data)>=1+2*22**2,"complete ambient input")
    G=[data[1+i*22:1+(i+1)*22] for i in range(22)]
    W=[data[1+22**2+i*22:1+22**2+(i+1)*22] for i in range(22)]
    need(G==tr(G),"symmetric trace Gram")
    for i in range(22): need(G[i][i]%2==0,"even trace diagonal")
    WW=mm(W,W)
    need(all(WW[i][j]-W[i][j]+6*(i==j)==0 for i in range(22) for j in range(22)),"omega satisfies w^2-w+6")
    GW=mm(G,W); WtG=mm(tr(W),G)
    need(all(GW[i][j]+WtG[i][j]==G[i][j] for i in range(22) for j in range(22)),"conjugate adjoint omega action")
    dG=det(G); need(dG==23**11,"ambient trace determinant")
    # Exact LDL; no floating point decisions are used.
    U=[[Q(int(i==j)) for j in range(22)] for i in range(22)]; diagonal=[]
    for i in range(22):
        di=Q(G[i][i])-sum(diagonal[k]*U[k][i]**2 for k in range(i))
        need(di>0,"ambient positive definiteness"); diagonal.append(di)
        for j in range(i+1,22): U[i][j]=(G[i][j]-sum(diagonal[k]*U[k][i]*U[k][j] for k in range(i)))/di
    records=[json.loads(line) for line in witness_path.read_text().splitlines() if line.strip()]
    need(sorted(r['rank'] for r in records)==[1,2,3,4],"exactly one witness for each rank")
    result=[]
    for rec in records:
        r=rec['rank']; vectors=rec['coordinates']; target=[6,16,27,48][r-1]
        need(len(vectors)==r and all(len(v)==22 and all(type(c) is int for c in v) for v in vectors),"integer J coordinates")
        H=[]
        for x in vectors:
            row=[]
            for y in vectors:
                s=bil(x,G,y); t=bil(mv(W,x),G,y)
                need((11*s+t)%23==0 and (s-2*t)%23==0,"integral Hermitian pairing")
                row.append(((11*s+t)//23,(s-2*t)//23))
            H.append(row)
        need(H==[[tuple(z) for z in row] for row in rec['gram']],"Gram agrees with independent left-action reconstruction")
        need(all(H[j][i]==conj(H[i][j]) for i in range(r) for j in range(r)),"Hermitian symmetry")
        df=det_field(H); need(df==(target,0),"exact desired Hermitian determinant")
        leading=[det_field([row[:i] for row in H[:i]]) for i in range(1,r+1)]
        need(all(b==0 and a>0 for a,b in leading),"positive Hermitian leading minors")
        # K columns x_1,omega*x_1,... . Positive trace determinant certifies
        # rank_Z(K)=2r, hence F independence and freeness over O_F.
        columns=[u for x in vectors for u in (x,mv(W,x))]
        TG=[[bil(x,G,y) for y in columns] for x in columns]
        td=det(TG); need(td==23**r*target**2,"independent trace determinant identity and rank2r")
        expected=[]
        for i in range(r):
            row0=[]; row1=[]
            for j in range(r):
                a,b=H[i][j]
                row0 += [2*a+b,a+12*b]
                row1 += [a-11*b,6*(2*a+b)]
            expected.extend([row0,row1])
        need(TG==expected,"all trace Gram entries agree with field trace")
        result.append({'rank':r,'hermitian_determinant':target,'trace_determinant':int(td),'norms':[H[i][i][0] for i in range(r)],'gram':H,'positive_leading_minors':[a for a,b in leading],'integer_coordinates':True,'O_linear_independence':True,'trace_rank':2*r,'trace_gram':TG})
    return {'status':'verified','claim':'upper bounds only; lower bounds require the separate complete enumerations','ambient_dimension':22,'ambient_trace_determinant':int(dG),'input_sha256':hashlib.sha256(input_path.read_bytes()).hexdigest(),'witnesses_sha256':hashlib.sha256(witness_path.read_bytes()).hexdigest(),'witnesses':result}

def main():
    root=Path(__file__).resolve().parent
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=root/'M_integer_data.txt');p.add_argument('--witnesses',type=Path,default=root/'witnesses.jsonl');p.add_argument('--output',type=Path,default=root/'verification.json');a=p.parse_args()
    started=datetime.datetime.now(datetime.timezone.utc).isoformat();t0=time.perf_counter();c0=resource.getrusage(resource.RUSAGE_SELF)
    out=verify(a.input,a.witnesses)
    c1=resource.getrusage(resource.RUSAGE_SELF);out['timing']={'started_utc':started,'wall_seconds':time.perf_counter()-t0,'user_cpu_seconds':c1.ru_utime-c0.ru_utime,'system_cpu_seconds':c1.ru_stime-c0.ru_stime};a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'determinants':[r['hermitian_determinant'] for r in out['witnesses']],'timing':out['timing']}))

if __name__=='__main__':main()
