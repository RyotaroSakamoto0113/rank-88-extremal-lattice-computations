#!/usr/bin/env python3
"""Independent Python integer/Fraction reconstruction of the tensor input.

No PARI library is used. Both small cyclotomic field calculations, the TeX
delta4 form, ideal47 membership, tensor entries, positivity, determinant and
the Phi115 cyclic-vector characteristic-polynomial certificate are exact.
"""
import datetime,hashlib,itertools,json,resource,time
from fractions import Fraction as Q
from pathlib import Path

def need(ok,s):
    if not ok:raise ValueError(s)

def read_matrix(p):
    t=list(map(int,p.read_text().split()));n,m=t[:2];need(len(t)==2+n*m,'matrix input size');return [t[2+i*m:2+(i+1)*m] for i in range(n)]

def tr(A):return list(map(list,zip(*A)))
def mm(A,B):return [[sum(x*y for x,y in zip(r,c)) for c in zip(*B)] for r in A]
def mv(A,v):return [sum(x*y for x,y in zip(r,v)) for r in A]
def bareiss(A,positive=False):
    B=[list(r) for r in A];n=len(B);sg=1;prev=1;piv=[]
    for k in range(n-1):
        if not B[k][k]:
            need(not positive,'positive leading principal minor');i=next((i for i in range(k+1,n) if B[i][k]),None)
            if i is None:return 0,piv
            B[k],B[i]=B[i],B[k];sg=-sg
        z=B[k][k];piv.append(z)
        if positive:need(z>0,'positive leading principal minor')
        for i in range(k+1,n):
            for j in range(k+1,n):
                t=B[i][j]*z-B[i][k]*B[k][j];need(t%prev==0,'Bareiss exact division');B[i][j]=t//prev
        for i in range(k+1,n):B[i][k]=0
        prev=z
    piv.append(B[-1][-1]);d=sg*B[-1][-1]
    if positive:need(d>0,'positive determinant')
    return d,piv

F0=(Q(0),Q(0));F1=(Q(1),Q(0));FW=(Q(0),Q(1))
def fa(x,y):return(x[0]+y[0],x[1]+y[1])
def fn(x):return(-x[0],-x[1])
def fm(x,y):a,b=x;c,d=y;return(a*c-6*b*d,a*d+b*c+b*d)
def fc(x):return(x[0]+x[1],-x[1])
def fi(x):n=fm(x,fc(x));need(n[1]==0 and n[0]!=0,'field inverse');c=fc(x);return(c[0]/n[0],c[1]/n[0])
def fs(n):return(Q(n),Q(0))

def ea(x,y):return[fa(a,b) for a,b in zip(x,y)]
def en(x):return[fn(a) for a in x]
def em(x,y):
    z=[F0]*7
    for i,a in enumerate(x):
        for j,b in enumerate(y):z[i+j]=fa(z[i+j],fm(a,b))
    for i in range(6,3,-1):
        for j in range(i-4,i):z[j]=fa(z[j],fn(z[i]))
    return z[:4]
def es(a):return[a,F0,F0,F0]
def ep(x,n):
    z=es(F1)
    while n:
        if n&1:z=em(z,x)
        x=em(x,x);n//=2
    return z
EZ=[F0,F1,F0,F0]
def ec(x):
    z=es(F0)
    for i,a in enumerate(x):z=ea(z,em(es(fc(a)),ep(EZ,(-i)%5)))
    return z
def ei(x):
    cols=[em(x,ep(EZ,i)) for i in range(4)];A=[list(row)+[F1 if i==0 else F0] for i,row in enumerate(zip(*cols))]
    for k in range(4):
        i=next(i for i in range(k,4) if A[i][k]!=F0);A[k],A[i]=A[i],A[k];v=fi(A[k][k]);A[k]=[fm(v,a) for a in A[k]]
        for i in range(4):
            if i!=k:
                v=A[i][k];A[i]=[fa(a,fn(fm(v,b))) for a,b in zip(A[i],A[k])]
    out=[row[-1] for row in A];need(em(x,out)==es(F1),'E inverse identity');return out
def et(x):
    z=fm(fs(4),x[0])
    for v in x[1:]:z=fa(z,fn(v))
    return z
def fd(A):
    out=F0;n=len(A)
    for p in itertools.permutations(range(n)):
        z=fs((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n)))
        for i in range(n):z=fm(z,A[i][p[i]])
        out=fa(out,z)
    return out

def ca(x,y):return[a+b for a,b in zip(x,y)]
def cn(x):return[-a for a in x]
def cm(x,y):
    n=len(x);z=[Q(0)]*(2*n-1)
    for i,a in enumerate(x):
        if a:
            for j,b in enumerate(y):
                if b:z[i+j]+=a*b
    for i in range(2*n-2,n-1,-1):
        if z[i]:
            for j in range(i-n,i):z[j]-=z[i]
    return z[:n]
def cp(x,n):
    z=[Q(1)]+[Q(0)]*(len(x)-1)
    while n:
        if n&1:z=cm(z,x)
        x=cm(x,x);n//=2
    return z
def ci(x):
    n=len(x);cols=[cm(x,[Q(int(i==j)) for i in range(n)]) for j in range(n)];A=[list(row)+[Q(int(i==0))] for i,row in enumerate(zip(*cols))]
    for k in range(n):
        i=next(i for i in range(k,n) if A[i][k]);A[k],A[i]=A[i],A[k];v=A[k][k];A[k]=[a/v for a in A[k]]
        for i in range(n):
            if i!=k:
                v=A[i][k]
                if v:A[i]=[a-v*b for a,b in zip(A[i],A[k])]
    return[row[-1] for row in A]

def main():
    root=Path(__file__).resolve().parent;t0=time.perf_counter();c0=resource.getrusage(resource.RUSAGE_SELF);started=datetime.datetime.now(datetime.timezone.utc).isoformat();timings=[]
    def stamp(name,begin):timings.append({'stage':name,'wall_seconds':time.perf_counter()-begin})
    inp=json.loads((root/'lattice.json').read_text());G88=inp['G'];M5=inp['M5'];M23=inp['M23'];read=lambda s:read_matrix(root/(s+'.txt'))
    JG=read('J_G');JW=read('J_W');JB=read('J_basis_power');JM=read('J_M23');P_a=read('E_P_a');P_b=read('E_P_b');A_a=read('E_M5_a');A_b=read('E_M5_b')
    t=time.perf_counter();P=[[fa(F1,fn(FW)),FW,F0,fs(-2)],[fs(-2),fs(-2),F0,fs(-3)],[F1,fs(-3),F0,fn(fa(F1,FW))],[fn(FW),F0,fs(-1),fn(fa(F1,FW))]]
    need([[P[j][i][0] for j in range(4)] for i in range(4)]==P_a and [[P[j][i][1] for j in range(4)] for i in range(4)]==P_b,'TeX E basis coefficients')
    need(fd(tr(P))==F1,'TeX basis determinant1');d=fa(fm(fs(2),FW),fs(-1));zinv=ep(EZ,4);diff=ea(EZ,en(zinv));sm=ea(EZ,zinv)
    delta=em(es(fs(Q(1,2))),ea(ea(em(es(fm(fs(11),d)),em(diff,sm)),em(es(fm(fs(-2),d)),diff)),ea(em(es(fs(-115)),sm),es(fs(115)))))
    need(ec(delta)==delta,'delta4 real');invdelta=ei(delta);H=[[et(em(invdelta,em(x,ec(y)))) for y in P] for x in P]
    need(fd(H)==fs(Q(1,529)),'E Hermitian determinant');Ha23=read('E_H4_a23');Hb23=read('E_H4_b23');need(all(H[i][j]==(Q(Ha23[i][j],23),Q(Hb23[i][j],23)) for i in range(4) for j in range(4)),'delta4 independently reproduces all h4 entries')
    for j in range(4):
        expected=em(EZ,P[j]);actual=es(F0)
        for i in range(4):actual=ea(actual,em(es((Q(A_a[i][j]),Q(A_b[i][j]))),P[i]))
        need(actual==expected,'E zeta5 multiplication in stated basis')
    stamp('independent_E_delta4_and_action',t)
    t=time.perf_counter();n=22;one=[Q(1)]+[Q(0)]*21;eta=[Q(0),Q(1)]+[Q(0)]*20;etainv=cp(eta,22);tt=ca(eta,etainv);omega=one[:]
    for k in range(1,23):
        if pow(k,11,23)==1:omega=ca(omega,cp(eta,k))
    need(ca(ca(cm(omega,omega),cn(omega)),[6*x for x in one])==[0]*22,'quadratic embedding into C23')
    B=tr(JB);need(abs(bareiss(JB)[0])==47,'ideal index47')
    for j,b in enumerate(B):
        need(sum(b[k]*pow(21,k,47) for k in range(22))%47==0,'ideal47 membership');need(mv(JB,[JW[i][j] for i in range(22)])==cm(omega,b),'J omega multiplication');need(mv(JB,[JM[i][j] for i in range(22)])==cm(eta,b),'J eta multiplication')
    coeff=[49,54,-112,-116,66,86,-2,-27,-6,3,1];aval=[Q(0)]*22
    for a in reversed(coeff):aval=ca(cm(aval,tt),[a*x for x in one])
    den=cm(cp(ca([2*x for x in one],cn(tt)),5),aval);weight=ci(den);need(cm(weight,den)==one,'J trace weight inverse')
    eta_neg=[cp(eta,(-i)%23) for i in range(22)];barB=[]
    for b in B:barB.append([sum(b[i]*eta_neg[i][j] for i in range(22)) for j in range(22)])
    weighted=[cm(weight,b) for b in B];direct=[]
    for x in weighted:
        row=[]
        for y in barB:
            z=cm(x,y);v=23*z[0]-sum(z);need(v.denominator==1,'J exact integral trace');row.append(int(v))
        direct.append(row)
    need(direct==JG,'original J trace form independently reconstructed');stamp('independent_C23_ideal_and_trace',t)
    t=time.perf_counter();WtG=mm(tr(JW),JG);J0=[[Q(11*JG[i][j]+WtG[i][j],23) for j in range(22)] for i in range(22)];J1=[[Q(JG[i][j]-2*WtG[i][j],23) for j in range(22)] for i in range(22)]
    need(all(x.denominator==1 for row in J0+J1 for x in row),'J Hermitian integral coefficients')
    for i in range(88):
        a,u=divmod(i,22)
        for j in range(88):
            b,v=divmod(j,22);p,q=H[a][b];need(G88[i][j]==(2*p+q)*J0[u][v]+(p-11*q)*J1[u][v],'tensor trace entry');need(M5[i][j]==A_a[a][b]*int(u==v)+A_b[a][b]*JW[u][v],'column zeta5 action');need(M23[i][j]==(JM[u][v] if a==b else 0),'column zeta23 action')
    need(G88==tr(G88) and all(G88[i][i]%2==0 for i in range(88)),'symmetric even integral tensor');detG,principal=bareiss(G88,positive=True);need(detG==1,'tensor unimodular');stamp('independent_tensor_entries_and_Bareiss_positivity',t)
    t=time.perf_counter();M=mm(M5,M23);v=[1]+[0]*87;cols=[]
    for k in range(89):cols.append(v);v=mv(M,v)
    K=tr(cols[:88]);detK,_=bareiss(K);need(detK!=0,'cyclic vector spans Q88')
    num=[0]*117;num[0]=1;num[1]=-1;num[115]=-1;num[116]=1;den=[0]*29;den[0]=1;den[5]=-1;den[23]=-1;den[28]=1;phi=[0]*89
    for i in range(116,27,-1):
        c=num[i];phi[i-28]=c
        for j in range(29):num[i-28+j]-=c*den[j]
    need(all(c==0 for c in num),'cyclotomic exact polynomial division');need(all(sum(phi[k]*cols[k][i] for k in range(89))==0 for i in range(88)),'Phi115 cyclic-vector annihilator')
    stamp('independent_Phi115_cyclic_certificate',t)
    (root/'cyclic_vector_certificate.json').write_text(json.dumps({'vector':[1]+[0]*87,'krylov_determinant':str(detK),'phi115_coefficients_ascending':phi,'leading_principal_minors':[str(z) for z in principal]},indent=2)+'\n')
    c1=resource.getrusage(resource.RUSAGE_SELF);out={'status':'verified','dimension':88,'G_determinant':1,'positive_definite':True,'even_integral':True,'E_basis_O_determinant':1,'E_delta4_Gram_matches_TeX':True,'J_ideal_norm':47,'J_trace_and_actions_reconstructed_from_field':True,'tensor_entries_exact':True,'column_actions_exact':True,'M5_times_M23_characteristic_polynomial':'Phi_115','cyclic_vector_determinant':str(detK),'Phi23_at21_mod43':sum(pow(21,k,43) for k in range(23))%43,'Phi23_at21_mod47':sum(pow(21,k,47) for k in range(23))%47,'input_sha256':hashlib.sha256((root/'lattice.json').read_bytes()).hexdigest(),'timing':{'started_utc':started,'wall_seconds':time.perf_counter()-t0,'user_cpu_seconds':c1.ru_utime-c0.ru_utime,'system_cpu_seconds':c1.ru_stime-c0.ru_stime,'stages':timings}}
    (root/'independent_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
