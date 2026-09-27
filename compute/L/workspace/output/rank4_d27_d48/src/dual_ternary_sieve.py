#!/usr/bin/env python3
"""Find and independently verify free ternary determinant <27 witnesses.

PARI is a witness search only. Failure to find a short primitive dual vector
never rejects a candidate. All rejections are checked with Python integers.
"""
import argparse, hashlib, itertools, json, math, pathlib, resource, subprocess, time

def add(x,y): return (x[0]+y[0],x[1]+y[1])
def neg(x): return (-x[0],-x[1])
def mul(x,y): return (x[0]*y[0]-6*x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def bar(x): return (x[0]+x[1],-x[1])
def det(H):
    n=len(H); z=(0,0)
    for p in itertools.permutations(range(n)):
        term=(1,0)
        for i in range(n): term=mul(term,H[i][p[i]])
        if sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))%2: term=neg(term)
        z=add(z,term)
    return z
def adj(H):
    n=len(H)
    return [[tuple((-1 if (i+j)%2 else 1)*a for a in det([[H[k][l] for l in range(n) if l!=i] for k in range(n) if k!=j])) for j in range(n)] for i in range(n)]
def norm(H,v):
    z=(0,0)
    for i in range(len(v)):
        for j in range(len(v)): z=add(z,mul(mul(v[i],H[i][j]),bar(v[j])))
    assert z[1]==0
    return z[0]
def gram(H,rows):
    return [[sumfield(mul(mul(u[i],H[i][j]),bar(v[j])) for i in range(4) for j in range(4)) for v in rows] for u in rows]
def sumfield(xs):
    z=(0,0)
    for x in xs:z=add(z,x)
    return z
def ideal_columns(v):
    return list(v)+[mul(x,(0,1)) for x in v]
def primitive_gcd(v):
    cols=ideal_columns(v); g=0
    for a,b in itertools.combinations(cols,2): g=math.gcd(g,a[0]*b[1]-a[1]*b[0])
    return abs(g)
def gpmat(H): return '['+';'.join(','.join(f'({a}+({b})*w)' for a,b in row) for row in H)+']'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('output_dir');ap.add_argument('--gp',default='/opt/homebrew/bin/gp');args=ap.parse_args()
    start=time.perf_counter();out=pathlib.Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
    targets=[json.loads(s) for s in pathlib.Path(args.input).read_text().splitlines() if s.strip()]
    script=out/'dual_search.gp';raw=out/'dual_search.stdout';err=out/'dual_search.stderr'
    header='''default(parisize,256000000);
w=Mod(t,t^2-t+6);
trf(z)=2*polcoef(lift(z),0)+polcoef(lift(z),1);
finddual(id,H)={my(A=matdet(H)*H^-1,G,R,V,v,c,F,B,q,found=0,hn,cpu);gettime();G=matrix(8,8,i,j,trf(w^(i>4)*A[(i-1)%4+1,(j-1)%4+1]*(1-w)^(j>4)));R=qfminim(G,52,,2);V=R[3];for(k=1,matsize(V)[2],v=V[,k];c=vector(4,i,v[i]+v[i+4]*w);B=matrix(2,8,i,j,polcoef(lift(c[(j-1)%4+1]*w^(j>4)),i-1));hn=mathnf(B);if(hn==matid(2),q=v~*G*v/2;print("W ",id," ",q," ",strjoin(vector(8,i,Str(v[i]))," ")," ",strjoin(vector(4,i,Str(hn[(i-1)%2+1,(i-1)\\2+1]))," "));found=1;break));cpu=gettime();print("S ",id," ",found," ",matsize(V)[2]," ",cpu);};
'''
    with script.open('w') as f:
        f.write(header)
        for t in targets: f.write(f"finddual({t['id']},{gpmat(t['H'])});\n")
        f.write('quit;\n')
    prepared=time.perf_counter();usage0=resource.getrusage(resource.RUSAGE_CHILDREN)
    with raw.open('w') as fo,err.open('w') as fe: proc=subprocess.run([args.gp,'-q',str(script)],stdout=fo,stderr=fe)
    search_end=time.perf_counter();usage1=resource.getrusage(resource.RUSAGE_CHILDREN)
    if proc.returncode or '***   at' in err.read_text(): raise RuntimeError(err.read_text()[:3000])
    byid={t['id']:t for t in targets};witnesses={};statuses={}
    for line in raw.read_text().splitlines():
        f=line.split()
        if f[0]=='W':
            ident,q,*nums=map(int,f[1:]);x=nums[:8];hn=nums[8:];assert hn==[1,0,0,1]
            t=byid[ident];H=t['H'];A=adj(H);v=[(x[i],x[i+4]) for i in range(4)]
            assert det(H)==(t['det'],0)
            assert all(sum((mul(H[i][k],A[k][j])[0] for k in range(4)))==(t['det'] if i==j else 0) and sum(mul(H[i][k],A[k][j])[1] for k in range(4))==0 for i in range(4) for j in range(4))
            assert 0<q<=26 and norm(A,v)==q and primitive_gcd(v)==1
            witnesses[ident]={'id':ident,'status':'dual_ternary_lt27','dual_coefficients':v,'kernel_equation':'sum(u_i * conjugate(lambda_i)) = 0','ternary_determinant':q,'coefficient_ideal_index':1,'ideal_Z_generators':ideal_columns(v),'pari_hnf':[[1,0],[0,1]],'adjugate':A}
            piv=next((i for i,z in enumerate(v) if z in ((1,0),(-1,0))),None)
            if piv is not None:
                rows=[]
                for i in range(4):
                    if i==piv:continue
                    row=[(0,0)]*4;row[i]=(1,0);row[piv]=neg(mul(bar(v[i]),v[piv]));rows.append(row)
                assert all(sumfield(mul(row[i],bar(v[i])) for i in range(4))==(0,0) for row in rows)
                K=gram(H,rows);assert det(K)==(q,0)
                witnesses[ident].update(kernel_basis=rows,kernel_Gram=K)
        elif f[0]=='S':
            ident,found,count,cpu=map(int,f[1:]);assert ident not in statuses;statuses[ident]={'id':ident,'found':bool(found),'short_pairs_found':count,'search_cpu_ms':cpu}
        else: raise ValueError(line)
    assert set(statuses)==set(byid) and {i for i,s in statuses.items() if s['found']}==set(witnesses)
    with (out/'dual_witnesses.jsonl').open('w') as f:
        for ident in byid:
            if ident in witnesses: f.write(json.dumps(witnesses[ident],separators=(',',':'))+'\n')
    with (out/'dual_survivors.jsonl').open('w') as f:
        for t in targets:
            if t['id'] not in witnesses: f.write(json.dumps(t,separators=(',',':'))+'\n')
    (out/'dual_candidate_timings.json').write_text(json.dumps(list(statuses.values()),indent=2)+'\n')
    end=time.perf_counter()
    summary={'input':str(pathlib.Path(args.input).resolve()),'input_sha256':hashlib.sha256(pathlib.Path(args.input).read_bytes()).hexdigest(),'total':len(targets),'rejected':len(witnesses),'survivors':len(targets)-len(witnesses),'seconds':{'prepare':prepared-start,'pari_search_wall':search_end-prepared,'pari_user_cpu':usage1.ru_utime-usage0.ru_utime,'pari_system_cpu':usage1.ru_stime-usage0.ru_stime,'independent_integer_verification_and_output':end-search_end,'total_wall':end-start},'proof_dependency':'PARI output is used only to propose witnesses. Every rejection has independently verified exact adjugate norm <27 and primitive coefficient ideal. Empty search output is not a mathematical assertion.'}
    (out/'dual_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))

if __name__=='__main__':main()
