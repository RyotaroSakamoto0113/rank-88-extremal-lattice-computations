#!/usr/bin/env python3
"""Reconstruct L0 from the original ideal, then export its order-115 actions.

PARI supplies exact number-field arithmetic; every integer identity is asserted.
The independent perfection verifier subsequently checks all exported matrices.
"""
from pathlib import Path
import argparse, json, subprocess, shutil

ROOT = Path(__file__).resolve().parents[1]
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args(); out = a.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    inp = ROOT/'compute/j0/run/input'
    code = r'''
default(parisizemax,2000000000);default(parisize,256000000);
need(b,s)={if(!b,error(s))};
flat(name,A)={my(f=fileopen(Str(OUT,"/",name,".txt"),"w"),sz=matsize(A));filewrite(f,Str(sz[1]," ",sz[2]));for(i=1,sz[1],filewrite(f,strjoin(vector(sz[2],j,Str(A[i,j]))," ")));fileclose(f)};
A23(t)=t^10+3*t^9-6*t^8-27*t^7-2*t^6+86*t^5+66*t^4-116*t^3-112*t^2+54*t+49;
main()={
my(eta=Mod(x,polcyclo(23)),al,t,aa,dd,wt,BC,B0,BC0,P0,Q0,g0,w0,m23,H0,H1,om,d,z,C,A5,A5a,A5b,H4,H4a,H4b,gg,m5,mm23,m115,gen);
al=1+sum(k=1,22,if(kronecker(k,23)==1,eta^k,0));t=eta+eta^-1;aa=A23(t);dd=(2-t)^5;wt=1/(2*dd*aa);
BC=matrix(22,22,i,j,polcoef(lift(B[j]),i-1,x));
B0=vector(22,j,sum(i=1,22,S[i,j]*B[i]));BC0=matrix(22,22,i,j,polcoef(lift(B0[j]),i-1,x));
P0=matrix(22,22,i,j,polcoef(lift((al-1)*B[j]),i-1,x));
need(mathnf(concat(2*BC,P0))==mathnf(BC0),"J0 = pbar J");
gen=sum(k=1,8,eta^[0,1,4,6,8,9,12,20][k]);
Q0=matrix(22,22,i,j,polcoef(lift(gen*eta^(j-1)),i-1,x));
need(mathnf(Q0)==mathnf(BC0),"displayed J0 generator");
g0=matrix(22,22,i,j,my(v=lift(wt*B0[i]*subst(lift(B0[j]),x,eta^-1)));23*polcoef(v,0,x)-subst(v,x,1));
w0=BC0^-1*matrix(22,22,i,j,polcoef(lift(al*B0[j]),i-1,x));
m23=BC0^-1*matrix(22,22,i,j,polcoef(lift(eta*B0[j]),i-1,x));
need(g0==G0&&w0==W0,"independent J0 trace reconstruction");
need(matdet(g0)==23^11&&denominator(m23)==1,"J0 determinant/action");
H0=(11*g0+w0~*g0)/23;H1=(g0-2*w0~*g0)/23;
need(denominator(H0)==1&&denominator(H1)==1,"Hermitian integrality");
om=Mod(y,y^2-y+6);d=2*om-1;z=Mod(x,polcyclo(5));
C=matrix(4,4,i,j,polcoef(lift([1-om+om*z-2*z^3,-2-2*z-3*z^3,1-3*z-(1+om)*z^3,-om-z^2-(1+om)*z^3][j]),i-1,x));
need(matdet(C)==1,"E basis");
A5=C^-1*matrix(4,4,i,j,polcoef(lift(z*z^(j-1)),i-1,x))*C;
A5a=matrix(4,4,i,j,polcoef(lift(A5[i,j]),0,y));A5b=matrix(4,4,i,j,polcoef(lift(A5[i,j]),1,y));
H4=[1,0,4/d,om/d;0,1,(om-1)/d,4/d;-4/d,om/d,1,0;(om-1)/d,-4/d,0,1];
H4a=matrix(4,4,i,j,polcoef(lift(H4[i,j]),0,y));H4b=matrix(4,4,i,j,polcoef(lift(H4[i,j]),1,y));
gg=matrix(88,88,i,j,my(a=(i-1)\22+1,b=(j-1)\22+1,u=(i-1)%22+1,v=(j-1)%22+1);(2*H4a[a,b]+H4b[a,b])*H0[u,v]+(H4a[a,b]-11*H4b[a,b])*H1[u,v]);
m5=matrix(88,88,i,j,my(a=(i-1)\22+1,b=(j-1)\22+1,u=(i-1)%22+1,v=(j-1)%22+1);A5a[a,b]*(u==v)+A5b[a,b]*w0[u,v]);
mm23=matrix(88,88,i,j,if((i-1)\22==(j-1)\22,m23[(i-1)%22+1,(j-1)%22+1],0));m115=m5*mm23;
need(gg==GL&&denominator(gg)==1&&matdet(gg)==1,"L0 tensor Gram");
for(i=1,88,need(gg[i,i]%2==0&&matdet(gg[1..i,1..i])>0,"even and positive"));
need(denominator(m5)==1&&m5^5==matid(88)&&mm23^23==matid(88),"integral actions");
need(m5*mm23==mm23*m5&&m5~*gg*m5==gg&&mm23~*gg*mm23==gg,"commuting isometries");
need(charpoly(m115)==polcyclo(115),"order115 characteristic");
flat("G88",gg);flat("M5",m5);flat("M23",mm23);flat("J0_basis_power",BC0);flat("J0_G",g0);flat("J0_W",w0);
print("L0_INPUT_VERIFIED J0_generator1 even1 determinant1 positive1 order1151");
};
main();
quit;
'''
    pre = 'OUT='+json.dumps(str(out))+';\n'
    pre += ''.join('read('+json.dumps(str(inp/n))+');\n' for n in ['J_ambient.gp','J0_exact_input.gp'])
    (out/'build.gp').write_text(pre+code)
    p = subprocess.run([shutil.which('gp'),'-fq',str(out/'build.gp')],capture_output=True,text=True)
    (out/'build.stdout').write_text(p.stdout); (out/'build.stderr').write_text(p.stderr)
    if p.returncode or 'L0_INPUT_VERIFIED' not in p.stdout: raise RuntimeError('L0 reconstruction failed; see '+str(out))
    subprocess.run([__import__('sys').executable, str(ROOT/'scripts/export_input.py'),str(out)],check=True)
    print(p.stdout)
if __name__=='__main__': main()
