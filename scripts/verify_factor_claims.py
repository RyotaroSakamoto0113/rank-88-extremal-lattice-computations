#!/usr/bin/env python3
"""Verify the exact factor-group assertions and their supporting inputs.

Requires a fresh run of compute/groups/reproduce.py and PARI/GP plus GAP.
"""
from pathlib import Path
import argparse,json,re,subprocess,shutil,time

def matrix(a):
    return '['+';'.join(','.join(str(x) for x in r) for r in a)+']'
def run(cmd,cwd,name):
    p=subprocess.run(list(map(str,cmd)),cwd=cwd,capture_output=True,text=True)
    (cwd/(name+'.stdout')).write_text(p.stdout)
    (cwd/(name+'.stderr')).write_text(p.stderr)
    if p.returncode or any(line.lstrip().startswith('Error,') or (line.lstrip().startswith('***') and 'Warning:' not in line) for line in (p.stdout+p.stderr).splitlines()):
        raise RuntimeError('Failed '+name+'; see saved logs')
    return p.stdout

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repository',type=Path,required=True)
    ap.add_argument('--groups',type=Path,required=True)
    ap.add_argument('--gap',required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();root=a.repository.resolve();old=a.groups.resolve();out=a.output.resolve()
    out.mkdir(parents=True,exist_ok=False);start=time.perf_counter()
    D=json.loads((root/'compute/groups/input.json').read_text())
    code='default(parisizemax,1000000000);\nneed(b,s)={if(!b,error(s))};\n'
    for tag,n in [('J',1012),('J0',506)]:
        code+=f'G={matrix(D[tag+"_G"])};V=qfminim(G,12,,2);need(V[1]=={n},"shell count");need(V[2]==12,"maximum enumerated norm");need(qfminim(G,11,,2)[1]==0,"minimum norm");need(qfminim(G,14,,2)[1]==V[1],"Hermitian norm7 absent");write("shells.g",Str("{tag}shell:=",Vec(V[3][,1]),";;"));print("{tag} norm12 count=",V[1]);\n'
    code+='BC='+matrix(D['J_basis_polynomials_columns'])+';S='+matrix(D['J0_basis_in_J'])+';G0='+matrix(D['J0_G'])+';W0='+matrix(D['J0_W'])+';GE='+matrix(D['G_E'])+';WE='+matrix(D['W_E'])+';\n'
    code+=r'''
et23=Mod(x,polcyclo(23));
c=sum(k=1,8,et23^[0,1,4,6,8,9,12,20][k]);ratio=c/subst(lift(c),x,et23^-1);BC0=BC*S;
CS=BC0^-1*matrix(22,22,i,j,polcoef(lift(ratio*sum(k=1,22,BC0[k,j]*et23^(1-k))),i-1,x));
need(denominator(CS)==1 && CS^2==matid(22) && CS~*G0*CS==G0 && CS*W0==(matid(22)-W0)*CS,"displayed J0 conjugation");
write("shells.g",Str("sigmaNatural:=",vector(22,i,Vec(CS[i,])),";;"));
om=Mod(y,y^2-y+6);z=Mod(x,polcyclo(5));
pp=[1-om+om*z-2*z^3,-2-2*z-3*z^3,1-3*z-(1+om)*z^3,-om-z^2-(1+om)*z^3];
CP=matrix(4,4,i,j,polcoef(lift(pp[j]),i-1,x));
CQ=matrix(4,4,i,j,polcoef(lift(subst(subst(liftall(pp[j]),y,1-om),x,z^-1)),i-1,x));
CE=CP^-1*CQ;A=matrix(4,4,i,j,polcoef(lift(CE[i,j]),0,y));BB=matrix(4,4,i,j,polcoef(lift(CE[i,j]),1,y));
TE=concat(concat(A,A+6*BB)~,concat(BB,-A)~)~;
need(denominator(TE)==1 && TE^2==matid(8) && TE~*GE*TE==GE && TE*WE==(matid(8)-WE)*TE,"natural E conjugation");
write("shells.g",Str("tauNatural:=",vector(8,i,Vec(TE[i,])),";;"));
al=1+sum(k=1,22,if(kronecker(k,23)==1,et23^k,0));
need(lift(al^2-al+6)==0,"quadratic embedding");
need(subst(lift(al),x,Mod(21,47))==Mod(34,47),"relative prime norm");
nf=nfinit(y^2-y+6);ppF=idealhnf(nf,2,om);qF=idealhnf(nf,47,om-34);
need(idealmul(nf,ppF,qF)==idealhnf(nf,8-3*om),"Steinitz ideal identity");
need(nfeltnorm(nf,8-3*om)==94,"norm94");
need(poldisc(polcyclo(23))==-(23^21) && poldisc(y^2-y+6)==-23,"relative discriminant support");
print("EXACT_IDEAL_AND_INVOLUTIONS_VERIFIED");quit;
'''
    (out/'supplement.gp').write_text(code)
    stdout=run([shutil.which('gp'),'-fq','supplement.gp'],out,'pari')
    if 'EXACT_IDEAL_AND_INVOLUTIONS_VERIFIED' not in stdout:raise RuntimeError('PARI completion missing')
    gap='check:=function(b,s) if not b then Error(s); fi; end;;\n'
    gap+='C:=Group('+repr(D['E_generators'])+');;\n'
    for tag in ['J','J0']:
        s=(old/(tag+'_aut.gp')).read_text()
        gens=[[[int(x)for x in row.split(',')]for row in m.strip('[]').split(';')]for m in re.findall(r'\[[^\[\]]*;[^\[\]]*\]',s)]
        if not gens:raise RuntimeError('Missing factor generators '+tag)
        gap+=tag+':=Group('+repr(gens)+');;\n'
        gap+=tag+'W:='+repr(D[tag+'_W'])+';;\n'
    gap+='Read("shells.g");;\n'
    gap+=r'''
check(Size(Center(C))=2,"E center");;
check(IsomorphismGroups(DerivedSubgroup(C),SL(2,5))<>fail,"E derived SL2(5)");;
actvec:=function(v,g) return g^-1*v; end;;
iso:=fail;;permact:=fail;;bij:=();;ss:=fail;;
for record in [[J,Jshell,JW,1012,12],[J0,J0shell,J0W,506,48]] do
  G:=record[1];;v:=record[2];;W:=record[3];;K:=DerivedSubgroup(G);;
  check(Size(K)=6072 and Size(Center(G))=2,"derived and center");;
  orb:=Orbit(G,v,actvec);;check(Length(orb)=record[4],"shell transitive");;
  check(Size(Stabilizer(G,v,actvec))=record[5],"shell stabilizer");;
  check(Size(Centralizer(G,W))=12144,"F-linear centralizer order");;
  if record[4]=1012 then check(RankMat(orb)=22,"faithful shell span"); fi;
  U:=PSL(2,23);;iso:=IsomorphismGroups(K,U);;check(iso<>fail,"PSL recognition");;
  syl:=SylowSubgroup(K,23);;sylows:=AsList(ConjugacyClassSubgroups(K,syl));;
  check(Length(sylows)=24,"24 Sylow23 subgroups");;
  permact:=ActionHomomorphism(K,sylows,OnPoints);;
  check(Size(Kernel(permact))=1,"faithful Sylow action");;
  points:=[];;
  for ss in sylows do
    fix:=Filtered([1..24],x->ForAll(GeneratorsOfGroup(Image(iso,ss)),g->x^g=x));;
    check(Length(fix)=1,"unique projective fixed point");;Add(points,fix[1]);;
  od;
  check(Set(points)=[1..24],"Sylow-projective bijection");;bij:=PermList(points);;
  check(ForAll(GeneratorsOfGroup(K),g->Image(permact,g)^bij=Image(iso,g)),"equivariant projective identification");;
  Print("ACTION_VERIFIED shell=",Length(orb)," stabilizer=",record[5]," Sylow23=24\n");
od;
K:=DerivedSubgroup(J0);;
check(sigmaNatural in J0 and Order(sigmaNatural)=2,"natural J0 involution in group");;
check(Size(Group(Concatenation(GeneratorsOfGroup(K),[sigmaNatural])))=12144,"natural PGL splitting");;
check(not sigmaNatural in Group(Concatenation(GeneratorsOfGroup(K),[-IdentityMat(22)])),"natural involution outer");;
extendedE:=Group(Concatenation(GeneratorsOfGroup(C),[tauNatural]));;
check(Size(extendedE)=480 and Order(tauNatural)=2,"natural E splitting");;
Print("REVISED_FACTOR_CLAIMS_VERIFIED\n");QUIT;
'''
    (out/'supplement.g').write_text(gap)
    gr=run([a.gap,'--print-gaproot'],out,'gaproot').strip().splitlines()[-1]
    (out/'gaproot').mkdir();(out/'gaproot/gap.ini').write_text('GAPInfo.Dependencies:=rec(NeededOtherPackages:=[["gapdoc",">= 1.2"],["smallgrp",">= 1.0"]]);\n')
    result=run([a.gap,'-l',str(out/'gaproot')+';'+gr,'-A','-q','-b','supplement.g'],out,'gap')
    if 'REVISED_FACTOR_CLAIMS_VERIFIED' not in result:raise RuntimeError('GAP completion missing')
    report={'status':'verified','J_minimum':12,'J_norm12_vectors':1012,'J0_minimum':12,'Hermitian_norm7_absent_both':True,'J0_norm12_vectors':506,'shell_stabilizers':[12,48],'factor_centers':[2,2],'derived_orders':[6072,6072],'Sylow23_counts':[24,24],'Sylow_actions_faithful_and_projectively_identified':True,'natural_E_and_displayed_J0_involutions_verified':True,'Steinitz_prime_norm_identity_verified':True,'seconds':time.perf_counter()-start,'scope':'Exact factor-group and arithmetic checks.'}
    (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
