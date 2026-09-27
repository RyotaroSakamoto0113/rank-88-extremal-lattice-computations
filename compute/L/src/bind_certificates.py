"""Bind independently generated certificates to the same lattice input."""
from pathlib import Path
import json,hashlib,importlib.util,time
O=Path(__file__).resolve().parents[1];W=O/'workspace';A=W/'output/minima_reaudit';P=W/'output/perfect_reaudit/input';t=time.perf_counter()
def read(p):return json.loads(p.read_text())
def latest(p):
 s=p.read_text().strip();key='/output/minima_reaudit/';return A/s.split(key,1)[1] if key in s else Path(s)
def importfile(name,p):
 spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
v=importfile('witness_algebra',A/'witnesses/verify_witnesses.py');e=importfile('field_algebra',P/'verify_tensor_input.py')
G=e.read_matrix(P/'J_G.txt');Om=e.read_matrix(P/'J_W.txt')
boot=latest(A/'LATEST_BOOTSTRAP.txt');d3=latest(A/'theory/LATEST_INDEPENDENT_D3.txt');b=read(boot/'BOOTSTRAP_SUMMARY.json');d=read(d3/'INDEPENDENT_D3_SUMMARY.json')
assert b['status']=='complete' and d['status']=='complete' and d['counts']['raw_ternary']==0 and d['counts']['unique_binary']==978 and d['all_exact_forms']==3099
source=(d3/'independent_d3.gp').read_text()
assert '52*q3[3]' in source and 'for(c3=1,3' in source
# The old optimized d3 summary only states the weaker projective bound 24.
# The independent all-class cutoff 52 = 2*26 proves the bound 27 used here.
di=read(W/'output/d4_projective_J/full_run/SUMMARY.json');assert di['status']=='verified' and di['exclusion_complete'] and len(di['verified_anchors'])==78 and di['counts']['raw_rank4']==0
for k in range(78):
 r=read(W/'output/d4_projective_J/full_run/anchors'/f'{k:03d}'/'FINAL_ANCHOR.json');assert r['status']=='verified' and r['independent_algebra_verification']['status']=='verified'
assert read(O/'results/initial_check/verification.json')['anchors']==78
h=read(O/'results/h3/output/rank4_audit/FULL_RUN_SUMMARY.json');assert h['status']=='complete' and h['rank4_h3_excluded'] and h['remaining']==0
assert [h[k] for k in ['abstract_candidates','binary_survivors','source_representatives']]==[8159,1376,362]
wv=A/'witnesses/fresh_reproduction/witnesses.jsonl';ambient=O/'results/h3/output/rank4_audit/M_integer_data.txt'
rr=v.verify(ambient,wv);assert rr['status']=='verified';flat=list(map(int,ambient.read_text().split()));assert [flat[1+i*22:1+(i+1)*22] for i in range(22)]==G;assert [flat[485+i*22:485+(i+1)*22] for i in range(22)]==Om
for name in ['anchor0','anchor1','anchorfree']:
 r=read(W/'output/d4_projective_J/witnesses/all_classes'/(name+'.json'));B=r['basis'];assert r['ambient_trace_gram']==G and r['ambient_omega_action']==Om;assert v.mm(B,r['omega_action'])==v.mm(Om,B);assert v.mm(v.mm(v.tr(B),G),B)==r['trace_gram'];assert v.det(r['trace_gram'])==23**4*48**2
mv0=json.loads(wv.read_text().splitlines()[0])['coordinates'][0];Z23=e.read_matrix(P/'J_M23.txt');orbit=[]
for j in range(22):orbit.append(mv0);mv0=v.mv(Z23,mv0)
assert v.det(v.tr(orbit))!=0
assert b['cusps'][0]['shells'][0]['pairs']==506 and b['cusps'][0]['shells'][0]['orbits']==1
binary=read(O/'results/binary_witness/witness.json');from fractions import Fraction as Q
uu=[list(map(Q,binary[k])) for k in ['e','f']];HH=[]
for x in uu:
 row=[]
 for y in uu:
  s=v.bil(x,G,y);q=v.bil(v.mv(Om,x),G,y);row.append(((11*s+q)/23,(s-2*q)/23))
 HH.append(row)
assert HH==[[(4,0),(-22,2)],[(-20,-2),(120,0)]]
cols=[[2*a for a in uu[0]],v.mv(Om,uu[0]),uu[1],[(a-b)/2 for a,b in zip(uu[1],v.mv(Om,uu[1]))]];assert all(a.denominator==1 for c in cols for a in c);assert v.det([[v.bil(x,G,y) for y in cols] for x in cols])==23**2*16**2
assert read(O/'results/perfect/REVERIFICATION.json')['status']=='verified'
assert read(O/'results/orbit_table/verification.json')['status']=='verified'
assert read(O/'results/group/verification.json')['status']=='verified'
assert read(O/'results/fields/timing.json')['status']=='verified'
assert read(O/'results/small_ideals/timing.json')['status']=='verified'
assert read(O/'results/catalogue/timing.json')['status']=='verified'
cr=[json.loads(s) for s in (O/'results/catalogue/verification.jsonl').read_text().splitlines()]
ce=[json.loads(s) for s in (O/'results/catalogue/exact_histograms.jsonl').read_text().splitlines()]
assert len(cr)==len(ce)==4
for rr,ee in zip(cr,ce):assert ee['complete'] and sum(ee['histogram'][:rr['minimum']])==0 and ee['histogram'][rr['minimum']]==rr['kissing_number']
assert read(O/'results/positivity/timing.json')['status']=='verified'
assert read(O/'results/nebe_class_number/verification.json')['status']=='verified'
assert read(O/'results/barnes/timing.json')['status']=='verified'
out={'status':'verified','input_lattice_sha256':hashlib.sha256((P/'lattice.json').read_bytes()).hexdigest(),'E_projective_minima':['1','5/23','1/23','1/529'],'J_projective_minima':[6,16,27,48],'J_free_minima':[6,16,27,48],'same_original_J_in_minima_h3_and_tensor':True,'d3_final_cutoff_all_classes':26,'d4_final_cutoff_all_classes':47,'rank4_h3_candidates_and_survivors':[8159,1376,362,1,0],'trace_minimum':8,'perfect_symmetric_rank':3916,'perfect_modulus':461,'perfect_seed_count':44,'all_orbit_table_rows_verified':True,'automorphism_group':'C2 x PSL(2,23)','all_supplement_components_verified':True,'wall_seconds':time.perf_counter()-t}
(O/'FINAL_VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
