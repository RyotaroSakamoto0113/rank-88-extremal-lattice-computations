from pathlib import Path
import importlib.util, json, hashlib, time, resource, datetime, re
from fractions import Fraction as Q
BASE=Path(__file__).resolve().parents[1];OUT=BASE/'results/additional';OUT.mkdir(exist_ok=True);PROJECT=BASE/'workspace'
INPUT=PROJECT/'output/perfect_reaudit/input'
spec=importlib.util.spec_from_file_location('exact',INPUT/'verify_tensor_input.py');e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
start=time.perf_counter();cpu=resource.getrusage(resource.RUSAGE_SELF)
read=lambda name:e.read_matrix(INPUT/(name+'.txt'))
Ha,Hb=read('E_H4_a23'),read('E_H4_b23')
H=[[(Q(a,23),Q(b,23)) for a,b in zip(ra,rb)] for ra,rb in zip(Ha,Hb)]
s=e.fa(e.fm(e.fs(2),e.FW),e.fs(-1))
U=[[e.fm(s,x) for x in row] for row in H]
assert all(a.denominator==b.denominator==1 for row in U for a,b in row)
assert e.fd(U)==e.F1
assert e.fd([row[:3] for row in H[:3]])==e.fs(Q(1,23))
# Derive the displayed binary witness from the original trace definition.
z=e.EZ;w=e.FW
sm=e.ea(z,e.ep(z,4));diff=e.ea(z,e.en(e.ep(z,4)))
delta=e.em(e.es(e.fs(Q(1,2))),e.ea(e.ea(e.em(e.es(e.fm(e.fs(11),s)),e.em(diff,sm)),e.em(e.es(e.fm(e.fs(-2),s)),diff)),e.ea(e.em(e.es(e.fs(-115)),sm),e.es(e.fs(115)))))
uu=e.ea(e.en(e.ea(z,e.ep(z,2))),e.em(e.es(w),e.ea(e.es(e.F1),e.ea(z,e.ep(z,3)))))
vv=e.ea(e.ea(e.es(e.fs(2)),e.en(e.ep(z,2))),e.ea(e.ep(z,3),e.em(e.es(w),z)))
Hw=[[e.et(e.em(e.ei(delta),e.em(x,e.ec(y)))) for y in (uu,vv)] for x in (uu,vv)]
assert Hw[0][0]==Hw[1][1]==e.F1 and e.fd(Hw)==e.fs(Q(5,23))
# Check the cosine expansion of the fixed input polynomial A23(t).
def padd(x,y):return [(x[i] if i<len(x) else 0)+(y[i] if i<len(y) else 0) for i in range(max(len(x),len(y)))]
cs=[[2],[0,1]]
for k in range(2,11):cs.append(padd([0]+cs[-1],[-a for a in cs[-2]]))
co=[13,-1,-4,-1,6,5,-5,0,4,3,1];poly=[13]
for k in range(1,11):poly=padd(poly,[co[k]*a for a in cs[k]])
assert poly==[49,54,-112,-116,66,86,-2,-27,-6,3,1]
certpath=PROJECT/'output/perfect_reaudit/runs/20260913-124726-split-u0rotrqc/certificate.json'
seeds=json.loads(certpath.read_text())['seeds'];G=json.loads((INPUT/'lattice.json').read_text())['G']
assert all(len(v)==88 and all(type(x)==int for x in v) for v in seeds)
norms=[sum(a*b for a,b in zip(v,e.mv(G,v))) for v in seeds]
assert set(norms)=={8}
(OUT/'norm8_witnesses.json').write_text(json.dumps({'basis':'p_i tensor b_j, i=1..4 then j=1..22; J_basis_power.txt specifies b_j','gram_sha256':hashlib.sha256((INPUT/'lattice.json').read_bytes()).hexdigest(),'vectors':seeds,'trace_norms':norms},indent=2)+'\n')
# Four-square identity, valid independently of the questionable general proposition.
Z=e.F0;I=e.F1;w=e.FW
D=[[s,Z,e.fs(-4),e.fa(w,e.fs(-1))],[Z,s,w,e.fs(-4)],[Z,Z,I,Z],[Z,Z,Z,I]]
product=[[e.F0 for j in range(4)] for i in range(4)]
for i in range(4):
 for j in range(4):
  for k in range(4):product[i][j]=e.fa(product[i][j],e.fm(D[k][i],e.fc(D[k][j])))
assert product==[[e.fm(e.fs(23),a) for a in row] for row in H]
assert e.fm(e.fd(D),e.fc(e.fd(D)))==e.fs(529)
report={'status':'verified','E_dual_equals_sqrt_minus23_E':True,'E_determinant':'1/529','E_free_upper_witnesses':['1','5/23','1/23','1/529'],'E_projective_lower_proof':'Paper Proposition 3.6: trace determinant integrality, projective binary Hermite bound, projective duality','J_weight_cosine_matches_polynomial':poly,'all_44_seed_norms_exactly_8':True,'four_square_identity_exact':True,'rank4_h_at_most2_excluded_by_d4_48':48>Q(529*2**4,4**4),'wall_seconds':time.perf_counter()-start,'user_cpu_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime-cpu.ru_utime,'system_cpu_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_stime-cpu.ru_stime}
(OUT/'additional_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
