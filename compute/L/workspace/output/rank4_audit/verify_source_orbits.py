import hashlib,json,pathlib,time
p=pathlib.Path(__file__).parent
source=(p/'independent_h2_abstract.py').read_text().split('start=time.monotonic()')[0]
scope={'__file__':str(p/'independent_h2_abstract.py')};exec(compile(source,str(p/'independent_h2_abstract.py'),'exec'),scope)
transform=scope['transform']
mul,add,scale,norm=map(scope.get,('mul','add','scale','norm'))
start=time.monotonic()
targets={r['id']:r for r in map(json.loads,(p/'independent_h3_targets.jsonl').read_text().splitlines())}
surv={r['id']for r in map(json.loads,(p/'independent_h3_binary_sieve.jsonl').read_text().splitlines())if r['status']=='survives'}
m=json.loads((p/'h3_source_orbit_map.json').read_text())
av=json.loads((p/'h3_source_orbit_automorphisms.json').read_text())
aut=[[[(v[i+4*j],v[i+4*j+16])for j in range(4)]for i in range(4)]for v in av]
reps=json.loads((p/'h3_source_representative_ids.json').read_text())
import itertools
for A in aut:
 det=(0,0)
 for per in itertools.permutations(range(4)):
  sign=(-1)**sum(per[i]>per[j]for i in range(4)for j in range(i+1,4))
  term=(1,0)
  for i in range(4):term=mul(term,A[i][per[i]])
  det=add(det,scale(sign,term))
 assert norm(det)==1

assert len(m)==len({v[0]for v in m})==1377
assert {v[0]for v in m}==surv
assert set(reps)=={v[1]for v in m}
for target,rep,ai in m:
 H=transform(aut[ai-1],targets[rep]['H'])
 assert H==[[tuple(x)for x in r]for r in targets[target]['H']]
outs=[targets[i]for i in reps]
(p/'independent_h3_representatives.jsonl').write_text('\n'.join(json.dumps(t)for t in outs)+'\n')
(p/'independent_h3_representatives.flat').write_text('\n'.join(str(t['id'])+' '+' '.join(str(x)for r in t['H']for e in r for x in e)for t in outs)+'\n')
status=dict(verified=True,unimodular_witness_matrices=len(aut),targets=len(m),representatives=len(reps),seconds=time.monotonic()-start,hashes={f:hashlib.sha256((p/f).read_bytes()).hexdigest()for f in ['independent_h3_targets.jsonl','independent_h3_binary_sieve.jsonl','h3_source_orbit_map.json','h3_source_orbit_automorphisms.json','independent_h3_representatives.jsonl']})
(p/'h3_source_orbit_verification.json').write_text(json.dumps(status,indent=2))
print(json.dumps(status,indent=2))
