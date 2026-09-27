from pathlib import Path
import json,time,csv
O=Path(__file__).resolve().parents[1];R=O/'results/orbit_table';t=time.perf_counter()
rows=[json.loads(s) for s in (R/'burnside_rows.jsonl').read_text().splitlines()];exact=[json.loads(s) for s in (O/'results/fast_orbit_all/stdout').read_text().splitlines()]
assert len(rows)==84 and len(exact)==45
used=[];ordinary=[[0]*22 for _ in range(3)]
for row in rows:
 if row['dimension']:
  r=exact[row['form_index']-1];used.append(row['form_index']);assert r['complete'] and r['dimension']==row['dimension'] and r['bound']==42
  assert r['histogram'][1::2]==[0]*21 and r['histogram'][2::2]==row['counts']
 else:assert row['counts']==[0]*21
 for d,v in enumerate(row['counts'],1):ordinary[row['cusp']-1][d]+=row['class_size']*v
assert sorted(used)==list(range(1,46))
for c in range(3):
 for d in range(22):assert ordinary[c][d]%12144==0;ordinary[c][d]//=12144
saturated=[r[:] for r in ordinary]
for d in range(1,22):
 for c in range(3):
  # The only possible nontrivial ideal norms are 2 and 3, since d/6<4.
  # Both primes split into the two nontrivial ideal classes of Cl(F)=C3.
  for m in [2,3]:
   if d%m==0:saturated[c][d]-=saturated[(c+1)%3][d//m]+saturated[(c+2)%3][d//m]
expected={6:[1,1,1],8:[3,4,3],9:[3,3,3],10:[12,9,12],11:[14,14,14],12:[43,45,43],13:[62,68,62],14:[158,160,158],15:[268,264,268],16:[533,541,533],17:[920,903,920],18:[1693,1661,1693],19:[2749,2748,2749],20:[4667,4722,4667],21:[7405,7499,7405]}
assert all([saturated[c][d] for c in range(3)]==v for d,v in expected.items())
with (R/'table.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['discriminant','O_F','p','pbar']);w.writerows([[d]+[saturated[c][d] for c in range(3)] for d in range(1,22)])
report={'status':'verified','method':'Burnside fixed-lattice counts; exact 128-bit enumeration with GMP preflight; remove ideal multiples of norms 2 and 3','conjugacy_classes':28,'cusps':3,'fixed_lattice_forms':45,'all_exact_histograms_agree_with_PARI':True,'ordinary_orbits':ordinary,'saturated_orbits':saturated,'table_rows_match_expected':True,'identity_signed_vector_counts':[sum(exact[k]['histogram']) for k in [0,15,30]],'wall_seconds':time.perf_counter()-t}
(R/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print('ORBIT_TABLE_VERIFIED 15 rows, 3 Steinitz classes')
