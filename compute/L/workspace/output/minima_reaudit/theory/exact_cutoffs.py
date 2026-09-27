#!/usr/bin/env python3
"""Exact finite-reduction cutoffs; no ambient minima are assumed."""
from fractions import Fraction as Q
import json
from pathlib import Path
C0=Q(30613,1029)
C=[C0,Q(46),Q(46)]
def rootfloor(q,n):
    r=0
    while (r+1)**n<=q:r+=1
    return r
result={
 'source':'https://www.math.rwth-aachen.de/homes/Oliver.Braun/documents/Braun.pdf',
 'C3_class_order_free_p_pbar':[str(c) for c in C],
 'gamma2_squared':'23/5',
 'rank2_projective_det_le15':{'line_det_bound':rootfloor(Q(23,5)*15,2)},
 'rank3':[],
 'rank4':[]
}
for t,B in [(0,26),(1,13),(2,13)]:
 pmax=rootfloor(C[t]*B,3)
 result['rank3'].append({'target_class':t,'det_bound':B,'line_det_bound':pmax,
 'rank2_det_bounds':{p:rootfloor(Q(23,5)*B*p,2) for p in range(1,pmax+1)}})
for t,B in [(0,47),(1,23),(2,23)]:
 mmax=rootfloor(529*B,4)
 result['rank4'].append({'target_class':t,'det_bound':B,'ordinary_norm_bound':mmax,
 'rank2_det_bounds_by_first_line_class':{a:{rho:rootfloor(C[(t-a)%3]*B*rho*rho,3) for rho in range(1,mmax+1)} for a in range(3)}})
result['rank4_union_D2_bound']={a:{rho:rootfloor(max(C[(t-a)%3]*B for t,B in [(0,47),(1,23),(2,23)])*rho*rho,3) for rho in range(1,13)} for a in range(3)}
result['rank4_rank3_det_bound_by_rank2_det']={B:{r:rootfloor(Q(23,5)*B*r,2) for r in range(1,68)} for B in [23,47]}
Path(__file__).with_name('exact_cutoffs.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='rank4_rank3_det_bound_by_rank2_det'},indent=2))
