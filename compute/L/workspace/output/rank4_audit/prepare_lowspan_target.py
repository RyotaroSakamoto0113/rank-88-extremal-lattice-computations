"""Verify an O-unimodular basis and prepare the last affine target."""
import itertools
import json
from pathlib import Path
from h2_candidate_witnesses import inner, mul, add

HERE = Path(__file__).resolve().parent
targets = [json.loads(line) for line in (HERE/'independent_h3_targets.jsonl').read_text().splitlines()]
source = next(t for t in targets if t['id'] == 11622)
H = source['H']
# Each row represents four coefficients a_i + b_i omega as [a_1,...,a_4,b_1,...,b_4].
C = [[0,0,0,1,0,0,0,0], [0,2,1,-1,0,0,0,1], [0,1,0,0,0,0,0,0], [1,0,0,-1,0,1,1,0]]
determinant = (0, 0)
for p in itertools.permutations(range(4)):
    value = ((-1)**sum(p[i] > p[j] for i in range(4) for j in range(i+1, 4)), 0)
    for i in range(4):
        value = mul(value, (C[i][p[i]], C[i][p[i]+4]))
    determinant = add(determinant, value)
assert determinant == (1, 0)
HC = [[list(inner(H, u, v)) for v in C] for u in C]
assert [HC[i][i][0] for i in range(4)] == [12, 14, 13, 13]
assert HC[0][1] == [-3, -4]
result = dict(source_id=source['id'], source_H=H, basis_rows=C, basis_determinant=determinant, H=HC)
(HERE/'h3_lowspan_rebased.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(dict(source_id=source['id'], basis_determinant=determinant, norms=[12,14,13,13])))
