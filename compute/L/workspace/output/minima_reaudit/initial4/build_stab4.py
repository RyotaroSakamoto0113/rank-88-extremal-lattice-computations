"""Extract the small stabilizers as direct GP matrices, with exact checks."""
from pathlib import Path
import hashlib
import json
import re
import time
from verify_initial4 import coords, det2, integral, mul, transpose

start_wall, start_cpu = time.perf_counter(), time.process_time()
out = Path(__file__).resolve().parent
source = out.parent / "witnesses" / "orbits"
stabs = [json.loads(line) for line in (source / "stabilizers.jsonl").read_text().splitlines()]
records = [json.loads(line) for line in (out / "records.jsonl").read_text().splitlines()]
g, w = json.loads((out / "ambient.json").read_text())
assert len(stabs) == len(records) == 78
wanted = {i for row in stabs for i in row["stabilizer_group_ids"]}
matrices = {}
count = 0
for match in re.finditer(r"\[([^\[\]]*)\]", (source / "group_matrices.gp").read_text()):
    if count in wanted:
        matrices[count] = [[int(x) for x in row.split(",")] for row in match.group(1).split(";")]
    count += 1
assert count == 12144 and set(matrices) == wanted
for a in matrices.values():
    assert len(a) == 22 and all(len(row) == 22 for row in a)
    assert mul(transpose(a), mul(g, a)) == g
    assert mul(w, a) == mul(a, w)
for row, rec in zip(stabs, records):
    rid, q, rho, cl, index, lam, indexcl, x, p, *_ = rec
    ids = row["stabilizer_group_ids"]
    assert rid == row["orbit_id"] and q == row["q"]
    assert len(ids) == len(set(ids)) and 0 in ids
    assert len(ids)*row["pair_orbit_size"] == count
    for gid in ids:
        a = matrices[gid]
        ax = [r[0] for r in mul(a, [[z] for z in x])]
        assert ax == x or ax == [-z for z in x]
        ap = coords(p, mul(a, p))
        assert integral(ap) and abs(det2(ap)) == 1

def gp_mat(a):
    return "["+";".join(",".join(map(str, row)) for row in a)+"]"

(out / "STAB4.gp").write_text("STAB4=["+",".join(
    "["+",".join(gp_mat(matrices[i]) for i in row["stabilizer_group_ids"])+"]"
    for row in stabs)+"];\n")
(out / "STAB4_IDS.gp").write_text("STAB4_IDS="+json.dumps([row["stabilizer_group_ids"] for row in stabs], separators=(",", ":"))+";\n")
(out / "stabilizers.jsonl").write_text("".join(json.dumps({**row, "matrices": [matrices[i] for i in row["stabilizer_group_ids"]]}, separators=(",", ":"))+"\n" for row in stabs))
result = {"records": len(stabs), "total_matrix_entries": sum(len(row["stabilizer_group_ids"]) for row in stabs),
          "distinct_matrices": len(matrices), "group_order": count,
          "cpu_seconds": time.process_time()-start_cpu, "wall_seconds": time.perf_counter()-start_wall,
          "sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [source / "stabilizers.jsonl", source / "group_matrices.gp", out / "STAB4.gp", out / "STAB4_IDS.gp", out / "build_stab4.py"]}}
(out / "stabilizers_verification.json").write_text(json.dumps(result, indent=2)+"\n")
print(json.dumps(result, indent=2))
