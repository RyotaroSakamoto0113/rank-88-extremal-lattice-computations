"""Independent rational/integer verification of GP-generated first-line data.

No algebraic number theory or normal-form library is used.  Saturation follows
from the gcd of all maximal minors being 1, together with equality of Q-spans.
The Steinitz class is checked using explicit p or pbar ideal bases, independently
of the quotient-eigenvalue determination.
"""
from fractions import Fraction as F
from functools import reduce
from math import gcd
from pathlib import Path
import collections
import hashlib
import json
import sys
import time


def transpose(a):
    return list(map(list, zip(*a)))


def mul(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def det2(a):
    return a[0][0]*a[1][1]-a[0][1]*a[1][0]


def content(a):
    return reduce(gcd, (abs(a[i][0]*a[j][1]-a[j][0]*a[i][1])
                        for i in range(22) for j in range(i+1, 22)), 0)


def coords(p, b):
    for i in range(22):
        for j in range(i+1, 22):
            a = [p[i], p[j]]
            d = det2(a)
            if d:
                inv = [[F(a[1][1], d), F(-a[0][1], d)],
                       [F(-a[1][0], d), F(a[0][0], d)]]
                c = mul(inv, [b[i], b[j]])
                assert mul(p, c) == b, "rational span differs"
                return c
    raise AssertionError("rank less than two")


def integral(a):
    return all(F(x).denominator == 1 for row in a for x in row)


def inside(p, v):
    return integral(coords(p, [[x] for x in v]))


def read_records(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def main():
    start_cpu, start_wall = time.process_time(), time.perf_counter()
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    g, w = json.loads((out / "ambient.json").read_text())
    records = read_records(out / "records.jsonl")
    original_gp = (out.parent / "witnesses" / "orbits" / "orbit_representatives.gp").read_text().strip()
    original = json.loads(original_gp.removeprefix("MINIMA_ORBIT_REPS=").removesuffix(";"))
    assert len(records) == 78
    assert g == transpose(g)
    assert mul(transpose(w), g) == [[g[i][j]-mul(g, w)[i][j] for j in range(22)] for i in range(22)]
    w2 = mul(w, w)
    assert all(w2[i][j]-w[i][j]+6*(i == j) == 0 for i in range(22) for j in range(22))
    stats = collections.Counter()
    nonfree = []
    for expected_id, rec in enumerate(records):
        (rid, q, rho, cl, index, lam, indexcl, x, p, k, c, y, b,
         vden, bc, wp) = rec
        assert rid == expected_id
        assert [rid, q, x] == [original[rid][0], original[rid][1], original[rid][3]]
        assert len(x) == 22 and len(p) == 22
        xc = [[z] for z in x]
        wx = [row[0] for row in mul(w, xc)]
        assert k == [[x[i], wx[i]] for i in range(22)]
        assert mul(transpose(xc), mul(g, xc))[0][0] == 2*q
        assert content(p) == 1
        assert coords(p, k) == c and integral(c)
        assert abs(det2(c)) == index == content(k)
        assert index in (1, 2)
        assert det2(mul(transpose(p), mul(g, p))) == 23*rho*rho
        assert q == index*rho and rho >= 6
        assert coords(p, mul(w, p)) == wp and integral(wp)
        if index == 1:
            assert (cl, lam, indexcl, vden) == (0, -1, -1, 1)
            assert b == k
        else:
            assert q == 12 and rho == 6 and vden == 2
            assert inside(p, y) and not inside(k, y) and inside(k, [2*z for z in y])
            wy = [row[0] for row in mul(w, [[z] for z in y])]
            assert inside(k, [wy[i]-lam*y[i] for i in range(22)])
            assert not inside(k, [wy[i]-(1-lam)*y[i] for i in range(22)])
            assert (indexcl, cl) == ((1, 2) if lam == 0 else (2, 1))
            # v=x/2.  Basis of p*v is [2v,omega*v]; basis of pbar*v
            # is [2v,(omega-1)*v].  Explicit equality fixes inverse-class direction.
            second = [F(wx[i]-(cl == 2)*x[i], 2) for i in range(22)]
            assert b == [[x[i], second[i]] for i in range(22)]
            nonfree.append({"id": rid, "q": q, "rho": rho,
                            "class": cl, "omega_quotient_eigenvalue": lam})
        assert integral(b)
        assert coords(p, b) == bc and integral(bc) and abs(det2(bc)) == 1
        stats[(q, rho, cl, index)] += 1
    result = {"verified_records": len(records), "nonfree": nonfree,
              "groups": [{"q": q, "rho": rho, "class": cl, "index": index, "count": n}
                         for (q, rho, cl, index), n in sorted(stats.items())],
              "cpu_seconds": time.process_time()-start_cpu,
              "wall_seconds": time.perf_counter()-start_wall,
              "sha256": {name: hashlib.sha256((out/name).read_bytes()).hexdigest()
                         for name in ["ambient.json", "records.jsonl", "FIRST4.gp",
                                      "build_initial4.gp", "verify_initial4.py"]}}
    (out / "verification.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
