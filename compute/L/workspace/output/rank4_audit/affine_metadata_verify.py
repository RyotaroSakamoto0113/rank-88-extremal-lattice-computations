"""Independent integer identities and Fraction enumeration for affine certificates."""
from fractions import Fraction as F
from math import isqrt


def num(x):
    return F(x) if '/' in str(x) else int(x)


def mat(x):
    return [[num(z) for z in row] for row in x]


def tr(A):
    return list(map(list, zip(*A)))


def mm(A, B):
    BT = tr(B)
    return [[sum(x*y for x, y in zip(row, col)) for col in BT] for row in A]


def mv(A, v):
    return [sum(x*y for x, y in zip(row, v)) for row in A]


def dot(u, v):
    return sum(x*y for x, y in zip(u, v))


def det(A):
    A = [list(row) for row in A]
    n, sign, previous = len(A), 1, 1
    for k in range(n-1):
        if not A[k][k]:
            j = next((i for i in range(k+1, n) if A[i][k]), None)
            if j is None:
                return 0
            A[k], A[j] = A[j], A[k]
            sign = -sign
        pivot = A[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                v = A[i][j]*pivot-A[i][k]*A[k][j]
                assert v % previous == 0
                A[i][j] = v//previous
            A[i][k] = 0
        previous = pivot
    return sign*A[-1][-1]


def solve(A, b):
    A = [[F(x) for x in row]+[F(z)] for row, z in zip(A, b)]
    n = len(A)
    for k in range(n):
        j = next(i for i in range(k, n) if A[i][k])
        A[k], A[j] = A[j], A[k]
        pivot = A[k][k]
        A[k] = [x/pivot for x in A[k]]
        for i in range(n):
            if i != k:
                scale = A[i][k]
                A[i] = [x-scale*y for x, y in zip(A[i], A[k])]
    return [row[-1] for row in A]


def verify_metadata(record, G):
    A, T, H = (mat(record[k]) for k in ('A', 'T', 'H'))
    rhs = list(map(num, record['rhs']))
    n, rows = len(T), len(A)
    assert abs(det(T)) == 1
    assert mm(A, T) == [[0]*(n-rows)+row for row in H]
    y = solve(H, rhs)
    assert record['consistent'] == all(x.denominator == 1 for x in y)
    if not record['consistent']:
        return
    K, B, Q = (mat(record[k]) for k in ('K', 'B', 'Q'))
    x0, c = ([num(x) for x in record[k]] for k in ('x0', 'center'))
    radius, bound = (num(record[k]) for k in ('radius', 'target_trace_norm'))
    assert all(isinstance(v, int) for row in K+B for v in row)
    assert all(isinstance(v, int) for v in x0)
    assert abs(det([a+b for a, b in zip(K, B)])) == 1
    assert mm(A, K) == [[0]*(n-rows) for _ in range(rows)]
    assert mm(A, B) == H
    assert mv(A, x0) == rhs
    assert mm(mm(tr(K), G), K) == Q
    assert mv(Q, c) == [-v for v in mv(mm(tr(K), G), x0)]
    assert radius == bound-dot(x0, mv(G, x0))+dot(c, mv(Q, c))


def enumerate_affine(Q, center, radius):
    """An independent exact LDL enumerator; retain both signs and the origin."""
    n = len(Q)
    Q = [[F(v) for v in row] for row in Q]
    center = list(map(F, center))
    diagonal, U = [], [[F(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        diagonal.append(Q[i][i]-sum(diagonal[k]*U[k][i]**2 for k in range(i)))
        assert diagonal[-1] > 0
        U[i][i] = F(1)
        for j in range(i+1, n):
            U[i][j] = (Q[i][j]-sum(diagonal[k]*U[k][i]*U[k][j] for k in range(i)))/diagonal[i]
    assert mm(mm(tr(U), [[diagonal[i] if i == j else 0 for j in range(n)] for i in range(n)]), U) == Q
    solutions, vector = [], [0]*n

    def recurse(i, remaining):
        if i < 0:
            solutions.append(tuple(vector))
            return
        offset = -center[i]+sum(U[i][j]*(vector[j]-center[j]) for j in range(i+1, n))
        radius2 = remaining/diagonal[i]
        a, b = offset.numerator, offset.denominator
        integer_radius = isqrt((radius2.numerator*b*b)//radius2.denominator)
        low = -((integer_radius+a)//b)
        high = (integer_radius-a)//b
        for z in range(low, high+1):
            left = remaining-diagonal[i]*(z+offset)**2
            assert left >= 0
            vector[i] = z
            recurse(i-1, left)

    if radius >= 0:
        recurse(n-1, F(radius))
    return set(solutions)
