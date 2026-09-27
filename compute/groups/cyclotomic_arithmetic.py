#!/usr/bin/env python3
"""Exact cyclotomic and quadratic arithmetic, using Python's standard library.

This does not use PARI, the search basis, or the saved Hermitian Gram.
It reconstructs h0 by summing the eleven Galois conjugates in Q(zeta_23).
"""
import itertools
import json
import time
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parent
N = 22
ONE = [Q(1)] + [Q(0)] * 21
ZERO = [Q(0)] * N
RESIDUES = [k for k in range(1, 23) if pow(k, 11, 23) == 1]

def add(x, y):
    return [a + b for a, b in zip(x, y)]

def scale(x, c):
    return [a * c for a in x]

def mul(x, y):
    z = [Q(0)] * 43
    for i, a in enumerate(x):
        if a:
            for j, b in enumerate(y):
                if b:
                    z[i + j] += a * b
    for i in range(42, 21, -1):
        for j in range(i - 22, i):
            z[j] -= z[i]
    return z[:22]

def power(x, n):
    z = ONE[:]
    while n:
        if n & 1:
            z = mul(z, x)
        x = mul(x, x)
        n //= 2
    return z

def inverse(x):
    cols = [mul(x, [Q(i == j) for i in range(N)]) for j in range(N)]
    a = [list(row) + [Q(i == 0)] for i, row in enumerate(zip(*cols))]
    for k in range(N):
        j = next(j for j in range(k, N) if a[j][k])
        a[k], a[j] = a[j], a[k]
        c = a[k][k]
        a[k] = [v / c for v in a[k]]
        for i in range(N):
            if i != k:
                c = a[i][k]
                if c:
                    a[i] = [u - c * v for u, v in zip(a[i], a[k])]
    y = [row[-1] for row in a]
    assert mul(x, y) == ONE
    return y

def sigma(x, k):
    z = [Q(0)] * 23
    for i, c in enumerate(x):
        z[(i * k) % 23] += c
    return [z[i] - z[22] for i in range(22)]

def gf2_rem(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a

def gf2_gcd(a, b):
    while b:
        a, b = b, gf2_rem(a, b)
    return a

def fa(x, y):
    return x[0] + y[0], x[1] + y[1]

def fm(x, y):
    a, b = x
    c, d = y
    return a * c - 6 * b * d, a * d + b * c + b * d

def det4(a):
    d = (Q(0), Q(0))
    for p in itertools.permutations(range(4)):
        s = (-1) ** sum(p[i] > p[j] for i in range(4) for j in range(i + 1, 4))
        v = (Q(s), Q(0))
        for i in range(4):
            v = fm(v, a[i][p[i]])
        d = fa(d, v)
    return d

def qrank(a):
    a = [[Q(v) for v in row] for row in a]
    r = 0
    for c in range(len(a[0])):
        i = next((i for i in range(r, len(a)) if a[i][c]), None)
        if i is None:
            continue
        a[r], a[i] = a[i], a[r]
        v = a[r][c]
        a[r] = [x / v for x in a[r]]
        for i in range(r + 1, len(a)):
            v = a[i][c]
            a[i] = [x - v * y for x, y in zip(a[i], a[r])]
        r += 1
    return r

