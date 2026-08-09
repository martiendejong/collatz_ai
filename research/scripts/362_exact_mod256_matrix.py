# 362: C1 from the 2026-08-09 analysis — the EXACT rational mod-256 transition
# matrix of the Syracuse macro-step (declared proof target of Obs 248/250).
# Structure: for odd r mod 256, 3n+1 = 3r+1 mod 256, so j = v2(3n+1) is
# DETERMINISTIC (= v2(3r+1)) whenever v2(3r+1) < 8; the successor
# (3n+1)/2^j mod 256 is then uniform over the 2^j lifts t: (3(r+t*256)+1)/2^j.
# The single exceptional class r* with 3r*+1 = 0 mod 256 has successor odd-part
# uniform over all 128 odd residues (unit part of a Haar-uniform 2-adic).
# Everything is an exact rational with power-of-2 denominator.
from fractions import Fraction
import numpy as np

M = 256
odds = [r for r in range(1, M, 2)]
idx = {r: i for i, r in enumerate(odds)}
n_st = len(odds)

def v2(x):
    c = 0
    while x % 2 == 0:
        x //= 2; c += 1
    return c

P = [[Fraction(0) for _ in range(n_st)] for _ in range(n_st)]
exceptional = []
for r in odds:
    j = v2(3*r + 1)
    if j >= 8:
        exceptional.append(r)
        for s in odds:
            P[idx[r]][idx[s]] = Fraction(1, n_st)
        continue
    w = Fraction(1, 2**j)
    for t in range(2**j):
        s = ((3*(r + t*M) + 1) // 2**j) % M
        P[idx[r]][idx[s]] += w

print(f"exact matrix built: {n_st} states, exceptional classes {exceptional}")
# row sums exactly 1?
assert all(sum(row) == 1 for row in P), "row sums != 1"
print("row sums exactly 1: OK")

# column sums (doubly stochastic <=> uniform stationary)
csums = [sum(P[i][jc] for i in range(n_st)) for jc in range(n_st)]
uniform_cols = all(c == 1 for c in csums)
print(f"doubly stochastic (all column sums exactly 1): {uniform_cols}")
if not uniform_cols:
    dev = sorted(set(csums))
    print(f"  distinct column sums: {[str(c) for c in dev[:6]]} ...")

# exact stationary via rational Gaussian elimination on (P^T - I) pi = 0, sum pi = 1
A_ = [[P[i][jc] - (1 if i == jc else 0) for i in range(n_st)] for jc in range(n_st)]
A_[0] = [Fraction(1) for _ in range(n_st)]  # replace one equation with normalization
b = [Fraction(0)]*n_st; b[0] = Fraction(1)
n = n_st
for col in range(n):
    piv = next(rw for rw in range(col, n) if A_[rw][col] != 0)
    A_[col], A_[piv] = A_[piv], A_[col]
    b[col], b[piv] = b[piv], b[col]
    inv = 1/A_[col][col]
    A_[col] = [x*inv for x in A_[col]]
    b[col] *= inv
    for rw in range(n):
        if rw != col and A_[rw][col] != 0:
            f = A_[rw][col]
            A_[rw] = [x - f*y for x, y in zip(A_[rw], A_[col])]
            b[rw] -= f*b[col]
pi = b
uni = Fraction(1, n_st)
exact_uniform = all(p == uni for p in pi)
print(f"exact stationary distribution uniform 1/128: {exact_uniform}")
if not exact_uniform:
    devs = [(odds[i], pi[i]) for i in range(n_st) if pi[i] != uni]
    print(f"  {len(devs)} deviating states, largest: "
          f"{max(devs, key=lambda x: abs(x[1]-uni))}")

# E[k0] under exact stationary: k0 = v2(3r+1) capped semantics — use exact j-distribution:
# for non-exceptional r: k0 = v2(3r+1) deterministic; exceptional r: E[j] = 8 + 1 = 9?
# (v2 >= 8, geometric beyond: E = 8 + 1 = 9)
Ek = Fraction(0)
for r in odds:
    j = v2(3*r + 1)
    contrib = Fraction(j) if j < 8 else Fraction(9)
    Ek += pi[idx[r]] * contrib
print(f"E[v2(3n+1)] under exact stationary = {Ek} = {float(Ek):.6f} (Obs 250 identity: 2)")

# spectral gap of the exact matrix (float eigenvalues of exact entries)
Pf = np.array([[float(x) for x in row] for row in P])
ev = np.linalg.eigvals(Pf)
ev = sorted(ev, key=lambda z: -abs(z))
print(f"eigenvalues: 1.000000 (Perron), |lam_2| = {abs(ev[1]):.6f}, |lam_3| = {abs(ev[2]):.6f}")
print(f"exact spectral gap at mod-256 level: {1-abs(ev[1]):.6f}")
