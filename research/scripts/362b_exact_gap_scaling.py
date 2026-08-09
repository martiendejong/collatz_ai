# 362b: exact spectral-gap scaling of the Haar-averaged Syracuse kernel at
# mod 2^m, m = 8..12. Same exact construction as 362 (deterministic j = v2(3r+1)
# below m, uniform over lifts; exceptional class -> exactly uniform row).
# Decides Obs 265's scaling question with EXACT matrices instead of sampled ones.
from fractions import Fraction
import numpy as np

def v2(x):
    c = 0
    while x % 2 == 0:
        x //= 2; c += 1
    return c

for m in range(8, 13):
    M = 2**m
    odds = list(range(1, M, 2))
    idx = {r: i for i, r in enumerate(odds)}
    n_st = len(odds)
    P = np.zeros((n_st, n_st))
    csum_exact_ok = True
    col_acc = [Fraction(0) for _ in range(n_st)]
    exceptional = 0
    for r in odds:
        j = v2(3*r + 1)
        if j >= m:
            exceptional += 1
            P[idx[r], :] += 1.0/n_st
            for s in range(n_st):
                col_acc[s] += Fraction(1, n_st)
            continue
        w = Fraction(1, 2**j)
        for t in range(2**j):
            s = ((3*(r + t*M) + 1) // 2**j) % M
            P[idx[r], idx[s]] += float(w)
            col_acc[idx[s]] += w
    doubly = all(c == 1 for c in col_acc)
    ev = np.linalg.eigvals(P)
    ev = sorted(ev, key=lambda z: -abs(z))
    l2 = abs(ev[1])
    print(f"mod 2^{m}: {n_st} states, exceptional {exceptional}, "
          f"doubly stochastic EXACT: {doubly}, |lam_2| = {l2:.6f} "
          f"(~{l2*2**m/3:.3f} x 3/2^m), gap = {1-l2:.6f}", flush=True)
