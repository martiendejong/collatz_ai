"""
172_log_correlation.py
=======================
IS THE CERTIFICATE FIELD LOG-CORRELATED? The two BRW fingerprints:

 F1: Var(log2 v) grows ~ linearly in k (per-level increments of a
     branching-random-walk-like field add variance).
 F2: Cov(log2 v(m), log2 v(m')) grows ~ linearly in the number of shared
     3-adic digits j = v3(m - m') (tree covariance structure), both for the
     raw field and for the residual field after removing the deterministic
     tower profile (strata means in (v3(m+1), v3(m+4))).

If both hold, the limit theorem lim gamma(k) = 1 sits in the domain of
log-correlated-field extreme value theory (Bramson, Fyodorov-Bouchaud),
and the k^-0.85 rate should be derivable from its universal corrections.
"""
import numpy as np
from math import log2

print("F1: variance of log2 v vs k")
print(f"{'k':>3} {'Var(log2 v)':>12} {'increment':>10}")
prev = None
for k in range(12, 18):
    try:
        v = np.load(f"certificate_k{k}.npy", mmap_mode='r')
    except FileNotFoundError:
        continue
    n = v.shape[0]
    idx = np.random.default_rng(0).integers(0, n, min(n, 2_000_000))
    lv = np.log2(np.asarray(v[np.sort(idx)], dtype=np.float64))
    var = float(lv.var())
    inc = var - prev if prev is not None else float('nan')
    print(f"{k:>3} {var:>12.4f} {inc:>10.4f}")
    prev = var
print()

k = 15
v = np.load(f"certificate_k{k}.npy").astype(np.float64)
v /= v.max()
N = 3**(k-1)
K3 = 3**k
lv = np.log2(v)
i_all = np.arange(N, dtype=np.int64)
m = (3*i_all + 2) % K3

def v3_arr(x, cap):
    val = np.zeros(len(x), dtype=np.int32)
    y = x.copy()
    for _ in range(cap):
        div = (y % 3 == 0) & (y > 0)
        if not div.any():
            break
        y = np.where(div, y // 3, y)
        val += div
    return val

w1 = np.minimum(v3_arr(m + 1, k), 8)
w4 = np.minimum(v3_arr(m + 4, k), 4)
# residual field: remove strata means
res = lv.copy()
for a in range(1, 9):
    for b in range(1, 5):
        sel = (w1 == a) & (w4 == b)
        if sel.any():
            res[sel] -= lv[sel].mean()

print(f"F2 (k={k}): correlation vs shared 3-adic depth j = v3(m-m')")
print(f"{'j':>3} {'corr raw':>9} {'corr resid':>11} {'cov resid':>10}")
rng = np.random.default_rng(1)
NP = 400_000
for j in range(1, 13):
    step = 3**(j-1)
    i1 = rng.integers(0, N, NP)
    t = rng.integers(1, 3**3, NP)
    t += (t % 3 == 0)                      # ensure 3 does not divide t
    i2 = (i1 + t*step) % N
    ok = ((i2 - i1) % (3*step) != 0)       # exact depth j (not deeper)
    a, b = lv[i1[ok]], lv[i2[ok]]
    ra, rb = res[i1[ok]], res[i2[ok]]
    craw = float(np.corrcoef(a, b)[0, 1])
    crs = float(np.corrcoef(ra, rb)[0, 1])
    cvr = float(np.cov(ra, rb)[0, 1])
    print(f"{j:>3} {craw:>9.4f} {crs:>11.4f} {cvr:>10.4f}")
print()
print("log-correlated verdict: covariance ~ linear in j (tree structure) --")
print("compare the resid-cov column against a straight line in j.")
