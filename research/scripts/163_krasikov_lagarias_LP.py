"""
163_krasikov_lagarias_LP.py
============================
FAITHFUL IMPLEMENTATION of the Krasikov-Lagarias system (arXiv:math/0205002,
Acta Arith. 109 (2003) 237-258), Proposition 2.1 + the LP family L_k^NT(lambda)
of (2.7)-(2.14), solved via nonlinear Perron (monotone homogeneous operator).

System I_k on classes [3^k] = {m mod 3^k : m == 2 mod 3}, alpha = log2(3):
  (D1) m == 2 (9):  phi^m(y) >= phi^{4m}(y-2) + phibar^{(4m-2)/3}(y+alpha-2)
  (D2) m == 5 (9):  phi^m(y) >= phi^{4m}(y-2)
  (D3) m == 8 (9):  phi^m(y) >= phi^{4m}(y-2) + phibar^{(2m-1)/3}(y+alpha-1)
  phibar_{k-1}^r(y) = min over the three lifts r, r+3^{k-1}, r+2*3^{k-1}.

Ansatz phi^m(y) >= c^m lambda^y  ==>  feasibility system (L1)-(L4):
  c^m <= c^{4m} L^-2  [+ cbar^{r1} L^{alpha-2}]  [+ cbar^{r3} L^{alpha-1}]
Theorem 2.2: feasible for lambda ==> pi_a(x) >> x^{log2 lambda}.

COMPACT INDEXING (exact): m = 3i+2, i in [0, N), N = 3^{k-1}:
  4m-map:        i -> (4i+2) mod N
  branch type:   i mod 3 = 0 -> D1;  1 -> D2;  2 -> D3
  D1 target:     i = 3s   -> r1 = 4s      mod N/3   (level k-1 compact)
  D3 target:     i = 3s+2 -> r3 = 2s+1    mod N/3
  cbar[r] = min(c[r], c[r+N/3], c[r+2N/3])

Feasibility test: F_lambda is monotone + positively homogeneous; the system
c <= F(c), c > 0 is feasible iff the nonlinear Perron growth rate rho(lambda)
of F is >= 1 (power iteration). gamma = log2(lambda*), lambda* = sup feasible.

CALIBRATION TARGETS (from the paper + history):
  k=2 -> 0.43 (Krasikov 1989), k=3 -> 0.48 (Wirsching), k=9 -> 0.81 (A-L 1995),
  k=11 -> 0.84 (Applegate computation, the 2003 record and current frontier).
THEN: k = 12..16 -- beyond the published frontier.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)

def make_maps(k):
    N = 3**(k-1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i + 2) % N
    s, r = np.divmod(i, 3)
    m1 = (r == 0)          # D1
    m3 = (r == 2)          # D3
    Nl = N // 3
    R1 = np.where(m1, (4*s) % Nl, 0).astype(np.int64)
    R3 = np.where(m3, (2*s + 1) % Nl, 0).astype(np.int64)
    return N, Nl, T4, m1, m3, R1, R3

def rho(lmb, maps, iters=500, v0=None):
    N, Nl, T4, m1, m3, R1, R3 = maps
    la2 = lmb**(-2.0)
    lb1 = lmb**(ALPHA - 2.0)
    lb3 = lmb**(ALPHA - 1.0)
    v = np.ones(N) if v0 is None else v0.copy()
    growth = 1.0
    for _ in range(iters):
        cbar = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:3*Nl])
        w = v[T4] * la2
        w[m1] += cbar[R1[m1]] * lb1
        w[m3] += cbar[R3[m3]] * lb3
        growth = w.max()
        v = w / growth
    return growth, v

def gamma_star(k, iters=500, tol=1e-4, verbose=False):
    maps = make_maps(k)
    lo, hi = 1.05, 1.95
    v = None
    while hi - lo > tol:
        mid = 0.5*(lo + hi)
        g, v = rho(mid, maps, iters=iters, v0=v)
        if g >= 1.0:
            lo = mid
        else:
            hi = mid
        if verbose:
            print(f"    lambda={mid:.5f} rho={g:.6f}", flush=True)
    return log2(lo)

if __name__ == "__main__":
    print("Krasikov-Lagarias LP exponents  pi(x) >> x^gamma,  gamma = log2 lambda*")
    print(f"{'k':>3} {'N=3^(k-1)':>10} {'gamma':>8}   calibration")
    calib = {2: 0.43, 3: 0.48, 9: 0.81, 11: 0.84}
    for k in [2, 3, 5, 7, 9, 10, 11]:
        g = gamma_star(k)
        note = f"(published ~{calib[k]})" if k in calib else ""
        print(f"{k:>3} {3**(k-1):>10,} {g:>8.4f}   {note}", flush=True)
    print()
    print("BEYOND THE PUBLISHED FRONTIER (k > 11):")
    for k in [12, 13, 14]:
        g = gamma_star(k, iters=400, tol=2e-4)
        print(f"{k:>3} {3**(k-1):>10,} {g:>8.4f}   <-- new territory", flush=True)
