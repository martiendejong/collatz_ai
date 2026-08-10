# 365: Route A step 1 — the q=1 rung of Lemma alpha (per-cell Chebyshev).
# Setting (script 311/312 convention): cells live on the cb index space s in Z/Nl;
# U = cell means of vbar3 = (v(s)+v(s+Nl)+v(s+2Nl))/3, W = cell means of
# G = vbar3 - cb. Obs 506 proved b_0 > 0 (q=0, main classes) via the pair form
# of Cov with (3b) as engine. The induction programme needs the same at q=1.
#
# Structure used here (derived from sibling-lag): all three lifts of s share
# class r = s mod 3; the tower transport is
#   rho*vbar3(s) = A*vbar3(sigma4(s)) + [r=0] B1*cb3(R1'(s)) + [r=2] B3*cb3(R3'(s))
# with sigma4(s) = (4s+2) mod Nl and cb3 = average of cb over Nl/3-translates.
# Parts:
#  (1) verify this tower identity exactly (machine precision) on cached vectors;
#  (2) q=1 pair decomposition: per parent cell d0 = s mod 3, sub-triples over
#      digit-1 d: U[d0,d], W[d0,d]; the three Chebyshev pairs per parent, signs
#      and margins across the (lambda, k) grid;
#  (3) the "(3b)-at-scale" candidate: association of the INPUT triples,
#      Cov_d(cell means of cb3, cell means of Gcb3) where Gcb3 is the gap of the
#      next tower level; margins across the grid.
import numpy as np
import os
from math import log2

ALPHA = log2(3.0)
CACHE = r"E:\projects\collatz\research\cache"

def load_or_make(lam, k, iters=1200):
    N = 3**(k-1)
    fn = os.path.join(CACHE, f"v_lam{lam:.2f}_k{k}.npy")
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_, r_ = np.divmod(i, 3)
    Nl = N//3
    m0, m2 = (r_ == 0), (r_ == 2)
    R1 = (4*s_) % Nl; R3 = (2*s_+1) % Nl
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    if os.path.exists(fn):
        v = np.load(fn).astype(np.float64)
        it = 300
    else:
        v = np.ones(N); it = iters
    rho = 1.0
    for _ in range(it):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        w[m2] += B3*cb[R3[m2]]
        w[m0] += B1*cb[R1[m0]]
        rho = float(w.max()); v = w/rho
    return v, rho, A, B1, B3

def analyse(lam, k):
    v, rho, A, B1, B3 = load_or_make(lam, k)
    N = v.size; Nl = N//3; Nll = Nl//3
    s = np.arange(Nl, dtype=np.int64)
    vb = (v[:Nl] + v[Nl:2*Nl] + v[2*Nl:])/3.0
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    G = vb - cb
    # tower level 2
    cb3 = (cb[:Nll] + cb[Nll:2*Nll] + cb[2*Nll:])/3.0
    cbb = np.minimum(np.minimum(cb[:Nll], cb[Nll:2*Nll]), cb[2*Nll:])
    Gc = cb3 - cbb
    # (1) tower identity check
    sig4 = (4*s+2) % Nl
    sl = s // 3
    R1p = (4*sl) % Nll; R3p = (2*sl+1) % Nll
    r = s % 3
    pred = A*vb[sig4]
    pred[r == 0] += B1*cb3[R1p[r == 0]]
    pred[r == 2] += B3*cb3[R3p[r == 2]]
    err = float(np.abs(rho*vb - pred).max()/np.abs(vb).max())
    # (2) q=1 sub-triples per parent
    # cell (d0, d): s mod 9 = d0 + 3*d?  digit-0 = s mod 3 (parent), digit-1 = (s//3) mod 3
    d0 = s % 3; d1 = (s // 3) % 3
    U = np.zeros((3, 3)); W = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            m = (d0 == a) & (d1 == b)
            U[a, b] = vb[m].mean()
            W[a, b] = G[m].mean()
    pair_info = []
    for a in range(3):
        u = U[a]; w = W[a]
        prods = []
        for (p, q_) in [(0, 1), (1, 2), (0, 2)]:
            prods.append((u[p]-u[q_])*(w[p]-w[q_]))
        cov = np.cov(u, w, bias=True)[0, 1]
        pair_info.append((cov, prods))
    # (3) input association at cb-level: sub-triples of cb3 and Gc over digit-1 of tower-2 index
    sl2 = np.arange(Nll)
    e0 = sl2 % 3
    uin = np.array([cb3[e0 == b].mean() for b in range(3)])
    win = np.array([Gc[e0 == b].mean() for b in range(3)])
    cov_in = np.cov(uin, win, bias=True)[0, 1]
    norm = np.std(uin)*np.std(win) + 1e-300
    return err, pair_info, cov_in/norm, U, W

print(f"{'lam':>5} {'k':>3} {'tower-id err':>12}  {'cov(d0=0)':>10} {'cov(d0=1)':>10} {'cov(d0=2)':>10}  {'neg pairs':>9}  {'input-assoc':>11}")
results = {}
for lam in [1.05, 1.30, 1.70, 2.00]:
    for k in [10, 12, 13, 14]:
        try:
            err, pinfo, ain, U, W = analyse(lam, k)
        except Exception as e:
            print(f"{lam:>5} {k:>3}  SKIP ({e})")
            continue
        covs = [p[0] for p in pinfo]
        negp = sum(1 for _, prods in pinfo for x in prods if x < 0)
        results[(lam, k)] = (covs, ain)
        print(f"{lam:>5} {k:>3} {err:>12.2e}  {covs[0]:>10.2e} {covs[1]:>10.2e} {covs[2]:>10.2e}  {negp:>9}  {ain:>11.4f}", flush=True)
print()
print("alle per-parent covarianties positief:",
      all(all(c > 0 for c in covs) for covs, _ in results.values()))
print("input-associatie (genormaliseerd) overal positief:",
      all(a > 0 for _, a in results.values()))
