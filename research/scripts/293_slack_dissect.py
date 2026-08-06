"""
293_slack_dissect.py
====================
Dissect the required inequality  S_ens2(t^2/lam, w) <= S_ens0(t, w)  (Obs 490)
into the two factors:
  weight step:    S_ens2(t^2/lam, w) / S_ens2(t, w)     (<= 1 by monotonicity, rigorous)
  ensemble step:  S_ens2(t, w)      / S_ens0(t, w)      (open — measure it)

Also test the cb-dominated linear regime for class 2:
  S(a, w) ~ a * D  with  D = E[x_{argmin y} - min x]  as a -> 0.
If class 2 sits in this regime at a = t^2/lam, then s2 ~ (t^2) * D2 * lam^{...}
and the comparison becomes quantitative.

Ensembles (per Obs 490):
  ens0: x = col_v2(4m), y = col_cb(4m)          [class-0 decomposition]
  ens2: x = col_v1(4m+3) = t*col_v0(..), y = col_cb(2m+1)   [class-2 decomposition]
NB: for ens2 we use x = col_v1 directly (weight t in the recursion), matching
S_ens2(a, b) with the ensemble (col_v1-pass, cb); s2 = S_ens2(t, lam*w) then, and via
homogeneity s2 = lam * S_ens2(t/lam, w). CAREFUL: two equivalent parameterizations:
  (i)  pass = col_v1, weights (t, lam*w)  ->  s2 = lam*S_i(t/lam, w)
  (ii) pass = col_v0, weights (t^2, lam*w) -> s2 = lam*S_ii(t^2/lam, w)
We use (i) here: S_i's pass ensemble col_v1 = t x col_v0 permuted, so
S_i(a, b) on col_v1 = S_ii(a*t, b) on col_v0. Both dissections reported.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)

def run(lam, k, niters):
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    N = 3**(k-1); Nl = N//3; Nl3 = Nl//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr == 0; m2 = r_arr == 2
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    del i, s_arr, r_arr
    v = np.ones(N)
    rho = 1.0
    for _ in range(niters):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w_ = A*v[T4]
        w_[m2] += B3*cb[R3[m2]]
        w_[m0] += B1*cb[R1[m0]]
        rho = float(w_.max())
        w_ /= rho
        v = w_
    del T4, R1, R3, m0, m2
    t = A/rho
    w = t*lam**ALPHA          # cb-weight of class 0
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    v0 = v[0::3]; v1 = v[1::3]; v2 = v[2::3]

    m = np.arange(Nl3)
    e = np.array([0, 1, 2])
    s_idx = m[:, None] + e[None, :]*Nl3

    def cols(vec, idx):
        return np.stack([vec[idx[:, 0]], vec[idx[:, 1]], vec[idx[:, 2]]], axis=1)

    # ensembles
    tgt0 = (4*s_idx) % Nl
    x0 = cols(v2, tgt0); y0 = cols(cb, tgt0)          # ens0
    tgt2p = (4*s_idx+3) % Nl
    tgt2c = (2*s_idx+1) % Nl
    x2 = cols(v1, tgt2p); y2 = cols(cb, tgt2c)        # ens2 (param (i): pass=col_v1)

    def S(a, b, x, y):
        mix = a*x + b*y
        return float((mix.min(axis=1) - a*x.min(axis=1) - b*y.min(axis=1)).mean())

    def D(x, y):
        # E[x at argmin y] - E[min x]  (slope of S at a=0+)
        jy = y.argmin(axis=1)
        xj = x[np.arange(len(x)), jy]
        return float((xj - x.min(axis=1)).mean())

    s0 = S(t, w, x0, y0)
    s2_direct = S(t, lam*w, x2, y2)                    # = s2 (class-2 slack)
    # homogeneity: s2 = lam * S(t/lam, w) on ens2
    s2_hom = lam * S(t/lam, w, x2, y2)
    # dissection at common weights (t, w):
    S2_tw = S(t, w, x2, y2)
    weight_factor = (lam * S(t/lam, w, x2, y2)) / (lam * S2_tw)   # = S(t/lam,w)/S(t,w)
    ens_factor = S2_tw / s0
    # linear-regime slopes
    D2 = D(x2, y2); D0 = D(x0, y0)
    lin2 = (t/lam)*D2 / S(t/lam, w, x2, y2)            # linear approx quality (should be >=1, ->1 in regime)
    return dict(t=t, w=w, s0=s0, s2=s2_direct, s2_hom=s2_hom, S2_tw=S2_tw,
                wf=weight_factor, ef=ens_factor, D2=D2, D0=D0, lin2=lin2)

ITERS = {6: 2000, 8: 1500, 10: 1000, 12: 500, 14: 200}
for lam in [1.05, 1.10, 1.30, 1.70, 2.00]:
    for k in [6, 8, 10, 12, 14]:
        r = run(lam, k, ITERS[k])
        need = r['s2']/r['s0']/lam    # must be < 1
        print(f"lam={lam:.2f} k={k:2d}: s2/(lam*s0)={need:.4f} | weight_f={r['wf']:.4f} "
              f"ens_f={r['ef']:.4f} prod={r['wf']*r['ef']:.4f} | lin2={r['lin2']:.4f} "
              f"D2/D0={r['D2']/r['D0']:.4f} hom_chk={r['s2_hom']/r['s2']:.6f}", flush=True)
print("DONE")
