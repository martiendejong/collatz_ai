"""
207_common_cause.py
===================
The common-cause test for the rich-rough coupling (Obs 410 -> proof
skeleton). Hypothesis: level and roughness are associated because both
are increasing functions of the SAME type variable -- the D3-share
along the ancestry -- whose trits are exactly i.i.d. (Freshness Lemma,
proved). D3 branches carry the largest weights (B3 > 1 > B1: they make
blocks rich) AND inject the most roughness (C_inj lifts: D3 = 2.94 vs
D1 = 1.28). If true, EPW closes the association:
shared cause + independent trits + monotonicity => Cov >= 0,
quantitatively via the type variance.

Test: per base triple (as in 206) compute (logL, logS) plus the
D3-share along the first G steps of the selected-feed type walk from
the triple's bases. Then:
  (i)  monotonicity: mean logL and mean logS per D3-share bucket;
  (ii) the common-cause kill: partial corr(logL, logS | D3-share) --
       residual correlation after regressing both on the share. A drop
       from 0.79 toward ~0 confirms the skeleton; a large residual
       means other shared structure carries part of the coupling.
Systems 15..17, lam = 1.70, G = 8 walk steps.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM = 1.70
A, B1, B3 = LAM ** -2.0, LAM ** (ALPHA - 2.0), LAM ** (ALPHA - 1.0)
G = 8


def perron_and_walk(k, n_iter=300):
    N = 3 ** (k - 1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    m1, m3 = (r == 0), (r == 2)
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        w2 = A * v[T4]
        w2[m1] += B1 * cb[R1[m1]]
        w2[m3] += B3 * cb[R3[m3]]
        v = w2 / w2.max()
    v = v / v.mean()
    # D3-share along the type walk, G steps. The type walk is
    # selection-independent (Type-rigidity Lemma), and each feed step
    # DESCENDS one level: the modulus shrinks by 3 per step (Freshness).
    pos = np.arange(Nl, dtype=np.int64)
    d3count = np.zeros(Nl, dtype=np.int64)
    alive = np.ones(Nl, dtype=bool)
    mod = Nl
    for _ in range(G):
        t = pos % 3
        d3count[(t == 2) & alive] += 1
        alive = alive & (t != 1)                  # D2 absorbs the walk
        sp = pos // 3
        mod //= 3
        nxt = np.where(t == 0, (4 * sp) % mod, (2 * sp + 1) % mod)
        pos = np.where(alive, nxt, pos)
    return v, d3count


print(f"common-cause test at lam = {LAM}  (G = {G} walk steps)", flush=True)
print("  sys  corr(L,S)  partial|share  mono(L)  mono(S)  "
      "corr(L,share) corr(S,share)", flush=True)
for kp1 in (15, 16, 17):
    v, d3 = perron_and_walk(kp1)
    Nl = v.size // 3
    mn = (v[:Nl] + v[Nl:2 * Nl] + v[2 * Nl:]) / 3.0
    Nl3 = Nl // 3
    M = np.stack([mn[:Nl3], mn[Nl3:2 * Nl3], mn[2 * Nl3:]])
    level = M.mean(axis=0)
    spread = M.std(axis=0)
    # triple-level D3-share: mean d3count over the triple's bases
    D = np.stack([d3[:Nl3], d3[Nl3:2 * Nl3], d3[2 * Nl3:]]).mean(axis=0)
    ok = spread > 0
    x = np.log(level[ok]); x -= x.mean()
    y = np.log(spread[ok]); y -= y.mean()
    z = D[ok].astype(np.float64); z -= z.mean()
    def c(a, b):
        return float(np.mean(a * b) / np.sqrt(np.mean(a * a) * np.mean(b * b)))
    cxy, cxz, cyz = c(x, y), c(x, z), c(y, z)
    partial = (cxy - cxz * cyz) / np.sqrt((1 - cxz ** 2) * (1 - cyz ** 2))
    # monotonicity: mean logL / logS per integer share bucket
    zz = D[ok]
    buckets = np.round(zz * 2) / 2
    ub = np.unique(buckets)
    mL = [float(np.mean(np.log(level[ok])[buckets == b])) for b in ub]
    mS = [float(np.mean(np.log(spread[ok])[buckets == b])) for b in ub]
    monoL = bool(np.all(np.diff(mL) > 0))
    monoS = bool(np.all(np.diff(mS) > 0))
    print(f"  {kp1:2d}  {cxy:.5f}    {partial:.5f}       {monoL}    "
          f"{monoS}    {cxz:.5f}      {cyz:.5f}", flush=True)
print("done", flush=True)
