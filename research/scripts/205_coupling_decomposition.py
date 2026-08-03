"""
205_coupling_decomposition.py
=============================
Decompose the selection-flatness coupling (Obs 408) into SCALE x SHAPE,
and run the desert test.

Per top-triple of the base space: the min field's deviation factors as
(local min/mean LEVEL ratio) x (mean-field deviation) x (shape
residual). Globally:

    f2 = scale x shape,
    scale = dev-weighted local q  = sqrt( sum (Cbar/Ubar)^2 dev_u^2
                                          / sum dev_u^2 )
    shape = f2 / scale            (the genuine flatness coupling beyond
                                   the trivial level effect)

Readings: if SCALE carries the creep, the f2-erosion is an offshoot of
the proven q-machinery (homogenization pushes the local min/mean level
ratio to 1) and the fork becomes "does shape stay below
threshold/lim-scale". If SHAPE carries it, the coupling itself erodes
-- one layer deeper still. Desert test: the same numbers per quartile
of the triple level Ubar (low quartile ~ desert-suppressed bases);
prediction (registered): the flat-selection effect concentrates in low
quartiles if desert suppression is the mechanism making lifts low AND
flat.  Systems 14..18, lam = 1.70.
"""
import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM = 1.70
A, B1, B3 = LAM ** -2.0, LAM ** (ALPHA - 2.0), LAM ** (ALPHA - 1.0)


def perron(k, n_iter=300):
    N = 3 ** (k - 1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    R1m = ((4 * s) % Nl)[r == 0]
    R3m = ((2 * s + 1) % Nl)[r == 2]
    m1, m3 = (r == 0), (r == 2)
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2 * Nl]), v[2 * Nl:])
        w2 = A * v[T4]
        w2[m1] += B1 * cb[R1m]
        w2[m3] += B3 * cb[R3m]
        v = w2 / w2.max()
    return v / v.mean()


print(f"coupling decomposition f2 = scale x shape at lam = {LAM}",
      flush=True)
print("  sys   f2       scale    shape    | per level-quartile "
      "(low->high): f2 / scale / shape", flush=True)
for kp1 in (14, 15, 16, 17, 18):
    v = perron(kp1)
    Nl = v.size // 3
    L0, L1, L2 = v[:Nl], v[Nl:2 * Nl], v[2 * Nl:]
    cb = np.minimum(np.minimum(L0, L1), L2)
    mn = (L0 + L1 + L2) / 3.0

    Nl3 = Nl // 3
    def trip(field):
        return np.stack([field[:Nl3], field[Nl3:2 * Nl3], field[2 * Nl3:]])
    Cb = trip(cb)
    Mn = trip(mn)
    Cbar = Cb.mean(axis=0)
    Ubar = Mn.mean(axis=0)
    dev_c = Cb - Cbar[None, :]
    dev_u = Mn - Ubar[None, :]
    w_u2 = (dev_u ** 2).sum(axis=0)          # per-triple mean-dev energy
    w_c2 = (dev_c ** 2).sum(axis=0)
    ratio = Cbar / Ubar                       # local q per triple

    def stats(mask):
        su = float(w_u2[mask].sum())
        sc = float(w_c2[mask].sum())
        f2 = (sc / su) ** 0.5
        scale = float(((ratio[mask] ** 2 * w_u2[mask]).sum() / su)) ** 0.5
        return f2, scale, f2 / scale

    allm = np.ones(Nl3, dtype=bool)
    f2a, sca, sha = stats(allm)
    qs = np.quantile(Ubar, [0.25, 0.5, 0.75])
    parts = []
    lo = -np.inf
    for hi in list(qs) + [np.inf]:
        m = (Ubar > lo) & (Ubar <= hi)
        parts.append(stats(m))
        lo = hi
    qtxt = "  ".join(f"{f:.3f}/{s:.3f}/{h:.3f}" for f, s, h in parts)
    print(f"  {kp1:2d}   {f2a:.5f}  {sca:.5f}  {sha:.5f}  |  {qtxt}",
          flush=True)
print("done", flush=True)
