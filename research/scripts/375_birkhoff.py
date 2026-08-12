# 375: Route A step 11 — Birkhoff contraction of the frozen Jacobian, measured
# in Hilbert's projective metric: d_H(x,y) = log max(x/y) + log max(y/x).
# Track d_H(J^m u, J^m w)/d_H-previous for positive starts: per-step contraction
# factor; uniformity across k => the spectral gap is Birkhoff-provable.
import numpy as np
import os
from math import log2

ALPHA = log2(3.0)
CACHE = r"E:\projects\collatz\research\cache"

def dH(x, y):
    r = x/y
    return float(np.log(r.max()) - np.log(r.min()))

for lam in [1.70, 2.00]:
    for k in [10, 12, 14]:
        N = 3**(k-1)
        i = np.arange(N, dtype=np.int64)
        T4 = (4*i+2) % N
        s_, r_ = np.divmod(i, 3)
        Nl = N//3
        m0, m2 = (r_ == 0), (r_ == 2)
        R1 = (4*s_) % Nl; R3 = (2*s_+1) % Nl
        A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
        v = np.load(os.path.join(CACHE, f"v_lam{lam:.2f}_k{k}.npy")).astype(np.float64)
        for _ in range(200):
            cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
            w = A*v[T4]
            w[m2] += B3*cb[R3[m2]]
            w[m0] += B1*cb[R1[m0]]
            v = w/w.max()
        pick = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]]).argmin(axis=0)*Nl + np.arange(Nl)
        def J(u):
            w = A*u[T4]
            uc = u[pick]
            w[m2] += B3*uc[R3[m2]]
            w[m0] += B1*uc[R1[m0]]
            return w/w.max()
        rng = np.random.default_rng(375)
        x = np.exp(0.5*rng.standard_normal(N))*v
        y = np.exp(0.5*rng.standard_normal(N))*v
        ds = [dH(x, y)]
        M = 6*k
        for m in range(M):
            x = J(x); y = J(y)
            ds.append(dH(x, y))
        # per-step factor in the settled regime (last third)
        tail = [ds[j+1]/ds[j] for j in range(2*M//3, M) if ds[j] > 1e-12]
        f = float(np.mean(tail)) if tail else float('nan')
        print(f"lam={lam} k={k:2d}: d_H {ds[0]:.2f} -> {ds[-1]:.2e}; "
              f"contractie/stap (staart) = {f:.4f}", flush=True)
