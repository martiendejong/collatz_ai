# 374: Route A step 10 — is the frozen-Jacobian spectral gap UNIFORM in k?
# If |lam2(J)|/rho stays flat below 1 across depths, the linear spectral-gap
# statement (the Psi-limit existence input) is empirically k-uniform.
import numpy as np
import os
from math import log2

ALPHA = log2(3.0)
CACHE = r"E:\projects\collatz\research\cache"

def rate(lam, k):
    N = 3**(k-1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_, r_ = np.divmod(i, 3)
    Nl = N//3
    m0, m2 = (r_ == 0), (r_ == 2)
    R1 = (4*s_) % Nl; R3 = (2*s_+1) % Nl
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    fn = os.path.join(CACHE, f"v_lam{lam:.2f}_k{k}.npy")
    v = np.load(fn).astype(np.float64) if os.path.exists(fn) else np.ones(N)
    rho = 1.0
    for _ in range(300 if os.path.exists(fn) else 1200):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        w[m2] += B3*cb[R3[m2]]
        w[m0] += B1*cb[R1[m0]]
        rho = float(w.max()); v = w/rho
    stack = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
    pick = stack.argmin(axis=0)*Nl + np.arange(Nl)
    rng = np.random.default_rng(374)
    u = rng.standard_normal(N)
    vn2 = float(v @ v)
    rs = []
    for it in range(260):
        u = u - (float(u @ v)/vn2)*v
        w = A*u[T4]
        ucb = u[pick]
        w[m2] += B3*ucb[R3[m2]]
        w[m0] += B1*ucb[R1[m0]]
        w /= rho
        g = float(np.linalg.norm(w)/np.linalg.norm(u))
        u = w
        if it > 210:
            rs.append(g)
    return float(np.mean(rs))

for lam in [1.70, 2.00]:
    out = []
    for k in [10, 12, 13, 14, 15, 16]:
        try:
            out.append((k, rate(lam, k)))
        except MemoryError:
            break
        print(f"lam={lam} k={k}: |lam2|/rho = {out[-1][1]:.4f}", flush=True)
    vals = [r for _, r in out]
    print(f"  reeks: {['%.4f' % x for x in vals]}  spreiding {max(vals)-min(vals):.4f}\n")
