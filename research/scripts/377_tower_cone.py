# 377: Route A step 13 — Birkhoff diameter of the ONE-LEVEL tower map.
# Exact transport (Obs 552): the level field solves g = (A/rho) g(sigma4 .) + f_h,
# f_h(s) = [s=0 mod 3] (B1/rho) h(R1'(s)) + [s=2 mod 3] (B3/rho) h(R3'(s)).
# L: positive input h (on Nll) -> positive output g (on Nl); coefficients are
# k-FREE. Measure the projective diameter of L(cone) by sampling: if finite and
# flat in k, the per-level Birkhoff factor is k-free by construction and the
# refined uniformity route works.
import numpy as np
import os
from math import log2

ALPHA = log2(3.0)
CACHE = r"E:\projects\collatz\research\cache"

for lam in [1.70, 2.00]:
    for k in [10, 12, 14]:
        N = 3**(k-1); Nl = N//3; Nll = Nl//3
        s = np.arange(Nl, dtype=np.int64)
        sig4 = (4*s+2) % Nl
        sl = s // 3
        R1p = (4*sl) % Nll; R3p = (2*sl+1) % Nll
        r = s % 3
        A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
        v = np.load(os.path.join(CACHE, f"v_lam{lam:.2f}_k{k}.npy")).astype(np.float64)
        rho = 1.0
        for _ in range(200):
            cbv = np.minimum(np.minimum(v[:N//3], v[N//3:2*(N//3)], ), v[2*(N//3):])
            i = np.arange(N); T4 = (4*i+2) % N
            w = A*v[T4]
            ri = i % 3; si = i//3
            w[ri == 2] += B3*cbv[((2*si+1) % (N//3))[ri == 2]]
            w[ri == 0] += B1*cbv[((4*si) % (N//3))[ri == 0]]
            rho = float(w.max()); v = w/rho
        t = A/rho
        def L(h):
            f = np.zeros(Nl)
            f[r == 0] = (B1/rho)*h[R1p[r == 0]]
            f[r == 2] = (B3/rho)*h[R3p[r == 2]]
            g = np.zeros(Nl)
            for _ in range(int(np.log(1e-14)/np.log(t)) + 1):
                g = t*g[sig4] + f
            return g
        rng = np.random.default_rng(377)
        outs = [L(np.exp(1.0*rng.standard_normal(Nll))) for _ in range(6)]
        dmax = 0.0
        for a_ in range(6):
            for b_ in range(a_+1, 6):
                q = outs[a_]/outs[b_]
                dmax = max(dmax, float(np.log(q.max()) - np.log(q.min())))
        print(f"lam={lam} k={k:2d}: diameter L(kegel) proxy = {dmax:.3f}  "
              f"-> per-niveau factor tanh(D/4) = {np.tanh(dmax/4):.4f}", flush=True)
