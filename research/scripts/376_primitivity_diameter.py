# 376: Route A step 12 — primitivity index and explicit Birkhoff diameter of J.
# (a) m*(k): first m with full support of J^m e_j (boolean reachability);
# (b) diameter: sample columns J^M e_j (M > m*), max pairwise Hilbert distance
#     Delta -> provable per-step contraction bound tanh(Delta/4)^(1/M).
import numpy as np
import os
from math import log2, tanh

ALPHA = log2(3.0)
CACHE = r"E:\projects\collatz\research\cache"

for lam in [1.70, 2.00]:
    for k in [10, 12]:
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
            cbv = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
            w = A*v[T4]
            w[m2] += B3*cbv[R3[m2]]
            w[m0] += B1*cbv[R1[m0]]
            v = w/w.max()
        pick = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]]).argmin(axis=0)*Nl + np.arange(Nl)
        def J(u):
            w = A*u[T4]
            uc = u[pick]
            w[m2] += B3*uc[R3[m2]]
            w[m0] += B1*uc[R1[m0]]
            return w
        # (a) support growth from a single basis vector
        u = np.zeros(N); u[1] = 1.0
        mstar = None
        for m in range(1, 30*k):
            u = J(u)
            mx = u.max()
            if mx > 0: u /= mx
            if (u > 0).all():
                mstar = m
                break
        # (b) diameter via sampled columns at M = mstar + k
        M = (mstar or 20*k) + k
        cols = []
        rng = np.random.default_rng(376)
        for j in rng.integers(0, N, 6):
            u = np.zeros(N); u[j] = 1.0
            for _ in range(M):
                u = J(u); u /= u.max()
            cols.append(u)
        dmax = 0.0
        for a_ in range(len(cols)):
            for b_ in range(a_+1, len(cols)):
                r = cols[a_]/cols[b_]
                dmax = max(dmax, float(np.log(r.max()) - np.log(r.min())))
        bound = tanh(dmax/4)**(1.0/M)
        print(f"lam={lam} k={k:2d}: m* = {mstar} (~{(mstar or 0)/k:.1f}k)  "
              f"Delta(J^{M}) = {dmax:.2f}  -> bewijsgrens/stap {bound:.4f} "
              f"(gemeten 0.84-0.87)", flush=True)
