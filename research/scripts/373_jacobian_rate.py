# 373: Route A step 9 — is the Psi-table Cauchy rate (~0.85) the subleading
# eigenvalue of the frozen Jacobian? Mechanism candidate: the constant-lift
# embedding (Lemma S1) maps depth k into k+1; the embedded vector relaxes to
# the true one under the deeper operator, governed by |lam2(J)|/rho. Then
# cross-depth table differences decay at that spectral rate.
import numpy as np
import os
from math import log2

ALPHA = log2(3.0)
CACHE = r"E:\projects\collatz\research\cache"

for lam in [1.70, 2.00]:
    k = 13
    N = 3**(k-1)
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_, r_ = np.divmod(i, 3)
    Nl = N//3
    m0, m2 = (r_ == 0), (r_ == 2)
    R1 = (4*s_) % Nl; R3 = (2*s_+1) % Nl
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    v = np.load(os.path.join(CACHE, f"v_lam{lam:.2f}_k{k}.npy")).astype(np.float64)
    rho = 1.0
    for _ in range(400):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        w[m2] += B3*cb[R3[m2]]
        w[m0] += B1*cb[R1[m0]]
        rho = float(w.max()); v = w/rho
    # frozen argmin per cell
    stack = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
    arg = stack.argmin(axis=0)
    pick = arg*Nl + np.arange(Nl)      # index of the winning lift per cell
    def J(u):
        ucb = u[pick]
        w = A*u[T4]
        w[m2] += B3*ucb[R3[m2]]
        w[m0] += B1*ucb[R1[m0]]
        return w
    # deflated power iteration: project out the Perron direction (weight v)
    rng = np.random.default_rng(373)
    u = rng.standard_normal(N)
    vn2 = float(v @ v)
    rates = []
    for it in range(300):
        u = u - (float(u @ v)/vn2)*v
        w = J(u)/rho
        g = float(np.linalg.norm(w)/np.linalg.norm(u))
        u = w
        if it > 250:
            rates.append(g)
    lam2 = float(np.mean(rates))
    print(f"lam={lam} k={k}: |lam2(J)|/rho = {lam2:.4f}   "
          f"(Psi-tabelrate ~0.85; kappa_deep 0.839; sqrt(c): "
          f"{np.sqrt(0.70 if lam == 1.70 else 0.835):.3f})", flush=True)
