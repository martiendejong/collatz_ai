# 371: Route A step 7 — taming the slack form factor via the kappa-linearised
# slack functional. Pointwise on parent-0 fibers (s = 0 mod 3):
#   rho*cb(s) = min_j [X_j + Y_j],  X_j = A*v(sigma4 s + j*Nl),
#                                   Y_j = B1*cb(R1'(s) + j*Nll)
#   slack(s) = rho*cb(s) - A*cb(sigma4 s) - B1*cbb(R1'(s))  (>= 0 exact)
# kappa model (Obs 533: E[min] = mean - kappa*sd):
#   slack(s) ~ kappa * D(s),  D(s) = sd(X) + sd(Y) - sd(X+Y)  (>= 0)
# D is a SECOND-MOMENT object: fully determined by the covariance structure,
# i.e. by the finite Psi table. If slack ~ kappa*D with stable kappa, then
# Lemma beta reduces to Psi-table facts. Tests:
#  (K1) pointwise correlation corr(slack, D) per fiber;
#  (K2) kappa = E[slack]/E[D]: stability across cells b, lambda, k;
#  (K3) shape reconstruction: does kappa*E[D|cell b] reproduce the stable
#       S1-triple shape [1, 1.1, ~3]?
import numpy as np
import os
from math import log2

ALPHA = log2(3.0)
CACHE = r"E:\projects\collatz\research\cache"

def load_or_make(lam, k, iters=1400):
    N = 3**(k-1)
    fn = os.path.join(CACHE, f"v_lam{lam:.2f}_k{k}.npy")
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_, r_ = np.divmod(i, 3)
    Nl = N//3
    m0, m2 = (r_ == 0), (r_ == 2)
    R1 = (4*s_) % Nl; R3 = (2*s_+1) % Nl
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    v = np.load(fn).astype(np.float64) if os.path.exists(fn) else np.ones(N)
    it = 400 if os.path.exists(fn) else iters
    rho = 1.0
    for _ in range(it):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        w[m2] += B3*cb[R3[m2]]
        w[m0] += B1*cb[R1[m0]]
        rho = float(w.max()); v = w/rho
    return v, rho, A, B1, B3

def sd3(a, b, c):
    m = (a+b+c)/3.0
    return np.sqrt(((a-m)**2 + (b-m)**2 + (c-m)**2)/3.0)

print(f"{'lam':>5} {'k':>3} {'corr(sl,D)':>10} {'kappa_glob':>10} {'kappa per b':>24} "
      f"{'S-shape':>20} {'kD-shape':>20}")
for lam in [1.05, 1.30, 1.70, 2.00]:
    for k in [12, 14]:
        v, rho, A, B1, B3 = load_or_make(lam, k)
        N = v.size; Nl = N//3; Nll = Nl//3
        s = np.arange(Nl, dtype=np.int64)
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        cbb = np.minimum(np.minimum(cb[:Nll], cb[Nll:2*Nll]), cb[2*Nll:])
        # parent-0 fibers
        f0 = s[s % 3 == 0]
        sig4 = (4*f0+2) % Nl
        r1 = (4*(f0//3)) % Nll
        X = [A*v[(sig4 + j*Nl) % N] for j in range(3)]
        Y = [B1*cb[(r1 + j*Nll) % Nl] for j in range(3)]
        # pointwise slack (exact identity)
        lhs = rho*cb[f0]
        slack = lhs - A*cb[sig4] - B1*cbb[r1]
        D = sd3(*X) + sd3(*Y) - sd3(X[0]+Y[0], X[1]+Y[1], X[2]+Y[2])
        m = (slack > 1e-300) | (D > 1e-300)
        corr = float(np.corrcoef(slack[m], D[m])[0, 1])
        kap_g = float(slack.mean()/D.mean())
        b_idx = (f0//3) % 3
        kaps = [float(slack[b_idx == b].mean()/D[b_idx == b].mean()) for b in range(3)]
        Sc = np.array([slack[b_idx == b].mean() for b in range(3)])
        Dc = np.array([D[b_idx == b].mean() for b in range(3)])
        print(f"{lam:>5} {k:>3} {corr:>10.4f} {kap_g:>10.4f} "
              f"{np.array2string(np.array(kaps), precision=3):>24} "
              f"{np.array2string(Sc/Sc[0], precision=3):>20} "
              f"{np.array2string(kap_g*Dc/Sc[0], precision=3):>20}", flush=True)
