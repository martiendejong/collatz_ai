# 380: Route A step 17 — evaluate the Doeblin quadratic form
# Q = sum w~_n w~_n' phi(n-n') explicitly: exact k-free weights (t^n, class
# pattern), phi measured from the log-fluctuations of cb3 along the sigma4
# orbit (the Psi machinery). Check sqrt(Q) vs the directly measured 0.60,
# stability k=12 vs 14, and a Cauchy-tail-inflated upper bound (< 1?).
import numpy as np
import os
from math import log2

ALPHA = log2(3.0)
CACHE = r"E:\projects\collatz\research\cache"

for k in [12, 14]:
    lam = 2.00
    N = 3**(k-1); Nl = N//3; Nll = Nl//3
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    v = np.load(os.path.join(CACHE, f"v_lam{lam:.2f}_k{k}.npy")).astype(np.float64)
    rho = 1.0
    i = np.arange(N); T4 = (4*i+2) % N
    ri = i % 3; si = i//3
    for _ in range(200):
        cbv = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        w[ri == 2] += B3*cbv[((2*si+1) % Nl)[ri == 2]]
        w[ri == 0] += B1*cbv[((4*si) % Nl)[ri == 0]]
        rho = float(w.max()); v = w/rho
    t = A/rho
    cbv = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    cb3 = (cbv[:Nll] + cbv[Nll:2*Nll] + cbv[2*Nll:])/3.0
    # weights along the backbone from a class-0 start: classes cycle 0->2->1
    nmax = min(int(np.log(1e-10)/np.log(t)), 60)
    cls_seq = [(0 + 2*n) % 3 for n in range(nmax)]   # class of sigma4^n s for s=0 mod 3
    wts = np.array([t**n*(B1 if c == 0 else (B3 if c == 2 else 0.0))
                    for n, c in enumerate(cls_seq)])
    wts = wts/wts.sum()
    # phi(d): correlation of centred log cb3 at feed-args separated by lag d
    lx = np.log(cb3)
    s0 = np.arange(0, Nl, 3)
    x = s0.copy()
    args = []
    for n in range(nmax):
        c = x % 3
        a_ = np.where(c == 0, (4*(x//3)) % Nll, np.where(c == 2, (2*(x//3)+1) % Nll, 0))
        args.append((c, a_))
        x = (4*x+2) % Nl
    fl = lx - lx.mean()
    var = float((fl**2).mean())
    phi = {}
    for d in range(nmax):
        vals = []
        for n in range(nmax - d):
            c1, a1 = args[n]; c2, a2 = args[n+d]
            m = (c1 != 1) & (c2 != 1)
            if m.sum() > 100:
                vals.append(float((fl[a1[m]]*fl[a2[m]]).mean())/var)
        phi[d] = float(np.mean(vals)) if vals else 0.0
    Q = 0.0
    for n in range(nmax):
        for m_ in range(nmax):
            if wts[n] > 0 and wts[m_] > 0:
                Q += wts[n]*wts[m_]*phi[abs(n-m_)]
    tail = 0.04/(1-0.85)   # Cauchy tail inflation on phi entries (Obs 560)
    Qhi = Q + tail*1.0     # crude: |dphi| <= tail on all entries, sum w=1
    print(f"lam={lam} k={k}: sqrt(Q) = {np.sqrt(Q):.4f}  (direct gemeten 0.60-0.62)  "
          f"Q = {Q:.4f}; met Cauchy-staartopslag Q_hi = {Qhi:.4f} "
          f"({'<' if Qhi < 1 else '>='} 1)", flush=True)
