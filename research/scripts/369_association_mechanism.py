# 369: Route A step 5 — the MECHANISM of the association rise (0.977 -> 0.996).
# Claim to test: by the gap-copy theorem applied at tower level 2, both input
# triples share their (1, t) part EXACTLY:
#   C  = C0  * (1, t, Rt)   with Rt = C[2]/C[0]
#   GC = GC0 * (1, t, xt)   with xt = GC[2]/GC[0]
# so the centred-triple correlation is a CLOSED FORM a = F(t, Rt, xt), and the
# entire misalignment is the single scalar m2 = xt - Rt: the level-2 margin,
# the same object as the certified (3b) criterion one tower level down.
# Predictions: (P1) C[1] = t*C[0] and GC[1] = t*GC[0] exactly;
# (P2) measured a equals F(t, Rt, xt) exactly; (P3) 1 - a ~ G(t,Rt) * m2^2
# (quadratic in the margin); (P4) the decay rate of 1-a per depth equals the
# SQUARE of the margin-decay rate (margin law, Obs 321: margin ~ sqrt(c)^k).
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
    return v, rho, A

def corr3(u, w):
    uc = u - u.mean(); wc = w - w.mean()
    return float(uc @ wc/np.sqrt((uc@uc)*(wc@wc)))

rows = {}
print(f"{'lam':>5} {'k':>3} {'copyC err':>9} {'copyGC err':>10} {'Rt':>8} {'xt':>8} "
      f"{'m2=xt-Rt':>9} {'a meten':>8} {'a=F(...)':>9} {'1-a':>9} {'Gm2^2/(1-a)':>11}")
for lam in [1.05, 1.30, 1.70, 2.00]:
    for k in [10, 12, 13, 14]:
        v, rho, A = load_or_make(lam, k)
        N = v.size; Nl = N//3; Nll = Nl//3
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        cb3 = (cb[:Nll] + cb[Nll:2*Nll] + cb[2*Nll:])/3.0
        cbb = np.minimum(np.minimum(cb[:Nll], cb[Nll:2*Nll]), cb[2*Nll:])
        Gc = cb3 - cbb
        x = np.arange(Nll)
        C = np.array([cb3[x % 3 == b].mean() for b in range(3)])
        GC = np.array([Gc[x % 3 == b].mean() for b in range(3)])
        t = A/rho
        # P1: copy identities
        ec = abs(C[1] - t*C[0])/C[0]
        eg = abs(GC[1] - t*GC[0])/GC[0]
        Rt = C[2]/C[0]; xt = GC[2]/GC[0]
        m2 = xt - Rt
        a_meas = corr3(C, GC)
        a_form = corr3(np.array([1, t, Rt]), np.array([1, t, xt]))
        # quadratic prefactor: 1 - F(t,Rt,Rt+m) ~ G*m^2 for small m: estimate G numerically
        eps = 1e-6
        G_pref = (1 - corr3(np.array([1, t, Rt]), np.array([1, t, Rt+eps])))/eps**2
        pred = G_pref*m2*m2
        rows.setdefault(lam, []).append((k, m2, 1-a_meas))
        print(f"{lam:>5} {k:>3} {ec:>9.1e} {eg:>10.1e} {Rt:>8.4f} {xt:>8.4f} "
              f"{m2:>9.5f} {a_meas:>8.5f} {a_form:>9.5f} {1-a_meas:>9.2e} "
              f"{pred/(1-a_meas):>11.4f}", flush=True)
print()
print("decay rates per depth step (P4):")
for lam, rr in rows.items():
    rr.sort()
    ks = [r[0] for r in rr]
    m2s = [r[1] for r in rr]
    das = [r[2] for r in rr]
    rm = [(m2s[i+1]/m2s[i])**(1/(ks[i+1]-ks[i])) for i in range(len(rr)-1)]
    ra = [(das[i+1]/das[i])**(1/(ks[i+1]-ks[i])) for i in range(len(rr)-1)]
    print(f"  lam={lam}: margin-rate {['%.3f' % r for r in rm]}  (1-a)-rate {['%.3f' % r for r in ra]}  "
          f"rate(1-a)/rate(m2)^2 {['%.3f' % (ra[i]/rm[i]**2) for i in range(len(rm))]}")
