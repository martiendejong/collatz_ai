"""Refine the v3=0 class: group Cov(w o arg_n, w o arg_m) by
(class at step n, class at step m, (m-n) mod 9). Does the scatter collapse?"""
import numpy as np
from math import log2
ALPHA = log2(3.0)

lam = 1.05; k = 10
A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
N = 3**(k-1); Nl = N//3
i = np.arange(N, dtype=np.int64)
T4 = (4*i+2) % N
s_arr, r_arr = np.divmod(i, 3)
m0 = r_arr==0; m2 = r_arr==2
R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
v = np.ones(N); rho = 1.0
for _ in range(1500):
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w_ = A*v[T4]; w_[m2] += B3*cb[R3[m2]]; w_[m0] += B1*cb[R1[m0]]
    rho = float(w_.max()); w_ /= rho; v = w_
t = A/rho
w = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
wc = w - w.mean(); varw = float((wc**2).mean())

nmax = 55
idx = i.copy(); coeffs = []; args = []; clss = []
for n in range(nmax):
    cls = idx % 3
    ar = np.where(cls == 0, R1[idx], np.where(cls == 2, R3[idx], 0)).astype(np.int64)
    coeffs.append(np.where(cls == 0, B1, np.where(cls == 2, B3, 0.0)))
    args.append(ar); clss.append(cls)
    idx = T4[idx]

groups = {}
for n in range(nmax):
    for m in range(n+1, min(n+28, nmax)):
        d = m - n
        if d % 3 == 0: continue  # only v3=0
        sel = (coeffs[n] > 0) & (coeffs[m] > 0)
        if sel.sum() < 100: continue
        # endpoint "R-type" per row varies; group rows by (cls_n, cls_m) pairs
        for ca in [0, 2]:
            for cb_ in [0, 2]:
                s2 = sel & (clss[n] == ca) & (clss[m] == cb_)
                if s2.sum() < 100: continue
                c = float((wc[args[n][s2]]*wc[args[m][s2]]).mean())/varw
                groups.setdefault((ca, cb_, d % 9), []).append(c)

print("verfijning v3=0: (cls_n, cls_m, d mod 9) -> mean +- sd (n)")
tight = 0; total = 0
for key in sorted(groups):
    g = groups[key]
    total += 1
    if np.std(g) < 0.05: tight += 1
    print(f"  {key}: {np.mean(g):+.4f} +- {np.std(g):.4f} (n={len(g)})")
print(f"\nklassen met sd < 0.05: {tight}/{total}")
