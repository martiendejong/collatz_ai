# 372: Route A step 8 — does the finite Psi table CONVERGE as k grows?
# The table (keys (cls_n, cls_m, d mod 27)) dictates all second moments, hence
# the sd-defect and beta. If entries are Cauchy in k with geometric rate, the
# table has a limit and the uniformity of (a_*, beta_*) reduces to that limit.
# Measure: entry-wise |Psi_k - Psi_{k+1}|, per-entry and max, and the rate.
import numpy as np
import os
from math import log2
from collections import defaultdict

ALPHA = log2(3.0)
CACHE = r"E:\projects\collatz\research\cache"

def psi_table(lam, k, MOD=27):
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
    it = 400 if os.path.exists(fn) else 1500
    rho = 1.0
    for _ in range(it):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]
        w[m2] += B3*cb[R3[m2]]
        w[m0] += B1*cb[R1[m0]]
        rho = float(w.max()); v = w/rho
    t = A/rho
    w = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    wc = w - w.mean(); varw = float((wc**2).mean())
    nmax = min(int(np.ceil(np.log(1e-10)/np.log(t))), 60)
    idx = i.copy(); coeffs = []; args = []; clss = []
    for n in range(nmax):
        cls = idx % 3
        ar = np.where(cls == 0, R1[idx], np.where(cls == 2, R3[idx], 0)).astype(np.int64)
        coeffs.append(np.where(cls == 0, B1, np.where(cls == 2, B3, 0.0)))
        args.append(ar); clss.append(cls)
        idx = T4[idx]
    acc = defaultdict(list)
    for n in range(0, nmax, 2):
        for m in range(n, nmax):
            d = m - n
            for ca in [0, 2]:
                for cb_ in [0, 2]:
                    s2 = (coeffs[n] > 0) & (coeffs[m] > 0) & (clss[n] == ca) & (clss[m] == cb_)
                    if s2.sum() < 80:
                        continue
                    c = float((wc[args[n][s2]]*wc[args[m][s2]]).mean())/varw
                    acc[(ca, cb_, d % MOD)].append(c)
    return {key: float(np.mean(vals)) for key, vals in acc.items()}

for lam in [1.70, 2.00]:
    tabs = {}
    for k in [10, 11, 12, 13, 14]:
        tabs[k] = psi_table(lam, k)
    ks = sorted(tabs)
    common = set(tabs[ks[0]])
    for k in ks[1:]:
        common &= set(tabs[k])
    common = sorted(common)
    print(f"lam={lam}: {len(common)} gemeenschappelijke tabel-entries")
    deltas = []
    for j in range(len(ks)-1):
        d = max(abs(tabs[ks[j+1]][key] - tabs[ks[j]][key]) for key in common)
        rms = np.sqrt(np.mean([(tabs[ks[j+1]][key] - tabs[ks[j]][key])**2 for key in common]))
        deltas.append((ks[j], ks[j+1], d, rms))
        print(f"  k={ks[j]}->{ks[j+1]}: max|dPsi| = {d:.5f}  rms = {rms:.5f}")
    rates = [deltas[j+1][3]/deltas[j][3] for j in range(len(deltas)-1)]
    print(f"  rms-verval per stap: {['%.3f' % r for r in rates]}")
    # largest entries for scale reference
    big = sorted(common, key=lambda key: -abs(tabs[14][key]))[:4]
    print(f"  grootste entries (k=14): " +
          ", ".join(f"{key}:{tabs[14][key]:+.4f}" for key in big), flush=True)
