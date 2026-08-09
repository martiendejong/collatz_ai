# 359: G2 from the 2026-08-09 analysis — refine the Psi structure-function
# modulus at lam=1.05 (Obs 532 left closure at 0.9727 with keys (cls,cls,d mod 27)
# and a TRUNCATED training window m <= n+30). Two changes:
#   (1) train on the FULL lag range (all m in [n, nmax), even n only — still
#       out-of-sample for odd n),
#   (2) sweep the modulus: d mod 9 / 27 / 81.
# Question: does Var_pred/Var_true close 0.9727 -> 1.000 at lam=1.05, and is the
# remaining defect a modulus issue (mod-81 fixes it) or a window issue (full
# training at mod-27 already fixes it)?
import numpy as np
from math import log2
from collections import defaultdict

ALPHA = log2(3.0)

def run(lam, k, MOD, full_window):
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    v = np.ones(N); rho = 1.0
    for _ in range(1500):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w_ = A*v[T4]; m2 = r_arr==2; m0 = r_arr==0
        w_[m2] += B3*cb[R3[m2]]; w_[m0] += B1*cb[R1[m0]]
        rho = float(w_.max()); w_ /= rho; v = w_
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
        mhi = nmax if full_window else min(n+30, nmax)
        for m in range(n, mhi):
            d = m - n
            for ca in [0, 2]:
                for cb_ in [0, 2]:
                    s2 = (coeffs[n] > 0) & (coeffs[m] > 0) & (clss[n] == ca) & (clss[m] == cb_)
                    if s2.sum() < 80: continue
                    c = float((wc[args[n][s2]]*wc[args[m][s2]]).mean())/varw
                    acc[(ca, cb_, d % MOD)].append(c)
    psi = {key: np.mean(vals) for key, vals in acc.items()}

    var_pred = 0.0; miss = 0
    for n in range(nmax):
        for m in range(nmax):
            ca_arr = clss[n]; cb_arr = clss[m]
            for ca in [0, 2]:
                for cb_ in [0, 2]:
                    kk = (ca, cb_, (m - n) % MOD) if m >= n else (cb_, ca, (n - m) % MOD)
                    if kk not in psi:
                        miss += 1; continue
                    s2 = (ca_arr == ca) & (cb_arr == cb_)
                    Ba = B1 if ca == 0 else B3
                    Bb = B1 if cb_ == 0 else B3
                    var_pred += (t**n)*(t**m)*Ba*Bb*float(s2.mean())*psi[kk]*varw/(rho**2)
    v_fluct = (1.0/rho)*sum((t**n)*coeffs[n]*wc[args[n]] for n in range(nmax))
    var_true = float((v_fluct**2).mean())
    return var_pred/var_true, len(psi), miss

for lam in [1.05, 1.70]:
    for k in [10]:
        base, nb, mb = run(lam, k, 27, full_window=False)
        print(f"lam={lam} k={k}  mod27 window30 (Obs 532 baseline): {base:.4f}  (klassen {nb}, gemist {mb})")
        for MOD in [9, 27, 81]:
            r, nc, ms = run(lam, k, MOD, full_window=True)
            print(f"lam={lam} k={k}  mod{MOD:<2} volledig venster:        {r:.4f}  (klassen {nc}, gemist {ms})")
    print()
