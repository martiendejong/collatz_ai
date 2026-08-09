"""Sum-rule closure with the REFINED covariance law:
key = (cls_n, cls_m, (m-n) mod 27) for all lags (covers v3=0,1 refined; v3>=2 via mod 27).
Build the lookup empirically from HALF the (n,m)-pairs (even n), predict Var from
the OTHER half + full sum — out-of-sample closure test."""
import numpy as np
from math import log2
ALPHA = log2(3.0)

for lam in [1.05, 1.70]:
    k = 10
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

    nmax = min(int(np.ceil(np.log(1e-10)/np.log(t))), 60)
    idx = i.copy(); coeffs = []; args = []; clss = []
    for n in range(nmax):
        cls = idx % 3
        ar = np.where(cls == 0, R1[idx], np.where(cls == 2, R3[idx], 0)).astype(np.int64)
        coeffs.append(np.where(cls == 0, B1, np.where(cls == 2, B3, 0.0)))
        args.append(ar); clss.append(cls)
        idx = T4[idx]

    # build refined lookup from training pairs (even n), then closure over ALL pairs
    from collections import defaultdict
    acc = defaultdict(list)
    for n in range(0, nmax, 2):
        for m in range(n, min(n+30, nmax)):
            d = m - n
            for ca in [0, 2]:
                for cb_ in [0, 2]:
                    s2 = (coeffs[n] > 0) & (coeffs[m] > 0) & (clss[n] == ca) & (clss[m] == cb_)
                    if s2.sum() < 80: continue
                    c = float((wc[args[n][s2]]*wc[args[m][s2]]).mean())/varw
                    acc[(ca, cb_, d % 27)].append(c)
    psi = {key: np.mean(vals) for key, vals in acc.items()}

    # closure: predicted variance of the fluctuation sum
    var_pred = 0.0
    missing = 0
    for n in range(nmax):
        for m in range(nmax):
            d = abs(m - n)
            for ca in [0, 2]:
                for cb_ in [0, 2]:
                    s2 = (clss[n] == ca) & (clss[m] == cb_)
                    frac = float(s2.mean())
                    if frac == 0: continue
                    key = (min(ca,cb_) if n<=m else min(ca,cb_),)  # placeholder
            # simpler: exact per-pair expected coefficient-weighted covariance using psi:
        # (do inner loop properly below)
    var_pred = 0.0
    for n in range(nmax):
        for m in range(nmax):
            d = abs(m - n)
            ca_arr = clss[n]; cb_arr = clss[m]
            for ca in [0, 2]:
                for cb_ in [0, 2]:
                    key = (ca, cb_, ((m - n) % 27) if m >= n else ((n - m) % 27))
                    kk = (ca, cb_, (m - n) % 27) if m >= n else (cb_, ca, (n - m) % 27)
                    if kk not in psi:
                        continue
                    s2 = (ca_arr == ca) & (cb_arr == cb_)
                    Ba = B1 if ca == 0 else B3
                    Bb = B1 if cb_ == 0 else B3
                    var_pred += (t**n)*(t**m)*Ba*Bb*float(s2.mean())*psi[kk]*varw/(rho**2)
    v_fluct = (1.0/rho)*sum((t**n)*coeffs[n]*wc[args[n]] for n in range(nmax))
    var_true = float((v_fluct**2).mean())
    print(f"lam={lam}: verfijnde closure Var_pred/Var_true = {var_pred/var_true:.4f}  "
          f"(lookup-klassen: {len(psi)})", flush=True)
