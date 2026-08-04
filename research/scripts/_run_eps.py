import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70

def make_operator(lam, eps):
    A  = lam ** -2.0
    B3_full = lam ** (ALPHA - 1.0)
    B1_full = lam ** (ALPHA - 2.0)
    return A, eps * B1_full, eps * B3_full

def perron_and_stats(k, lam, eps, n_iter=250):
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    A, B1, B3 = make_operator(lam, eps)
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        if eps > 0:
            w[m2] += B3 * cb[R3[m2]]
            w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    v /= v.mean()
    cv = float(np.std(v) / np.mean(v))
    cb2 = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    w2  = A * v[T4]
    if eps > 0:
        w2[m2] += B3 * cb2[R3[m2]]
        w2[m0] += B1 * cb2[R1[m0]]
    rho = float(w2.sum() / v.sum())
    stack = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
    sel   = stack.argmin(axis=0).astype(np.int64)
    tgt2  = R3 + sel[R3] * Nl
    tgt0  = R1 + sel[R1] * Nl
    def PW(x):
        m = (x[:Nl] + x[Nl:2*Nl] + x[2*Nl:]) / 3.0
        y = x.copy(); y[:Nl] -= m; y[Nl:2*Nl] -= m; y[2*Nl:] -= m
        return y
    rng2 = np.random.default_rng(1)
    d    = PW(rng2.standard_normal(N)); d /= np.linalg.norm(d)
    rates = []
    for _ in range(180):
        y = A * d[T4]
        if eps > 0:
            y[m2] += B3 * d[tgt2[m2]]
            y[m0] += B1 * d[tgt0[m0]]
        y = PW(y); nrm = np.linalg.norm(y); rates.append(nrm); d = y/(nrm+1e-300)
    sw_rho = float(np.exp(np.mean(np.log(np.array(rates[-60:])+1e-300)))) / rho
    if N >= 27:
        Nl3 = Nl // 3
        M_mat = np.stack([v[:Nl], v[Nl:2*Nl], v[2*Nl:]])
        level  = M_mat.mean(axis=0)[:Nl3]
        spread = M_mat.std(axis=0)[:Nl3]
        ok = spread > 1e-12
        corr_ls = float(np.corrcoef(np.log(level[ok]+1e-12), np.log(spread[ok]+1e-12))[0,1]) if ok.sum()>10 else float("nan")
    else:
        corr_ls = float("nan")
    return rho, cv, sw_rho, corr_ls

print("Methode 8: eps-perturbatie VOLLEDIG (lam=1.70, k=12)")
print("=" * 65)
print("  eps    rho      CV(v)    sw/rho   corr(logL,logS)")
k = 12
for eps in np.linspace(0.0, 1.0, 11):
    rho, cv, sw_rho, corr_ls = perron_and_stats(k, LAM, eps)
    corr_str = "nan" if np.isnan(corr_ls) else f"{corr_ls:.5f}"
    print(f"  {eps:.2f}   {rho:.5f}  {cv:.5f}  {sw_rho:.5f}  {corr_str}")

_, cv_hi, _, _ = perron_and_stats(k, LAM, 0.98)
_, cv_lo, _, _ = perron_and_stats(k, LAM, 0.92)
print(f"  d(CV)/d(eps) bij eps~1: {(cv_hi-cv_lo)/0.06:.4f}")
print("done")
