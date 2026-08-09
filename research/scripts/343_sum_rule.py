"""
Hierarchical sum rule test (k=10):
v(i) = (1/rho) sum_n t^n B(cls_n(i)) w(arg_n(i)).
(a) Measure the covariance matrix C(n,m) = Cov(w o arg_n, w o arg_m) for n,m <= nmax
    and test the valuation law: C(n,m) ~ phi_w(v3(n-m)) * Var-scale.
(b) Reconstruct Var(v) from the phi-approximation and compare with true Var(v):
    the quantitative closure of "valuation sums replace mixing".
"""
import numpy as np
from math import log2
CACHE = "E:/projects/collatz/research/cache"
ALPHA = log2(3.0)

def v3(n):
    if n == 0: return 99
    c = 0
    while n % 3 == 0:
        n //= 3; c += 1
    return c

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
    wc = w - w.mean()
    varw = float((wc**2).mean())

    nmax = min(int(np.ceil(np.log(1e-10)/np.log(t))), 60)
    # per step n: coefficient field b_n(i) in {B1, 0, B3}, index field arg_n(i)
    idx = i.copy()
    coeffs = []; args = []
    for n in range(nmax):
        cls = idx % 3
        b = np.where(cls == 0, B1, np.where(cls == 2, B3, 0.0))
        s_of = idx // 3
        ar = np.where(cls == 0, R1[idx], np.where(cls == 2, R3[idx], 0)).astype(np.int64)
        coeffs.append(b); args.append(ar)
        idx = T4[idx]

    # (a) covariances grouped by v3(n-m): sample pairs
    groups = {}
    for n in range(nmax):
        for m in range(n, min(n+28, nmax)):
            sel = (coeffs[n] > 0) & (coeffs[m] > 0)
            if sel.sum() < 100: continue
            cnm = float((wc[args[n][sel]]*wc[args[m][sel]]).mean())/varw
            d = m - n
            groups.setdefault(v3(d) if d > 0 else 'zero', []).append(cnm)
    print(f"\nlam={lam} k=10: genormaliseerde Cov per v3(n-m):")
    for key in sorted(groups, key=lambda x: (isinstance(x, str), x)):
        g = groups[key]
        print(f"  v3={key}: mean {np.mean(g):+.4f} sd {np.std(g):.4f} (n={len(g)})")

    # (b) closure: Var(v) predicted from phi-approx vs true
    phi = {key: np.mean(g) for key, g in groups.items()}
    var_pred = 0.0
    for n in range(nmax):
        for m in range(nmax):
            d = abs(m - n)
            key = 'zero' if d == 0 else v3(d)
            if key not in phi: continue
            # mean coefficient product over i (only both-active rows)
            cp = float((coeffs[n]*coeffs[m]).mean())
            var_pred += (t**n)*(t**m)*cp*phi[key]*varw/(rho**2)
    # true Var(v) minus the mean-structure part: Var(v) computed on centered w-sum:
    # v - E-structure: use directly Var(v_centered by class means)? simplest: total var of v around its mean has also mean-coefficient cross-terms; instead verify on the centered field:
    v_fluct = (1.0/rho)*sum((t**n)*coeffs[n]*wc[args[n]] for n in range(nmax))
    var_true = float((v_fluct**2).mean())
    print(f"  somregel-closure: Var_pred/Var_true = {var_pred/var_true:.4f}")
