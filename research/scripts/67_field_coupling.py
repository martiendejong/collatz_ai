"""R631-640: EFFECTIVE FIELD COUPLING. Regress Delta_P(j) on same-base
neighbor-scale differences D_{P-1}(j), D_{P+1}(j) (and P+-2): the coefficients
the PROFILE recurrence actually reflects. Does a < c live here?"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
k = 13
C = np.load("certificates/cert_k13.npy").astype(np.float64)
C /= C.mean()
N = 3 ** (k - 1)
rng = np.random.default_rng(631)

for P in (4, 5, 6, 7):
    d = 2 * 3 ** P
    js = rng.integers(0, N - 20 * 3 ** P, 100000)
    def D(scale_off):
        return C[(js + scale_off) % N] - C[js]
    LHS = D(d)
    Xdn1 = D(d // 3); Xup1 = D(3 * d)
    Xdn2 = D(d // 9); Xup2 = D(9 * d)
    X = np.column_stack([Xdn1, Xup1, Xdn2, Xup2, np.ones(len(js))])
    coef, *_ = np.linalg.lstsq(X, LHS, rcond=None)
    pred = X @ coef
    r2 = 1 - ((LHS - pred)**2).sum()/((LHS - LHS.mean())**2).sum()
    a1, c1, a2, c2 = coef[:4]
    # amplitude-normalized (per unit CV)
    sa1, sc1 = a1 * Xdn1.std()/LHS.std(), c1 * Xup1.std()/LHS.std()
    print(f"P={P}: a1={a1:+.4f} c1={c1:+.4f} (norm a*={sa1:+.4f} c*={sc1:+.4f}) "
          f"a2={a2:+.4f} c2={c2:+.4f}  R2={r2:.3f}  a*<c*: {sa1 < sc1}")
