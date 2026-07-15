"""R621-630: RECONSTRUCT (a,c) from the signed channels. Regress LHS (= Delta_P)
on the pure cross-scale source differences: X_dn = branch-base difference at
scale P-1, X_up = transport-base difference at scale P+1. The regression
coefficients, scaled by the CV-ratios between scales, give the effective
lattice (a, c). Target: (0.47, 0.52) at k=13-ish; a < c must emerge."""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
k = 13
C = np.load("certificates/cert_k13.npy").astype(np.float64)
C /= C.mean()
N = 3 ** (k - 1); M = 3 ** k; Mc = 3 ** (k - 1); M3 = 3 ** (k - 2)
lam = 1.818; A = math.log2(3)
W0 = lam**-2; W2 = lam**(A-2); W8 = lam**(A-1)
rng = np.random.default_rng(621)

def cls(m): return ((m % M) - 2) // 3
def tri_min(b): return np.minimum(np.minimum(C[b % M3], C[b % M3 + M3]), C[b % M3 + 2*M3])

for P in (4, 5, 6):
    d = 2 * 3 ** P
    js = rng.integers(0, N, 80000)
    m = 3 * js + 2
    LHS = C[cls(m + d)] - C[cls(m)]
    # up-source: transport 3delta-part at image
    y = 4 * m
    X_up = C[cls(y + 3*d)] - C[cls(y)]        # scale P+1 difference at transport base
    # dn-source: branch base difference at scale P-1 (offset d/3 at t)
    mod9 = m % 9
    X_dn = np.zeros_like(LHS)
    for mv, mul, sub in ((2, 4, 2), (8, 2, 1)):
        sel = mod9 == mv
        mm = m[sel]
        t = ((mul * mm - sub) // 3) % Mc
        b = (t - 2) // 3
        db = d // 3
        X_dn[sel] = tri_min(((t + db) % Mc - 2) // 3) - tri_min(b)
    X = np.column_stack([X_dn, X_up])
    coef, *_ = np.linalg.lstsq(np.column_stack([X, np.ones(len(LHS))]), LHS, rcond=None)
    a_r, c_r = coef[0], coef[1]
    pred = X @ coef[:2] + coef[2]
    r2 = 1 - ((LHS - pred)**2).sum()/((LHS - LHS.mean())**2).sum()
    # scale to CV-ratio convention: a_eff = a_r * CV_{P-1}/CV_P etc. -- report both raw and CV-scaled
    s_dn = X_dn.std()/LHS.std(); s_up = X_up.std()/LHS.std()
    print(f"P={P}: raw regression a={a_r:.4f} c={c_r:.4f}  (R2 {r2:.3f})"
          f"  amplitude-scaled: a*={a_r*s_dn:.4f} c*={c_r*s_up:.4f}  a<c: {a_r*s_dn < c_r*s_up}")
