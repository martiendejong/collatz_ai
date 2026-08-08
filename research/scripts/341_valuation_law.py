"""Valuation law test: ACF(n) along the T4-cycle grouped by v3(n) should match
the prefix-variance fractions (digits preserved by T4^n = v3(n)+1, up to shift).
Predict: mean ACF over lags with v3(n)=j ~ monotone increasing in j, tracking
C(j)/Var (prefix fraction)."""
import numpy as np
from math import log2
CACHE = "E:/projects/collatz/research/cache"

def v3(n):
    v = 0
    while n % 3 == 0:
        n //= 3; v += 1
    return v

for lam, k in [(1.05, 12), (1.70, 12), (2.00, 12)]:
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    N = v.size
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    F = np.log2(v); F -= F.mean()
    var = float((F**2).mean())
    # prefix fractions C(p)/var
    Cfr = []
    prev = None
    for p in range(k-1):
        M = 3**(p+1)
        cm = F.reshape(N//M, M).mean(axis=0)
        Cfr.append(float((cm**2).mean())/var)
    # ACF up to lag 90, grouped by v3
    idx = i.copy()
    groups = {}
    for n in range(1, 91):
        idx = T4[idx]
        if n % 3 != 0 and n > 6:  # skip most non-divisible for speed? keep all under 91
            pass
        acf = float((F*F[idx]).mean())/var
        groups.setdefault(v3(n), []).append(acf)
    print(f"\nlam={lam}: prefix-fracties C(p)/Var:", " ".join(f"p{p}:{c:.4f}" for p, c in enumerate(Cfr[:6])))
    for j in sorted(groups):
        g = groups[j]
        print(f"  v3(n)={j}: mean ACF {np.mean(g):+.4f} (sd {np.std(g):.4f}, n={len(g)}) "
              f"| voorspelde bovengrens C({j})/Var = {Cfr[j] if j < len(Cfr) else float('nan'):.4f}")
