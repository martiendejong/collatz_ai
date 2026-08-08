"""
The fundamental mixing function: autocorrelation of the log-field along the
single T4-cycle. ACF(n) = Corr(F(i), F(T4^n i)) over the whole cycle.
This controls the cross-terms in inc_p(K w) — the key to the contraction proof.
Expect: 3-periodic structure (class routing) x geometric decay (mixing).
"""
import numpy as np
from math import log2
CACHE = "E:/projects/collatz/research/cache"
ALPHA = log2(3.0)

for lam, k in [(1.05, 12), (1.70, 12), (2.00, 12)]:
    import os
    fn = f"{CACHE}/v_lam{lam:.2f}_k{k}.npy"
    if os.path.exists(fn):
        v = np.load(fn)
    else:
        A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
        N = 3**(k-1); Nl = N//3
        i = np.arange(N, dtype=np.int64)
        T4 = (4*i+2) % N
        s_arr, r_arr = np.divmod(i, 3)
        m0 = r_arr==0; m2 = r_arr==2
        R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
        v = np.ones(N)
        for _ in range(900):
            cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
            w = A*v[T4]; w[m2] += B3*cb[R3[m2]]; w[m0] += B1*cb[R1[m0]]
            v = w/w.max()
        np.save(fn, v)
    N = v.size
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    F = np.log2(v); F -= F.mean()
    var = float((F**2).mean())
    # ACF along the T4 map: Corr(F(i), F(T4^n i)) — computed via repeated permutation
    idx = i.copy()
    print(f"\n=== lam={lam} k={k} ===  Var={var:.5f}")
    line = []
    acfs = []
    for n in range(1, 31):
        idx = T4[idx]
        acf = float((F*F[idx]).mean())/var
        acfs.append(acf)
        line.append(f"{n}:{acf:+.4f}")
        if n % 10 == 0:
            print("  " + " ".join(line)); line = []
    if line: print("  " + " ".join(line))
    # decay of |ACF| at multiples of 3 (same-class lags)
    sub3 = [abs(acfs[n-1]) for n in [3,6,9,12,15,18,21,24,27,30]]
    rats = [sub3[j+1]/sub3[j] for j in range(len(sub3)-1) if sub3[j] > 1e-6]
    print("  |ACF| bij lags 3,6,...,30:", " ".join(f"{x:.4f}" for x in sub3))
    print("  vervalratio's:", " ".join(f"{r:.3f}" for r in rats))
