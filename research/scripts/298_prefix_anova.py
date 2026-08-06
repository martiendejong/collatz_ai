import numpy as np
from math import log2

ALPHA = log2(3.0)

for lam in [1.05, 1.70]:
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    k = 13
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr==0; m2 = r_arr==2
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    del s_arr, r_arr
    v = np.ones(N)
    for _ in range(400):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w_ = A*v[T4]
        w_[m2] += B3*cb[R3[m2]]
        w_[m0] += B1*cb[R1[m0]]
        v = w_ / w_.max()
    F = np.log2(v); F -= F.mean()
    total = F.var()

    # C(p) = Var(E[F | i mod 3^(p+1)]) : condition on the p+1 least-significant digits
    print(f"lam={lam}: total Var = {total:.5f}")
    prev = 0.0
    incs = []
    for p in range(k-1):
        M = 3**(p+1)
        cm = F.reshape(N//M, M).mean(axis=0)
        C = (cm**2).mean()          # variance of conditional mean (F centered)
        inc = C - prev
        incs.append(inc)
        prev = C
    print("  C(p):     ", " ".join(f"{x:.5f}" for x in np.cumsum(incs)))
    print("  increments:", " ".join(f"{x:.2e}" for x in incs))
    rats = [incs[p+1]/incs[p] if incs[p] > 1e-12 else float('nan') for p in range(len(incs)-1)]
    print("  inc ratios:", " ".join(f"{r:.3f}" for r in rats))
    print(f"  coverage C(max)/total = {prev/total:.4f}")
