import numpy as np
from math import log2

ALPHA = log2(3.0)
lam = 1.70
A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
ITERS = {10: 1200, 11: 1000, 12: 600, 13: 400, 14: 250, 15: 150}

print(f"lam={lam}: per-digit prefix-increment ratios, stability across k")
for k in [10, 11, 12, 13, 14, 15]:
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr==0; m2 = r_arr==2
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    del s_arr, r_arr
    v = np.ones(N)
    for _ in range(ITERS[k]):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w_ = A*v[T4]
        w_[m2] += B3*cb[R3[m2]]
        w_[m0] += B1*cb[R1[m0]]
        v = w_ / w_.max()
    F = np.log2(v); F -= F.mean()
    incs = []
    prev = 0.0
    for p in range(k-1):
        M = 3**(p+1)
        cm = F.reshape(N//M, M).mean(axis=0)
        C = (cm**2).mean()
        incs.append(C - prev); prev = C
    rats = [incs[p+1]/incs[p] for p in range(len(incs)-1)]
    # print ratios aligned by digit position p
    print(f"k={k:2d}: " + " ".join(f"{r:.3f}" for r in rats))
print()
print("column p is the ratio inc(p+1)/inc(p); read down a column for k-stability")
