import numpy as np
from math import log2

ALPHA = log2(3.0)
lam = 1.05
A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
k = 12
N = 3**(k-1); Nl = N//3
i = np.arange(N, dtype=np.int64)
T4 = (4*i+2) % N
s_arr, r_arr = np.divmod(i, 3)
m0 = r_arr==0; m2 = r_arr==2
R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
del s_arr, r_arr

v = np.ones(N)
res = {}
it_done = 0
for target in [300, 600, 1200, 2400]:
    for _ in range(target - it_done):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w_ = A*v[T4]
        w_[m2] += B3*cb[R3[m2]]
        w_[m0] += B1*cb[R1[m0]]
        v = w_ / w_.max()
    it_done = target
    F = np.log2(v); F -= F.mean()
    incs = []
    prev = 0.0
    for p in range(k-1):
        M = 3**(p+1)
        cm = F.reshape(N//M, M).mean(axis=0)
        C = (cm**2).mean()
        incs.append(C - prev); prev = C
    res[target] = np.array(incs)
    rats = [incs[p+1]/incs[p] for p in range(len(incs)-1)]
    print(f"iters={target:4d}: deep incs: " + " ".join(f"{x:.3e}" for x in incs[6:]) +
          " | ratios: " + " ".join(f"{r:.3f}" for r in rats[5:]), flush=True)
# stability check
for a, b in [(300, 2400), (1200, 2400)]:
    rel = np.abs(res[a]-res[b])/res[2400]
    print(f"max rel diff iters {a} vs 2400: {rel.max():.2e}")
