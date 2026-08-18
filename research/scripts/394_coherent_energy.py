# 394: the coherent-energy instrument — Parseval splits the Walsh spectrum of
# the parity function into a Haar noise floor (|corr| ~ 2^{-j/2}) and coherent
# structure above it. If max-coefficients are amplitudes (~sqrt(c)), the
# coherent ENERGY should decay at c itself. Exact per depth via full FWHT.
import numpy as np
import math

def fwht_int(a):
    h = 1; n = len(a)
    while h < n:
        a = a.reshape(-1, 2, h)
        x = a[:, 0, :].copy()
        a[:, 0, :] = x + a[:, 1, :]
        a[:, 1, :] = x - a[:, 1, :]
        a = a.reshape(n); h *= 2
    return a

print(f"{'j':>3} {'E_coh (>4x ruis)':>16} {'#coherent':>9} {'max^2':>10}")
E = {}
for j in range(12, 25):
    Nj = 1 << (j+1)
    n = np.arange(Nj, dtype=np.int64)
    for _ in range(j):
        odd = (n & 1).astype(bool)
        n = np.where(odd, (3*n + 1) >> 1, n >> 1)
    s = (1 - 2*(n & 1)[1::2]).astype(np.int64)
    del n
    W = fwht_int(s)
    del s
    L = 1 << j
    corr2 = (W.astype(np.float64)/L)**2
    del W
    corr2[0] = 0.0   # DC eruit
    thr = (4.0/math.sqrt(L))**2
    mask = corr2 > thr
    Ec = float(corr2[mask].sum())
    cnt = int(mask.sum())
    mx2 = float(corr2.max())
    E[j] = Ec
    print(f"{j:>3} {Ec:>16.6f} {cnt:>9} {mx2:>10.6f}", flush=True)
    del corr2
ks = sorted(E)
r = [E[ks[i+1]]/E[ks[i]] for i in range(len(ks)-1)]
print("\nE_coh-ratio's per stap:", ["%.3f" % x for x in r])
tail = r[-6:]
gm = math.exp(sum(math.log(x) for x in tail)/len(tail))
print(f"staart-vervalvoet energie: {gm:.4f}   (kandidaten: c(2) = 0.835 | rate(gammabar) = 0.917 | max^2-rate)")
