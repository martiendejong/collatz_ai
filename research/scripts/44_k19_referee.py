"""R330-333: REFEREE TEST on cert_k19 (387,420,489 entries, memmap chunked).
Prediction from the factorization law (R323-329): residual CV at k=19 should be
~ 0.150 * 0.874 = 0.131, and the residual shape should correlate 0.999 with g."""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

def theory_stationary(j):
    M = 3 ** j
    states = [r for r in range(M) if r % 3 != 0]
    idx = {r: i for i, r in enumerate(states)}
    P = np.zeros((len(states), len(states)))
    inv2 = pow((M + 1) // 2, 1, M)
    for r in states:
        b = (3 * r + 1) % M
        p = 0.5; x = b
        for w in range(1, 100):
            x = (x * inv2) % M
            P[idx[r], idx[x]] += p
            p *= 0.5
    P /= P.sum(1, keepdims=True)
    v = np.ones(len(states)) / len(states)
    for _ in range(4000): v = v @ P
    return dict(zip(states, v))

j = 7; M = 3 ** j
th = theory_stationary(j)
coset = np.array(sorted(s for s in th if s % 3 == 2))
t = np.array([th[s] for s in coset]); t /= t.sum()

def residual(k, path):
    C = np.load(path, mmap_mode="r")
    N = 3 ** (k - 1)
    B = M // 3
    sums = np.zeros(B); cnts = np.zeros(B)
    step = 3 ** 12
    for lo in range(0, N, step):
        hi = min(lo + step, N)
        ev = np.asarray(C[lo:hi], dtype=np.float64)
        ii = (np.arange(lo, hi, dtype=np.int64)) % B
        sums += np.bincount(ii, weights=ev, minlength=B)
        cnts += np.bincount(ii, minlength=B)
    bm = sums / cnts
    e = bm[(coset - 2) // 3]; e /= e.sum()
    return e / t

r13 = residual(13, "certificates/cert_k13.npy")
r19 = residual(19, "certificates/cert_k19.npy")
cv19 = r19.std()/r19.mean()
print(f"k=19 residual CV = {cv19:.4f}   [prediction 0.131 = 0.150 x 0.874]")
print(f"corr(res_k19, res_k13 shape g) = {np.corrcoef(r19, r13)[0,1]:.4f}   [prediction ~0.999]")
# lambda sequence: CVs 13,15,17,19
print(f"lambda ratio k17->k19: {cv19/0.1502:.4f}   [prediction 0.874]")
