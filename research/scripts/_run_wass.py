import numpy as np
from math import log2

ALPHA = log2(3.0)
LAM   = 1.70
A  = LAM ** -2.0
B1 = LAM ** (ALPHA - 2.0)
B3 = LAM ** (ALPHA - 1.0)

def perron(k, n_iter=300):
    N  = 3 ** (k - 1)
    i  = np.arange(N, dtype=np.int64)
    T4 = (4 * i + 2) % N
    s, r = np.divmod(i, 3)
    Nl = N // 3
    m0, m2 = (r == 0), (r == 2)
    R1 = (4 * s) % Nl
    R3 = (2 * s + 1) % Nl
    v = np.ones(N, dtype=np.float64)
    for _ in range(n_iter):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w  = A * v[T4]
        w[m2] += B3 * cb[R3[m2]]
        w[m0] += B1 * cb[R1[m0]]
        v = w / w.max()
    v /= v.mean()
    return v

print("Wasserstein Parts B+C  (lam=1.70)")
print("="*60)
print()
print("(B) KL-divergentie per k:")
kls = []
for k in (10, 11, 12, 13):
    v = perron(k)
    mu = v / v.sum()
    kl = float(np.sum(mu * np.log(mu * len(mu) + 1e-300)))
    kls.append(kl)
    print(f"  k={k}: KL(v^k || uniform) = {kl:.6f} nats")
print()
print("  KL ratios (k+1 / k):")
for i in range(len(kls)-1):
    print(f"  k={10+i}->{11+i}: {kls[i+1]/kls[i]:.5f}")

print()
print("(C) Gecoarsende W1 krimp:")
for k in (10, 11, 12):
    v_k  = perron(k)
    v_k1 = perron(k+1)
    N = v_k.size
    v_c = v_k1[:3*N].reshape(N, 3).mean(axis=1)
    v_c /= v_c.mean()
    diff = float(np.abs(np.sort(v_k) - np.sort(v_c)).mean())
    pw   = float(np.abs(v_k - v_c).max())
    print(f"  k={k}: sorted-W1={diff:.6f}  max-pointwise={pw:.6f}")

print()
print("done")
