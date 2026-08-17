# 387: the halving law decomposed in the T-frame (accelerated map).
# Terras: parity p_j = n_j XOR f_j(n_0..n_{j-1}) (triangular bijection).
# Measure per T-step j: max linear correlation (single bits + XOR pairs).
# If the T-frame law is clean, the Syracuse pairs (Obs 573) fall out of it.
import numpy as np
rng = np.random.default_rng(387)
M = 400000
n0 = (rng.integers(0, 2**62, M, dtype=np.int64) << 1) | 1
B = 22
bits = np.array([((n0 >> b) & 1) for b in range(B)], dtype=np.int8)
n = n0.copy()
print("T-stap j | max|corr| lineair | dominante voorspeller")
for j in range(1, 17):
    odd = (n & 1).astype(bool)
    n = np.where(odd, 3*n + 1, n) >> 1
    p = (n & 1).astype(np.float64)
    pc = p - p.mean(); sd = pc.std()
    best, who = 0.0, ""
    for a in range(B):
        xc = bits[a] - bits[a].mean()
        c = abs(float((xc*pc).mean())/(xc.std()*sd + 1e-12))
        if c > best: best, who = c, f"bit {a}"
    for a in range(B):
        for b2 in range(a+1, B):
            x = bits[a] ^ bits[b2]
            xc = x - x.mean()
            c = abs(float((xc*pc).mean())/(xc.std()*sd + 1e-12))
            if c > best: best, who = c, f"bit{a}^bit{b2}"
    print(f"{j:>8} | {best:.4f}          | {who}", flush=True)
