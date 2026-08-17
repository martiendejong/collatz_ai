# 391: test of the AND-tree mechanism for the halving law (Obs 588).
# Hypothesis: the nonlinearity f_j is an AND-tree growing one level per two
# T-steps; then the best DEGREE-2 predictor (one product term allowed) should
# recover a factor ~2 over the best linear predictor at every j.
import numpy as np
rng = np.random.default_rng(391)
M = 200000
n0 = (rng.integers(0, 2**62, M, dtype=np.int64) << 1) | 1
B = 20
bits = np.array([((n0 >> b) & 1) for b in range(B)], dtype=np.int8)
n = n0.copy()
print("j | best lineair | best +1 AND-term | ratio")
for j in range(1, 15):
    odd = (n & 1).astype(bool)
    n = np.where(odd, 3*n + 1, n) >> 1
    p = (n & 1).astype(np.int8)
    pc = p - p.mean(); sd = pc.std()
    bestL = 0.0
    for a in range(B):
        for b2 in range(a, B):
            x = bits[a] if b2 == a else (bits[a] ^ bits[b2])
            xc = x - x.mean()
            c = abs(float((xc*pc).mean())/(xc.std()*sd + 1e-12))
            bestL = max(bestL, c)
    bestQ = 0.0
    # degree-2: single bit XOR (product of two nearby bits)
    for a in range(B):
        for g1 in range(B):
            for g2 in range(g1+1, min(g1+5, B)):
                x = bits[a] ^ (bits[g1] & bits[g2])
                xc = x - x.mean()
                c = abs(float((xc*pc).mean())/(xc.std()*sd + 1e-12))
                bestQ = max(bestQ, c)
    print(f"{j:>2} | {bestL:.4f}      | {bestQ:.4f}          | {bestQ/max(bestL,1e-9):.2f}", flush=True)
