"""R381-394: MULTI-FEATURE g-FIELD REGRESSION. Features per class r (mod 3^7,
r==2 mod 3): spring-blockage of backward rungs w=1..6, caste digits, roulette
one-step-source composition. Target: the universal g field. Report R^2 buildup."""
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
        b = (3 * r + 1) % M; p = 0.5; x = b
        for w in range(1, 100):
            x = (x * inv2) % M
            P[idx[r], idx[x]] += p; p *= 0.5
    P /= P.sum(1, keepdims=True)
    v = np.ones(len(states)) / len(states)
    for _ in range(4000): v = v @ P
    return dict(zip(states, v))

j = 7; M = 3 ** j
th = theory_stationary(j)
coset = np.array(sorted(s for s in th if s % 3 == 2))
t = np.array([th[s] for s in coset]); t /= t.sum()
C = np.load("certificates/cert_k13.npy", mmap_mode="r")
N = 3 ** 12
ev = np.asarray(C[:N], dtype=np.float64)
B = M // 3
ii = np.arange(N, dtype=np.int64) % B
bm = np.bincount(ii, weights=ev, minlength=B) / np.bincount(ii, minlength=B)
e = bm[(coset - 2) // 3]; e /= e.sum()
g = e / t; g = g / g.mean() - 1.0

feats, names = [], []
# backward rung structure: for class r, rung w preimage p_w = (r*2^w - 1)/3 exists iff
# r*2^w == 1 mod 3 (w parity by caste); blocked iff p_w == 0 mod 3 (mod 3^2 info of r)
for w in range(1, 7):
    col = []
    for r in coset:
        x = (r * (1 << w) - 1)
        col.append(1.0 if (x % 3 == 0 and (x // 3) % 3 == 0) else 0.0)
    col = np.array(col)
    if col.std() > 0: feats.append(col - col.mean()); names.append(f"rung{w}_spring")
# 3-adic digits of r (one-hot digits 1..3)
for pos in range(1, 5):
    dig = (coset // 3**pos) % 3
    for val in (0, 1):
        col = (dig == val).astype(float)
        feats.append(col - col.mean()); names.append(f"digit{pos}={val}")
# log-theory (does g correlate with the roulette weight itself?)
lt = np.log(t); feats.append(lt - lt.mean()); names.append("log_theory")

X = np.array(feats).T
print(f"features: {len(names)}; single-feature correlations with g:")
for i, nm in enumerate(names):
    r = np.corrcoef(g, X[:, i])[0, 1]
    if abs(r) > 0.25: print(f"  {nm:>14}: r = {r:+.3f}")
coef, *_ = np.linalg.lstsq(X, g, rcond=None)
pred = X @ coef
r2 = 1 - ((g - pred)**2).sum()/(g**2).sum()
print(f"\nfull linear model R^2 = {r2:.4f}")
# interaction: rung1_spring x log_theory
inter = feats[0] * (lt - lt.mean())
X2 = np.column_stack([X, inter - inter.mean()])
coef2, *_ = np.linalg.lstsq(X2, g, rcond=None)
r22 = 1 - ((g - X2 @ coef2)**2).sum()/(g**2).sum()
print(f"+ rung1xlogtheory interaction: R^2 = {r22:.4f}")
