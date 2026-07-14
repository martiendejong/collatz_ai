"""R334-340: IDENTIFYING THE g-FIELD.
Candidates for the universal correction field g (residual shape):
(A) subdominant LEFT eigenvectors of the roulette transition matrix at mod 3^7
    (the slowest-decaying perturbation modes of the chain);
(B) simple local functions: f1(r) = 1[first preimage is spring],
    f2(r) = v3-type offsets, f3(r) = log of class value proxy.
Correlate g (= res_k13 centered) against each; report R^2, and check whether
lambda's decay ratio matches |eig2(P)| in the appropriate power."""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

def build_P(j):
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
    return states, P

j = 7; M = 3 ** j
states, P = build_P(j)
sidx = {s: i for i, s in enumerate(states)}
w_, vl = np.linalg.eig(P.T)
order = np.argsort(-np.abs(w_))
print("roulette chain spectrum |lambda_i|:", [f"{abs(w_[order[i]]):.4f}" for i in range(6)])

# stationary + g field
pi = np.real(vl[:, order[0]]); pi /= pi.sum()
C = np.load("certificates/cert_k13.npy", mmap_mode="r")
k = 13; N = 3 ** (k - 1)
ev = np.asarray(C[:N], dtype=np.float64)
B = M // 3
ii = np.arange(N, dtype=np.int64) % B
bm = np.bincount(ii, weights=ev, minlength=B) / np.bincount(ii, minlength=B)
coset = np.array(sorted(s for s in states if s % 3 == 2))
t = np.array([pi[sidx[s]] for s in coset]); t /= t.sum()
e = bm[(coset - 2) // 3]; e /= e.sum()
g = e / t; g = g / g.mean() - 1.0   # centered residual shape

print(f"\ng-field: {len(g)} classes, std {g.std():.4f}")
# (A) project onto top eigenmodes (restricted to coset, real parts)
print("correlation of g with subdominant left-eigenmodes:")
for i in range(1, 8):
    mode = np.real(vl[:, order[i]])
    mode_c = np.array([mode[sidx[s]] for s in coset])
    mode_c = mode_c / np.array([pi[sidx[s]] for s in coset])  # density perturbation shape
    mode_c -= mode_c.mean()
    if mode_c.std() < 1e-12: continue
    r = np.corrcoef(g, mode_c)[0, 1]
    print(f"  mode {i} (|lam|={abs(w_[order[i]]):.4f}): r = {r:+.4f}")

# combined R^2 of top-5 modes (least squares)
X = []
for i in range(1, 12):
    mode = np.real(vl[:, order[i]]); mo = np.array([mode[sidx[s]] for s in coset])
    mo = mo / np.array([pi[sidx[s]] for s in coset]); mo -= mo.mean()
    if mo.std() > 1e-12: X.append(mo)
    mode = np.imag(vl[:, order[i]]); mo = np.array([mode[sidx[s]] for s in coset])
    mo = mo / np.array([pi[sidx[s]] for s in coset]); mo -= mo.mean()
    if mo.std() > 1e-12: X.append(mo)
X = np.array(X).T
coef, res_, *_ = np.linalg.lstsq(X, g, rcond=None)
pred = X @ coef
r2 = 1 - ((g - pred)**2).sum()/ (g**2).sum()
print(f"\ncombined R^2 of top eigenmode space (dim {X.shape[1]}): {r2:.4f}")

# (B) simple local predictors
f_spring = np.array([1.0 if ((2*s - 1) % 3 == 0 and ((2*s - 1)//3) % 3 == 0) else 0.0 for s in coset])
f_spring -= f_spring.mean()
r_sp = np.corrcoef(g, f_spring)[0, 1] if f_spring.std() > 0 else 0
print(f"corr(g, first-preimage-is-spring) = {r_sp:+.4f}")
# lambda decay vs |lam2|: lambda ratio per 2 digits measured 0.874-0.887
print(f"\n|lam2(P)| = {abs(w_[order[1]]):.4f}; |lam2|^2 = {abs(w_[order[1]])**2:.4f} "
      f"[measured lambda ratio/2digits ~ 0.874-0.887]")
