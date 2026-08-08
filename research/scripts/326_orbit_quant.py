"""
Quantitative orbit law: richness deviation delta(L) vs orbit length L of u->2u+1
mod 3^j, at j=2,3,4 (Z/9, Z/27, Z/81) and lambda = 1.05, 1.70, 1.90.
Candidate laws: delta = A + B/L  (ergodic averaging)  vs  delta = A + B*c^L (recirculation).
Fit quality compared per (lambda, modulus); cross-validated.
"""
import numpy as np

CACHE = "E:/projects/collatz/research/cache"

def orbits_mod(M):
    p = [(2*u+1) % M for u in range(M)]
    seen = set(); out = []
    for s0 in range(M):
        if s0 in seen: continue
        c = [s0]; seen.add(s0); x = p[s0]
        while x != s0:
            c.append(x); seen.add(x); x = p[x]
        out.append(c)
    return out

for lam in [1.05, 1.70, 1.90]:
    k = 13
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    vbar = (v[:Nl]+v[Nl:2*Nl]+v[2*Nl:])/3.0
    G = vbar - cb
    print(f"\n=== lam={lam} ===")
    for j, M in [(2, 9), (3, 27), (4, 81)]:
        su = np.arange(Nl) % M
        rich_u = np.array([G[su == u].sum()/vbar[su == u].sum() for u in range(M)])
        orbs = orbits_mod(M)
        Ls, ds, ws = [], [], []
        for o in orbs:
            Ls.append(len(o)); ds.append(rich_u[o].mean()); ws.append(len(o))
        Ls = np.array(Ls, float); ds = np.array(ds); ws = np.array(ws, float)
        # fit A + B/L (weighted by orbit size)
        X1 = np.column_stack([np.ones_like(Ls), 1.0/Ls])
        W = np.diag(ws)
        beta1 = np.linalg.solve(X1.T@W@X1, X1.T@W@ds)
        r1 = ds - X1@beta1
        sse1 = float((ws*r1**2).sum())
        # fit A + B*c^L with c from our measurements (lambda-dependent contraction)
        cmap = {1.05: 0.41, 1.70: 0.70, 1.90: 0.78}
        cl = cmap[lam]**Ls
        X2 = np.column_stack([np.ones_like(Ls), cl])
        beta2 = np.linalg.solve(X2.T@W@X2, X2.T@W@ds)
        r2 = ds - X2@beta2
        sse2 = float((ws*r2**2).sum())
        tot = float((ws*(ds-ds.mean())**2).sum())
        pairs = " ".join(f"L={int(L)}:{d:.5f}" for L, d in sorted(zip(Ls, ds)))
        print(f" mod {M}: {pairs}")
        print(f"   1/L-fit: A={beta1[0]:.5f} B={beta1[1]:.5f} R2={1-sse1/tot:.4f} | "
              f"c^L-fit: R2={1-sse2/tot:.4f}")
