"""Preregistered test: mod 243 (j=5). Predictions (frozen before this run):
A unchanged per lambda: ~0.00043 / 0.0300 / 0.0519
B(243) = 1.44 * B(81):  ~0.00025 / 0.0299 / 0.0434
Candidate law: B proportional to 3^(j/3)."""
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
    M = 243
    su = np.arange(Nl) % M
    rich_u = np.array([G[su == u].sum()/vbar[su == u].sum() for u in range(M)])
    orbs = orbits_mod(M)
    Ls, ds, ws = [], [], []
    for o in orbs:
        Ls.append(len(o)); ds.append(rich_u[o].mean()); ws.append(len(o))
    Ls = np.array(Ls, float); ds = np.array(ds); ws = np.array(ws, float)
    X1 = np.column_stack([np.ones_like(Ls), 1.0/Ls])
    W = np.diag(ws)
    beta = np.linalg.solve(X1.T@W@X1, X1.T@W@ds)
    r = ds - X1@beta
    tot = float((ws*(ds-ds.mean())**2).sum())
    R2 = 1 - float((ws*r**2).sum())/tot
    pairs = " ".join(f"L={int(L)}:{d:.5f}" for L, d in sorted(zip(Ls, ds)))
    print(f"lam={lam} mod 243: {pairs}")
    print(f"  A={beta[0]:.5f} B={beta[1]:.5f} R2={R2:.4f}")
