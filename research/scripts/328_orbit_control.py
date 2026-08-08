"""Control: random partitions with the same size distribution as the true orbits.
CLT gives 1/L concentration for any subsets; the DYNAMICAL content of the law is the
sign-definite enrichment (B > 0). Prediction: control B ~ 0 within noise."""
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

rng = np.random.default_rng(11)
for lam in [1.70, 1.90]:
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
    sizes = sorted([len(o) for o in orbs])
    def fitAB(partition):
        Ls, ds, ws = [], [], []
        for o in partition:
            Ls.append(len(o)); ds.append(rich_u[list(o)].mean()); ws.append(len(o))
        Ls = np.array(Ls, float); ds = np.array(ds); ws = np.array(ws, float)
        X = np.column_stack([np.ones_like(Ls), 1.0/Ls])
        W = np.diag(ws)
        beta = np.linalg.solve(X.T@W@X, X.T@W@ds)
        return beta
    bt = fitAB(orbs)
    # 200 random partitions with identical size distribution
    Bs = []
    for trial in range(200):
        perm = rng.permutation(M)
        parts = []; pos = 0
        for s in sizes:
            parts.append(perm[pos:pos+s]); pos += s
        Bs.append(fitAB(parts)[1])
    Bs = np.array(Bs)
    p_emp = (Bs >= bt[1]).mean()
    print(f"lam={lam}: true B = {bt[1]:.5f} | control B: mean {Bs.mean():+.5f} sd {Bs.std():.5f} "
          f"| z = {(bt[1]-Bs.mean())/Bs.std():.1f} | empirical p = {p_emp:.3f}")
