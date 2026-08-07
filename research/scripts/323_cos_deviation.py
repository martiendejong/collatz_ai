"""Sharp test: 1 - cos(cb-profile, vbar-profile) per digit and lambda,
benchmarked against gamma^2 (naive rotation scale) and gamma-spread^2.
If 1-cos << gamma^2: exact-cancellation mechanism (H1 at digit level holds
beyond naive scaling) -> provable-shaped invariant."""
import numpy as np
from math import log2
CACHE = "E:/projects/collatz/research/cache"

for lam in [1.05, 1.30, 1.70, 1.90]:
    k = 13
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    vbar = (v[:Nl]+v[Nl:2*Nl]+v[2*Nl:])/3.0
    G = vbar - cb
    gam = G.mean()/vbar.mean()
    j = np.arange(Nl, dtype=np.int64)
    P = k-2
    digs = np.empty((P, Nl), dtype=np.int8)
    x = j.copy()
    for p in range(P):
        digs[p] = x % 3; x //= 3
    print(f"\nlam={lam}: gamma={gam:.5f}  gamma^2={gam**2:.2e}")
    print("p | 1-cos(cb,vbar) | 1-cos(G,vbar) | gamma-profile-ratio spread")
    for p in range(min(P, 6)):
        pc = np.array([cb[digs[p]==d].mean() for d in range(3)]); pc -= pc.mean()
        pv = np.array([vbar[digs[p]==d].mean() for d in range(3)]); pv -= pv.mean()
        pg = np.array([G[digs[p]==d].mean() for d in range(3)]); pg -= pg.mean()
        def onemcos(a, b):
            na, nb = np.linalg.norm(a), np.linalg.norm(b)
            return 1 - a@b/(na*nb) if na>0 and nb>0 else np.nan
        # per-digit relative gap ratios (richness profile)
        gd = np.array([G[digs[p]==d].mean()/vbar[digs[p]==d].mean() for d in range(3)])
        spread = gd.max()-gd.min()
        print(f"{p} | {onemcos(pc,pv):.2e} | {onemcos(pg,pv):.2e} | {spread:.2e}")
