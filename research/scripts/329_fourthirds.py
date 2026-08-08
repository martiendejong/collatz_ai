"""
Preregistered: fixed-point excess law delta(1,j) - A = D*(4/3)^j.
Predictions for mod 729 (frozen): delta(1)-A = 0.000227 / 0.0355 / 0.0512
for lam = 1.05 / 1.70 / 1.90.
Also: (b) AR coefficient of richness along the orbit ordering;
(c) sd of the richness field per level (growth rate).
"""
import numpy as np
from math import log2
CACHE = "E:/projects/collatz/research/cache"
ALPHA = log2(3.0)

def get_v(lam, k):
    import os
    fn = f"{CACHE}/v_lam{lam:.2f}_k{k}.npy"
    if os.path.exists(fn):
        return np.load(fn)
    A = lam**-2; B1 = lam**(ALPHA-2); B3 = lam**(ALPHA-1)
    N = 3**(k-1); Nl = N//3
    i = np.arange(N, dtype=np.int64)
    T4 = (4*i+2) % N
    s_arr, r_arr = np.divmod(i, 3)
    m0 = r_arr==0; m2 = r_arr==2
    R1 = (4*s_arr) % Nl; R3 = (2*s_arr+1) % Nl
    v = np.ones(N)
    for _ in range(250):
        cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
        w = A*v[T4]; w[m2] += B3*cb[R3[m2]]; w[m0] += B1*cb[R1[m0]]
        v = w/w.max()
    np.save(fn, v)
    return v

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

for lam, k in [(1.05, 14), (1.70, 14), (1.90, 14)]:
    v = get_v(lam, k)
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    vbar = (v[:Nl]+v[Nl:2*Nl]+v[2*Nl:])/3.0
    G = vbar - cb
    print(f"\n=== lam={lam} (k={k}) ===")
    # (a) fixed-point excess across moduli, incl. 729
    row = []
    for j, M in [(2,9),(3,27),(4,81),(5,243),(6,729)]:
        su = np.arange(Nl) % M
        rich_u = np.array([G[su == u].sum()/vbar[su == u].sum() for u in range(M)])
        orbs = orbits_mod(M)
        # A = weighted mean over the LONGEST orbit (homogenized reference)
        big = max(orbs, key=len)
        Aref = rich_u[big].mean()
        d1 = rich_u[M-1]   # fixed point u = -1 mod M
        row.append((M, d1 - Aref))
    print("  delta(1) - A per modulus:", " ".join(f"m{M}:{x:.5f}" for M, x in row))
    rats = [row[i+1][1]/row[i][1] for i in range(len(row)-1)]
    print("  groeiratio's:", " ".join(f"{r:.3f}" for r in rats), " (4/3 = 1.333)")
    # (b) AR coefficient along the longest orbit at mod 243
    M = 243
    su = np.arange(Nl) % M
    rich_u = np.array([G[su == u].sum()/vbar[su == u].sum() for u in range(M)])
    orbs = orbits_mod(M)
    big = max(orbs, key=len)
    seq = rich_u[big] - rich_u[big].mean()
    ar1 = float(np.corrcoef(seq[1:], seq[:-1])[0,1])
    ar2 = float(np.corrcoef(seq[2:], seq[:-2])[0,1])
    print(f"  AR langs orbit (L={len(big)}): lag1={ar1:+.3f} lag2={ar2:+.3f}")
    # (c) sd of richness field per level
    sds = []
    for j, M2 in [(2,9),(3,27),(4,81),(5,243),(6,729)]:
        su2 = np.arange(Nl) % M2
        r2 = np.array([G[su2 == u].sum()/vbar[su2 == u].sum() for u in range(M2)])
        sds.append(r2.std())
    print("  sd(rich)/niveau:", " ".join(f"{s:.5f}" for s in sds),
          "| ratio's:", " ".join(f"{sds[i+1]/sds[i]:.3f}" for i in range(len(sds)-1)))
