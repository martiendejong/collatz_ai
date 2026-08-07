"""Cross-validation at digit-3 (Z/27): cb-map u->2u+1 mod 27 orbits are
18-cycle, 6-cycle, 2-cycle (8,17), fixed (26).
PREREGISTERED prediction: richness ordering fixed > 2-cycle > 6-cycle > 18-cycle."""
import numpy as np
CACHE = "E:/projects/collatz/research/cache"

# orbit decomposition of u -> 2u+1 mod 27
p = [(2*u+1) % 27 for u in range(27)]
seen = set(); orbits = []
for s0 in range(27):
    if s0 in seen: continue
    c = [s0]; seen.add(s0); x = p[s0]
    while x != s0:
        c.append(x); seen.add(x); x = p[x]
    orbits.append(c)
orbits.sort(key=len)
print("orbits:", [(len(o), o[:4]) for o in orbits])

for lam in [1.05, 1.70, 1.90]:
    k = 13
    v = np.load(f"{CACHE}/v_lam{lam:.2f}_k{k}.npy")
    N = v.size; Nl = N//3
    cb = np.minimum(np.minimum(v[:Nl], v[Nl:2*Nl]), v[2*Nl:])
    vbar = (v[:Nl]+v[Nl:2*Nl]+v[2*Nl:])/3.0
    G = vbar - cb
    su = np.arange(Nl) % 27
    rich = np.array([G[su == u].mean()/vbar[su == u].mean() for u in range(27)])
    print(f"\nlam={lam}:")
    for o in orbits:
        vals = rich[o]
        print(f"  len={len(o):2d}: richness = {vals.mean():.5f} (sd {vals.std():.5f})")
    print(f"  global: {G.mean()/vbar.mean():.5f}")
