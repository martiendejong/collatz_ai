"""Obs 617: the eigen-visit bridge deepened (S5 candidate 7 / agenda item 2).

Claim under test (R285 measured r=0.985 at mod 27): the K-L certificate
weights ARE the forward visit measure under time reversal. Two levels:
  (a) coarse: aggregate cert_k13 mass per residue class mod 81 and 243
      (classes = 2 mod 3) vs measured orbit-visit frequencies.
  (b) pointwise: P(orbit hits v) for every v = 2 mod 3, 5 <= v < 6561,
      vs the certificate value at v's exact class (idx = (v-2)/3),
      including the turnstile prefactors of Obs 616
      (d(341)/d(85) = 1.582, rungs 13/53 near-50/50).
"""
import math
import numpy as np
from collections import Counter

def v2(x):
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return v

C = np.load(r"E:\projects\collatz\research\certificates\cert_k13.npy")
print(f"cert_k13: {C.shape[0]:,} classes (m = 3i+2 mod 3^13)")

# ---- forward measurement: full orbit walks, no memo shortcuts ----
NSTARTS = 300_000
starts = range(1_000_001, 1_000_001 + 2 * NSTARTS, 2)
POINTMAX = 6561  # 3^8
vis81 = Counter(); vis243 = Counter(); hits = Counter()
for n0 in starts:
    n = n0
    while n != 1:
        if n % 3 == 2:
            vis81[n % 81] += 1
            vis243[n % 243] += 1
            if n < POINTMAX:
                hits[n] += 1
        m = 3 * n + 1
        n = m >> v2(m)

def cert_mass(mod):
    """Aggregate certificate mass per residue class mod `mod`."""
    idx = np.arange(C.shape[0], dtype=np.int64)
    m = (3 * idx + 2) % mod
    mass = np.zeros(mod)
    np.add.at(mass, m, C)
    return {r: mass[r] for r in range(mod) if r % 3 == 2}

def corr(xs, ys):
    lx = np.log(np.array(xs)); ly = np.log(np.array(ys))
    return np.corrcoef(lx, ly)[0, 1]

print("\n(a) coarse bridge (log-log correlation cert-mass vs visit-freq):")
for mod, vis in ((81, vis81), (243, vis243)):
    cm = cert_mass(mod)
    common = [r for r in cm if vis.get(r, 0) > 0]
    r = corr([cm[c] for c in common], [vis[c] for c in common])
    print(f"  mod {mod:>3}: {len(common)} classes, r(log) = {r:.4f}")

print("\n(b) pointwise: P(hit v) vs cert value at v's own class:")
vs = sorted(v for v in hits if hits[v] >= 200)
pv = [hits[v] / NSTARTS for v in vs]
cv = [float(C[(v - 2) // 3]) for v in vs]
r = corr(pv, cv)
print(f"  {len(vs)} values v (>=200 hits), corr(log P, log C) = {r:.4f}")
# power-law fit log P = a*log C + b
a, b = np.polyfit(np.log(cv), np.log(pv), 1)
print(f"  fit: P ~ C^{a:.3f}")
print("  key ratios (measured Obs 616 vs certificate):")
for v1, v2_ in ((341, 85 * 4 + 1), (341, 5), (53, 5), (3413, 53)):
    pass
for (va, vb) in ((341, 5), (53, 5), (3413, 53), (21845, 341)):
    ca, cb = float(C[(va - 2) // 3]), float(C[(vb - 2) // 3])
    pa = hits.get(va, 0) / NSTARTS
    pb = hits.get(vb, 0) / NSTARTS
    print(f"    P({va})/P({vb}) = {pa / pb if pb else float('nan'):8.4f}   "
          f"C({va})/C({vb}) = {ca / cb:8.4f}   "
          f"(C-ratio)^{a:.2f} = {(ca / cb) ** a:8.4f}")
