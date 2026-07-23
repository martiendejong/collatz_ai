"""
168_cycle_backbone.py
======================
THE LOCAL MODEL AT -4: single-cycle backbone + self-loop at -1.

Theory (derived on paper, to be verified here):
  In u-coordinates (m = -u mod 3^k, u in 1+3Z_3), the lambda^-2 backbone
  m -> 4m is the single N-cycle u_j = 4^j (4 topologically generates 1+3Z_3;
  ord_{3^k}(4) = 3^{k-1} = N). Information flows v(j) <- v(j+1).
  Special algebra:
    * m = -1 is D3 (m == 8 mod 9) with (2m-1)/3 = -1: A SELF-LOOP through
      the advanced coefficient lambda^(alpha-1) > 1. Prediction: global max.
    * 4*(-1) = -4: the cycle successor of the peak; m = -4 is D2 (pure decay,
      no side-feed). Prediction: global min sits at -4 BECAUSE it is the
      pure-decay successor of the peak (v(-1) large does NOT feed v(-4);
      flow is the other way).
Verify:
  V1: the x4-map on compact indices is a single N-cycle.
  V2: argmax v = class -1, argmin v = class -4, at every k available.
  V3: the profile v(j) along the cycle: decay/recovery lengths; extract the
      effective recovery scale and compare with 1-gamma(k) ~ k^-0.85.
"""
import numpy as np
from math import log2

def compact_of_m(m, k):
    m = m % 3**k
    assert m % 3 == 2
    return (m - 2) // 3

print("V1: single-cycle check (small k):")
for k in [3, 5, 7]:
    N = 3**(k-1)
    seen = np.zeros(N, bool)
    i = 0; L = 0
    while not seen[i]:
        seen[i] = True
        i = (4*i + 2) % N
        L += 1
    print(f"  k={k}: cycle through i=0 has length {L} of N={N}: "
          f"{'SINGLE CYCLE' if L == N else 'MULTIPLE CYCLES'}")
print()

print("V2: peak and trough positions:")
for k in [12, 13, 14, 15, 16, 17]:
    try:
        v = np.load(f"certificate_k{k}.npy").astype(np.float64)
    except FileNotFoundError:
        continue
    N = 3**(k-1)
    im1 = compact_of_m(-1, k)
    im4 = compact_of_m(-4, k)
    amax, amin = int(np.argmax(v)), int(np.argmin(v))
    print(f"  k={k}: argmax={'m=-1' if amax==im1 else f'i={amax} (m={3*amax+2})'}"
          f"  argmin={'m=-4' if amin==im4 else f'i={amin} (m={3*amin+2})'}"
          f"  v(-1)/v(-4) = {v[im1]/v[im4]:.1f}")
print()

print("V3: profile along the x4-cycle from m=-1 (k=15):")
k = 15
v = np.load(f"certificate_k{k}.npy").astype(np.float64); v /= v.max()
N = 3**(k-1)
M = 3**k
# walk the cycle in the DEPENDENCE direction: position j has m_j = -4^j
u = 1
prof = []
for j in range(0, 61):
    m = (-u) % M
    prof.append(v[compact_of_m(m, k)])
    u = (4*u) % M
print("  j : v(j)  (j=0 is m=-1, j=1 is m=-4, ...)")
for j in range(0, 31, 1):
    bar = '#' * max(1, int(48 + 8*log2(prof[j])))
    print(f"  {j:>3}: {prof[j]:.3e}  {bar}")
rat = [prof[j]/prof[j+1] for j in range(30)]
print("  ratios v(j)/v(j+1):", np.round(rat[:12], 3))
print()
# recovery length: how many steps until v returns to median level?
med = float(np.median(v))
back = next(j for j in range(1, 61) if prof[j] >= med)
print(f"  median v = {med:.3e}; recovery to median at j = {back}")
print()
print("V3b: recovery length vs k (steps from the trough back to median):")
for k in [12, 13, 14, 15, 16, 17]:
    try:
        v = np.load(f"certificate_k{k}.npy").astype(np.float64); v /= v.max()
    except FileNotFoundError:
        continue
    M = 3**k
    med = float(np.median(v))
    u = 1; rec = None
    for j in range(0, 200):
        m = (-u) % M
        val = v[(m - 2)//3]
        if j >= 1 and val >= med:
            rec = j; break
        u = (4*u) % M
    print(f"  k={k}: recovery length = {rec}")
